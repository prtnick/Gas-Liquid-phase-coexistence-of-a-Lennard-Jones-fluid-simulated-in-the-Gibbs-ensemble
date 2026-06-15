import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np


# ============================================================
# Input / output
# ============================================================

MEASUREMENTS_FILE = Path("measurements.dat")
MU_FILE = Path("mu_measurements.dat")

OUTPUT_DIR = Path("plots")
OUTPUT_DIR.mkdir(exist_ok=True)

PRESSURE_BLOCK_SIZE = 50
EQUILIBRATION_FRACTION = 0.3

TEMPERATURE = 0.70
MU_BLOCK_SIZE = 10

# ============================================================
# Reading utilities
# ============================================================

def read_table(path, required_columns):
    if not path.exists():
        raise FileNotFoundError(f"Could not find input file: {path}")

    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df.columns = df.columns.str.strip()

    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")

    return df


def block_average(df, columns, block_size):
    blocks = df.index // block_size
    return df.groupby(blocks, as_index=False)[columns].mean()


def mean_and_uncertainty(series):
    n = series.count()
    mean = series.mean()
    uncertainty = series.std(ddof=1) / n**0.5 if n > 1 else float("nan")
    return mean, uncertainty

def chemical_potential_block_analysis(df, weight_column, temperature, block_size):
    data = df[["t", weight_column]].copy()
    data = data.replace([np.inf, -np.inf], np.nan)
    data = data.dropna()
    data = data[data[weight_column] > 0.0]
    data = data.reset_index(drop=True)

    if data.empty:
        return float("nan"), float("nan"), pd.DataFrame(columns=["t", "mu"])

    # Central estimate: logarithm of the full equilibrated average weight
    mu_central = -temperature * np.log(data[weight_column].mean())

    # Block uncertainty: compute one mu estimate per block
    data["block"] = data.index // block_size

    block_data = data.groupby("block").agg(
        t=("t", "mean"),
        W=(weight_column, "mean")
    ).reset_index(drop=True)

    block_data = block_data[block_data["W"] > 0.0]
    block_data["mu"] = -temperature * np.log(block_data["W"])

    n_blocks = len(block_data)

    if n_blocks > 1:
        mu_uncertainty = block_data["mu"].std(ddof=1) / np.sqrt(n_blocks)
    else:
        mu_uncertainty = float("nan")

    return mu_central, mu_uncertainty, block_data

# ============================================================
# Read measurements
# ============================================================

meas_cols = ["t", "E_tot", "n_gas", "n_liq", "V_gas", "V_liq", "P_gas", "P_liq"]
widom_cols = ["t", "W_gas", "W_liq"]

df = read_table(MEASUREMENTS_FILE, meas_cols)
df_mu = read_table(MU_FILE, widom_cols)


# ============================================================
# Derived quantities
# ============================================================

df["rho_gas"] = df["n_gas"] / df["V_gas"]
df["rho_liq"] = df["n_liq"] / df["V_liq"]
df_pressure_blocks = block_average(
    df,
    columns=["t", "P_gas", "P_liq"],
    block_size=PRESSURE_BLOCK_SIZE
)

mc_steps = max(df["t"].max(), df_mu["t"].max())
equilibration_cutoff = EQUILIBRATION_FRACTION * mc_steps
df_eq = df[df["t"] >= equilibration_cutoff]
df_mu_eq = df_mu[df_mu["t"] >= equilibration_cutoff]

if df_eq.empty:
    raise ValueError("No pressure measurements remain after the equilibration cutoff.")

if df_mu_eq.empty:
    raise ValueError("No chemical-potential measurements remain after the equilibration cutoff.")

P_gas_avg, P_gas_unc = mean_and_uncertainty(df_eq["P_gas"])
P_liq_avg, P_liq_unc = mean_and_uncertainty(df_eq["P_liq"])

mu_gas_avg, mu_gas_unc, mu_gas_blocks = chemical_potential_block_analysis(
    df_mu_eq,
    "W_gas",
    TEMPERATURE,
    MU_BLOCK_SIZE
)
mu_liq_avg, mu_liq_unc, mu_liq_blocks = chemical_potential_block_analysis(
    df_mu_eq,
    "W_liq",
    TEMPERATURE,
    MU_BLOCK_SIZE
)


# ============================================================
# Plot style
# ============================================================

plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "legend.fontsize": 11,
    "figure.figsize": (10, 5),
    "lines.linewidth": 1.0,
})


# ============================================================
# Generic gas-liquid plotter
# ============================================================

def plot_gas_liquid(
    t,
    gas,
    liq,
    ylabel,
    title,
    gas_label,
    liq_label,
    output_name
):
    plt.figure()

    # Liquid first, I need it to be behind the gas curve, otherwise fluctuations cover gas data
    plt.plot(t, liq, label=liq_label, linewidth=0.9, alpha=0.75, zorder=1, rasterized=True)
    plt.plot(t, gas, label=gas_label, linewidth=0.9, alpha=1.00, zorder=2, rasterized=True)

    plt.xlabel("MC step t")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / output_name, dpi=300)
    plt.close()


# ============================================================
# Plots
# No total quantities are plotted: no E_tot, no n_tot, no V_tot, no rho_tot.
# ============================================================

plot_gas_liquid(
    t=df["t"],
    gas=df["n_gas"],
    liq=df["n_liq"],
    ylabel="Number of particles",
    title="Particle numbers vs MC step",
    gas_label=r"$n_{\mathrm{gas}}$",
    liq_label=r"$n_{\mathrm{liq}}$",
    output_name="01_particle_numbers_gas_liquid.png"
)

plot_gas_liquid(
    t=df["t"],
    gas=df["V_gas"],
    liq=df["V_liq"],
    ylabel="Volume",
    title="Volumes vs MC step",
    gas_label=r"$V_{\mathrm{gas}}$",
    liq_label=r"$V_{\mathrm{liq}}$",
    output_name="02_volumes_gas_liquid.png"
)

plot_gas_liquid(
    t=df["t"],
    gas=df["rho_gas"],
    liq=df["rho_liq"],
    ylabel=r"Density $\rho$",
    title="Densities vs MC step",
    gas_label=r"$\rho_{\mathrm{gas}}$",
    liq_label=r"$\rho_{\mathrm{liq}}$",
    output_name="03_densities_gas_liquid.png"
)

plot_gas_liquid(
    t=df_pressure_blocks["t"],
    gas=df_pressure_blocks["P_gas"],
    liq=df_pressure_blocks["P_liq"],
    ylabel="Pressure",
    title=f"Pressures vs MC step, averaged over {PRESSURE_BLOCK_SIZE}-point blocks",
    gas_label=rf"$P_{{\mathrm{{gas}}}}$, eq avg = {P_gas_avg:.6f} +/- {P_gas_unc:.6f}",
    liq_label=rf"$P_{{\mathrm{{liq}}}}$, eq avg = {P_liq_avg:.6f} +/- {P_liq_unc:.6f}",
    output_name="04_pressures_gas_liquid.png",
)

plot_gas_liquid(
    t=mu_gas_blocks["t"],
    gas=mu_gas_blocks["mu"],
    liq=mu_liq_blocks["mu"],
    ylabel=r"Chemical potential $\mu$",
    title=r"Chemical potential from block-averaged Widom weights",
    gas_label=rf"$\mu_{{\mathrm{{gas}}}}$, eq avg = {mu_gas_avg:.6f} +/- {mu_gas_unc:.6f}",
    liq_label=rf"$\mu_{{\mathrm{{liq}}}}$, eq avg = {mu_liq_avg:.6f} +/- {mu_liq_unc:.6f}",
    output_name="05_mu_gas_liquid.png",
)

# ============================================================
# Minimal diagnostics
# ============================================================

print(f"Read measurements from: {MEASUREMENTS_FILE}")
print(f"Read chemical potentials from: {MU_FILE}")
print(f"Saved plots in: {OUTPUT_DIR.resolve()}")

print("\nFinal measured values:")
print(f"n_gas   = {df['n_gas'].iloc[-1]}")
print(f"n_liq   = {df['n_liq'].iloc[-1]}")
print(f"V_gas   = {df['V_gas'].iloc[-1]:.6f}")
print(f"V_liq   = {df['V_liq'].iloc[-1]:.6f}")
print(f"rho_gas = {df['rho_gas'].iloc[-1]:.6f}")
print(f"rho_liq = {df['rho_liq'].iloc[-1]:.6f}")
print(f"P_gas   = {df['P_gas'].iloc[-1]:.6f}")
print(f"P_liq   = {df['P_liq'].iloc[-1]:.6f}")


print("\nEquilibrated averages and uncertainties:")
print(f"Equilibration cutoff: t >= {equilibration_cutoff:.0f} ({EQUILIBRATION_FRACTION:.1%} of mc_steps)")
print(f"<P_gas>  = {P_gas_avg:.6f} +/- {P_gas_unc:.6f}")
print(f"<P_liq>  = {P_liq_avg:.6f} +/- {P_liq_unc:.6f}")
print(f"<mu_gas> = {mu_gas_avg:.6f} +/- {mu_gas_unc:.6f}")
print(f"<mu_liq> = {mu_liq_avg:.6f} +/- {mu_liq_unc:.6f}")

