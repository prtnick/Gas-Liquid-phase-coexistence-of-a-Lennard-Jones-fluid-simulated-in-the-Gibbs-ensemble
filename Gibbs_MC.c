#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 513

const int    mc_steps        = 20000000;
const int    output_steps    = 10000;
const int    measures_steps  = 1000;
const int    mu_measure_steps = 100000;
const double overall_density = 0.12;
const double delta           = 0.1;
const double delta_V         = 0.005;
const double r_cut           = 2.5; 
const double Temperature     = 0.70;
const double beta            = 1.0 / Temperature;

const double ratio_displacement= 100;
const double ratio_volumechange = 1;
const double ratio_transfer = 2000;

typedef struct {
    int n;
    double energy; 
    double virial; 
    double r[N][NDIM];
    double box[NDIM];
    const char* input_file;
    const char* label;
} Box;

Box gas = { .n = 0, .energy=0, .input_file = "gas.dat",    .label = "gas" };
Box liq = { .n = 0, .energy=0, .input_file = "liquid.dat", .label = "liq" };

typedef struct {
   double V_i; double V_f; 
   double boxn;
   double E_i; double E_f; 
   double Virial_i; double Virial_f; 
   double scale_factor; 
   double r[N][NDIM];
} Volume_Move; 

Volume_Move Volume_Move_g; 
Volume_Move Volume_Move_l; 

typedef struct {
    double energy; 
    double virial; 
} Energy_Virial;

Energy_Virial particle_Energy_Virial_at_position(const Box* b, const double pos[NDIM], int skip_index);
Energy_Virial particle_Energy_Virial(const Box* b, int i);

double E_tot=0 ;
double n_tot=0; 
double V_tot;
double e_cut;
double dV = 0.0; // value will be assigned in the main

char selected_target_box[4]= "xxx";
int attempts_to_liq = 0;
int attempts_to_gas = 0;
int accepted_to_liq = 0;
int accepted_to_gas = 0;

/*FUNCTIONS*/

double box_volume(const Box* b)
{
    double volume = 1.0;

    for(int d = 0; d < NDIM; ++d){
        volume *= b->box[d];
    }

    return volume;
}

void read_box(Box* b)
{
    FILE* fp = fopen(b->input_file, "r");
    if(fp == NULL){
        fprintf(stderr, "Error: could not open input file %s\n", b->input_file);
        exit(EXIT_FAILURE);
    }

    if(fscanf(fp, "%d\n", &b->n) != 1){
        fprintf(stderr, "Error: could not read number of particles from %s\n", b->input_file);
        exit(EXIT_FAILURE);
    }

    if(b->n < 0 || b->n > N){
        fprintf(stderr, "Error: file %s contains n = %d, but maximum allowed is N = %d and n must be > 0\n",
                b->input_file, b->n, N);
        exit(EXIT_FAILURE);
    }

    for(int d = 0; d < NDIM; ++d){
        double dmin, dmax;

        if(fscanf(fp, "%lf %lf\n", &dmin, &dmax) != 2){
            fprintf(stderr, "Error: could not read box size from %s\n", b->input_file);
            exit(EXIT_FAILURE);
        }

        b->box[d] = fabs(dmax - dmin);
    }

    for(int n = 0; n < b->n; ++n){
        for(int d = 0; d < NDIM; ++d){
            if(fscanf(fp, "%lf", &b->r[n][d]) != 1){
                fprintf(stderr, "Error: could not read particle coordinates from %s\n", b->input_file);
                exit(EXIT_FAILURE);
            }
        }

        double diameter;

        if(fscanf(fp, "%lf\n", &diameter) != 1){
            fprintf(stderr, "Error: could not read particle diameter from %s\n", b->input_file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(fp);
}

void read_data(void)
{
    read_box(&gas);
    read_box(&liq);

    n_tot = gas.n + liq.n;
}

void write_box(const Box* b, int step)
{
    char filename[128];
    snprintf(filename, sizeof(filename), "coords_%s_step%07d.dat", b->label, step);

    FILE* fp = fopen(filename, "w");
    if(fp == NULL){
        fprintf(stderr, "Error: could not write output file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "%d\n", b->n);

    for(int d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n", 0.0, b->box[d]);
    }

    for(int n = 0; n < b->n; ++n){
        for(int d = 0; d < NDIM; ++d){
            fprintf(fp, "%f\t", b->r[n][d]);
        }

        fprintf(fp, "%lf\n", 1.0);
    }

    fclose(fp);
}

void write_data(int step)
{
    write_box(&gas, step);
    write_box(&liq, step);
}

void rescale_box(Box* b, double scale_factor)
{
    for(int n = 0; n < b->n; ++n){
        for(int d = 0; d < NDIM; ++d){
            b->r[n][d] *= scale_factor;
        }
    }

    for(int d = 0; d < NDIM; ++d){
        b->box[d] *= scale_factor;
    }
}

void set_overall_density(void)
{
    double current_volume = box_volume(&gas) + box_volume(&liq);
    double target_volume  = ((double)n_tot) / overall_density;
    double scale_factor   = pow(target_volume / current_volume, 1.0 / NDIM);

    rescale_box(&gas, scale_factor);
    rescale_box(&liq, scale_factor);
}

double distance_squared_pbc(const Box* b, const double a[NDIM], const double c[NDIM])
{
    double r2 = 0.0;

    for(int d = 0; d < NDIM; ++d){
        double dr = a[d] - c[d];

        dr -= b->box[d] * round(dr / b->box[d]);
        r2 += dr * dr;
    }

    return r2;
}

Energy_Virial particle_Energy_Virial_at_position(const Box* b, const double pos[NDIM], int skip_index) 
{   // pos is just a point, we need this for particle_transfer,  
    // because there will be a new proposed position that is of course not a position of a particle in either of the 2 boxes
    // in that case it will be the skip_index value to be useless and it can be set to b->n   
    // of course I'll want the new particle to interact with all the other particles in the box
    Energy_Virial info; 
    info.energy=0.0; 
    info.virial=0.0; 

    double r_cut2 = r_cut * r_cut;                                                   
    for(int j = 0; j < b->n; ++j){                
        if(j == skip_index){
            continue;
        }

        double r2 = distance_squared_pbc(b, pos, b->r[j]);

        if(r2 < r_cut2){
            double inv_r2  = 1.0 / r2;
            double inv_r6  = inv_r2 * inv_r2 * inv_r2;
            double inv_r12 = inv_r6 * inv_r6;

            info.energy += 4.0 * (inv_r12 - inv_r6) - e_cut;
            info.virial += 24.0 * inv_r6 * (2.0 * inv_r6 - 1.0);
        }
    }

    return info;
}

Energy_Virial particle_Energy_Virial(const Box* b, int i)
{
     return particle_Energy_Virial_at_position(b, b->r[i], i);
}

Energy_Virial box_Energy_Virial(const Box* b)
{
    Energy_Virial info;
    info.energy = 0.0;
    info.virial = 0.0;

    for(int i = 0; i < b->n; ++i){
        Energy_Virial particle_info = particle_Energy_Virial(b, i);

        info.energy += particle_info.energy;
        info.virial += particle_info.virial;
    }

    info.energy *= 0.5;
    info.virial *= 0.5;

    return info;
}

double pressure_measurement(const Box* b)
{
    double volume = box_volume(b);
    double density = (double)b->n / volume;

    return density / beta + b->virial / (3.0 * volume);
}

double mu_measurement(const Box* b)
{
    if(b->n <= 0){
    return NAN;
    }

    int ntest = 200000;
    double sum_boltz = 0.0;
    double volume = box_volume(b); 

    for(int i = 0; i < ntest; ++i){
        double r_test[NDIM];

        for(int d = 0; d < NDIM; ++d){
            r_test[d] = b->box[d] * dsfmt_genrand();
        }

        Energy_Virial info = particle_Energy_Virial_at_position(b, r_test, -1);
        sum_boltz += (volume*exp(-beta * info.energy))/(b->n + 1);
    }

    
    

    //chemical potential is given by mu_id + mu_excess (mu_id= kTln(rho) by less than an unimportant constant given by the thermal wavelenght)

    return (1.0 / beta) *( - log(sum_boltz / ntest) );
}

void check_box(const Box* b)
{
    if(b->n == 0){
        fprintf(stderr, "Error: number of particles in %s box is zero.\n", b->label);
        exit(EXIT_FAILURE);
    }

    for(int d = 0; d < NDIM; ++d){
        assert(r_cut <= 0.5 * b->box[d]);
    }
}

/*MOVES*/
//1. Displace particle within the same box

int displacement(void)
{
// choose the box where to displace a random particle
    Box* b = NULL;
    double u=dsfmt_genrand(); 
    
    if(u < 0.5){
        b = &gas;
    } else if(u < 1.0){
            b = &liq;
        } else {
                fprintf(stderr, "Error: invalid random number in displacement function: u = %f\n", u);
                return -1;
            }
// you cannot displace anything if there are no particles (and sometimes it happens)

    if(b->n <= 0){
            printf(" the Box %s had 0 particles in the displacement function", b->label);
            return 0;
        }
            

// choose the particle to be displaced 
    int n=(int)((b->n)*dsfmt_genrand());
    // save the old position
    double old_pos[NDIM]; 
    for(int d=0; d<NDIM; d++){
        old_pos[d]=b->r[n][d]; 
    }
    Energy_Virial Old_particle_Energy_Virial=particle_Energy_Virial(b,n);

    //make a trial move 
    for(int d=0; d<NDIM; d++){
        double shift=(dsfmt_genrand() - 0.5) * 2.0 * delta;
        b->r[n][d]+=shift; 
        //respect PBC, keep particles in box
        if(b->r[n][d]<0){b->r[n][d]+= b->box[d];}
        if(b->r[n][d]>=b->box[d]){b->r[n][d]-= b->box[d];}
    }
    
    // test energy 
    Energy_Virial Trial_particle_Energy_Virial=particle_Energy_Virial(b,n); 
    //Monte Carlo move 
    double dE = Trial_particle_Energy_Virial.energy - Old_particle_Energy_Virial.energy;
    double dVirial = Trial_particle_Energy_Virial.virial - Old_particle_Energy_Virial.virial;
    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){
        //move accepted (also update energy of the box)
        b->energy += dE; 
        b->virial += dVirial;
        E_tot+=dE; //important because energy isn't conserved 
        return 1;
    }
    else{ //move rejected (revert)
        for(int d = 0; d < NDIM; ++d) b->r[n][d] = old_pos[d];
        return 0; 
    }
}    

//2. Volume change

void proposed_volume_move(Volume_Move*m, const Box*b){   // the const notation is just to signal that I am not modifying the old box (of course, it is the whole point of the function) I am just reading from it
    Box trial=*b; // I make a box that is identical to the original one, but now I modify it 
    for(int d=0; d<NDIM; d++){
         trial.box[d]=m->boxn;
    }
    for(int n=0; n<trial.n; n++){
        for(int d=0; d<NDIM; d++){
            trial.r[n][d]=b->r[n][d]*(m->scale_factor); 
        }
    }  // here I have created a modified trial-box on which to test the new energy 

    for(int n=0; n<trial.n; n++){
        for(int d=0; d<NDIM; d++){
            m->r[n][d]=trial.r[n][d]; 
        } //these modified positions will be instead be put in the Volume_Move struct, they are conceptually different, this data is saved, and will update the true boxes if the move is accepted
    }
    
    // intialise energy and virial 
   m->E_f=0.0; 
   m->Virial_f=0.0; 

   for(int n=0; n<trial.n; n++){
    Energy_Virial Trial_Energy_Virial= particle_Energy_Virial_at_position(&trial,trial.r[n],n);
   m->E_f += Trial_Energy_Virial.energy;
   m->Virial_f += Trial_Energy_Virial.virial; 
   }
   m->E_f *=0.5; m->Virial_f *=0.5; // don't count couple twice 
}

int change_volume(void)
{   Volume_Move_g.V_i = box_volume(&gas); 
    Volume_Move_l.V_i=V_tot-Volume_Move_g.V_i;
    Volume_Move_g.V_f=Volume_Move_g.V_i+ dV*((dsfmt_genrand()-0.5)*2);
    if(Volume_Move_g.V_f >=V_tot || Volume_Move_g.V_f<=0){
    return 0;
    }
    Volume_Move_l.V_f=V_tot-Volume_Move_g.V_f;  
    if(pow(Volume_Move_g.V_f, 1.0 / NDIM) < 2.0 * r_cut){
    return 0;
    }
    if(pow(Volume_Move_l.V_f, 1.0 / NDIM) < 2.0 * r_cut){
        return 0;
    }
    Volume_Move_g.boxn=pow((Volume_Move_g.V_f),1.0/3.0);
    Volume_Move_l.boxn=pow((Volume_Move_l.V_f),1.0/3.0);    
    Volume_Move_g.scale_factor=pow((Volume_Move_g.V_f/Volume_Move_g.V_i),1.0/3.0);
    Volume_Move_l.scale_factor=pow((Volume_Move_l.V_f/Volume_Move_l.V_i),1.0/3.0);
    Volume_Move_g.E_i=gas.energy; 
    Volume_Move_l.E_i=liq.energy; 
    proposed_volume_move(&Volume_Move_g,&gas);
    proposed_volume_move(&Volume_Move_l,&liq);
    double deltaE_g=Volume_Move_g.E_f-Volume_Move_g.E_i;
    double deltaE_l=Volume_Move_l.E_f-Volume_Move_l.E_i;
    double deltaV_g=Volume_Move_g.V_f-Volume_Move_g.V_i; double deltaV_l=-deltaV_g;
    double log_acc=-beta*deltaE_g-beta*deltaE_l+gas.n*log((Volume_Move_g.V_i+deltaV_g)/Volume_Move_g.V_i)+liq.n*log((Volume_Move_l.V_i+deltaV_l)/Volume_Move_l.V_i);
    double u=dsfmt_genrand();
    if (u <= 0.0) u = 1e-16; // in order not to have -inf (the consequence is I always refute moves that have probability 10^-16 which is not a problem)
    if(log_acc>=0 || log(u)<log_acc){
        gas.energy=Volume_Move_g.E_f;
        liq.energy=Volume_Move_l.E_f;
        gas.virial=Volume_Move_g.Virial_f; 
        liq.virial=Volume_Move_l.Virial_f; 
        E_tot=gas.energy+liq.energy; 
        for(int d=0; d<NDIM; d++){
            gas.box[d]=Volume_Move_g.boxn; 
            liq.box[d]=Volume_Move_l.boxn;
        }
        for(int n=0; n<gas.n; n++){
            for(int d=0; d<NDIM; d++){
                gas.r[n][d]=Volume_Move_g.r[n][d]; 
            }
        } 
        for(int n=0; n<liq.n; n++){
            for(int d=0; d<NDIM; d++){
                liq.r[n][d]=Volume_Move_l.r[n][d];            }
        } 
        
        return 1;
    } else return 0;
    
}

//3. Particle Transfer

int particle_transfer(void)
{  
    // choose the source box (and the target box))
    Box* source = NULL;
    Box* target = NULL;
    double u=dsfmt_genrand(); 
    if(u < 0.5){
        source = &gas;
        target = &liq;
        strcpy(selected_target_box, "liq");
    } else if(u < 1.0){
            source = &liq;
            target = &gas;
            strcpy(selected_target_box, "gas");
        } else {
                fprintf(stderr, "Error: invalid random number in displacement function: u = %f\n", u);
                strcpy(selected_target_box, "xxx");
                return -1;
            }

    // check on the validity of n_gas and n_liq 
    if(source->n <= 0){
        perror("\n ATTENTION : one box has 0 particles \n");
        return 0;
    }
    if(target->n >= N){
        perror("\n ATTENTION : risking overflow of the target box ");
        return 0; 
    }        

    // choose a random particle in the source box 
    int n=(int)((source->n)*dsfmt_genrand());
    // calculate its energy and Virial 

    Energy_Virial Source_Energy_Virial=particle_Energy_Virial(source,n); 
    double DeltaE_source = -Source_Energy_Virial.energy; 
    double DeltaVirial_source = -Source_Energy_Virial.virial; 

    // identify a new target position where to send the particle 
    double r_target[NDIM]; 
    for(int d=0;d<NDIM;d++){
        r_target[d]=(target->box[d])*dsfmt_genrand();
    }
    // identify what energy would be associated with the new particle
    // again E_final= energy old + energy of the target particle, so the delta is just the target particle energy
    
    Energy_Virial Target_Energy_Virial=particle_Energy_Virial_at_position(target,r_target,-1);
    double DeltaE_target = Target_Energy_Virial.energy;
    double DeltaVirial_target = Target_Energy_Virial.virial;

    double V_source=1;
    for(int d=0; d<NDIM; d++){V_source*=source->box[d];}
    double V_target=1;
    for(int d=0; d<NDIM; d++){V_target*=target->box[d];}
    double log_acc=-beta*DeltaE_source-beta*DeltaE_target-log(((target->n+1)*V_source)/(V_target*source->n));
    // use again u as a different random number for a different use 
    u=dsfmt_genrand(); 
    if (u <= 0.0) u = 1e-16; // in order not to have -inf (the consequence is I always refute moves that have probability 10^-16 which is not a problem)
    if(log_acc>=0 || log(u)<log_acc){//acceptance 

        //modify source box
        source->n -= 1;
        for(int idx = n; idx<source->n; idx ++){
            for(int d=0; d<NDIM; d++){
                source->r[idx][d]=source->r[idx+1][d];
            }
        }
        source->energy+=DeltaE_source;
        source->virial+=DeltaVirial_source; 

        //modify target box
        if(target->n >= N){
        fprintf(stderr, "Error: target box overflow\n");
        exit(EXIT_FAILURE);
}
        for(int d=0; d<NDIM; d++){
            target->r[target->n][d]=r_target[d]; // the particle is added on the n+1 th index for simplicity
        }
        target->n += 1;   if(target->n >=N){perror("\n target box overflow \n");}
        target->energy+=DeltaE_target;
        E_tot= source->energy + target->energy;
        target->virial+=DeltaVirial_target; 
        return 1;
    }
    else return 0; 
}


int main(int argc, char* argv[]){
    assert(delta > 0.0);
    e_cut = 4.0 * (pow(1.0 / r_cut, 12.0) - pow(1.0 / r_cut, 6.0)); //did you not make this a simulation variable on purpose?

    read_data();
    set_overall_density();  // Rescale to correct density so particles interact from start
    printf("nr of particles %lf\n", n_tot);
    check_box(&gas);
    check_box(&liq);
    printf("\n ngas=%d , nliq=%d\n", gas.n ,liq.n );
    

    size_t seed = time(NULL);
    dsfmt_seed(seed);
    printf("box volume: gas = %lf, liq= %lf \n", box_volume(&gas), box_volume(&liq));
    dV = delta_V * (box_volume(&gas) + box_volume(&liq)); 

    // intialise energy and virial of the 2 boxes 
    Energy_Virial init; 
    init = box_Energy_Virial(&gas); 
    gas.energy = init.energy; 
    gas.virial = init.virial; 
    init = box_Energy_Virial(&liq); 
    liq.energy = init.energy; 
    liq.virial = init.virial; 


    E_tot=gas.energy+liq.energy;
    V_tot=box_volume(&gas)+box_volume(&liq);

    printf("Starting gas volume:    %f\n", box_volume(&gas));
    printf("Starting liquid volume: %f\n", box_volume(&liq));
    printf("Starting total volume:  %f\n", V_tot);
    printf("Starting total energy:  %lf\n", E_tot);
    printf("Starting seed:          %lu\n", seed);

    FILE* fp1 = fopen("measurements.dat", "w");
    if(fp1 == NULL){
        fprintf(stderr, "Error: could not write measurements.dat\n");
        exit(EXIT_FAILURE);
    }

    FILE* fp_mu = fopen("mu_measurements.dat", "w");
    if(fp_mu == NULL){
        fprintf(stderr, "Error: could not write mu_measurements.dat\n");
        exit(EXIT_FAILURE);
    }

    // set how likely each move is to be selected 
 
    double total_ratio = ratio_displacement + ratio_transfer + ratio_volumechange;
    double displacment_divide = ratio_displacement; 
    double changevolume_divide = ratio_displacement + ratio_volumechange; 
    double transfer_divide = total_ratio;

    fprintf(fp1, "t \t E_tot \t\t n_gas \t n_liq \t V_gas \t\t V_liq \t P_gas \t P_liq \n"); 
    fprintf(fp_mu, "t \t mu_gas \t mu_liq \n");
 
    int accepted_displacement = 0;  int tot_disp = 0;
    int accepted_volume = 0;        int tot_vol_ch = 0; 
    int accepted_transfer = 0;      int tot_transf = 0;
    int vol_print_count=0; 
    const int vol_print_count_threshold = 5;

    for(int step=0; step < mc_steps; step++){
        
        //select a random step
        double  mc_move_extraction = dsfmt_genrand()*total_ratio;
        if(mc_move_extraction<0){perror("invalid value of mc_move_extraction \n");}

        if(mc_move_extraction < displacment_divide){
            tot_disp++;
            accepted_displacement+= displacement(); 
        }   else if(mc_move_extraction < changevolume_divide){
            tot_vol_ch++;
            accepted_volume +=  change_volume(); 
            }   else if(mc_move_extraction <= transfer_divide){
                tot_transf++;

                int transfer_output= particle_transfer();
                    if(strcmp(selected_target_box, "gas") == 0){
                    attempts_to_gas ++;  accepted_to_gas += transfer_output;
                    } else  if(strcmp(selected_target_box, "liq") == 0){attempts_to_liq ++; accepted_to_liq+= transfer_output; }
                    accepted_transfer += transfer_output; 

                } else perror("there is some mistake in the implementation of the mc step, dsfmt_genrand can't generate numbers greater than 1 \n");
        

        if( step % measures_steps == 0) fprintf(fp1, "%d \t %lf \t %d \t %d \t %lf \t %lf \t %lf \t %lf \n", step, E_tot, gas.n, liq.n, gas.box[0]*gas.box[1]*gas.box[2], liq.box[0]*liq.box[1]*liq.box[2], pressure_measurement(&gas), pressure_measurement(&liq)); 

       if(step % mu_measure_steps == 0){
            fprintf(fp_mu, "%d \t %lf \t %lf \n", step, mu_measurement(&gas), mu_measurement(&liq));
        }  

        if(step % output_steps == 0){
            vol_print_count++; 
            if(step > 0){
                printf("Step %d. Move acceptance: displacement acceptance = %lf. \t transfer acceptance = %lf. \t", step, (double)accepted_displacement/tot_disp, (double)accepted_transfer/tot_transf);
                accepted_displacement=accepted_transfer=tot_disp=tot_transf=0; 
            }
            if(vol_print_count == vol_print_count_threshold){
                printf(" change volume acceptance = %lf. \n",(double)accepted_volume/tot_vol_ch );
                vol_print_count=0; accepted_volume=tot_vol_ch=0; 
            } 
            printf("n_gas = %d \t n_liq = %d \n \n", gas.n, liq.n);
           // write_data(step);
        }
    }

    printf("\n attempts to gas = %d, attempts to liq = %d,  effective transfer to gas  = %d, effective transfer to liq = %d \n", attempts_to_gas, attempts_to_liq, accepted_to_gas, accepted_to_liq);
    printf(" ratio_gas = %lf , ratio_liq = %lf \n", (double)accepted_to_gas/attempts_to_gas , (double)accepted_to_liq/attempts_to_liq ); 
    write_data(mc_steps);
    fclose(fp1);
    fclose(fp_mu);
    return 0;
}
