#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 513

const int    mc_steps        = 100;
const int    output_steps    = 10;
const double overall_density = 0.2;
const double delta           = 0.1; //displacements
const double delta_V         = 0.05; //volume changes
const double r_cut           = 0.5; //2.5;
const double Temperature     = 2.0;
const double beta            = 1.0 / Temperature;

/*DEFINING STRUCTURES (BOX (liquid + gas) + 1 FOR EVERY MOVE(liquid + gas))*/
typedef struct {
    int n;
    //int n_initial;
    double energy; 
    double r[2*N][NDIM]; //extra space for added particles
    double box[NDIM];
    const char* input_file;
    const char* label;
} Box;

Box gas = { .n=0, .energy=0, .input_file = "gas.dat",    .label = "gas" };
Box liq = { .n=0, .energy=0, .input_file = "liquid.dat", .label = "liq" };

typedef struct {
   double V_i; double V_f; 
   double boxn;
   double E_i; double E_f; 
   double scale_factor; 
   double r[2*N][NDIM];
} Volume_Move; 
Volume_Move Volume_Move_g; 
Volume_Move Volume_Move_l;


double E_tot = 0 ;
double n_tot = 0; 
double V_tot;
double e_cut;
double dV = 0.0;
double diameter;

/*FUNCTIONS*/
void read_box(Box* b)
{
    FILE* fp = fopen(b->input_file, "r");
    if(fp == NULL){
        fprintf(stderr, "Error: could not open input file %s\n", b->input_file);
        exit(EXIT_FAILURE);
    }

    if(fscanf(fp, "%d\n", &b->n) != 1){ //!? n seems to stay zero //n_initial
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

        if(fscanf(fp, "%lf\n", &diameter) != 1){ //so it reads already and then it immediately checks this constraint?
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

    //gas.n = gas.n_initial;
    //liq.n = liq.n_initial;
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


double box_volume(const Box* b)
{
    double volume = 1.0;

    for(int d = 0; d < NDIM; ++d){
        volume *= b->box[d];
    }

    return volume;
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

double particle_energy_at_position(const Box* b, const double pos[NDIM], int skip_index) 
{   // pos is just a point, we need this for particle_transfer,  
    // because there will be a new proposed position that is of course not a position of a particle in either of the 2 boxes
    // in that case it will be the skip_index value to be useless and it can be set to b->n   
    // of course I'll want the new particle to interact with all the other particles in the box
    double en = 0.0;                             
    double r_cut2 = r_cut * r_cut;                                                   
    for(int j = 0; j < b->n; ++j){                
        if(j == skip_index){
            continue;
        }

        double r2 = distance_squared_pbc(b, pos, b->r[j]);

        if(r2 < r_cut2 && r2>0){
            double inv_r2  = 1.0 / r2;
            double inv_r6  = inv_r2 * inv_r2 * inv_r2;
            double inv_r12 = inv_r6 * inv_r6;

            en += 4.0 * (inv_r12 - inv_r6) - e_cut;
        }
        if(r2==0){
            printf("error: division by 0");
        }
    }

    return en;
}

double particle_energy(const Box* b, int i)
{
    assert(i >= 0 && i < b->n);

     return particle_energy_at_position(b, b->r[i], i);
}

double box_energy(const Box* b)
{
    double en = 0.0;

    for(int i = 0; i < b->n; ++i){
        en += particle_energy(b, i);
    }

    return 0.5 * en;
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

void update_positions(Box *target, int t, const double r_added[NDIM], Box *source, int s){ 
    //note: this function assumes particle number is already updated
    double updated_source[2*N][NDIM];
    double updated_target[2*N][NDIM];
    int source_n = source->n;
    int target_n = target->n;

    //delete s from source box
    for(int n = 0; n < s; n++){
        for(int d = 0; d < NDIM; d++){
            updated_source[n][d] = source->r[n][d]; //copy until line l-1
        }
    }
    for(int l=s; l < source_n; l++){ //change to l < ((target->n)-1) for non-updated n
        for(int d = 0; d < NDIM; d++){
            updated_source[l][d] = source->r[l+1][d]; //proceed to copy the positions in lines l+1 until n-1. Note we skipped l
        }
    }

    //add t to target box
    for(int n=0; n < t; n++){
        for(int d = 0; d < NDIM; d++){
            updated_target[n][d]=target->r[n][d]; //copy until line t-1
        }
    }
    for(int l=t; l < target_n; l++){ //change to l < ((target->n)+1) for non-updated n
        for(int d = 0; d < NDIM; d++){
            if(l == t) updated_target[l][d] = r_added[d];
            else updated_target[l][d] = target->r[l-1][d];
        }
    }

//put updated values in Box
    for(int n=0; n < source_n; n++){ 
        for(int d=0; d<NDIM; d++){
            source->r[n][d]=updated_source[n][d]; 
        }
    } 
        for(int n=0; n < target_n; n++){
            for(int d=0; d<NDIM; d++){
                target->r[n][d]=updated_target[n][d];            
            }
        } 
}




/*MOVES*/
//1. Displacement
int displacement(void){
    Box* b = NULL;
    double u = dsfmt_genrand();

    //choose box randomly
    if(u < 0.5){
        b = &gas;
    } else {
        b = &liq;
    }

    if(b->n <= 0){
        return 0;
    }

    int n = (int)(b->n * dsfmt_genrand()); //select random particle

    double old_pos[NDIM];
    for(int d = 0; d < NDIM; d++){
        old_pos[d] = b->r[n][d];
    }

    double old_particle_energy = particle_energy(b, n); //negative

    for(int d = 0; d < NDIM; d++){
        double shift = (dsfmt_genrand() - 0.5) * 2.0 * delta;
        b->r[n][d] += shift;
        //pbc keep particles in box
        if(b->r[n][d] < 0.0){
            b->r[n][d] += b->box[d];
        }
        if(b->r[n][d] >= b->box[d]){
            b->r[n][d] -= b->box[d];
        }
    }

    double new_particle_energy = particle_energy(b, n); //prob: super small e-310
    double dE = new_particle_energy - old_particle_energy; //small - negative

    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){
        b->energy += dE;
        E_tot += dE;
        return 1;
    } else {
        for(int d = 0; d < NDIM; d++){
            b->r[n][d] = old_pos[d];
        }
        return 0;
    }
}

//2. Volume change
void proposed_volume_move(Volume_Move*m, const Box*b){   // the const notation is just to signal that I am not modifying the old box (of course, it is the whole point of the function) I am just reading from it
    Box trial=*b; // I make a box that is identical to the original one, but now I modify it

    for(int d=0; d<NDIM; d++){
         trial.box[d] = m->boxn; //m->boxn is the proposed new box length
    }
    for(int n=0; n<trial.n; n++){
        for(int d=0; d<NDIM; d++){
            trial.r[n][d] = b->r[n][d]*(m->scale_factor); 
        }
    }  // here I have created a modified trial-box on which to test the new energy (by means of the new particle positions)

    for(int n=0; n<trial.n; n++){
        for(int d=0; d<NDIM; d++){
            m->r[n][d] = trial.r[n][d]; 
        } //these modified positions will be instead be put in the Volume_Move struct, they are conceptually different, this data is saved, and will update the true boxes if the move is accepted
    }
   m->E_f=0; 
   for(int n=0; n<trial.n; n++){
   m->E_f+= particle_energy_at_position(&trial,trial.r[n],n); 
   }
   m->E_f *=0.5; // don't count couple twice 
   trial.energy=m->E_f; //just for the sake of clarity 
}

int change_volume(void){
    //set initial information in g_i and l_i
    Volume_Move_g.V_i = box_volume(&gas); //old volume of gas
    Volume_Move_l.V_i = V_tot - Volume_Move_g.V_i; 
    Volume_Move_g.V_f = Volume_Move_g.V_i+ dV*((dsfmt_genrand()-0.5)*2); //new volume of gas
    if(Volume_Move_g.V_f >= V_tot || Volume_Move_g.V_f<=0){ //smart
    return 0;
    }
    Volume_Move_l.V_f = V_tot - Volume_Move_g.V_f; //(if new volume is not accepted this doesnt happen)
    if(pow(Volume_Move_g.V_f, 1.0 / NDIM) < 2.0 * r_cut){ //why do we include this? side < 2r_cut so the box shouldn´t be bigger than the cutoff radius?
    return 0;
    }
    if(pow(Volume_Move_l.V_f, 1.0 / NDIM) < 2.0 * r_cut){
        return 0;
    }
    Volume_Move_g.boxn = pow((Volume_Move_g.V_f),1.0/3.0);
    Volume_Move_l.boxn = pow((Volume_Move_l.V_f),1.0/3.0);    
    Volume_Move_g.scale_factor = pow((Volume_Move_g.V_f/Volume_Move_g.V_i),1.0/3.0); //side_new/side_old
    Volume_Move_l.scale_factor = pow((Volume_Move_l.V_f/Volume_Move_l.V_i),1.0/3.0);
    Volume_Move_g.E_i = gas.energy; //for myself: remeber gas.energy = (box->gas).energy 
    Volume_Move_l.E_i = liq.energy; 
    //generate potentail new volumes
    proposed_volume_move(&Volume_Move_g,&gas); 
    proposed_volume_move(&Volume_Move_l,&liq); 
    //accept or reject based on energy calculation (remember E_f calculated in proposed_particle_volume())
    double deltaE_g = Volume_Move_g.E_f-Volume_Move_g.E_i;
    double deltaE_l = Volume_Move_l.E_f-Volume_Move_l.E_i;
    double deltaV_g = Volume_Move_g.V_f-Volume_Move_g.V_i; double deltaV_l=-deltaV_g;
    double log_acc = -beta * deltaE_g - beta*deltaE_l + gas.n*log((Volume_Move_g.V_i + deltaV_g)/Volume_Move_g.V_i) + liq.n*log((Volume_Move_l.V_i + deltaV_l)/Volume_Move_l.V_i);
    double u = dsfmt_genrand();
    if (u <= 0.0) u = 1e-16; // in order not to have -inf (the consequence is we always refute moves that have probability 10^-16 which is not a problem)
    if(log_acc>=0 || log(u)<log_acc){ //accept
        gas.energy=Volume_Move_g.E_f; //update energies
        liq.energy=Volume_Move_l.E_f;
        for(int d=0; d<NDIM; d++){
            gas.box[d]=Volume_Move_g.boxn; //update box dimensions (NDIM is same for liq and gas)
            liq.box[d]=Volume_Move_l.boxn;
        }
        for(int n=0; n<gas.n; n++){ //update particle positions. note: nr of particles may be different so we need 2 seperate loops
            for(int d=0; d<NDIM; d++){
                gas.r[n][d]=Volume_Move_g.r[n][d]; 
            }
        } 
        for(int n=0; n<liq.n; n++){
            for(int d=0; d<NDIM; d++){
                liq.r[n][d]=Volume_Move_l.r[n][d];            
            }
        } 
        
        return 1;
    } else return 0;
    
}

//3. transfer partcile from one box to the other
int particle_transfer(void){
    Box*target = NULL;
    Box*source = NULL;
    double u = dsfmt_genrand();

    //choose source and target randomly (50/50)
    if(u < 0.5){
        target = &gas;
        source = &liq;
    } else {
        target = &liq;
        source = &gas;
    }

    //energy of removing a random particle:
    int s = (int)(dsfmt_genrand() * source->n); //particle removed from source (if accepted)
    double E_removed_particle = particle_energy_at_position(source, source->r[s], s); 

    //energy of adding a particle at random position:
    int t = (int)(dsfmt_genrand() * target->n); //particle added to target (if accepted)
    double r_added[NDIM];
    for(int d = 0; d < NDIM; d++){
        double side = target->box[d];
        double random_coordinate = dsfmt_genrand() * side; //i think theres no need for pbc now
        r_added[d]=random_coordinate;
    } 
    //check for overlap (commented out for now)
    // for(int n=0; n < (target->n)+1; n++){
    // if(n!=t){ //skip self
    //     double dr[NDIM]=0;
    //     for(int d=0; d<NDIM; d++){
    //         dr[d]=(m->r_target[t][d]-m->r_target[n][d]);
    //         dr[d]-= (target->box)[d]*round(dr[d]/(target->box)[d]); //compute shortest distance
    //     }
    //     if(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]<diameter*diameter){ // delete the move
    //         for(int d=0;d<NDIM;d++){
    //             m->r_target[t][d] = (target->r)[t][d]; //delete the changes
    //         } 
    //         return 0;
    //     }    
    //     else continue; 
    // }
    double E_added_particle = particle_energy_at_position(target, r_added, t); //prob: this becomes huge e7 probably due to overlap since i inserted a fcc file hahha

    //calculate energy difference
    double deltaE = E_added_particle - E_removed_particle; 
    double V_s = pow(source->box[0],3); //volume source
    double V_t = pow(target->box[0],3); //volume target 
    double configE = ( ( (target->n)+1 ) * V_s)  / ((source->n)* V_t);

    if(deltaE < 0.0 || dsfmt_genrand() < exp(-beta * (deltaE) - log(configE))){
    //accept: update energies
        target->energy += E_added_particle; 
        source->energy -= E_removed_particle;
    //update positions
        target->n += 1;
        source->n -= 1;
        update_positions(target, t, r_added, source, s); //this function assumes updated particle nr n
        return 1;
    } else {
        return 0;
    }
}


int main(int argc, char* argv[]){
    assert(delta > 0.0);
    printf("hoi!!!!!!!!!!!!!!!1\n");
    e_cut = 4.0 * (pow(1.0 / r_cut, 12.0) - pow(1.0 / r_cut, 6.0)); //did you not make this a simulation variable on purpose?

    read_data();
    printf("nr of particles %lf\n", n_tot);
    check_box(&gas);
    check_box(&liq);

    size_t seed = time(NULL);
    dsfmt_seed(seed);

    dV = delta_V * (box_volume(&gas) + box_volume(&liq)); 
    gas.energy = box_energy(&gas);
    liq.energy = box_energy(&liq);
    E_tot=gas.energy+liq.energy;
    V_tot=box_volume(&gas)+box_volume(&liq);

    printf("Starting gas volume:    %f\n", box_volume(&gas));
    printf("Starting liquid volume: %f\n", box_volume(&liq));
    printf("Starting total volume:  %f\n", V_tot);
    printf("Starting total energy:  %lf\n", E_tot);
    printf("Starting seed:          %lu\n", seed);

    FILE* fp = fopen("measurements.dat", "w");
    if(fp == NULL){
        fprintf(stderr, "Error: could not write measurements.dat\n");
        exit(EXIT_FAILURE);
    }
    

    //ratios
    double ratio_displacement= 10;
    double ratio_volumechange = 100;
    double ratio_transfer = 1; 
    double total_ratio = ratio_displacement + ratio_transfer + ratio_volumechange;

    int accepted = 0; //maybe keep track seperately of each move
    for(int step=0; step < mc_steps; step++){
        double m = dsfmt_genrand() * total_ratio; //which move

        // if(m <= ratio_transfer) {
        //     accepted += particle_transfer();
        // }
        // if(m > ratio_displacement && m < ratio_volumechange){
        //     accepted += displacement();
        // }
        // if(m >= ratio_volumechange){
        //     accepted += change_volume(); //Is it correct this one randomly chooses a box to shrink inside the function?
        // }

        accepted += particle_transfer();

        if(step % output_steps == 0){
            printf("Step %d. Move acceptance: %lf.\n", step, (double)accepted);
            printf("n_gas = %d \t n_liq = %d \n", gas.n, liq.n);
            accepted = 0;
            //write_data(step);
        }
    }
    write_data(mc_steps);

    fclose(fp);
    return 0;
}
