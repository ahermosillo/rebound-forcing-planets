/**
 *
 * Planetary migration forces
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

double* tau_a;     /**< Migration timescale in years for all particles */
double* delta_a;     /**< Migration timescale in years for all particles */
double tmax;

void migration_forces(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
	 char* filename = "simulationarchive.bin";
    // Setup constants
    r->integrator    = REB_INTEGRATOR_MERCURIUS;
    //r->integrator    = REB_INTEGRATOR_IAS15;
    r->dt = 0.25*2.*M_PI;        // in year/(2*pi)
    r->additional_forces = migration_forces;     //Set function pointer to add dissipative forces.
    //r->heartbeat = heartbeat;  
    r->force_is_velocity_dependent = 1;
    tmax = 1.0*M_PI; //e7*2.*M_PI;    // in year/(2*pi)

	 reb_simulationarchive_automate_interval(r,filename,1e4*2.*M_PI);

    // Initial conditions
    // Parameters are those of Lee & Peale 2002, Figure 4. 
    struct reb_particle star = {0};
    star.m  = 1.0;            // This is a sub-solar mass star
    reb_add(r, star); 
    
	 //current planetary radii set to 1 hill radius. Change later!
	 struct reb_particle pN = reb_tools_orbit_to_particle(r->G, r->particles[0], 5.15E-05, 2.31036615E+01, 1.12142690E-02, 1.76797500E+00*M_PI/180.0, 2.65646853E+02*M_PI/180.0, 1.31794310E+02*M_PI/180.0, 2.67767281E+02*M_PI/180.0);
	 double rhill = pow(pN.m/(3.*star.m),1./3.)*30.1;
	 pN.r = rhill;
	 struct reb_particle pU = reb_tools_orbit_to_particle(r->G, r->particles[0], 4.36E-05, 1.62294119E+01, 4.44055860E-02, 7.72556000E-01*M_PI/180.0, 9.65413180E+01*M_PI/180.0, 7.39898210E+01*M_PI/180.0, 1.42955717E+02*M_PI/180.0);
	 rhill = pow(pU.m/(3.*star.m),1./3.)*19.2;
	 pU.r = rhill;
	 struct reb_particle pS = reb_tools_orbit_to_particle(r->G, r->particles[0], 2.86E-04, 8.78201719E+00, 5.57232190E-02, 2.48524000E+00*M_PI/180.0, 3.36013862E+02*M_PI/180.0, 1.13642811E+02*M_PI/180.0, 3.20346750E+02*M_PI/180.0);
	 rhill = pow(pS.m/(3.*star.m),1./3.)*9.4;
	 pS.r = rhill;
	 struct reb_particle pJ = reb_tools_orbit_to_particle(r->G, r->particles[0], 9.54E-04, 5.40426700E+00, 4.87750000E-02, 1.30500000E+00*M_PI/180.0, 2.75066000E+02*M_PI/180.0, 1.00492000E+02*M_PI/180.0, 1.88180000E+01*M_PI/180.0);
	 rhill = pow(pJ.m/(3.*star.m),1./3.)*5.2;
	 pJ.r = rhill;

    

	 reb_add(r, pJ); 
	 reb_add(r, pS); 
	 reb_add(r, pU); 
	 reb_add(r, pN); 

	 r->N_active = r->N;
    
	 reb_move_to_com(r);  

	// add test particles anytime after this	 

    tau_a = calloc(sizeof(double),r->N_active);
    delta_a = calloc(sizeof(double),r->N_active);

	 tau_a[1] = 2e6*2.*M_PI;
	 tau_a[2] = 2e6*2.*M_PI;
	 tau_a[3] = 2e6*2.*M_PI;
	 tau_a[4] = 2e6*2.*M_PI;

	 delta_a[1] = -0.2;
	 delta_a[2] = 0.8;
	 delta_a[3] = 3.0;
	 delta_a[4] = 7.0;

    tau_a[2] = 2.*M_PI*20000.0;    // Migration timescale of planet 2 is 20000 years.



    reb_integrate(r, tmax);
}

void migration_forces(struct reb_simulation* r){
    const int N = r->N_active;
	 const double time = r->t;
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = r->particles[0]; // calculate migration forces with respect to center of mass;
    for(int i=1;i<N;i++){
        if (tau_a[i]!=0){
            struct reb_particle* p = &(particles[i]);
				struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);
            const double dvx = p->vx - com.vx;
            const double dvy = p->vy - com.vy;
            const double dvz = p->vz - com.vz;
				//const double dx = p->x - com.x;
            //const double dy = p->y - com.y;
            //const double dz = p->z - com.z;
				//const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
				//const double r_1 = 1.0/sqrt(dx*dx + dy*dy + dz*dz);
				//const double helr = sqrt(dx*dx + dy*dy + dz*dz);



				//const double semi2 = 1.0/(2.0*r_1-v2/(1.+p->m));

				//printf("\n\n\n");
				//printf("%f %f\n", semi2,o.a);
				//printf("%f\n", helr);
				//printf("%f %f %f\n", dvx, dvy, dvz);
				//printf("%f %f %f\n", dx,dy,dz);

                
				p->ax += 0.5*delta_a[i]/tau_a[i]*exp(-time/tau_a[i])*dvx/o.a; 
				p->ay += 0.5*delta_a[i]/tau_a[i]*exp(-time/tau_a[i])*dvy/o.a; 
				p->az += 0.5*delta_a[i]/tau_a[i]*exp(-time/tau_a[i])*dvz/o.a; 
            }
        
       // com = reb_get_com_of_pair(com,particles[i]);
    }
}

//void heartbeat(struct reb_simulation* r){
//    if(reb_output_check(r, 20.*M_PI)){
//        reb_output_timing(r, tmax);
//    }
//    if(reb_output_check(r, 40.)){
//        reb_integrator_synchronize(r);
//        reb_output_orbits(r,"orbits.txt");
//        reb_move_to_com(r); 
//    }
//}
