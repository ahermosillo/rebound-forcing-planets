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

double* tau_a;      /**< Migration timescale in years for all particles */
double* tau_e;      /**< Migration timescale in years for all particles */
double* delta_a;    /**< Migration timescale in years for all particles */
double tmax;


void migration_forces(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);
void reb_check_remove(struct reb_simulation* r);

int main(int argc, char* argv[]){

   //initializing the simulation variable so it's always accessible outside if statements
   struct reb_simulation* r = reb_create_simulation();
   reb_free_simulation(r);
   r = NULL;
   
   
   //reading in the basic integration parameters
   FILE *myfile;
   double myvariable;
   char* intfile = "integration-parameters.txt";
   myfile=fopen(intfile, "r");
   fscanf(myfile,"%lf",&myvariable);
   double dt = myvariable;
   fscanf(myfile,"%lf",&myvariable);
   double tmax = 2.0*M_PI*myvariable;
   fscanf(myfile,"%lf",&myvariable);
   double tout = myvariable;
   fscanf(myfile,"%lf",&myvariable);
   int trial_sim = myvariable;
   fclose(myfile);


	int ndtout = (int)floor(tout/dt);
	int highresndtout = (int)floor(1.e3/dt);

	//convert the timestep to the correct units
	dt = 2.0*M_PI*dt;

   char* migfile = "planet-migration-parameters.txt";
   tau_a = calloc(sizeof(double),5);
   delta_a = calloc(sizeof(double),5);
   tau_e = calloc(sizeof(double),5);
   myfile=fopen(migfile, "r");
   for(int i = 0; i <= 4; i++){
      fscanf(myfile,"%lf",&myvariable);
      tau_a[i] = 2.0*M_PI*myvariable;
      fscanf(myfile,"%lf",&myvariable);
      delta_a[i] = myvariable;
      fscanf(myfile,"%lf",&myvariable);
      tau_e[i] = 2.0*M_PI*myvariable;
   }
   fclose(myfile);

   int continued_sim = 0;

   char* filename = "simulationarchive.bin";
   // Trying to open a SimulationArchive file
   struct reb_simulationarchive* sa = reb_open_simulationarchive(filename);
   if (sa==NULL){
      printf("Cannot open simulation archive file. Not a continuation of previous simulation\n");
   }
   else{
      // Get a simulation from the file (if possible, otherwise NULL is returned)
      printf("Reading in from a simulation archive to do a continuation of previous simulation\n");
      r = reb_create_simulation_from_simulationarchive(sa,-1);
		if(r!=NULL){
         continued_sim = 1;
         printf("Simulation to be continued from time %f \n", r->t);
      }
   }
   reb_close_simulationarchive(sa);
   
   char* binplfile = "planet-ic.bin";
   if (r==NULL){
      //////**************************************************************** ////////////    
      //////This is all code for either reading in a binary ic file and adding test particles
      //////or reading in planet positions from a text file and skipping test particles
      //////to just see what the planets do
      //////**************************************************************** ////////////    

      r = reb_create_simulation_from_binary(binplfile);      
   }
   // Check if we were successful
   if (r==NULL){
      printf("No binary initial conditions found. Testing planet initial conditions from text file.\n");
      printf("There will be no test particles\n");
      //trial_sim = 1;
      r = reb_create_simulation();
      //r->integrator = REB_INTEGRATOR_MERCURIUS;
      //r->dt = dt;       
      //r->t = 0.;      
      //r->additional_forces = migration_forces;      
      //r->force_is_velocity_dependent = 1;
		//r->heartbeat = heartbeat;

      char* txticfile = "planet-ic.txt";
     
      double* a0 = calloc(sizeof(double),4);
      double* e0 = calloc(sizeof(double),4);
      double* inc0 = calloc(sizeof(double),4);
      double* node0 = calloc(sizeof(double),4);
      double* w0 = calloc(sizeof(double),4);
      double* f0 = calloc(sizeof(double),4);
      double* rpl = calloc(sizeof(double),4);
      double* masspl = calloc(sizeof(double),4);

      myfile=fopen(txticfile, "r");
      for(int i = 0; i < 4; i++){
         fscanf(myfile,"%lf",&myvariable);
         masspl[i] = myvariable;
         fscanf(myfile,"%lf",&myvariable);
         rpl[i] = myvariable;
         fscanf(myfile,"%lf",&myvariable);
         a0[i] = myvariable;
         fscanf(myfile,"%lf",&myvariable);
         e0[i] = myvariable;
         fscanf(myfile,"%lf",&myvariable);
         inc0[i] = (M_PI/180.0)*myvariable;
         fscanf(myfile,"%lf",&myvariable);
         node0[i] = (M_PI/180.0)*myvariable;
         fscanf(myfile,"%lf",&myvariable);
         w0[i] = (M_PI/180.0)*myvariable;
         fscanf(myfile,"%lf",&myvariable);
         double M0 = (M_PI/180.0)*myvariable;
			f0[i] = reb_tools_M_to_f(e0[i], M0);
      }
      fclose(myfile);

      struct reb_particle star = {0};
      star.m  = 1.0;        
      reb_add(r, star); 
   
      struct reb_particle pJ = reb_tools_orbit_to_particle(r->G, r->particles[0],masspl[0],a0[0],e0[0],inc0[0],node0[0],w0[0],f0[0]);
      pJ.r = rpl[0];
      pJ.hash = 1;
      struct reb_particle pS = reb_tools_orbit_to_particle(r->G, r->particles[0],masspl[1],a0[1],e0[1],inc0[1],node0[1],w0[1],f0[1]);
      pS.r = rpl[1];
      pS.hash = 2;
      struct reb_particle pU = reb_tools_orbit_to_particle(r->G, r->particles[0],masspl[2],a0[2],e0[2],inc0[2],node0[2],w0[2],f0[2]);
      pU.r = rpl[2];
      pU.hash = 3;
      struct reb_particle pN = reb_tools_orbit_to_particle(r->G, r->particles[0],masspl[3],a0[3],e0[3],inc0[3],node0[3],w0[3],f0[3]);
      pN.r = rpl[3];
      pN.hash = 4;
      reb_add(r, pJ); 
      reb_add(r, pS); 
      reb_add(r, pU); 
      reb_add(r, pN);
     
      r->N_active = r->N;
      reb_move_to_com(r);  
      //write these initial conditions to a binary file so they can be used again reproducibly
      reb_simulationarchive_snapshot(r,binplfile);
   }
   else if (trial_sim > 0){
      printf("Loaded simulation state from binary file. Continuing simulation \n");
      //r->integrator = REB_INTEGRATOR_MERCURIUS;
      //r->dt = dt;       
      //r->additional_forces = migration_forces;      
      //r->force_is_velocity_dependent = 1;
		//r->heartbeat = heartbeat;
      //reb_move_to_com(r);
   }
   else if (continued_sim == 0){
      printf("Loaded planet initial conditions from binary file. Adding test particles \n");
      //r->integrator = REB_INTEGRATOR_MERCURIUS;
      //r->dt = dt;       
      //r->additional_forces = migration_forces;      
      //r->force_is_velocity_dependent = 1;
      //r->N_active = r->N;
      //reb_move_to_com(r);
		//r->heartbeat = heartbeat;
 
      char* tpfile = "tp-parameters.txt";
      myfile=fopen(tpfile, "r");
      int myintvariable;
      fscanf(myfile,"%d",&myintvariable);
      int N_tp = myintvariable;
      printf("ntp: %d",N_tp);
      fscanf(myfile,"%lf",&myvariable);
      double tpamin = myvariable;
      fscanf(myfile,"%lf",&myvariable);
      double tpamax = myvariable;
      fscanf(myfile,"%lf",&myvariable);
      double tpde = myvariable;
      fscanf(myfile,"%lf",&myvariable);
      double tpdi = myvariable;
      fscanf(myfile,"%d",&myintvariable);
      int new_seed = myintvariable;
      fclose(myfile);      

      r->rand_seed = new_seed;      
      int k=101;
      printf("\nAdding %i test particles from  a=%.16f - %.16f with de=%.16f, di=%.16f.\n\n",N_tp,tpamin,tpamax,tpde,tpdi);
      for( int j = r->N; j < N_tp + r->N_active; j++){
         double a = reb_random_uniform(r,tpamin,tpamax);
         double e = reb_random_rayleigh(r, tpde);
         double inc = reb_random_rayleigh(r, tpdi);
         double Omega = reb_random_uniform(r, 0.,2.*M_PI);
         double omega = reb_random_uniform(r, 0.,2.*M_PI);
         double f = reb_random_uniform(r, 0.,2.*M_PI);
         struct reb_particle pt = reb_tools_orbit_to_particle(r->G, r->particles[0], 0., a, e, inc, Omega, omega, f);
         pt.hash = k;
         k+=1;
         reb_add(r, pt);
      }

      //////**************************************************************** ////////////     
      ////// End of just setting up initial conditions!
      ////////**************************************************************** ////////////  
   }
       
   printf("nactive = %i\n",r->N_active);
   printf("ntotal = %i\n",r->N);
   r->integrator = REB_INTEGRATOR_MERCURIUS;
   r->testparticle_type = 0;
   r->dt = dt;        
	r->heartbeat = heartbeat;
   r->additional_forces = migration_forces;      
   r->force_is_velocity_dependent = 1;
   printf("integrating to %.16f saving every %.16f time units\n",tmax,tout);
	reb_simulationarchive_automate_step(r,filename,ndtout);
   //reb_simulationarchive_automate_interval(r,filename,tout);
   reb_move_to_com(r);
   reb_integrate(r, tmax);
	reb_check_remove(r);
	reb_simulationarchive_snapshot(r,filename);

   reb_free_simulation(r);
   r = NULL;
   tout = 1e3*2.*M_PI;
   printf("integrating 10 Myr from the final state at %.16f saving every %.16f time units\n",tmax,tout);
   struct reb_simulationarchive* sae = reb_open_simulationarchive(filename);
   r = reb_create_simulation_from_simulationarchive(sae,-1);
   reb_close_simulationarchive(sae);
   reb_move_to_com(r);
   r->integrator = REB_INTEGRATOR_MERCURIUS;
	r->testparticle_type = 0;
   r->dt = dt;      
   r->additional_forces = migration_forces;      
   r->force_is_velocity_dependent = 1;
	r->heartbeat = heartbeat;

	printf("nactive = %i\n",r->N_active);
   printf("ntotal = %i\n",r->N);

   char* endfilename = "high-res-end-archive.bin";
	reb_simulationarchive_automate_step(r,endfilename,highresndtout);
   //reb_simulationarchive_automate_interval(r,endfilename,tout);
   double newtmax = r->t + 1.e7*2.*M_PI;	
   reb_integrate(r, newtmax);

}

void migration_forces(struct reb_simulation* r){
   const int N = r->N_active;
   const double time = r->t;
   struct reb_particle* const particles = r->particles;
   struct reb_particle com = reb_get_com(r); // calculate migration forces with respect to the center of mass;
   for(int i=1;i<N;i++){
      if (tau_a[i]!=0){
         struct reb_particle* p = &(particles[i]);
         struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);
         const double dvx = p->vx - com.vx;
         const double dvy = p->vy - com.vy;
         const double dvz = p->vz - com.vz;
         p->ax += 0.5*delta_a[i]/tau_a[i]*exp(-time/tau_a[i])*dvx/o.a;
         p->ay += 0.5*delta_a[i]/tau_a[i]*exp(-time/tau_a[i])*dvy/o.a; 
         p->az += 0.5*delta_a[i]/tau_a[i]*exp(-time/tau_a[i])*dvz/o.a; 
      }
   }
}




void reb_check_remove(struct reb_simulation* r){
   const int N = r->N;
   struct reb_particle com = r->particles[0];
   for (int i=N-1;i>=r->N_active;i+=-1){
      struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);
		double q = o.a*(1.-o.e);
      if(q < 5.2 || o.a > 500. || o.a < 0.05){
         printf("removing particle %u at time %e \n", r->particles[i].hash, r->t);			
         reb_remove_by_hash(r, r->particles[i].hash, 1);
      }
   }
}




void heartbeat(struct reb_simulation* r){
   if(reb_output_check(r, 1000.*r->dt)){
      reb_check_remove(r);
	}
}
