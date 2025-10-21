/**
 * example code to evolve bodies orbital elements with a desired functional form
 * the modifications are found in directory srcDev for mercurius and whfast integrators
 * the modifications are explained in Hermosillo Ruiz et. al. 2025 
 ** "Forcing Planets to Evolve: Interactions Between Uranus and Neptune at Late Stages of Dynamical Evolution"
 * compile with make -f MakefileDev
 * run with ./rebound
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "rebound.h"

// parameters to use in migration forces function
double tout;
double* tau_e;
double* tau_a;     // Migration timescale in years for all particles
double* tau_i; 
double* tau_w;
double* tau_Om;
double* final_a;     // migration change in semimajor axis for planets  
double* final_e;
double* final_i;
double* final_w;
double* final_Om;
double* idot;
double* edot;
double* adot;
double* wdot;
double* Omdot;
double tmax;
double* a0;
double* e0;
double* inc0;
double* node0;
double* w0;
char* outputfile = "sim-output.txt";

void migration_forces(struct reb_simulation* sim);
void heartbeat(struct reb_simulation* sim);

int main(int argc, char* argv[]){
	//initializing the simulation variable so it's always accessible outside if statements
	struct reb_simulation* sim = reb_create_simulation();
	reb_free_simulation(sim);
    sim = NULL;

	//reading in the basic integration parameters
	FILE *myfile;
	double myvariable;
	char* intfile = "integration-parameters.txt";
	myfile=fopen(intfile, "r");
	fscanf(myfile,"%lf",&myvariable);
    double dt = 2.0*M_PI*myvariable;
	fscanf(myfile,"%lf",&myvariable);
	double tmax = 2.0*M_PI*myvariable;
	fscanf(myfile,"%lf",&myvariable);
	tout = 2.0*M_PI*myvariable;
	fscanf(myfile,"%lf",&myvariable);
	int trial_sim = myvariable;
	fscanf(myfile,"%lf",&myvariable);
	int numP = myvariable;
	fclose(myfile);

	char* migfile = "planet-migration-parameters.txt";
	tau_a = calloc(sizeof(double),5);
	final_a = calloc(sizeof(double),5);
	final_e = calloc(sizeof(double),5);
	final_i = calloc(sizeof(double),5);
	final_w = calloc(sizeof(double),5);
	final_Om = calloc(sizeof(double),5);
	tau_e = calloc(sizeof(double),5);
	tau_i = calloc(sizeof(double),5);
	tau_w = calloc(sizeof(double),5);
	tau_Om = calloc(sizeof(double),5);
    idot = calloc(sizeof(double), 5);
    edot = calloc(sizeof(double), 5);
    adot = calloc(sizeof(double), 5);
	wdot = calloc(sizeof(double),5);
	Omdot = calloc(sizeof(double),5);
	myfile=fopen(migfile, "r");
	for(int i = 1; i < 5; i++){
		fscanf(myfile,"%lf",&myvariable);
		tau_a[i] = 2.0*M_PI*myvariable;
		fscanf(myfile,"%lf",&myvariable);
		tau_e[i] = 2.0*M_PI*myvariable;
		fscanf(myfile,"%lf",&myvariable);
		tau_i[i] = 2.0*M_PI*myvariable;
		fscanf(myfile,"%lf",&myvariable);
		tau_w[i] = 2.0*M_PI*myvariable;
		fscanf(myfile,"%lf",&myvariable);
		tau_Om[i] = 2.0*M_PI*myvariable;
		fscanf(myfile,"%lf",&myvariable);
		final_a[i] = myvariable;
		fscanf(myfile,"%lf",&myvariable);
		final_e[i] = myvariable;
		fscanf(myfile,"%lf",&myvariable);
		final_i[i] = (M_PI/180.0)*myvariable;
		fscanf(myfile,"%lf",&myvariable);
		final_w[i] = (M_PI/180.0)*myvariable;
		fscanf(myfile,"%lf",&myvariable);
		final_Om[i] = (M_PI/180.0)*myvariable;
		printf("%f; %f; %f; %f; %f\n",tau_a[i],tau_e[i],tau_i[i],tau_w[i],tau_Om[i]);
		printf("%f; %f; %f; %f; %f; %f\n",final_a[i],final_e[i],final_i[i],final_w[i],final_Om[i],M_PI);
        idot[i] = 0.0;
        edot[i] = 0.0;
        adot[i] = 0.0;
		wdot[i] = 0.0;
		Omdot[i] = 0.0;
	}
	fclose(myfile);

	char* filename = "simulationarchive.bin";
	// Trying to open a SimulationArchive file
	struct reb_simulationarchive* sa = reb_open_simulationarchive(filename);
	if (sa==NULL){
		printf("simulationarchive.bin doesn't exist\n");
	}
	else{
		// Get a simulation from the file (if possible, otherwise NULL is returned)
		printf("Reading in from a simulation archive to do a continuation of previous simulation\n");
		sim = reb_create_simulation_from_simulationarchive(sa,-1);
	}
	reb_close_simulationarchive(sa); 
	
	char* binplfile = "planet-ic.bin";
	if (sim==NULL){
		//////**************************************************************** ////////////		
		//////This is all code for either reading in a binary ic file and adding test particles
		//////or reading in planet positions from a text file and skipping test particles
		//////to just see what the planets do
		//////**************************************************************** ////////////		
		sim = reb_create_simulation_from_binary(binplfile);		
	}
	// Check if we were successful
	if (sim==NULL){
		printf("No binary initial conditions found. Testing planet initial conditions from text file.\n");
	       
		sim = reb_create_simulation();
		sim->integrator = REB_INTEGRATOR_MERCURIUS;
		sim->dt = dt;        
		sim->t = 0.;        
		sim->additional_forces = migration_forces;     
		sim->force_is_velocity_dependent = 1;

		char* txticfile = "planet-ic.txt";
		
		a0 = calloc(sizeof(double),5);
		e0 = calloc(sizeof(double),5);
		inc0 = calloc(sizeof(double),5);	
		node0 = calloc(sizeof(double),5);
		w0 = calloc(sizeof(double),5);
		double* M0 = calloc(sizeof(double),5);
		double* rpl = calloc(sizeof(double),5);
		double* masspl = calloc(sizeof(double),5);

		myfile=fopen(txticfile, "r");
		for(int i = 1; i < 5; i++){
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
			w0[i] = (M_PI/180.0)*myvariable;
			fscanf(myfile,"%lf",&myvariable);
			node0[i] = (M_PI/180.0)*myvariable;
			fscanf(myfile,"%lf",&myvariable);
			M0[i] = (M_PI/180.0)*myvariable;
		}
		fclose(myfile);

		struct reb_particle star = {0};
		star.m  = 1.0;          
		reb_add(sim, star); 
	
      	for(int i = 1; i <= numP; i++){
			struct reb_particle pl = reb_tools_orbit_to_particle(sim->G, sim->particles[0],masspl[i],a0[i],e0[i],inc0[i],node0[i],w0[i],M0[i]);
			pl.r = rpl[i];
			pl.hash = i;
			reb_add(sim, pl);
         }
		
		sim->N_active = sim->N;
		reb_move_to_com(sim);  
		//write these initial conditions to a binary file so they can be used again reproducibly
        // commenting since sim archive not working with the changes we made in src code
		//reb_simulationarchive_snapshot(r,binplfile);
	}
	if (trial_sim < 1){
		printf("Loaded planet initial conditions from txt file. Adding test particles \n");
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

		sim->rand_seed = new_seed;
		int k=101;
		printf("\nAdding %i test particles from  a=%.16f - %.16f with de=%.16f, di=%.16f.\n\n",N_tp,tpamin,tpamax,tpde,tpdi);
		for( int j = sim->N; j < N_tp + sim->N_active; j++){
			double a    = reb_random_uniform(sim,tpamin,tpamax);
			double e    = reb_random_rayleigh(sim, tpde);
			double inc  = reb_random_rayleigh(sim, tpdi);
			double Omega = reb_random_uniform(sim, 0.,2.*M_PI);
			double omega = reb_random_uniform(sim, 0.,2.*M_PI);
			double f 	= reb_random_uniform(sim, 0.,2.*M_PI);
			struct reb_particle pt = reb_tools_orbit_to_particle(sim->G, sim->particles[0], 0., a, e, inc, Omega, omega, f);
			pt.hash = k;
			k+=1;
			reb_add(sim, pt);
		}

	//////**************************************************************** ////////////		
	////// End of just setting up initial conditions!
	////////**************************************************************** ////////////	
	}
			
	printf("nactive = %i\n",sim->N_active);
	printf("ntotal = %i\n",sim->N);
	sim->integrator = REB_INTEGRATOR_MERCURIUS;
	sim->dt = dt;        
	sim->additional_forces = migration_forces; 
    sim->heartbeat = heartbeat;     
	sim->force_is_velocity_dependent = 1;
	printf("integrating to %.16f saving every %.16f time units to %s\n",tmax,tout,outputfile);
	//reb_simulationarchive_automate_interval(r,filename,tout); // sim archive does not work with modifications
	reb_integrate(sim, tmax);
}

void migration_forces(struct reb_simulation* sim){
    const int N = sim->N_active;
    struct reb_particle* const particles = sim->particles;
    struct reb_particle hel = particles[0]; // calculate migration forces with respect to star
    // struct reb_particle com = reb_get_com(sim); // or with respect to com
    for(int i=1;i<N;i++){
		struct reb_particle* p = &(particles[i]);
		p->vusrx = 0.0; 
		p->vusry = 0.0;
		p->vusrz = 0.0;
		p->ausrx = 0.0;
		p->ausry = 0.0;
		p->ausrz = 0.0;
		adot[i] = 0.0;
		idot[i] = 0.0;
		edot[i] = 0.0;
		wdot[i] = 0.0;
		Omdot[i] = 0.0;

		struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, particles[i],hel);
		const double msum = p->m + hel.m;

		double semi = o.a;
		double e = o.e;
		double f = o.f;
		double inc = o.inc;
		double node = o.Omega;
		double peri = o.omega;
		
		double delta_a=final_a[i]-a0[i];
		double de=final_e[i]-e0[i];
		double di=final_i[i]-inc0[i];
		double dw = final_w[i]-w0[i];
		double dnode = final_Om[i] - node0[i];
		
		// here define the functional form for each orbital element

		if(tau_a[i] != 0.0){
			adot[i] = delta_a*exp(-sim->t/tau_a[i])/(tau_a[i]);
		}
		if(tau_e[i]!=0.0){
			edot[i] = de*exp(-sim->t/tau_e[i])/(tau_e[i]);
		}
		if(tau_i[i]!=0.0){
			idot[i] = di*2.0*M_PI*cos(2.0*M_PI*sim->t/tau_i[i])/(tau_i[i]);
		}
		if(tau_w[i] != 0.0){
			wdot[i] = dw*2.0*M_PI*cos(2.0*M_PI*sim->t/tau_w[i])/(tau_w[i]);
		}
		if(tau_Om[i] != 0.0){
			Omdot[i] = dnode*2.0*M_PI*cos(2.0*M_PI*sim->t/tau_Om[i])/(tau_Om[i]);
		}
		
		const double r = semi * (1.0-e*e)/(1.0+e*cos(f));
		const double rfdot = sqrt(msum)*pow(semi,(-1.0/2))*(1.0+e*cos(f))*pow(1.0-e*e,-1.0/2);
		const double drda = (1.0-e*e)/(1.0+e*cos(f));
		const double drde = -2.0*semi*e/(1.0+e*cos(f))- semi*(1.0-e*e)*cos(f)/pow(1.0+e*cos(f),2.0);
		const double drdf = semi*(1.0-e*e)*e*sin(f)/pow(1.0+e*cos(f),2.0);
		const double drfdotda = -rfdot/(2.0*semi);
		const double drfdotde = rfdot*(e+cos(f))/((1.0-e*e)*(1.0+e*cos(f)));
		// commenting since not forcing true anomoly to evolve
		// const double drfdotdf=-semi*sqrt(msum/pow(semi,3.0))*e*sin(f)/sqrt(1.0-e*e);
//
// *****************    assuming keplerian rdot, use these (only f is time dependent)

//                 fdot_f = 0.0
//                 rdot = sqrt(msum/semi**3.0)*semi*e*sin(f)/sqrt(1.0-e*e)
//                 dfdotda = 0.0
//                 dfdotde = 0.0
//                 dfdotdf = 0.0
//                 drdotda = (-1.0/2.0)*rdot/(semi)
//                 drdotde = rdot/(e*(1.0-e*e))
//                 drdotdf = sqrt(msum/semi**3.0)*semi*e*cos(f)/sqrt(1.0-e*e)
//
//*****************  now a,e,f are time dependent
//
		const double fdot_f = rfdot/r;
		const double rdot = drda*adot[i]+drde*edot[i]+drdf*fdot_f; 
		const double dfdotda = -(3.0/2.0)*sqrt(msum/pow(semi,5.0))*(pow(1.0+e*cos(f),2.0))/pow(1.0-e*e,3.0/2.0);
		const double dfdotde = sqrt(msum/pow(semi,3.0))*(2.0*(1.0+e*cos(f))*cos(f)/pow(1.0-e*e,3.0/2.0)\
		+(3.0)*(pow(1.0+e*cos(f),2.0))*e/(pow(1.0-e*e,5.0/2.0)));
		// const double dfdotdf = -sqrt(msum/pow(semi,3.0))*2.0*(1.0+e*cos(f))*e*sin(f)/(pow(1.0-e*e,3.0/2.0));
		const double drdotda = (-2.0*e/(1.0+e*cos(f))-(1.0-e*e)*cos(f)/pow(1.0+e*cos(f),2.0))*edot[i]\
		+((1.0-e*e)*e*sin(f)/(pow(1.0+e*cos(f),2.0)))*fdot_f\
		+(semi*(1.0-e*e)*e*sin(f)/pow(1.0+e*cos(f),2.0))*dfdotda;
		const double drdotde = (-2.0*e/(1.0+e*cos(f))- (1.0-e*e)*cos(f)/pow(1.0+e*cos(f),2.0))*adot[i]\
		+(-2.0*semi/(1.0+e*cos(f))+4.0*semi*e*cos(f)/pow(1.0+e*cos(f),2.0)\
		+(2.0*semi*(1.0-e*e)*pow(cos(f),2.0))/pow(1.0+e*cos(f),3.0))*edot[i]\
		+(semi*(1.0-3.0*e*e)*sin(f)/pow(1.0+e*cos(f),2.0)\
		-2.0*semi*(e-e*e*e)*sin(f)*cos(f)/pow(1.0+e*cos(f),3.0))*fdot_f\
		+(semi*(1.0-e*e)*e*sin(f)/pow(1.0+e*cos(f),2.0))*dfdotde;
		// commenting since not forcing true anomoly to evolve
		// const double drdotdf = ((1.0-e*e)*e*sin(f)/pow(1.0+e*cos(f),2.0))*adotovera[i]*semi\
		// +(-2.0*semi*e*e*sin(f)/pow(1.0+e*cos(f),2.0)\
		// +semi*(1.0-e*e)*sin(f)/pow(1.0+e*cos(f),2.0)\
		// -2.0*semi*(1.0-e*e)*cos(f)*e*sin(f)/pow(1.0+e*cos(f),3.0))*edotovere[i]*e\
		// +(semi*(1.0-e*e)*e*cos(f)/pow(1.0+e*cos(f),2.0)\
		// +(2.0*semi*(1.0-e*e)*e*e*pow(sin(f),2.0))/pow(1.0+e*cos(f),3.0))*fdot_f\
		// +(semi*(1.0-e*e)*e*sin(f)/pow(1.0+e*cos(f),2.0))*dfdotdf;
// **************************

		const double dxda = drda*cos(node)*cos(peri+f)-drda*cos(inc)*sin(node)*sin(peri+f);
		const double dxde = drde*cos(node)*cos(peri+f)- drde*cos(inc)*sin(node)*sin(peri+f);
		const double dxdi = r*sin(inc)*sin(node)*sin(peri+f);
		const double dxdOm = -r*sin(node)*cos(peri+f)-r*cos(inc)*cos(node)*sin(peri+f);
		const double dxdw =-r*cos(node)*sin(peri+f)-r*cos(inc)*sin(node)*cos(peri+f);
		// commenting since not forcing true anomoly to evolve		
		// const double dxdf = drdf*cos(node)*cos(peri+f)-r*cos(node)*sin(peri+f)\
		// -drdf*cos(inc)*sin(node)*sin(peri+f)-r*cos(inc)*sin(node)*cos(peri+f);
		
		const double dyda = drda*sin(node)*cos(peri+f)+drda*cos(inc)*cos(node)*sin(peri+f) ;
		const double dyde=drde*sin(node)*cos(peri+f)+drde*cos(inc)*cos(node)*sin(peri+f);
		const double dydi = -r*sin(inc)*cos(node)*sin(peri+f);
		const double dydOm=r*cos(node)*cos(peri+f)-r*cos(inc)*sin(node)*sin(peri+f);
		const double dydw= -r*sin(node)*sin(peri+f)+r*cos(inc)*cos(node)*cos(peri+f);
		// commenting since not forcing true anomoly to evolve
		// const double dydf = drdf*sin(node)*cos(peri+f)-r*sin(node)*sin(peri+f)\
		// +drdf*cos(inc)*cos(node)*sin(peri+f)+r*cos(inc)*cos(node)*cos(peri+f);

		const double dzda = drda*sin(inc)*sin(peri+f);
		const double dzde = drde*sin(inc)*sin(peri+f);
		const double dzdi = r*cos(inc)*sin(peri+f);
		const double dzdOm = 0.0;
		const double dzdw = r*sin(inc)*cos(peri+f);
		// commenting since not forcing true anomoly to evolve
		// const double dzdf= drdf*sin(inc)*sin(peri+f)+ r*sin(inc)*cos(peri+f);

		const double dxdotda = drdotda*cos(node)*cos(peri+f)\
		-drda*sin(node)*Omdot[i]*cos(peri+f)\
		-(drda*wdot[i]+drfdotda)*cos(node)*sin(peri+f)\
		-drdotda*cos(inc)*sin(node)*sin(peri+f)\
		+drda*sin(inc)*idot[i]*sin(node)*sin(peri+f)\
		-drda*cos(inc)*cos(node)*Omdot[i]*sin(peri+f)\
		-(drda*wdot[i]+drfdotda)*cos(inc)*sin(node)*cos(peri+f);

		const double dxdotde = drdotde*cos(node)*cos(peri+f)\
		-drde*sin(node)*Omdot[i]*cos(peri+f)\
		-(drde*wdot[i]+drfdotde)*cos(node)*sin(peri+f)\
		-drdotde*cos(inc)*sin(node)*sin(peri+f)\
		+drde*sin(inc)*idot[i]*sin(node)*sin(peri+f)\
		-drde*cos(inc)*cos(node)*Omdot[i]*sin(peri+f)\
		-(drde*wdot[i]+drfdotde)*cos(inc)*sin(node)*cos(peri+f);
		
		const double dxdotdi = rdot*sin(inc)*sin(node)*sin(peri+f)\
		+r*cos(inc)*idot[i]*sin(node)*sin(peri+f)\
		+r*sin(inc)*cos(node)*Omdot[i]*sin(peri+f)+(r*wdot[i]\
		+rfdot)*sin(inc)*sin(node)*cos(peri+f);

		const double dxdotdOm = -rdot*sin(node)*cos(peri+f)\
		-r*cos(node)*Omdot[i]*cos(peri+f)\
		+(r*wdot[i]+rfdot)*sin(node)*sin(peri+f)\
		-rdot*cos(inc)*cos(node)*sin(peri+f)\
		+r*sin(inc)*idot[i]*cos(node)*sin(peri+f)\
		+r*cos(inc)*sin(node)*Omdot[i]*sin(peri+f)\
		-(r*wdot[i]+rfdot)*cos(inc)*cos(node)*cos(peri+f);

		const double dxdotdw = -rdot*cos(node)*sin(peri+f)\
		+r*sin(node)*Omdot[i]*sin(peri+f)\
		-(r*wdot[i]+rfdot)*cos(node)*cos(peri+f)\
		-rdot*cos(inc)*sin(node)*cos(peri+f)\
		+r*sin(inc)*idot[i]*sin(node)*cos(peri+f)\
		-r*cos(inc)*cos(node)*Omdot[i]*cos(peri+f)\
		+(r*wdot[i]+rfdot)*cos(inc)*sin(node)*sin(peri+f);

		// commenting since not forcing true anomoly to evolve

		// const double dxdotdf = -rdot*cos(node)*sin(peri+f)\
		// +r*sin(node)*OmdotoverOm[i]*node*sin(peri+f)\
		// -(r*wdotoverw[i]*peri+rfdot)*cos(node)*cos(peri+f)\
		// -rdot*cos(inc)*sin(node)*cos(peri+f)\
		// +r*sin(inc)*idotoveri[i]*inc*sin(node)*cos(peri+f)\
		// -r*cos(inc)*cos(node)*OmdotoverOm[i]*node*cos(peri+f)\
		// +(r*wdotoverw[i]*peri+rfdot)*cos(inc)*sin(node)*sin(peri+f);
		// +drdotdf*cos(node)*cos(peri+f)\
		// -drdf*sin(node)*OmdotoverOm[i]*node*cos(peri+f)\
		// -(drdf*wdotoverw[i]*peri+drfdotdf)*cos(node)*sin(peri+f)\
		// -drdotdf*cos(inc)*sin(node)*sin(peri+f) \
		// +drdf*sin(inc)*idotoveri[i]*inc*sin(node)*sin(peri+f)\
		// -drdf*cos(inc)*cos(node)*OmdotoverOm[i]*node*sin(peri+f)\
		// -(drdf*wdotoverw[i]*peri+drfdotdf)*cos(inc)*sin(node)*cos(peri+f);              

		const double dydotda = drdotda*sin(node)*cos(peri+f)\
		+drda*cos(node)*Omdot[i]*cos(peri+f)\
		-(drda*wdot[i]+drfdotda)*sin(node)*sin(peri+f)\
		+drdotda*cos(inc)*cos(node)*sin(peri+f)\
		-drda*sin(inc)*idot[i]*cos(node)*sin(peri+f)\
		-drda*cos(inc)*sin(node)*Omdot[i]*sin(peri+f)\
		+(drda*wdot[i]+drfdotda)*cos(inc)*cos(node)*cos(peri+f);

		const double dydotde=drdotde*sin(node)*cos(peri+f)\
		+drde*cos(node)*Omdot[i]*cos(peri+f)\
		-(drde*wdot[i]+drfdotde)*sin(node)*sin(peri+f)\
		+drdotde*cos(inc)*cos(node)*sin(peri+f)\
		-drde*sin(inc)*idot[i]*cos(node)*sin(peri+f)\
		-drde*cos(inc)*sin(node)*Omdot[i]*sin(peri+f)\
		+(drde*wdot[i]+drfdotde)*cos(inc)*cos(node)*cos(peri+f);

		const double dydotdi=-rdot*sin(inc)*cos(node)*sin(peri+f)\
		-r*cos(inc)*idot[i]*cos(node)*sin(peri+f)\
		+r*sin(inc)*sin(node)*Omdot[i]*sin(peri+f)\
		-(r*wdot[i]+rfdot)*sin(inc)*cos(node)*cos(peri+f);

		const double dydotdOm=rdot*cos(node)*cos(peri+f)\
		-r*sin(node)*Omdot[i]*cos(peri+f)\
		-(r*wdot[i]*peri+rfdot)*cos(node)*sin(peri+f)\
		-rdot*cos(inc)*sin(node)*sin(peri+f)\
		+r*sin(inc)*idot[i]*sin(node)*sin(peri+f)\
		-r*cos(inc)*cos(node)*Omdot[i]*sin(peri+f)\
		-(r*wdot[i]+rfdot)*cos(inc)*sin(node)*cos(peri+f);

		const double dydotdw=-rdot*sin(node)*sin(peri+f)\
		-r*cos(node)*Omdot[i]*sin(peri+f)\
		-(r*wdot[i]+rfdot)*sin(node)*cos(peri+f)\
		+rdot*cos(inc)*cos(node)*cos(peri+f)\
		-r*sin(inc)*idot[i]*cos(node)*cos(peri+f)\
		-r*cos(inc)*sin(node)*Omdot[i]*cos(peri+f)\
		-(r*wdot[i]+rfdot)*cos(inc)*cos(node)*sin(peri+f);

		// commenting since not forcing true anomoly to evolve

		// const double dydotdf= -rdot*sin(node)*sin(peri+f)\
		// -r*cos(node)*OmdotoverOm[i]*node*sin(peri+f)\
		// -(r*wdotoverw[i]*peri+rfdot)*sin(node)*cos(peri+f)\
		// +rdot*cos(inc)*cos(node)*cos(peri+f)\
		// -r*sin(inc)*idotoveri[i]*inc*cos(node)*cos(peri+f)\
		// -r*cos(inc)*sin(node)*OmdotoverOm[i]*node*cos(peri+f)\
		// -(r*wdotoverw[i]*peri+rfdot)*cos(inc)*cos(node)*sin(peri+f)\
		// +drdotdf*sin(node)*cos(peri+f)\
		// +drdf*cos(node)*OmdotoverOm[i]*node*cos(peri+f)\
		// -(drdf*wdotoverw[i]*peri+drfdotdf)*sin(node)*sin(peri+f)\
		// +drdotdf*cos(inc)*cos(node)*sin(peri+f)\
		// -drdf*sin(inc)*idotoveri[i]*inc*cos(node)*sin(peri+f)\
		// -drdf*cos(inc)*sin(node)*OmdotoverOm[i]*node*sin(peri+f)\
		// +(drdf*wdotoverw[i]*peri+drfdotdf)*cos(inc)*cos(node)*cos(peri+f);

		const double dzdotda = drdotda*sin(inc)*sin(peri+f)\
		+drda*cos(inc)*idot[i]*sin(peri+f)\
		+(drda*wdot[i]+drfdotda)*sin(inc)*cos(peri+f);

		const double dzdotde = drdotde*sin(inc)*sin(peri+f)\
		+drde*cos(inc)*idot[i]*sin(peri+f)\
		+(drde*wdot[i]+drfdotde)*sin(inc)*cos(peri+f);

		const double dzdotdi = rdot*cos(inc)*sin(peri+f)\
		-r*sin(inc)*idot[i]*sin(peri+f)\
		+(r*wdot[i]+rfdot)*cos(inc)*cos(peri+f);

		const double dzdotdOm=0.0;

		const double dzdotdw=rdot*sin(inc)*cos(peri+f)\
		+r*cos(inc)*idot[i]*cos(peri+f)\
		-(r*wdot[i]+rfdot)*sin(inc)*sin(peri+f);

		// commenting since not forcing true anomoly to evolve

		// const double dzdotdf=rdot*sin(inc)*cos(peri+f)\
		// +r*cos(inc)*idotoveri[i]*inc*cos(peri+f)\
		// -(r*wdotoverw[i]*peri+rfdot)*sin(inc)*sin(peri+f)\
		// +drdotdf*sin(inc)*sin(peri+f)\
		// +drdf*cos(inc)*idotoveri[i]*inc*sin(peri+f)\
		// +(drdf*wdotoverw[i]*peri+drfdotdf)*sin(inc)*cos(peri+f);

		if (adot[i] != 0.0){
			p->ausrx+= dxdotda*adot[i];
			p->ausry+= dydotda*adot[i];
			p->ausrz+= dzdotda*adot[i];
			p->vusrx+= dxda*adot[i];
			p->vusry+= dyda*adot[i];
			p->vusrz+= dzda*adot[i];
		}

		if (edot[i] != 0.0){
			p->ausrx+= dxdotde*edot[i];
			p->ausry+= dydotde*edot[i];
			p->ausrz+= dzdotde*edot[i];
			p->vusrx+= dxde*edot[i];
			p->vusry+= dyde*edot[i];
			p->vusrz+= dzde*edot[i];
		}

		if (idot[i]!=0.0) {
			p->ausrx+= dxdotdi*idot[i];
			p->ausry+= dydotdi*idot[i];
			p->ausrz+= dzdotdi*idot[i];
			p->vusrx+= dxdi*idot[i];
			p->vusry+= dydi*idot[i];
			p->vusrz+= dzdi*idot[i];
		}

		if (wdot[i] != 0.0) {
			p->ausrx+= dxdotdw*wdot[i];
			p->ausry+= dydotdw*wdot[i];
			p->ausrz+= dzdotdw*wdot[i];
			p->vusrx+= dxdw*wdot[i];
			p->vusry+= dydw*wdot[i];
			p->vusrz+= dzdw*wdot[i];
		}

		if (Omdot[i]!=0.0){ 
			p->ausrx+= dxdotdOm*Omdot[i];
			p->ausry+= dydotdOm*Omdot[i];
			p->ausrz+= dzdotdOm*Omdot[i];
			p->vusrx+= dxdOm*Omdot[i];
			p->vusry+= dydOm*Omdot[i];
			p->vusrz+= dzdOm*Omdot[i];
		}
    }
}


void reb_output_all_orbits(struct reb_simulation* sim, char* filename){
    const int N = sim->N;
#ifdef MPI
    char filename_mpi[1024];
    sprintf(filename_mpi,"%s_%d",filename,sim->mpi_id);
    FILE* of = fopen(filename_mpi,"a");
#else // MPI
    FILE* of = fopen(filename,"a");
#endif // MPI
	struct reb_particle* const particles = sim->particles;
	struct reb_particle hel = particles[0];
    for (int i=1;i<N;i++){
		// struct reb_particle* p = &(particles[i]);
        struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[i],hel);
		fprintf(of,"%u\t%20.16e\t%20.16e\t%20.16e\t%20.16e\t%20.16e\t%20.16e\t%20.16e\n",i,sim->t,o.a,o.e,o.inc,o.omega,o.Omega,o.f);
    }
    fclose(of);
}


void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, tout)){ // tout is already in 2*pi
		reb_integrator_synchronize(sim);
      	reb_output_all_orbits(sim,outputfile);
     	reb_move_to_com(sim); 
    }
}
