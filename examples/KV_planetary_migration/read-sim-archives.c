/**
 * Simulation Archive
 *
 * This example shows how to use the Simulation Archive.
 * We integrate a two planet system forward in time using
 * the WHFast integrator. The simulation can be interrupted
 * at any time. On the next run, the program will try to reload
 * the latest data from the Simulation Archive. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "rebound.h"

void output_orbits_with_phi(struct reb_simulation* r, char* filename, double pr, double qr, int planet);

int main(int argc, char* argv[]) {
    char* filename = "simulationarchive.bin";
    char* filename2 = "high-res-end-archive.bin";

    // Trying to open a SimulationArchive file
    struct reb_simulationarchive* sa = reb_open_simulationarchive(filename);
    if (sa==NULL){
        printf("Can not open file.\n");
    }


	int nsnaps = sa->nblobs;
	double pr = 2.;
	double qr = 1.;
	int planet = 4;
	
	for( int j = 0; j < nsnaps; j++)
		{
		struct reb_simulation* r = reb_create_simulation_from_simulationarchive(sa,j);
		output_orbits_with_phi(r,"c-all-orbits.txt",pr,qr,planet);
    	reb_free_simulation(r);
		}

    struct reb_simulationarchive* sa2 = reb_open_simulationarchive(filename2);
    if (sa2==NULL){
        printf("Can not open file.\n");
    }


	nsnaps = sa2->nblobs;
	
	for( int j = 0; j < nsnaps; j++)
		{
		struct reb_simulation* r = reb_create_simulation_from_simulationarchive(sa2,j);
		output_orbits_with_phi(r,"c-end-orbits.txt",pr,qr,planet);
    	reb_free_simulation(r);
		}



}

void output_orbits_with_phi(struct reb_simulation* r, char* filename, double pr, double qr, int planet){
    const int N = r->N;
    FILE* of = fopen(filename,"a");
	 struct reb_particle com = reb_get_com(r);
    double time = r->t/(2.*M_PI);
	 double lpl = 0.;
    for (int i=1;i<N;i++){
        struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);
		  double lam = o.Omega + o.omega + o.M;
		  double wt = o.Omega + o.omega;
		  if(i==planet)
			  	{
				lpl = lam;
				}
		  double resangle = pr*lam - qr*lpl - wt;
		  while(resangle < 0){resangle = resangle+2.0*M_PI;}
        while(resangle > 2.0*M_PI){resangle = resangle-2.0*M_PI;}
        fprintf(of,"%u\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",r->particles[i].hash,time,o.a,o.e,o.inc,o.Omega,o.omega,o.M,resangle,lam,lpl);
    }
    fclose(of);
}
