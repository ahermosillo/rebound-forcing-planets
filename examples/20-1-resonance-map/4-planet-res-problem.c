/**
 * 4 giant planet resonance mapping simulation
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

double fprev[5000];
struct reb_particle com;

double pr; 
double qr; 
double a0res;


int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
    // Setup constants
	r->integrator   = REB_INTEGRATOR_MERCURIUS;
	r->collision_resolve_keep_sorted = 1;
	r->collision = REB_COLLISION_DIRECT;
	r->dt           = 0.25; //3e-3*2.*M_PI;
	r->heartbeat    = heartbeat;
	r->t = 0.;


	FILE *myfile;
	double myvariable;
	double qperi;
	double inc;
	double node;
	double ncycles;
	int N_tp; // =360;
	int N_tpr; // = 360;
	myfile=fopen("resonance-parameters.txt", "r");
	fscanf(myfile,"%lf",&myvariable);
	pr = myvariable;
	fscanf(myfile,"%lf",&myvariable);
	qr = myvariable;
	fscanf(myfile,"%lf",&myvariable);
	qperi = myvariable;
	
	fscanf(myfile,"%lf",&myvariable);
	inc = myvariable*M_PI/180.0;
	
	fscanf(myfile,"%lf",&myvariable);
	node = myvariable*M_PI/180.0;
	
	fscanf(myfile,"%lf",&myvariable);
	ncycles = myvariable;
	fscanf(myfile,"%lf",&myvariable);
	N_tpr = myvariable;
	fscanf(myfile,"%lf",&myvariable);
	N_tp = myvariable;

	printf("\n res particles %d\n",N_tpr);
	printf("\n particles %d\n",N_tp);
	
	fclose(myfile);
	
	double phimax = 2.*M_PI/qr;
	double phi0 = 0.;
	
	double dphi = phimax/N_tpr;

	
	if(qr > 2){
	    phi0 = phimax;
	    phimax = phimax+phimax;
	}
	




	//add the star
	struct reb_particle star = {0};
	star.m = 1.0; 
	reb_add(r, star);
	
	//add Planets from Horizons data for JD2459200.500000000
	//from orbit it's G, primary,mass, a, e, inc, node, arg peri, true anomaly
	//though adding here in cartesian coords
	struct reb_particle jup = {0};
    jup.x = 2.951175666074504E+00;
    jup.y = -4.160109295473768E+00;
    jup.z = -4.874924149414127E-02;
    jup.vx = 6.070224750969053E-03*365.25/(2.0*M_PI);
    jup.vy = 4.727007450896939E-03*365.25/(2.0*M_PI);
    jup.vz = -1.554455127043464E-04*365.25/(2.0*M_PI);
    jup.m = 9.545823E-04;
    double rhill =  5.2*pow(jup.m/(3.*star.m),1./3.);
    jup.r = rhill;
    reb_add(r, jup);
    
    struct reb_particle sat = {0};
    sat.x = 5.424825699753957E+00;
    sat.y = -8.387710923072389E+00;
    sat.z = -7.007721520557159E-02;
    sat.vx = 4.381055765854487E-03*365.25/(2.0*M_PI);
    sat.vy = 3.020425159647322E-03*365.25/(2.0*M_PI);
    sat.vz = -2.268231137361701E-04*365.25/(2.0*M_PI);
    sat.m = 2.8580206E-04;
    rhill =  9.2*pow(sat.m/(3.*star.m),1./3.);
    sat.r = rhill;
    reb_add(r, sat);

    struct reb_particle ura = {0};
    ura.x = 1.538750382669396E+01;
    ura.y = 1.241842514052818E+01;
    ura.z = -1.532205483223961E-01;
    ura.vx = -2.493455589404909E-03*365.25/(2.0*M_PI);
    ura.vy = 2.881947144487664E-03*365.25/(2.0*M_PI);
    ura.vz = 4.285848895278896E-05*365.25/(2.0*M_PI);
    ura.m = 4.365738E-05;
    rhill =  19.2*pow(ura.m/(3.*star.m),1./3.);
    ura.r = rhill;
    reb_add(r, ura);
    
    struct reb_particle nep = {0};
    nep.x = 2.945247244787225E+01;
    nep.y = -5.278545608764948E+00;
    nep.z = -5.701377712423963E-01;
    nep.vx = 5.393961892265245E-04*365.25/(2.0*M_PI);
    nep.vy = 3.114249566999407E-03*365.25/(2.0*M_PI);
    nep.vz = -7.650435608534635E-05*365.25/(2.0*M_PI);
    nep.m = 5.150244E-05;
    rhill =  30.1*pow(nep.m/(3.*star.m),1./3.);
    nep.r = rhill;
    reb_add(r, nep);
	
	
	
	double tmax = ncycles*pr*2.*M_PI*164.79132; 

	//struct reb_particle p = reb_tools_orbit_to_particle(r->G, r->particles[0], mubar, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	//double rhill =  pow(p.m/(3.*star.m),1./3.);
    //set the collision radius to 1 hill sphere to save time
	//p.r = rhill;
	//reb_add(r, p);

	r->N_active = r->N;
    
	reb_move_to_com(r);
	com = reb_get_com(r);
	
	struct reb_orbit on = reb_tools_particle_to_orbit(r->G, r->particles[4],com);
	//calculate the exact p:q resonant semimajor axis
	a0res = on.a*pow((pr/qr),(2./3.));
	printf("\n a0res =%.16f.\n\n",a0res);

	double ecc = 1. - qperi/a0res;
    printf("\n eccentricity=%.16f.\n\n",ecc);



	double count = 0.;
	double count3 = 0.;
	for( int j = r->N_active; j < N_tp + r->N_active; j++)
		{
		if(count > N_tpr && qr == 3)
			{
			count = count+1.0;
			double deltaa = 0.8/( (float)N_tp - (float)N_tpr);
			double atp = (a0res - 0.4) + deltaa*(float)count3;
			count3 = count3+1; 
			double necc = 1. - qperi/(atp);
			double phi = 2.98; //M_PI;
			double argp = phi-node;
			struct reb_particle ptr = reb_tools_orbit_to_particle(r->G,com,0., atp, necc, inc, node, argp, 0.0);
		    	ptr.hash = 99+j;
		    	fprev[j] = 3.0;
		    	reb_add(r, ptr);
			}
		else if(count > N_tpr)
			{
			count = count+1.0;
			double rand = 0.01;
			if(a0res > 150.)
				{rand = reb_random_uniform(r, -1,1);}
			else
				{rand = reb_random_uniform(r, -0.5,0.5);}
			double atp = a0res + rand; 
			double necc = 1. - qperi/(atp);
			double phi = reb_random_uniform(r, 0., phimax);
			double argp = phi-node;
			struct reb_particle ptr = reb_tools_orbit_to_particle(r->G,com,0., atp, necc, inc, node, argp, 0.0);
		    ptr.hash = 99+j;
		    fprev[j] = 3.0;
		    reb_add(r, ptr);
			}
		else
		    {
		    double phi = phi0 + count*dphi; 
		    count = count+1.0;
		    double atp = a0res;
			double argp = phi-node;
		    struct reb_particle pt = reb_tools_orbit_to_particle(r->G,com,0., atp, ecc, inc, node, argp, 0.0);
		    pt.hash = 99+j;
		    fprev[j] = 3.0;
		    reb_add(r, pt);
		    }

    	}

	const int N = r->N;
	FILE* of = fopen("sections.txt","a");

	//output the initial conditions
	for (int i=N-1;i>=r->N_active;i+=-1)
		{
		struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);
		on = reb_tools_particle_to_orbit(r->G, r->particles[4],com);
	    double lnep = on.M + on.Omega + on.omega;
        while(lnep < 0){lnep = lnep+2.0*M_PI;}
        while(lnep > 2.0*M_PI){lnep = lnep-2.0*M_PI;}		
		double lkbo = o.M + o.Omega + o.omega;
		double resangle = pr*lkbo - qr*lnep - (pr-qr)*(o.Omega + o.omega);
		double psi = lkbo-lnep;
        while(psi < 0){psi = psi+2.0*M_PI;}
        while(psi > 2.0*M_PI){psi = psi-2.0*M_PI;}		
        while(resangle < 0){resangle = resangle+2.0*M_PI;}
        while(resangle > 2.0*M_PI){resangle = resangle-2.0*M_PI;}		
		fprintf(of,"%u\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",r->particles[i].hash,r->t,o.a,psi,
				o.e,o.inc,o.Omega,o.omega,o.M,lnep,resangle);
		}
	fclose(of);
	
 	//run the integration	
    printf("\n integrate to t=%.16f.\n\n",tmax);
	reb_integrate(r, tmax);
	}

void heartbeat(struct reb_simulation* r)
	{
	struct reb_orbit on = reb_tools_particle_to_orbit(r->G, r->particles[4],com);
	double lnep = on.M + on.Omega + on.omega;
    while(lnep < 0){lnep = lnep+2.0*M_PI;}
    while(lnep > 2.0*M_PI){lnep = lnep-2.0*M_PI;}
	FILE* of = fopen("sections.txt","a"); 
	double time = r->t/(2.*M_PI);
	const int N = r->N;
		
	for (int i=N-1;i>=r->N_active;i+=-1)
		{
		struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);
		double MA = o.M;
		if(MA < 0){MA = MA+2.0*M_PI;}
		double dar = fabs(o.a - a0res);
		if(MA < 1.0 && fprev[i] > 5.0)
			{
		    double lkbo = o.M + o.Omega + o.omega;
		    double resangle = pr*lkbo - qr*lnep - (pr-qr)*(o.Omega + o.omega);
		    double psi = lkbo-lnep;
            while(psi < 0){psi = psi+2.0*M_PI;}
            while(psi > 2.0*M_PI){psi = psi-2.0*M_PI;}		
            while(resangle < 0){resangle = resangle+2.0*M_PI;}
            while(resangle > 2.0*M_PI){resangle = resangle-2.0*M_PI;}		
		    fprintf(of,"%u\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",r->particles[i].hash,time,o.a,psi,
				o.e,o.inc,o.Omega,o.omega,o.M,lnep,resangle);
			}
			
		//remove particles that are too far from the resonance
		if(dar > 10.)
			{
			printf("removing particle %u\n",r->particles[i].hash);
			reb_remove_by_hash(r, r->particles[i].hash, 1);
			}
 		}
    	
		for (int i=r->N-1;i>=r->N_active;i+=-1)
			{
			struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);
			double MA = o.M;
			if(MA < 0){MA = MA+2.0*M_PI;}
			fprev[i] = MA;
			}

	fclose(of);
	}

