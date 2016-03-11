//Initialize velocities of all atoms
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 200 //number of atoms

int main(){

    const double m=40.0*1.66053892173; //mass of a single atom of argon (1E-27kg)
    const double T0=300.0; //initial temperature given (K)
    const double k=1.3806488E-6; //boltzmann constant(1E-17(m^2)(kg)(s^-2)(K^-1))

	int i=0;
	double vel[N][2]; //velocity of each atom

	//generate velocitie by random number
	double L=0.005; //range of random numbers (1E5 m/s)
	srand(time(NULL));

	for(i=0;i<N;i++){
		vel[i][0]=L*(rand()/((double)RAND_MAX));
		vel[i][1]=L*(rand()/((double)RAND_MAX));
	}

	//zero momentum
	double avg_vel[2]; //average velocity

	for(i=0;i<N;i++){
		avg_vel[0]+=vel[i][0]/(double)N;
		avg_vel[1]+=vel[i][1]/(double)N;
	}
	for(i=0;i<N;i++){
		vel[i][0]-=avg_vel[0];
		vel[i][1]-=avg_vel[1];
	}

	//Initialize temperature
	double Ek=0.0; //Total kinetic energy
    double T=0.0;
	double mod_T=0.0;
	for(i=0;i<N;i++){
		Ek+=0.5*m*(vel[i][0]*vel[i][0]+vel[i][1]*vel[i][1]);
	}

	T=Ek/N/k; //Calculate temperature now
	mod_T=sqrt(T0/T);

	for(i=0;i<N;i++){
		vel[i][0]*=mod_T;
		vel[i][1]*=mod_T;
	}

	//Print
	avg_vel[0]=0.0;
	avg_vel[1]=0.0;

	for(i=0;i<N;i++){
		avg_vel[0]+=vel[i][0]/(double)N;
		avg_vel[1]+=vel[i][1]/(double)N;
	}

	Ek=0.0;
	T=0.0;
	for(i=0;i<N;i++){
		Ek+=0.5*m*(vel[i][0]*vel[i][0]+vel[i][1]*vel[i][1]);
	}
	T=Ek/N/k; //Calculate temperature now
	//Verify the initialized temperature and velocity
	printf("Total temperature (K) : %f\n", T);
	printf("Average velocity (m/s) : %f\t%f\n", avg_vel[0], avg_vel[1]); //should be zeros
	FILE *pFile;
	pFile=fopen("velocity.txt","w");
	for(i=0;i<N;i++){
			fprintf(pFile,"%f\t%f\n", vel[i][0], vel[i][1]);
	}
	fclose(pFile);

	printf("Done!\n");

	return 0;

}
