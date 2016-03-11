//2D molecular dynamics solver of argons
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NP 200 //number of atoms
#define NT 1000 //# of steps
const double m=40.0*1.66053892173; //mass of each atom (1E-27kg)
const double dt = 1.0;//1E-15(s)
const double L = 40.0;//1E-10(m)
const double s = 3.4;//1E-10(m)
const double esp = 0.01*1.602E-2;//1E-17(J)
const double kb=1.3806488E-6; //boltzmann constant(1E-17(m^2)(kg)(s^-2)(K^-1))
const double a=12.0;//cut-off distance, 1E-10(m)
//constants of predictor-corrector method
const double a0=3.0/16.0;
const double a1=251.0/360.0;
const double a2=1.0;
const double a3=11.0/18.0;
const double a4=1.0/6.0;
const double a5=1.0/60.0;

//Loading initial positions and velocities
void Load(double *x, double *y, double *vx, double *vy){

	int i=0;
	FILE *sheet;
	sheet = fopen("coords.txt", "r");
	if (!sheet) {
		printf("Coords loading failed!\n");
		exit(1);
	}

	for(i=0;i<NP;i++){
		fscanf(sheet, "%lE\t%lE\n", &x[i], &y[i]);
	}

	fclose(sheet);

	sheet = fopen("velocity.txt", "r");
	if (!sheet) {
		printf("Velocity loading failed!\n");
		exit(1);
	}

	for(i=0;i<NP;i++){
		fscanf(sheet, "%lE\t%lE\n", &vx[i], &vy[i]);
	}

	fclose(sheet);
}

void CalcTemp(int k, double *T, double vx[NP], double vy[NP]){
	int i=0.0;
	double Ek=0.0;
	for(i=0;i<NP;i++){
		Ek+=0.5*m*(vx[i]*vx[i] + vy[i]*vy[i]);
	}
	T[k]=Ek/NP/kb; //Calculate temperature now
}

void CalcMomentum(int k, double *P, double *H, double x[NP], double y[NP], double vx[NP], double vy[NP]){
	int i=0;
	double vxsum=0.0, vysum=0.0, rotsum=0.0;

	for(i=0;i<NP;i++){
		vxsum += vx[i];
		vysum += vy[i];
		rotsum += vy[i]*x[i] - vx[i]*y[i];
	}

	P[k] = m*sqrt(vxsum*vxsum + vysum*vysum);
	H[k] = m*rotsum;
}

void CalcEnergy(int k, double *U, double x[NP], double y[NP]){
	int i=0, j=0;
	double rx=0.0, ry=0.0;
	double r=0.0;
	double sr=0.0;
	double kr=0.0;

	for(i=0;i<NP;i++){
		for(j=i+1;j<NP;j++){
			rx = x[i]-x[j];
			ry = y[i]-y[j];
			//Periodic BC
			if(rx<-L/2.0){rx+=L;}
			else if(rx>L/2.0){rx-=L;}
			if(ry<-L/2.0){ry+=L;}
			else if(ry>L/2.0){ry-=L;}

			r=sqrt(rx*rx + ry*ry);
			if(r<a){
			sr = s/r;
			U[k] += 4*esp*(pow(sr,12)-pow(sr,6));
			}
		}
	}
}

void CalcForce(double *Fx, double *Fy, double x[NP], double y[NP]){
	int i=0, j=0;
	double rx=0.0, ry=0.0;
	double r=0.0;
	double sr=0.0;
	double kr=0.0;

	for(i=0;i<NP;i++){
		Fx[i]=0.0;
		Fy[i]=0.0;
	}

	for(i=0;i<NP;i++){
		for(j=i+1;j<NP;j++){
			rx = x[i]-x[j];
			ry = y[i]-y[j];
			//Periodic BC
			if(rx<-L/2.0){rx+=L;}
			else if(rx>L/2.0){rx-=L;}
			if(ry<-L/2.0){ry+=L;}
			else if(ry>L/2.0){ry-=L;}

			r=sqrt(rx*rx + ry*ry);
			if(r<a){
			sr = s/r;
			kr = 24*esp*(2*pow(sr,12)-pow(sr,6))/(r*r);
			Fx[i]+=kr*rx;
			Fy[i]+=kr*ry;
			Fx[j]-=kr*rx;
			Fy[j]-=kr*ry;
			}
		}
	}
}

void Move(double *x, double *vx, double *x2, double *x3, double *x4, double *x5, double *y, double *vy, double *y2, double *y3, double *y4, double *y5, double *Fx, double *Fy){

    double corx[NP], cory[NP]; //correctors in x & y direction
    int i=0;

    for(i=0;i<NP;i++){
        x2[i]=Fx[i]/m;
        y2[i]=Fy[i]/m;
    }
    //predictor
    for(i=0;i<NP;i++){
        x[i] = x[i] + vx[i]*dt + x2[i]*dt*dt/2.0 + x3[i]*dt*dt*dt/6.0 + x4[i]*dt*dt*dt*dt/24.0 + x5[i]*dt*dt*dt*dt*dt/120.0;
        vx[i] = vx[i] + x2[i]*dt + x3[i]*dt*dt/2.0 + x4[i]*dt*dt*dt/6.0 + x5[i]*dt*dt*dt*dt/24.0;
        x2[i] = x2[i] + x3[i]*dt + x4[i]*dt*dt/2.0 + x5[i]*dt*dt*dt/6.0;
        x3[i] = x3[i] + x4[i]*dt + x5[i]*dt*dt/2.0;
        x4[i] = x4[i] + x5[i]*dt;
        y[i] = y[i] + vy[i]*dt + y2[i]*dt*dt/2.0 + y3[i]*dt*dt*dt/6.0 + y4[i]*dt*dt*dt*dt/24.0 + y5[i]*dt*dt*dt*dt*dt/120.0;
        vy[i] = vy[i] + y2[i]*dt + y3[i]*dt*dt/2.0 + y4[i]*dt*dt*dt/6.0 + y5[i]*dt*dt*dt*dt/24.0;
        y2[i] = y2[i] + y3[i]*dt + y4[i]*dt*dt/2.0 + y5[i]*dt*dt*dt/6.0;
        y3[i] = y3[i] + y4[i]*dt + y5[i]*dt*dt/2.0;
        y4[i] = y4[i] + y5[i]*dt;
    }
    //corrector
    CalcForce(Fx, Fy, x, y);
    for(i=0;i<NP;i++){
        corx[i] = (Fx[i]/m - x2[i])*dt*dt/2.0;
        cory[i] = (Fy[i]/m - y2[i])*dt*dt/2.0;
    }

    //renew
    for(i=0;i<NP;i++){
        x[i] += a0*corx[i];
        vx[i] += a1*corx[i]/dt;
        x2[i] += a2*corx[i]*2.0/dt/dt;
        x3[i] += a3*corx[i]*6.0/dt/dt/dt;
        x4[i] += a4*corx[i]*24.0/dt/dt/dt/dt;
        x5[i] += a5*corx[i]*120.0/dt/dt/dt/dt/dt;
        y[i] += a0*cory[i];
        vy[i] += a1*cory[i]/dt;
        y2[i] += a2*cory[i]*2.0/dt/dt;
        y3[i] += a3*cory[i]*6.0/dt/dt/dt;
        y4[i] += a4*cory[i]*24.0/dt/dt/dt/dt;
        y5[i] += a5*cory[i]*120.0/dt/dt/dt/dt/dt;
    }

    //Periodic boundary conditions
    for(i=0;i<NP;i++){
        while(x[i]>0.5*L){x[i]-=L;}
        while(x[i]<-0.5*L){x[i]+=L;}
        while(y[i]>0.5*L){y[i]-=L;}
        while(y[i]<-0.5*L){y[i]+=L;}
    }

}

void ExportResults(double T[NP], double U[NP], double P[NP], double H[NP]){

	FILE *pFile;
    int i=0;

	//export temperature
	pFile=fopen("Temperature.txt","w");
	for(i=0;i<NT+1;i++){
			fprintf(pFile,"%f\n", T[i]);
	}
	fclose(pFile);
	//export potential energy
	pFile=fopen("Potential.txt","w");
	for(i=0;i<NT+1;i++){
			fprintf(pFile,"%lE\n", U[i]);
	}
	fclose(pFile);
	//export momentum
	pFile=fopen("Momentum.txt","w");
	for(i=0;i<NT+1;i++){
			fprintf(pFile,"%lE\n", P[i]);
	}
	fclose(pFile);
	//export angular momentum
	pFile=fopen("AngMomentum.txt","w");
	for(i=0;i<NT+1;i++){
			fprintf(pFile,"%lE\n", H[i]);
	}

}

int main(){

	int i=0, k=0;
	double x[NP], vx[NP], x2[NP], x3[NP], x4[NP], x5[NP], corx[NP];
	double y[NP], vy[NP], y2[NP], y3[NP], y4[NP], y5[NP], cory[NP];
	double Fx[NP], Fy[NP];  //force on each atom
	double T[NT+1]; //Temperature
	double U[NT+1]; //potential energy
	double P[NT+1], H[NT+1];  //momentum and angular momentum

	Load(x, y, vx, vy);

	for(i=0;i<NP;i++){
		x3[i]=0.0;
		y3[i]=0.0;
		x4[i]=0.0;
		y4[i]=0.0;
		x5[i]=0.0;
		y5[i]=0.0;
	}

	for(k=0;k<NT+1;k++){
		printf("k=%d\n",k);
		CalcTemp(k, T, vx, vy);
		CalcMomentum(k, P, H, x, y, vx, vy);
		CalcEnergy(k, U, x, y);
		CalcForce(Fx, Fy, x, y);
        Move(x, vx, x2, x3, x4, x5, y, vy, y2, y3, y4, y5, Fx, Fy);
	}

	ExportResults(T, U, P, H);

	printf("Done!\n");
	return 0;
}
