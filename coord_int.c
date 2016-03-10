//Initialize positions of th atoms
//The atoms are BCC in this program
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 10 //number of units per side
#define NY 10 //number of units per side

void Calc();


int main(){

	double coord[2*NX*NY][2]; //coordinates of each atom
	int i=0, j=0, m=0;
	int NP=2*NX*NY;
	double L=40.0;//1E-10(m)
	const double a=4.0; //lattice const

	//first two atoms
	coord[0][0]=0.0;
	coord[0][1]=0.0;
	coord[1][0]=0.5*a;
	coord[1][1]=0.5*a;

	//generate coords of other atoms
	for(j=0;j<NY;j++){
		for(i=0;i<NX;i++){
			for(m=0;m<2;m++){
				coord[m+2*(i+NX*j)][0]=coord[m][0]+a*(double)i;
				coord[m+2*(i+NX*j)][1]=coord[m][1]+a*(double)j;
			}
		}
	}

	//center box at (0,0,0)
	for(i=0;i<NP;i++){
		coord[i][0]-=0.5*L;
		coord[i][1]-=0.5*L;
	}

	//export coords of each atoms
	FILE *pFile;
	pFile=fopen("coords.txt","w");
	for(i=0;i<NP;i++){
			fprintf(pFile,"%lE\t%lE\n", coord[i][0], coord[i][1]);
	}
	fclose(pFile);

	printf("Done!\n");
	return 0;

}
