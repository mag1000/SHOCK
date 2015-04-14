#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "unistd.h"
#include "FUNCTIONS.h"
#define FS 1.25

int main(int argc, char *argv[])
{
	int MeshPoints[5];
	MeshPoints[0]=8;	MeshPoints[1]=16;	MeshPoints[2]=32;	MeshPoints[3]=64;	MeshPoints[4]=128;
	int mesh;
	float dx;
	int i;
	float f1,f2,f_abl;
	f_abl;
	FILE * file;
	char output[300];
	for(mesh=0;mesh<5;mesh++)
	{
		sprintf(output,"Mesh_%d.dat",mesh);
		file=fopen(output,"w");
		fprintf(file,"variables = \"x\" \"y\"\n");
		fprintf(file,"zone t=\"%s\", i= %d, f=point \n",output,MeshPoints[mesh]);		
		dx=(float)2.0*M_PI/(MeshPoints[mesh]-2);
		for(i=0;i<MeshPoints[mesh];i++)
		{
			if(i==0)
			{
				f1=sin((i)*dx);
				f2=sin((i+1)*dx);
				f_abl=(f2-f1)/(1.0*dx);
			}else if(i==MeshPoints[mesh]-1)
			{
				f1=sin((i-1)*dx);
				f2=sin((i)*dx);
				f_abl=(f2-f1)/(1.0*dx);
			}
			else
			{
				f1=sin((i-1)*dx);
				f2=sin((i+1)*dx);
				f_abl=(f2-f1)/(2.0*dx);
			}

			fprintf(file," %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",(i*dx),f_abl,f_abl,f_abl,f_abl,f_abl,f_abl,f_abl,f_abl,f_abl,f_abl,f_abl,f_abl,f_abl,f_abl);
		}
		fclose(file);
	}

}

