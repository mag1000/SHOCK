#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"
#include "unistd.h"
#include "limits.h"
#include <inttypes.h>
#include <stdint.h>

#define CG(cmd) if(cmd)cg_error_exit( );


int main(int argc, char *argv[])
{
	char *inputFile = NULL;
	char *outputFile = NULL;
	int inputflag;
	int outputflag;	
	int c;

	opterr = 0;

	while ((c = getopt (argc, argv, "i:o:")) != -1)
	{
		switch (c)
		{
			case 'i':
				inputflag=1;
				inputFile = optarg;
				break;
			case 'o':
				outputflag=1;
				outputFile = optarg;
				break;				
			case '?':
				if (optopt == 'o')
					fprintf (stderr, "Option -%o requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%o'.\n", optopt);
				else
					fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
			return 1;
		}
	}
	
	if(inputflag!=1)
	{
		printf("ERROR: Kein Inputfile angegeben (-i 'filename')\n");
		return 1;
	}
	else if(outputflag!=1)
	{
		printf("ERROR: Kein Outputfile angegeben (-o 'filename')\n");
		return 1;
	}
	
	int index_file_in;
	printf("Opening source: %s\n",inputFile);
	CG( cg_open(inputFile,CG_MODE_READ,&index_file_in));


	
	int index_file_out;
	printf("Opening destination: %s\n",outputFile);
	CG( cg_open(outputFile,CG_MODE_MODIFY,&index_file_out));	
	
	
	char zonename[33];
	cgsize_t zonesize[3][3];

	printf("Reading zone...\n");
	CG( cg_zone_read(index_file_out,1,1,zonename,zonesize[0]));
	
	printf("Generating mesh for zone %zux%zux%zu\n",zonesize[0][0],zonesize[0][1],zonesize[0][2]);
	
	cgsize_t irmin_in[3];
	cgsize_t irmax_in[3];
	
	irmin_in[0]=1;
	irmin_in[1]=1;
	irmin_in[2]=1;
	irmax_in[0]=zonesize[0][0];
	irmax_in[1]=zonesize[0][1];
	irmax_in[2]=1;		

	float *buffer_out;
	float *buffer_in;
	long long int iMeshPoints,jMeshPoints,kMeshPoints;
	iMeshPoints=zonesize[0][0];
	jMeshPoints=zonesize[0][1];
	kMeshPoints=zonesize[0][2];
	long long int size;
	size=iMeshPoints*jMeshPoints*kMeshPoints;
	printf("size: %lldx%lldx%lld = %lld (INT_MAX:%d LONG_MAX:%lu)\n",iMeshPoints,jMeshPoints,kMeshPoints,size,INT_MAX,LONG_MAX);

	printf("Allocating memory...\n");
	buffer_in=(float*)calloc(iMeshPoints*jMeshPoints,sizeof(float));
	buffer_out=(float*)calloc(size,sizeof(float));

	int index_coord;
	long long i,j,k,ijk_out,ijk_in;	
	////////////////////////////////COORDINATE X ////////////////////////////////////////////
	printf("Transfering CoordinateX...\n");
	cg_coord_read(index_file_in,1,1,"CoordinateX",RealSingle,irmin_in,irmax_in,buffer_in);
	for(i=0;i<iMeshPoints;i++)
	{
		printf("%d/100 \r",(int)((i*100)/(iMeshPoints-1)));fflush(stdout);
		for(j=0;j<jMeshPoints;j++)
		{
			for(k=0;k<kMeshPoints;k++)
			{
				ijk_in=0*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ijk_out=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				buffer_out[ijk_out]=buffer_in[ijk_in];
			}
		}
	}
	printf("\n");
	printf("Writing CoordinateX...\n");
	cg_coord_write(index_file_out,1,1,RealSingle,"CoordinateX",buffer_out,&index_coord);
	
	////////////////////////////////COORDINATE Y ////////////////////////////////////////////
	printf("Transfering CoordinateY...\n");
	cg_coord_read(index_file_in,1,1,"CoordinateY",RealSingle,irmin_in,irmax_in,buffer_in);
	for(i=0;i<iMeshPoints;i++)
	{
		printf("%d/100 \r",(int)((i*100)/iMeshPoints));fflush(stdout);
		for(j=0;j<jMeshPoints;j++)
		{
			for(k=0;k<kMeshPoints;k++)
			{
				ijk_in=0*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ijk_out=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				buffer_out[ijk_out]=buffer_in[ijk_in];
			}
		}
	}
	printf("\n");
	printf("Writing CoordinateY...\n");
	cg_coord_write(index_file_out,1,1,RealSingle,"CoordinateY",buffer_out,&index_coord);	
	
	free(buffer_in);
	
		
	////////////////////////////////COORDINATE Z ////////////////////////////////////////////
	float CoordinateZ;
	for(i=0;i<iMeshPoints;i++)
	{
		printf("%d/100 \r",(int)((i*100)/iMeshPoints));fflush(stdout);
		for(j=0;j<jMeshPoints;j++)
		{
			for(k=0;k<kMeshPoints;k++)
			{
				CoordinateZ=0.1*k/(kMeshPoints-1);
				ijk_in=0*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ijk_out=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				buffer_out[ijk_out]=CoordinateZ;
			}
		}
	}
	printf("\n");
	printf("Writing CoordinateZ...\n");
	cg_coord_write(index_file_out,1,1,RealSingle,"CoordinateZ",buffer_out,&index_coord);

	free(buffer_out);
	
	CG (cg_close(index_file_out));
	CG (cg_close(index_file_in));
}
