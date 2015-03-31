#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"
#include "unistd.h"

#define CG(cmd) if(cmd)cg_error_exit( );

void copyFine2Coarse(cgsize_t irmax_coarse[2],int *translationmatrix,float *buffer_coarse,float *buffer_fine)
{
	int i,j,k,ijk;
	for(i=0;i<irmax_coarse[0];i++)
	{
		for(j=0;j<irmax_coarse[1];j++)
		{
			ijk=j*irmax_coarse[0]+i;
			buffer_coarse[ijk]=buffer_fine[translationmatrix[ijk]];
		}
	}
}


int main(int argc, char *argv[])
{
	char file_fine[500];
	char file_coarse[500];	
	
	char *ivalue = NULL;
	int inputflag;
	int c;

	opterr = 0;

	while ((c = getopt (argc, argv, "i:")) != -1)
	{
		switch (c)
		{
			case 'i':
				inputflag=1;
				ivalue = optarg;
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
		printf("ERROR: Keine Datei für das grobe Gitter angegeben (-i 'filename')\n");
		return 1;
	}
	else
	{
		strcpy(file_coarse,ivalue);
	}	
	

	sprintf(file_fine,"/home/mag/WORK2/NACA3D_M065a4/2d_NACA0012_yMin40_500000_8192x512.cgns");
	//sprintf(file_coarse,"/home/mag/DISSERTATION/3D-GKS/512x128x64/3d_NACA0012_yMin40_500000_512x128x64.cgns");
	printf("Daten werden aus hochaufgelöster Gitter-Datei: %s in gröbere aufgelösteres Gitter: %s kopiert\n",file_fine,file_coarse);


//	Öffnen beider Dateien
	int index_file_fine,index_file_coarse;
	CG( cg_open(file_fine,CG_MODE_READ,&index_file_fine));
	CG( cg_open(file_coarse,CG_MODE_MODIFY,&index_file_coarse));

	int index_base=1;
	int index_zone=1;
	int index_sol=1;
	

	
	char zonename[33];
	cgsize_t zonesize_fine[2][2];
	cgsize_t zonesize_coarse[2][2];	
	printf("Reading zone of fine...\n");
	CG( cg_zone_read(index_file_fine,index_base,index_zone,zonename,zonesize_fine[0]));
	printf("Reading zone of coarse...\n");
	CG( cg_zone_read(index_file_coarse,index_base,index_zone,zonename,zonesize_coarse[0]));
	
	cgsize_t irmin_coarse[2];
	cgsize_t irmax_coarse[2];	
	cgsize_t irmin_fine[2];
	cgsize_t irmax_fine[2];
	
	irmin_fine[0]=1;
	irmin_fine[1]=1;
	irmax_fine[0]=zonesize_fine[0][0];
	irmax_fine[1]=zonesize_fine[0][1];

	irmin_coarse[0]=1;
	irmin_coarse[1]=1;
	irmax_coarse[0]=zonesize_coarse[0][0];
	irmax_coarse[1]=zonesize_coarse[0][1];
	
	printf("Fine grid: %zux%zu\n",irmax_fine[0],irmax_fine[1]);
	printf("Coarse grid: %zux%zu\n",irmax_coarse[0],irmax_coarse[1]);
	
	float *CoordinateX_coarse;
	float *CoordinateY_coarse;
	float *CoordinateX_fine;
	float *CoordinateY_fine;	


	CoordinateX_coarse=(float*)calloc(irmax_coarse[0]*irmax_coarse[1],sizeof(float));
	CoordinateY_coarse=(float*)calloc(irmax_coarse[0]*irmax_coarse[1],sizeof(float));
	CoordinateX_fine=(float*)calloc(irmax_fine[0]*irmax_fine[1],sizeof(float));
	CoordinateY_fine=(float*)calloc(irmax_fine[0]*irmax_fine[1],sizeof(float));	

	printf("Creating TranslationMatrix\n");
	printf("Reading CoordinateX&Y of fine...\n");
	CG( cg_coord_read(index_file_fine,index_base,index_zone,"CoordinateX",RealSingle,irmin_fine,irmax_fine,CoordinateX_fine));
	CG( cg_coord_read(index_file_fine,index_base,index_zone,"CoordinateY",RealSingle,irmin_fine,irmax_fine,CoordinateY_fine));
	printf("Reading CoordinateX&Y of coarse...\n");
	CG( cg_coord_read(index_file_coarse,index_base,index_zone,"CoordinateX",RealSingle,irmin_coarse,irmax_coarse,CoordinateX_coarse));
	CG( cg_coord_read(index_file_coarse,index_base,index_zone,"CoordinateY",RealSingle,irmin_coarse,irmax_coarse,CoordinateY_coarse));

	int i_fine,j_fine,ijk_fine;
	int i_coarse,j_coarse,ijk_coarse;
	int ijk_correspondig;
	int *translationmatrix;
	translationmatrix=(int*)calloc(irmax_coarse[0]*irmax_coarse[1],sizeof(int));
	float distance;
	float distance_smallest;
	printf("Searching for corresponding meshpoints...\n");

	int i_faktor=zonesize_fine[0][0]/zonesize_coarse[0][0];
	int j_faktor=zonesize_fine[0][1]/zonesize_coarse[0][1];
	
	int i_start_fine,i_end_fine;		
	int j_start_fine,j_end_fine;
	
	int searchWindow=100;
	printf("The search-window is set to %d\n",searchWindow);
	for(i_coarse=0;i_coarse<irmax_coarse[0];i_coarse++)
	{
		for(j_coarse=0;j_coarse<irmax_coarse[1];j_coarse++)
		{
			ijk_coarse=j_coarse*zonesize_coarse[0][0]+i_coarse;
			distance_smallest=9999999999.;
			
			
			i_start_fine=i_coarse*i_faktor-searchWindow;
			i_end_fine=i_coarse*i_faktor+searchWindow;
			j_start_fine=j_coarse*j_faktor-searchWindow;
			j_end_fine=j_coarse*j_faktor+searchWindow;
			
			if(i_start_fine<0){i_start_fine=0;}
			if(j_start_fine<0){j_start_fine=0;}
			if(i_end_fine>=irmax_fine[0]){i_end_fine=irmax_fine[0];}
			if(j_end_fine>=irmax_fine[1]){j_end_fine=irmax_fine[1];}
			
			for(i_fine=i_start_fine;i_fine<i_end_fine;i_fine++)
			{
				for(j_fine=j_start_fine;j_fine<j_end_fine;j_fine++)
				{
					ijk_fine=j_fine*zonesize_fine[0][0]+i_fine;
				
					distance=sqrt(
					(CoordinateX_fine[ijk_fine]-CoordinateX_coarse[ijk_coarse])*(CoordinateX_fine[ijk_fine]-CoordinateX_coarse[ijk_coarse])
					+(CoordinateY_fine[ijk_fine]-CoordinateY_coarse[ijk_coarse])*(CoordinateY_fine[ijk_fine]-CoordinateY_coarse[ijk_coarse])
					);
					if(distance_smallest>distance)
					{
						ijk_correspondig=ijk_fine;
						distance_smallest=distance;
					}
				}
			}
			//printf("distance_smallest: %g\n",distance_smallest);
			translationmatrix[ijk_coarse]=ijk_correspondig;
		}
	}
	
	printf("Free memory of CoordinateX&Y-arrays\n");
	free(CoordinateX_coarse);
	free(CoordinateY_coarse);
	free(CoordinateX_fine);
	free(CoordinateY_fine);
	
	
	/*
	int ijk0;
	int k_coarse;
	printf("Determing other k-Planes...\n");
	for(i_coarse=0;i_coarse<irmax_coarse[0];i_coarse++)
	{
		for(j_coarse=0;j_coarse<irmax_coarse[1];j_coarse++)
		{
			ijk0=0*zonesize_coarse[0][1]*zonesize_coarse[0][0]+j_coarse*zonesize_coarse[0][0]+i_coarse;
			for(k_coarse=1;k_coarse<irmax_coarse[2];k_coarse++)
			{
				ijk_coarse=k_coarse*zonesize_coarse[0][1]*zonesize_coarse[0][0]+j_coarse*zonesize_coarse[0][0]+i_coarse;
				translationmatrix[ijk_coarse]=translationmatrix[ijk0]+(k_faktor*k_coarse)*zonesize_fine[0][1]*zonesize_fine[0][0];
			}
		}
	}
	*/

	char iterationstr[ 255 ];
	int int_actualIteration=0;
	snprintf( iterationstr,255,"%i",int_actualIteration );
	CG( cg_goto( index_file_coarse,index_base,"end" ) );
	CG( cg_descriptor_write( "Iterations",iterationstr ) );
	printf("Writing BasiterativeData...\n");
	int samples=1;
	float* timearray;
	timearray = malloc( ( samples )*sizeof( float ) );
	timearray[0]=0.0;
	CG( cg_biter_write( index_file_coarse,index_base,"BaseIterativeData",samples ) );
	CG( cg_goto( index_file_coarse,index_base,"BaseIterativeData",0,"end" ) );
	CG( cg_array_write( "TimeValues",RealSingle,1,( cgsize_t[ ] ){ samples },timearray ) );

	printf("Writing ZoneIterativeData...\n");
	char solname[ 33 ];
	snprintf( solname,33,"Sol_%04d",0 );
	CG( cg_ziter_write( index_file_coarse,index_base,1,"ZoneIterativeData" ) );
	CG( cg_goto( index_file_coarse,index_base,"Zone_t",1,"ZoneIterativeData_t",1,"end" ) );
	CG( cg_array_write( "FlowSolutionPointers",Character,2,(cgsize_t[ ]){ 32,samples },solname ) );
	printf("Allocating memory for buffer-array\n");
	float *buffer_coarse;
	float *buffer_fine;
	buffer_coarse=(float*)calloc(irmax_coarse[0]*irmax_coarse[1],sizeof(float));
	buffer_fine=(float*)calloc(irmax_fine[0]*irmax_fine[1],sizeof(float));

	int n;
	CG( cg_sol_write( index_file_coarse,index_base,index_zone,solname,Vertex,&n ) );
	int field;
	int nfields;
	int index_field;
	char fieldname[ 33 ];
	DataType_t datatype;
	CG( cg_nfields( index_file_fine,index_base,index_zone,index_sol,&nfields ) );
	printf("Es gibt %d Variablen.\n",nfields);
	nfields=5;
	printf("Es werden %d Variablen übertragen.\n",nfields);
	for( field = 1; field<=nfields; field++ ) 
	{
		cg_field_info(index_file_fine,index_base,index_zone,index_sol, field, &datatype, fieldname);
		printf("%s...\n",fieldname);
		cg_field_read(index_file_fine,index_base,index_zone,index_sol,fieldname,RealSingle,irmin_fine,irmax_fine,buffer_fine);
		copyFine2Coarse(irmax_coarse,translationmatrix,buffer_coarse,buffer_fine);
		cg_field_write(index_file_coarse,index_base,index_zone,index_sol,RealSingle,fieldname,buffer_coarse,&index_field);
	}
	
	
	printf("Closing files...\n");
	CG( cg_close(index_file_fine));
	CG( cg_close(index_file_coarse));
}
