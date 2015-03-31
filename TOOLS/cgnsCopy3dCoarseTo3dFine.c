#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"
#include "unistd.h"

#define CG(cmd) if(cmd)cg_error_exit( );

/*void copycoarse2fine(cgsize_t irmax_fine[3],int *translationmatrix,float *buffer_fine,float *buffer_coarse)
{
	int i,j,k,ijk;
	for(i=0;i<irmax_fine[0];i++)
	{
		for(j=0;j<irmax_fine[1];j++)
		{
			for(k=0;k<irmax_fine[2];k++)
			{
				ijk=k*irmax_fine[1]*irmax_fine[0]+j*irmax_fine[0]+i;
				buffer_fine[ijk]=buffer_coarse[translationmatrix[ijk]];
			}
		}
	}
}*/

void copycoarse2fine(int k_offset,int k_faktor,cgsize_t irmax_fine[3], float *buffer_fine,float *buffer_coarse,
int *topLeftIJK,int *topRightIJK,int *bottomLeftIJK,int *bottomRightIJK,
float *topLeftWeighting,float *topRightWeighting,float *bottomLeftWeighting,float *bottomRightWeighting
)
{
	int i,j,k,ijk,k2;
	int ijk0;
	for(i=0;i<irmax_fine[0];i++)
	{
		for(j=0;j<irmax_fine[1];j++)
		{
			ijk0=0*irmax_fine[1]*irmax_fine[0]+j*irmax_fine[0]+i;
			for(k=0;k<irmax_fine[2];k++)
			{
				ijk=k*irmax_fine[1]*irmax_fine[0]+j*irmax_fine[0]+i;
				k2=(int)k/k_faktor*k_offset;
				buffer_fine[ijk]=
				(buffer_coarse[topLeftIJK[ijk0]+k2]*topLeftWeighting[ijk0]+
				buffer_coarse[topRightIJK[ijk0]+k2]*topRightWeighting[ijk0]+
				buffer_coarse[bottomLeftIJK[ijk0]+k2]*bottomLeftWeighting[ijk0]+
				buffer_coarse[bottomRightIJK[ijk0]+k2]*bottomRightWeighting[ijk0])/
				(topLeftWeighting[ijk0]+topRightWeighting[ijk0]+bottomLeftWeighting[ijk0]+bottomRightWeighting[ijk0]);
			}
		}
	}
}


int main(int argc, char *argv[])
{
        char *ovalue = NULL;
        char *ivalue = NULL;
	char file_coarse[500];
	char file_fine[500];	
        int c;

        opterr = 0;

        while ((c = getopt (argc, argv, "f:c:")) != -1)
        {
                switch (c)
                {
                        case 'c':
				strcpy(file_coarse,optarg);
                                break;
                        case 'f':
				strcpy(file_fine,optarg);
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



	printf("Daten werden aus grober Gitter-Datei: %s in hochaufgelöstes Gitter: %s kopiert\n",file_coarse,file_fine);


//	Öffnen beider Dateien
	int index_file_coarse,index_file_fine;
	CG( cg_open(file_coarse,CG_MODE_READ,&index_file_coarse));
	CG( cg_open(file_fine,CG_MODE_MODIFY,&index_file_fine));

	int index_base=1;
	int index_zone=1;
	int index_sol=1;
	

	
	char zonename[33];
	cgsize_t zonesize_coarse[3][3];
	cgsize_t zonesize_fine[3][3];	
	printf("Reading zone of coarse...\n");
	CG( cg_zone_read(index_file_coarse,index_base,index_zone,zonename,zonesize_coarse[0]));
	printf("Reading zone of fine...\n");
	CG( cg_zone_read(index_file_fine,index_base,index_zone,zonename,zonesize_fine[0]));
	
	cgsize_t irmin_fine[3];
	cgsize_t irmax_fine[3];	
	cgsize_t irmin_coarse[3];
	cgsize_t irmax_coarse[3];
	
	irmin_coarse[0]=1;
	irmin_coarse[1]=1;
	irmin_coarse[2]=1;
	irmax_coarse[0]=zonesize_coarse[0][0];
	irmax_coarse[1]=zonesize_coarse[0][1];
	irmax_coarse[2]=zonesize_coarse[0][2];		

	irmin_fine[0]=1;
	irmin_fine[1]=1;
	irmin_fine[2]=1;
	irmax_fine[0]=zonesize_fine[0][0];
	irmax_fine[1]=zonesize_fine[0][1];
	irmax_fine[2]=zonesize_fine[0][2];	
	
	printf("coarse grid: %zux%zux%zu\n",irmax_coarse[0],irmax_coarse[1],irmax_coarse[2]);
	printf("fine grid: %zux%zux%zu\n",irmax_fine[0],irmax_fine[1],irmax_fine[2]);
	
	float *CoordinateX_fine;
	float *CoordinateY_fine;
	float *CoordinateX_coarse;
	float *CoordinateY_coarse;	


	CoordinateX_fine=(float*)calloc(irmax_fine[0]*irmax_fine[1]*irmax_fine[2],sizeof(float));
	CoordinateY_fine=(float*)calloc(irmax_fine[0]*irmax_fine[1]*irmax_fine[2],sizeof(float));
	CoordinateX_coarse=(float*)calloc(irmax_coarse[0]*irmax_coarse[1]*irmax_coarse[2],sizeof(float));
	CoordinateY_coarse=(float*)calloc(irmax_coarse[0]*irmax_coarse[1]*irmax_coarse[2],sizeof(float));	

	printf("Creating TranslationMatrix\n");
	printf("Reading CoordinateX&Y of coarse...\n");
	CG( cg_coord_read(index_file_coarse,index_base,index_zone,"CoordinateX",RealSingle,irmin_coarse,irmax_coarse,CoordinateX_coarse));
	CG( cg_coord_read(index_file_coarse,index_base,index_zone,"CoordinateY",RealSingle,irmin_coarse,irmax_coarse,CoordinateY_coarse));
	printf("Reading CoordinateX&Y of fine...\n");
	CG( cg_coord_read(index_file_fine,index_base,index_zone,"CoordinateX",RealSingle,irmin_fine,irmax_fine,CoordinateX_fine));
	CG( cg_coord_read(index_file_fine,index_base,index_zone,"CoordinateY",RealSingle,irmin_fine,irmax_fine,CoordinateY_fine));

	int i_coarse,j_coarse,k,ijk_coarse;
	int i_fine,j_fine,ijk_fine;
	//int *translationmatrix;
	//translationmatrix=(int*)calloc(irmax_fine[0]*irmax_fine[1]*irmax_fine[2],sizeof(int));
	int *topLeftIJK=(int*)calloc(irmax_fine[0]*irmax_fine[1],sizeof(int));
	int *topRightIJK=(int*)calloc(irmax_fine[0]*irmax_fine[1],sizeof(int));
	int *bottomLeftIJK=(int*)calloc(irmax_fine[0]*irmax_fine[1],sizeof(int));
	int *bottomRightIJK=(int*)calloc(irmax_fine[0]*irmax_fine[1],sizeof(int));
	
	float *topLeftWeighting=(float*)calloc(irmax_fine[0]*irmax_fine[1],sizeof(float));
	float *topRightWeighting=(float*)calloc(irmax_fine[0]*irmax_fine[1],sizeof(float));
	float *bottomLeftWeighting=(float*)calloc(irmax_fine[0]*irmax_fine[1],sizeof(float));
	float *bottomRightWeighting=(float*)calloc(irmax_fine[0]*irmax_fine[1],sizeof(float));	
	
	float distance;
	float distance_smallest;
	printf("Searching for corresponding meshpoints within k=0 plane...\n");
	k=0;

	int i_faktor=(int)zonesize_fine[0][0]/(int)zonesize_coarse[0][0];
	int j_faktor=(int)zonesize_fine[0][1]/(int)zonesize_coarse[0][1];
	int k_faktor=(int)zonesize_fine[0][2]/(int)zonesize_coarse[0][2];

	int i_start_coarse,i_end_coarse;		
	int j_start_coarse,j_end_coarse;

	int i_found;
	int j_found;	
	
	int i,i_situation;
	int j,j_situation;
	int ijk;
	

	j_start_coarse=1;
	j_end_coarse=30;
	for(j_fine=0;j_fine<irmax_fine[1];j_fine++)
	{
		i_start_coarse=1;
		i_end_coarse=30;
		for(i_fine=0;i_fine<irmax_fine[0];i_fine++)
		{
			ijk_fine=k*zonesize_fine[0][1]*zonesize_fine[0][0]+j_fine*zonesize_fine[0][0]+i_fine;
			distance_smallest=9999999999.;
			
			if(i_start_coarse<1){i_start_coarse=1;}
			if(j_start_coarse<1){j_start_coarse=1;}
			if(i_end_coarse>=irmax_coarse[0]-1){i_end_coarse=irmax_coarse[0]-1;}
			if(j_end_coarse>=irmax_coarse[1]-1){j_end_coarse=irmax_coarse[1]-1;}

			for(i_coarse=i_start_coarse;i_coarse<i_end_coarse;i_coarse++)
			{
				for(j_coarse=j_start_coarse;j_coarse<j_end_coarse;j_coarse++)
				{
					ijk_coarse=k*zonesize_coarse[0][1]*zonesize_coarse[0][0]+j_coarse*zonesize_coarse[0][0]+i_coarse;
				
					distance=sqrt(
					(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])*(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])
					+(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])*(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])
					);
					if(distance_smallest>distance)
					{
						i_found=i_coarse;
						j_found=j_coarse;
						
						distance_smallest=distance;
					}
				}
			}
						//if((i_fine==600)&&(j_fine<50)){printf("ifound:%d jfound:%d \n",i_found,j_found);printf("i:%d-%d | j:%d-%d     --- distance:%f\n",i_start_coarse,i_end_coarse,j_start_coarse,j_end_coarse,distance_smallest);}
			i_start_coarse=i_found-15;
			i_end_coarse=i_found+15;
			j_start_coarse=j_found-15;
			j_end_coarse=j_found+15;			
			
			//printf("ifound:%d jfound:%d \n",i_found,j_found);

			distance_smallest=999999.;
			for(i=-1;i<2;i=i+2){for(j=-1;j<2;j=j+2){
				i_coarse=i_found+i;	j_coarse=j_found+j;
				ijk_coarse=k*zonesize_coarse[0][1]*zonesize_coarse[0][0]+j_coarse*zonesize_coarse[0][0]+i_coarse;			
			
				distance=sqrt(
				(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])*(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])
				+(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])*(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])
				);
				if(distance<distance_smallest)
				{
					distance_smallest=distance;
					if((i==1)&&(j==1)){i_situation=1;j_situation=1;} 	//case1: i,j_situation enthalten die offsets von topright bezueglich ij_nearest
					if((i==-1)&&(j==1)){i_situation=0;j_situation=1;}	//case2
					if((i==-1)&&(j==-1)){i_situation=0;j_situation=0;}	//case3
					if((i==1)&&(j==-1)){i_situation=1;j_situation=0;}	//case4
				}
			}}

			
			// topRight
			ijk_coarse=k*zonesize_coarse[0][1]*zonesize_coarse[0][0]+(j_found+j_situation)*zonesize_coarse[0][0]+(i_found+i_situation);
			distance=sqrt(
			(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])*(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])
			+(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])*(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])
			);			
			topRightIJK[ijk_fine]=ijk_coarse;	topRightWeighting[ijk_fine]=1.0/(distance+1.e-20);
			
			// topLeft
			ijk_coarse=k*zonesize_coarse[0][1]*zonesize_coarse[0][0]+(j_found+j_situation)*zonesize_coarse[0][0]+(i_found+i_situation-1);
			distance=sqrt(
			(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])*(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])
			+(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])*(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])
			);			
			topLeftIJK[ijk_fine]=ijk_coarse;	topLeftWeighting[ijk_fine]=1.0/(distance+1.e-20);
			
			// bottomRight
			ijk_coarse=k*zonesize_coarse[0][1]*zonesize_coarse[0][0]+(j_found+j_situation-1)*zonesize_coarse[0][0]+(i_found+i_situation);
			distance=sqrt(
			(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])*(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])
			+(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])*(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])
			);			
			bottomRightIJK[ijk_fine]=ijk_coarse;	bottomRightWeighting[ijk_fine]=1.0/(distance+1.e-20);
			
			// bottomLeft
			ijk_coarse=k*zonesize_coarse[0][1]*zonesize_coarse[0][0]+(j_found+j_situation-1)*zonesize_coarse[0][0]+(i_found+i_situation-1);
			distance=sqrt(
			(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])*(CoordinateX_coarse[ijk_coarse]-CoordinateX_fine[ijk_fine])
			+(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])*(CoordinateY_coarse[ijk_coarse]-CoordinateY_fine[ijk_fine])
			);			
			bottomLeftIJK[ijk_fine]=ijk_coarse;	bottomLeftWeighting[ijk_fine]=1.0/(distance+1.e-20);
					
					
			
		}
	}
	
	printf("Free memory of CoordinateX&Y-arrays\n");
	free(CoordinateX_fine);
	free(CoordinateY_fine);
	free(CoordinateX_coarse);
	free(CoordinateY_coarse);
	
	int ijk0;
	int k_fine;

	char iterationstr[ 255 ];
	int int_actualIteration=0;
	snprintf( iterationstr,255,"%i",int_actualIteration );
	CG( cg_goto( index_file_fine,index_base,"end" ) );
	CG( cg_descriptor_write( "Iterations",iterationstr ) );

	int samples=1;
	float* timearray;
	timearray = malloc( ( samples )*sizeof( float ) );
	timearray[0]=0.0;
	CG( cg_biter_write( index_file_fine,index_base,"BaseIterativeData",samples ) );
	CG( cg_goto( index_file_fine,index_base,"BaseIterativeData",0,"end" ) );
	CG( cg_array_write( "TimeValues",RealSingle,1,( cgsize_t[ ] ){ samples },timearray ) );

	char solname[ 33 ];
	snprintf( solname,33,"Sol_%04d",0 );
	CG( cg_ziter_write( index_file_fine,index_base,1,"ZoneIterativeData" ) );
	CG( cg_goto( index_file_fine,index_base,"Zone_t",1,"ZoneIterativeData_t",1,"end" ) );
	CG( cg_array_write( "FlowSolutionPointers",Character,2,(cgsize_t[ ]){ 32,samples },solname ) );	
	printf("Allocating memory for buffer-array\n");
	float *buffer_fine;
	float *buffer_coarse;
	buffer_fine=(float*)calloc(irmax_fine[0]*irmax_fine[1]*irmax_fine[2],sizeof(float));
	buffer_coarse=(float*)calloc(irmax_coarse[0]*irmax_coarse[1]*irmax_coarse[2],sizeof(float));
	
	int n;
	CG( cg_sol_write( index_file_fine,index_base,index_zone,solname,Vertex,&n ) );
	int field;
	int nfields;
	int index_field;
	char fieldname[ 33 ];
	int k_offset=zonesize_coarse[0][1]*zonesize_coarse[0][0];
	DataType_t datatype;
	CG( cg_nfields( index_file_coarse,index_base,index_zone,index_sol,&nfields ) );
	printf("Es gibt %d Variablen.\n",nfields);
	nfields=5;
	printf("Es werden %d Variablen übertragen.\n",nfields);
	for( field = 1; field<=nfields; field++ ) 
	{
		cg_field_info(index_file_coarse,index_base,index_zone,index_sol, field, &datatype, fieldname);
		printf("%s...\n",fieldname);
		cg_field_read(index_file_coarse,index_base,index_zone,index_sol,fieldname,RealSingle,irmin_coarse,irmax_coarse,buffer_coarse);
		copycoarse2fine(k_offset,k_faktor,irmax_fine,buffer_fine,buffer_coarse,
			topLeftIJK,topRightIJK,bottomLeftIJK,bottomRightIJK,
			topLeftWeighting,topRightWeighting,bottomLeftWeighting,bottomRightWeighting);
		cg_field_write(index_file_fine,index_base,index_zone,index_sol,RealSingle,fieldname,buffer_fine,&index_field);
	}
	
	
	printf("Closing files...\n");
	CG( cg_close(index_file_coarse));
	CG( cg_close(index_file_fine));
}
