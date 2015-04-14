#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"
#include "unistd.h"
#include <inttypes.h>
#include <stdint.h>

#define CG(cmd) if(cmd)cg_error_exit( );


void copyGrid1ToGrid2(int k_offset,float k_faktor,cgsize_t irmax_out[3], cgsize_t irmax_in[3], float *buffer_out,float *buffer_in,
int *topLeftIJK,int *topRightIJK,int *bottomLeftIJK,int *bottomRightIJK,
float *topLeftWeighting,float *topRightWeighting,float *bottomLeftWeighting,float *bottomRightWeighting
)
{
	long long int i,j,k,ijk,k_in,k_up,k_low;
	int ijk0;
	int kMeshPoints_in=irmax_in[2];
	int kMeshPoints_out=irmax_out[2];
	long long int irmax_out0=irmax_out[0];
	long long int irmax_out1=irmax_out[1];
	long long int irmax_out2=irmax_out[2];
	float delta_Z_in=0.1/(kMeshPoints_in-1);
	float delta_Z_out=0.1/(kMeshPoints_out-1);
	float z_in,z_out;
	float Weighting_low, Weighting_up;
	
	for(k=0;k<irmax_out[2];k++)
	{	
		z_out=delta_Z_out*k;
		z_in=0;
		k_in=0;
		while (z_in<z_out)
		{
			k_in++;
			z_in=delta_Z_in*k_in;			
		}
		Weighting_low=1.0-(z_out-delta_Z_in*(k_in-1))/delta_Z_in;
		Weighting_up=1.0-Weighting_low;
		if(k==0)
		{
			k_in=1;
			Weighting_low=1.0;
			Weighting_up=0.0;
		}
		if(k==irmax_out[2]-1)
		{
			k_in=irmax_in[2]-1;
			Weighting_low=0.0;
			Weighting_up=1.0;
		}
		
		for(i=0;i<irmax_out[0];i++)
		{
			for(j=0;j<irmax_out[1];j++)
			{
				ijk0=j*irmax_out[0]+i;

				ijk=k*irmax_out1*irmax_out0+j*irmax_out0+i;
				k_low=(k_in-1)*k_offset;
				k_up=k_in*k_offset;
				
				buffer_out[ijk]=
				Weighting_low*(
				(buffer_in[topLeftIJK[ijk0]+k_low]*topLeftWeighting[ijk0]+
				buffer_in[topRightIJK[ijk0]+k_low]*topRightWeighting[ijk0]+
				buffer_in[bottomLeftIJK[ijk0]+k_low]*bottomLeftWeighting[ijk0]+
				buffer_in[bottomRightIJK[ijk0]+k_low]*bottomRightWeighting[ijk0]))
				+
				Weighting_up*(
				(buffer_in[topLeftIJK[ijk0]+k_up]*topLeftWeighting[ijk0]+
				buffer_in[topRightIJK[ijk0]+k_up]*topRightWeighting[ijk0]+
				buffer_in[bottomLeftIJK[ijk0]+k_up]*bottomLeftWeighting[ijk0]+
				buffer_in[bottomRightIJK[ijk0]+k_up]*bottomRightWeighting[ijk0]))
				;
				
			}
		}
	}
}


int main(int argc, char *argv[])
{
        char *ovalue = NULL;
        char *ivalue = NULL;
	char file_in[500];
	char file_out[500];	
        int c;

        opterr = 0;

        while ((c = getopt (argc, argv, "i:o:")) != -1)
        {
                switch (c)
                {
                        case 'i':
				strcpy(file_in,optarg);
                                break;
                        case 'o':
				strcpy(file_out,optarg);
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



	printf("Daten werden aus Gitter-Datei: %s in Gitter: %s kopiert\n",file_in,file_out);


//	Öffnen beider Dateien
	int index_file_in,index_file_out;
	CG( cg_open(file_in,CG_MODE_READ,&index_file_in));
	CG( cg_open(file_out,CG_MODE_MODIFY,&index_file_out));

	int index_base=1;
	int index_zone=1;
	int index_sol=1;
	

	
	char zonename[33];
	cgsize_t zonesize_in[3][3];
	cgsize_t zonesize_out[3][3];	
	printf("Reading zone of grid1...\n");
	CG( cg_zone_read(index_file_in,index_base,index_zone,zonename,zonesize_in[0]));
	printf("Reading zone of grid2...\n");
	CG( cg_zone_read(index_file_out,index_base,index_zone,zonename,zonesize_out[0]));
	
	cgsize_t irmin_out[3];
	cgsize_t irmax_out[3];	
	cgsize_t irmin_in[3];
	cgsize_t irmax_in[3];
	
	irmin_in[0]=1;
	irmin_in[1]=1;
	irmin_in[2]=1;
	irmax_in[0]=zonesize_in[0][0];
	irmax_in[1]=zonesize_in[0][1];
	irmax_in[2]=1;		

	irmin_out[0]=1;
	irmin_out[1]=1;
	irmin_out[2]=1;
	irmax_out[0]=zonesize_out[0][0];
	irmax_out[1]=zonesize_out[0][1];
	irmax_out[2]=1;	
	
	printf("grid1 grid: %zux%zux%zu\n",irmax_in[0],irmax_in[1],irmax_in[2]);
	printf("grid2 grid: %zux%zux%zu\n",irmax_out[0],irmax_out[1],irmax_out[2]);
	
	float *CoordinateX_out;
	float *CoordinateY_out;
	float *CoordinateX_in;
	float *CoordinateY_in;	


	CoordinateX_out=(float*)calloc(irmax_out[0]*irmax_out[1],sizeof(float));
	CoordinateY_out=(float*)calloc(irmax_out[0]*irmax_out[1],sizeof(float));
	CoordinateX_in=(float*)calloc(irmax_in[0]*irmax_in[1],sizeof(float));
	CoordinateY_in=(float*)calloc(irmax_in[0]*irmax_in[1],sizeof(float));	

	printf("Creating TranslationMatrix\n");
	printf("Reading CoordinateX&Y of grid1...\n");
	CG( cg_coord_read(index_file_in,index_base,index_zone,"CoordinateX",RealSingle,irmin_in,irmax_in,CoordinateX_in));
	CG( cg_coord_read(index_file_in,index_base,index_zone,"CoordinateY",RealSingle,irmin_in,irmax_in,CoordinateY_in));
	printf("Reading CoordinateX&Y of grid2...\n");
	CG( cg_coord_read(index_file_out,index_base,index_zone,"CoordinateX",RealSingle,irmin_out,irmax_out,CoordinateX_out));
	CG( cg_coord_read(index_file_out,index_base,index_zone,"CoordinateY",RealSingle,irmin_out,irmax_out,CoordinateY_out));
	
	irmax_in[2]=zonesize_in[0][2];
	irmax_out[2]=zonesize_out[0][2];
	
	int i_in,j_in,k,ijk_in;
	int i_out,j_out,ijk_out;
	//int *translationmatrix;
	//translationmatrix=(int*)calloc(irmax_out[0]*irmax_out[1]*irmax_out[2],sizeof(int));
	int *topLeftIJK=(int*)calloc(irmax_out[0]*irmax_out[1],sizeof(int));
	int *topRightIJK=(int*)calloc(irmax_out[0]*irmax_out[1],sizeof(int));
	int *bottomLeftIJK=(int*)calloc(irmax_out[0]*irmax_out[1],sizeof(int));
	int *bottomRightIJK=(int*)calloc(irmax_out[0]*irmax_out[1],sizeof(int));
	
	float *topLeftWeighting=(float*)calloc(irmax_out[0]*irmax_out[1],sizeof(float));
	float *topRightWeighting=(float*)calloc(irmax_out[0]*irmax_out[1],sizeof(float));
	float *bottomLeftWeighting=(float*)calloc(irmax_out[0]*irmax_out[1],sizeof(float));
	float *bottomRightWeighting=(float*)calloc(irmax_out[0]*irmax_out[1],sizeof(float));	
	
	float distance;
	float distance_smallest;
	printf("Searching for corresponding meshpoints within k=0 plane...\n");
	k=0;

	float i_faktor;
	float j_faktor;
	float k_faktor;
	i_faktor=(float)zonesize_out[0][0]/zonesize_in[0][0];	
	j_faktor=(float)zonesize_out[0][1]/zonesize_in[0][1];	
	k_faktor=(float)zonesize_out[0][2]/zonesize_in[0][2];

	printf("Verhältnisse out/in i-j-k: %f %f %f\n",i_faktor,j_faktor,k_faktor);

	int i_start_in,i_end_in;		
	int j_start_in,j_end_in;

	int i_found;
	int j_found;	
	
	int i,i_situation;
	int j,j_situation;
	int ijk;
	
	int flag_corners;
	int flag_tr;
	int flag_tl;
	int flag_br;
	int flag_bl;
	int search_size;	
	float area_tr;
	float area_tl;
	float area_br;
	float area_bl;
	
	int i_delta=15;
	int j_delta=15;
	j_start_in=1;
	j_end_in=30;
	int found_nall=0;
	int found_all=0;
	
	for(j_out=0;j_out<irmax_out[1];j_out++)
	{
		i_start_in=1;
		i_end_in=30;
		for(i_out=0;i_out<irmax_out[0];i_out++)
		{
			ijk_out=k*zonesize_out[0][1]*zonesize_out[0][0]+j_out*zonesize_out[0][0]+i_out;
			distance_smallest=9999999999.;
			
			if(i_start_in<1){i_start_in=1;}
			if(j_start_in<1){j_start_in=1;}
			if(i_end_in>=irmax_in[0]-1){i_end_in=irmax_in[0]-1;}
			if(j_end_in>=irmax_in[1]-1){j_end_in=irmax_in[1]-1;}

//			Finde den Punkt in 'in' der am nächsten zum aktuellen Punkt von 'out' ist

			for(i_in=i_start_in;i_in<i_end_in;i_in++)
			{
				for(j_in=j_start_in;j_in<j_end_in;j_in++)
				{
					ijk_in=k*zonesize_in[0][1]*zonesize_in[0][0]+j_in*zonesize_in[0][0]+i_in;
				
					distance=sqrt(
					(CoordinateX_in[ijk_in]-CoordinateX_out[ijk_out])*(CoordinateX_in[ijk_in]-CoordinateX_out[ijk_out])
					+(CoordinateY_in[ijk_in]-CoordinateY_out[ijk_out])*(CoordinateY_in[ijk_in]-CoordinateY_out[ijk_out])
					);
					if(distance_smallest>distance)
					{
						i_found=i_in;
						j_found=j_in;
						
						distance_smallest=distance;
					}
				}
			}
			//if((i_out==100)&&(j_out<50)){printf("ifound:%d jfound:%d \n",i_found,j_found);printf("i:%d-%d | j:%d-%d     --- distance:%f\n",i_start_in,i_end_in,j_start_in,j_end_in,distance_smallest);}
			
			//Festlegung des Such-Intervalls für nächsten 'out'Punkt
			i_delta=(int)(15.+10./(1.+i_found-i_start_in))/i_faktor;
			j_delta=(int)(15.+10./(1.+j_found-j_start_in))/j_faktor;
			
			//printf("i_delta: %d - j_delta: %d\n",i_delta,j_delta);
			
			i_start_in=i_found-i_delta;
			i_end_in=i_found+i_delta;
			j_start_in=j_found-j_delta;
			j_end_in=j_found+j_delta;
			
			//i_start_in=i_found-15;
			//i_end_in=i_found+15;
			//j_start_in=j_found-15;
			//j_end_in=j_found+15;			
			
			//printf("ifound:%d jfound:%d \n",i_found,j_found);

//			Überprüfe wo der 'in' Punkt bezogen zum 'out' Punkt liegt (tr,tl,br,bl)

			//Punkte liegen quasi übereinander
			if(distance_smallest<1.0e-10)
			{
			ijk_in=j_found*zonesize_in[0][0]+i_found;
			topLeftIJK[ijk_out]=ijk_in;
			topLeftWeighting[ijk_out]=1.0;
			
			topRightWeighting[ijk_out]=0.0;
			bottomLeftWeighting[ijk_out]=0.0;
			bottomRightWeighting[ijk_out]=0.0;
			}
			else
			{
			flag_corners=0;
			flag_tr=0;	flag_tl=0;	flag_br=0;	flag_bl=0;
			area_tr=0.0;	area_tl=0.0;	area_br=0.0;	area_bl=0.0;

			distance_smallest=999999.;
			for(i=-1;i<2;i=i+2){for(j=-1;j<2;j=j+2){
				i_in=i_found+i;	j_in=j_found+j;
				ijk_in=j_in*zonesize_in[0][0]+i_in;			
			
				distance=sqrt(
				(CoordinateX_in[ijk_in]-CoordinateX_out[ijk_out])*(CoordinateX_in[ijk_in]-CoordinateX_out[ijk_out])
				+(CoordinateY_in[ijk_in]-CoordinateY_out[ijk_out])*(CoordinateY_in[ijk_in]-CoordinateY_out[ijk_out])
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
			ijk_in=(j_found+j_situation)*zonesize_in[0][0]+(i_found+i_situation);
			area_tr=fabs((CoordinateX_in[ijk_in]-CoordinateX_out[ijk_out])*(CoordinateY_in[ijk_in]-CoordinateY_out[ijk_out]));
			topRightIJK[ijk_out]=ijk_in;
			
			// topLeft
			ijk_in=(j_found+j_situation)*zonesize_in[0][0]+(i_found+i_situation-1);
			area_tl=fabs((CoordinateX_in[ijk_in]-CoordinateX_out[ijk_out])*(CoordinateY_in[ijk_in]-CoordinateY_out[ijk_out]));
			topLeftIJK[ijk_out]=ijk_in;
			
			// bottomRight
			ijk_in=(j_found+j_situation-1)*zonesize_in[0][0]+(i_found+i_situation);
			area_br=fabs((CoordinateX_in[ijk_in]-CoordinateX_out[ijk_out])*(CoordinateY_in[ijk_in]-CoordinateY_out[ijk_out]));
			bottomRightIJK[ijk_out]=ijk_in;	
			
			// bottomLeft
			ijk_in=(j_found+j_situation-1)*zonesize_in[0][0]+(i_found+i_situation-1);
			area_bl=fabs((CoordinateX_in[ijk_in]-CoordinateX_out[ijk_out])*(CoordinateY_in[ijk_in]-CoordinateY_out[ijk_out]));
			bottomLeftIJK[ijk_out]=ijk_in;

			
			topLeftWeighting[ijk_out]=area_br/(area_br+area_bl+area_tr+area_tl);
			topRightWeighting[ijk_out]=area_bl/(area_br+area_bl+area_tr+area_tl);
			bottomLeftWeighting[ijk_out]=area_tr/(area_br+area_bl+area_tr+area_tl);
			bottomRightWeighting[ijk_out]=area_tl/(area_br+area_bl+area_tr+area_tl);
			}
			/*printf("Gefundene Ecken für Position x:%f y:%f\n bl:x:%f y:%f\n  br:x:%f y:%f\n  tl:x:%f y:%f\n  tr:x:%f y:%f\n ",
			CoordinateX_out[ijk_out],CoordinateY_out[ijk_out],
			CoordinateX_in[bottomLeftIJK[ijk_out]],CoordinateY_in[bottomLeftIJK[ijk_out]],
			CoordinateX_in[bottomRightIJK[ijk_out]],CoordinateY_in[bottomRightIJK[ijk_out]],
			CoordinateX_in[topLeftIJK[ijk_out]],CoordinateY_in[topLeftIJK[ijk_out]],
			CoordinateX_in[topRightIJK[ijk_out]],CoordinateY_in[topRightIJK[ijk_out]]);*/
		}
	}
	printf("Free memory of CoordinateX&Y-arrays\n");
	free(CoordinateX_out);
	free(CoordinateY_out);
	free(CoordinateX_in);
	free(CoordinateY_in);
	
	int ijk0;
	int k_out;

	char iterationstr[ 255 ];
	int int_actualIteration=0;
	snprintf( iterationstr,255,"%i",int_actualIteration );
	CG( cg_goto( index_file_out,index_base,"end" ) );
	CG( cg_descriptor_write( "Iterations",iterationstr ) );

	int samples=1;
	float* timearray;
	timearray = malloc( ( samples )*sizeof( float ) );
	timearray[0]=0.0;
	CG( cg_biter_write( index_file_out,index_base,"BaseIterativeData",samples ) );
	CG( cg_goto( index_file_out,index_base,"BaseIterativeData",0,"end" ) );
	CG( cg_array_write( "TimeValues",RealSingle,1,( cgsize_t[ ] ){ samples },timearray ) );

	char solname[ 33 ];
	snprintf( solname,33,"Sol_%04d",0 );
	CG( cg_ziter_write( index_file_out,index_base,1,"ZoneIterativeData" ) );
	CG( cg_goto( index_file_out,index_base,"Zone_t",1,"ZoneIterativeData_t",1,"end" ) );
	CG( cg_array_write( "FlowSolutionPointers",Character,2,(cgsize_t[ ]){ 32,samples },solname ) );	
	float *buffer_out;
	float *buffer_in;
	long long int irmax_out0=irmax_out[0];
	long long int irmax_out1=irmax_out[1];
	long long int irmax_out2=irmax_out[2];
	long long int size_out=irmax_out0*irmax_out1*irmax_out2;
	printf("Allocating memory for buffer-array (size_out:%lld)\n",size_out);
	buffer_out=(float*)calloc(size_out,sizeof(float));
	buffer_in=(float*)calloc(irmax_in[0]*irmax_in[1]*irmax_in[2],sizeof(float));
	
	int n;
	CG( cg_sol_write( index_file_out,index_base,index_zone,solname,Vertex,&n ) );
	int field;
	int nfields;
	int index_field;
	char fieldname[ 33 ];
	int k_offset=zonesize_in[0][1]*zonesize_in[0][0];
	DataType_t datatype;
	CG( cg_nfields( index_file_in,index_base,index_zone,index_sol,&nfields ) );
	printf("Es gibt %d Variablen.\n",nfields);
	nfields=5;
	printf("Es werden %d Variablen übertragen.\n",nfields);
	for( field = 1; field<=nfields; field++ ) 
	{
		cg_field_info(index_file_in,index_base,index_zone,index_sol, field, &datatype, fieldname);
		printf("%s...\n",fieldname);
		printf("reading...\n");
		cg_field_read(index_file_in,index_base,index_zone,index_sol,fieldname,RealSingle,irmin_in,irmax_in,buffer_in);
		printf("copying...\n");
		copyGrid1ToGrid2(k_offset,k_faktor,irmax_out,irmax_in,buffer_out,buffer_in,
			topLeftIJK,topRightIJK,bottomLeftIJK,bottomRightIJK,
			topLeftWeighting,topRightWeighting,bottomLeftWeighting,bottomRightWeighting);
		printf("writing...\n");
		cg_field_write(index_file_out,index_base,index_zone,index_sol,RealSingle,fieldname,buffer_out,&index_field);
	}
	
	
	printf("Closing files...\n");
	CG( cg_close(index_file_in));
	CG( cg_close(index_file_out));
}

