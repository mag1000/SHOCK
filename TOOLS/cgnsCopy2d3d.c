#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"
#include "unistd.h"

#define CG(cmd) if(cmd)cg_error_exit( );


int main(int argc, char *argv[])
{
	int index_file;
	char input_file[500];
	char output_file[500];	
	//sprintf(input_file,"/home/mag/CLUSTER_JUQUEEN2/NACA3D/2d_NACA0012_yMin20_500000_4096x256.cgns");
	//sprintf(output_file,"/home/mag/CLUSTER_JUQUEEN2/NACA3D/3d_NACA0012_yMin20_500000_4096x512x256_2D.cgns");
	
	sprintf(input_file,"/home/mag/WORK2/NACA3D_M065a4/2d_NACA0012_yMin40_500000_2048x256.cgns");
	sprintf(output_file,"/home/mag/WORK2/NACA3D_M065a4/3d_NACA0012_yMin40_500000_2048x256x128.cgns");
	printf("2D-Daten werden in 3D-Daten kopiert...\n\n%s ---> %s\n\n\n",input_file,output_file);



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////2D - DATEI//////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
	CG( cg_open(input_file,CG_MODE_READ,&index_file));

	int index_base=1;
	char baseIterName[100];	
	int number_solutions;
	CG( cg_biter_read(index_file, index_base, baseIterName, &number_solutions));	

	float *time,*time_new;
	time=(float *)calloc(number_solutions, sizeof(float));
	CG( cg_goto(index_file,index_base,"BaseIterativeData_t",1,"end"));
	CG( cg_array_read_as(1,RealSingle,time));
	time_new=(float *)calloc(1, sizeof(float));
	time_new[0]=time[number_solutions-1];

	char solname[ 32 ];
	GridLocation_t gl;
	int index_zone=1;
	CG( cg_sol_info( index_file,index_base,index_zone,number_solutions,solname,&gl ) );
	printf("Ergebnisse namens %s des Zeitschrittes %e werden gelesen...\n",solname,time_new[0]);
		
	cgsize_t zonesize[2][2];
	char zonename[32];
	CG( cg_zone_read(index_file,index_base,index_zone,zonename,zonesize[0]));
	
	cgsize_t irmin[2];
	cgsize_t irmax[2];
	irmin[0]=1;
	irmin[1]=1;	
	irmax[0]=zonesize[0][0];
	irmax[1]=zonesize[0][1];
	int buffer=zonesize[0][0]*zonesize[0][1];
	float *u,*v,*p,*rho;
	u=(float *)calloc(buffer, sizeof(float));
	v=(float *)calloc(buffer, sizeof(float));
	p=(float *)calloc(buffer, sizeof(float));
	rho=(float *)calloc(buffer, sizeof(float));
	printf("Lade Ergebnisse aus Zone %s(%zux%zu)...\n",zonename,irmax[0],irmax[1]);
	CG( cg_field_read(index_file,index_base,index_zone,number_solutions,"VelocityX",
	RealSingle,irmin,irmax,u));
	CG( cg_field_read(index_file,index_base,index_zone,number_solutions,"VelocityY",
	RealSingle,irmin,irmax,v));
	CG( cg_field_read(index_file,index_base,index_zone,number_solutions,"Pressure",
	RealSingle,irmin,irmax,p));
	CG( cg_field_read(index_file,index_base,index_zone,number_solutions,"Density",
	RealSingle,irmin,irmax,rho));
	
	CG( cg_close(index_file));

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////2D - DATEI ENDE/////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////	
	
	
	
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////3D - DATEI//////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
	CG( cg_open(output_file,CG_MODE_MODIFY,&index_file));

	cgsize_t zonesize_3d[3][3];
	CG( cg_zone_read(index_file,index_base,index_zone,zonename,zonesize_3d[0]));
	
	cgsize_t irmin_3d[3];
	cgsize_t irmax_3d[3];
	irmin_3d[0]=1;
	irmin_3d[1]=1;
	irmin_3d[2]=1;		
	irmax_3d[0]=zonesize_3d[0][0];
	irmax_3d[1]=zonesize_3d[0][1];
	irmax_3d[2]=zonesize_3d[0][2];
	buffer=zonesize_3d[0][0]*zonesize_3d[0][1]*zonesize_3d[0][2];
	printf("Anzahl Gitterpunkte in 3D-Gitter: %d\n",buffer);
	float *u_3d,*v_3d,*w_3d,*p_3d,*rho_3d;
	u_3d=(float *)calloc(buffer, sizeof(float));
	v_3d=(float *)calloc(buffer, sizeof(float));
	w_3d=(float *)calloc(buffer, sizeof(float));
	p_3d=(float *)calloc(buffer, sizeof(float));
	rho_3d=(float *)calloc(buffer, sizeof(float));
	printf("Schreibe Ergebnisse in Zone %s(%zux%zux%zu)...\n",zonename,irmax_3d[0],irmax_3d[1],irmax_3d[2]);
	
	int imax=irmax_3d[0];
	int jmax=irmax_3d[1];
	int kmax=irmax_3d[2];	
	int i,j,k,ijk,ijk_3d;
	int index_solution;
	int index_field;		
	for(i=0;i<imax;i++)
	{
		for(j=0;j<jmax;j++)
		{
			for(k=0;k<kmax;k++)
			{
				ijk_3d=k*jmax*imax+j*imax+i;
				ijk=j*imax+i;

				u_3d[ijk_3d]=u[ijk];
				v_3d[ijk_3d]=v[ijk];				
				p_3d[ijk_3d]=p[ijk];
				rho_3d[ijk_3d]=rho[ijk];
				if(j<20)
				{
					w_3d[ijk_3d]=0.01*sin(k*5.0/kmax*2.0*M_PI);
				}
				else
				{
					w_3d[ijk_3d]=0.0;
				}
			}
		}
	}

	CG( cg_sol_write( index_file,index_base,index_zone,solname,Vertex,&index_solution ) );
	CG( cg_field_write(index_file,index_base,index_zone,index_solution,RealSingle,"VelocityX",u_3d,&index_field)) ;
	CG( cg_field_write(index_file,index_base,index_zone,index_solution,RealSingle,"VelocityY",v_3d,&index_field)) ;
	CG( cg_field_write(index_file,index_base,index_zone,index_solution,RealSingle,"VelocityZ",w_3d,&index_field)) ;
	CG( cg_field_write(index_file,index_base,index_zone,index_solution,RealSingle,"Density",rho_3d,&index_field)) ;
	CG( cg_field_write(index_file,index_base,index_zone,index_solution,RealSingle,"Pressure",p_3d,&index_field)) ;



	
	CG( cg_close(index_file));
	
	return 0;
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////3D - DATEI ENDE/////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////		
}

