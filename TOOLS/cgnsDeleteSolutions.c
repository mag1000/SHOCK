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
	char file[500];
	char command[500];
	strcpy(file,argv[1]);
	int i;
	int index_base=1;
	int index_zone;
	int number_zones;
	int index_solution;
	int number_solutions;
	int index_flow;
	int index_field;
	int nsteps;
	int int_meshDimensions;
	CG( cg_open(file,CG_MODE_MODIFY,&index_file));
	CG (cg_cell_dim(index_file, index_base, &int_meshDimensions));
	CG( cg_nzones(index_file,index_base,&number_zones));
	

	int buffer;
	float *u,*v,*w,*p,*rho,*mach,*gradrho,*Lambda2;
	cgsize_t idata[2],nuse;

	char zonename[32];
	char baseIterName[100];
	GridLocation_t gl;
	char solname[ 32 ];
	char solname_first[ 32 ];

	cgsize_t zonesize[int_meshDimensions][int_meshDimensions];

	cgsize_t irmin[int_meshDimensions];
	cgsize_t irmax[int_meshDimensions];

	/* lower range index */
	irmin[0]=1;
	irmin[1]=1;
	if (int_meshDimensions==3){irmin[2]=1;}

	float *time,*time_new;




	float * x,*y,*z;
	int nsols;
	CG( cg_biter_read(index_file, index_base, baseIterName, &number_solutions));
	time=(float *)calloc(number_solutions, sizeof(float));
	CG( cg_goto(index_file,index_base,"BaseIterativeData_t",1,"end"));
	CG( cg_array_read_as(1,RealSingle,time));
	
	time_new=(float *)calloc(1, sizeof(float));
	time_new[0]=time[number_solutions-1];
	printf("last timestep: %e\n",time_new[0]);
	for(index_zone=1;index_zone<=number_zones;index_zone++)
	{
		CG( cg_zone_read(index_file,index_base,index_zone,zonename,zonesize[0]));				
		irmax[0]=zonesize[0][0];
		irmax[1]=zonesize[0][1];
		if (int_meshDimensions==3){irmax[2]=zonesize[0][2];}
		
		if (int_meshDimensions==2){buffer=zonesize[0][0]*zonesize[0][1];}
		if (int_meshDimensions==3){buffer=zonesize[0][0]*zonesize[0][1]*zonesize[0][2];}		

		
		u=(float *)calloc(buffer, sizeof(float));
		v=(float *)calloc(buffer, sizeof(float));
		w=(float *)calloc(buffer, sizeof(float));
		p=(float *)calloc(buffer, sizeof(float));
		rho=(float *)calloc(buffer, sizeof(float));
		gradrho=(float *)calloc(buffer, sizeof(float));
		mach=(float *)calloc(buffer, sizeof(float));
		Lambda2=(float *)calloc(buffer, sizeof(float));

		printf("Loading last solution of zone %s\n",zonename);
		cg_field_read(index_file,index_base,index_zone,number_solutions,"VelocityX",
		RealSingle,irmin,irmax,u);
		cg_field_read(index_file,index_base,index_zone,number_solutions,"VelocityY",
		RealSingle,irmin,irmax,v);
		cg_field_read(index_file,index_base,index_zone,number_solutions,"VelocityZ",
		RealSingle,irmin,irmax,w);
		cg_field_read(index_file,index_base,index_zone,number_solutions,"Pressure",
		RealSingle,irmin,irmax,p);
		cg_field_read(index_file,index_base,index_zone,number_solutions,"Density",
		RealSingle,irmin,irmax,rho);
		cg_field_read(index_file,index_base,index_zone,number_solutions,"DensityGradient",
		RealSingle,irmin,irmax,gradrho);
		cg_field_read(index_file,index_base,index_zone,number_solutions,"MachNumber",
		RealSingle,irmin,irmax,mach);
		cg_field_read(index_file,index_base,index_zone,number_solutions,"Lambda2",
		RealSingle,irmin,irmax,Lambda2);
		
		CG( cg_nsols( index_file,index_base,index_zone,&nsols ) );
		printf("Deleting %d solutions for zone %s\n",nsols-1,zonename);
		CG( cg_goto( index_file,index_base,"Zone_t",index_zone,"end" ) );
		for( index_solution = 1; index_solution<=nsols; index_solution++ )
		{
			CG( cg_sol_info( index_file,index_base,index_zone,1,solname,&gl ) );
			CG( cg_delete_node( solname ) );
			if (index_solution==1)
			{
				strcpy(solname_first,solname);
			}
		}
		
		printf("Writing last solution for zone %s, ",zonename);
		CG( cg_sol_write(index_file,index_base,index_zone,solname_first,Vertex,&index_flow) );
		printf("solution Index %d\n",index_flow);
				
		CG( cg_field_write(index_file,index_base,index_zone,index_flow,
		RealSingle,"VelocityX",u,&index_field) );
		CG( cg_field_write(index_file,index_base,index_zone,index_flow,
		RealSingle,"VelocityY",v,&index_field) );
		CG( cg_field_write(index_file,index_base,index_zone,index_flow,
		RealSingle,"VelocityZ",w,&index_field) );
		CG( cg_field_write(index_file,index_base,index_zone,index_flow,
		RealSingle,"Density",rho,&index_field) );
		CG( cg_field_write(index_file,index_base,index_zone,index_flow,
		RealSingle,"Pressure",p,&index_field) );
		CG( cg_field_write(index_file,index_base,index_zone,index_flow,
		RealSingle,"DensityGradient",gradrho,&index_field) );
		CG( cg_field_write(index_file,index_base,index_zone,index_flow,
		RealSingle,"MachNumber",mach,&index_field) );
		CG( cg_field_write(index_file,index_base,index_zone,index_flow,
		RealSingle,"Lambda2",Lambda2,&index_field) );
		
		printf("Deleting ZoneIterativeData of zone %s\n",zonename);
		CG( cg_delete_node( "ZoneIterativeData" ) );
		
		printf("Writing new ZoneIterativeData of zone %s\n",zonename);
		/* create ZoneIterativeData */
		CG( cg_ziter_write(index_file,index_base,index_zone,"ZoneIterativeData") );
		/* go to ZoneIterativeData level and give info telling which */
		/* flow solution corresponds with which time (solname(1) corresponds */
		/* with time(1), solname(2) with time(2), and solname(3) with time(3)) */
		CG( cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneIterativeData_t",1,"end") );
		idata[0]=32;
		idata[1]=1;
		CG( cg_array_write("FlowSolutionPointers",Character,2,idata,solname_first) );

		free(u);
		free(v);
		free(w);
		free(p);
		free(rho);
		free(Lambda2);
		free(mach);
		free(gradrho);	
	}
	
	printf("Deleting BaseIterativeData\n");
	CG( cg_goto( index_file,index_base,"end" ) );
	CG( cg_delete_node( "BaseIterativeData" ) );		
	
	printf("Writing new BaseIterativeData\n");
	nsteps=1;
	CG( cg_biter_write(index_file,index_base,"BaseIterativeData",nsteps) );
	/* go to BaseIterativeData level and write time values */
	CG( cg_goto(index_file,index_base,"BaseIterativeData",0,"end") );
	nuse=1;
	CG( cg_array_write("TimeValues",RealSingle,1,&nuse,time_new) );


	CG( cg_close(index_file) );
	
	sprintf(command,"cgnscompress %s\n",file);
	system(command);
}



