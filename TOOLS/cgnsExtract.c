#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"
#include "unistd.h"

#define CG(cmd) if(cmd)cg_error_exit( );
#define NUMBERPOSITIONS 11

/*
* Check if a file exist using fopen() function
* return 1 if the file exist otherwise return 0
*/
int cfileexists(const char * filename){
    /* try to open file to read */
    FILE *file;
    if (file = fopen(filename, "r")){
        fclose(file);
        return 1;
    }
    return 0;
}

int main(int argc, char *argv[])
{
	int c;
	int fileflag=0;
	opterr = 0;
	char actual_file[500];
	while ((c = getopt (argc, argv, "f:")) != -1)
	{
		switch (c)
		{
			case 'f':
				fileflag+=1;
				strcpy(actual_file,optarg);
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
	if(fileflag!=1)
	{
		printf("ERROR: Kein Datei angegeben (-f 'file.cgns')!\n");
		errorflag=1;
	}

	int file,index_file_out,index_file_in;
	char path[500];
	char destination_file[500];
	
	float p_test=0.0;
	
	
	int index_base=1;
	int index_zone,index_zone_1;
	int number_zones;
	int index_solution;
	int number_solutions;
	int index_field;
	int i;
	int numberFiles;

strcpy(destination_file,"/home/mag/WORK/Post_BAC_Clean/Total_TimeAverage_Film_3d_bac_2048x256x128_VG_Re1E6_periodic_ct_16384CPU_sdrc.cgns");
	
	float * u,*v,*w,*p,*rho,*mach,*gradrho,*Lambda2;
	cgsize_t idata[2],nuse;

	char zonename[32];
	char baseIterName[100];
	
	cgsize_t zonesize[3][3];

	cgsize_t irmin[3];
	cgsize_t irmax[3];

		/* lower range index */
	irmin[0]=1;
	irmin[1]=1;
	irmin[2]=1;
	
	int step,startStep,endStep,stepSize;
	float *time;
	float *time_0;
	float time_factor;
	time_factor=0.000008*30000./10.*0.08/(0.76*sqrt(287.*300.*1.4));

	float x_extr[NUMBERPOSITIONS],y_extr[NUMBERPOSITIONS],z_extr[NUMBERPOSITIONS];
	float *p_extr[NUMBERPOSITIONS];
	int position[NUMBERPOSITIONS];
	float last_p[NUMBERPOSITIONS];
	float distance[NUMBERPOSITIONS];

		

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	startStep=455;	endStep=1325;	stepSize=30;
	
	x_extr[0]=0.0;	y_extr[0]=0.1;	z_extr[0]=0.01;	distance[0]=9999.;
	x_extr[1]=0.1;	y_extr[1]=0.1;	z_extr[1]=0.01;	distance[1]=9999.;
	x_extr[2]=0.2;	y_extr[2]=0.1;	z_extr[2]=0.01;	distance[2]=9999.;
	x_extr[3]=0.3;	y_extr[3]=0.1;	z_extr[3]=0.01;	distance[3]=9999.;
	x_extr[4]=0.4;	y_extr[4]=0.1;	z_extr[4]=0.01;	distance[4]=9999.;
	x_extr[5]=0.5;	y_extr[5]=0.1;	z_extr[5]=0.01;	distance[5]=9999.;
	x_extr[6]=0.6;	y_extr[6]=0.1;	z_extr[6]=0.01;	distance[6]=9999.;
	x_extr[7]=0.7;	y_extr[7]=0.1;	z_extr[7]=0.01;	distance[7]=9999.;
	x_extr[8]=0.8;	y_extr[8]=0.1;	z_extr[8]=0.01;	distance[8]=9999.;
	x_extr[9]=0.9;	y_extr[9]=0.1;	z_extr[9]=0.01;	distance[9]=9999.;
	x_extr[10]=1.0;	y_extr[10]=0.1;	z_extr[10]=0.01;distance[10]=9999.;	
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	
	step=startStep;
	

	int pos;
	int buffer;


	numberFiles=1;
	
	CG( cg_open(destination_file,CG_MODE_MODIFY,&index_file_out));
		
	for (step=startStep;step<=endStep;step+=stepSize)
	{
		sprintf(path,"/home/mag/NAS/RESULTS/VG/BAC_pur/16384CPU/%dk-%dk/",(step-30),step);
		sprintf(actual_file,"%sFilm_3d_bac_2048x256x128_VG_Re1E6_periodic_ct_16384CPU_sdrc.cgns",path);
		if(cfileexists(actual_file))
		{
		numberFiles++;
		printf("Opening %s...",actual_file);
		printf("Datei: %d/%d\n",(step-startStep)/stepSize+1,(endStep-startStep)/stepSize+1);
		CG( cg_open(actual_file,CG_MODE_READ,&index_file_in));
		CG( cg_nzones(index_file_in,index_base,&number_zones));
		if (step==startStep)
		{
			float * x,*y,*z;
			for(index_zone=1;index_zone<=number_zones;index_zone++)
			{
				CG( cg_zone_read(index_file_in,index_base,index_zone,zonename,zonesize[0]));				
				irmax[0]=zonesize[0][0];
				irmax[1]=zonesize[0][1];
				irmax[2]=zonesize[0][2];

				buffer=zonesize[0][0]*zonesize[0][1]*zonesize[0][2];
				
				x=(float *)calloc(buffer, sizeof(float));
				y=(float *)calloc(buffer, sizeof(float));
				z=(float *)calloc(buffer, sizeof(float));

				CG( cg_coord_read(index_file_in,index_base,index_zone,"CoordinateX",RealSingle,irmin,irmax,x));
				CG( cg_coord_read(index_file_in,index_base,index_zone,"CoordinateY",RealSingle,irmin,irmax,y));
				CG( cg_coord_read(index_file_in,index_base,index_zone,"CoordinateZ",RealSingle,irmin,irmax,z));
				for(pos=0;pos<NUMBERPOSITIONS;pos++)
				{
					for(i=0;i<buffer;i++)
					{
						if(sqrt(pow((x[i]-x_extr[pos]),2.)+pow((y[i]-y_extr[pos]),2.)+pow((z[i]-z_extr[pos]),2.))<distance[pos])
						{
							distance[pos]=sqrt(pow((x[i]-x_extr[pos]),2.)+pow((y[i]-y_extr[pos]),2.)+pow((z[i]-z_extr[pos]),2.));
							position[pos]=i;
							index_zone_1=index_zone;
						}
					}
					printf("actual result: zone %d: x=%g y=%g z=%g\n",index_zone_1,x[position[pos]],y[position[pos]],z[position[pos]]);
				}
		
				free(x);
				free(y);
				free(z);
			}
		}
		

		
		index_zone=index_zone_1;

		CG( cg_biter_read(index_file_in, index_base, baseIterName, &number_solutions));
		
		if (step==startStep)
		{
			time_0=(float *)calloc((number_solutions*((endStep-startStep)/stepSize+1)), sizeof(float));
			for(pos=0;pos<NUMBERPOSITIONS;pos++)
			{
				p_extr[pos]=(float *)calloc((number_solutions*((endStep-startStep)/stepSize+1)), sizeof(float));
			}
		}
		time=(float *)calloc(number_solutions, sizeof(float));
		CG( cg_goto(index_file_in,index_base,"BaseIterativeData_t",1,"end"));
		CG( cg_array_read_as(1,RealSingle,time));
		for(index_solution=1;index_solution<=number_solutions;index_solution++)
		{
			printf("Solution %d/%d \r",index_solution,number_solutions);
			fflush(stdout);
			CG( cg_zone_read(index_file_in,index_base,index_zone,zonename,zonesize[0]));				
			irmax[0]=zonesize[0][0];
			irmax[1]=zonesize[0][1];
			irmax[2]=zonesize[0][2];

			buffer=zonesize[0][0]*zonesize[0][1]*zonesize[0][2];
			//u=(float *)calloc(buffer, sizeof(float));
			//v=(float *)calloc(buffer, sizeof(float));
			//w=(float *)calloc(buffer, sizeof(float));
			p=(float *)calloc(buffer, sizeof(float));
			//rho=(float *)calloc(buffer, sizeof(float));
			//gradrho=(float *)calloc(buffer, sizeof(float));
			//mach=(float *)calloc(buffer, sizeof(float));
			//Lambda2=(float *)calloc(buffer, sizeof(float));
			
	
			//printf("Reading File %d/%d...\n",step,endStep);
			/*cg_field_read(index_file_in,index_base,index_zone,index_solution,"VelocityX",
				RealSingle,irmin,irmax,u);*/
			/*cg_field_read(index_file_in,index_base,index_zone,index_solution,"VelocityY",
				RealSingle,irmin,irmax,v);*/
			/*cg_field_read(index_file_in,index_base,index_zone,index_solution,"VelocityZ",
				RealSingle,irmin,irmax,w);*/
			cg_field_read(index_file_in,index_base,index_zone,index_solution,"Pressure",
				RealSingle,irmin,irmax,p);
			/*cg_field_read(index_file_in,index_base,index_zone,index_solution,"Density",
				RealSingle,irmin,irmax,rho);*/
			/*cg_field_read(index_file_in,index_base,index_zone,index_solution,"DensityGradient",
				RealSingle,irmin,irmax,gradrho);*/
			/*cg_field_read(index_file_in,index_base,index_zone,index_solution,"MachNumber",
				RealSingle,irmin,irmax,mach);*/
			/*cg_field_read(index_file_in,index_base,index_zone,index_solution,"Lambda2",
				RealSingle,irmin,irmax,Lambda2);*/

			time_0[number_solutions*((step-startStep)/stepSize)+(index_solution-1)]=time[index_solution-1];
			for(pos=0;pos<NUMBERPOSITIONS;pos++)
			{
				p_extr[pos][number_solutions*((step-startStep)/stepSize)+(index_solution-1)]=p[position[pos]];
				last_p[pos]=p[position[pos]];
			}
	
			//printf("Druck %g\n",p_1[number_solutions*((step-startStep)/stepSize)+(index_solution-1)]);
			
			//free(u);
			//free(v);
			//free(w);
			free(p);
			//free(rho);
			//free(Lambda2);
			//free(mach);
			//free(gradrho);			
			time_0[number_solutions*((step-startStep)/stepSize)+(index_solution-1)]=(number_solutions*((step-startStep)/stepSize)+(index_solution-1))*time_factor;
		}
		free(time);
		
		CG( cg_close(index_file_in));
		}
		else
		{
			printf("%s and not found!\n",actual_file);
			for(index_solution=1;index_solution<=number_solutions;index_solution++)
			{
				time_0[number_solutions*((step-startStep)/stepSize)+(index_solution-1)]=(number_solutions*((step-startStep)/stepSize)+(index_solution-1))*time_factor;
			for(pos=0;pos<NUMBERPOSITIONS;pos++)
			{
				p_extr[pos][number_solutions*((step-startStep)/stepSize)+(index_solution-1)]=last_p[pos];
			}
			}
			
		}


	}
	
	CG( cg_close(index_file_out));
	
//	Zeit ist mit offsets gespeichert. fehler
	
	

	FILE * file0;
	sprintf(actual_file,"PressureHistory_ExtractFromFilm.dat");
	file0=fopen(actual_file,"w");
	fprintf(file0,"TITLE = \"PressureHistory von FilmDateien\"\n");
	fprintf(file0,"VARIABLES = \"Time (S)\" \"Pressure\"\n");



	for(pos=0;pos<NUMBERPOSITIONS;pos++)
	{
		fprintf(file0,"ZONE T=\"Point%d(x=%f, y=%f, z=%f)\", F=POINT, I=%d, DT=(SINGLE)\n",
				1,
				x_extr[pos],
				y_extr[pos],
				z_extr[pos],
				(number_solutions*((endStep-startStep)/stepSize+1)));
		for(i=0;i<(number_solutions*((endStep-startStep)/stepSize+1));i++)
		{
			fprintf(file0,"%le %le\n",time_0[i],p_extr[pos][i]);
		}
	}


	fclose(file0);
	
}



