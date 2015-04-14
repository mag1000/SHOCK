#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"
#include "unistd.h"
#include "FUNCTIONS.h"

#define CG(cmd) if(cmd)cg_error_exit( );

int main(int argc, char *argv[])
{
	char file[500];
	FILE* fileOut;
	int fileflag=0;		
	int c;
	int flag_convergence=0;
	int flag_pressureExtraction=0;
	int optionflag=0;
	opterr = 0;
	char zonename[100];
	strcpy(zonename,"Zone");
	float ReynoldsNumber=500000;
	float MachNumber=0.65;	
	
	float x_extr=0.5;
	float y_extr=-0.5;

	while ((c = getopt (argc, argv, "cf:r:m:z:p:")) != -1)
	{
		switch (c)
		{
			case 'p':
				flag_pressureExtraction=1;
				x_extr=atof(optarg);				
				optionflag++;
				break;			
			case 'c':
				flag_convergence=1;
				optionflag++;
				break;
			case 'm':
				MachNumber=atof(optarg);
				break;			
			case 'r':
				ReynoldsNumber=atof(optarg);
				break;
			case 'f':
				fileflag+=1;
				strcpy(file,optarg);
				break;						
			case 'z':
				strcpy(zonename,optarg);	
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
	printf("Starting postTransientData....\n");
	printf("Die Reynolds-Nummer (y+-Berechnung, option -r (FLOAT)) ist: %g\n",ReynoldsNumber);
	printf("Die Mach-Nummer (option -m (FLOAT)) ist: %g\n",MachNumber);			
	printf("Dieses Tool dient der Auswertung von transienten cgns-Dateien\n");
	int errorflag=0;
	if(fileflag!=1)
	{
		printf("ERROR: Kein Datei angegeben (-f 'file.cgns')!\n");
		errorflag=1;
	}
	
	if(optionflag<1)
	{
		printf("ERROR: Keine Auswerteoption angegeben: convergence: -c !\n");
		errorflag=1;
	}

	if(errorflag==1)
	{
		return 1;
	}
	
	
	printf("Öffne Datei %s im read mode\n",file);
	int index_file;
	int index_base=1;
	int NumberSolutions=1;
	int index_zone=1;
	int cell_dim;	
	char zonename_tmp[33];	
	CG( cg_open(file,CG_MODE_READ,&index_file));
	CG( cg_cell_dim(index_file,index_base,&cell_dim));
	
	
	cgsize_t zonesize[cell_dim][cell_dim];	
	
	CG( cg_zone_read(index_file,index_base,index_zone,zonename_tmp,zonesize[0]));

	cgsize_t irmin[cell_dim];
	cgsize_t irmax[cell_dim];
	int buffer;
	int iMeshPoints;
	int jMeshPoints;
	int kMeshPoints;		
	int i,j,k,ijk;
	irmin[0]=1;
	irmin[1]=1;	
	irmax[0]=zonesize[0][0];
	irmax[1]=zonesize[0][1];

	iMeshPoints=irmax[0];
	jMeshPoints=irmax[1];
	kMeshPoints=1;	
	
	if(cell_dim==3)
	{
		irmin[2]=1;	
		irmax[2]=zonesize[0][2];		

		kMeshPoints=irmax[2];
	}
	
	
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////convergence //////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////	
	if((flag_convergence==1)&&(optionflag==1))
	{
		irmax[1]=2;
		jMeshPoints=2;
		float *x_tmp;
		x_tmp=(float *)calloc(iMeshPoints*jMeshPoints*kMeshPoints, sizeof(float));
		CG( cg_coord_read(index_file,index_base,index_zone,"CoordinateX",
		RealSingle,irmin,irmax,x_tmp));
		int ijk,i1jk;		
		for(i=1;i<iMeshPoints-1;i++)
		{
			ijk=0*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
			i1jk=0*jMeshPoints*iMeshPoints+0*iMeshPoints+i-1;	
			if((x_tmp[ijk]<=1.1)&&(x_tmp[i1jk]>1.1))
			{
				irmin[0]=i;
			}
			if((x_tmp[ijk]>=1.1)&&(x_tmp[i1jk]<1.1))
			{
				irmax[0]=i;
			}
		}
		printf("Bereichseingrenzung auf Profiloberflaeche moeglich: i_start:%zu, i_ende:%zu\n",irmin[0],irmax[0]);
		iMeshPoints=irmax[0]-irmin[0]+1;
		free(x_tmp);
	}
	
	
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////pressureExtraction //////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////	
	if((flag_pressureExtraction==1)&&(optionflag==1))
	{
		printf("Suche Gitterpunkte nahe x=%f, y=%f\n",x_extr,y_extr);
		float distance=9999999.;
		float distance_tmp;
		irmax[2]=1;
		kMeshPoints=1;
		float *x_tmp,*y_tmp;
		x_tmp=(float *)calloc(iMeshPoints*jMeshPoints*kMeshPoints, sizeof(float));
		y_tmp=(float *)calloc(iMeshPoints*jMeshPoints*kMeshPoints, sizeof(float));		
		CG( cg_coord_read(index_file,index_base,index_zone,"CoordinateX",RealSingle,irmin,irmax,x_tmp));
		CG( cg_coord_read(index_file,index_base,index_zone,"CoordinateY",RealSingle,irmin,irmax,y_tmp));		
		int ijk,ijk_extr;	
		for(i=0;i<iMeshPoints;i++)
		{
		for(j=0;j<jMeshPoints;j++)
		{		
			ijk=0*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
			distance_tmp=sqrt((x_tmp[ijk]-x_extr)*(x_tmp[ijk]-x_extr)+(y_tmp[ijk]-y_extr)*(y_tmp[ijk]-y_extr));
			if(distance_tmp<distance)
			{
				distance=distance_tmp;
				irmin[0]=i+1;
				irmax[0]=i+1;
				irmin[1]=j+1;
				irmax[1]=j+1;
				ijk_extr=ijk;
			}
		}
		}
		printf("Druckwertextraktion bei : x=%f, y=%f\n",x_tmp[ijk_extr],y_tmp[ijk_extr]);
		iMeshPoints=1;
		jMeshPoints=1;
		free(x_tmp);
		free(y_tmp);		
	}	
	
	buffer=iMeshPoints*jMeshPoints*kMeshPoints;
	printf("Lade %dD-Ergebnisse aus Zone %s (%dx%dx%d)...\n",cell_dim,zonename_tmp,iMeshPoints,jMeshPoints,kMeshPoints);


	float *x,*y,*z;
	float *u_transient,*v_transient,*w_transient,*p_transient,*rho_transient;

	if((flag_pressureExtraction==1)&&(optionflag==1))
	{
		p_transient=(float *)calloc(1, sizeof(float));	
	}
	else
	{
		x=(float *)calloc(buffer, sizeof(float));
		y=(float *)calloc(buffer, sizeof(float));
		z=(float *)calloc(buffer, sizeof(float));
	
		u_transient=(float *)calloc(buffer, sizeof(float));
		v_transient=(float *)calloc(buffer, sizeof(float));
		w_transient=(float *)calloc(buffer, sizeof(float));
		p_transient=(float *)calloc(buffer, sizeof(float));
		rho_transient=(float *)calloc(buffer, sizeof(float));
	
	
		printf("Lade Gitter...\n");
	
		CG( cg_coord_read(index_file,index_base,index_zone,"CoordinateX",
		RealSingle,irmin,irmax,x));
		CG( cg_coord_read(index_file,index_base,index_zone,"CoordinateY",
		RealSingle,irmin,irmax,y));
		if(cell_dim==3)
		{
		CG( cg_coord_read(index_file,index_base,index_zone,"CoordinateZ",
		RealSingle,irmin,irmax,z));	
		}	
	}	


	
	CG( cg_nsols(index_file,index_base,1,&NumberSolutions ) );
	int index_solution;
	int index_sol;
	int index_field;
	//char solname[NumberSolutions*32+1];
	//char solnameArray[NumberSolutions][32];
	
	//strcpy(solname,"");
	
	float *time;
	time=(float *)calloc(NumberSolutions, sizeof(float));
	printf("Lade Zeit-array...\n");
	CG( cg_goto(index_file,1,"BaseIterativeData_t",1,"end"));
	CG( cg_array_read_as(1,RealSingle,time));
	
	float *p_extracted;
	p_extracted=(float *)calloc(NumberSolutions, sizeof(float));


	int field,nfields;
	char fieldname[ 33 ];
	DataType_t datatype;
	printf("Starte mit transientem postprocessing mit %d solutions\n",NumberSolutions);
	for(index_solution=1;index_solution<=NumberSolutions;index_solution++)
	{
		printf("Verarbeite Solution: %d von %d bei Zeit: %f\n",index_solution,NumberSolutions,time[index_solution-1]);
		if((flag_pressureExtraction==1)&&(optionflag==1))
		{
			CG( cg_field_read(index_file,index_base,index_zone,index_solution,"Pressure",
			RealSingle,irmin,irmax,p_transient));
			p_extracted[index_solution-1]=p_transient[0];
		}
		else
		{
			CG( cg_field_read(index_file,index_base,index_zone,index_solution,"VelocityX",
			RealSingle,irmin,irmax,u_transient));
			CG( cg_field_read(index_file,index_base,index_zone,index_solution,"VelocityY",
			RealSingle,irmin,irmax,v_transient));
			CG( cg_field_read(index_file,index_base,index_zone,index_solution,"VelocityZ",
			RealSingle,irmin,irmax,w_transient));	
			CG( cg_field_read(index_file,index_base,index_zone,index_solution,"Pressure",
			RealSingle,irmin,irmax,p_transient));
			CG( cg_field_read(index_file,index_base,index_zone,index_solution,"Density",
			RealSingle,irmin,irmax,rho_transient));
		}
		
		
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////convergence //////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		if(flag_convergence==1)
		{
			//LE:LeadingEdge
			//TE:TrailingEdge
			//SS:MidOfSuctionSide
			//PS:MidOfPressureSide
			float cp_LE=0.0;
			float cp_TE=0.0;
			float cp_SS=0.0;
			float cp_PS=0.0;
			float cf_LE=0.0;
			float cf_TE=0.0;
			float cf_SS=0.0;
			float cf_PS=0.0;
		

			float distance_LE=999999.;
			float distance_TE=999999.;
			float distance_PS=999999.;
			float distance_SS=999999.;
			float distance_LE_tmp;
			float distance_TE_tmp;
			float distance_PS_tmp;
			float distance_SS_tmp;
			int i_LE;
			int i_TE;
			int i_PS;
			int i_SS;
			for(i=0;i<iMeshPoints;i++)
			{
				ijk=0*jMeshPoints*iMeshPoints+0*iMeshPoints+i;

				distance_LE_tmp=sqrt((x[ijk]-0.0)*(x[ijk]-0.0)+(y[ijk]-0.0)*(y[ijk]-0.0));
				if(distance_LE_tmp<distance_LE)
				{
					i_LE=i;
					distance_LE=distance_LE_tmp;
				}
				
				distance_TE_tmp=sqrt((x[ijk]-1.0)*(x[ijk]-1.0)+(y[ijk]-0.0)*(y[ijk]-0.0));
				if(distance_TE_tmp<distance_TE)
				{
					i_TE=i;
					distance_TE=distance_TE_tmp;
				}
				
				distance_SS_tmp=sqrt((x[ijk]-0.5)*(x[ijk]-0.5)+(y[ijk]-0.0001)*(y[ijk]-0.0001));
				if(distance_SS_tmp<distance_SS)
				{
					i_SS=i;
					distance_SS=distance_SS_tmp;
				}
				
				distance_PS_tmp=sqrt((x[ijk]-0.5)*(x[ijk]-0.5)+(y[ijk]+0.0001)*(y[ijk]+0.0001));
				if(distance_PS_tmp<distance_PS)
				{
					i_PS=i;
					distance_PS=distance_PS_tmp;
				}												

			}

			int ij1k,ij0k,i1jk;
			int iPlus1jk,iMinus1jk;
			float vorzeichen;
			
			j=0;

			for(k=0;k<kMeshPoints;k++)
			{
				i=i_LE;
				ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
				ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;								
				cp_LE+=(p_transient[ijk]-1.)/(0.5*1.4*MachNumber*MachNumber)/kMeshPoints;
				cf_LE+=calcCF(iMeshPoints, iMeshPoints, iMeshPoints, i, 0, k, u_transient, v_transient, w_transient, p_transient, rho_transient, x, y, ReynoldsNumber)/kMeshPoints;
				if((k==kMeshPoints-1)&&(index_solution==NumberSolutions))printf("LE: x=%g y=%g\n",x[ijk],y[ijk]);
				
				i=i_TE;
				ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
				ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;								
				cp_TE+=(p_transient[ijk]-1.)/(0.5*1.4*MachNumber*MachNumber)/kMeshPoints;
				cf_TE+=calcCF(iMeshPoints, iMeshPoints, iMeshPoints, i, 0, k, u_transient, v_transient, w_transient, p_transient, rho_transient, x, y, ReynoldsNumber)/kMeshPoints;
				if((k==kMeshPoints-1)&&(index_solution==NumberSolutions))printf("TE: x=%g y=%g\n",x[ijk],y[ijk]);
				
				i=i_SS;
				ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
				ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;								
				cp_SS+=(p_transient[ijk]-1.)/(0.5*1.4*MachNumber*MachNumber)/kMeshPoints;
				cf_SS+=calcCF(iMeshPoints, iMeshPoints, iMeshPoints, i, 0, k, u_transient, v_transient, w_transient, p_transient, rho_transient, x, y, ReynoldsNumber)/kMeshPoints;
				if((k==kMeshPoints-1)&&(index_solution==NumberSolutions))printf("SS: x=%g y=%g\n",x[ijk],y[ijk]);
				
				i=i_PS;
				ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
				ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;								
				cp_PS+=(p_transient[ijk]-1.)/(0.5*1.4*MachNumber*MachNumber)/kMeshPoints;
				cf_PS+=calcCF(iMeshPoints, iMeshPoints, iMeshPoints, i, 0, k, u_transient, v_transient, w_transient, p_transient, rho_transient, x, y, ReynoldsNumber)/kMeshPoints;
				if((k==kMeshPoints-1)&&(index_solution==NumberSolutions))printf("PS: x=%g y=%g\n",x[ijk],y[ijk]);
			}
				
				
			if(index_solution==1)
			{
				char output[300];
				sprintf(zonename,"convergence");
				sprintf(output,"convergence.dat");
				fileOut=fopen(output,"w");
				fprintf(fileOut,"variables = \"t [s]\" \"c<sub>p,LE</sub> [-]\" \"c<sub>p,TE</sub> [-]\" \"c<sub>p,SS</sub> [-]\" \"c<sub>p,PS</sub> [-]\" \"c<sub>f,LE</sub> [-]\" \"c<sub>f,TE</sub> [-]\" \"c<sub>f,SS</sub> [-]\" \"c<sub>f,PS</sub> [-]\" \n");
				fprintf(fileOut,"zone t=\"%s\", i= %d, f=point \n",zonename,NumberSolutions);
			}


			fprintf(fileOut," %le %le %le %le %le %le %le %le %le\n",time[index_solution-1],cp_LE,cp_TE,cp_SS,cp_PS,cf_LE,cf_TE,cf_SS,cf_PS);

		}/////////////////////////////////ENDE convergence//////////////////////////////////////

			
	}///////////////////// ENDE SOLUTION-LOOP
	if(flag_convergence==1)
	{
		fclose(fileOut);
	}
	

	FILE * file0;
	char output2[300];
	sprintf(output2,"PressureHistory_Extracted_x%.1fy%.1f.dat",x_extr,y_extr);
	file0=fopen(output2,"w");
	fprintf(file0,"TITLE = \"PressureHistory extrahiert\"\n");
	fprintf(file0,"VARIABLES = \"t [s]\" \"PH_0_x%f_y%f_z%f\"\n",x_extr,y_extr,0.0);
	fprintf(file0,"ZONE T=\"PressureHistory_Extracted_x%.1fy%.1f.dat\", F=BLOCK, I=%d, DT=(SINGLE)\n",x_extr,y_extr,NumberSolutions);
	for(i=0;i<NumberSolutions;i++)
	{
		fprintf(file0,"%le\n",time[i]);
	}
	for(i=0;i<NumberSolutions;i++)
	{
		fprintf(file0,"%le\n",p_extracted[i]);
	}	
	fclose(file0);	


	printf("...fertig. Schließe Datei.\n");
	CG( cg_close(index_file));	
	


	return 0;
}///////////////////// ENDE Main
	


