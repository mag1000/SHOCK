#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"
#include "unistd.h"
#include <stdbool.h>

#define CG(cmd) if(cmd)cg_error_exit( );

double CalcLambda2();
double CalcQ();
double CalcPFluctuation();
double CalcMa();
double CalcGradRho();
double CalcGradRho2D();
long long int iMeshPoints,jMeshPoints,kMeshPoints;
long long int buffer;
float *buffer_tmp;
float *CoordinateX,*CoordinateY,*CoordinateZ,*VelocityX,*VelocityY,*VelocityZ,*Pressure,*Density;
float *Pressure_timeAverage; 
long long int i,j,k,ijk,ijk_buffer;
long long int i_help,j_help,k_help;
char varname[32];	
float mach;

char *replace_str(char *str, char *orig, char *rep)
{
  static char buffer[4096];
  char *p;

  if(!(p = strstr(str, orig)))  // Is 'orig' even in 'str'?
    return str;

  strncpy(buffer, str, p-str); // Copy characters from 'str' start to 'orig' st$
  buffer[p-str] = '\0';

  sprintf(buffer+(p-str), "%s%s", rep, p+strlen(orig));

  return buffer;
}

int main(int argc, char *argv[])
{
	char file_in[500];
	char file_out[500];	
	

	int inputflag=0;
	int c;

	int Solution_Start=1;
	int Solution_End=0;	
	int mflag;	
	int reducedflag=0;
	int endflag=0;
	int startflag=0;

	opterr = 0;

	while ((c = getopt (argc, argv, "i:s:e:m:r")) != -1)
	{
		switch (c)
		{
			case 'i':
				inputflag=1;
				strcpy(file_in,optarg);
				break;
			case 's':
				Solution_Start=atoi(optarg);
				startflag=1;
				break;				
			case 'e':
				Solution_End=atoi(optarg);
				endflag=1;				
				break;								
			case 'm':
				mflag=1;
				mach = atof(optarg);
				break;	
			case 'r':
				reducedflag=1;
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
	
	printf("Starting cgnsVisualizeData....\n");	
	printf("Die Output Datei ist eine cgns-Datei, die automatisch in eine szplt-Datei konvertiert wird.\n");
	
	printf("Folgende Optionen sind verfügbar: \n i: Inputfile (gefordert ) \n m: Mach (gefordert) \n s: start at solution (optional) \n e: end at solution (optional) \n r: reduced (optional)\n");	
	
	if(inputflag!=1)
	{
		printf("ERROR: Kein Inputfile angegeben (-i 'filename')\n");
		return 1;
	}
	
	if(mflag!=1)
	{
		printf("ERROR: Machzahl muss angegeben werden! (-m MACH)\n");
		return 1;
	}
	else
	{
		printf("Machzahl: %f\n",mach);
	}		
	
	if(startflag!=1)
	{
		printf("Keine start_solution angegeben (-s START, 0=letzte). Die erste solution wird verwendet.\n");
		Solution_Start=1;
	}
	
	if(endflag!=1)
	{
		printf("Keine end_solution angegeben (-e ENDE, 0=letzte). Die letzte solution wird verwendet.\n");
		Solution_End=0;
	}	

	if(reducedflag==1)
	{
		printf("Reduzierter Export (Nur Ma, Q und DensityGradient)\n");
	}


	printf("Öffne Datei %s im read mode\n",file_in);
	int index_file_in;
	int index_file_out;
	int index_base=1;
	CG( cg_open(file_in,CG_MODE_READ,&index_file_in));
	int cell_dim;
	CG( cg_cell_dim(index_file_in,index_base,&cell_dim));
	char solname[32];
	cgsize_t zonesize[cell_dim][cell_dim];	

	int number_solutions=1;
	char zonename[33];
	int nzones;
	CG( cg_nzones(index_file_in,index_base, &nzones));
	int index_zone;
	
	int NumberSolutions;
	CG( cg_nsols(index_file_in,1,1,&NumberSolutions ) );
	if(Solution_End==0){Solution_End=NumberSolutions;}
	float time[NumberSolutions];
	printf("Reading time-array 1...\n");
	CG( cg_goto(index_file_in,1,"BaseIterativeData_t",1,"end"));
	CG( cg_array_read_as(1,RealSingle,&time));	
	
	strcpy(file_out,replace_str(file_in, ".cgns", "_Ma-DensityGradient-L2-Q.cgns"));	
	printf("Öffne Datei %s im write mode\n",file_out);

	CG( cg_open(file_out,CG_MODE_WRITE,&index_file_out));
	
	printf("Writing base...\n");
	CG( cg_base_write(index_file_out,"Base",cell_dim,cell_dim,&index_base));	
	
	index_zone=1;
	CG( cg_zone_read(index_file_in,index_base,index_zone,zonename,zonesize[0]));

	cgsize_t irmin[cell_dim];
	cgsize_t irmax[cell_dim];
	
	if (cell_dim==2)
	{
		irmin[0]=1;
		irmin[1]=1;	
		irmax[0]=zonesize[0][0];
		irmax[1]=zonesize[0][1];
		printf("Lade Ergebnisse aus Zone %s (%zux%zu)...\n",zonename,irmax[0],irmax[1]);
		iMeshPoints=irmax[0];
		jMeshPoints=irmax[1];
		kMeshPoints=0;
	}
	else
	{
		irmin[0]=1;
		irmin[1]=1;	
		irmin[2]=1;	
		irmax[0]=zonesize[0][0];
		irmax[1]=zonesize[0][1];
		irmax[2]=zonesize[0][2];
		printf("Lade Ergebnisse aus Zone %s (%zux%zux%zu)...\n",zonename,irmax[0],irmax[1],irmax[2]);
		iMeshPoints=irmax[0];
		jMeshPoints=irmax[1];
		kMeshPoints=irmax[2];
	}
	
	printf("Writing zone...\n");
	CG( cg_zone_write(index_file_out,index_base,zonename,*zonesize,Structured,&index_zone));

	if(cell_dim==2)
	{
		buffer=irmax[0]*irmax[1];
	}
	else
	{
		buffer=irmax[0]*irmax[1]*irmax[2];
	}
	
	CoordinateX=(float *)calloc(buffer, sizeof(float));
	CoordinateY=(float *)calloc(buffer, sizeof(float));
	CoordinateZ=(float *)calloc(buffer, sizeof(float));
	VelocityX=(float *)calloc(buffer, sizeof(float));
	VelocityY=(float *)calloc(buffer, sizeof(float));
	VelocityZ=(float *)calloc(buffer, sizeof(float));
	Pressure=(float *)calloc(buffer, sizeof(float));
	Density=(float *)calloc(buffer, sizeof(float));	
	buffer_tmp=(float*)calloc(buffer,sizeof(float));


	//////////////////////////////////////////////////////////////////////////
	printf("Copy coordinates for zone:%d\n",index_zone);
	
	int index_coord;
	printf("Transfering CoordinateX ...\n");
	cg_coord_read(index_file_in,1,index_zone,"CoordinateX",RealSingle,irmin,irmax,CoordinateX);
	cg_coord_write(index_file_out,1,index_zone,RealSingle,"CoordinateX",CoordinateX,&index_coord);
	
	printf("Transfering CoordinateY ...\n");
	cg_coord_read(index_file_in,1,index_zone,"CoordinateY",RealSingle,irmin,irmax,CoordinateY);
	cg_coord_write(index_file_out,1,index_zone,RealSingle,"CoordinateY",CoordinateY,&index_coord);
	
	if(cell_dim==3)
	{
		printf("Transfering CoordinateZ ...\n");
		cg_coord_read(index_file_in,1,index_zone,"CoordinateZ",RealSingle,irmin,irmax,CoordinateZ);
		cg_coord_write(index_file_out,1,index_zone,RealSingle,"CoordinateZ",CoordinateZ,&index_coord);	
	}
	//////////////////////////////////////////////////////////////////////////
	
	
	///////////////////// Lade timeAverage-Ergebnisse
	if(reducedflag==0)
	{
	printf("Lese timeAverage Werte für die Erstellung von PFluctuation (Instationäre Druckschwankung).\n");
	char file_timeAverage[500];
	int index_file_timeAverage;
	strcpy(file_timeAverage,replace_str(file_in, ".cgns","_VarianzTimeAverage.cgns"));
	CG( cg_open(file_timeAverage,CG_MODE_READ,&index_file_timeAverage));
	
	
	CG( cg_zone_read(index_file_timeAverage,index_base,index_zone,zonename,zonesize[0]));
	
	irmin[0]=1;
	irmin[1]=1;	
	irmax[0]=zonesize[0][0]; 
	irmax[1]=zonesize[0][1];
	
	iMeshPoints=irmax[0]-irmin[0]+1;
	jMeshPoints=irmax[1]-irmin[1]+1;
	
	if(cell_dim==3)
	{
		irmin[2]=1;
		irmax[2]=zonesize[0][2];
		kMeshPoints=irmax[2]-irmin[2]+1;
	}
	else
	{
		kMeshPoints=1;
	}
	
	Pressure_timeAverage=(float *)calloc(buffer, sizeof(float));
	
	number_solutions=1;	
	CG( cg_field_read(index_file_timeAverage,index_base,index_zone,number_solutions,"Pressure_timeAverage",
	RealSingle,irmin,irmax,Pressure_timeAverage));
	printf("fertig.\n");
	}
	///////////////////// timeAverage fertig			
	

	//////////////////////////////////////////////////////////////////////////
	printf("Copy solutions for zone:%d\n",index_zone);
	char *solnames;
	solnames = malloc( 32* (Solution_End-Solution_Start+1)  );
	int index_solution;
	int index_solution_out;	
	int index_field;
	GridLocation_t location; //Vertex, CellCenter, IFaceCenter, JFaceCenter, and KFaceCenter. 
	printf("Kopiere Solution:%d bis %d\n",Solution_Start,Solution_End);	
	for( index_solution = Solution_Start; index_solution<=Solution_End; index_solution++ ) 
	{
	printf("Solution:%d\n",index_solution);	
	
	printf("Lade Ergebnisse\n");
	CG( cg_field_read(index_file_in,index_base,index_zone,index_solution,"VelocityX",RealSingle,irmin,irmax,VelocityX));
	CG( cg_field_read(index_file_in,index_base,index_zone,index_solution,"VelocityY",RealSingle,irmin,irmax,VelocityY));
	CG( cg_field_read(index_file_in,index_base,index_zone,index_solution,"VelocityZ",RealSingle,irmin,irmax,VelocityZ));	
	CG( cg_field_read(index_file_in,index_base,index_zone,index_solution,"Pressure",RealSingle,irmin,irmax,Pressure));
	CG( cg_field_read(index_file_in,index_base,index_zone,index_solution,"Density",RealSingle,irmin,irmax,Density));

	CG( cg_sol_info(index_file_in, 1, index_zone, index_solution, solname, &location));		
	CG( cg_sol_write( index_file_out,1,index_zone,solname,location,&index_solution_out) );		

	memcpy( solnames+(index_solution-Solution_Start)*32,solname,32 );
	
	transferMa();
	cg_field_write(index_file_out,1,index_zone,index_solution_out,RealSingle,"M",buffer_tmp,&index_field);	
	
	if(cell_dim==3)
	{	
		transferGradRho();
		cg_field_write(index_file_out,1,index_zone,index_solution_out,RealSingle,"DensityGradient",buffer_tmp,&index_field);
	}
	if(cell_dim==2)
	{
		transferGradRho2D();
		cg_field_write(index_file_out,1,index_zone,index_solution_out,RealSingle,"DensityGradient",buffer_tmp,&index_field);		
	}				
	
	if((reducedflag==0)&&(cell_dim==3))
	{
		transferL2();
		cg_field_write(index_file_out,1,index_zone,index_solution_out,RealSingle,"Lambda2",buffer_tmp,&index_field);		
	}
	
	if(cell_dim==3)
	{
		transferQ();
		cg_field_write(index_file_out,1,index_zone,index_solution_out,RealSingle,"Q",buffer_tmp,&index_field);		
	}
	
	if(reducedflag==0)
	{
		transferPFluctuation();
		cg_field_write(index_file_out,1,index_zone,index_solution_out,RealSingle,"p'",buffer_tmp,&index_field);		
	}		

	}
	//////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////
	printf("Copy zoneiter for zone:%d\n",index_zone);
	CG( cg_ziter_write( index_file_out,1,1,"ZoneIterativeData"));
	CG( cg_goto(index_file_out,1,"Zone_t",1,"ZoneIterativeData_t",1,"end"));
	CG( cg_array_write( "FlowSolutionPointers",Character,2,(cgsize_t[ ]){ 32,(Solution_End-Solution_Start+1) },solnames));	
	//////////////////////////////////////////////////////////////////////////	
	
	
	
	//////////////////////////////////////////////////////////////////////////
	printf("Writing new baseiter\n");
	CG( cg_biter_write(index_file_out,1,"BaseIterativeData",Solution_End-Solution_Start+1));
	CG( cg_goto(index_file_out,1,"BaseIterativeData_t",1,"end"));
	float time_out[Solution_End-Solution_Start+1];
	
	int t;
	for(t=Solution_Start;t<Solution_End;t++)
	{
		time_out[t-Solution_Start]=time[t-1];
	}
	cgsize_t nuse=Solution_End-Solution_Start+1;
	CG( cg_array_write("TimeValues",RealSingle,1,&nuse,&time_out));
	//////////////////////////////////////////////////////////////////////////
	
	
	CG (cg_close(index_file_out));
	CG (cg_close(index_file_in));		
	
	
	printf("Erzeuge szplt-Datei\n");
	char command[500];

	sprintf(command,"cp /home/mag/TOOLS/TECPLOT/convert2szplt.* ./");
	system(command);	
	
	strcpy(file_out,replace_str(file_out, ".cgns", ".szplt"));	
	sprintf(command,"sed -i s/SAMPLEFILENAME/%s/g convert2szplt.mcr",file_out);
	system(command);
	
	strcpy(file_out,replace_str(file_out, ".szplt", ".cgns"));	
	sprintf(command,"sed -i s/SAMPLEFILENAME/%s/g convert2szplt*",file_out);
	system(command);	
	
	char Varlist[500];
	int var,varMax=3;
	strcpy(Varlist,"");
	if ((reducedflag==1) && (cell_dim==2))
	{
		varMax=2;
	}
	if ((reducedflag==0) && (cell_dim==3))
	{
		varMax=5;
	}
	for( var = 0; var<varMax; var++ )
	{
		sprintf(Varlist,"%s%d,",Varlist,var);
	}	
	sprintf(command,"sed -i s/SAMPLEVARLIST/%s/g convert2szplt*",Varlist);
	system(command);
	
	char Zonelist[500];
	strcpy(Zonelist,"");
	for( index_solution = Solution_Start; index_solution<=Solution_End; index_solution++ )
	{
		sprintf(Zonelist,"%s%d,",Zonelist,index_solution-Solution_Start);
	}
	sprintf(command,"sed -i s/SAMPLEZONELIST/%s/g convert2szplt*",Zonelist);
	system(command);	
		
	sprintf(command,"/usr/local/Tecplot2014R2/bin/tec360 -b -p convert2szplt.mcr convert2szplt.lay");	
	printf("Executing:\n%s\n",command);
	system(command);
	
	//Aufräumarbeiten
	/*printf("Aufräumarbeiten\n");
	strcpy(file_out,replace_str(file_out, ".szplt", ".cgns"));	
	sprintf(command,"rm -fr %s convert2szplt* batch.log",file_out);
	printf("Executing:\n%s\n",command);
	system(command);*/
	
	
	return 0;
}

int transferL2()
{
	printf("Uebertrage Lamda2...\n");
	for(k_help=0;k_help<kMeshPoints;k_help++)
	{
		k=k_help;
		if (k_help==0){k=1;}
		if (k_help==kMeshPoints-1){k=k_help-1;}
		for(j_help=0;j_help<jMeshPoints;j_help++)
		{
			j=j_help;
			if (j_help==0){j=1;}
			if (j_help==jMeshPoints-1){j=j_help-1;}
			for(i_help=0;i_help<iMeshPoints;i_help++)
			{
				i=i_help;
				if (i_help==0){i=1;}
				if (i_help==iMeshPoints-1){i=i_help-1;}
				ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ijk_buffer=k_help*jMeshPoints*iMeshPoints+j_help*iMeshPoints+i_help;
				buffer_tmp[ijk_buffer]=CalcLambda2();
			}
		}
	}
}

int transferQ()
{
	printf("Uebertrage Q...\n");
	for(k_help=0;k_help<kMeshPoints;k_help++)
	{
		k=k_help;
		if (k_help==0){k=1;}
		if (k_help==kMeshPoints-1){k=k_help-1;}
		for(j_help=0;j_help<jMeshPoints;j_help++)
		{
			j=j_help;
			if (j_help==0){j=1;}
			if (j_help==jMeshPoints-1){j=j_help-1;}
			for(i_help=0;i_help<iMeshPoints;i_help++)
			{
				i=i_help;
				if (i_help==0){i=1;}
				if (i_help==iMeshPoints-1){i=i_help-1;}
				ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ijk_buffer=k_help*jMeshPoints*iMeshPoints+j_help*iMeshPoints+i_help;
				buffer_tmp[ijk_buffer]=CalcQ();
			}
		}
	}
}

int transferMa()
{
	printf("Uebertrage M\n");
	for(k_help=0;k_help<kMeshPoints;k_help++)
	{
		k=k_help;
		for(j_help=0;j_help<jMeshPoints;j_help++)
		{
			j=j_help;
			for(i_help=0;i_help<iMeshPoints;i_help++)
			{
				i=i_help;
				ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ijk_buffer=k_help*jMeshPoints*iMeshPoints+j_help*iMeshPoints+i_help;
				buffer_tmp[ijk_buffer]=CalcMa();
			}
		}
	}
}

int transferPFluctuation()
{
	printf("Uebertrage PFluctuation...\n");
	for(k_help=0;k_help<kMeshPoints;k_help++)
	{
		k=k_help;
		for(j_help=0;j_help<jMeshPoints;j_help++)
		{
			j=j_help;
			for(i_help=0;i_help<iMeshPoints;i_help++)
			{
				i=i_help;
				ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ijk_buffer=k_help*jMeshPoints*iMeshPoints+j_help*iMeshPoints+i_help;
				buffer_tmp[ijk_buffer]=CalcPFluctuation();				
			}
		}
	}
}

int transferGradRho()
{
	printf("Uebertrage GradRho...\n");
	for(k_help=0;k_help<kMeshPoints;k_help++)
	{
		k=k_help;
		if (k_help==0){k=1;}
		if (k_help==kMeshPoints-1){k=k_help-1;}
		for(j_help=0;j_help<jMeshPoints;j_help++)
		{
			j=j_help;
			if (j_help==0){j=1;}
			if (j_help==jMeshPoints-1){j=j_help-1;}
			for(i_help=0;i_help<iMeshPoints;i_help++)
			{
				i=i_help;
				if (i_help==0){i=1;}
				if (i_help==iMeshPoints-1){i=i_help-1;}
				ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
				ijk_buffer=k_help*jMeshPoints*iMeshPoints+j_help*iMeshPoints+i_help;
				buffer_tmp[ijk_buffer]=CalcGradRho();	
			}
		}
	}
}

int transferGradRho2D()
{
	printf("Uebertrage GradRho...\n");
	k=0;
	for(j_help=0;j_help<jMeshPoints;j_help++)
	{
		j=j_help;
		if (j_help==0){j=1;}
		if (j_help==jMeshPoints-1){j=j_help-1;}
		for(i_help=0;i_help<iMeshPoints;i_help++)
		{
			i=i_help;
			if (i_help==0){i=1;}
			if (i_help==iMeshPoints-1){i=i_help-1;}
			ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
			ijk_buffer=k_help*jMeshPoints*iMeshPoints+j_help*iMeshPoints+i_help;
			buffer_tmp[ijk_buffer]=CalcGradRho2D();	
		}
	}
}

double CalcMa()
{
	double Ma;
	Ma=(sqrt(VelocityX[ijk]*VelocityX[ijk]+VelocityY[ijk]*VelocityY[ijk]+VelocityZ[ijk]*VelocityZ[ijk]))/sqrt(Pressure[ijk]/Density[ijk])*mach;
	return Ma;
}

double CalcPFluctuation()
{
	double PFluctuation;
	PFluctuation=Pressure[ijk]-Pressure_timeAverage[ijk];
	return PFluctuation;
}

double CalcGradRho()
{
	long long int iPlus1jk,iMinus1jk;
	long long int ijPlus1k,ijMinus1k;
	long long int ijkPlus1,ijkMinus1;
	
	double jacobian;
	double xi_x,eta_x,zeta_x;
	double xi_y,eta_y,zeta_y;
	double xi_z,eta_z,zeta_z;
	double x_xi,x_eta,x_zeta;
	double y_xi,y_eta,y_zeta;
	double z_xi,z_eta,z_zeta;	
	
	double GradRho;

	iPlus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i+1;
	iMinus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;
	ijPlus1k=k*jMeshPoints*iMeshPoints+(j+1)*iMeshPoints+i;
	ijMinus1k=k*jMeshPoints*iMeshPoints+(j-1)*iMeshPoints+i;
	ijkPlus1=(k+1)*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
	ijkMinus1=(k-1)*jMeshPoints*iMeshPoints+j*iMeshPoints+i;


	x_xi=(-0.5*CoordinateX[iMinus1jk]+0.5*CoordinateX[iPlus1jk]);
	y_xi=(-0.5*CoordinateY[iMinus1jk]+0.5*CoordinateY[iPlus1jk]);
	z_xi=(-0.5*CoordinateZ[iMinus1jk]+0.5*CoordinateZ[iPlus1jk]);
	x_eta=(-0.5*CoordinateX[ijMinus1k]+0.5*CoordinateX[ijPlus1k]);
	y_eta=(-0.5*CoordinateY[ijMinus1k]+0.5*CoordinateY[ijPlus1k]);
	z_eta=(-0.5*CoordinateZ[ijMinus1k]+0.5*CoordinateZ[ijPlus1k]);
	x_zeta=(-0.5*CoordinateX[ijkMinus1]+0.5*CoordinateX[ijkPlus1]);
	y_zeta=(-0.5*CoordinateY[ijkMinus1]+0.5*CoordinateY[ijkPlus1]);
	z_zeta=(-0.5*CoordinateZ[ijkMinus1]+0.5*CoordinateZ[ijkPlus1]);


	jacobian=x_xi*y_eta*z_zeta+x_eta*y_zeta*z_xi+x_zeta*y_xi*z_eta
		-x_zeta*y_eta*z_xi-x_eta*y_xi*z_zeta-x_xi*y_zeta*z_eta;


	xi_x=(y_eta*z_zeta-y_zeta*z_eta)/jacobian;
	xi_y=(x_zeta*z_eta-x_eta*z_zeta)/jacobian;
	xi_z=(x_eta*y_zeta-x_zeta*y_eta)/jacobian;
	eta_x=(y_zeta*z_xi-y_xi*z_zeta)/jacobian;
	eta_y=(x_xi*z_zeta-x_zeta*z_xi)/jacobian;
	eta_z=(x_zeta*y_xi-x_xi*y_zeta)/jacobian;
	zeta_x=(y_xi*z_eta-y_eta*z_xi)/jacobian;
	zeta_y=(x_eta*z_xi-x_xi*z_eta)/jacobian;
	zeta_z=(x_xi*y_eta-x_eta*y_xi)/jacobian;
	
	GradRho=
		sqrt(
		pow((0.5*Density[iPlus1jk]-0.5*Density[iMinus1jk])*xi_x,2.0)+
		pow((0.5*Density[ijPlus1k]-0.5*Density[ijMinus1k])*eta_x,2.0)+
		pow((0.5*Density[ijkPlus1]-0.5*Density[ijkMinus1])*zeta_x,2.0)+
		pow((0.5*Density[iPlus1jk]-0.5*Density[iMinus1jk])*xi_y,2.0)+
		pow((0.5*Density[ijPlus1k]-0.5*Density[ijMinus1k])*eta_y,2.0)+
		pow((0.5*Density[ijkPlus1]-0.5*Density[ijkMinus1])*zeta_y,2.0)+
		pow((0.5*Density[iPlus1jk]-0.5*Density[iMinus1jk])*xi_z,2.0)+
		pow((0.5*Density[ijPlus1k]-0.5*Density[ijMinus1k])*eta_z,2.0)+
		pow((0.5*Density[ijkPlus1]-0.5*Density[ijkMinus1])*zeta_z,2.0)
		);
	return GradRho;
}

double CalcGradRho2D()
{
	long long int iPlus1jk,iMinus1jk;
	long long int ijPlus1k,ijMinus1k;
	
	double jacobian;
	double xi_x,eta_x,zeta_x;
	double xi_y,eta_y,zeta_y;
	double xi_z,eta_z,zeta_z;
	double x_xi,x_eta,x_zeta;
	double y_xi,y_eta,y_zeta;
	double z_xi,z_eta,z_zeta;	
	
	double GradRho;

	iPlus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i+1;
	iMinus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;
	ijPlus1k=k*jMeshPoints*iMeshPoints+(j+1)*iMeshPoints+i;
	ijMinus1k=k*jMeshPoints*iMeshPoints+(j-1)*iMeshPoints+i;


	x_xi=(-0.5*CoordinateX[iMinus1jk]+0.5*CoordinateX[iPlus1jk]);
	y_xi=(-0.5*CoordinateY[iMinus1jk]+0.5*CoordinateY[iPlus1jk]);
	z_xi=0.0;
	x_eta=(-0.5*CoordinateX[ijMinus1k]+0.5*CoordinateX[ijPlus1k]);
	y_eta=(-0.5*CoordinateY[ijMinus1k]+0.5*CoordinateY[ijPlus1k]);
	z_eta=0.0;
	x_zeta=0.0;
	y_zeta=0.0;
	z_zeta=1.0;


	jacobian=x_xi*y_eta*z_zeta+x_eta*y_zeta*z_xi+x_zeta*y_xi*z_eta
		-x_zeta*y_eta*z_xi-x_eta*y_xi*z_zeta-x_xi*y_zeta*z_eta;


	xi_x=(y_eta*z_zeta-y_zeta*z_eta)/jacobian;
	xi_y=(x_zeta*z_eta-x_eta*z_zeta)/jacobian;
	xi_z=(x_eta*y_zeta-x_zeta*y_eta)/jacobian;
	eta_x=(y_zeta*z_xi-y_xi*z_zeta)/jacobian;
	eta_y=(x_xi*z_zeta-x_zeta*z_xi)/jacobian;
	eta_z=(x_zeta*y_xi-x_xi*y_zeta)/jacobian;
	zeta_x=(y_xi*z_eta-y_eta*z_xi)/jacobian;
	zeta_y=(x_eta*z_xi-x_xi*z_eta)/jacobian;
	zeta_z=(x_xi*y_eta-x_eta*y_xi)/jacobian;
	
	GradRho=
		sqrt(
		pow((0.5*Density[iPlus1jk]-0.5*Density[iMinus1jk])*xi_x,2.0)+
		pow((0.5*Density[ijPlus1k]-0.5*Density[ijMinus1k])*eta_x,2.0)+
		pow((0.5*Density[iPlus1jk]-0.5*Density[iMinus1jk])*xi_y,2.0)+
		pow((0.5*Density[ijPlus1k]-0.5*Density[ijMinus1k])*eta_y,2.0)
		);
	return GradRho;
}


double CalcLambda2()
{
	long long int iPlus1jk,iMinus1jk;
	long long int ijPlus1k,ijMinus1k;
	long long int ijkPlus1,ijkMinus1;
	
	double jacobian;

	double p,q,rrr,al,yy1,yy2,yy3,lam1,lam2,lam3,lam,lamax,lamin;
	double xi_x,eta_x,zeta_x;
	double xi_y,eta_y,zeta_y;
	double xi_z,eta_z,zeta_z;
	double x_xi,x_eta,x_zeta;
	double y_xi,y_eta,y_zeta;
	double z_xi,z_eta,z_zeta;	

	double udx,udy,udz;
	double vdx,vdy,vdz;
	double wdx,wdy,wdz;

	double u_xi,u_eta,u_zeta;
	double v_xi,v_eta,v_zeta;
	double w_xi,w_eta,w_zeta;

	double om12,om13,om23,oo11,oo12,oo13,oo22,oo23,oo33;
	double s11,s12,s13,s22,s23,s33,ss11,ss12,ss13,ss22,ss23,ss33;
	double a11,a12,a13,a22,a23,a33,f1,f2,f3;

	lam =0.;

	/****************************************************************************/
	iPlus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i+1;
	iMinus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;
	ijPlus1k=k*jMeshPoints*iMeshPoints+(j+1)*iMeshPoints+i;
	ijMinus1k=k*jMeshPoints*iMeshPoints+(j-1)*iMeshPoints+i;
	ijkPlus1=(k+1)*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
	ijkMinus1=(k-1)*jMeshPoints*iMeshPoints+j*iMeshPoints+i;


	x_xi=(-0.5*CoordinateX[iMinus1jk]+0.5*CoordinateX[iPlus1jk]);
	y_xi=(-0.5*CoordinateY[iMinus1jk]+0.5*CoordinateY[iPlus1jk]);
	z_xi=(-0.5*CoordinateZ[iMinus1jk]+0.5*CoordinateZ[iPlus1jk]);
	x_eta=(-0.5*CoordinateX[ijMinus1k]+0.5*CoordinateX[ijPlus1k]);
	y_eta=(-0.5*CoordinateY[ijMinus1k]+0.5*CoordinateY[ijPlus1k]);
	z_eta=(-0.5*CoordinateZ[ijMinus1k]+0.5*CoordinateZ[ijPlus1k]);
	x_zeta=(-0.5*CoordinateX[ijkMinus1]+0.5*CoordinateX[ijkPlus1]);
	y_zeta=(-0.5*CoordinateY[ijkMinus1]+0.5*CoordinateY[ijkPlus1]);
	z_zeta=(-0.5*CoordinateZ[ijkMinus1]+0.5*CoordinateZ[ijkPlus1]);


	jacobian=x_xi*y_eta*z_zeta+x_eta*y_zeta*z_xi+x_zeta*y_xi*z_eta
		-x_zeta*y_eta*z_xi-x_eta*y_xi*z_zeta-x_xi*y_zeta*z_eta;


	xi_x=(y_eta*z_zeta-y_zeta*z_eta)/jacobian;
	xi_y=(x_zeta*z_eta-x_eta*z_zeta)/jacobian;
	xi_z=(x_eta*y_zeta-x_zeta*y_eta)/jacobian;
	eta_x=(y_zeta*z_xi-y_xi*z_zeta)/jacobian;
	eta_y=(x_xi*z_zeta-x_zeta*z_xi)/jacobian;
	eta_z=(x_zeta*y_xi-x_xi*y_zeta)/jacobian;
	zeta_x=(y_xi*z_eta-y_eta*z_xi)/jacobian;
	zeta_y=(x_eta*z_xi-x_xi*z_eta)/jacobian;
	zeta_z=(x_xi*y_eta-x_eta*y_xi)/jacobian;

	u_xi=(-0.5*VelocityX[iMinus1jk]+0.5*VelocityX[iPlus1jk]);
	v_xi=(-0.5*VelocityY[iMinus1jk]+0.5*VelocityY[iPlus1jk]);
	w_xi=(-0.5*VelocityZ[iMinus1jk]+0.5*VelocityZ[iPlus1jk]);
	u_eta=(-0.5*VelocityX[ijMinus1k]+0.5*VelocityX[ijPlus1k]);
	v_eta=(-0.5*VelocityY[ijMinus1k]+0.5*VelocityY[ijPlus1k]);
	w_eta=(-0.5*VelocityZ[ijMinus1k]+0.5*VelocityZ[ijPlus1k]);
	u_zeta=(-0.5*VelocityX[ijkMinus1]+0.5*VelocityX[ijkPlus1]);
	v_zeta=(-0.5*VelocityY[ijkMinus1]+0.5*VelocityY[ijkPlus1]);
	w_zeta=(-0.5*VelocityZ[ijkMinus1]+0.5*VelocityZ[ijkPlus1]);

	udx=u_xi*xi_x+u_eta*eta_x+u_zeta*zeta_x;
	udy=u_xi*xi_y+u_eta*eta_y+u_zeta*zeta_y;
	udz=u_xi*xi_z+u_eta*eta_z+u_zeta*zeta_z;
	vdx=v_xi*xi_x+v_eta*eta_x+v_zeta*zeta_x;
	vdy=v_xi*xi_y+v_eta*eta_y+v_zeta*zeta_y;
	vdz=v_xi*xi_z+v_eta*eta_z+v_zeta*zeta_z;
	wdx=w_xi*xi_x+w_eta*eta_x+w_zeta*zeta_x;
	wdy=w_xi*xi_y+w_eta*eta_y+w_zeta*zeta_y;
	wdz=w_xi*xi_z+w_eta*eta_z+w_zeta*zeta_z;

	/****************************************************************************/

	om12 = 0.5 * (udy-vdx);
	om13 = 0.5 * (udz-wdx);
	om23 = 0.5 * (vdz-wdy);

	/****************************************************************************/

	s11 = 0.5 * (udx);
	s12 = 0.5 * (udy+vdx);
	s13 = 0.5 * (udz+wdx);
	s22 = 0.5 * (vdy);
	s23 = 0.5 * (vdz+wdy);
	s33 = 0.5 * (wdz);

	/*****************************************************************************/

	oo11 =-om12*om12-om13*om13;
	oo12 =-om13*om23;
	oo13 =om12*om23;
	oo22 =-om12*om12-om23*om23;
	oo23 =-om12*om13;
	oo33 =-om13*om13-om23*om23;

	/*****************************************************************************/

	ss11 =s11*s11+s12*s12+s13*s13;
	ss12 =s11*s12+s12*s22+s13*s23;
	ss13 =s11*s13+s12*s23+s13*s33;
	ss22 =s12*s12+s22*s22+s23*s23;
	ss23 =s12*s13+s22*s23+s23*s33;
	ss33 =s13*s13+s23*s23+s33*s33;

	/******************************************************************************/

	a11 = ss11+oo11;
	a12 = ss12+oo12;
	a13 = ss13+oo13;
	a22 = ss22+oo22;
	a23 = ss23+oo23;
	a33 = ss33+oo33;

	/******************************************************************************/

	f1 = -a11-a22-a33;
	f2 = -a23*a23-a12*a12-a13*a13+a11*a33+a11*a22+a22*a33;
	f3 = -a11*a22*a33+a11*a23*a23+a22*a13*a13+a33*a12*a12-2.*a12*a13*a23;

	/******************************************************************************/

	p = -f1*f1/3.+f2;    q = 2.*f1*f1*f1/27.-f1*f2/3.+f3;

	if(p < 0.){rrr=-p;}
	else      {rrr=p; }
	rrr = sqrt(rrr/3.);

	al = acos(q/(2.*rrr*rrr*rrr));

	yy1 =  -2.*rrr*cos(al/3.);
	yy2 = 2.*rrr*cos(M_PI/3.-al/3.);
	yy3 = 2.*rrr*cos(M_PI/3.+al/3.);

	/**********************************************************************/
	lam1 = yy1 - f1/3;    lam2 = yy2 - f1/3;    lam3 = yy3 - f1/3;
	/*****************************************************/
	/*****************************************************/
	/****************************************************/

	if(lam2>=lam1)
	{
		lamax=lam2;
	}
	else
	{
		lamax=lam1;
	}
	if(lam3>=lamax)
	{
		lamax=lam3;
	}

	if(lam2<=lam1)
	{
		lamin=lam2;
	}
	else
	{
		lamin=lam1;
	}
	if(lam3<=lamin)
	{
		lamin=lam3;
	}

	if(lam1==lam2 || lam2==lam3)
	{
		lam = lam2;
	}

	if(lam1>lamin && lam1<lamax)
	{
		lam = lam1;
	}
	else
	{
		if(lam2>lamin && lam2<lamax)
		{
			lam = lam2;
		}
		else
		{
			if(lam3>lamin && lam3<lamax)
			{
				lam = lam3;
			}
		}
	}

	if(lam >= 0.)
	{
		lam =0.;
	}

	return lam;
}


double CalcQ()
{
	long long int iPlus1jk,iMinus1jk;
	long long int ijPlus1k,ijMinus1k;
	long long int ijkPlus1,ijkMinus1;
	
	double q;

	double xi_x,eta_x,zeta_x;
	double xi_y,eta_y,zeta_y;
	double xi_z,eta_z,zeta_z;
	double x_xi,x_eta,x_zeta;
	double y_xi,y_eta,y_zeta;
	double z_xi,z_eta,z_zeta;
	double jacobian;	

	double udx,udy,udz;
	double vdx,vdy,vdz;
	double wdx,wdy,wdz;

	double u_xi,u_eta,u_zeta;
	double v_xi,v_eta,v_zeta;
	double w_xi,w_eta,w_zeta;


	/****************************************************************************/
	iPlus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i+1;
	iMinus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;
	ijPlus1k=k*jMeshPoints*iMeshPoints+(j+1)*iMeshPoints+i;
	ijMinus1k=k*jMeshPoints*iMeshPoints+(j-1)*iMeshPoints+i;
	ijkPlus1=(k+1)*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
	ijkMinus1=(k-1)*jMeshPoints*iMeshPoints+j*iMeshPoints+i;


	x_xi=(-0.5*CoordinateX[iMinus1jk]+0.5*CoordinateX[iPlus1jk]);
	y_xi=(-0.5*CoordinateY[iMinus1jk]+0.5*CoordinateY[iPlus1jk]);
	z_xi=(-0.5*CoordinateZ[iMinus1jk]+0.5*CoordinateZ[iPlus1jk]);
	x_eta=(-0.5*CoordinateX[ijMinus1k]+0.5*CoordinateX[ijPlus1k]);
	y_eta=(-0.5*CoordinateY[ijMinus1k]+0.5*CoordinateY[ijPlus1k]);
	z_eta=(-0.5*CoordinateZ[ijMinus1k]+0.5*CoordinateZ[ijPlus1k]);
	x_zeta=(-0.5*CoordinateX[ijkMinus1]+0.5*CoordinateX[ijkPlus1]);
	y_zeta=(-0.5*CoordinateY[ijkMinus1]+0.5*CoordinateY[ijkPlus1]);
	z_zeta=(-0.5*CoordinateZ[ijkMinus1]+0.5*CoordinateZ[ijkPlus1]);


	jacobian=x_xi*y_eta*z_zeta+x_eta*y_zeta*z_xi+x_zeta*y_xi*z_eta
		-x_zeta*y_eta*z_xi-x_eta*y_xi*z_zeta-x_xi*y_zeta*z_eta;


	xi_x=(y_eta*z_zeta-y_zeta*z_eta)/jacobian;
	xi_y=(x_zeta*z_eta-x_eta*z_zeta)/jacobian;
	xi_z=(x_eta*y_zeta-x_zeta*y_eta)/jacobian;
	eta_x=(y_zeta*z_xi-y_xi*z_zeta)/jacobian;
	eta_y=(x_xi*z_zeta-x_zeta*z_xi)/jacobian;
	eta_z=(x_zeta*y_xi-x_xi*y_zeta)/jacobian;
	zeta_x=(y_xi*z_eta-y_eta*z_xi)/jacobian;
	zeta_y=(x_eta*z_xi-x_xi*z_eta)/jacobian;
	zeta_z=(x_xi*y_eta-x_eta*y_xi)/jacobian;

	u_xi=(-0.5*VelocityX[iMinus1jk]+0.5*VelocityX[iPlus1jk]);
	v_xi=(-0.5*VelocityY[iMinus1jk]+0.5*VelocityY[iPlus1jk]);
	w_xi=(-0.5*VelocityZ[iMinus1jk]+0.5*VelocityZ[iPlus1jk]);
	u_eta=(-0.5*VelocityX[ijMinus1k]+0.5*VelocityX[ijPlus1k]);
	v_eta=(-0.5*VelocityY[ijMinus1k]+0.5*VelocityY[ijPlus1k]);
	w_eta=(-0.5*VelocityZ[ijMinus1k]+0.5*VelocityZ[ijPlus1k]);
	u_zeta=(-0.5*VelocityX[ijkMinus1]+0.5*VelocityX[ijkPlus1]);
	v_zeta=(-0.5*VelocityY[ijkMinus1]+0.5*VelocityY[ijkPlus1]);
	w_zeta=(-0.5*VelocityZ[ijkMinus1]+0.5*VelocityZ[ijkPlus1]);

	udx=u_xi*xi_x+u_eta*eta_x+u_zeta*zeta_x;
	udy=u_xi*xi_y+u_eta*eta_y+u_zeta*zeta_y;
	udz=u_xi*xi_z+u_eta*eta_z+u_zeta*zeta_z;
	vdx=v_xi*xi_x+v_eta*eta_x+v_zeta*zeta_x;
	vdy=v_xi*xi_y+v_eta*eta_y+v_zeta*zeta_y;
	vdz=v_xi*xi_z+v_eta*eta_z+v_zeta*zeta_z;
	wdx=w_xi*xi_x+w_eta*eta_x+w_zeta*zeta_x;
	wdy=w_xi*xi_y+w_eta*eta_y+w_zeta*zeta_y;
	wdz=w_xi*xi_z+w_eta*eta_z+w_zeta*zeta_z;

	/*
	Diese Formel für Q wird in der Literatur ebenfalls angegeben ist jedoch veraltet.
	Eigene Vergleiche haben gezeigt, dass diese Variante Schmutz enthaelt, den die neuere Variante nicht hat.
	Daher wird die Variante unten stehend verwendet.
	q=udx*vdy+udx*wdz+vdy*wdz
	 -udy*vdx-udz*wdx-vdz*wdy;
	*/
	q=-0.5*(udx*udx+vdy*vdy+wdz*wdz)
	 -udy*vdx-udz*wdx-vdz*wdy;	 

	return q;
}
