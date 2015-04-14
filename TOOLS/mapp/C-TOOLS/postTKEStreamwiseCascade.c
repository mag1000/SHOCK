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
	char file_transient[500];
	char file_timeAverage3D[500];

	int fileflag=0;		
	int c;
	int optionflag=0;
	opterr = 0;
	char zonename[100];
	strcpy(zonename,"Zone");
	float Re_deltaStern_vorgabe[2];

	int time_samples, time_stepwidth;
	float ReynoldsNumber=500000.;
	float MachNumber=0.65;	
	while ((c = getopt (argc, argv, "m:f:g:x:y:s:w:r:z:")) != -1)
	{
		switch (c)
		{
			case 'm':
				MachNumber=atof(optarg);
				break;	
			case 'r':
				ReynoldsNumber=atof(optarg);
				break;					
			case 'f':
				fileflag+=1;
				strcpy(file_transient,optarg);
				break;
			case 'g':
				fileflag+=1;
				strcpy(file_timeAverage3D,optarg);
				break;
			case 'x':
				optionflag+=1;
				Re_deltaStern_vorgabe[0]=atof(optarg);
				break;
			case 'y':
				optionflag+=1;
				Re_deltaStern_vorgabe[1]=atof(optarg);
				break;						
			case 's':
				optionflag+=1;
				time_samples=atoi(optarg);
				break;
			case 'w':
				optionflag+=1;
				time_stepwidth=atoi(optarg);
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
	printf("Starting postTKEStreamwiseCascade....\n");
	printf("Die Reynolds-Nummer (y+-Berechnung, option -r (FLOAT)) ist: %g\n",ReynoldsNumber);
	printf("Die Mach-Nummer (option -m (FLOAT)) ist: %g\n",MachNumber);				
	printf("Dieses Tool dient der Erstellung von anströmungsparallelen TKE unter Verwendung einer zeitlich gemittelten und einer transienten Lösung.\n");
	int errorflag=0;
	if(fileflag!=2)
	{
		printf("ERROR: Keine Dateien angegeben (-f 'file_transient.cgns' -g 'file_timeAverage3D.cgns')!\n");
		errorflag=1;
	}
	
	if(optionflag!=4)
	{
		printf("ERROR: Keine Reynolszahlen Re_thetaStern_start/-ende zur Positionsbestimmung angegeben: -x (FLOAT) -y (FLOAT )!\n");
		printf("ERROR: Kein Zeitfenster angegeben (-s SAMPLES(INT) -w SCHRITTWEITE(INT))! Ausgehend vom letzten Sample der Datei und diesen Angaben ist das Zeitfenster definiert.\n");
		errorflag=1;
	}

	if(errorflag==1)
	{
		return 1;
	}
	
	printf("Die turbulente kinetische Energie wird in der ebene (zeta=const.) zwischen Re_deltaStern_vorgabe_start=%f und Re_deltaStern_vorgabe_ende=%f für die lezten %d Samples mit einer Schrittweite von %d ausgewertet.\n ",
	Re_deltaStern_vorgabe[0],Re_deltaStern_vorgabe[1],time_samples,time_stepwidth);
	
	int index_file_timeAverage3D;
	int index_base=1;
	int number_solutions=1;
	int index_zone=1;
	int cell_dim;	
	char zonename_tmp[33];	
	int buffer;
	int iMeshPoints=1;
	int jMeshPoints=1;
	int kMeshPoints=1;		
	int i,k,ijk,ijk2,j;
	float *x,*y,*z;
	float *u_timeAverage,*v_timeAverage,*w_timeAverage,*p_timeAverage,*rho_timeAverage;
		
	cell_dim=3;
	cgsize_t zonesize3D[cell_dim][cell_dim];	
	cgsize_t irmin3D[cell_dim];
	cgsize_t irmax3D[cell_dim];
	CG( cg_open(file_timeAverage3D,CG_MODE_READ,&index_file_timeAverage3D));
	CG( cg_zone_read(index_file_timeAverage3D,index_base,index_zone,zonename_tmp,zonesize3D[0]));


	int i_vorgabe[2];		
	iMeshPoints=zonesize3D[0][0];
	float distance[2];	distance[0]=999999.;	distance[1]=999999.;
	int i2;
	for(i2=0;i2<2;i2++)
	{
	if(Re_deltaStern_vorgabe[i2]>1.0)
	{
	///////////////////////////////Re_detaStern-SUCHE-START//////////////////////////		
	float *Delta;
	Delta=readDelta(iMeshPoints);
	float *Re_deltaStern;
	Re_deltaStern=calloc(sizeof(float),iMeshPoints);		
	printf("Beginne Suche von Re_detaStern_vorgabe: %f\n",Re_deltaStern_vorgabe[i2]);		
	
	for(i=0;i<iMeshPoints;i++)
	{
		Re_deltaStern[i]=Delta[i*6+4];
	}

	for(i=3./4.*iMeshPoints;i<iMeshPoints;i++)
	{
		if(Delta[i*6+0]>1.1)
		{
			break;
		}

		if(fabs(Re_deltaStern[i]-Re_deltaStern_vorgabe[i2])<distance[i2])
		{
			distance[i2]=fabs(Re_deltaStern[i]-Re_deltaStern_vorgabe[i2]);
			i_vorgabe[i2]=i;
		}
	}
	printf("Die ermittelte Position ist: x=%f, Re_deltaStern=%f ||| Vorgabe Re_deltaStern_vorgabe=%f\n",
	Delta[i_vorgabe[i2]*6+0],Re_deltaStern[i_vorgabe[i2]],Re_deltaStern_vorgabe[i2]);	
	///////////////////////////////Re_detaStern-SUCHE-ENDE//////////////////////////		
	}
	else
	{
	///////////////////////////////x-SUCHE-START//////////////////////////	
	float *x_search;
	x_search=(float *)calloc(iMeshPoints, sizeof(float));
	irmin3D[0]=1;		irmin3D[1]=1;	irmin3D[2]=1;	
	irmax3D[0]=iMeshPoints;	irmax3D[1]=1;	irmax3D[2]=1;	
	CG( cg_coord_read(index_file_timeAverage3D,index_base,index_zone,"CoordinateX",
	RealSingle,irmin3D,irmax3D,x_search));	


	float x_vorgabe;
	x_vorgabe=Re_deltaStern_vorgabe[i2];	
	printf("Beginne Suche von x_vorgabe auf Profiloberseite: %f\n",x_vorgabe);		
	for(i=1./2.*iMeshPoints;i<iMeshPoints;i++)
	{
		ijk=0*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
		if(fabs(x_search[ijk]-x_vorgabe)<distance[i2])
		{
			distance[i2]=fabs(x_search[ijk]-x_vorgabe);
			i_vorgabe[i2]=i;
		}
	}
	free(x_search);
	printf("Auswertung bei x %g \n ",x_search[0*jMeshPoints*iMeshPoints+0*iMeshPoints+i_vorgabe[i2]]);
	///////////////////////////////x-SUCHE-SUCHE-ENDE//////////////////////////	
	}
	}



	
	irmin3D[0]=i_vorgabe[0]+1;
	irmin3D[1]=1;	
	irmin3D[2]=1;	
	irmax3D[0]=i_vorgabe[1]+1;
	irmax3D[1]=zonesize3D[0][1];
	irmax3D[2]=1;
	
	printf("Lade Ergebnisse aus Zone %s (%zux%zux%zu)...\n",zonename_tmp,irmax3D[0]-irmin3D[0]+1,irmax3D[1]-irmin3D[1]+1,irmax3D[2]-irmin3D[2]+1);
	iMeshPoints=irmax3D[0]-irmin3D[0]+1;
	jMeshPoints=irmax3D[1]-irmin3D[1]+1;
	kMeshPoints=irmax3D[2]-irmin3D[2]+1;
	
	buffer=iMeshPoints*jMeshPoints*kMeshPoints;

	float *TKE;
	float *u,*v,*w;

	x=(float *)calloc(buffer, sizeof(float));
	y=(float *)calloc(buffer, sizeof(float));
	z=(float *)calloc(buffer, sizeof(float));
	
	u_timeAverage=(float *)calloc(buffer, sizeof(float));
	v_timeAverage=(float *)calloc(buffer, sizeof(float));
	w_timeAverage=(float *)calloc(buffer, sizeof(float));
	p_timeAverage=(float *)calloc(buffer, sizeof(float));
	rho_timeAverage=(float *)calloc(buffer, sizeof(float));		
	
	TKE=(float *)calloc(buffer, sizeof(float));
	
	u=(float *)calloc(buffer, sizeof(float));
	v=(float *)calloc(buffer, sizeof(float));
	w=(float *)calloc(buffer, sizeof(float));
	
	printf("Lesen timeAverage Werte.\n");
	CG( cg_coord_read(index_file_timeAverage3D,index_base,index_zone,"CoordinateX",
	RealSingle,irmin3D,irmax3D,x));
	CG( cg_coord_read(index_file_timeAverage3D,index_base,index_zone,"CoordinateY",
	RealSingle,irmin3D,irmax3D,y));
	CG( cg_coord_read(index_file_timeAverage3D,index_base,index_zone,"CoordinateZ",
	RealSingle,irmin3D,irmax3D,z));

	CG( cg_field_read(index_file_timeAverage3D,index_base,index_zone,number_solutions,"VelocityX_timeAverage",
	RealSingle,irmin3D,irmax3D,u_timeAverage));
	CG( cg_field_read(index_file_timeAverage3D,index_base,index_zone,number_solutions,"VelocityY_timeAverage",
	RealSingle,irmin3D,irmax3D,v_timeAverage));
	CG( cg_field_read(index_file_timeAverage3D,index_base,index_zone,number_solutions,"VelocityZ_timeAverage",
	RealSingle,irmin3D,irmax3D,w_timeAverage));	
	CG( cg_field_read(index_file_timeAverage3D,index_base,index_zone,number_solutions,"Pressure_timeAverage",
	RealSingle,irmin3D,irmax3D,p_timeAverage));
	CG( cg_field_read(index_file_timeAverage3D,index_base,index_zone,number_solutions,"Density_timeAverage",
	RealSingle,irmin3D,irmax3D,rho_timeAverage));		
	printf("fertig.\n");
	
	printf("Öffne Datei %s im read mode\n",file_transient);
	int index_file_transient;	
	CG( cg_open(file_transient,CG_MODE_READ,&index_file_transient));

	int solution;
	int time_start;
	int time_end;
	int NumberSolutions;
	CG( cg_nsols(index_file_transient,1,1,&NumberSolutions ) );
	time_end=NumberSolutions;
	time_start=time_end-time_stepwidth*(time_samples-1);
	if(time_start<0)
	{
		printf("ERROR: time_start ist kleiner 0.\n");
		return 1;
	}
	printf("Auswertefenster: x:%g->x:%g, start:%d - ende:%d - samples:%d.\n",x[0],x[iMeshPoints-1],time_start,time_end,time_stepwidth);
	
	int n=0;
	
	for(solution=time_start;solution<=time_end;solution=solution+time_stepwidth)
	{
		n++;
		printf("Reading solution %d.\n",solution);
		CG( cg_field_read(index_file_transient,index_base,index_zone,solution,"VelocityX",
		RealSingle,irmin3D,irmax3D,u));
		CG( cg_field_read(index_file_transient,index_base,index_zone,solution,"VelocityY",
		RealSingle,irmin3D,irmax3D,v));
		CG( cg_field_read(index_file_transient,index_base,index_zone,solution,"VelocityZ",
		RealSingle,irmin3D,irmax3D,w));
		
		for(k=0;k<buffer;k++)
		{
			
			TKE[k]=0.5*(
			(u_timeAverage[k]-u[k])*(u_timeAverage[k]-u[k])+
			(v_timeAverage[k]-v[k])*(v_timeAverage[k]-v[k])+
			(w_timeAverage[k]-w[k])*(w_timeAverage[k]-w[k]));
		}
		
		printf("Schreibe Datei %d von %d.\n",n,time_samples);
		FILE * file1;
		char output_TKE[300];
		float yPlus,deltaY,viscLength;
		
		viscLength=calcViscLength(iMeshPoints, jMeshPoints, kMeshPoints, (int)((iMeshPoints-1)/2), 0, (kMeshPoints-1), u_timeAverage, v_timeAverage, w_timeAverage, p_timeAverage, rho_timeAverage, x, y, ReynoldsNumber);
		
		sprintf(output_TKE,"TKEStreamwiseCascade_%d.dat",n);
		file1=fopen(output_TKE,"w");
		fprintf(file1,"TITLE = \"TKE-Streamwise-Cascade\"\n");
		fprintf(file1,"variables = \"x/c [-]\"");

		for(j=0;j<jMeshPoints;j++)
		{
			ijk=0*iMeshPoints*jMeshPoints+j*iMeshPoints+0;
			deltaY=sqrt((x[ijk]-x[0])*(x[ijk]-x[0])+(y[ijk]-y[0])*(y[ijk]-y[0]));
			yPlus=deltaY/viscLength;
			fprintf(file1," \"TKE_%d_x%le_y%le_yPlus%le\"",j,x[ijk],y[ijk],yPlus);
		}
		fprintf(file1,"\n");

		fprintf(file1,"ZONE T=\"%s\", F=BLOCK, I=%d, DT=(SINGLE) \n",zonename,iMeshPoints);
		for(i=0;i<iMeshPoints;i++)
		{
			ijk=0*iMeshPoints*jMeshPoints+0*iMeshPoints+i;
			fprintf(file1," %le\n",x[ijk]);
		}
	
		k=0;
		for(j=0;j<jMeshPoints;j++)
		{
		for(i=0;i<iMeshPoints;i++)
		{
			ijk=k*iMeshPoints*jMeshPoints+j*iMeshPoints+i;
			fprintf(file1," %le\n",TKE[ijk]);
		}
		}
		fclose(file1);			
	}

	printf("...fertig. Schließe Datei.\n");
	CG( cg_close(index_file_transient));
	

		
	


	return 0;
}


