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
	
	int fileflag=0;		
	int c;
	int flag_aerodynamicCoefficient=0;
	int flag_uPlus=0;	
	int flag_wallUnits=0;	
	int flag_flucQuant=0;
	int flag_boundaryExtraction=0;
	int flag_gs;
	int optionflag=0;
	opterr = 0;
	char zonename[100];
	char zonename_add[100];
	strcpy(zonename,"Zone");
	strcpy(zonename_add,"add");
	float ReynoldsNumber=500000;
	float MachNumber=0.65;	
	float AoA=0;
	float extractionPosition=0.;

	while ((c = getopt (argc, argv, "f:a:z:ur:m:gwqb:")) != -1)
	{
		switch (c)
		{
			case 'b':
				flag_boundaryExtraction=1;
				extractionPosition=atof(optarg);				
				optionflag=1;				
				break;		
			case 'w':
				flag_wallUnits=1;
				optionflag=1;				
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
			case 'g':
				optionflag=1;
				flag_gs=1;
				break;				
			case 'u':
				optionflag=1;
				flag_uPlus=1;
				break;				
			case 'q':
				optionflag=1;
				flag_flucQuant=1;
				break;								
			case 'a':
				optionflag=1;
				flag_aerodynamicCoefficient=1;
				AoA=atof(optarg);				
				break;	
			case 'z':
				strcpy(zonename_add,optarg);	
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
	printf("Starting postTimeAverage....\n");
	printf("Die Reynolds-Nummer (y+-Berechnung, option -r (FLOAT)) ist: %g\n",ReynoldsNumber);
	printf("Die Mach-Nummer (option -m (FLOAT)) ist: %g\n",MachNumber);			
	printf("Dieses Tool dient der Auswertung von zeitlich gemittelten cgns-Dateien\n");
	int errorflag=0;
	if(fileflag!=1)
	{
		printf("ERROR: Kein Datei angegeben (-f 'file.cgns')!\n");
		errorflag=1;
	}
	
	if(optionflag!=1)
	{
		printf("ERROR: Keine Auswerteoption angegeben: vorticity: -v | aerodynmic Coefficients -a | turbulent kinetic energy -t | uPlus profiles -u (FLOAT) [Re_deltaStern(displacement thickness)_vorgabe] | delta99,deltaStern-plots-g | wallunits -w | fluctuating Quantities -q (FLOAT) [Re_deltaStern(displacement thickness)_vorgabe]\n");
		errorflag=1;
	}

	if(errorflag==1)
	{
		return 1;
	}
	
	
	
	
	int index_base=1;
	int number_solutions=1;
	int index_zone=1;

	char zonename_tmp[33];	

	int iMeshPoints;
	int jMeshPoints;
	int kMeshPoints;		
	int i,j,k,ijk;
		

	float *x,*y;	
	float *u_timeAverage2D,*v_timeAverage2D,*w_timeAverage2D,*p_timeAverage2D,*rho_timeAverage2D;
	float *u_varianz2D,*v_varianz2D,*w_varianz2D,*p_varianz2D,*rho_varianz2D;			

	char file2D[500];

	strcpy(file2D,file);
	strcpy(file2D,replace_str(file2D, "3d","2d"));				
	
	printf("Öffne Datei %s im read mode\n",file2D);
	int index_file2D;

	int cell_dim2D;	

	CG( cg_open(file2D,CG_MODE_READ,&index_file2D));
	CG( cg_cell_dim(index_file2D,index_base,&cell_dim2D));

	cgsize_t zonesize2D[cell_dim2D][cell_dim2D];	

	CG( cg_zone_read(index_file2D,index_base,index_zone,zonename_tmp,zonesize2D[0]));

	cgsize_t irmin2D[cell_dim2D];
	cgsize_t irmax2D[cell_dim2D];

	int buffer2D;

	irmin2D[0]=1;
	irmin2D[1]=1;	
	irmax2D[0]=zonesize2D[0][0];
	irmax2D[1]=zonesize2D[0][1];
	
	iMeshPoints=irmax2D[0];
	jMeshPoints=irmax2D[1];

	buffer2D=zonesize2D[0][0]*zonesize2D[0][1];	
	printf("Lade 2D-Ergebnisse aus Zone %s (%dx%d)...\n",zonename_tmp,iMeshPoints,jMeshPoints);
	
	
	x=(float *)calloc(buffer2D, sizeof(float));
	y=(float *)calloc(buffer2D, sizeof(float));	
	u_timeAverage2D=(float *)calloc(buffer2D, sizeof(float));
	v_timeAverage2D=(float *)calloc(buffer2D, sizeof(float));
	w_timeAverage2D=(float *)calloc(buffer2D, sizeof(float));
	p_timeAverage2D=(float *)calloc(buffer2D, sizeof(float));
	rho_timeAverage2D=(float *)calloc(buffer2D, sizeof(float));				
	u_varianz2D=(float *)calloc(buffer2D, sizeof(float));
	v_varianz2D=(float *)calloc(buffer2D, sizeof(float));
	w_varianz2D=(float *)calloc(buffer2D, sizeof(float));
	p_varianz2D=(float *)calloc(buffer2D, sizeof(float));
	rho_varianz2D=(float *)calloc(buffer2D, sizeof(float));		
	
	
	CG( cg_coord_read(index_file2D,index_base,index_zone,"CoordinateX",
	RealSingle,irmin2D,irmax2D,x));
	CG( cg_coord_read(index_file2D,index_base,index_zone,"CoordinateY",
	RealSingle,irmin2D,irmax2D,y));	
	
	CG( cg_field_read(index_file2D,index_base,index_zone,number_solutions,"VelocityX_timeAverage",
	RealSingle,irmin2D,irmax2D,u_timeAverage2D));
	CG( cg_field_read(index_file2D,index_base,index_zone,number_solutions,"VelocityY_timeAverage",
	RealSingle,irmin2D,irmax2D,v_timeAverage2D));
	CG( cg_field_read(index_file2D,index_base,index_zone,number_solutions,"VelocityZ_timeAverage",
	RealSingle,irmin2D,irmax2D,w_timeAverage2D));	
	CG( cg_field_read(index_file2D,index_base,index_zone,number_solutions,"Pressure_timeAverage",
	RealSingle,irmin2D,irmax2D,p_timeAverage2D));
	CG( cg_field_read(index_file2D,index_base,index_zone,number_solutions,"Density_timeAverage",
	RealSingle,irmin2D,irmax2D,rho_timeAverage2D));		
	
	printf("Lade Varianzgroessen.\n");

	
	CG( cg_field_read(index_file2D,index_base,index_zone,number_solutions,"VelocityX_varianz",
	RealSingle,irmin2D,irmax2D,u_varianz2D));
	CG( cg_field_read(index_file2D,index_base,index_zone,number_solutions,"VelocityY_varianz",
	RealSingle,irmin2D,irmax2D,v_varianz2D));
	CG( cg_field_read(index_file2D,index_base,index_zone,number_solutions,"VelocityZ_varianz",
	RealSingle,irmin2D,irmax2D,w_varianz2D));	
	CG( cg_field_read(index_file2D,index_base,index_zone,number_solutions,"Pressure_varianz",
	RealSingle,irmin2D,irmax2D,p_varianz2D));
	CG( cg_field_read(index_file2D,index_base,index_zone,number_solutions,"Density_varianz",
	RealSingle,irmin2D,irmax2D,rho_varianz2D));	
	printf("...fertig. Schließe 2D-Datei.\n");
	CG( cg_close(index_file2D));		

	printf("Öffne Datei %s im read mode\n",file);		
	int index_file;	
	CG( cg_open(file,CG_MODE_READ,&index_file));
	int cell_dim3D;			
	CG( cg_cell_dim(index_file,index_base,&cell_dim3D));
	if(cell_dim3D==2)
	{
		printf("ERROR: Auswertung nur mit 3D-Dateien möglich.\n");
		return 0;
	}
	cgsize_t zonesize[cell_dim3D][cell_dim3D];	

	CG( cg_zone_read(index_file,index_base,index_zone,zonename_tmp,zonesize[0]));
	kMeshPoints=zonesize[0][2];				
	cgsize_t irmin[cell_dim3D];
	cgsize_t irmax[cell_dim3D];

	
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////GRENZSCHICHTDICKEN//////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	if(flag_gs==1)
	{
		printf("Bestimme Grenzschichtdicken.\n");
		printf("Achtung!!! Grenzschichtdickenbestimmung mit 3D-Datei.\n");
		int buffer;
		irmin[0]=1;
		irmin[1]=1;	
		irmin[2]=1;	
		irmax[0]=zonesize[0][0];
		irmax[1]=zonesize[0][1];
		irmax[2]=zonesize[0][2];
		buffer=zonesize[0][0]*zonesize[0][1]*zonesize[0][2];	
		printf("Lade 3D-Ergebnisse aus Zone %s (%dx%dx%d)...\n",zonename_tmp,iMeshPoints,jMeshPoints,kMeshPoints);

		float *x3D,*y3D,*z3D;	
		float *u_timeAverage,*v_timeAverage,*w_timeAverage;

		x3D=(float *)calloc(buffer, sizeof(float));
		y3D=(float *)calloc(buffer, sizeof(float));
		z3D=(float *)calloc(buffer, sizeof(float));		
		
		u_timeAverage=(float *)calloc(buffer, sizeof(float));
		v_timeAverage=(float *)calloc(buffer, sizeof(float));
		w_timeAverage=(float *)calloc(buffer, sizeof(float));
		
		CG( cg_coord_read(index_file,index_base,index_zone,"CoordinateX",
		RealSingle,irmin,irmax,x3D));
		CG( cg_coord_read(index_file,index_base,index_zone,"CoordinateY",
		RealSingle,irmin,irmax,y3D));
		CG( cg_coord_read(index_file,index_base,index_zone,"CoordinateZ",
		RealSingle,irmin,irmax,z3D));						
	
		CG( cg_field_read(index_file,index_base,index_zone,number_solutions,"VelocityX_timeAverage",
		RealSingle,irmin,irmax,u_timeAverage));
		CG( cg_field_read(index_file,index_base,index_zone,number_solutions,"VelocityY_timeAverage",
		RealSingle,irmin,irmax,v_timeAverage));
		CG( cg_field_read(index_file,index_base,index_zone,number_solutions,"VelocityZ_timeAverage",
		RealSingle,irmin,irmax,w_timeAverage));	

		printf("...fertig. Schließe 3D-Datei.\n");
		CG( cg_close(index_file));		
		
		
		
		
		
		printf("Grenzschichtdicken...\n");
		float *delta99_verlauf=(float*)calloc(iMeshPoints,sizeof(float));;
		float *delta99_j=(float*)calloc(iMeshPoints,sizeof(float));;
		float *deltaStern_verlauf=(float*)calloc(iMeshPoints,sizeof(float));;
		float *theta_verlauf=(float*)calloc(iMeshPoints,sizeof(float));;
		float *Re_deltaStern_verlauf=(float*)calloc(iMeshPoints,sizeof(float));;
		float *DeltaCalc;

		printf("Beginne Grenzschichtdickenverlaufs-Erstellung\n");
		for(i=0;i<iMeshPoints;i++)
		{
			DeltaCalc=calcDelta(i, iMeshPoints, jMeshPoints, kMeshPoints, x3D, y3D, z3D, u_timeAverage, v_timeAverage, w_timeAverage);

			delta99_verlauf[i]=DeltaCalc[0];
			delta99_j[i]=DeltaCalc[3];
			deltaStern_verlauf[i]=DeltaCalc[1];
			theta_verlauf[i]=DeltaCalc[2];
			Re_deltaStern_verlauf[i]=ReynoldsNumber*deltaStern_verlauf[i];
			
			free(DeltaCalc);
		}
		printf("Beginne Grenzschichtdickenverlaufs-Export\n");
		FILE * file3;
		char output_gs[300];
		sprintf(zonename,"%s",zonename_add);
		sprintf(output_gs,"Grenzschichtdickenverlaeufe.dat");
		file3=fopen(output_gs,"w");
		fprintf(file3,"variables = \"x/c [-]\" \"<greek>d</greek><sub>99</sub> [-]\" \"<greek>d</greek><sup>*</sup> [-]\" \"<greek>Q</greek> [-]\" \"Re<sub><greek>d</greek>*</sub> [-]\" \"<greek>d</greek><sub>99,j</sub> [-]\"\n");
		fprintf(file3,"zone t=\"%s\", i= %d, f=point \n",zonename,iMeshPoints);
		for(i=0;i<iMeshPoints;i++)
		{
			ijk=0*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
			fprintf(file3," %le %le %le %le %le %le\n",x3D[ijk],delta99_verlauf[i],deltaStern_verlauf[i],theta_verlauf[i],Re_deltaStern_verlauf[i],delta99_j[i]);

		}
		fclose(file3);	

	}

	
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////aerodynamic coefficients//////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	if(flag_aerodynamicCoefficient==1)
	{
		printf("aerodynamic coefficients...\n");	
		float *cp_verlauf=(float*)calloc(iMeshPoints,sizeof(float));;
		float *cl_verlauf=(float*)calloc(iMeshPoints,sizeof(float));;		
		float *cf_verlauf=(float*)calloc(iMeshPoints,sizeof(float));;
		float *cd_verlauf=(float*)calloc(iMeshPoints,sizeof(float));;						
		float cp_mean=0.0;
		float cl_mean=0.0;		
		float cf_mean=0.0;
		float cd_mean=0.0;			
		float alpha;			
		float deltaX;
		float totalX=0.0;
		float cl_withoutAoA;
		float cd_withoutAoA;
		int ij1k,ij0k,i1jk;
		int iPlus1jk,iMinus1jk;
		float vorzeichen;
		printf("Beginne Verlaufs-Erstellung der aerodynamischen Koeffizienten\n");
		j=1;
		for(i=1;i<iMeshPoints;i++)
		{
			k=0;
			ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
			i1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;						
			alpha=atan((y[ijk]-y[i1jk])/(x[ijk]-x[i1jk]));

			ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
			ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
			ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;								
			cp_verlauf[i]=(p_timeAverage2D[ijk]-1.)/(0.5*1.4*MachNumber*MachNumber);
			cf_verlauf[i]=calcCF(iMeshPoints, iMeshPoints, iMeshPoints, i, 1, k, u_timeAverage2D, v_timeAverage2D, w_timeAverage2D, p_timeAverage2D, rho_timeAverage2D, x, y, ReynoldsNumber);

			k=0;
			iMinus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;			
			if(y[iMinus1jk]>=0){vorzeichen=1.0;}
			else{vorzeichen=-1.0;}

			cl_withoutAoA=sin(alpha)*cf_verlauf[i]-vorzeichen*cos(alpha)*(cp_verlauf[i]+1./(0.5*1.4*MachNumber*MachNumber));			
			cd_withoutAoA=cos(alpha)*cf_verlauf[i]+vorzeichen*sin(alpha)*(cp_verlauf[i]+1./(0.5*1.4*MachNumber*MachNumber));

			cl_verlauf[i]=cos(AoA)*cl_withoutAoA+sin(AoA)*cd_withoutAoA;
			cd_verlauf[i]=sin(AoA)*cl_withoutAoA+cos(AoA)*cd_withoutAoA;
			//printf("i: %d cf=%g, cp=%g, cl=%g, cd=%g \n",i,cp_verlauf[i],cf_verlauf[i],cl_verlauf[i],cd_verlauf[i]);
			k=0;
			ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;			
			iPlus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i+1;
			iMinus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;						
			if(
			((x[iPlus1jk]<=1.0)&&(y[ijk]>0.0))
			||
		   	((x[iMinus1jk]<=1.0)&&(y[ijk]<0.0))
		   	)
			{
				//printf("i:%d x:%f y:%f\n",i,x[ijk],y[ijk]);
				deltaX=0.5*sqrt((x[iPlus1jk]-x[iMinus1jk])*(x[iPlus1jk]-x[iMinus1jk])+(y[iPlus1jk]-y[iMinus1jk])*(y[iPlus1jk]-y[iMinus1jk]));
				totalX+=deltaX;
				cp_mean+=cp_verlauf[i]*deltaX;
				cf_mean+=cf_verlauf[i]*deltaX;
				cl_mean+=cl_verlauf[i]*deltaX;
				cd_mean+=cd_verlauf[i]*deltaX;						
			}
			else if((x[ijk]==1.0)&&(y[ijk]>0.0))
			{
				//printf("i:%d x:%f y:%f\n",i,x[ijk],y[ijk]);
				deltaX=sqrt((x[ijk]-x[iMinus1jk])*(x[ijk]-x[iMinus1jk])+(y[ijk]-y[iMinus1jk])*(y[ijk]-y[iMinus1jk]));
				totalX+=deltaX;
				cp_mean+=cp_verlauf[i]*deltaX/2.;
				cf_mean+=cf_verlauf[i]*deltaX/2.;
				cl_mean+=cl_verlauf[i]*deltaX/2.;
				cd_mean+=cd_verlauf[i]*deltaX/2.;						
			}			
			else if((x[ijk]==1.0)&&(y[ijk]<0.0))
			{
				//printf("i:%d x:%f y:%f\n",i,x[ijk],y[ijk]);
				deltaX=sqrt((x[ijk]-x[iPlus1jk])*(x[ijk]-x[iPlus1jk])+(y[ijk]-y[iPlus1jk])*(y[ijk]-y[iPlus1jk]));
				totalX+=deltaX;
				cp_mean+=cp_verlauf[i]*deltaX/2.;
				cf_mean+=cf_verlauf[i]*deltaX/2.;
				cl_mean+=cl_verlauf[i]*deltaX/2.;
				cd_mean+=cd_verlauf[i]*deltaX/2.;
			}			
		}
		cp_mean=cp_mean/totalX;		
		printf("cp-mean: %f\n",cp_mean);				
		cf_mean=cf_mean/totalX;		
		printf("cf-mean: %f\n",cf_mean);
		cl_mean=cl_mean/totalX;		
		printf("cl-mean: %f\n",cl_mean);
		cd_mean=cd_mean/totalX;		
		printf("cd-mean: %f\n",cd_mean);						
		printf("Beginne aerodynamischenKoeffizienten-verlaufs-Export\n");
		FILE * file4;
		char output[300];
		sprintf(zonename,"%s",zonename_add);
		sprintf(output,"aerodynmicCoefficients.dat");
		file4=fopen(output,"w");
		fprintf(file4,"variables = \"x/c [-]\" \"c<sub>p</sub> [-]\"  \"c<sub>p,mean</sub> [-]\"  \"c<sub>f</sub> [-]\"  \"c<sub>f,mean</sub> [-]\"  \"c<sub>l</sub> [-]\"  \"c<sub>l,mean</sub> [-]\"  \"c<sub>d</sub> [-]\"  \"c<sub>d,mean</sub> [-]\"\n");
		fprintf(file4,"zone t=\"%s\", i= %d, f=point \n",zonename,iMeshPoints);
		for(i=0;i<iMeshPoints;i++)
		{
			ijk=0*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
			fprintf(file4," %le %le %le %le %le %le %le %le %le\n",x[ijk],cp_verlauf[i],cp_mean,cf_verlauf[i],cf_mean,cl_verlauf[i],cl_mean,cd_verlauf[i],cd_mean);

		}
		fclose(file4);	

	}
	
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////Wall Units//////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	if(flag_wallUnits==1)
	{
		printf("Wall Units...\n");		
		float *deltaXPlus_verlauf=(float*)calloc(iMeshPoints,sizeof(float));;
		float *deltaYPlus_verlauf=(float*)calloc(iMeshPoints,sizeof(float));;
		float *deltaZPlus_verlauf=(float*)calloc(iMeshPoints,sizeof(float));;
		int ij1k,ij0k,i1jk;
		float deltaX,deltaY,deltaZ,viscLength;
		ijk=0*jMeshPoints*iMeshPoints+0*iMeshPoints+0;

		deltaZ=0.1/(float)(kMeshPoints-1.);		
		printf("Beginne WallUnits-verlaufs-Erstellung\n");
		j=0;
		for(i=0;i<iMeshPoints;i++)
		{
			k=0;
			ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;			
			if (i==0){i1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i+1;}
			else{i1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;}
			deltaX=sqrt((x[ijk]-x[i1jk])*(x[ijk]-x[i1jk])+(y[ijk]-y[i1jk])*(y[ijk]-y[i1jk]));

			ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
			ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;
			deltaY=sqrt((x[ij1k]-x[ij0k])*(x[ij1k]-x[ij0k])+(y[ij1k]-y[ij0k])*(y[ij1k]-y[ij0k]));
			viscLength=calcViscLength(iMeshPoints, jMeshPoints, kMeshPoints, i, 0, k, u_timeAverage2D, v_timeAverage2D, w_timeAverage2D, p_timeAverage2D, rho_timeAverage2D, x, y, ReynoldsNumber);
			deltaYPlus_verlauf[i]=deltaY/viscLength;
			deltaXPlus_verlauf[i]=deltaX/viscLength;
			deltaZPlus_verlauf[i]=deltaZ/viscLength;
		}
		printf("Beginne WallUnits-verlaufs-Export\n");
		FILE * file6;
		char output_gs[300];
		sprintf(zonename,"%s",zonename_add);
		sprintf(output_gs,"WallUnitsverlaeufe.dat");
		file6=fopen(output_gs,"w");
		fprintf(file6,"variables = \"x/c [-]\" \"<greek>D</greek>x<sup>+</sup> [-]\" \"<greek>D</greek>y<sup>+</sup> [-]\" \"<greek>D</greek>z<sup>+</sup> [-]\"\n");
		fprintf(file6,"zone t=\"%s\", i= %d, f=point \n",zonename,iMeshPoints);
		for(i=0;i<iMeshPoints;i++)
		{
			ijk=0*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
			fprintf(file6," %le %le %le %le\n",x[ijk],deltaXPlus_verlauf[i],deltaYPlus_verlauf[i],deltaZPlus_verlauf[i]);

		}
		fclose(file6);	

	}	


	/////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////UPLUS//////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////	
	if(flag_uPlus==1)
	{
		printf("UPLUS...\n");


		///////////////////////////////x-SUCHE-START//////////////////////////	
		float x_vorgabe;
		x_vorgabe=0.815;

		int i_vorgabe;		
		float distance=99999.;		
		printf("Beginne Suche von x_vorgabe auf Profiloberseite: %f\n",x_vorgabe);		
		for(i=1./2.*iMeshPoints;i<iMeshPoints;i++)
		{
			ijk=0*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
			if(fabs(x[ijk]-x_vorgabe)<distance)
			{
				distance=fabs(x[ijk]-x_vorgabe);
				i_vorgabe=i;
			}
		}
		printf("Auswertung bei x,y %g,%g \n ",x[0*jMeshPoints*iMeshPoints+0*iMeshPoints+i_vorgabe],y[0*jMeshPoints*iMeshPoints+0*iMeshPoints+i_vorgabe]);
		///////////////////////////////x-SUCHE-SUCHE-ENDE//////////////////////////	

		float yPlus,uPlus,uxPlus,uPlusVD;
		float M_tau,Pr_tau,H,R,D;
		FILE * file2;
		char output_uPlus[300];
		sprintf(zonename,"%s",zonename_add);
		sprintf(output_uPlus,"MeanVelocityProfile.dat");
		file2=fopen(output_uPlus,"w");
		fprintf(file2,"variables = \"y<sup>+</sup>\" \"u<sup>+</sup> [-]\" \"u<sub>x</sub><sup>+</sup> [-]\"\n");
		fprintf(file2,"zone t=\"%s\", i= %d, f=point \n",zonename,jMeshPoints-1);
		//i=position_xi;
		int ij0k,ij1k;
		int iPlus1jk,iMinus1jk;
		int ijPlus1k,ijMinus1k;	

		i=i_vorgabe;
		k=0;
		fprintf(file2," %le %le %le\n",0.,0.,0.);
		
		float uTau;
		float deltaY;
		float deltaY_w;
		float mue, dpdx, uP, lambda, uPlusAdvPreGra;
		ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
		ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
		ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;
				
		j=0;
		iPlus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i+1;
		iMinus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;
		ijPlus1k=k*jMeshPoints*iMeshPoints+(j+1)*iMeshPoints+i;
		ijMinus1k=k*jMeshPoints*iMeshPoints+(j-1)*iMeshPoints+i;			

		uTau=calcUTau(iMeshPoints, iMeshPoints, iMeshPoints, i, 1, k, u_timeAverage2D, v_timeAverage2D, w_timeAverage2D, p_timeAverage2D, rho_timeAverage2D, x, y, ReynoldsNumber);
		ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
		ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;		
		float viscLength=calcViscLength(iMeshPoints, jMeshPoints, kMeshPoints, i, 1, k, u_timeAverage2D, v_timeAverage2D, w_timeAverage2D, p_timeAverage2D, rho_timeAverage2D, x, y, ReynoldsNumber);
		float x_xi,y_xi,x_eta,y_eta;
		float jacobian;
		float xi_x,xi_y,eta_x,eta_y;
		float dp_xi,dp_eta;
		deltaY_w=sqrt(pow((x[ij1k]-x[ij0k]),2.)+pow((y[ij1k]-y[ij0k]),2.));
		yPlus=deltaY_w/2./viscLength;
		for(j=1;j<jMeshPoints-1;j++)
		{
			ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;

			ijMinus1k=k*jMeshPoints*iMeshPoints+(j-1)*iMeshPoints+i;			
			ijPlus1k=k*jMeshPoints*iMeshPoints+(j+1)*iMeshPoints+i;
			iPlus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i+1;
			iMinus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;

			deltaY=sqrt(pow((x[ijk]-x[ijMinus1k]),2.)+pow((y[ijk]-y[ijMinus1k]),2.))/viscLength;
			yPlus=yPlus+deltaY;
			//printf("x:%f, y: %f, deltay:%f, yPlus:%f,u:%f,v:%f,w:%f,p:%f,rho:%f,\n",x[ijk],y[ijk],sqrt(pow((x[ijk]-x[ij0k]),2.)+pow((y[ijk]-y[ij0k]),2.)),yPlus,u_timeAverage2D[ijk],v_timeAverage2D[ijk],w_timeAverage2D[ijk],p_timeAverage2D[ijk],rho_timeAverage2D[ijk]);			

			uPlus=sqrt(u_timeAverage2D[ijk]*u_timeAverage2D[ijk]+v_timeAverage2D[ijk]*v_timeAverage2D[ijk]+w_timeAverage2D[ijk]*w_timeAverage2D[ijk])/uTau;
			
			mue=((1.0+110.4/300.0)*pow((p_timeAverage2D[ijk]/rho_timeAverage2D[ijk]),1.5))/(p_timeAverage2D[ijk]/rho_timeAverage2D[ijk]+110.4/300.0);
			
			x_xi=(-0.5*x[iMinus1jk]+0.5*x[iPlus1jk]);	y_xi=(-0.5*y[iMinus1jk]+0.5*y[iPlus1jk]);
			x_eta=(-0.5*x[ijMinus1k]+0.5*x[ijPlus1k]);	y_eta=(-0.5*y[ijMinus1k]+0.5*y[ijPlus1k]);
			jacobian=x_xi*y_eta-x_eta*y_xi;
			xi_x=(y_eta)/jacobian;	xi_y=(-x_eta)/jacobian;
			eta_x=(-y_xi)/jacobian;	eta_y=(x_xi)/jacobian;

			dp_xi=(0.5*p_timeAverage2D[iPlus1jk]-0.5*p_timeAverage2D[iMinus1jk]);
			dp_eta=(0.5*p_timeAverage2D[ijPlus1k]-0.5*p_timeAverage2D[ijMinus1k]);
			
			dpdx=dp_xi*xi_x+dp_eta*eta_x;
			//dpdx=(p_timeAverage2D[iPlus1jk]-p_timeAverage2D[iMinus1jk])/sqrt((x[iPlus1jk]-x[iMinus1jk])*(x[iPlus1jk]-x[iMinus1jk])+(y[iPlus1jk]-y[iMinus1jk])*(y[iPlus1jk]-y[iMinus1jk]));
			if(dpdx<0)
			 dpdx=0.0;

			uP=pow(((mue/rho_timeAverage2D[ijk]/rho_timeAverage2D[ijk]*dpdx)/ReynoldsNumber/MachNumber/1.4),1./3.);
			lambda=pow((uP/uTau),3.0);
			uPlusAdvPreGra=1./0.41*(log(yPlus)-2.*log((sqrt(1.+lambda*yPlus)+1.)/2.)+2.*(sqrt(1.+lambda*yPlus)-1))+5.24;
			//printf("dpdx:%e uplus%e\n",dpdx,uPlusAdvPreGra);
			
			M_tau=uTau/sqrt(p_timeAverage2D[ij0k]/rho_timeAverage2D[ij0k])*MachNumber;
			Pr_tau=0.9;
			R=M_tau*sqrt((1.4-1.)*Pr_tau/2.);
			H=0.;//adiabat
			D=sqrt(1.+R*R*H*H);
			uPlusVD=1./R*(asin(R/D*(uPlus+H))-asin(R*H/D));
			uxPlus=u_timeAverage2D[ijk]/uTau;
			fprintf(file2," %le %le %le\n",yPlus,uPlus,uPlusAdvPreGra);
		
		}
		fclose(file2);	

		
	}


	/////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////FluctuatingQuantities////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////	
	if(flag_flucQuant==1)
	{
		float x_vorgabe;
		x_vorgabe=0.8;
		printf("FluctuatingQuantities at x=%f...\n",x_vorgabe);
		///////////////////////////////x-SUCHE-START//////////////////////////	


		int i_vorgabe;		
		float distance=99999.;		
		printf("Beginne Suche von x_vorgabe auf Profiloberseite: %f\n",x_vorgabe);		
		for(i=1./2.*iMeshPoints;i<iMeshPoints;i++)
		{
			ijk=0*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
			if(fabs(x[ijk]-x_vorgabe)<distance)
			{
				distance=fabs(x[ijk]-x_vorgabe);
				i_vorgabe=i;
			}
		}
		printf("Auswertung bei x,y %g,%g \n ",x[0*jMeshPoints*iMeshPoints+0*iMeshPoints+i_vorgabe],y[0*jMeshPoints*iMeshPoints+0*iMeshPoints+i_vorgabe]);
		///////////////////////////////x-SUCHE-SUCHE-ENDE//////////////////////////			


		float yPlus,uu,vv,ww,uv,tke,utau;
		FILE * filefq;
		char output_fq[300];
		sprintf(zonename,"%s",zonename_add);
		sprintf(output_fq,"FluctuatingQuantities.dat");
		filefq=fopen(output_fq,"w");
		fprintf(filefq,"variables = \"y/c [-]\" \"y<sup>+</sup>\" \"<greek>D</greek>/<greek>d</greek><sub>99</sub>\" \"<u'u'>\" \"<v'v'>\" \"<w'w'>\" \"<u'v'>\"  \"<u'u'>/k\" \"<v'v'>/k\" \"<w'w'>/k\" \"<u'v'>/k\" \"<u'u'>/u<sub><greek>t</greek></sub><sup>2</sup>\" \"<v'v'>/u<sub><greek>t</greek></sub><sup>2</sup>\" \"<w'w'>/u<sub><greek>t</greek></sub><sup>2</sup>\" \"<u'v'>/u<sub><greek>t</greek></sub><sup>2</sup>\" \"TKE\" \"u<sub><greek>t</greek></sub>\"\n");
		fprintf(filefq,"zone t=\"%s\", i= %d, f=point \n",zonename,jMeshPoints-2);

		int ij0k,ij1k;
		int iPlus1jk,iMinus1jk;
		int ijPlus1k,ijMinus1k;


		float y_mesh;
		float y_delta,delta99;
		float *Delta;
		float deltaY_w;
		float deltaY,deltaY2;		
		Delta=readDelta(iMeshPoints);
		fprintf(filefq," %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
		i=i_vorgabe;
		k=0;
		
		ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
		ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;
		float viscLength=calcViscLength(iMeshPoints, jMeshPoints, kMeshPoints, i, 1, k, u_timeAverage2D, v_timeAverage2D, w_timeAverage2D, p_timeAverage2D, rho_timeAverage2D, x, y, ReynoldsNumber);
		
		deltaY_w=sqrt(pow((x[ij1k]-x[ij0k]),2.)+pow((y[ij1k]-y[ij0k]),2.));
		yPlus=deltaY_w/2./viscLength;
		
		delta99=0.0226689; //Wert von M6 bei 0.8
		y_delta=deltaY_w/2./delta99;
				
		for(j=1;j<jMeshPoints-1;j++)
		{
			iPlus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i+1;
			iMinus1jk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i-1;
			ijPlus1k=k*jMeshPoints*iMeshPoints+(j+1)*iMeshPoints+i;
			ijMinus1k=k*jMeshPoints*iMeshPoints+(j-1)*iMeshPoints+i;
			
			ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
			ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
			ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;
			
			deltaY=sqrt(pow((x[ijk]-x[ijMinus1k]),2.)+pow((y[ijk]-y[ijMinus1k]),2.))/viscLength;
			yPlus=yPlus+deltaY;
			
			utau=calcUTau(iMeshPoints, iMeshPoints, iMeshPoints, i, 1, k, u_timeAverage2D, v_timeAverage2D, w_timeAverage2D, p_timeAverage2D, rho_timeAverage2D, x, y, ReynoldsNumber);

			delta99=Delta[i*6+1];
			
			deltaY2=sqrt(pow((x[ijk]-x[ijMinus1k]),2.)+pow((y[ijk]-y[ijMinus1k]),2.))/delta99;
			y_delta=y_delta+deltaY2;
			
			uu=u_varianz2D[ijk];
			vv=v_varianz2D[ijk];
			ww=w_varianz2D[ijk];
			uv=sqrt(fabs(u_varianz2D[ijk]))*sqrt(fabs(v_varianz2D[ijk]));
			
			tke=0.5*(u_varianz2D[ijk]+v_varianz2D[ijk]+w_varianz2D[ijk]);
			y_mesh=y[ijk];
	
			fprintf(filefq," %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",y_mesh,yPlus,y_delta,uu,vv,ww,uv,uu/tke,vv/tke,ww/tke,uv/tke,uu/utau/utau,vv/utau/utau,ww/utau/utau,uv/utau/utau,tke,utau);
		
		}
		fclose(filefq);
		
		
		printf("FluctuatingQuantities along airfoil...\n");

		sprintf(zonename,"%s",zonename_add);
		sprintf(output_fq,"FluctuatingQuantities_x.dat");
		filefq=fopen(output_fq,"w");
		fprintf(filefq,"variables = \"x/c [-]\" \"TKE<sub>MAX</sub>\" \"u<sub><greek>t</greek></sub>\"\n");
		fprintf(filefq,"zone t=\"%s\", i= %d, f=point \n",zonename,iMeshPoints);

		float tke_max;
		float tke_tmp;
		int j_max;
		k=0;
		for(i=0;i<iMeshPoints;i++)
		{		
		j_max=(int)Delta[i*6+5];
		tke_max=0.0;
		ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
		ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;
		for(j=0;j<j_max;j++)
		{
			ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;

			utau=calcUTau(iMeshPoints, iMeshPoints, iMeshPoints, i, 1, k, u_timeAverage2D, v_timeAverage2D, w_timeAverage2D, p_timeAverage2D, rho_timeAverage2D, x, y, ReynoldsNumber);
			tke_tmp=0.5*(u_varianz2D[ijk]+v_varianz2D[ijk]+w_varianz2D[ijk]);
			if (tke_tmp>tke_max)
			{
				tke_max=tke_tmp;
			}
	
		}
		fprintf(filefq," %le %le %le\n",x[ij0k],tke_max,utau);		
		}
		fclose(filefq);				
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////BoundaryExtration////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////	
	if(flag_boundaryExtraction==1)
	{
		printf("Boundary Extraction...\n");
		float *Delta;
		Delta=readDelta(iMeshPoints);
		int iStart,iEnde;
		float x_vorgabe;
		if (extractionPosition<0.)
		{
			iStart=0;
			iEnde=1./2.*iMeshPoints;
			x_vorgabe=fabs(extractionPosition);
			printf("Beginne Suche von x_vorgabe auf Profilunterseite: %f\n",x_vorgabe);		
		}
		else
		{
			iStart=1./2.*iMeshPoints;
			iEnde=iMeshPoints;
			x_vorgabe=fabs(extractionPosition);
			printf("Beginne Suche von x_vorgabe auf Profiloberseite: %f\n",x_vorgabe);					
		}
		///////////////////////////////x-SUCHE-START//////////////////////////	

		int i_vorgabe;		
		float distance=99999.;		

		for(i=iStart;i<iEnde;i++)
		{
			ijk=0*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
			if(fabs(x[ijk]-x_vorgabe)<distance)
			{
				distance=fabs(x[ijk]-x_vorgabe);
				i_vorgabe=i;
			}
		}
		i=i_vorgabe;
		printf("Auswertung bei x,y %g,%g \n ",x[0*jMeshPoints*iMeshPoints+0*iMeshPoints+i],y[0*jMeshPoints*iMeshPoints+0*iMeshPoints+i]);

		///////////////////////////////x-SUCHE-SUCHE-ENDE//////////////////////////	

		float Velocity;
		FILE * file7;
		char output_eB[300];
		float re_deltaStern;
		sprintf(zonename,"%s",zonename_add);
		sprintf(output_eB,"BoundaryLayer_%.2f.dat",extractionPosition);
		file7=fopen(output_eB,"w");
		fprintf(file7,"variables = \"y\" \"u\" \"deltaStern\" \"delta99\" \"Re_deltaStern\" \"Re_deltaStern\"\n");
		fprintf(file7,"zone t=\"%s\", i= %d, f=point \n",zonename,jMeshPoints);
		//i=position_xi;
		int ij0k,ij1k;

		k=0;
		float deltaStern=Delta[i*6+2];
		float delta99=Delta[i*6+1];				
		
		float deltaY;
		float deltaY_w;
		ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
		ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
		ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;
				
		float relDelta;
		float VelocityOld=0.00001;
		float VelocityFreeStream;
		int j_break;
		for(j=0;j<jMeshPoints-1;j++)
		{
			ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
			Velocity=sqrt(u_timeAverage2D[ijk]*u_timeAverage2D[ijk]+v_timeAverage2D[ijk]*v_timeAverage2D[ijk]);
			relDelta=fabs(Velocity-VelocityOld)/VelocityOld*100.;
			VelocityOld=Velocity;
			if(relDelta<0.5)
			{
				VelocityFreeStream=Velocity;
				j_break=j;
				re_deltaStern=ReynoldsNumber*rho_timeAverage2D[ijk]*VelocityFreeStream*deltaStern;
				break;
			}
		}

		printf("VelocityFreeStream=%g bei j=%d\n",VelocityFreeStream,j_break);
		
		fprintf(file7," %le %le %le %le %le %le\n",0.,0.,deltaStern,delta99,re_deltaStern,VelocityFreeStream);
		for(j=0;j<j_break;j++)
		{
			ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
			ij0k=k*jMeshPoints*iMeshPoints+0*iMeshPoints+i;
			ij1k=k*jMeshPoints*iMeshPoints+1*iMeshPoints+i;
		
			if (j==0)
			{
				deltaY=sqrt(pow((x[ij1k]-x[ij0k]),2.)+pow((y[ij1k]-y[ij0k]),2.))/2.;
			}
			else
			{
				deltaY_w=sqrt(pow((x[ij1k]-x[ij0k]),2.)+pow((y[ij1k]-y[ij0k]),2.))/2.;
				deltaY=sqrt(pow((x[ijk]-x[ij0k]),2.)+pow((y[ijk]-y[ij0k]),2.))+deltaY_w;
			}

			Velocity=sqrt(u_timeAverage2D[ijk]*u_timeAverage2D[ijk]+v_timeAverage2D[ijk]*v_timeAverage2D[ijk]);
			fprintf(file7," %le %le %le %le %le %le\n",deltaY,Velocity/VelocityFreeStream,deltaStern,delta99,re_deltaStern,VelocityFreeStream);
		
		}
		for(j=j_break;j<jMeshPoints;j++)
		{
			ijk=k*jMeshPoints*iMeshPoints+j*iMeshPoints+i;
			deltaY_w=sqrt(pow((x[ij1k]-x[ij0k]),2.)+pow((y[ij1k]-y[ij0k]),2.))/2.;
			deltaY=sqrt(pow((x[ijk]-x[ij0k]),2.)+pow((y[ijk]-y[ij0k]),2.))+deltaY_w;

			Velocity=VelocityFreeStream;
			fprintf(file7," %le %le %le %le %le %le\n",deltaY,Velocity/VelocityFreeStream,deltaStern,delta99,re_deltaStern,VelocityFreeStream);
		
		}		
		fclose(file7);	

		
	}
	
	return 0;	

}
