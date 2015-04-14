#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "unistd.h"
#include "FUNCTIONS.h"
#define FS 1.25

//calcHebelarm
float cH(float f, float fm1, float hebelarm)
{
	float result;
	result=fm1+hebelarm*(f-fm1);
	return result;
}


float calcGCI(float r, float p, float f1, float f2, int version)
{
	float epsilon;
	
	float rp;
	rp=pow(r,p);
	float gci;
	if (version == 0){epsilon=fabs((f1-f2)/f1)*100.;gci=fabs((rp*FS*epsilon)/(1.0-rp));}	//grob
	if (version == 1){epsilon=fabs((f1-f2)/f2)*100.;gci=fabs((FS*epsilon)/(1.0-rp));}	//fein
	
	return gci;

}

int main(int argc, char *argv[])
{

////////////////////////////////////////////////////////////////////////////////
//////////////7 EINSTELLUNG DER PARAMETER FÜR DIE GCI AUSWERTUNG ///////////////
////////////////////////////////////////////////////////////////////////////////
	char folders[5][500];
	int iMeshPoints[5];
	int jMeshPoints[5];
	int kMeshPoints[5];
	int version;//[5][3];
	
	strcpy(folders[0],"/home/mag/WORK3/3D-GKS/512x128x64/Run2/");	iMeshPoints[0]=512;	jMeshPoints[0]=128;	kMeshPoints[0]=64;
	strcpy(folders[1],"/home/mag/WORK3/3D-GKS/1024x128x128/Run2/");	iMeshPoints[1]=1024;	jMeshPoints[1]=128;	kMeshPoints[1]=128;
	strcpy(folders[2],"/home/mag/WORK3/3D-GKS/2048x256x128/Run2/");	iMeshPoints[2]=2048;	jMeshPoints[2]=256;	kMeshPoints[2]=128;
	strcpy(folders[3],"/home/mag/WORK3/3D-GKS/4096x256x256/Run2/");	iMeshPoints[3]=4096;	jMeshPoints[3]=256;	kMeshPoints[3]=256;
	strcpy(folders[4],"/home/mag/WORK3/3D-GKS/8192x512x256/Run1/");	iMeshPoints[4]=8192;	jMeshPoints[4]=512;	kMeshPoints[4]=256;

	
	//float r[5][3];
	//int reference[5][3];
//	r[0][0]=16.0;	r[0][1]=3.62;	r[0][2]=4.0;	reference[0][0]=4;	reference[0][1]=4;	reference[0][2]=4;
//	r[1][0]=8.0;	r[1][1]=3.57;	r[1][2]=2.0;	reference[1][0]=4;	reference[1][1]=4;	reference[1][2]=4;
//	r[2][0]=4.0;	r[2][1]=2.99;	r[2][2]=2.0;	reference[2][0]=4;	reference[2][1]=4;	reference[2][2]=4;
//	r[3][0]=2.0;	r[3][1]=2.95;	r[3][2]=2.0;	reference[3][0]=4;	reference[3][1]=4;	reference[3][2]=2;
//	r[4][0]=2.0;	r[4][1]=2.95;	r[4][2]=2.0;	reference[4][0]=3;	reference[4][1]=3;	reference[4][2]=2;

	//Die j-Auflösung basiert auf dem j-wer für d99 an hinterkante:
	//	512	1024	2048	4096	8192
	//j:	64	66	79	79	227
	int MeshPoints[5][3];
	MeshPoints[0][0]=512;	MeshPoints[0][1]=65;	MeshPoints[0][2]=64;
	MeshPoints[1][0]=1024;	MeshPoints[1][1]=65;	MeshPoints[1][2]=128;
	MeshPoints[2][0]=2048;	MeshPoints[2][1]=79;	MeshPoints[2][2]=128;
	MeshPoints[3][0]=4096;	MeshPoints[3][1]=79;	MeshPoints[3][2]=256;
	MeshPoints[4][0]=8192;	MeshPoints[4][1]=227;	MeshPoints[4][2]=256;	
	
	float p_gci=4.0; // Konservative Abschätzung der effektiven Ordnung : WENO=5, ZENTRAL=6, RK=4  -> konservativ: p=4
	float r_gci;
	char extension[4][10];
	sprintf(extension[0],"_x");
	sprintf(extension[1],"_y");
	sprintf(extension[2],"_z");
	sprintf(extension[3],"_total");



//test:
/*
strcpy(folders[0],"/home/mag/workspace/mapp/C-TOOLS/");	iMeshPoints[0]=8;	jMeshPoints[0]=8;	kMeshPoints[0]=8;
strcpy(folders[1],"/home/mag/workspace/mapp/C-TOOLS/");	iMeshPoints[1]=16;	jMeshPoints[1]=16;	kMeshPoints[1]=16;
strcpy(folders[2],"/home/mag/workspace/mapp/C-TOOLS/");	iMeshPoints[2]=32;	jMeshPoints[2]=32;	kMeshPoints[2]=32;
strcpy(folders[3],"/home/mag/workspace/mapp/C-TOOLS/");	iMeshPoints[3]=64;	jMeshPoints[3]=64;	kMeshPoints[3]=64;
strcpy(folders[4],"/home/mag/workspace/mapp/C-TOOLS/");	iMeshPoints[4]=128;	jMeshPoints[4]=128;	kMeshPoints[4]=128;
MeshPoints[0][0]=8;	MeshPoints[0][1]=8;	MeshPoints[0][2]=8;
MeshPoints[1][0]=16;	MeshPoints[1][1]=16;	MeshPoints[1][2]=16;
MeshPoints[2][0]=32;	MeshPoints[2][1]=32;	MeshPoints[2][2]=32;
MeshPoints[3][0]=64;	MeshPoints[3][1]=64;	MeshPoints[3][2]=64;
MeshPoints[4][0]=128;	MeshPoints[4][1]=128;	MeshPoints[4][2]=128;
p_gci=2.0;	
*/
	


	int c;
	opterr = 0;
	
	int flag_flucQuant=0;
	int flag_All=0;	
	int optionflag=0;
	int flag_aeroCoeff=0;



	while ((c = getopt (argc, argv, "aqc")) != -1)
	{
		switch (c)
		{
			case 'a':
				optionflag++;
				flag_All=1;
				break;						
			case 'c':
				optionflag++;
				flag_aeroCoeff=1;				
			case 'q':
				optionflag++;
				flag_flucQuant=1;
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
	printf("Starting postGCI...\n");
	printf("Dieses Tool dient der weiteren Auswertung von Ergebnissen, um den Grid Convergence Index (GCI) zu erstellen.\n");
	int errorflag=0;
	
	if(optionflag<1)
	{
		printf("ERROR: Keine Auswerteoption angegeben.\nAll: -a | Fluctuating Quantities: -q | \n");
		errorflag=1;
	}

	if(errorflag==1)
	{
		return 1;
	}
	

	
	int mesh;
	int otherMesh;
	for(mesh=0;mesh<5;mesh++)
	{
		printf("The used folders for mesh %dx%dx%d is: %s\n",iMeshPoints[mesh],jMeshPoints[mesh],kMeshPoints[mesh],folders[mesh]);
	}
	
	if((flag_flucQuant==1)||(flag_All==1))
	{
		float gci_uu_tmp;
		float gci_vv_tmp;
		float gci_ww_tmp;
		float gci_uv_tmp;
		float gci_tke_tmp;

		float *gci_uu;
		float *gci_vv;
		float *gci_ww;
		float *gci_uv;
		float *gci_tke;
		
		float *gci_uu_total;
		float *gci_vv_total;
		float *gci_ww_total;
		float *gci_uv_total;
		float *gci_tke_total;		
		
		char dummy[500];
		float **yMesh,**yPlus,**uu,**vv,**ww,**uv,**tke;
		float value_dummy;
		int i,ijk;
		int status;
		yPlus=calloc(5,sizeof(float*));
		yMesh=calloc(5,sizeof(float*));
		uu=calloc(5,sizeof(float*));
		vv=calloc(5,sizeof(float*));
		ww=calloc(5,sizeof(float*));
		uv=calloc(5,sizeof(float*));
		tke=calloc(5,sizeof(float*));	
		
		///////////////// Loading data for gci calculation
		
		for(mesh=0;mesh<5;mesh++)
		{

			yPlus[mesh]=calloc(jMeshPoints[mesh],sizeof(float));
			yMesh[mesh]=calloc(jMeshPoints[mesh],sizeof(float));
			uu[mesh]=calloc(jMeshPoints[mesh],sizeof(float));
			vv[mesh]=calloc(jMeshPoints[mesh],sizeof(float));
			ww[mesh]=calloc(jMeshPoints[mesh],sizeof(float));
			uv[mesh]=calloc(jMeshPoints[mesh],sizeof(float));
			tke[mesh]=calloc(jMeshPoints[mesh],sizeof(float));			

			FILE * filefq;
			char input_fq[300];
			sprintf(input_fq,"%sFluctuatingQuantities.dat",folders[mesh]);
			//sprintf(input_fq,"%sMesh_%d.dat",folders[mesh],mesh);
			printf("Beginne FluctuatingQuantities-Import von Datei: %s\n",input_fq);		
			filefq=fopen(input_fq,"r");
			fgets(dummy,499,filefq);
			fgets(dummy,499,filefq);
			for(i=0;i<jMeshPoints[mesh];i++)
			{
				fscanf(filefq," %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",&yMesh[mesh][i],&yPlus[mesh][i],&uu[mesh][i],&vv[mesh][i],&ww[mesh][i],&uv[mesh][i],
				&value_dummy,&value_dummy,&value_dummy,&value_dummy,
				&value_dummy,&value_dummy,&value_dummy,&value_dummy,&tke[mesh][i],&value_dummy);

			}
			fclose(filefq);	
			printf("Import fertig\n");
		}
		
		///////////////// Preparing Outputfile

		FILE * file_outfq;
		char outputfq[300];
		sprintf(outputfq,"GCI_Verlaeufe_FluctuatingQuantities.dat");
		printf("Outputfile:%s\n",outputfq);
		file_outfq=fopen(outputfq,"w");
		char zonename[300];

		fprintf(file_outfq,"variables = \"y\" \"GCI(<u'u'>)\" \"GCI(<v'v'>)\" \"GCI(<w'w'>)\" \"GCI(<u'v'>)\" \"GCI(TKE)\"\n");

		int ii,ii2,direction;
		float hebelarm;
		float distance;
		float distance_tmp;
		int comparisons;
		for(mesh=0;mesh<5;mesh++)
		{
			gci_uu_tmp=0.;
			gci_vv_tmp=0.;
			gci_ww_tmp=0.;
			gci_uv_tmp=0.;
			gci_tke_tmp=0.;	
		
		
			//gci_*_total stores data averaged over all compared otherMeshes: gci_total= gci_x+gci_y+gci_z 
			gci_uu_total=calloc(jMeshPoints[mesh],sizeof(float));
			gci_vv_total=calloc(jMeshPoints[mesh],sizeof(float));
			gci_ww_total=calloc(jMeshPoints[mesh],sizeof(float));
			gci_uv_total=calloc(jMeshPoints[mesh],sizeof(float));
			gci_tke_total=calloc(jMeshPoints[mesh],sizeof(float));		
		
			for(direction=0;direction<3;direction++)
			{			

				//gci_* stores data averaged over all compared otherMeshes - this gci-variable depends on the direction
				gci_uu=calloc(jMeshPoints[mesh],sizeof(float));
				gci_vv=calloc(jMeshPoints[mesh],sizeof(float));
				gci_ww=calloc(jMeshPoints[mesh],sizeof(float));
				gci_uv=calloc(jMeshPoints[mesh],sizeof(float));
				gci_tke=calloc(jMeshPoints[mesh],sizeof(float));			
			
				sprintf(zonename,"GCI%s FluctuatingQuantities %dx%dx%d",extension[direction],iMeshPoints[mesh],jMeshPoints[mesh],kMeshPoints[mesh]);
				fprintf(file_outfq,"zone t=\"%s\", i= %d, f=point \n",zonename,jMeshPoints[mesh]-3);
			
				comparisons=0;
				for(otherMesh=0;otherMesh<5;otherMesh++)
				{
					if(MeshPoints[mesh][direction]!=MeshPoints[otherMesh][direction])
					{
						comparisons++;
						if(MeshPoints[mesh][direction]<MeshPoints[otherMesh][direction])//grob
						{
							version=0;
							r_gci=(float)MeshPoints[otherMesh][direction]/MeshPoints[mesh][direction];					
						} 
						if(MeshPoints[mesh][direction]>MeshPoints[otherMesh][direction])//fein
						{
							version=1;
							r_gci=(float)MeshPoints[mesh][direction]/MeshPoints[otherMesh][direction];					
						} 
					
							for(i=0;i<jMeshPoints[mesh]-2;i++)
						{
							distance=9999999.9;
							for(ii2=1;ii2<jMeshPoints[otherMesh]-3;ii2++)
							{
								distance_tmp=fabs(yMesh[mesh][i]-yMesh[otherMesh][ii2]);
								if (distance_tmp<distance)
								{
									distance=distance_tmp;
									if(yMesh[mesh][i]>yMesh[otherMesh][ii2])
									{
										ii=ii2+1;
									}
									else
									{
										ii=ii2;
									}
								}
							}
							hebelarm=(yMesh[mesh][i]-yMesh[otherMesh][ii-1])/(yMesh[otherMesh][ii]-yMesh[otherMesh][ii-1]);

							if((yMesh[mesh][i]-cH(yMesh[otherMesh][ii],yMesh[otherMesh][ii-1],hebelarm))/yMesh[mesh][i]>0.1)
							{
								printf(RED "error bei mesh=%d: nicht richtige interpolation: mesh:%f ref:%f\n" RESET,
								mesh,yMesh[mesh][i],cH(yMesh[otherMesh][ii],yMesh[otherMesh][ii-1],hebelarm));
							}

							gci_uu[i]+=calcGCI(r_gci, p_gci, cH(uu[otherMesh][ii],uu[otherMesh][ii-1],hebelarm), uu[mesh][i],version);
							gci_vv[i]+=calcGCI(r_gci, p_gci, cH(vv[otherMesh][ii],vv[otherMesh][ii-1],hebelarm), vv[mesh][i],version); 
							gci_ww[i]+=calcGCI(r_gci, p_gci, cH(ww[otherMesh][ii],ww[otherMesh][ii-1],hebelarm), ww[mesh][i],version);
							gci_uv[i]+=calcGCI(r_gci, p_gci, cH(uv[otherMesh][ii],uv[otherMesh][ii-1],hebelarm), uv[mesh][i],version);
							gci_tke[i]+=calcGCI(r_gci, p_gci, cH(tke[otherMesh][ii],tke[otherMesh][ii-1],hebelarm), tke[mesh][i],version);

							//ii--;
						}
						
						i=0;
						gci_tke_tmp=0.0;
						do
						{
							gci_tke_tmp+=gci_tke[i];
							i++;
						}while(yMesh[mesh][i]<100.);
		
						gci_tke_tmp=gci_tke_tmp/i;						
						
						printf("\tMesh: %d (%d points) - OtherMesh %d (%d points) - Direction:%d - Version: %d - r_gci: %f - gci(tke): %f\n",
						mesh,MeshPoints[mesh][direction],otherMesh,MeshPoints[otherMesh][direction],direction,version,r_gci,gci_tke_tmp);

					}
				}
				
				printf("Durchgeführte Vergleiche: %d\n",comparisons);
			
				for(i=0;i<jMeshPoints[mesh]-2;i++)
				{
					fprintf(file_outfq," %f %f %f %f %f %f\n",
					yMesh[mesh][i],gci_uu[i]/comparisons,gci_vv[i]/comparisons,gci_ww[i]/comparisons,gci_uv[i]/comparisons,gci_tke[i]/comparisons);
				
					gci_uu_total[i]+=gci_uu[i]/comparisons;
					gci_vv_total[i]+=gci_vv[i]/comparisons;
					gci_ww_total[i]+=gci_ww[i]/comparisons;
					gci_uv_total[i]+=gci_uv[i]/comparisons;
					gci_tke_total[i]+=gci_tke[i]/comparisons;
				}
			
			

				free(gci_uu);
				free(gci_vv);
				free(gci_ww);
				free(gci_uv);
				free(gci_tke);

			}
		
			direction=3;
			sprintf(zonename,"GCI%s %dx%dx%d",extension[direction],iMeshPoints[mesh],jMeshPoints[mesh],kMeshPoints[mesh]);
			fprintf(file_outfq,"zone t=\"%s\", i= %d, f=point \n",zonename,jMeshPoints[mesh]-3);
			for(i=0;i<jMeshPoints[mesh]-2;i++)
			{
				fprintf(file_outfq," %f %f %f %f %f %f\n",yMesh[mesh][i],gci_uu_total[i],gci_vv_total[i],gci_ww_total[i],gci_uv_total[i],gci_tke_total[i]);
			}
		
			i=0;
			gci_uu_tmp=0.0;
			gci_vv_tmp=0.0;
			gci_ww_tmp=0.0;
			gci_uv_tmp=0.0;
			gci_tke_tmp=0.0;
			do
			{
				gci_uu_tmp+=gci_uu_total[i];
				gci_vv_tmp+=gci_vv_total[i];
				gci_ww_tmp+=gci_ww_total[i];
				gci_uv_tmp+=gci_uv_total[i];
				gci_tke_tmp+=gci_tke_total[i];
			
			
				i++;
			}while(yMesh[mesh][i]<0.04);
		
			gci_uu_tmp=gci_uu_tmp/i;
			gci_vv_tmp=gci_vv_tmp/i;
			gci_ww_tmp=gci_ww_tmp/i;
			gci_uv_tmp=gci_uv_tmp/i;
			gci_tke_tmp=gci_tke_tmp/i;
			
			printf(GREEN "Mittelungen für %dx%dx%d 0<yMesh<%g : uu:%g, vv:%g, ww:%g, uv:%g, tke:%g\n" RESET,iMeshPoints[mesh],jMeshPoints[mesh],kMeshPoints[mesh],yMesh[mesh][i],gci_uu_tmp,gci_vv_tmp,gci_ww_tmp,gci_uv_tmp,gci_tke_tmp);			
		
			free(gci_uu_total);
			free(gci_vv_total);
			free(gci_ww_total);
			free(gci_uv_total);
			free(gci_tke_total);		
		}
		
		fclose(file_outfq);
		

	}
	
	if((flag_aeroCoeff==1)||(flag_All==1))
	{
		float gci_cp_tmp;
		float gci_cf_tmp;
		float gci_cl_tmp;
		float gci_cd_tmp;

		float *gci_cp;
		float *gci_cf;
		float *gci_cl;
		float *gci_cd;
		
		float *gci_cp_total;
		float *gci_cf_total;
		float *gci_cl_total;
		float *gci_cd_total;
		
		char dummy[500];
		float **xMesh,**cp,**cf,**cl,**cd;
		float value_dummy;
		int i,ijk;
		int status;
		xMesh=calloc(5,sizeof(float*));
		cp=calloc(5,sizeof(float*));
		cf=calloc(5,sizeof(float*));
		cl=calloc(5,sizeof(float*));
		cd=calloc(5,sizeof(float*));
		
		///////////////// Loading data for gci calculation
		
		for(mesh=0;mesh<5;mesh++)
		{

			xMesh[mesh]=calloc(iMeshPoints[mesh],sizeof(float));
			cp[mesh]=calloc(iMeshPoints[mesh],sizeof(float));
			cf[mesh]=calloc(iMeshPoints[mesh],sizeof(float));
			cl[mesh]=calloc(iMeshPoints[mesh],sizeof(float));
			cd[mesh]=calloc(iMeshPoints[mesh],sizeof(float));

			FILE * fileac;
			char input_ac[300];
			sprintf(input_ac,"%saerodynmicCoefficients.dat",folders[mesh]);
			//sprintf(input_ac,"%sMesh_%d.dat",folders[mesh],mesh);
			printf("Beginne AerodynamicCoefficients-Import von Datei: %s\n",input_ac);		
			fileac=fopen(input_ac,"r");
			fgets(dummy,499,fileac);
			fgets(dummy,499,fileac);
			for(i=0;i<iMeshPoints[mesh];i++)
			{
				fscanf(fileac," %f %f %f %f %f %f %f %f %f\n",
				&xMesh[mesh][i],&cp[mesh][i],&value_dummy,&cf[mesh][i],&value_dummy,&cl[mesh][i],&value_dummy,&cd[mesh][i],&value_dummy);

			}
			fclose(fileac);	
			printf("Import fertig\n");
		}
		
		///////////////// Preparing Outputfile

		FILE * file_outac;
		char outputac[300];
		sprintf(outputac,"GCI_Verlaeufe_aerodynamicCoefficients.dat");
		printf("Outputfile:%s\n",outputac);
		file_outac=fopen(outputac,"w");
		char zonename[300];

		fprintf(file_outac,"variables = \"x\" \"GCI(c<sub>p</sub>)\" \"GCI(c<sub>f</sub>)\" \"GCI(c<sub>l</sub>)\" \"GCI(c<sub>d</sub>)\"\n");

		int ii,ii2,direction;
		float hebelarm;
		float distance;
		float distance_tmp;
		int comparisons;
		int counter;
		for(mesh=0;mesh<5;mesh++)
		{
			gci_cp_tmp=0.;
			gci_cf_tmp=0.;
			gci_cl_tmp=0.;
			gci_cd_tmp=0.;
		
		
			//gci_*_total stores data averaged over all compared otherMeshes: gci_total= gci_x+gci_y+gci_z 
			gci_cp_total=calloc(iMeshPoints[mesh],sizeof(float));
			gci_cf_total=calloc(iMeshPoints[mesh],sizeof(float));
			gci_cl_total=calloc(iMeshPoints[mesh],sizeof(float));
			gci_cd_total=calloc(iMeshPoints[mesh],sizeof(float));
		
			for(direction=0;direction<3;direction++)
			{			

				//gci_* stores data averaged over all compared otherMeshes - this gci-variable depends on the direction
				gci_cp=calloc(iMeshPoints[mesh],sizeof(float));
				gci_cf=calloc(iMeshPoints[mesh],sizeof(float));
				gci_cl=calloc(iMeshPoints[mesh],sizeof(float));
				gci_cd=calloc(iMeshPoints[mesh],sizeof(float));
			
				sprintf(zonename,"GCI%s AerodynamicCoefficients %dx%dx%d",extension[direction],iMeshPoints[mesh],jMeshPoints[mesh],kMeshPoints[mesh]);
				fprintf(file_outac,"zone t=\"%s\", i= %d, f=point \n",zonename,iMeshPoints[mesh]-4);
			
				comparisons=0;
				for(otherMesh=0;otherMesh<5;otherMesh++)
				{
					if(MeshPoints[mesh][direction]!=MeshPoints[otherMesh][direction])
					{
						comparisons++;
						if(MeshPoints[mesh][direction]<MeshPoints[otherMesh][direction])//grob
						{
							version=0;
							r_gci=(float)MeshPoints[otherMesh][direction]/MeshPoints[mesh][direction];					
						} 
						if(MeshPoints[mesh][direction]>MeshPoints[otherMesh][direction])//fein
						{
							version=1;
							r_gci=(float)MeshPoints[mesh][direction]/MeshPoints[otherMesh][direction];					
						} 
					
						for(i=1;i<iMeshPoints[mesh]-2;i++)
						{
							distance=9999999.9;
							for(ii2=1;ii2<iMeshPoints[otherMesh]-3;ii2++)
							{
								distance_tmp=fabs(xMesh[mesh][i]-xMesh[otherMesh][ii2]);
								if ((distance_tmp<distance)&&(cl[mesh][i]*cl[otherMesh][ii2]>0.0))
								{
									distance=distance_tmp;
									if(xMesh[mesh][i]>xMesh[otherMesh][ii2])
									{
										ii=ii2+1;
									}
									else
									{
										ii=ii2;
									}
								}
							}
							hebelarm=(xMesh[mesh][i]-xMesh[otherMesh][ii-1])/(xMesh[otherMesh][ii]-xMesh[otherMesh][ii-1]);

							if((xMesh[mesh][i]-cH(xMesh[otherMesh][ii],xMesh[otherMesh][ii-1],hebelarm))/xMesh[mesh][i]>0.1)
							{
								printf(RED "error bei mesh=%d: nicht richtige interpolation: mesh:%f ref:%f\n" RESET,
								mesh,xMesh[mesh][i],cH(xMesh[otherMesh][ii],xMesh[otherMesh][ii-1],hebelarm));
							}
							
							/*printf("Vergleiche bei: mesh:x=%f,cp=%f - ref:x=%f,cp=%f mit hebelarm:x=%f,cp=%f ->gci:%f\n",
							xMesh[mesh][i],cp[mesh][i],xMesh[otherMesh][ii],cp[otherMesh][ii],
							cH(xMesh[otherMesh][ii],xMesh[otherMesh][ii-1],hebelarm),
							cH(cp[otherMesh][ii],cp[otherMesh][ii-1],hebelarm),
							calcGCI(r_gci, p_gci, cH(cp[otherMesh][ii],cp[otherMesh][ii-1],hebelarm), cp[mesh][i],version));*/
							
							
							gci_cp[i]+=calcGCI(r_gci, p_gci, cH(cp[otherMesh][ii],cp[otherMesh][ii-1],hebelarm), cp[mesh][i],version);
							gci_cf[i]+=calcGCI(r_gci, p_gci, cH(cf[otherMesh][ii],cf[otherMesh][ii-1],hebelarm), cf[mesh][i],version); 
							gci_cl[i]+=calcGCI(r_gci, p_gci, cH(cl[otherMesh][ii],cl[otherMesh][ii-1],hebelarm), cl[mesh][i],version);
							gci_cd[i]+=calcGCI(r_gci, p_gci, cH(cd[otherMesh][ii],cd[otherMesh][ii-1],hebelarm), cd[mesh][i],version);
						}
						
						counter=0;
						gci_cp_tmp=0.0;
						for(i=0;i<iMeshPoints[mesh];i++)
						{
							if(xMesh[mesh][i]<=1)
							{
								gci_cp_tmp+=gci_cp[i];
								counter++;
							}
						}
						gci_cp_tmp=gci_cp_tmp/counter;						
						
						printf("\tMesh: %d (%d points) - OtherMesh %d (%d points) - Direction:%d - Version: %d - r_gci: %f - gci(cp): %f\n",
						mesh,MeshPoints[mesh][direction],otherMesh,MeshPoints[otherMesh][direction],direction,version,r_gci,gci_cp_tmp);

					}
				}
				
				printf("Durchgeführte Vergleiche: %d\n",comparisons);
			
				for(i=1;i<iMeshPoints[mesh]-2;i++)
				{
					fprintf(file_outac," %f %f %f %f %f\n",
					xMesh[mesh][i],gci_cp[i]/comparisons,gci_cf[i]/comparisons,gci_cl[i]/comparisons,gci_cd[i]/comparisons);
				
					gci_cp_total[i]+=gci_cp[i]/comparisons;
					gci_cf_total[i]+=gci_cf[i]/comparisons;
					gci_cl_total[i]+=gci_cl[i]/comparisons;
					gci_cd_total[i]+=gci_cd[i]/comparisons;
				}
			
			

				free(gci_cp);
				free(gci_cf);
				free(gci_cl);
				free(gci_cd);

			}
		
			direction=3;
			sprintf(zonename,"GCI%s %dx%dx%d",extension[direction],iMeshPoints[mesh],jMeshPoints[mesh],kMeshPoints[mesh]);
			fprintf(file_outac,"zone t=\"%s\", i= %d, f=point \n",zonename,iMeshPoints[mesh]-4);
			for(i=0;i<iMeshPoints[mesh]-2;i++)
			{
				fprintf(file_outac," %f %f %f %f %f\n",xMesh[mesh][i],gci_cp_total[i],gci_cf_total[i],gci_cl_total[i],gci_cd_total[i]);
			}
		
			counter=0;
			gci_cp_tmp=0.0;
			gci_cf_tmp=0.0;
			gci_cl_tmp=0.0;
			gci_cd_tmp=0.0;
			for(i=0;i<iMeshPoints[mesh];i++)
			{
				if(xMesh[mesh][i]<=1)
				{
					gci_cp_tmp+=gci_cp_total[i];
					gci_cf_tmp+=gci_cf_total[i];
					gci_cl_tmp+=gci_cl_total[i];
					gci_cd_tmp+=gci_cd_total[i];
					counter++;
				}				
			}
		
			gci_cp_tmp=gci_cp_tmp/counter;
			gci_cf_tmp=gci_cf_tmp/counter;
			gci_cl_tmp=gci_cl_tmp/counter;
			gci_cd_tmp=gci_cd_tmp/counter;
			
			printf(GREEN "Mittelungen für %dx%dx%d 0<xMesh<1 : cp:%g, cf:%g, cl:%g, cd:%g\n" RESET,
			iMeshPoints[mesh],jMeshPoints[mesh],kMeshPoints[mesh],gci_cp_tmp,gci_cf_tmp,gci_cl_tmp,gci_cd_tmp);			
		
			free(gci_cp_total);
			free(gci_cf_total);
			free(gci_cl_total);
			free(gci_cd_total);
		}
		
		fclose(file_outac);
		
	}	

}


























