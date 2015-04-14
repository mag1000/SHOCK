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
		float p_uu_tmp;
		float p_vv_tmp;
		float p_ww_tmp;
		float p_uv_tmp;
		float p_tke_tmp;

		float *p_uu;
		float *p_vv;
		float *p_ww;
		float *p_uv;
		float *p_tke;
		
		
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
		sprintf(outputfq,"GCI_p_Verlaeufe_FluctuatingQuantities.dat");
		printf("Outputfile:%s\n",outputfq);
		file_outfq=fopen(outputfq,"w");
		char zonename[300];

		fprintf(file_outfq,"variables = \"y\" \"p(<u'u'>)\" \"p(<v'v'>)\" \"p(<w'w'>)\" \"p(<u'v'>)\" \"p(TKE)\"\n");

		int ii,ii2,direction;
		float hebelarm;
		float distance;
		float distance_tmp;
		int comparisons;
		for(mesh=1;mesh<4;mesh++)
		{
			p_uu_tmp=0.;
			p_vv_tmp=0.;
			p_ww_tmp=0.;
			p_uv_tmp=0.;
			p_tke_tmp=0.;	
		
		
			p_uu=calloc(jMeshPoints[mesh],sizeof(float));
			p_vv=calloc(jMeshPoints[mesh],sizeof(float));
			p_ww=calloc(jMeshPoints[mesh],sizeof(float));
			p_uv=calloc(jMeshPoints[mesh],sizeof(float));
			p_tke=calloc(jMeshPoints[mesh],sizeof(float));			
		

		
			direction=3;
			sprintf(zonename,"p_%d",mesh-1);
			fprintf(file_outfq,"zone t=\"%s\", i= %d, f=point \n",zonename,jMeshPoints[mesh]-3);
		

			r_gci=4;
			for(i=0;i<jMeshPoints[mesh]-2;i++)
			{
				otherMesh=mesh-1;
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

				p_uu[i]=log(fabs(cH(uu[otherMesh][ii],uu[otherMesh][ii-1],hebelarm)-uu[mesh][i]));
				p_vv[i]=log(fabs(cH(vv[otherMesh][ii],vv[otherMesh][ii-1],hebelarm)-vv[mesh][i]));
				p_ww[i]=log(fabs(cH(ww[otherMesh][ii],ww[otherMesh][ii-1],hebelarm)-ww[mesh][i]));
				p_uv[i]=log(fabs(cH(uv[otherMesh][ii],uv[otherMesh][ii-1],hebelarm)-uv[mesh][i]));
				p_tke[i]=log(fabs(cH(tke[otherMesh][ii],tke[otherMesh][ii-1],hebelarm)-tke[mesh][i]));
				
				
				
				
				
				otherMesh=mesh+1;
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

				p_uu[i]-=log(fabs(cH(uu[otherMesh][ii],uu[otherMesh][ii-1],hebelarm)-uu[mesh][i]));
				p_vv[i]-=log(fabs(cH(vv[otherMesh][ii],vv[otherMesh][ii-1],hebelarm)-vv[mesh][i]));
				p_ww[i]-=log(fabs(cH(ww[otherMesh][ii],ww[otherMesh][ii-1],hebelarm)-ww[mesh][i]));
				p_uv[i]-=log(fabs(cH(uv[otherMesh][ii],uv[otherMesh][ii-1],hebelarm)-uv[mesh][i]));
				p_tke[i]-=log(fabs(cH(tke[otherMesh][ii],tke[otherMesh][ii-1],hebelarm)-tke[mesh][i]));				


				fprintf(file_outfq," %f %f %f %f %f %f\n",yMesh[mesh][i],p_uu[i],p_vv[i],p_ww[i],p_uv[i],p_tke[i]);
			}
			
			free(p_uu);
			free(p_vv);
			free(p_ww);
			free(p_uv);
			free(p_tke);

	
		}
		
		fclose(file_outfq);
		

	}
}


























