#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include <dirent.h>
#include <getopt.h>


int main(int argc, char *argv[])
{
	printf("Merging Velocity Histories from folder \"/VelocityHistories\" \n");
	char command[500];
	int j,k;
	char temp[256];
	int j_zahl,k_zahl;
	FILE *inputfile;
	FILE *outputfile;
	int actualSampleNumber;
	size_t len = 0;
	double time_dummy;
	double * time_out;
	char out_file[300];	
	int t,t2,f,p,i,total_t;	
	int totalNumberSamples;
	int ende;
	int start;
	char actualFile[500];
	char *dummy1;
	char *dummy2;
	char *dummy3;
	char **header;
	int samplesUser;
	int NoProbes;
	int c;
	int sflag=0;
	int nflag=0;
	int zflag=0;
	int jflag=0;
	int jump=1;
	char zonename[256];	
	char filenames[256][256];
	DIR *d;
	struct dirent *dir;
	int count = 0;
	
	
	d = opendir("./VelocityHistories/");
	if (d)
	{
		while ((dir = readdir(d)) != NULL)
		{
			if((strcmp(dir->d_name,".")!=0)&&(strcmp(dir->d_name,"..")!=0))
			{
			//printf("%s\n", dir->d_name);
			strcpy(filenames[count],dir->d_name);
			count++;
			}
		}

		closedir(d);
	}
	
//	Sortiere nach Groeße
	
	for ( k = 0; k < count; k++ )
	{
	    for ( j = 0; j < k; j++ )
	    {
		sscanf(filenames[j], "VelocityHistory_Points_Iteration%d.dat", &j_zahl);
		sscanf(filenames[k], "VelocityHistory_Points_Iteration%d.dat", &k_zahl);
		if ( k_zahl < j_zahl )
		{
		    strcpy(temp,filenames[k]);
		    strcpy(filenames[k],filenames[j]);
		    strcpy(filenames[j],temp);
		}
	    }
	}

	printf("The following files are found...\n");
	for ( k = 0; k < count; k++ )
	{
	printf("%s\r\n",filenames[k]);
	}
	

	opterr = 0;
	while ((c = getopt (argc, argv, "s:n:z:j:")) != -1)
	{
		switch (c)
		{
			case 'j':
				jump = atoi(optarg);
				jflag=1;
				break;		
			case 's':
				samplesUser = atoi(optarg);
				sflag=1;
				break;
			case 'n':
				NoProbes = atoi(optarg);
				nflag=1;
				break;				
			case 'z':
				strcpy(zonename,optarg);
				zflag=1;
				break;								
			return 1;
		}
	}
	if(zflag!=1)	
	{
		strcpy(zonename,"Merged_Velocity_History");
		printf("No zonename is set (-z 'zonename'). '%s' is used.\n",zonename);
	}
	else
	{
		printf("Zonename is set to %s.\n",zonename);
	}
	

	
	if(nflag!=1)
	{
		printf("ERROR: The needed Number of Probes is not set (-n X)!\n");
		return(1);
	}
	
	if(sflag!=1)
	{
		printf("No samplesUser set (-s X). All samples are taken.\n");
		samplesUser=0;
	}

	printf("The Velocity histories of %d probes between \n%s\nand \n%s are used.\n",NoProbes,filenames[0],filenames[count-1]);

	double * VelocityX[NoProbes];	
	double   VelocityX_avrg[NoProbes];		
	double * VelocityY[NoProbes];	
	double   VelocityY_avrg[NoProbes];		
	double * VelocityZ[NoProbes];	
	double   VelocityZ_avrg[NoProbes];	

	header = (char **)calloc(NoProbes, sizeof(char*));
	sscanf(filenames[count-1], "VelocityHistory_Points_Iteration%d.dat", &ende);
	sscanf(filenames[0], "VelocityHistory_Points_Iteration%d.dat", &start);
	
	sprintf(actualFile,"VelocityHistories/%s",filenames[0]);
	inputfile=fopen(actualFile,"r");
	dummy1=NULL;
	dummy2=NULL;
	dummy3=NULL;
	getline(&dummy1,&len,inputfile);
	free(dummy1);
	getline(&dummy2,&len,inputfile);	
	free(dummy2);
	getline(&dummy3,&len,inputfile);
	sscanf(dummy3, "%*s %*s %*s %*s %*s %*[^=]=%d,",&actualSampleNumber);
	totalNumberSamples=ende-start+actualSampleNumber;
	printf("Samples found: %d\n",totalNumberSamples);
	printf("Samples reduced: %d\n",totalNumberSamples/jump);
	if(samplesUser==0)
		samplesUser=totalNumberSamples/jump;
		
	printf("Samples taken: %d\n",samplesUser);
	free(dummy3);
	fclose(inputfile);

	total_t=0;

	strcpy(out_file,"VelocityHistory_Points_All.dat");


	time_out=(double *)calloc(totalNumberSamples,sizeof(double));




	for(p=0;p<NoProbes;p++){VelocityX[p]=(double *)calloc(totalNumberSamples,sizeof(double));}
	for(p=0;p<NoProbes;p++){VelocityY[p]=(double *)calloc(totalNumberSamples,sizeof(double));}
	for(p=0;p<NoProbes;p++){VelocityZ[p]=(double *)calloc(totalNumberSamples,sizeof(double));}


	for(f=0;f<count;f++)
	{
		//sscanf(filenames[f], "VelocityHistory_Points_Iteration%d.dat", &ende);
		

		//else{start=0;}
		//actualSampleNumber=ende-start;

		sprintf(actualFile,"VelocityHistories/%s",filenames[f]);
		printf("reading %s...\n",actualFile);

		inputfile=fopen(actualFile,"r");
		dummy1=NULL;
		dummy2=NULL;
		getline(&dummy1,&len,inputfile);
		getline(&dummy2,&len,inputfile);
		
 		//printf("d1 %s\n",dummy1);
 		//printf("d2 %s\n",dummy2);	
 				

		for(p=0;p<NoProbes;p++)
		{
			getline(&header[p],&len,inputfile);
			sscanf(header[p], "%*s %*s %*s %*s %*s %*[^=]=%d,",&actualSampleNumber);
			if(p==0)
				printf("enthält %d Samples.\n",actualSampleNumber);
			for(t=0;t<actualSampleNumber;t++)
			{
				if(p==0)
				{
					fscanf(inputfile,"%le %le %le %le\n",&time_out[total_t+t],&VelocityX[p][total_t+t],&VelocityY[p][total_t+t],&VelocityZ[p][total_t+t]);
				}
				else
				{
					fscanf(inputfile,"%le %le %le %le\n",&time_dummy,&VelocityX[p][total_t+t],&VelocityY[p][total_t+t],&VelocityZ[p][total_t+t]);
				}
			}
		}		
		total_t=total_t+actualSampleNumber;

		fclose(inputfile);
		//printf("file:%d-%d\n",f,total_t);
		free(dummy1);
		free(dummy2);
	}

	printf("writing ...\n");
	outputfile=fopen(out_file,"w");
	fprintf(outputfile,"TITLE = \"%s\"\n",zonename);
	fprintf(outputfile,"VARIABLES = \"t [s]\"");
	float probeX;
	float probeY;
	float probeZ;
	for(p=0;p<NoProbes;p++)
	{
		//printf("%s\n",header[p]);
		sscanf(header[p], "%*s %*[^=]=%*[^=]=%f,",&probeX);
		sscanf(header[p], "%*s %*s %*[^=]=%f,",&probeY);
		sscanf(header[p], "%*s %*s %*s %*[^=]=%f,",&probeZ);
		//printf("(%s)-(%s):%f,%f,%f\n",dummy1,dummy2,probeX,probeY,probeZ);
		fprintf(outputfile," \"VelocityX_Probe_%d(x=%.2f,y=%.2f,z=%.2f)\"",p,probeX,probeY,probeZ);
	}
	for(p=0;p<NoProbes;p++)
	{
		//printf("%s\n",header[p]);
		sscanf(header[p], "%*s %*[^=]=%*[^=]=%f,",&probeX);
		sscanf(header[p], "%*s %*s %*[^=]=%f,",&probeY);
		sscanf(header[p], "%*s %*s %*s %*[^=]=%f,",&probeZ);
		//printf("(%s)-(%s):%f,%f,%f\n",dummy1,dummy2,probeX,probeY,probeZ);
		fprintf(outputfile," \"VelocityY_Probe_%d(x=%.2f,y=%.2f,z=%.2f)\"",p,probeX,probeY,probeZ);
	}
	for(p=0;p<NoProbes;p++)
	{
		//printf("%s\n",header[p]);
		sscanf(header[p], "%*s %*[^=]=%*[^=]=%f,",&probeX);
		sscanf(header[p], "%*s %*s %*[^=]=%f,",&probeY);
		sscanf(header[p], "%*s %*s %*s %*[^=]=%f,",&probeZ);
		//printf("(%s)-(%s):%f,%f,%f\n",dummy1,dummy2,probeX,probeY,probeZ);
		fprintf(outputfile," \"VelocityZ_Probe_%d(x=%.2f,y=%.2f,z=%.2f)\"",p,probeX,probeY,probeZ);
	}		
	fprintf(outputfile,"\nZONE T=\"%s\", F=BLOCK, I=%d, DT=(SINGLE)\n",zonename,samplesUser);

	printf("writing time...\n");
	int indexUser;
	for(indexUser=0;indexUser<samplesUser;indexUser++)
	{
	total_t=sflag*(totalNumberSamples-jump*samplesUser-1)+indexUser*jump;
	fprintf(outputfile,"%le\n",time_out[total_t]);
	}
	printf("writing Velocity histories...\n");
	for(p=0;p<NoProbes;p++)
	{
		for(indexUser=0;indexUser<samplesUser;indexUser++)
		{
			total_t=sflag*(totalNumberSamples-jump*samplesUser-1)+indexUser*jump;
			fprintf(outputfile,"%le\n",VelocityX[p][total_t]);
		}
	}
	for(p=0;p<NoProbes;p++)
	{
		for(indexUser=0;indexUser<samplesUser;indexUser++)
		{
			total_t=sflag*(totalNumberSamples-jump*samplesUser-1)+indexUser*jump;
			fprintf(outputfile,"%le\n",VelocityY[p][total_t]);
		}
	}
	for(p=0;p<NoProbes;p++)
	{
		for(indexUser=0;indexUser<samplesUser;indexUser++)
		{
			total_t=sflag*(totalNumberSamples-jump*samplesUser-1)+indexUser*jump;
			fprintf(outputfile,"%le\n",VelocityZ[p][total_t]);
		}
	}		
	fclose(outputfile);
	
	sprintf(command,"/usr/local/Tecplot2015R1/bin/preplot %s",out_file);
	printf("Executing /usr/local/Tecplot2015R1/bin/preplot...\n");
	printf("%s\n",command);		
	system(command);
	
	
	
//////////////////////////////////////////////////////////	TKE
	printf("averaging velocites ...\n");
	for(p=0;p<NoProbes;p++)
	{
		VelocityX_avrg[p]=0.0;
		VelocityY_avrg[p]=0.0;
		VelocityZ_avrg[p]=0.0;
		for(total_t=0;total_t<totalNumberSamples;total_t++)
		{
			VelocityX_avrg[p]+=VelocityX[p][t]/totalNumberSamples;
			VelocityY_avrg[p]+=VelocityY[p][t]/totalNumberSamples;
			VelocityZ_avrg[p]+=VelocityZ[p][t]/totalNumberSamples;
		}
	}

	strcpy(out_file,"TKEHistory_Points_All.dat");	
	printf("writing ...\n");
	outputfile=fopen(out_file,"w");
	fprintf(outputfile,"TITLE = \"%s\"\n",zonename);
	fprintf(outputfile,"VARIABLES = \"t [s]\"");
	for(p=0;p<NoProbes;p++)
	{
		//printf("%s\n",header[p]);
		sscanf(header[p], "%*s %*[^=]=%*[^=]=%f,",&probeX);
		sscanf(header[p], "%*s %*s %*[^=]=%f,",&probeY);
		sscanf(header[p], "%*s %*s %*s %*[^=]=%f,",&probeZ);
		//printf("(%s)-(%s):%f,%f,%f\n",dummy1,dummy2,probeX,probeY,probeZ);
		fprintf(outputfile," \"TKE_%d_x%le_y%le_z%le\"",p,probeX,probeY,probeZ);
	}
		
	fprintf(outputfile,"\nZONE T=\"%s\", F=BLOCK, I=%d, DT=(SINGLE)\n",zonename,samplesUser);

	printf("writing time...\n");
	for(indexUser=0;indexUser<samplesUser;indexUser++)
	{
	total_t=sflag*(totalNumberSamples-jump*samplesUser-1)+indexUser*jump;
	fprintf(outputfile,"%le\n",time_out[total_t]);
	}
	printf("writing Velocity histories...\n");
	double tke;
	for(p=0;p<NoProbes;p++)
	{
		for(indexUser=0;indexUser<samplesUser;indexUser++)
		{
			total_t=sflag*(totalNumberSamples-jump*samplesUser-1)+indexUser*jump;
			tke=0.5*(
			(VelocityX[p][total_t]-VelocityX_avrg[p])*(VelocityX[p][total_t]-VelocityX_avrg[p])+
			(VelocityY[p][total_t]-VelocityY_avrg[p])*(VelocityY[p][total_t]-VelocityY_avrg[p])+
			(VelocityZ[p][total_t]-VelocityZ_avrg[p])*(VelocityZ[p][total_t]-VelocityZ_avrg[p])
			);
			fprintf(outputfile,"%le\n",tke);
		}
	}
	fclose(outputfile);
	
	/*
	sprintf(command,"sed -i 's/I=%d/I=%d/g' %s",actualSampleNumber,totalNumberSamples,out_file);
	printf("Executing sed...\n");
	printf("%s\n",command);
	system(command);*/
	sprintf(command,"/usr/local/Tecplot2015R1/bin/preplot %s",out_file);
	printf("Executing /usr/local/Tecplot2015R1/bin/preplot...\n");
	printf("%s\n",command);		
	system(command);	
}
