#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"
#include "unistd.h"

#define CG(cmd) if(cmd)cg_error_exit( );

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
	
	char *ivalue = NULL;
	int inputflag;

	int fileflag=0;
	int c;
	int optionflag=0;
	char extension[100];
	strcpy(extension,".cgns");
	int NumberSamples=0;
	int flag_samples;
	opterr = 0;

	while ((c = getopt (argc, argv, "f:s:")) != -1)
	{
		switch (c)
		{
			case 'f':
				fileflag+=1;
				strcpy(file_in,optarg);
				break;
			case 's':
				flag_samples=1;
				NumberSamples=atoi(optarg);
				break;				
								
			case '?':
				printf("Geben Sie bei Bedarf eine Option an: numer of Samples (-s (INT))\n");
				if (optopt == 'o')
					fprintf (stderr, "Option -%o requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%o'.\n", optopt);
				else
					fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
			return 1;
		}
	}
			
	printf("Starting cgnsTimeAverage....\n");
	int errorflag=0;
	if(fileflag!=1)
	{
		printf("ERROR: Keine Datei angegeben (-f 'filename.cgns')!\n");
		errorflag=1;		
	}
	
	if(errorflag==1)
	{
		return 1;
	}


	printf("Öffne Datei %s im read mode\n",file_in);
	int index_file_in;
	int index_file_out;
	int index_base=1;
	CG( cg_open(file_in,CG_MODE_READ,&index_file_in));
	int cell_dim;
	CG( cg_cell_dim(index_file_in,index_base,&cell_dim));

	cgsize_t zonesize[cell_dim][cell_dim];	
	
	int number_solutions=1;
	int index_zone=1;
	char zonename[33];

	CG( cg_zone_read(index_file_in,index_base,index_zone,zonename,zonesize[0]));

	cgsize_t irmin[cell_dim];
	cgsize_t irmax[cell_dim];
	
	long long int iMeshPoints;
	long long int jMeshPoints;
	long long int kMeshPoints;		
	long long int i,j,k,ijk;
	float *buffer;
	float *buffer_avrg;
	float *buffer_varianz;
		
	if (cell_dim==2)
	{
		irmin[0]=1;
		irmin[1]=1;	
		irmax[0]=zonesize[0][0];
		irmax[1]=zonesize[0][1];
		printf("Lade Ergebnisse aus Zone %s (%zux%zu)...\n",zonename,irmax[0],irmax[1]);
		iMeshPoints=irmax[0];
		jMeshPoints=irmax[1];
		kMeshPoints=1;
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
	
	int NumberSolutions;
	CG( cg_nsols(index_file_in,1,1,&NumberSolutions ) );
	int SolutionStart;
	int SolutionEnd;
	if(NumberSamples>=1)
	{
		SolutionStart=NumberSolutions-NumberSamples+1;
		if(SolutionStart<1)
		{
			SolutionStart=1;
			NumberSamples=NumberSolutions;
		}
	}
	else
	{
		printf("Illegal value for Numer Samples. Have to be >= 1.\n");
		return 1;
	}

	SolutionEnd=NumberSolutions;
	
	int NumberSolutionsEff=SolutionEnd-SolutionStart+1;

	char replaceString[50];
	sprintf(replaceString,"_VarianzTimeAverage_%dSamples.cgns",NumberSamples);
	strcpy(extension,replace_str(extension, ".cgns",replaceString));	
		
	strcpy(file_out,replace_str(file_in, ".cgns",extension));
	printf("Öffne Datei %s im write mode\n",file_out);

	CG( cg_open(file_out,CG_MODE_WRITE,&index_file_out));
	
	printf("Writing base...\n");
	CG( cg_base_write(index_file_out,"Base",cell_dim,cell_dim,&index_base));
	
	printf("Writing zone...\n");
	CG( cg_zone_write(index_file_out,index_base,zonename,*zonesize,Structured,&index_zone));

	if(cell_dim==2)
	{
		buffer=(float*)calloc(irmax[0]*irmax[1],sizeof(float));
	}
	else
	{
		buffer=(float*)calloc(irmax[0]*irmax[1]*irmax[2],sizeof(float));
	}


	int index_coord;
	printf("Transfering CoordinateX ...\n");
	cg_coord_read(index_file_in,1,1,"CoordinateX",RealSingle,irmin,irmax,buffer);
	cg_coord_write(index_file_out,1,1,RealSingle,"CoordinateX",buffer,&index_coord);
	
	printf("Transfering CoordinateY ...\n");
	cg_coord_read(index_file_in,1,1,"CoordinateY",RealSingle,irmin,irmax,buffer);
	cg_coord_write(index_file_out,1,1,RealSingle,"CoordinateY",buffer,&index_coord);
	
	if(cell_dim==3)
	{
		printf("Transfering CoordinateZ ...\n");
		cg_coord_read(index_file_in,1,1,"CoordinateZ",RealSingle,irmin,irmax,buffer);
		cg_coord_write(index_file_out,1,1,RealSingle,"CoordinateZ",buffer,&index_coord);	
	}	
	
	free(buffer);
	

	float time[NumberSolutions];
	printf("Reading time-array 1...\n");
	CG( cg_goto(index_file_in,1,"BaseIterativeData_t",1,"end"));
	CG( cg_array_read_as(1,RealSingle,&time));
		
	printf("Das Auswertefenster beginne bei Lösung %d(%f s) und endet bei %d(%f s) - %d von %d Lösungen werden betrachtet. Physikalischer Zeitraum: %f s\n",SolutionStart,time[SolutionStart-1],SolutionEnd,time[SolutionEnd-1],NumberSolutionsEff,NumberSolutions,(time[SolutionEnd-1]-time[SolutionStart-1]));
	

	int index_solution;
	int index_field;
	//char solname[NumberSolutions*32+1];
	//char solnameArray[NumberSolutions][32];
	
	//strcpy(solname,"");

	GridLocation_t location; //Vertex, CellCenter, IFaceCenter, JFaceCenter, and KFaceCenter. 

	char solname[32];

	CG( cg_sol_info(index_file_in, 1, 1, 1, solname, &location));		
	sprintf(solname,"TimeAverage_%dSamples",NumberSamples);
	CG( cg_sol_write( index_file_out,1,1,solname,location,&index_solution) );		

	int field,nfields;
	char fieldname[ 33 ];
	char fieldname_var[ 33 ];
	char fieldname_time[ 33 ];
	DataType_t datatype;
	CG( cg_nfields( index_file_in,1,1,index_solution,&nfields ) );
	for( field = 1; field<=nfields; field++ ) 
	{
		if(cell_dim==2)
		{
			buffer=(float*)calloc(irmax[0]*irmax[1],sizeof(float));
			buffer_avrg=(float*)calloc(irmax[0]*irmax[1],sizeof(float));
			buffer_varianz=(float*)calloc(irmax[0]*irmax[1],sizeof(float));	
		}
		else
		{
			buffer=(float*)calloc(irmax[0]*irmax[1]*irmax[2],sizeof(float));
			buffer_avrg=(float*)calloc(irmax[0]*irmax[1]*irmax[2],sizeof(float));
			buffer_varianz=(float*)calloc(irmax[0]*irmax[1]*irmax[2],sizeof(float));			
		}
		cg_field_info(index_file_in,1,1,1, field, &datatype, fieldname);
		printf("############## Variable: %s ############\n",fieldname);
		
		printf("Beginne mit Auswertung von %d Lösungen.\n",NumberSolutionsEff);
		for(index_solution=SolutionStart;index_solution<=SolutionEnd;index_solution++)
		{
			printf("Reading   Solution %d of %d\r",index_solution,NumberSolutions);
			fflush(stdout);
			cg_field_read(index_file_in,1,1,index_solution,fieldname,RealSingle,irmin,irmax,buffer);
			printf("Averaging Solution %d of %d\r",index_solution,NumberSolutions);
			fflush(stdout);

			for(i=0;i<iMeshPoints*jMeshPoints*kMeshPoints;i++)
			{
				buffer_avrg[i]+=buffer[i];
				buffer_varianz[i]+=(buffer[i]*buffer[i]);
			}			
		}

		for(i=0;i<iMeshPoints*jMeshPoints*kMeshPoints;i++)
		{

			buffer_avrg[i]=buffer_avrg[i]/NumberSolutionsEff;
			buffer_varianz[i]=1./NumberSolutionsEff*(buffer_varianz[i]-1./NumberSolutionsEff*(buffer_avrg[i]*NumberSolutionsEff)*(buffer_avrg[i]*NumberSolutionsEff));
		}			
		printf("\n");

		sprintf(fieldname_time,"%s_timeAverage",fieldname);
		cg_field_write(index_file_out,1,1,1,RealSingle,fieldname_time,buffer_avrg,&index_field);

		sprintf(fieldname_var,"%s_varianz",fieldname);
		cg_field_write(index_file_out,1,1,1,RealSingle,fieldname_var,buffer_varianz,&index_field);	
		
		free(buffer);
		free(buffer_avrg);
		free(buffer_varianz);				
		
	}
	
	
	printf("Writing new baseiter,zoneiter...\n");
	CG( cg_biter_write(index_file_out,1,"TimeIterValues",1));
	CG( cg_goto(index_file_out,1,"BaseIterativeData_t",1,"end"));
	float time_out=(time[SolutionEnd-1]-time[SolutionStart-1]);
	cgsize_t nuse=1;
	CG( cg_array_write("TimeValues",RealSingle,1,&nuse,&time_out));
	cgsize_t idata[2];
	idata[0]=32;
	idata[1]=1;
	CG( cg_ziter_write(index_file_out,1,1,"ZoneIterativeData"));
	CG( cg_goto(index_file_out,1,"Zone_t",1,"ZoneIterativeData_t",1,"end"));
	CG( cg_array_write("FlowSolutionPointers",Character,2,idata,solname));	
	CG (cg_close(index_file_out));
	CG (cg_close(index_file_in));
}

