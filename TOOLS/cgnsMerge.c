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
	char file_1[500];
	char file_2[500];	
	
	int fileflag=0;		
	int c;
	
	opterr = 0;

	while ((c = getopt (argc, argv, "f:g:")) != -1)
	{
		switch (c)
		{
			case 'f':
				fileflag+=1;
				strcpy(file_1,optarg);
				break;						
			case 'g':
				fileflag+=1;
				strcpy(file_2,optarg);
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
	printf("Starting cgnsMerge....\n");
	printf("cgnsMerge führt zwei cgnsDateien mit identischer Zonengröße aber unterschiedlichen Lösungen (Ergebnisse aus zwei aufeinanderfolgenden Simulationen) zu einer neuen Datei zusammen.\n");
	printf("Dabei werden BoundaryConditions und Interfaces ignoriert, so dass die erstellte Datei zum Weiterrechnen unbrauchbar ist.\n");
	int errorflag=0;
	if(fileflag!=2)
	{
		printf("ERROR: Keine zwei Datein zum Zusammenführen angegeben (-f 'file1.cgns' -g 'file2.cgns')!\n");
		errorflag=1;
	}

	if(errorflag==1)
	{
		return 1;
	}
	
	char file_merge[500];
	strcpy(file_merge,replace_str(file_1, ".cgns","_merged.cgns"));

	printf("Die Dateien %s (1.) und %s (2.) werden in die Datei %s zusammengeführt.\n",
	file_1,file_2,file_merge);
	

//	Öffnen beider Dateien
	int index_file_1;
	int index_file_2;
	int index_file_merge;
	int index_base=1;
	int index_zone=1;
	int index_sol=1;
	char zonename[33];

	printf("Oeffne %s im read mode...\n",file_1);	
	CG( cg_open(file_1,CG_MODE_READ,&index_file_1));
	int cell_dim;
	CG( cg_cell_dim(index_file_1,index_base,&cell_dim));
	cgsize_t zonesize[cell_dim][cell_dim];	
	printf("Reading zone...\n");	
	CG( cg_zone_read(index_file_1,index_base,index_zone,zonename,zonesize[0]));
	
	
	printf("Oeffne %s im read mode...\n",file_2);	
	CG( cg_open(file_2,CG_MODE_READ,&index_file_2));


	printf("Oeffne %s im write mode...\n",file_merge);	
	CG( cg_open(file_merge,CG_MODE_WRITE,&index_file_merge));
	printf("Writing base...\n");
	CG( cg_base_write(index_file_merge,"Base",cell_dim,cell_dim,&index_base));
	
	printf("Writing zone...\n");
	CG( cg_zone_write(index_file_merge,index_base,zonename,*zonesize,Structured,&index_zone));
	

	cgsize_t irmin[cell_dim];
	cgsize_t irmax[cell_dim];	
	
	irmin[0]=1;
	irmin[1]=1;
	irmax[0]=zonesize[0][0];
	irmax[1]=zonesize[0][1];
	if(cell_dim==3)
	{
		irmin[2]=1;
		irmax[2]=zonesize[0][2];		
	}	

	float *buffer;

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
	cg_coord_read(index_file_1,1,1,"CoordinateX",RealSingle,irmin,irmax,buffer);
	cg_coord_write(index_file_merge,1,1,RealSingle,"CoordinateX",buffer,&index_coord);
	
	printf("Transfering CoordinateY ...\n");
	cg_coord_read(index_file_1,1,1,"CoordinateY",RealSingle,irmin,irmax,buffer);
	cg_coord_write(index_file_merge,1,1,RealSingle,"CoordinateY",buffer,&index_coord);
	
	if(cell_dim==3)
	{
		printf("Transfering CoordinateZ ...\n");
		cg_coord_read(index_file_1,1,1,"CoordinateZ",RealSingle,irmin,irmax,buffer);
		cg_coord_write(index_file_merge,1,1,RealSingle,"CoordinateZ",buffer,&index_coord);	
	}
	
	
	int NumberSolutions_1;
	int NumberSolutions_2;
	CG( cg_nsols(index_file_1,index_base,index_zone,&NumberSolutions_1 ) );
	CG( cg_nsols(index_file_2,index_base,index_zone,&NumberSolutions_2 ) );
	int NumberSolutions_merge;
	NumberSolutions_merge=NumberSolutions_1+NumberSolutions_2-1; //-1 weil die letzte Solution von 1 gleich der ersten von 2 ist

	int index_solution;
	int index_field;
	//char solname[NumberSolutions*32+1];
	//char solnameArray[NumberSolutions][32];
	
	//strcpy(solname,"");
	
	float time_1[NumberSolutions_1];
	float time_2[NumberSolutions_2];
	float time_merge[NumberSolutions_merge];

	printf("Reading time-array 1...\n");
	CG( cg_goto(index_file_1,1,"BaseIterativeData_t",1,"end"));
	CG( cg_array_read_as(1,RealSingle,&time_1));
	
	printf("Reading time-array 2...\n");
	CG( cg_goto(index_file_2,1,"BaseIterativeData_t",1,"end"));
	CG( cg_array_read_as(1,RealSingle,&time_2));
	

	for(index_solution=0;index_solution<NumberSolutions_merge;index_solution++)
	{
		if(index_solution<NumberSolutions_1)
		{
			time_merge[index_solution]=time_1[index_solution];
		}
		else
		{
			time_merge[index_solution]=time_2[index_solution-NumberSolutions_1+1];
		}
	}
	
	
	char *solnames;
	solnames = malloc( 32* NumberSolutions_merge  );
	char solname[32];
	int field,nfields;
	char fieldname[ 33 ];
	GridLocation_t location; //Vertex, CellCenter, IFaceCenter, JFaceCenter, and KFaceCenter. 
	CG( cg_sol_info(index_file_1, 1, 1, 1, solname, &location));		
	DataType_t datatype;
	printf("Start with transfering %d solutions\n",NumberSolutions_merge);
	for(index_solution=1;index_solution<=NumberSolutions_merge;index_solution++)
	{

		
		snprintf( solname,33,"Sol_%04d",index_solution-1 );
		memcpy( solnames+(index_solution-1)*32,solname,32 );
		printf("Transfering Solution %s (%d of %d) at time: %f\n",solname,index_solution,NumberSolutions_merge,time_merge[index_solution-1]);
		CG( cg_sol_write( index_file_merge,1,1,solname,location,&index_sol) );
		if(index_solution<=NumberSolutions_1)
		{
			printf("Datei 1.\n");
			CG( cg_nfields( index_file_1,1,1,index_solution,&nfields ) );
			for( field = 1; field<=nfields; field++ ) 
			{
				CG( cg_field_info(index_file_1,1,1,index_solution, field, &datatype, fieldname) );
				printf("%s...\n",fieldname);
				CG( cg_field_read(index_file_1,1,1,index_solution,fieldname,RealSingle,irmin,irmax,buffer) );
				CG( cg_field_write(index_file_merge,1,1,index_solution,RealSingle,fieldname,buffer,&index_field) );
			}
		}
		else
		{
			printf("Datei 2.\n");
			CG( cg_nfields( index_file_2,1,1,index_solution-NumberSolutions_1+1,&nfields ) );
			for( field = 1; field<=nfields; field++ ) 
			{
				CG( cg_field_info(index_file_2,1,1,index_solution-NumberSolutions_1+1, field, &datatype, fieldname) );
				printf("%s...\n",fieldname);
				CG( cg_field_read(index_file_2,1,1,index_solution-NumberSolutions_1+1,fieldname,RealSingle,irmin,irmax,buffer) );
				CG( cg_field_write(index_file_merge,1,1,index_solution,RealSingle,fieldname,buffer,&index_field) );
			}
		}

	}

	printf("Writing new baseiter,zoneiter...\n");
	CG( cg_biter_write(index_file_merge,1,"TimeIterValues",NumberSolutions_merge));
	CG( cg_goto(index_file_merge,1,"BaseIterativeData_t",1,"end"));
	cgsize_t nuse=NumberSolutions_merge;
	CG( cg_array_write("TimeValues",RealSingle,1,&nuse,&time_merge));
	cgsize_t idata[2];
	idata[0]=32;
	idata[1]=NumberSolutions_merge;
	CG( cg_ziter_write(index_file_merge,1,1,"ZoneIterativeData"));
	CG( cg_goto(index_file_merge,1,"Zone_t",1,"ZoneIterativeData_t",1,"end"));
	CG( cg_array_write("FlowSolutionPointers",Character,2,idata,solnames));			
	
	
	CG (cg_close(index_file_merge));
	CG (cg_close(index_file_1));
	CG (cg_close(index_file_2));
}

