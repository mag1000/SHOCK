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

void average(cgsize_t zonesize_3D[3][3],float *buffer_2D,float *buffer_3D)
{
	long long int i,j,k;
	long long int ijk2D,ijk3D;
	for(i=0;i<zonesize_3D[0][0];i++)
	{
		for(j=0;j<zonesize_3D[0][1];j++)
		{
			ijk2D=0*zonesize_3D[0][1]*zonesize_3D[0][0]+j*zonesize_3D[0][0]+i;
			buffer_2D[ijk2D]=0.0;
			for(k=0;k<zonesize_3D[0][2];k++)
			{
				ijk3D=k*zonesize_3D[0][1]*zonesize_3D[0][0]+j*zonesize_3D[0][0]+i;
				buffer_2D[ijk2D]+=buffer_3D[ijk3D]/zonesize_3D[0][2];
			}
		}
	}
}

int main(int argc, char *argv[])
{
	char *ovalue = NULL;
	char *ivalue = NULL;
	int inputflag;
	char file[500];
	int c;

	opterr = 0;

	while ((c = getopt (argc, argv, "i:")) != -1)
	{
		switch (c)
		{
			case 'i':
				inputflag=1;
				ivalue = optarg;
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
	
	if(inputflag!=1)
	{
		printf("ERROR: Kein Inputfile angegeben (-i 'filename')\n");
		return 1;
	}
	else
	{
		strcpy(file,ivalue);
	}	
	
	int index_file_src;
	printf("Opening source: %s\n",file);
	CG( cg_open(file,CG_MODE_READ,&index_file_src));

	strcpy(file,replace_str(file, "3d", "2d"));
	strcpy(file,replace_str(file, "3D", "2D"));
	if (strcmp(file,ivalue)==0)
	{
		printf("Error: Source-Datei und Destination-Datei sind identisch. Source Datei muss 3d oder 3D im Namen enthalten!\n");
		return 0;
	}
	
	int index_file_dest;
	printf("Opening destination: %s\n",file);
	CG( cg_open(file,CG_MODE_WRITE,&index_file_dest));	
	
	char basename[33];
	sprintf(basename,"Base");
	int index_base;
	CG( cg_base_write(index_file_dest,basename,2,2,&index_base));
	
	char zonename[33];
	cgsize_t zonesize_2D[3][2];
	cgsize_t zonesize_3D[3][3];
   	int index_zone;

	printf("Reading zone...\n");
	CG( cg_zone_read(index_file_src,1,1,zonename,zonesize_3D[0]));
	

	/* vertex size */
	   zonesize_2D[0][0]=zonesize_3D[0][0];
	   zonesize_2D[0][1]=zonesize_3D[0][1];

	/* cell size */
	   zonesize_2D[1][0]=zonesize_3D[1][0];
	   zonesize_2D[1][1]=zonesize_3D[1][1];

	/* boundary vertex size (always zero for structured grids) */
	   zonesize_2D[2][0]=0;
	   zonesize_2D[2][1]=0;

	printf("Writing zone: %s...\n",zonename);
	CG( cg_zone_write(index_file_dest,1,zonename,*zonesize_2D,Structured,&index_zone));
	
	cgsize_t irmin_3D[3];
	cgsize_t irmax_3D[3];	
	cgsize_t irmin_2D[3];
	cgsize_t irmax_2D[3];
	
	irmin_2D[0]=1;
	irmin_2D[1]=1;
	irmax_2D[0]=zonesize_3D[0][0];
	irmax_2D[1]=zonesize_3D[0][1];	

	irmin_3D[0]=1;
	irmin_3D[1]=1;
	irmin_3D[2]=1;
	irmax_3D[0]=zonesize_3D[0][0];
	irmax_3D[1]=zonesize_3D[0][1];
	irmax_3D[2]=zonesize_3D[0][2];	
	
	float *buffer_3D;
	float *buffer_2D;

	buffer_2D=(float*)calloc(irmax_2D[0]*irmax_2D[1],sizeof(float));
	buffer_3D=(float*)calloc(irmax_3D[0]*irmax_3D[1]*irmax_3D[2],sizeof(float));

	int index_coord;
	printf("Transfering CoordinateX (%dx%d)...\n",(int)irmax_2D[0],(int)irmax_2D[1]);
	cg_coord_read(index_file_src,1,1,"CoordinateX",RealSingle,irmin_3D,irmax_3D,buffer_3D);
	average(zonesize_3D,buffer_2D,buffer_3D);
	cg_coord_write(index_file_dest,1,1,RealSingle,"CoordinateX",buffer_2D,&index_coord);
	
	printf("Transfering CoordinateY (%dx%d)...\n",(int)irmax_2D[0],(int)irmax_2D[1]);
	cg_coord_read(index_file_src,1,1,"CoordinateY",RealSingle,irmin_3D,irmax_3D,buffer_3D);
	average(zonesize_3D,buffer_2D,buffer_3D);
	cg_coord_write(index_file_dest,1,1,RealSingle,"CoordinateY",buffer_2D,&index_coord);
	
	int NumberSolutions;
	CG( cg_nsols(index_file_src,1,1,&NumberSolutions ) );
	int index_solution;
	int index_sol;
	int index_field;
	//char solname[NumberSolutions*32+1];
	//char solnameArray[NumberSolutions][32];
	
	//strcpy(solname,"");
	long long int i;
	
	float time[NumberSolutions];
	printf("Reading time-array...\n");
	CG( cg_goto(index_file_src,1,"BaseIterativeData_t",1,"end"));
	CG( cg_array_read_as(1,RealSingle,&time));

	GridLocation_t location; //Vertex, CellCenter, IFaceCenter, JFaceCenter, and KFaceCenter. 
	
	for(i=0;i<NumberSolutions;i++)
	{
		//sprintf(solnameArray[i],"FlowSolution%d",i);
		//sprintf(solname,"%s%-32s",solname,solnameArray[i]);
		//time[i]=L/(M*sqrt(1.4*287.*300.))*TAU*ITERBETWSAMPLES*i;
	}	
	char solname[32];
	char* solnames=NULL;
	solnames = malloc( 32*NumberSolutions );
	int field,nfields;
	char fieldname[ 33 ];
	DataType_t datatype;
	printf("Start with transfereing %d solutions\n",NumberSolutions);
	for(index_solution=1;index_solution<=NumberSolutions;index_solution++)
	{
		CG( cg_sol_info(index_file_src, 1, 1, index_solution, solname, &location));
		memcpy( solnames+(index_solution-1)*32,solname,32 );
		
		printf("Transforming Solution:\"%s\" (%d of %d) at time: %f\n",solname,index_solution,NumberSolutions,time[index_solution]);
		CG( cg_sol_write( index_file_dest,1,1,solname,location,&index_sol) );
		
		CG( cg_nfields( index_file_src,1,1,index_solution,&nfields ) );
		for( field = 1; field<=nfields; field++ ) 
		{
			cg_field_info(index_file_src,1,1,index_solution, field, &datatype, fieldname);
			printf("%s...\n",fieldname);
			cg_field_read(index_file_src,1,1,index_solution,fieldname,RealSingle,irmin_3D,irmax_3D,buffer_3D);
			average(zonesize_3D,buffer_2D,buffer_3D);
			cg_field_write(index_file_dest,1,1,index_sol,RealSingle,fieldname,buffer_2D,&index_field);
		}
	}
	
	
	printf("Check Solnames: %s\n",solnames);
	printf("Writing new baseiter,zoneiter...\n");
	CG( cg_biter_write(index_file_dest,1,"TimeIterValues",NumberSolutions));
	CG( cg_goto(index_file_dest,1,"BaseIterativeData_t",1,"end"));
	cgsize_t nuse=NumberSolutions;
	CG( cg_array_write("TimeValues",RealSingle,1,&nuse,&time));
	cgsize_t idata[2];
	idata[0]=32;
	idata[1]=NumberSolutions;
	CG( cg_ziter_write(index_file_dest,1,1,"ZoneIterativeData"));
	CG( cg_goto(index_file_dest,1,"Zone_t",1,"ZoneIterativeData_t",1,"end"));
	CG( cg_array_write("FlowSolutionPointers",Character,2,idata,solnames));	
	CG (cg_close(index_file_dest));
	CG (cg_close(index_file_src));
}
