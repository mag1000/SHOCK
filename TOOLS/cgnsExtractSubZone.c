#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"
#include "unistd.h"
#include <stdbool.h>

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

	char extension[100];

	int Solution_Start=1;
	int Solution_End=0;	
	int flag_solution=0;
	
	long long int iStart=0;
	long long int iEnd=0;
	long long int jStart=0;
	long long int jEnd=0;
	long long int kStart=0;
	long long int kEnd=0;

	opterr = 0;

	while ((c = getopt (argc, argv, "f:s:e:i:x:j:y:k:z:")) != -1)
	{
		switch (c)
		{
			case 'f':
				fileflag+=1;
				strcpy(file_in,optarg);
				break;
			case 's':
				Solution_Start=atoi(optarg);
				break;				
			case 'e':
				Solution_End=atoi(optarg);
				break;
			case 'i':
				iStart=atoi(optarg);
				break;
			case 'x':
				iEnd=atoi(optarg);
				break;
			case 'j':
				jStart=atoi(optarg);
				break;
			case 'y':
				jEnd=atoi(optarg);
				break;
			case 'k':
				kStart=atoi(optarg);
				break;
			case 'z':
				kEnd=atoi(optarg);
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
	printf("Geben Sie bei Bedarf eine Option an: (s)tart at solution 'x' - (e)nd at solution 'x' !\n");	
	printf("Starting cgnsExtractSubZone....\n");
	int errorflag=0;
	if(fileflag!=1)
	{
		printf("ERROR: Keine Datei angegeben (-f 'filename.cgns')!\n");
		errorflag=1;		
	}

	if(((iEnd-iStart)*(jEnd-jStart)*(kEnd-kStart))<2)
	{
		printf("ERROR: Ungültiger Subzone-Bereich!\n");
		printf("Geben Sie den zu extrahierenden Bereich an:\n -i (iStart) -x (iEnd)\n -j (jStart) -y (jEnd)\n -k (kStart) -z (kEnd)\n");			
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
	char solname[32];
	cgsize_t zonesize[3][cell_dim];
	cgsize_t zonesize_SubZone[3][cell_dim];		
	
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
	
	sprintf(extension,"_i%lld-%lld_j%lld-%lld_k%lld-%lld_Sol%d-%d.cgns",iStart,iEnd,jStart,jEnd,kStart,kEnd,Solution_Start,Solution_End);
	strcpy(file_out,replace_str(file_in, ".cgns",extension));
	printf("Öffne Datei %s im write mode\n",file_out);

	CG( cg_open(file_out,CG_MODE_WRITE,&index_file_out));
	
	printf("Writing base...\n");
	CG( cg_base_write(index_file_out,"Base",cell_dim,cell_dim,&index_base));	
	
	for(index_zone=1;index_zone<=nzones;index_zone++)
	{
	
	CG( cg_zone_read(index_file_in,index_base,index_zone,zonename,zonesize[0]));

	cgsize_t irmin[cell_dim];
	cgsize_t irmax[cell_dim];
	
	long long int i,j,k,ijk;
	float *buffer;
	long long int iMeshPoints;
	long long int jMeshPoints;
	long long int kMeshPoints;
		
	if (cell_dim==2)
	{
		irmin[0]=iStart;
		irmin[1]=jStart;	
		irmax[0]=iEnd;
		irmax[1]=jEnd;
		printf("Lade Ergebnisse aus SubZone (%zux%zu)-(%zux%zu) von Zone:%s (%zux%zu)...\n",irmin[0],irmin[1],irmax[0],irmax[1],zonename,zonesize[0][0],zonesize[0][1]);
		
		iMeshPoints=irmax[0]-irmin[0]+1;
		jMeshPoints=irmax[1]-irmin[1]+1;
		
		zonesize_SubZone[0][0]=iMeshPoints;
		zonesize_SubZone[0][1]=jMeshPoints;
		zonesize_SubZone[1][0]=iMeshPoints-1;
		zonesize_SubZone[1][1]=jMeshPoints-1;
		zonesize_SubZone[2][0]=0;
		zonesize_SubZone[2][1]=0;				
	}
	else
	{
		irmin[0]=iStart;
		irmin[1]=jStart;	
		irmin[2]=kStart;	
		irmax[0]=iEnd;
		irmax[1]=jEnd;
		irmax[2]=kEnd;
		printf("Lade Ergebnisse aus SubZone (%zux%zux%zu)-(%zux%zux%zu) von Zone:%s (%zux%zux%zu)...\n",irmin[0],irmin[1],irmin[2],irmax[0],irmax[1],irmax[2],zonename,zonesize[0][0],zonesize[0][1],zonesize[0][2]);
		
		iMeshPoints=irmax[0]-irmin[0]+1;
		jMeshPoints=irmax[1]-irmin[1]+1;
		kMeshPoints=irmax[2]-irmin[2]+1;	
		
		zonesize_SubZone[0][0]=iMeshPoints;
		zonesize_SubZone[0][1]=jMeshPoints;
		zonesize_SubZone[0][2]=kMeshPoints;		
		zonesize_SubZone[1][0]=iMeshPoints-1;
		zonesize_SubZone[1][1]=jMeshPoints-1;
		zonesize_SubZone[1][2]=kMeshPoints-1;		
		zonesize_SubZone[2][0]=0;
		zonesize_SubZone[2][1]=0;									
		zonesize_SubZone[2][2]=0;			
	}
	
	printf("Writing zone...\n");
	CG( cg_zone_write(index_file_out,index_base,zonename,*zonesize_SubZone,Structured,&index_zone));

	if(cell_dim==2)
	{
		buffer=(float*)calloc(iMeshPoints*jMeshPoints,sizeof(float));
	}
	else
	{
		buffer=(float*)calloc(iMeshPoints*jMeshPoints*kMeshPoints,sizeof(float));
	}


	//////////////////////////////////////////////////////////////////////////
	printf("Copy coordinates for zone:%d\n",index_zone);
	
	int index_coord;
	printf("Transfering CoordinateX ...\n");
	cg_coord_read(index_file_in,1,index_zone,"CoordinateX",RealSingle,irmin,irmax,buffer);
	cg_coord_write(index_file_out,1,index_zone,RealSingle,"CoordinateX",buffer,&index_coord);
	
	printf("Transfering CoordinateY ...\n");
	cg_coord_read(index_file_in,1,index_zone,"CoordinateY",RealSingle,irmin,irmax,buffer);
	cg_coord_write(index_file_out,1,index_zone,RealSingle,"CoordinateY",buffer,&index_coord);
	
	if(cell_dim==3)
	{
		printf("Transfering CoordinateZ ...\n");
		cg_coord_read(index_file_in,1,index_zone,"CoordinateZ",RealSingle,irmin,irmax,buffer);
		cg_coord_write(index_file_out,1,index_zone,RealSingle,"CoordinateZ",buffer,&index_coord);	
	}
	//////////////////////////////////////////////////////////////////////////
	
			
	

	//////////////////////////////////////////////////////////////////////////
	printf("Copy solutions for zone:%d\n",index_zone);
	char *solnames;
	solnames = malloc( 32* (Solution_End-Solution_Start+1)  );
	int index_solution;
	int index_solution_out;	
	int index_field;
	int field,nfields;
	char fieldname[ 33 ];
	DataType_t datatype;
	GridLocation_t location; //Vertex, CellCenter, IFaceCenter, JFaceCenter, and KFaceCenter. 
	printf("Kopiere Solution:%d bis %d\n",Solution_Start,Solution_End);	
	for( index_solution = Solution_Start; index_solution<=Solution_End; index_solution++ ) 
	{
	printf("Solution:%d\n",index_solution);	

	CG( cg_sol_info(index_file_in, 1, index_zone, index_solution, solname, &location));		
	CG( cg_sol_write( index_file_out,1,index_zone,solname,location,&index_solution_out) );		

	memcpy( solnames+(index_solution-Solution_Start)*32,solname,32 );
	
	CG( cg_nfields( index_file_in,1,1,index_solution,&nfields ) );
	for( field = 1; field<=nfields; field++ ) 
	{
		cg_field_info(index_file_in,1,1,1, field, &datatype, fieldname);
		printf("Variable: %s\n",fieldname);
		
		cg_field_read(index_file_in,1,index_zone,index_solution,fieldname,RealSingle,irmin,irmax,buffer);
		cg_field_write(index_file_out,1,index_zone,index_solution_out,RealSingle,fieldname,buffer,&index_field);

	}
	}
	free(buffer);
	//////////////////////////////////////////////////////////////////////////
	
	
	
	//////////////////////////////////////////////////////////////////////////
	printf("Copy zoneiter for zone:%d\n",index_zone);
	CG( cg_ziter_write( index_file_out,1,1,"ZoneIterativeData"));
	CG( cg_goto(index_file_out,1,"Zone_t",1,"ZoneIterativeData_t",1,"end"));
	CG( cg_array_write( "FlowSolutionPointers",Character,2,(cgsize_t[ ]){ 32,(Solution_End-Solution_Start+1) },solnames));	
	//////////////////////////////////////////////////////////////////////////	
	
	////////// Copies within zone
	}
	
	
	
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

}
