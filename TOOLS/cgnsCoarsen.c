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


void coarse(long long int coarsefactors[3],cgsize_t zonesize_coarse[3][3],float *buffer_coarse,float *buffer_fine)
{
	long long int i,j,k;
	long long int i_fine,j_fine,k_fine;
	long long int ijk_fine,ijk_coarse;
	long long int zonesize_fine_i=zonesize_coarse[0][0]*coarsefactors[0];
	long long int zonesize_fine_j=zonesize_coarse[0][1]*coarsefactors[1];
	long long int zonesize_fine_k=zonesize_coarse[0][2]*coarsefactors[2];

	for(i=0;i<zonesize_coarse[0][0];i++)
	{
		i_fine=coarsefactors[0]*i;
		for(j=0;j<zonesize_coarse[0][1];j++)
		{
			j_fine=coarsefactors[1]*j;
			for(k=0;k<zonesize_coarse[0][2];k++)
			{
				k_fine=coarsefactors[2]*k;
				ijk_fine=k_fine*zonesize_fine_j*zonesize_fine_i+j_fine*zonesize_fine_i+i_fine;
				ijk_coarse=k*zonesize_coarse[0][1]*zonesize_coarse[0][0]+j*zonesize_coarse[0][0]+i;				
				buffer_coarse[ijk_coarse]=buffer_fine[ijk_fine];
			}
		}
	}
}


int main(int argc, char *argv[])
{
	char file_fine[500];
	char file_coarse[500];	
	
	char *ivalue = NULL;
	int inputflag;
	int ijkflag=0;
	int fileflag=0;		
	int c;
	
	long long int coarsefactors[3];


	opterr = 0;

	while ((c = getopt (argc, argv, "f:i:j:k:")) != -1)
	{
		switch (c)
		{
			case 'f':
				fileflag+=1;
				strcpy(file_fine,optarg);
				break;						
			case 'i':
				ijkflag+=1;
				coarsefactors[0] = atoi(optarg);
				break;				
			case 'j':
				ijkflag+=1;
				coarsefactors[1] = atoi(optarg);
				break;	
			case 'k':
				ijkflag+=1;
				coarsefactors[2] = atoi(optarg);
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
	printf("Starting cgnsCoarsen....\n");
	int errorflag=0;
	if(fileflag!=1)
	{
		printf("ERROR: Keine Datei angegeben (-f 'filename.cgns')!\n");
		errorflag=1;		
	}

	if(ijkflag!=3)
	{
		printf("ERROR: Geben Sie die Vergroeberungsfaktoren in i-/j-/k-Richtung an (z.B. -i 2 -j 2 -k 1 )!\n");
		errorflag=1;
	}
	
	
	if(errorflag==1)
	{
		return 1;
	}

	printf("Das Gitter %s wird mit folgenden Faktoren vergroebert: %lld %lld %lld .\n",
	file_fine,coarsefactors[0],coarsefactors[1],coarsefactors[2]);
	

//	Öffnen beider Dateien
	int index_file_fine,index_file_coarse;
	int index_base=1;
	int index_zone=1;
	int index_sol=1;
	char zonename[33];
	cgsize_t zonesize_fine[3][3];
	cgsize_t zonesize_coarse[3][3];	

	printf("Opening %s in read mode...\n",file_fine);	
	CG( cg_open(file_fine,CG_MODE_READ,&index_file_fine));
	printf("Reading zone...\n");	
	CG( cg_zone_read(index_file_fine,index_base,index_zone,zonename,zonesize_fine[0]));

	long long int iMeshPoints=zonesize_fine[0][0];
	long long int jMeshPoints=zonesize_fine[0][1];
	long long int kMeshPoints=zonesize_fine[0][2];
	
	long long int iMeshPointsCoarse=(long long int)iMeshPoints/coarsefactors[0];
	long long int jMeshPointsCoarse=(long long int)jMeshPoints/coarsefactors[1];
	long long int kMeshPointsCoarse=(long long int)kMeshPoints/coarsefactors[2];

	char oldname[500];
	char newname[500];
	sprintf(oldname,"%lldx%lldx%lld",iMeshPoints,jMeshPoints,kMeshPoints);
	sprintf(newname,"%lldx%lldx%lld-%lldx%lldx%lld",iMeshPoints,jMeshPoints,kMeshPoints,iMeshPointsCoarse,jMeshPointsCoarse,kMeshPointsCoarse);
	strcpy(file_coarse,replace_str(file_fine, oldname,newname));
	strcpy(file_coarse,replace_str(file_coarse,".cgns","_coarsen.cgns"));
	printf("Opening %s in write mode...\n",file_coarse);	

	CG( cg_open(file_coarse,CG_MODE_WRITE,&index_file_coarse));
	
	printf("Writing base...\n");
	CG( cg_base_write(index_file_coarse,"Base",3,3,&index_base));
	
	printf("Writing zone...\n");
	zonesize_coarse[0][0]=iMeshPointsCoarse;	zonesize_coarse[0][1]=jMeshPointsCoarse;	zonesize_coarse[0][2]=kMeshPointsCoarse;
	zonesize_coarse[1][0]=iMeshPointsCoarse-1;	zonesize_coarse[1][1]=jMeshPointsCoarse-1;	zonesize_coarse[1][2]=kMeshPointsCoarse-1;
	zonesize_coarse[2][0]=0;			zonesize_coarse[2][1]=0;			zonesize_coarse[2][2]=0;		
	CG( cg_zone_write(index_file_coarse,index_base,zonename,*zonesize_coarse,Structured,&index_zone));
	

	cgsize_t irmin_fine[3];
	cgsize_t irmax_fine[3];	
	cgsize_t irmin_coarse[3];
	cgsize_t irmax_coarse[3];
	
	irmin_coarse[0]=1;
	irmin_coarse[1]=1;
	irmin_coarse[2]=1;
	irmax_coarse[0]=zonesize_coarse[0][0];
	irmax_coarse[1]=zonesize_coarse[0][1];
	irmax_coarse[2]=zonesize_coarse[0][2];		

	irmin_fine[0]=1;
	irmin_fine[1]=1;
	irmin_fine[2]=1;
	irmax_fine[0]=zonesize_fine[0][0];
	irmax_fine[1]=zonesize_fine[0][1];
	irmax_fine[2]=zonesize_fine[0][2];	
	
	float *buffer_fine;
	float *buffer_coarse;

	buffer_coarse=(float*)calloc(irmax_coarse[0]*irmax_coarse[1]*irmax_coarse[2],sizeof(float));
	buffer_fine=(float*)calloc(irmax_fine[0]*irmax_fine[1]*irmax_fine[2],sizeof(float));

	int index_coord;
	printf("Transfering CoordinateX ...\n");
	cg_coord_read(index_file_fine,1,1,"CoordinateX",RealSingle,irmin_fine,irmax_fine,buffer_fine);
	coarse(coarsefactors,zonesize_coarse,buffer_coarse,buffer_fine);
	cg_coord_write(index_file_coarse,1,1,RealSingle,"CoordinateX",buffer_coarse,&index_coord);
	
	printf("Transfering CoordinateY ...\n");
	cg_coord_read(index_file_fine,1,1,"CoordinateY",RealSingle,irmin_fine,irmax_fine,buffer_fine);
	coarse(coarsefactors,zonesize_coarse,buffer_coarse,buffer_fine);
	cg_coord_write(index_file_coarse,1,1,RealSingle,"CoordinateY",buffer_coarse,&index_coord);
	
	printf("Transfering CoordinateZ ...\n");
	cg_coord_read(index_file_fine,1,1,"CoordinateZ",RealSingle,irmin_fine,irmax_fine,buffer_fine);
	coarse(coarsefactors,zonesize_coarse,buffer_coarse,buffer_fine);
	cg_coord_write(index_file_coarse,1,1,RealSingle,"CoordinateZ",buffer_coarse,&index_coord);	

	
	int NumberSolutions;
	CG( cg_nsols(index_file_fine,1,1,&NumberSolutions ) );

	int index_solution;
	int index_field;
	//char solname[NumberSolutions*32+1];
	//char solnameArray[NumberSolutions][32];
	
	//strcpy(solname,"");
	long long int i;
	
	float time[NumberSolutions];
	printf("Reading time-array...\n");
	CG( cg_goto(index_file_fine,1,"BaseIterativeData_t",1,"end"));
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
		CG( cg_sol_info(index_file_fine, 1, 1, index_solution, solname, &location));
		memcpy( solnames+(index_solution-1)*32,solname,32 );
		
		printf("Transforming Solution:\"%s\" (%d of %d) at time: %f\n",solname,index_solution,NumberSolutions,time[index_solution]);
		CG( cg_sol_write( index_file_coarse,1,1,solname,location,&index_sol) );
		
		CG( cg_nfields( index_file_fine,1,1,index_solution,&nfields ) );
		for( field = 1; field<=nfields; field++ ) 
		{
			cg_field_info(index_file_fine,1,1,index_solution, field, &datatype, fieldname);
			printf("%s...\n",fieldname);
			cg_field_read(index_file_fine,1,1,index_solution,fieldname,RealSingle,irmin_fine,irmax_fine,buffer_fine);
			coarse(coarsefactors,zonesize_coarse,buffer_coarse,buffer_fine);
			cg_field_write(index_file_coarse,1,1,index_sol,RealSingle,fieldname,buffer_coarse,&index_field);
		}
	}
	
	
	printf("Check Solnames: %s\n",solnames);
	printf("Writing new baseiter,zoneiter...\n");
	CG( cg_biter_write(index_file_coarse,1,"TimeIterValues",NumberSolutions));
	CG( cg_goto(index_file_coarse,1,"BaseIterativeData_t",1,"end"));
	cgsize_t nuse=NumberSolutions;
	CG( cg_array_write("TimeValues",RealSingle,1,&nuse,&time));
	cgsize_t idata[2];
	idata[0]=32;
	idata[1]=NumberSolutions;
	CG( cg_ziter_write(index_file_coarse,1,1,"ZoneIterativeData"));
	CG( cg_goto(index_file_coarse,1,"Zone_t",1,"ZoneIterativeData_t",1,"end"));
	CG( cg_array_write("FlowSolutionPointers",Character,2,idata,solnames));	
	CG (cg_close(index_file_coarse));
	CG (cg_close(index_file_fine));
}

