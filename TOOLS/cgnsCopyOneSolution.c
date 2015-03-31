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

	int Solution=0;
	int flag_solution=0;

	opterr = 0;
	int flag_reduce,flag_reduce_i,flag_reduce_j,flag_reduce_k;
	int outflag=0;
	flag_reduce=0;
	int iMeshPoints_reduced;
	int jMeshPoints_reduced;
	int kMeshPoints_reduced;
	while ((c = getopt (argc, argv, "f:tvs:e:ri:j:k:o:")) != -1)
	{
		switch (c)
		{
			case 'f':
				fileflag+=1;
				strcpy(file_in,optarg);
				break;
			case 'o':
				outflag=1;
				strcpy(file_out,optarg);
				break;
			case 's':
				flag_solution=1;
				Solution=atoi(optarg);
				break;				
			case 'r':
				flag_reduce=1;
				break;
			case 'i':
				flag_reduce_i=1;
				iMeshPoints_reduced=atoi(optarg);
				break;
			case 'j':
				flag_reduce_j=1;
				jMeshPoints_reduced=atoi(optarg);
				break;
			case 'k':
				flag_reduce_k=1;
				kMeshPoints_reduced=atoi(optarg);
				break;								
			case '?':
				printf("Geben Sie bei Bedarf eine Option an: (s)tart at solution 'x' - (e)nd at solution 'x' !\n");
				if (optopt == 'o')
					fprintf (stderr, "Option -%o requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%o'.\n", optopt);
				else
					fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
			return 1;
		}
	}
	
	printf("Starting cgnsCopyOneSolution....\n");
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
	
	if(flag_reduce==1)
	{
		if((flag_reduce_k==1)&&(flag_reduce_k==1)&&(flag_reduce_k==1))
		{
			printf("Reduziere Größe der Kopie.\n");
		}
		else
		{
			printf("ERROR: Bei Verwendung von reduzierter Größe sind die Angaben der Gitterpunkte in i,j,k Richtung verlangt (-i INT -j INT -k INT).\n");
			return 1;		
		}

	}


	printf("Öffne Datei %s im read mode\n",file_in);
	int index_file_in;
	int index_file_out;
	int index_base=1;
	CG( cg_open(file_in,CG_MODE_READ,&index_file_in));
	int cell_dim;
	CG( cg_cell_dim(index_file_in,index_base,&cell_dim));
	char solname[32];
	cgsize_t zonesize[cell_dim][cell_dim];	
	
	int number_solutions=1;
	char zonename[33];
	int nzones;
	CG( cg_nzones(index_file_in,index_base, &nzones));
	int index_zone;
	
	int NumberSolutions;
	CG( cg_nsols(index_file_in,1,1,&NumberSolutions ) );
	if(flag_solution==0){Solution=NumberSolutions;}
	float time[NumberSolutions];
	printf("Reading time-array 1...\n");
	CG( cg_goto(index_file_in,1,"BaseIterativeData_t",1,"end"));
	CG( cg_array_read_as(1,RealSingle,&time));	
	
	if(outflag==0)
	{
		sprintf(extension,"_Sol%d.cgns",Solution);
		strcpy(file_out,replace_str(file_in, ".cgns",extension));
	}
	
	if(flag_reduce==1)
	{
		CG( cg_zone_read(index_file_in,index_base,1,zonename,zonesize[0]));

		char gridChar[50];
		sprintf(gridChar,"%zux%zux%zu",zonesize[0][0],zonesize[0][1],zonesize[0][2]);
		char gridChar_reduced[50];
		sprintf(gridChar_reduced,"%zux%zux%zu-%dx%dx%d",zonesize[0][0],zonesize[0][1],zonesize[0][2],iMeshPoints_reduced,jMeshPoints_reduced,kMeshPoints_reduced);

		strcpy(file_out,replace_str(file_out, gridChar,gridChar_reduced));

		zonesize[0][0]=iMeshPoints_reduced;
		zonesize[0][1]=jMeshPoints_reduced;
		zonesize[0][2]=kMeshPoints_reduced;
		zonesize[1][0]=iMeshPoints_reduced-1;
		zonesize[1][1]=jMeshPoints_reduced-1;
		zonesize[1][2]=kMeshPoints_reduced-1;			
		
	}
	
	printf("Öffne Datei %s im write mode\n",file_out);

	CG( cg_open(file_out,CG_MODE_WRITE,&index_file_out));
	
	printf("Writing base...\n");
	CG( cg_base_write(index_file_out,"Base",cell_dim,cell_dim,&index_base));	
	
	for(index_zone=1;index_zone<=nzones;index_zone++)
	{
	
	if(flag_reduce!=1)
	{
	CG( cg_zone_read(index_file_in,index_base,index_zone,zonename,zonesize[0]));
	}

	cgsize_t irmin[cell_dim];
	cgsize_t irmax[cell_dim];
	
	long long int iMeshPoints;
	long long int jMeshPoints;
	long long int kMeshPoints;		
	long long int i,j,k,ijk;
	float *buffer;
		
	if (cell_dim==2)
	{
		irmin[0]=1;
		irmin[1]=1;	
		irmax[0]=zonesize[0][0];
		irmax[1]=zonesize[0][1];
		printf("Lade Ergebnisse aus Zone %s (%zux%zu)...\n",zonename,irmax[0],irmax[1]);
		iMeshPoints=irmax[0];
		jMeshPoints=irmax[1];
		kMeshPoints=0;
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
	int index_solution;
	int index_field;
	GridLocation_t location; //Vertex, CellCenter, IFaceCenter, JFaceCenter, and KFaceCenter. 
	CG( cg_sol_info(index_file_in, 1, index_zone, 1, solname, &location));		
	CG( cg_sol_write( index_file_out,1,index_zone,solname,location,&index_solution) );		

	int field,nfields;
	char fieldname[ 33 ];
	DataType_t datatype;
	CG( cg_nfields( index_file_in,1,1,index_solution,&nfields ) );
	for( field = 1; field<=nfields; field++ ) 
	{
		cg_field_info(index_file_in,1,1,1, field, &datatype, fieldname);
		printf("Variable: %s\n",fieldname);
		
		index_solution=Solution;

		cg_field_read(index_file_in,1,index_zone,index_solution,fieldname,RealSingle,irmin,irmax,buffer);
		cg_field_write(index_file_out,1,index_zone,1,RealSingle,fieldname,buffer,&index_field);

	}
	free(buffer);
	//////////////////////////////////////////////////////////////////////////
	

	if(flag_reduce==1)
	{
	printf("Because flag_reduced==!, Skip: Copy interfaces for zone:%d\n",index_zone);
	}
	else
	{
	
	//////////////////////////////////////////////////////////////////////////
	printf("Copy interfaces for zone:%d\n",index_zone);
	int n1to1;
	int interface;
	char connectname[33];
	char donorname[33];
	cgsize_t range[2][cell_dim];
	cgsize_t donor_range[2][cell_dim];
	int transform[cell_dim];	
	int index_interface;
	float RotationCenter[cell_dim]; 
	float RotationAngle[cell_dim];
	float Translation[cell_dim];
	CG( cg_n1to1(index_file_in,1,index_zone,&n1to1));
	for(index_interface=1;index_interface<=n1to1;index_interface++) 
	{
		CG( cg_1to1_read(index_file_in,1,index_zone,index_interface, connectname,donorname,range[0],donor_range[0],transform));
		printf("connectname: %s\n",connectname);
		CG( cg_1to1_write(index_file_out,1,index_zone, connectname, donorname, range[0], donor_range[0],transform, &interface));

		bool periodic = !cg_1to1_periodic_read(index_file_in,1,index_zone,index_interface, RotationCenter, RotationAngle, Translation);
		if(periodic)
		{
			CG( cg_1to1_periodic_write(index_file_out,1,index_zone,index_interface, RotationCenter, RotationAngle, Translation));
			printf("Periodic BC detected and copied.\n");
		}
		
	}
	
	//////////////////////////////////////////////////////////////////////////

	}
	
	
	//////////////////////////////////////////////////////////////////////////
	printf("Copy boundary conditions for zone:%d\n",index_zone);
	int nbc,bc;
	CG( cg_nbocos( index_file_in,1,index_zone,&nbc ) );
	int BC;
	for( bc = 1; bc<=nbc; bc++ ) 
	{
		PointSetType_t bcPointSetType;
		int NormalIndexVektor[ cell_dim ];
		cgsize_t NormalListSize;
		int NumberDatasets;
		char bcName[ 33 ];
		char bcFamName[ 33 ];
		cgsize_t bcNumberPoints;
		int bcndata,bcnorm[ cell_dim ];
		DataType_t bcDataType;
		BCType_t bcType;
		cgsize_t bcPoints[ 2 ][ cell_dim ];
		void *NormalList;

		CG( cg_boco_info( index_file_in,1,index_zone,bc,bcName,&bcType,&bcPointSetType,&bcNumberPoints,NormalIndexVektor,&NormalListSize,&bcDataType,&NumberDatasets ) );
		CG( cg_boco_read(index_file_in,1,index_zone,bc,bcPoints [0], NormalList));
		CG( cg_goto( index_file_in,1,zonename,0,"ZoneBC",0,bcName,0,"end" ) );
		CG( cg_famname_read( bcFamName ) );
					
		printf("bc: %s - family-name: %s\n",bcName,bcFamName);
		CG( cg_boco_write( index_file_out,1,index_zone,bcName,bcType,bcPointSetType,bcNumberPoints,bcPoints [0],&BC ) );
		CG( cg_goto( index_file_out,1,zonename,0,"ZoneBC",0,bcName,0,"end" ) );
		CG( cg_famname_write( bcFamName ) );		
	}
	//////////////////////////////////////////////////////////////////////////




	//////////////////////////////////////////////////////////////////////////
	printf("Copy zoneiter for zone:%d\n",index_zone);
	cgsize_t idata[2];
	idata[0]=32;
	idata[1]=1;
	CG( cg_ziter_write(index_file_out,1,index_zone,"ZoneIterativeData"));
	CG( cg_goto(index_file_out,1,"Zone_t",index_zone,"ZoneIterativeData_t",1,"end"));
	CG( cg_array_write("FlowSolutionPointers",Character,2,idata,solname));			
	//////////////////////////////////////////////////////////////////////////	
	
	////////// Copies within zone
	}
	
	
	
	//////////////////////////////////////////////////////////////////////////
	printf("Writing new baseiter\n");
	CG( cg_biter_write(index_file_out,1,"BaseIterativeData",1));
	CG( cg_goto(index_file_out,1,"BaseIterativeData_t",1,"end"));
	float time_out=time[Solution-1];
	cgsize_t nuse=1;
	CG( cg_array_write("TimeValues",RealSingle,1,&nuse,&time_out));
	//////////////////////////////////////////////////////////////////////////
	
	
	
	//////////////////////////////////////////////////////////////////////////
	printf("Copy descriptors\n");
	CG( cg_goto( index_file_in,1,"end" ) );
	int desc,ndesc;
	CG( cg_ndescriptors( &ndesc ) );
	for( desc = 1; desc<=ndesc; desc++ )
	{
		char* dtext,dname[ 33 ];
		CG( cg_goto( index_file_in,1,"end" ) );
		CG( cg_descriptor_read( desc,dname,&dtext ) );
		CG( cg_goto( index_file_out,1,"end" ) );
		CG( cg_descriptor_write( dname,dtext ) );
	}
	
	
	//////////////////////////////////////////////////////////////////////////
	printf("Copy Families for boundary conditions\n");
	int nfamilies;
	int index_family;
	char FamilyName[33];
	char FamBCName[33];
	int nFamBC, nGeo;
	int Fam;
	int BC=1;
	BCType_t BCType_Fam;
	CG( cg_nfamilies(index_file_in,1, &nfamilies));
	for(index_family=1;index_family<=nfamilies;index_family++)
	{
		CG( cg_family_read(index_file_in,1,index_family,FamilyName,&nFamBC,&nGeo));
		if (strcmp(FamilyName,"Unspecified")!=0)
		{
		CG( cg_family_write(index_file_out,1,FamilyName, &Fam)); 
		printf("Family: %s, index. %d\n",FamilyName,index_family);
		CG( cg_fambc_read(index_file_in,1,index_family,1,FamBCName, &BCType_Fam));
		CG( cg_fambc_write(index_file_out,1,Fam, FamBCName,BCType_Fam, &BC));
		}


	}
	//////////////////////////////////////////////////////////////////////////	
	
	CG (cg_close(index_file_out));
	CG (cg_close(index_file_in));		

}
