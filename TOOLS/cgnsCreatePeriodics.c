#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"
#include "unistd.h"
#define CG(cmd) if(cmd)cg_error_exit( );
int main(int argc, char *argv[])
{
	char *ivalue = NULL;
	int inputflag;
	float translation_value;
	char file[500];
	int c;
	int start_solution,startflag;
	int end_solution,endflag;
	opterr = 0;


	while ((c = getopt (argc, argv, "i:t:")) != -1)
	{
		switch (c)
		{
			case 'i':
				inputflag=1;
				ivalue = optarg;
				break;
			case 't':
				translation_value = atof(optarg);
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

	printf("Erzeuge periodische Randbedingungen für %s mit %f Verschiebung\n",file,translation_value);

	char donorname[33];
	char interfacename[33];
	int index_file,index_base,index_zone,index_interface;
	cgsize_t* rangeOfInterface;
	cgsize_t* donorRangeOfInterface;
	cgsize_t zonesize[3][3];
      	int transform[]={1,2,3};
	float rotationCenter[]={0,0,0};
	float rotationAngle[]={0,0,0};
	float translation[]={0,0,0};

		
	index_base=1;
	index_zone=1;
	
	CG( cg_open(file,CG_MODE_MODIFY,&index_file));

	CG( cg_zone_read(index_file,index_base,index_zone,donorname, zonesize[0]));
	
	rangeOfInterface= (cgsize_t *)calloc(2*3, sizeof(cgsize_t));
	donorRangeOfInterface= (cgsize_t *)calloc(2*3, sizeof(cgsize_t));
	
//WRITING INTERFACE1
	sprintf(interfacename,"PERIODIC1");
	/* lower range index */
	rangeOfInterface[0]=1;
	rangeOfInterface[1]=1;
	rangeOfInterface[2]=1;		
	rangeOfInterface[3+0]=zonesize[0][0];
	rangeOfInterface[3+1]=zonesize[0][1];
	rangeOfInterface[3+2]=1;	

	donorRangeOfInterface[0]=1;
	donorRangeOfInterface[1]=1;
	donorRangeOfInterface[2]=zonesize[0][2];
	donorRangeOfInterface[3+0]=zonesize[0][0];
	donorRangeOfInterface[3+1]=zonesize[0][1];
	donorRangeOfInterface[3+2]=zonesize[0][2];	
			
	translation[2]=-translation_value;
	//		CREATE INTERFACE1
	CG( cg_1to1_write(index_file,index_base,index_zone,interfacename,donorname,rangeOfInterface,donorRangeOfInterface,transform,&index_interface)) ;
	//		CREATE PERIODICITY1
	CG( cg_1to1_periodic_write(index_file,index_base,index_zone,index_interface,rotationCenter,rotationAngle,translation));	



	
//WRITING INTERFACE2
	sprintf(interfacename,"PERIODIC2");
	/* lower range index */
	donorRangeOfInterface[0]=1;
	donorRangeOfInterface[1]=1;
	donorRangeOfInterface[2]=1;		
	donorRangeOfInterface[3+0]=zonesize[0][0];
	donorRangeOfInterface[3+1]=zonesize[0][1];
	donorRangeOfInterface[3+2]=1;	

	rangeOfInterface[0]=1;
	rangeOfInterface[1]=1;
	rangeOfInterface[2]=zonesize[0][2];
	rangeOfInterface[3+0]=zonesize[0][0];
	rangeOfInterface[3+1]=zonesize[0][1];
	rangeOfInterface[3+2]=zonesize[0][2];	
			
	translation[2]=translation_value;
	//		CREATE INTERFACE1
	CG( cg_1to1_write(index_file,index_base,index_zone,interfacename,donorname,rangeOfInterface,donorRangeOfInterface,transform,&index_interface)) ;
	//		CREATE PERIODICITY1
	CG( cg_1to1_periodic_write(index_file,index_base,index_zone,index_interface,rotationCenter,rotationAngle,translation));		
	
	CG( cg_close(index_file));
}
