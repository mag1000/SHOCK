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
	
	int index_file,nsols;

	printf("Opening source: %s\n",file);
	CG( cg_open(file,CG_MODE_MODIFY,&index_file));
	CG( cg_nsols(index_file,1,1,&nsols ) );
	
	printf("Es gibt %d solutions\n",nsols);
	
	CG( cg_goto(index_file,1,"Zone_t",1,"ZoneIterativeData_t",1,"end"));
	char *solutionarray = malloc(nsols*32*sizeof(char));
	int i;
	char solname[32];
	for(i=0;i<nsols;i++)
	{
	sprintf(solname,"FlowSolution%d",i);
	memcpy(solutionarray+i*32,solname,32);
	}
	CG( cg_array_write( "FlowSolutionPointers",Character,2,(cgsize_t[ ]){ 32,nsols},solutionarray ) );
	
	CG( cg_array_read( 1,solutionarray ) );
	printf("flowsolutionpointer:\n %s\n",solutionarray);
}










