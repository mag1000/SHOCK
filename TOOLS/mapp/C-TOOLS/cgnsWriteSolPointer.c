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
	
	int start_value=0;

	while ((c = getopt (argc, argv, "i:s:")) != -1)
	{
		switch (c)
		{
			case 'i':
				inputflag=1;
				ivalue = optarg;
				break;
			case 's':
				start_value = atoi(optarg);
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
	printf("Der Index der ersten Solution lautet: %d (Startwert durch Option -s (INT))\n",start_value);
	int index_file,nsols;

	printf("Opening source: %s\n",file);
	CG( cg_open(file,CG_MODE_MODIFY,&index_file));
	CG( cg_nsols(index_file,1,1,&nsols ) );
	
	printf("Es gibt %d solutions. Pointers werden geschrieben...\n",nsols);



	int sol;

	char *solnames;
	solnames = malloc( 32* nsols  );
	for( sol = 0; sol<nsols; sol++ ) {

	char solname[33];
	if( snprintf( solname,33,"Sol_%04d",sol+start_value )>33 ) {
	printf( "(EE) Could not write solution #%i because consecutive numbering exceeded string length\n",sol );
	exit( 1 );
	}

	memcpy( solnames+sol*32,solname,32 );
	}
	CG( cg_ziter_write( index_file,1,1,"ZoneIterativeData"));
	CG( cg_goto(index_file,1,"Zone_t",1,"ZoneIterativeData_t",1,"end"));
	CG( cg_array_write( "FlowSolutionPointers",Character,2,(cgsize_t[ ]){ 32,nsols },solnames));
}

