#include "Load.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pcgnslib.h>

#define PHYSDIM 3

#define ORIGREF "SDDC_OREF"
#define CG(cmd) if(cmd)cg_error_exit( );
#define DIM MESHDIMENSIONS
#define DEBUG 0
#define TM 1

#if DIM==2
#	define DIM3(...)
#else
#	define DIM3(...) __VA_ARGS__
#endif

float ts;

#if TM==1
#	define TM_START() ts=(float)MPI_Wtime( );
#else
#	define TM_START()
#endif

#if TM==1
#	define TM_END(msg) if(pnt_config->MPI_rank==0){printf(msg " took %g seconds\n",(float)MPI_Wtime( )-ts );}
#else
#	define TM_END(msg)
#endif

static void assignIF( struct strct_configuration* pnt_config,int i,bool j,int index) {
	int *const sidemap[ ][ 2 ]= {
		{ &pnt_config->InterfaceNeighbourLeft,&pnt_config->InterfaceNeighbourRight },
		{ &pnt_config->InterfaceNeighbourBottom,&pnt_config->InterfaceNeighbourTop },
		{ &pnt_config->InterfaceNeighbourBehind,&pnt_config->InterfaceNeighbourInFront }
	};

	*sidemap[ i ][ j ]= index;
}

static bool assignBC( struct strct_configuration* pnt_config,int i,bool j,BCType_t bct ) {
	char farfieldstring[ ]= "BCFarfield";
	char inflowstring[ ]= "BCInflow";
	char outflowstring[ ]= "BCOutflow";
	char outflowsubstring[ ]= "BCOutflowSubsonic";
	char wallivstring[ ]= "BCWallInviscid";
	char wallvstring[ ]= "BCWallViscous";
	char inflowsupstring[ ]= "BCInflowSupersonic";
	char inflowsubstring[ ]= "BCInflowSubsonic";
	char wallvsisotring[ ]= "BCWallViscousIsothermal";

	char (*const sidemap[ ][ 2 ])[ 30 ]= {
		{ &pnt_config->BC_Left,&pnt_config->BC_Right },
		{ &pnt_config->BC_Bottom,&pnt_config->BC_Top }
		DIM3( ,{ &pnt_config->BC_Behind,&pnt_config->BC_InFront } )
	};

	switch( bct ) {
	case BCFarfield:
		strcpy( *sidemap[ i ][ j ],farfieldstring );
		break;
	case BCInflow:
		strcpy( *sidemap[ i ][ j ],inflowstring );
		break;
	case BCOutflow:
		strcpy( *sidemap[ i ][ j ],outflowstring );
		break;
	case BCOutflowSubsonic:
		strcpy( *sidemap[ i ][ j ],outflowsubstring );
		break;
	case BCWallInviscid:
		strcpy( *sidemap[ i ][ j ],wallivstring );
		break;
	case BCWallViscous:
		strcpy( *sidemap[ i ][ j ],wallvstring );
		break;
	case BCInflowSupersonic:
		strcpy( *sidemap[ i ][ j ],inflowsupstring );
		break;
	case BCInflowSubsonic:
		strcpy( *sidemap[ i ][ j ],inflowsubstring );
		break;
	case BCWallViscousIsothermal:
		strcpy( *sidemap[ i ][ j ],wallvsisotring );
		break;
	default:
		printf( "(EE) Unsupported BC with numerical representation %i\n",bct );
		return false;
	}

	return true;
}

/** Helper struct for interface data
 *
 * Holds all information which needs to be passed for connectivity
 * information. Usually, this would only be rank and transform. For
 * simplicity, the members are statically allocated.
 */
typedef struct {
	unsigned neighbour; ///< Zoneid of the neighbour at that interface
	cgsize_t drange[ 2 ][ DIM ]; ///< donor range in the neighbour
	cgsize_t range[ 2 ][ DIM ]; ///< range of interface
	int transform[ DIM ]; ///< Transformation values
	float rc[ PHYSDIM ],ra[ PHYSDIM ],tr[ PHYSDIM ]; ///< Auxiliary data
} interface;

void* normalizeBytes( struct strct_configuration* pnt_config, void* mem,size_t size,int count ) {

	if (pnt_config->flag_swapDivisionFile==1) //force to stop swapping
	{
		return mem;
	}
	if (pnt_config->flag_swapDivisionFile==2) //force to swap
	{
		char swap;



		size_t i;
		int k;
		for( k = 0; k<count; k++ ) {
			char* buf = ( (char*)mem )+k*size;
			for( i = 0; i<(int)( size/2 ); i++ ) {
				swap = buf[ i ];
				buf[ i ]= buf[ size-1-i ];
				buf[ size-1-i ]= swap;
			}
		}

		return mem;	
	}
		
	//automatic ("scans" architecture)
	
	const int one = 1;
	if( *( (char*)(&one) )==1 )
	{
		return mem;
	}
	else
	{

		char swap;
		size_t i;
		int k;
		for( k = 0; k<count; k++ ) {
			char* buf = ( (char*)mem )+k*size;
			for( i = 0; i<(int)( size/2 ); i++ ) {
				swap = buf[ i ];
				buf[ i ]= buf[ size-1-i ];
				buf[ size-1-i ]= swap;
			}
		}

		return mem;
	}
}

/** Determine domain and offset in the grid
 *
 * Based upon the rank and the division file, determines the name of the
 * zone and the region within it which the given rank is supposed to
 * calculate for.
 *
 * @param Filename or path to the division file (currently CGNS)
 *
 * @param Rank for which the region is to be sought
 *
 * @param[out,caller-allocates] Name of the zone which the rank is to
 * calculate
 *
 * @param[out,caller-allocates] Lower and higher corner of the region
 * inside of the zone which the rank is to calculate. That is
 * 2*Dimension elements.
 *
 * @param[out,caller-allocates] Neighbouring ranks at the lower-x,
 * higher-x, lower-y, higher-y, lower-z, higher-z border (last two only
 * for 3 dimensions).
 *
 * @param[out] Number of interfaces which have a rank other than 0.
 *
 * @return The CG Error code corresponding to the last error or 0 if no
 * error occured.
 */
static void loadDivision( struct strct_configuration* pnt_config,char* const filename,unsigned divzone,char zonename[ 33 ],cgsize_t corners[ 2 ][ DIM ],interface neighbours[ DIM ][ 2 ],int* ni,const int proccount ) {

	MPI_File divfile;
	MPI_Status status;
	MPI_Offset filesize;

	size_t ifblocklen = sizeof( int )+DIM*sizeof( int )+2*DIM*sizeof( cgsize_t )*2+PHYSDIM*sizeof( float )*3;
	size_t datalen = 32+sizeof( cgsize_t )*DIM*2+ifblocklen*2*DIM;

	MPI_File_open( MPI_COMM_WORLD,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,&divfile );
	MPI_File_get_size( divfile,&filesize ); /* in bytes */

	if( filesize!=( proccount*datalen ) ) {
		printf( "(EE) Division file '%s' has size %i but should be processors (%i) times data (%i)\n",filename,(int)filesize,proccount,(int)datalen );
		exit( 1 );
	}

	char* dataset = malloc( datalen );
	MPI_File_set_view( divfile,( divzone-1 )*datalen,MPI_CHAR,MPI_CHAR,"native",MPI_INFO_NULL );
	MPI_File_read( divfile,dataset,datalen,MPI_CHAR,&status );

	MPI_File_close( &divfile );

/** Untangle 'em */
	memcpy( zonename,dataset,32 );
	zonename[ 32 ]= '\0';

	cgsize_t(* incorners)[ DIM ]= (cgsize_t(*)[ DIM ])( dataset+32 );
	normalizeBytes( pnt_config,memcpy( corners[ 0 ],incorners,2*DIM*sizeof( cgsize_t ) ),sizeof( cgsize_t ),2*DIM );

	int dir;
	*ni = 0;
	for( dir = 0; dir<DIM; dir++ ) {
		int higher;
		for( higher = 0; higher<2; higher++ ) {
			int index = dir*2+higher;

			int* neighbour= (int*)( dataset+32+DIM*2*sizeof( cgsize_t )+index*ifblocklen );
			int* transform= neighbour+1;
			cgsize_t(* range)[ DIM ]= (cgsize_t(*)[ DIM ])( transform+DIM ),
				(* drange)[ DIM ]= range+2;
			float* rc = (float*)( drange+2 ),
				* ra = rc+PHYSDIM,
				* tr = ra+PHYSDIM;

			if( *neighbour ) {
				( *ni )++;
				neighbours[ dir ][ higher ].neighbour = *(int*)normalizeBytes( pnt_config,neighbour,sizeof( int ),1 );
				normalizeBytes( pnt_config,memcpy( &neighbours[ dir ][ higher ].range,range,2*DIM*sizeof( cgsize_t ) ),sizeof( cgsize_t ),2*DIM );
				normalizeBytes( pnt_config,memcpy( &neighbours[ dir ][ higher ].drange,drange,2*DIM*sizeof( cgsize_t ) ),sizeof( cgsize_t ),2*DIM );
				normalizeBytes( pnt_config,memcpy( &neighbours[ dir ][ higher ].transform,transform,DIM*sizeof( int ) ),sizeof( int ),DIM );
				normalizeBytes( pnt_config,memcpy( &neighbours[ dir ][ higher ].rc,rc,PHYSDIM*sizeof( float ) ),sizeof( float ),PHYSDIM );
				normalizeBytes( pnt_config,memcpy( &neighbours[ dir ][ higher ].ra,ra,PHYSDIM*sizeof( float ) ),sizeof( float ),PHYSDIM );
				normalizeBytes( pnt_config,memcpy( &neighbours[ dir ][ higher ].tr,tr,PHYSDIM*sizeof( float ) ),sizeof( float ),PHYSDIM );
			}
		}
	}

	free( dataset );
}

bool loadFile( struct strct_configuration* pnt_config,vta* x,vta* y,vta* z,vta* u,vta* v,vta* w,vta* rho,vta* p ) {
/** Check that celldimension is correct */
	int file,physdim,celldim,base = 1;
	char basename[ 33 ];
	char* zonename = pnt_config->Zonename;
	interface neighbours[ DIM ][ 2 ] = { { { 0 } } };
	cgsize_t corners[ 2 ][ DIM ],zonesize[ 3 ][ DIM ];
	int ni;
/** Prepare most following steps. I.e.: Analyze division file and figure
 * out related data such as the id of the original zone (which is
 * refered to by name by the division file). */
	TM_START( )
#if (CGNS_VERSION<=3210)
	CG( cgp_pio_mode(CGP_COLLECTIVE));
#else
	CG( cgp_pio_mode(CGP_COLLECTIVE, MPI_INFO_NULL));
#endif
	CG( cgp_open( pnt_config->chr_MeshPath,CG_MODE_READ,&file ) );
	TM_END( "SHOCK: Opening meshfile" )

	TM_START ( )

	if(pnt_config->MPI_rank==0)
	{
		if(pnt_config->flag_swapDivisionFile==1)
		{
			printf("SHOCK: Forced: No Swapping!\n");
		}
		else if(pnt_config->flag_swapDivisionFile==2)
		{
			printf("SHOCK: Forced: Swapping!\n");
		}
		else
		{
		const int one = 1;
		if( *( (char*)(&one) )==1 )
		{
			printf("SHOCK: LittleEndian-System: No Swapping.\n");
		}
		if( *( (char*)(&one) )!=1 )
		{
			printf("SHOCK: BigEndian-System: DivisionFile is swapped.\n");
		}
		}
	}
	loadDivision( pnt_config,pnt_config->chr_DivisionPath,pnt_config->MPI_rank+1,zonename,corners,neighbours,&ni,pnt_config->MPI_size );
	TM_END ( "SHOCK: Loading divisionfile" )
	CG( cg_base_read( file,base,basename,&celldim,&physdim ) );

	if( celldim!=DIM ) {
		printf( "(EE) Mesh file dimension %i does not match compile time dimension %i\n",celldim,DIM );
		exit( 1 );
	}

/** Load number of iterations, if available */
	int ndesc,desc;

	CG( cg_goto( file,base,"end" ) );
	CG( cg_ndescriptors( &ndesc ) );
	for( desc = 1; desc<=ndesc; desc++ ) {
		char* dtext,dname[ 33 ];
		CG( cg_descriptor_read( desc,dname,&dtext ) );
		if( !strcmp( dname,"Iterations" ) ) {
			pnt_config->int_StartIteration = atol( dtext );

			cg_free( dtext );
			break;
		} else
			cg_free( dtext );
	}
/** Open division file, figure out which zone we're part of and where in
 * that zone our domain lies. */
	int nzones,zone;

	CG( cg_nzones( file,base,&nzones ) );

	for( zone = 1; zone<=nzones; zone++ ) {
		char tzname[ 33 ];
		CG( cg_zone_read( file,base,zone,tzname,zonesize[ 0 ] ) );
		if( !strcmp( tzname,zonename ) )
			break;
	}

	if( zone>nzones ) {
		printf( "(EE) Could not find associated zone '%s' in the mesh file\n",zonename );
		exit( 1 );
	}
	int i;
/** ATTENTION: If we would not assume that pnt_config is initially
 * "clean", we would have to free the pointed-to-pointers of a few
 * pointer-to-pointer arrays like so:
 *
 * for( i = 0; i<pnt_config->NumberInterfaces; i++ ) {
 * 	free( pnt_config->TransformMatrixOfInterface[ i ] );
 * 	free( pnt_config->DonorRangeOfInterface[ i ] );
 * 	free( pnt_config->Translation[ i ] );
 * 	free( pnt_config->RotationCenter[ i ] );
 * 	free( pnt_config->RotationAngle[ i ] );
 * }
 *
 * Also, the pnt_config entries should then be intialized to NULL on
 * startup and be realloc'd rather than malloced. For now, we assume a
 * clean start and alloc everything without checking.
 */

/** Allocation of members in pnt_config based upon the number of
 * neighbours (rather than always 2*DIM) and, of course, the
 * celldimension. */
	pnt_config->MPI_rankNeighbours = malloc( ni*sizeof( unsigned ) );

	pnt_config->TransformMatrixOfInterface = malloc( ni*sizeof( int* ) );
	pnt_config->DonorRangeOfInterface = malloc( ni*sizeof( cgsize_t* ) );
	pnt_config->RangeOfInterface = malloc( ni*sizeof( cgsize_t* ) );
	pnt_config->Translation = malloc( ni*sizeof( float* ) );
	pnt_config->RotationCenter = malloc( ni*sizeof( float* ) );
	pnt_config->RotationAngle = malloc( ni*sizeof( float* ) );

/** Load information about neighbours into pnt_config. */
	pnt_config->NumberInterfaces = ni;
	for( i = 0; i<DIM; i++ ) {
		pnt_config->zonesize[ i ]= corners[ 1 ][ i ]-corners[ 0 ][ i ]+1;
		pnt_config->offset[ i ]= corners[ 0 ][ i ]-1;
	}

	int j,ii = 0;

	for( i = 0; i<3; i++ )
		for( j = 0; j<2; j++ ) {
			assignIF( pnt_config,i,j,neighbours[ i ][ j ].neighbour?ii:NO_NEIGHBOUR );

			if( ( i<2 || DIM==3 )&& neighbours[ i ][ j ].neighbour ) {
				pnt_config->TransformMatrixOfInterface[ ii ] = malloc( DIM*sizeof( int ) );
				pnt_config->DonorRangeOfInterface[ ii ]= malloc( 2*DIM*sizeof( cgsize_t ) );
				pnt_config->RangeOfInterface[ ii ]= malloc( 2*DIM*sizeof( cgsize_t ) );
				pnt_config->Translation[ ii ]= malloc( PHYSDIM*sizeof( float ) );
				pnt_config->RotationCenter[ ii ]= malloc( PHYSDIM*sizeof( float ) );
				pnt_config->RotationAngle[ ii ]= malloc( PHYSDIM*sizeof( float ) );

				pnt_config->MPI_rankNeighbours[ ii ]= neighbours[ i ][ j ].neighbour-1;

				memcpy( pnt_config->TransformMatrixOfInterface[ ii ],neighbours[ i ][ j ].transform,DIM*sizeof( int ) );
				memcpy( pnt_config->DonorRangeOfInterface[ ii ],neighbours[ i ][ j ].drange,2*DIM*sizeof( cgsize_t ) );
				memcpy( pnt_config->RangeOfInterface[ ii ],neighbours[ i ][ j ].range,2*DIM*sizeof( cgsize_t ) );
				memcpy( pnt_config->Translation[ ii ],neighbours[ i ][ j ].tr,PHYSDIM*sizeof( float ) );
				memcpy( pnt_config->RotationCenter[ ii ],neighbours[ i ][ j ].rc,PHYSDIM*sizeof( float ) );
				memcpy( pnt_config->RotationAngle[ ii ],neighbours[ i ][ j ].ra,PHYSDIM*sizeof( float ) );

				ii++;
			}
		}
/**
 * Load boundary conditions into pnt_config. We assume that all boundary
 * conditions satisfy PointSetType = PointRange and NormalListSize = 0. */
#if DEBUG
	printf( "--- loadFile STATUS ---\n"
		"Dimension:\t\t3\n"
		"Rank:\t\t\t%i\n"
		"Zone name:\t\t%s\n"
		"Given offset:\t\t[%i,%i" DIM3( ",%i" ) "]\n"
		"Size:\t\t\t[%i,%i" DIM3( ",%i" ) "]\n"
		"ID in mesh:\t\t%i\n"
		"Neigbour: Zone, Transform, Donor-Range, Range, Translation, Rotation Center, Rotation Angle\n",
		pnt_config->MPI_rank,
		zonename,
		(int)pnt_config->offset[ 0 ],(int)pnt_config->offset[ 1 ] DIM3( ,(int)pnt_config->offset[ 2 ] ),
		(int)pnt_config->zonesize[ 0 ],(int)pnt_config->zonesize[ 1 ] DIM3( ,(int)pnt_config->zonesize[ 2 ] ),
		zone
	);

	char sidenames[ ]={ 'i','j','k' };
	char endnames[ ][ 5 ]={ "low ","high" };
#endif

	for( i = 0; i<DIM; i++ )
		for( j = 0; j<2; j++ ) {
#if DEBUG
			printf( " %c-%s:\t",sidenames[ i ],endnames[ j ] );
#endif

			if( !neighbours[ i ][ j ].neighbour ) {
				int nbc,bc;
				CG( cg_nbocos( file,base,zone,&nbc ) );

				for( bc = 1; bc<=nbc; bc++ ) {
					PointSetType_t bcpst;
					char bcname[ 33 ];
					cgsize_t bcnpnts,bcnls;
					int bcndata,bcnorm[ DIM ];
					DataType_t bcdt;
					BCType_t bctype;
					cgsize_t bcpnts[ 2 ][ DIM ];

					CG( cg_boco_info( file,base,zone,bc,bcname,&bctype,&bcpst,&bcnpnts,bcnorm,&bcnls,&bcdt,&bcndata ) );
					CG( cg_boco_read( file,base,zone,bc,bcpnts[ 0 ],NULL ) );

// Correct orientation
					bool found = bcpnts[ 0 ][ i ]==bcpnts[ 1 ][ i ];
// Correct plane along normal
					found = found &&( bcpnts[ 0 ][ i ]==1 )==( j==0 );
// Intersects in tangent plane
					int dir;
					if( found )
						for( dir =( i+1 )%DIM; dir!=i; dir = ( dir + 1 )%DIM ) {
							cgsize_t low = 1+pnt_config->offset[ dir ];
							cgsize_t high = 1+pnt_config->zonesize[ dir ]+pnt_config->offset[ dir ];

							found = found &&( ( low>bcpnts[ 0 ][ dir ]&& low<bcpnts[ 1 ][ dir ] )||
							( high>bcpnts[ 0 ][ dir ]&& high<bcpnts[ 1 ][ dir ] )||
							( bcpnts[ 0 ][ dir ]>low && bcpnts[ 0 ][ dir ]<high )||
							( bcpnts[ 1 ][ dir ]>low && bcpnts[ 1 ][ dir ]<high ) );

						}

					if( found ) {
						BCType_t realbctype = bctype;

						if( bctype==FamilySpecified ) {
							char famname[ 33 ];
							CG( cg_goto( file,base,zonename,0,"ZoneBC",0,bcname,0,"end" ) );
							CG( cg_famname_read( famname ) );

							int nfams,fam;
							CG( cg_nfamilies( file,base,&nfams ) );

							for( fam = 1; fam<=nfams; fam++ ) {
								char tfamname[ 33 ];
								int nfbc,nfg;
								CG( cg_family_read( file,base,fam,tfamname,&nfbc,&nfg ) );
								if( !strcmp( famname,tfamname ) )
									break;
							}

							if( fam>nfams ) {
								printf( "(EE) Family '%s' could not be found in file\n",famname );
								exit( 1 );
							}

							char fambcname[ 33 ];
							CG( cg_fambc_read( file,base,fam,1,fambcname,&realbctype ) );
						}

						if( !assignBC( pnt_config,i,j,realbctype ) )
							exit( 1 );
						else {
#if DEBUG
							char (*const sidemap[ ][ 2 ])[ 30 ]= {
								{ &pnt_config->BC_Left,&pnt_config->BC_Right },
								{ &pnt_config->BC_Bottom,&pnt_config->BC_Top }
								DIM3( ,{ &pnt_config->BC_Behind,&pnt_config->BC_InFront } )
							};
							printf( "BC %s\n",*sidemap[ i ][ j ] );
#endif
							break;
						}
					}
				}

				if( bc>nbc ) {
					printf( "(EE) BC for side %i (%s end) could not be found\n",i,j?"higher":"lower" );
					exit( 1 );
				}
			} else {
#if DEBUG
				printf( "%i (%i,%i" DIM3( ",%i" )"), ",
					neighbours[ i ][ j ].neighbour,
					neighbours[ i ][ j ].transform[ 0 ],
					neighbours[ i ][ j ].transform[ 1 ]
					DIM3( ,neighbours[ i ][ j ].transform[ 2 ] )
				);
				printf( "[%i,%i" DIM3( ",%i" )"]-[%i,%i" DIM3( ",%i" )"], ",
					(int)( neighbours[ i ][ j ].drange[ 0 ][ 0 ] ),
					(int)( neighbours[ i ][ j ].drange[ 0 ][ 1 ] )
					DIM3( ,(int)( neighbours[ i ][ j ].drange[ 0 ][ 2 ] ) ),
					(int)( neighbours[ i ][ j ].drange[ 1 ][ 0 ] ),
					(int)( neighbours[ i ][ j ].drange[ 1 ][ 1 ] )
					DIM3( ,(int)( neighbours[ i ][ j ].drange[ 1 ][ 2 ] ) )
				);
				printf( "[%i,%i" DIM3( ",%i" )"]-[%i,%i" DIM3( ",%i" )"], ",
					(int)( neighbours[ i ][ j ].range[ 0 ][ 0 ] ),
					(int)( neighbours[ i ][ j ].range[ 0 ][ 1 ] )
					DIM3( ,(int)( neighbours[ i ][ j ].range[ 0 ][ 2 ] ) ),
					(int)( neighbours[ i ][ j ].range[ 1 ][ 0 ] ),
					(int)( neighbours[ i ][ j ].range[ 1 ][ 1 ] )
					DIM3( ,(int)( neighbours[ i ][ j ].range[ 1 ][ 2 ] ) )
				);
				printf( "(%g,%g" DIM3( ",%g" )"), ",
					neighbours[ i ][ j ].tr[ 0 ],
					neighbours[ i ][ j ].tr[ 1 ]
					DIM3( ,neighbours[ i ][ j ].tr[ 2 ] )
				);
				printf( "(%g,%g" DIM3( ",%g" )"), ",
					neighbours[ i ][ j ].rc[ 0 ],
					neighbours[ i ][ j ].rc[ 1 ]
					DIM3( ,neighbours[ i ][ j ].rc[ 2 ] )
				);
				printf( "(%g,%g" DIM3( ",%g" )")\n",
					neighbours[ i ][ j ].ra[ 0 ],
					neighbours[ i ][ j ].ra[ 1 ]
					DIM3( ,neighbours[ i ][ j ].ra[ 2 ] )
				);
#endif
			}
		}
	TM_END ( "SHOCK: Loading meshConfigs" )
	TM_START ( )
	const long long pointcount = pnt_config->zonesize[ 0 ]*pnt_config->zonesize[ 1 ] DIM3( *pnt_config->zonesize[ 2 ] );

/**
 * Load grid coordinates into VTAs for later processing. */

	int ncoords,coord;
	CG( cg_ncoords( file,base,zone,&ncoords ) );

	for( coord = 1; coord<=ncoords; coord++ ) {
		DataType_t dt;
		char coordname[ 33 ];
		vta* target = NULL;

		CG( cg_coord_info( file,base,zone,coord,&dt,coordname ) );

		if( !strcmp( coordname,"CoordinateX" ) )
			target = x;
		else if( !strcmp( coordname,"CoordinateY" ) )
			target = y;
#if DIM==3
		else if( !strcmp( coordname,"CoordinateZ" ) )
			target = z;
#endif

		if( target ) {
			target->ptr = malloc( pointcount*( dt==RealDouble?sizeof( float ):sizeof( float ) ) );
			target->dt = dt;
			#if DIM==2
				CG( cgp_coord_read_data( file,base,zone,coord,corners[ 0 ],corners[ 1 ],target->ptr ) );
			#endif
			#if ( HDF5_HAVE_MULTI_DATASETS_ == 0 && DIM == 3 )
				CG( cgp_coord_read_data( file,base,zone,coord,corners[ 0 ],corners[ 1 ],target->ptr ) );
			#endif
		}
	}
	#if ( HDF5_HAVE_MULTI_DATASETS_ == 1 && DIM == 3 )
		int coordVec[3];
		CG( cgp_coord_multi_read_data( file,base,zone,coordVec,corners[ 0 ],corners[ 1 ],x,y,z ) );
	#endif

/** Load solution data from last timestep into vta. If there is a
 * BaseIterativeData node, obtain the number of steps from it, find the
 * pointer to the according FlowSolution and read that FlowSolution. If
 * there is no BaseIterativeData, take the first FlowSolution node. If
 * there is no such node, self-destruct. */
	char bitername[ 33 ];
	int nsteps,nsols,sol = 1;
	CG( cg_nsols( file,base,zone,&nsols ) );

/** Find the index of the solution field or default to solution node 1.
 */
	if( !cg_biter_read( file,base,bitername,&nsteps ) ) {
		char zitername[ 33 ];
		int narrays,array;

		CG( cg_goto( file,base,"BaseIterativeData",0,"end" ) );
		CG( cg_narrays( &narrays ) );

		for( array = 1; array<=narrays; array++ ) {
			char arrayname[ 33 ];
			int dim;
			DataType_t dt;
			cgsize_t dimv[ 12 ];

			CG( cg_array_info( array,arrayname,&dt,&dim,dimv ) );

			void* timearray;

			if( !strcmp( arrayname,"TimeValues" ) ) {
				timearray = malloc( dimv[ 0 ]*( dt==RealDouble?sizeof( float ):sizeof( float ) ) );

				CG( cg_array_read( array,timearray ) );

				if( dt==RealDouble )
					pnt_config->flt_time_dim =( (float*)timearray )[ dimv[ 0 ]-1 ];
				else
					pnt_config->flt_time_dim =( (float*)timearray )[ dimv[ 0 ]-1 ];

				free( timearray );
				break;
			}
		}

		CG( cg_ziter_read( file,base,zone,zitername ) );
		CG( cg_goto( file,base,zonename,0,zitername,0,"end" ) );
		CG( cg_narrays( &narrays ) );
		for( array = 1; array<=narrays; array++ ) {
			char arrname[ 33 ];
			DataType_t dt;
			int dd;
			cgsize_t dv[ 12 ];

			CG( cg_array_info( array,arrname,&dt,&dd,dv ) );
			// A bit of redundancy despite any assumptions...
			if( !strcmp( arrname,"FlowSolutionPointers" )&& dt==Character && dd==2 && dv[ 0 ]==32 && dv[ 1 ]==nsteps ) {
				char* fsp = malloc( 32*nsteps );

				CG( cg_array_read( array,fsp ) );

				for( sol = 1; sol<=nsols; sol++ ) {
					char solname[ 33 ];
					GridLocation_t gl;

					CG( cg_sol_info( file,base,zone,sol,solname,&gl ) );

					if( !strncmp( solname,fsp+32*( nsteps-1 ),32 ) )
						break;
				}

				free( fsp );
			}
		}
	} else
		pnt_config->flt_time_dim = 0;
	TM_END ( "SHOCK: Loading mesh" )
	TM_START ( )
	if(( sol<=nsols )&&(abs(pnt_config->int_initializeType)==1)) {
/** Load solution data into vta for later processing. */
		int nfields,field,fieldVec[5];
		CG( cg_nfields( file,base,zone,sol,&nfields ) );
		for( field = 1; field<=nfields; field++ ) {
			DataType_t dt;
			char fieldname[ 33 ];
			vta* target = NULL;
			CG( cg_field_info( file,base,zone,sol,field,&dt,fieldname ) );

			// TODO: CGNS Conventions or not? u,v,w, etc are not.
			if( !strcmp( fieldname,"VelocityX" ) )
				{target = u;fieldVec[0]=field;}
			else if( !strcmp( fieldname,"VelocityY" ) )
				{target = v;fieldVec[1]=field;}
			else if( !strcmp( fieldname,"VelocityZ" ) )
				{target = w;fieldVec[2]=field;}
			else if( !strcmp( fieldname,"Density" ) )
				{target = rho;fieldVec[3]=field;}
			else if( !strcmp( fieldname,"Pressure" ) )
				{target = p;fieldVec[4]=field;}

			if( target ) {
				target->ptr = malloc( pointcount*( dt==RealDouble?sizeof( float ):sizeof( float ) ) );
				target->dt = dt;
				#if DIM==2
						if(pnt_config->MPI_rank==0){printf("SHOCK:  Using command 'cgp_field_read_data'\n");}
						CG( cgp_field_read_data( file,base,zone,sol,field,corners[ 0 ],corners[ 1 ],target->ptr) );
				#endif
				#if ( HDF5_HAVE_MULTI_DATASETS_ == 0 && DIM == 3 )
						if(pnt_config->MPI_rank==0){printf("SHOCK:  Using command 'cgp_field_read_data'\n");}
						CG( cgp_field_read_data( file,base,zone,sol,field,corners[ 0 ],corners[ 1 ],target->ptr) );
				#endif
			}
		}

		#if ( HDF5_HAVE_MULTI_DATASETS_ == 1 && DIM == 3 )
			if(pnt_config->MPI_rank==0){printf("SHOCK:  Using command 'cgp_field_multi_read_data'\n");}
			CG( cgp_field_multi_read_data( file,base,zone,sol,fieldVec,corners[ 0 ],corners[ 1 ],5,u->ptr,v->ptr,w->ptr,rho->ptr,p->ptr) );
		#endif

	}
	TM_END ( "SHOCK: Loading results" )
	TM_START ( )
	cgp_close( file );
	TM_END ( "SHOCK: Closing meshfile" )
	TM_START ( )
	return 0;
}
