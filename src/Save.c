#include "Save.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <pcgnslib.h>

#define CG(cmd) if(cmd)cg_error_exit( );
#define DIM MESHDIMENSIONS
#define PHYSDIM 3
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

void saveFile( struct strct_configuration* pnt_config,struct strct_Film* pnt_Film ) {

/** Cleaning prepares the file such that all solution nodes, all
 * ZoneIterativeData and BaseIterativeData are deleted */

	int nprocs;
	MPI_Comm_size( MPI_COMM_WORLD,&nprocs );
	bool worldleader = nprocs==1 || pnt_config->MPI_rank==0;

	char* zonename = pnt_config->Zonename;
/** Make the outfile the infile, if it is not already. At this point, we
 * assume that the outfile is the infile and do not do anything. */

	int file,base = 1;
	int zone,nzones;
	if( worldleader ) {
		TM_START( )
		CG( cg_open( pnt_config->chr_MeshPath,CG_MODE_MODIFY,&file ) );
		TM_END( "SHOCK: Opening savefile - rank0 - preparing Savefile" )
		if( abs(pnt_config->int_initializeType)!=1 ) {
			CG( cg_nzones( file,base,&nzones ) );

			for( zone = 1; zone<=nzones; zone++ ) {
				int nsols,sol;
				CG( cg_nsols( file,base,zone,&nsols ) );
				CG( cg_goto( file,base,"Zone_t",zone,"end" ) );
				for( sol = 1; sol<=nsols; sol++ ) {
					GridLocation_t gl;
					char solname[ 33 ];
					CG( cg_sol_info( file,base,zone,1,solname,&gl ) );
					CG( cg_delete_node( solname ) );
				}
				cg_delete_node( "ZoneIterativeData" );
			}

			CG( cg_goto( file,base,"end" ) );
			cg_delete_node( "BaseIterativeData" );
		}
/** Adjust number of iterations */
		CG( cg_goto( file,base,"end" ) );
		char iterationstr[ 255 ];
		snprintf( iterationstr,255,"%i",pnt_config->int_actualIteration-1 );
		int ndesc,desc;

		CG( cg_ndescriptors( &ndesc ) );

		for( desc = 1; desc<=ndesc; desc++ ) {
			char* dtext,dname[ 33 ];
			CG( cg_descriptor_read( desc,dname,&dtext ) );
			if( !strcmp( dname,"Iterations" ) ) {
				CG( cg_delete_node( "Iterations" ) );
				cg_free( dtext );
				break;
			} else
				cg_free( dtext );
		}

		CG( cg_descriptor_write( "Iterations",iterationstr ) );

		int nsteps; ///< Number of steps already in the file
		char bitername[ 33 ];
		float* timearray;

		if( cg_biter_read( file,base,bitername,&nsteps ) ) {
			nsteps = 0;
			CG( cg_biter_write( file,base,"BaseIterativeData",pnt_config->int_Samples ) );

			timearray = malloc( pnt_config->int_Samples*sizeof( float ) );

			CG( cg_goto( file,base,"BaseIterativeData",0,"end" ) );
		} else {
			CG( cg_goto( file,base,"BaseIterativeData",0,"end" ) );
			int array,narrays;
			CG( cg_narrays( &narrays ) );

			timearray = malloc( ( pnt_config->int_Samples+nsteps )*sizeof( float ) );

			for( array = 1; array<=narrays; array++ ) {
				char arrayname[ 33 ];
				int dim;
				DataType_t dt;
				cgsize_t dimv[ 12 ];

				CG( cg_array_info( array,arrayname,&dt,&dim,dimv ) );

				if( !strcmp( arrayname,"TimeValues" ) ) {
					if( dt==RealDouble ) {
						float* buffer = malloc( nsteps*sizeof( float ) );

						CG( cg_array_read( array,buffer ) );

						int i;
						for( i = 0; i<nsteps; i++ )
							timearray[ i ]= buffer[ i ];

						free( buffer );
					} else
						CG( cg_array_read( array,timearray ) );
					break;
				}
			}

			CG( cg_goto( file,base,"end" ) );
			CG( cg_delete_node( bitername ) );

			CG( cg_biter_write( file,base,"BaseIterativeData",nsteps+pnt_config->int_Samples ) );
			CG( cg_goto( file,base,"BaseIterativeData",0,"end" ) );
		}

		int i;

		for( i = 0; i<pnt_config->int_Samples; i++ )
			timearray[ i+nsteps ]= pnt_Film->flt_time_dim[ i ];

		CG( cg_array_write( "TimeValues",RealSingle,1,( cgsize_t[ ] ){ nsteps+pnt_config->int_Samples },timearray ) );

		free( timearray );

		CG( cg_nzones( file,base,&nzones ) );

/** Responsibility of creating all the rest in the zones but the actual
 * field data, which will be written by each process individually. */
		for( zone = 1; zone<=nzones; zone++ ) {
			char zonename[ 33 ];
			cgsize_t zonesize[ 3 ][ DIM ];

			CG( cg_zone_read( file,base,zone,zonename,zonesize[ 0 ] ) );

			int nsols;
			CG( cg_nsols( file,base,zone,&nsols ) );

			char zitername[ 33 ];

			int prepends = 0;
			char* solnames = NULL;

			if( cg_ziter_read( file,base,zone,zitername ) ) {
				CG( cg_ziter_write( file,base,zone,"ZoneIterativeData" ) );

				solnames = malloc( 32*( pnt_config->int_Samples ) );
				CG( cg_goto( file,base,zonename,0,"ZoneIterativeData",0,"end" ) );
			} else {
				CG( cg_goto( file,base,zonename,0,"ZoneIterativeData",0,"end" ) );
				int narrays,array;

				CG( cg_narrays( &narrays ) );
				for( array = 1; array<=narrays; array++ ) {
					char arrayname[ 33 ];
					int dim;
					DataType_t dt;
					cgsize_t dimv[ 12 ];

					CG( cg_array_info( array,arrayname,&dt,&dim,dimv ) );

					if( !strcmp( arrayname,"FlowSolutionPointers" ) ) {
						prepends = dimv[ 1 ];
						solnames = malloc( 32*( dimv[ 1 ]+pnt_config->int_Samples ) );

						CG( cg_array_read( array,solnames ) );

						CG( cg_delete_node( arrayname ) );
						break;
					}
				}

				if( array>narrays )
					solnames = malloc( 32*( pnt_config->int_Samples ) );
			}

			int sol;
			int n_check;
			if (pnt_config->flag_ReducedExport==0){n_check = 8;}
			else{n_check = 5;}

			for( sol = 0; sol<pnt_config->int_Samples; sol++ ) {
				int n;
				char solname[ 33 ];

				if( snprintf( solname,33,"Sol_%04d",sol+nsols )>33 ) {
					printf( "(EE) Could not write solution #%i because consecutive numbering exceeded string length\n",sol );
					exit( 1 );
				}

				memcpy( solnames+prepends*32+sol*32,solname,32 );

				CG( cg_sol_write( file,base,zone,solname,Vertex,&n ) );
				// Just a quick check whether things worked as expected
				if( n!=nsols+1+sol ) {
					printf( "(EE) Order of solutions did not emerge as expected appended %ith solution to %i solutions, got solution id %i\n",sol+1,nsols,n );
					exit( 1 );
				}

				CG( cg_field_write( file,base,zone,nsols+1+sol,RealSingle,"VelocityX",NULL,&n ) );
				CG( cg_field_write( file,base,zone,nsols+1+sol,RealSingle,"VelocityY",NULL,&n ) );
				CG( cg_field_write( file,base,zone,nsols+1+sol,RealSingle,"VelocityZ",NULL,&n ) );
				CG( cg_field_write( file,base,zone,nsols+1+sol,RealSingle,"Pressure",NULL,&n ) );
				CG( cg_field_write( file,base,zone,nsols+1+sol,RealSingle,"Density",NULL,&n ) );
				if(pnt_config->flag_ReducedExport==0)
				{
					CG( cg_field_write( file,base,zone,nsols+1+sol,RealSingle,"DensityGradient",NULL,&n ) );
					CG( cg_field_write( file,base,zone,nsols+1+sol,RealSingle,"MachNumber",NULL,&n ) );
					CG( cg_field_write( file,base,zone,nsols+1+sol,RealSingle,"Lambda2",NULL,&n ) );
				}

				// Just a quick check whether things worked as expected
				if( n!=n_check ) {
					printf( "(EE) Order of solution fields did not emerge as expected\n" );
					exit( 1 );
				}
			}

			CG( cg_array_write( "FlowSolutionPointers",Character,2,(cgsize_t[ ]){ 32,prepends+pnt_config->int_Samples },solnames ) );

			free( solnames );
		}

		CG( cg_close( file ) );
	}

	MPI_Barrier( MPI_COMM_WORLD );
	TM_START( )
	CG( cgp_open( pnt_config->chr_MeshPath,CG_MODE_MODIFY,&file ) );
//#if (CGNS_VERSION<=3210)
//	CG( cgp_pio_mode(CGP_COLLECTIVE));
//#else
//	CG( cgp_pio_mode(CGP_COLLECTIVE, MPI_INFO_NULL));
//#endif
//	CG( cgp_queue_set( 1 ) );
	TM_END( "SHOCK: Opening savefile" )

	CG( cg_nzones( file,base,&nzones ) );
	for( zone = 1; zone<=nzones; zone++ ) {
		char tzname[ 33 ];
		cgsize_t zonesize[ 3 ][ DIM ];
		CG( cg_zone_read( file,base,zone,tzname,zonesize[ 0 ] ) );
		if( !strcmp( tzname,zonename ) )
			break;
	}
	if( zone>nzones ) {
		printf( "(EE) Could not find associated zone '%s' in the mesh file\n",zonename );
		exit( 1 );
	}

	CG( cg_goto( file,base,zonename,0,"ZoneIterativeData",0,"end" ) );
	int narrays,array;
	char* solnames;
	int tcount;

	CG( cg_narrays( &narrays ) );
	for( array = 1; array<=narrays; array++ ) {
		char arrayname[ 33 ];
		int dim;
		DataType_t dt;
		cgsize_t dimv[ 12 ];

		CG( cg_array_info( array,arrayname,&dt,&dim,dimv ) );

		if( !strcmp( arrayname,"FlowSolutionPointers" ) ) {
			tcount = dimv[ 1 ];
			solnames = malloc( 32*tcount );
			CG( cg_array_read( array,solnames ) );
			break;
		}
	}

	if( array>narrays ) {
		printf( "(EE) Can't find solution pointers for zone '%s'\n",zonename );
		exit( 1 );
	}

	int nsols,sol;
	CG( cg_nsols( file,base,zone,&nsols ) );

	const int imax = pnt_config->zonesize[ 0 ];
	const int jmax = pnt_config->zonesize[ 1 ];
#if DIM==3
	const int kmax = pnt_config->zonesize[ 2 ];
#endif

	size_t bufsize = sizeof( float )*imax*jmax DIM3( *kmax );

	float* u = malloc	( bufsize ),
		* v = malloc( bufsize ),
		* w = malloc( bufsize ),
		* p = malloc( bufsize ),
		* rho = malloc( bufsize ),
		* gradrho = malloc( bufsize ),
		* ma = malloc( bufsize ),
		* lambda = malloc( bufsize );

	for( sol = 1; sol<=nsols; sol++ ) {
		int t,tindex;
		char solname[ 33 ];
		GridLocation_t gl;
		CG( cg_sol_info( file,base,zone,sol,solname,&gl ) );

		for( t = tcount-pnt_config->int_Samples; t<tcount; t++ )
			if( !strncmp( solnames+32*t,solname,32 ) )
				break;

		// Solution not listed at all or has already been there
		if( t>=tcount )
			continue;

		tindex = t-( tcount-pnt_config->int_Samples );

/* Write data from film into suitable memory without ghost cells. Map
 * of memory in film: ...[t][i][j][k]. Map of memory as required by as
 * required by CGNS Midlevel: [k][j][i] for each t (this is [i][j][k] in
 * Fotran memory mapping). */

		int i;
		for( i = 0; i<pnt_config->zonesize[ 0 ]; i++ ) {
			int j;
			for( j = 0; j<pnt_config->zonesize[ 1 ]; j++ ) {
#if DIM==3
				int k;
				for( k = 0; k<pnt_config->zonesize[ 2 ]; k++ ) {
					int memindex = tindex*imax*jmax*kmax+i*jmax*kmax+j*kmax+k;
					int cgindex = k*imax*jmax+j*imax+i;
#else
				{
					int memindex = tindex*imax*jmax+i*jmax+j;
					int cgindex = j*imax+i;
#endif
					u[ cgindex ]= pnt_Film->u[ memindex ];
					v[ cgindex ]= pnt_Film->v[ memindex ];
					w[ cgindex ]= pnt_Film->w[ memindex ];
					rho[ cgindex ]= pnt_Film->rho[ memindex ];
					p[ cgindex ]= pnt_Film->p[ memindex ];
					gradrho[ cgindex ]= pnt_Film->gradRho[ memindex ];
					ma[ cgindex ]= pnt_Film->MachNumber[ memindex ];
					lambda[ cgindex ]= pnt_Film->Lambda2[ memindex ];
// Just for correct syntax highlighting
#if DIM==3
				}
#else
				}
#endif
			}
		}

		cgsize_t from[ ]={ pnt_config->offset[ 0 ]+1,pnt_config->offset[ 1 ]+1 DIM3( ,pnt_config->offset[ 2 ]+1 ) };
		cgsize_t to[ ]={ from[ 0 ]+imax-1,from[ 1 ]+jmax-1 DIM3( ,from[ 2 ]+kmax-1 ) };
		int field,nfields;
		if (pnt_config->flag_ReducedExport==0){nfields = 8;}
		else{nfields = 5;}
		int *fieldVec;
		fieldVec=malloc(nfields*(sizeof(int)));

		for( field = 1; field<=nfields; field++ ) {
			char fieldname[ 33 ];
			DataType_t dt;
			CG( cg_field_info( file,base,zone,sol,field,&dt,fieldname ) );

			#if DIM==2
				if((pnt_config->MPI_rank==0)&&(field==1)&&(sol==1)){printf("SHOCK:  using cgp_field_write_data (2D).\n");}
				if( !strcmp( fieldname,"VelocityX" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,u ) );fieldVec[0]=field;
				} else if( !strcmp( fieldname,"VelocityY" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,v ) );fieldVec[1]=field;
				} else if( !strcmp( fieldname,"VelocityZ" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,w ) );fieldVec[2]=field;
				} else if( !strcmp( fieldname,"Density" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,rho ) );fieldVec[3]=field;
				} else if( !strcmp( fieldname,"Pressure" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,p ) );fieldVec[4]=field;
				} else if( !strcmp( fieldname,"DensityGradient" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,gradrho ) );fieldVec[5]=field;
				} else if( !strcmp( fieldname,"MachNumber" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,ma ) );fieldVec[6]=field;
				} else if( !strcmp( fieldname,"Lambda2" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,lambda ) );fieldVec[7]=field;
				}
			#endif

			#if ( HDF5_HAVE_MULTI_DATASETS_ == 0 && DIM == 3 )
				if((pnt_config->MPI_rank==0)&&(field==1)&&(sol==1)){printf("SHOCK:  using cgp_field_write_data (3D-Hdf5<1.8.14).\n");}
				if( !strcmp( fieldname,"VelocityX" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,u ) );fieldVec[0]=field;
				} else if( !strcmp( fieldname,"VelocityY" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,v ) );fieldVec[1]=field;
				} else if( !strcmp( fieldname,"VelocityZ" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,w ) );fieldVec[2]=field;
				} else if( !strcmp( fieldname,"Density" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,rho ) );fieldVec[3]=field;
				} else if( !strcmp( fieldname,"Pressure" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,p ) );fieldVec[4]=field;
				} else if( !strcmp( fieldname,"DensityGradient" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,gradrho ) );fieldVec[5]=field;
				} else if( !strcmp( fieldname,"MachNumber" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,ma ) );fieldVec[6]=field;
				} else if( !strcmp( fieldname,"Lambda2" ) ) {
					CG( cgp_field_write_data( file,base,zone,sol,field,from,to,lambda ) );fieldVec[7]=field;
				}
			#endif

		}
		#if ( HDF5_HAVE_MULTI_DATASETS_ == 1 && DIM == 3 )
			if((pnt_config->MPI_rank==0)&&(field==1)){printf("SHOCK:  using cgp_field_multi_write_data.\n");}
			if(nfields==5)
			{	CG( cgp_field_multi_write_data( file,base,zone,sol,fieldVec,from,to,nfields,u,v,w,rho,p ) );}
			if(nfields==8)
			{	CG( cgp_field_multi_write_data( file,base,zone,sol,fieldVec,from,to,nfields,u,v,w,rho,p,gradrho,ma,lambda ) );}
		#endif
	}
//	TM_START( )
//	CG(cgp_queue_flush()); flush is not recommended in improved pcgns-I/O
//	TM_END( "cgp_queue_flush" )

	free( u );
	free( v );
	free( w );
	free( rho );
	free( p );
	free( gradrho );
	free( ma );
	free( lambda );

	CG( cgp_close( file ) );
}
