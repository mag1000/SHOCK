#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "cgnslib.h"

#define MESHDIM 3  // CAN BE VARIED
#define PHYSDIM 3  // HAS TO BE 3

int main(int argc, char *argv[])
{
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//	Division-File Creator for Only 3D C-Grid (1 domain)
//	Sinve Version 3 With decomposition in k-direction
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

	char path[300];
	char zonename[32];
	sprintf(path,"3d_NACA0012_yMin20_500000_8192x1024x1024.sddc");
//	sprintf(zonename,"dom-2");
	sprintf(zonename,"blk-2");
	int len_zonename=5;
	int i_size_mesh=8192;
	int j_size_mesh=1024;
	int k_size_mesh=1024;
	
	
	int i_size_CPU=32;
	int j_size_CPU=64;
	int k_size_CPU=64;
	
	int i_parts=i_size_mesh/i_size_CPU;
	int j_parts=j_size_mesh/j_size_CPU;
	int k_parts=k_size_mesh/k_size_CPU;
	
	int i_points_connected=1024;		// number of grid points which are connected downstream of trailing edge
	

	
	printf("Das %dD-Gitter (%dx%dx%d) wird für %d(%dx%dx%d) CPU aufgeteilt.\nDie zerteilten Gitter sind %dx%dx%d groß.\n",
	MESHDIM,i_size_mesh,j_size_mesh,k_size_mesh,(i_parts*j_parts*k_parts),i_parts,j_parts,k_parts,i_size_CPU,j_size_CPU,k_size_CPU);
	
	FILE* divfile;
	divfile=fopen(path,"w");
	int i,j,k,rank;

	int index,index_limit;
	size_t ifblocklen = sizeof( int )+MESHDIM*sizeof( int )+2*MESHDIM*sizeof( cgsize_t )*2+PHYSDIM*sizeof( float )*3;
	size_t datalen = 32+sizeof( cgsize_t )*MESHDIM*2+ifblocklen*2*MESHDIM;
	
	printf( "Writing division file, %i bytes per zone, %i bytes per interface, sizeof(cgsize_t)=%zu\n",(int)datalen,(int)ifblocklen,sizeof(cgsize_t) );

/** Adresses in dataset are arranged as following, assuming
 * big-endianess when serializing to and from file:
 * 
 * zonename(32*char)
 * cornerlo(DIM*cgsize_t),
 * cornerhi(DIM*cgsize_t),
 * 2*DIM*neighbour(
 * 	zoneid(int),
 * 	transform(DIM*int),
 * 	rcs(PHYSDIM*float),
 * 	ras(PHYSDIM*float),
 * 	trs(PHYSDIM*float)
 * )
 */
	char* dataset = calloc( 1,datalen );
	for( k = 0; k<k_parts; k++ ) 
	{
	for( j = 0; j<j_parts; j++ ) 
	{
	for( i = 0; i<i_parts; i++ ) 
	{
		rank=k*i_parts*j_parts+j*i_parts+i;
		memcpy( dataset,zonename,len_zonename );

		cgsize_t(* corners)[ MESHDIM ]= (cgsize_t(*)[ MESHDIM ])( dataset+32 );
		corners[0][0]=i*i_size_CPU+1;
		corners[0][1]=j*j_size_CPU+1;
		corners[0][2]=k*k_size_CPU+1;
		
		corners[1][0]=(i+1)*i_size_CPU;
		corners[1][1]=(j+1)*j_size_CPU;
		corners[1][2]=(k+1)*k_size_CPU;

		
		index_limit=6;
		for( index = 0; index<index_limit; index++ )
		{
			int* neighbour= (int*)( dataset+32+MESHDIM*2*sizeof( cgsize_t )+index*ifblocklen );
			int* transform= neighbour+1;
			cgsize_t(* outranges)[ MESHDIM ]= (cgsize_t(*)[ MESHDIM ])( transform+MESHDIM ),
				(* outdranges)[ MESHDIM ]= outranges+2;
			float* rcs = (float*)( outdranges+2 ),
				* ras = rcs+PHYSDIM,
				* trs = ras+PHYSDIM;
				
			ras[0]=0.;
			ras[1]=0.;
			ras[2]=0.;			
						
			rcs[0]=0.;
			rcs[1]=0.;
			rcs[2]=0.;						

			trs[0]=0.;
			trs[1]=0.;
			trs[2]=0.;

			switch( index ) 
			{
			case 0: //	linker Nachbar
				if((rank%i_parts)==0)// linker oder rechter Rand
				{
					*neighbour=0;
				}
				else
				{
					*neighbour=rank-1+1;
				}
				transform[0]=1;
				transform[1]=2;
				transform[2]=3;
				outranges[0][0]=1;
				outranges[0][1]=1;
				outranges[0][2]=1;
				outranges[1][0]=1;
				outranges[1][1]=j_size_CPU;
				outranges[1][2]=k_size_CPU;
				outdranges[0][0]=i_size_CPU;
				outdranges[0][1]=1;
				outdranges[0][2]=1;
				outdranges[1][0]=i_size_CPU;
				outdranges[1][1]=j_size_CPU;
				outdranges[1][2]=k_size_CPU;
				break;
			case 1: //	rechter Nachbar
				if((rank%i_parts)==(i_parts-1))
				{
					*neighbour=0;
				}
				else
				{
					*neighbour=rank+1+1;
				}
				transform[0]=1;
				transform[1]=2;
				transform[2]=3;
				outranges[0][0]=i_size_CPU;
				outranges[0][1]=1;
				outranges[0][2]=1;
				outranges[1][0]=i_size_CPU;
				outranges[1][1]=j_size_CPU;
				outranges[1][2]=k_size_CPU;
				outdranges[0][0]=1;
				outdranges[0][1]=1;
				outdranges[0][2]=1;
				outdranges[1][0]=1;
				outdranges[1][1]=j_size_CPU;
				outdranges[1][2]=k_size_CPU;
				break;
			case 2: //	unterer Nachbar
				if(j==0)
				{
				if(
				(i>=(i_points_connected/i_size_CPU))
				&&
				(i<(i_parts-(i_points_connected/i_size_CPU)))
				)
				 //	Abfrage wegen Verbindungen hinter Profil im wake
				{
				*neighbour=0;	
				}			
				else
				{
				//printf( "Selfcon i=%i j=%i\n",i,j );
				*neighbour=i_parts-i+k*i_parts*j_parts;
				transform[0]=-1;
				transform[1]=-2;
				transform[2]=3;
				
				outranges[0][0]=1;
				outranges[0][1]=1;
				outranges[0][2]=1;
				outranges[1][0]=i_size_CPU;
				outranges[1][1]=1;
				outranges[1][2]=k_size_CPU;

				outdranges[0][0]=1;
				outdranges[0][1]=1;
				outdranges[0][2]=1;
				outdranges[1][0]=i_size_CPU;
				outdranges[1][1]=1;
				outdranges[1][2]=k_size_CPU;			
				}
				}
				else //	normaler unterer Nachbar
				{
				*neighbour=rank-i_parts+1;
				transform[0]=1;
				transform[1]=2;
				transform[2]=3;	
				
				outranges[0][0]=1;
				outranges[0][1]=1;
				outranges[0][2]=1;
				outranges[1][0]=i_size_CPU;
				outranges[1][1]=1;
				outranges[1][2]=k_size_CPU;
									
				outdranges[0][0]=1;
				outdranges[0][1]=j_size_CPU;
				outdranges[0][2]=1;
				outdranges[1][0]=i_size_CPU;
				outdranges[1][1]=j_size_CPU;
				outdranges[1][2]=k_size_CPU;
				}												
				break;
			case 3: //	oberer Nachbar
				if(j==j_parts-1) //	Abfrage oberer Rand
				{
				*neighbour=0;
				}
				else
				{
				*neighbour=rank+i_parts+1;
				}
				transform[0]=1;
				transform[1]=2;
				transform[2]=3;	
				
				outranges[0][0]=1;
				outranges[0][1]=j_size_CPU;
				outranges[0][2]=1;
				outranges[1][0]=i_size_CPU;
				outranges[1][1]=j_size_CPU;
				outranges[1][2]=k_size_CPU;
									
				outdranges[0][0]=1;
				outdranges[0][1]=1;
				outdranges[0][2]=1;
				outdranges[1][0]=i_size_CPU;
				outdranges[1][1]=1;
				outdranges[1][2]=k_size_CPU;
				break;
			case 4: //	hinterer Nachbar
				if(k==0)
				{
				*neighbour=rank+i_parts*j_parts*(k_parts-1)+1;
				trs[2]=-(0.1+0.1/(float)(k_size_mesh-1));
				}
				else
				{
				*neighbour=rank-i_parts*j_parts+1;
				}
				transform[0]=1;
				transform[1]=2;
				transform[2]=3;	
				
				outranges[0][0]=1;
				outranges[0][1]=1;
				outranges[0][2]=1;
				outranges[1][0]=i_size_CPU;
				outranges[1][1]=j_size_CPU;
				outranges[1][2]=1;
									
				outdranges[0][0]=1;
				outdranges[0][1]=1;
				outdranges[0][2]=k_size_CPU;
				outdranges[1][0]=i_size_CPU;
				outdranges[1][1]=j_size_CPU;
				outdranges[1][2]=k_size_CPU;
				
				break;
			case 5: //	vorderer Nachbar
				if(k==k_parts-1)
				{
				*neighbour=rank-i_parts*j_parts*(k_parts-1)+1;
				trs[2]=(0.1+0.1/(float)(k_size_mesh-1));
				}
				else
				{
				*neighbour=rank+i_parts*j_parts+1;
				}			
				transform[0]=1;
				transform[1]=2;
				transform[2]=3;	

				outranges[0][0]=1;
				outranges[0][1]=1;
				outranges[0][2]=k_size_CPU;
				outranges[1][0]=i_size_CPU;
				outranges[1][1]=j_size_CPU;
				outranges[1][2]=k_size_CPU;
									
				outdranges[0][0]=1;
				outdranges[0][1]=1;
				outdranges[0][2]=1;
				outdranges[1][0]=i_size_CPU;
				outdranges[1][1]=j_size_CPU;
				outdranges[1][2]=1;	

				break;
			}

			if((*neighbour<1)||(*neighbour>(i_parts*j_parts*k_parts)))
				{*neighbour=0;}
//if (index==3)
	//	printf("rank %d top%d\n",rank,*neighbour);
			
				
		}

		fwrite( dataset,datalen,1,divfile );

	}
	}
	}
}


