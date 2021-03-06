//#define _ISOC99_SOURCE  1

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "float.h"
#include "string.h"
#include "cgnslib.h"
//make and remove directories
#include <unistd.h>
#include <sys/types.h>  /* Linux/UNIX */
#include <sys/stat.h>   /* Linux/UNIX */

#include "SHOCK.h"
#include "Functions.h"
#include "Import.h"
#include "Export.h"
#include "WENO.h"
#include "ZD.h"
#include "BC.h"
#include "ManufacturedSolution.h"

int i,j,k,ijk;

void postprocessLoad(
	vta* x,vta* y,vta* z,vta* u,vta* v,vta* w,vta* rho,vta* p,
	struct strct_configuration * pnt_config,
	struct strct_mesh * pnt_mesh,
	struct strct_U * pnt_U_lastStep,
	struct strct_U * pnt_U_RK,
	struct strct_Flux * pnt_Flux,
	struct strct_Flux * pnt_Flux_PlusHalf,
	struct strct_Flux * pnt_Q,
	struct strct_Flux * pnt_Q_sum,
	struct strct_Film * pnt_Film,
	struct strct_U * pnt_U_backup1,
	struct strct_U * pnt_U_backup2)
{
	int i2,j2,k2,ijk2;

	pnt_config->int_meshDimensions=MESHDIMENSIONS;

	if(abs(pnt_config->int_initializeType)!=1)
	{
		pnt_config->time_dim=0;
		pnt_config->int_StartIteration=0;
	}

//	Da in Load.c StartIteration gesetzt wurde (Weiterrechnen),
//	muessen nun folgende zwei Variablen angepasst werden
	pnt_config->int_actualIteration=pnt_config->int_StartIteration;
	pnt_config->int_EndIteration=pnt_config->int_StartIteration+pnt_config->int_TotalIterations;


	pnt_config->int_iMeshPoints=pnt_config->zonesize[0];
	pnt_config->int_jMeshPoints=pnt_config->zonesize[1];
	pnt_config->int_iMeshPointsGhostCells=pnt_config->int_iMeshPoints+(pnt_config->int_SpaceOrder+1);
	pnt_config->int_jMeshPointsGhostCells=pnt_config->int_jMeshPoints+(pnt_config->int_SpaceOrder+1);

	if(MESHDIMENSIONS==2)
	{
		pnt_config->int_kMeshPoints=1;
		pnt_config->int_kMeshPointsGhostCells=1;
	}
	else
	{
		pnt_config->int_kMeshPoints=pnt_config->zonesize[2];
		pnt_config->int_kMeshPointsGhostCells=pnt_config->zonesize[2]+(pnt_config->int_SpaceOrder+1);
	}


	pnt_config->int_iStartGhosts=0;
	pnt_config->int_jStartGhosts=0;
	pnt_config->int_kStartGhosts=0;
	pnt_config->int_iEndGhosts=pnt_config->int_iStartGhosts+pnt_config->int_iMeshPoints+pnt_config->int_SpaceOrder;
	pnt_config->int_jEndGhosts=pnt_config->int_jStartGhosts+pnt_config->int_jMeshPoints+pnt_config->int_SpaceOrder;
	pnt_config->int_kEndGhosts=pnt_config->int_kStartGhosts+pnt_config->int_kMeshPoints+pnt_config->int_SpaceOrder;

	pnt_config->int_iStartReal=(pnt_config->int_SpaceOrder-1)/2+1;
	pnt_config->int_jStartReal=(pnt_config->int_SpaceOrder-1)/2+1;
	pnt_config->int_kStartReal=(pnt_config->int_SpaceOrder-1)/2+1;
	pnt_config->int_iEndReal=pnt_config->int_iStartReal+pnt_config->int_iMeshPoints-1;
	pnt_config->int_jEndReal=pnt_config->int_jStartReal+pnt_config->int_jMeshPoints-1;
	pnt_config->int_kEndReal=pnt_config->int_kStartReal+pnt_config->int_kMeshPoints-1;

	pnt_config->int_iMid=(int)(0.5*(pnt_config->int_iEndGhosts+pnt_config->int_iStartGhosts));
	pnt_config->int_jMid=(int)(0.5*(pnt_config->int_jEndGhosts+pnt_config->int_jStartGhosts));
	pnt_config->int_kMid=(int)(0.5*(pnt_config->int_kEndGhosts+pnt_config->int_kStartGhosts));

	pnt_config->ijkMid=pnt_config->int_iMid*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_jMid*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_kMid;

	if(MESHDIMENSIONS==2)
	{
		pnt_config->int_kStartReal=0;
		pnt_config->int_kEndReal=0;

		pnt_config->int_kStartGhosts=0;
		pnt_config->int_kEndGhosts=0;

		pnt_config->int_kMid=0;

	}

	pnt_config->int_kStartReal_original=pnt_config->int_kStartReal;
	pnt_config->int_kEndReal_original=pnt_config->int_kEndReal;

//	Extrapolate
	pnt_config->int_iStartReal_extrapolate=pnt_config->int_iStartReal+1;
	pnt_config->int_jStartReal_extrapolate=pnt_config->int_jStartReal+1;
	pnt_config->int_kStartReal_extrapolate=pnt_config->int_kStartReal+1;
	pnt_config->int_iEndReal_extrapolate=pnt_config->int_iEndReal+1;
	pnt_config->int_jEndReal_extrapolate=pnt_config->int_jEndReal+1;
	pnt_config->int_kEndReal_extrapolate=pnt_config->int_kEndReal+1;

	pnt_config->int_iStartGhosts_extrapolate=pnt_config->int_iStartGhosts;
	pnt_config->int_jStartGhosts_extrapolate=pnt_config->int_jStartGhosts;
	pnt_config->int_kStartGhosts_extrapolate=pnt_config->int_kStartGhosts;
	pnt_config->int_iEndGhosts_extrapolate=pnt_config->int_iEndGhosts+2;
	pnt_config->int_jEndGhosts_extrapolate=pnt_config->int_jEndGhosts+2;
	pnt_config->int_kEndGhosts_extrapolate=pnt_config->int_kEndGhosts+2;

	pnt_config->int_iMeshPointsGhostCells_extrapolate=pnt_config->int_iMeshPointsGhostCells+2;
	pnt_config->int_jMeshPointsGhostCells_extrapolate=pnt_config->int_jMeshPointsGhostCells+2;
	pnt_config->int_kMeshPointsGhostCells_extrapolate=pnt_config->int_kMeshPointsGhostCells+2;

	if(MESHDIMENSIONS==2)
	{
		pnt_config->int_kStartReal_extrapolate=0;
		pnt_config->int_kEndReal_extrapolate=0;

		pnt_config->int_kStartGhosts_extrapolate=0;
		pnt_config->int_kEndGhosts_extrapolate=0;

		pnt_config->int_kMeshPointsGhostCells_extrapolate=1;

	}


//  Speicherallokierung für lokales Gamma (Prandtl)
    pnt_config->Gamma=(FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
                pnt_config->Gamma[ijk]=1.0/((pnt_config->gammaNumber-1.0)*pow(pnt_config->machNumber,2.0)*pnt_config->reynoldsNumber*pnt_config->prandtlNumber);
            }
        }
    }


	//		Transfer Interface-Informations
//		getInterfaceInformations(
//				pnt_config);

	//	pnt_config->bufferSendFlowLeft = (FLT *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints,sizeof(FLT));
	//	pnt_config->bufferSendFlowRight= (FLT *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints,sizeof(FLT));
	//	pnt_config->bufferSendFlowBottom = (FLT *)calloc(5*pnt_config->int_iMeshPoints * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPoints,sizeof(FLT));
	//	pnt_config->bufferSendFlowTop= (FLT *)calloc(5*pnt_config->int_iMeshPoints * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPoints,sizeof(FLT));
	//	pnt_config->bufferRecieveFlowLeft = (FLT *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints,sizeof(FLT));
	//	pnt_config->bufferRecieveFlowRight= (FLT *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints,sizeof(FLT));
	//	pnt_config->bufferRecieveFlowBottom = (FLT *)calloc(5*pnt_config->int_iMeshPoints * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPoints,sizeof(FLT));
	//	pnt_config->bufferRecieveFlowTop= (FLT *)calloc(5*pnt_config->int_iMeshPoints * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPoints,sizeof(FLT));

		pnt_config->bufferSendMeshLeft = (FLT *)calloc(13*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
		pnt_config->bufferSendMeshRight = (FLT *)calloc(13*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
		pnt_config->bufferSendMeshBottom = (FLT *)calloc(13*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
		pnt_config->bufferSendMeshTop = (FLT *)calloc(13*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
		pnt_config->bufferRecieveMeshLeft = (FLT *)calloc(13*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
		pnt_config->bufferRecieveMeshRight = (FLT *)calloc(13*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
		pnt_config->bufferRecieveMeshBottom = (FLT *)calloc(13*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
		pnt_config->bufferRecieveMeshTop = (FLT *)calloc(13*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));

		if(MESHDIMENSIONS==2)
		{
			pnt_config->bufferSendFlowWithGhostsLeft = (FLT *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferSendFlowWithGhostsRight= (FLT *)calloc(4*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferSendFlowWithGhostsBottom = (FLT *)calloc(4*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferSendFlowWithGhostsTop= (FLT *)calloc(4*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));

			pnt_config->bufferRecieveFlowWithGhostsLeft = (FLT *)calloc(4*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferRecieveFlowWithGhostsRight= (FLT *)calloc(4*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferRecieveFlowWithGhostsBottom = (FLT *)calloc(4*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferRecieveFlowWithGhostsTop= (FLT *)calloc(4*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
		}
		if(MESHDIMENSIONS==3)
		{
			pnt_config->bufferSendFlowWithGhostsLeft = (FLT *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferSendFlowWithGhostsRight= (FLT *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferSendFlowWithGhostsBottom = (FLT *)calloc(5*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferSendFlowWithGhostsTop= (FLT *)calloc(5*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));

			pnt_config->bufferRecieveFlowWithGhostsLeft = (FLT *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferRecieveFlowWithGhostsRight= (FLT *)calloc(5*((pnt_config->int_SpaceOrder+1)/2) * pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferRecieveFlowWithGhostsBottom = (FLT *)calloc(5*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));
			pnt_config->bufferRecieveFlowWithGhostsTop= (FLT *)calloc(5*pnt_config->int_iMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2)*pnt_config->int_kMeshPointsGhostCells,sizeof(FLT));

			pnt_config->bufferSendFlowWithGhostsBehind = (FLT *)calloc(5*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(FLT));
			pnt_config->bufferSendFlowWithGhostsInFront= (FLT *)calloc(5*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(FLT));
			pnt_config->bufferRecieveFlowWithGhostsBehind = (FLT *)calloc(5*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(FLT));
			pnt_config->bufferRecieveFlowWithGhostsInFront= (FLT *)calloc(5*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(FLT));

			pnt_config->bufferSendMeshBehind = (FLT *)calloc(13*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(FLT));
			pnt_config->bufferSendMeshInFront = (FLT *)calloc(13*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(FLT));
			pnt_config->bufferRecieveMeshBehind = (FLT *)calloc(13*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(FLT));
			pnt_config->bufferRecieveMeshInFront= (FLT *)calloc(13*pnt_config->int_iMeshPointsGhostCells * pnt_config->int_jMeshPointsGhostCells * ((pnt_config->int_SpaceOrder+1)/2),sizeof(FLT));
		}

		/* get neighbour */
//		pnt_config->MPI_rankNeighbours= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_tag= (long *)calloc(pnt_config->NumberInterfaces, sizeof(long));
//		getNeighbour(
//				pnt_config);


		int MPI_tag_factor[6];
		MPI_tag_factor[0]=1;
		MPI_tag_factor[1]=3;
		MPI_tag_factor[2]=5;
		MPI_tag_factor[3]=7;
		MPI_tag_factor[4]=11;
		MPI_tag_factor[5]=13;

		for(i=0;i<pnt_config->NumberInterfaces;i++)
		{
			int flag,l;
			void *v;
			int value_MPI_TAG_UB;
			MPI_Comm_get_attr( MPI_COMM_WORLD, MPI_TAG_UB, &v, &flag );
			value_MPI_TAG_UB = *(int*)v;
	//		printf("MAX TAG: %d\n",value_MPI_TAG_UB);

			if(pnt_config->MPI_rank<=pnt_config->MPI_rankNeighbours[i])
			{
				pnt_config->MPI_tag[i]=0;
				for(k=0;k<2*pnt_config->int_meshDimensions;k++)
		        {
					pnt_config->MPI_tag[i]+=(pnt_config->RangeOfInterface[i][k]*MPI_tag_factor[k]);
		        }
				while(pnt_config->MPI_tag[i]>value_MPI_TAG_UB)
				{
					pnt_config->MPI_tag[i]=pnt_config->MPI_tag[i]-value_MPI_TAG_UB;
				}
			}
			else
			{
				pnt_config->MPI_tag[i]=0;
				for(k=0;k<2*pnt_config->int_meshDimensions;k++)
		        {
					pnt_config->MPI_tag[i]+=(pnt_config->DonorRangeOfInterface[i][k]*MPI_tag_factor[k]);
		        }
				while(pnt_config->MPI_tag[i]>value_MPI_TAG_UB)
				{
					pnt_config->MPI_tag[i]=pnt_config->MPI_tag[i]-value_MPI_TAG_UB;
				}
			}
			for(l=0;l<i;l++)
			{
				if ((pnt_config->MPI_tag[i]==pnt_config->MPI_tag[l])&&(pnt_config->MPI_rankNeighbours[i]==pnt_config->MPI_rankNeighbours[l]))
				{
					printf("ERROR: Zwei Interfaces mit identischer tag entdeckt:");
					printf("tag1/2: rank %d - partner %d - TAG: %ld\n",pnt_config->MPI_rank,pnt_config->MPI_rankNeighbours[l],pnt_config->MPI_tag[l]);
				}
			}
		}

		pnt_config->MPI_intTransformation_IMax= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_JMax= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_KMax= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_IMax_Mesh= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_JMax_Mesh= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_KMax_Mesh= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_flag_I0_I= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_flag_I0_J= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_flag_I0_K= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_flag_J0_I= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_flag_J0_J= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_flag_J0_K= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_flag_K0_I= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_flag_K0_J= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_flag_K0_K= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_Offset_I= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_Offset_J= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_Offset_K= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_Offset_I_Ghosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_Offset_J_Ghosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransformation_Offset_K_Ghosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));

		pnt_config->MPI_dblTransformation_xi_x= (FLT *)calloc(pnt_config->NumberInterfaces, sizeof(FLT));
		pnt_config->MPI_dblTransformation_xi_y= (FLT *)calloc(pnt_config->NumberInterfaces, sizeof(FLT));
		pnt_config->MPI_dblTransformation_xi_z= (FLT *)calloc(pnt_config->NumberInterfaces, sizeof(FLT));
		pnt_config->MPI_dblTransformation_eta_x= (FLT *)calloc(pnt_config->NumberInterfaces, sizeof(FLT));
		pnt_config->MPI_dblTransformation_eta_y= (FLT *)calloc(pnt_config->NumberInterfaces, sizeof(FLT));
		pnt_config->MPI_dblTransformation_eta_z= (FLT *)calloc(pnt_config->NumberInterfaces, sizeof(FLT));
		pnt_config->MPI_dblTransformation_zeta_x= (FLT *)calloc(pnt_config->NumberInterfaces, sizeof(FLT));
		pnt_config->MPI_dblTransformation_zeta_y= (FLT *)calloc(pnt_config->NumberInterfaces, sizeof(FLT));
		pnt_config->MPI_dblTransformation_zeta_z= (FLT *)calloc(pnt_config->NumberInterfaces, sizeof(FLT));

		//		pnt_config->MPI_charNeighbours= (char **)calloc(pnt_config->NumberInterfaces, sizeof(char*));
		for(i=0;i<pnt_config->NumberInterfaces;i++)
		{
		//			pnt_config->RankNeighbours[i]= (int *)calloc(pnt_config->int_meshDimensions+1, sizeof(int));
		//			pnt_config->MPI_charNeighbours[i]= (char *)calloc(30, sizeof(char));

			pnt_config->MPI_dblTransformation_xi_x[i]=1.0;
			pnt_config->MPI_dblTransformation_xi_y[i]=1.0;
			pnt_config->MPI_dblTransformation_xi_z[i]=1.0;
			pnt_config->MPI_dblTransformation_eta_x[i]=1.0;
			pnt_config->MPI_dblTransformation_eta_y[i]=1.0;
			pnt_config->MPI_dblTransformation_eta_z[i]=1.0;
			pnt_config->MPI_dblTransformation_zeta_x[i]=1.0;
			pnt_config->MPI_dblTransformation_zeta_y[i]=1.0;
			pnt_config->MPI_dblTransformation_zeta_z[i]=1.0;
		}

		pnt_config->MPI_intTransferSizeMesh=(int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intTransferSizeFlow_WithGhosts=(int *)calloc(pnt_config->NumberInterfaces, sizeof(int));

		pnt_config->MPI_SendBufferMesh= (FLT **)calloc(pnt_config->NumberInterfaces, sizeof(FLT *));
		pnt_config->MPI_RecieveBufferMesh= (FLT **)calloc(pnt_config->NumberInterfaces, sizeof(FLT *));
		pnt_config->MPI_SendBufferFlowWithGhosts= (FLT **)calloc(pnt_config->NumberInterfaces, sizeof(FLT *));
		pnt_config->MPI_RecieveBufferFlowWithGhosts= (FLT **)calloc(pnt_config->NumberInterfaces, sizeof(FLT *));

		pnt_config->MPI_intIStartSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intIEndSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intJStartSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intJEndSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intKStartSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intKEndSend= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));

		pnt_config->MPI_intIStartRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intIEndRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intJStartRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intJEndRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intKStartRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intKEndRecieve= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));

		pnt_config->MPI_intIStartSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intIEndSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intJStartSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intJEndSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intKStartSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intKEndSend_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));

		pnt_config->MPI_intIStartRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intIEndRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intJStartRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intJEndRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intKStartRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));
		pnt_config->MPI_intKEndRecieve_WithGhosts= (int *)calloc(pnt_config->NumberInterfaces, sizeof(int));



		//		Teile den Interfaces, an denen der einzelne Rank beteiligt ist, die Informationen zu,
		//		ob dieses Interface aus seiner Sicht links, rechts... ist.
		for(i=0;i<pnt_config->NumberInterfaces;i++)
		{
				if(i==pnt_config->InterfaceNeighbourLeft)
				{
	//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshI;
	//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshI;
	//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowI;
	//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowI;
	//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidI;
	//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidI;

					pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshLeft;
					pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshLeft;
					pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsLeft;
					pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsLeft;

				   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iStartReal;
				   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iStartReal+((pnt_config->int_SpaceOrder+1)/2-1);
				   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jStartReal;
				   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jEndReal;
				   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kStartReal;
				   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kEndReal;

				   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iStartGhosts;
				   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iStartReal-1;
				   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jStartReal;
				   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jEndReal;
				   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kStartReal;
				   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kEndReal;


				   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iStartReal;
				   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iStartReal+((pnt_config->int_SpaceOrder+1)/2-1);
				   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jStartGhosts;
				   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jEndGhosts;
				   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kStartGhosts;
				   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kEndGhosts;

				   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iStartGhosts;
				   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iStartReal-1;
				   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jStartGhosts;
				   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jEndGhosts;
				   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kStartGhosts;
				   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kEndGhosts;
				   if(MESHDIMENSIONS==2)
				   {
					   pnt_config->MPI_intKStartSend_WithGhosts[i]=0;
					   pnt_config->MPI_intKEndSend_WithGhosts[i]=0;
					   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=0;
					   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=0;

					   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
							   (pnt_config->int_SpaceOrder+1)/2*4*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
				   }
				   else
				   {
					   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
						   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
				   }
				   pnt_config->MPI_intTransferSizeMesh[i]=
						   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;



				}
					if(i==pnt_config->InterfaceNeighbourRight)
				{
	//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshI;
	//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshI;
	//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowI;
	//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowI;
	//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidI;
	//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidI;

					pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshRight;
					pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshRight;
					pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsRight;
					pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsRight;

				   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
				   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iEndReal;
				   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jStartReal;
				   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jEndReal;
				   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kStartReal;
				   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kEndReal;

				   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iEndReal+1;
				   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iEndGhosts;
				   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jStartReal;
				   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jEndReal;
				   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kStartReal;
				   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kEndReal;


				   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
				   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iEndReal;
				   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jStartGhosts;
				   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jEndGhosts;
				   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kStartGhosts;
				   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kEndGhosts;

				   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iEndReal+1;
				   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iEndGhosts;
				   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jStartGhosts;
				   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jEndGhosts;
				   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kStartGhosts;
				   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kEndGhosts;
				   if(MESHDIMENSIONS==2)
				   {
					   pnt_config->MPI_intKStartSend_WithGhosts[i]=0;
					   pnt_config->MPI_intKEndSend_WithGhosts[i]=0;
					   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=0;
					   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=0;
					   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
					   					   (pnt_config->int_SpaceOrder+1)/2*4*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
				   }
				   else
				   {
					   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
					   					   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
				   }

				   pnt_config->MPI_intTransferSizeMesh[i]=
						   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;

				}
					if(i==pnt_config->InterfaceNeighbourBottom)
				{
	//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshJ;
	//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshJ;
	//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowJ;
	//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowJ;
	//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidJ;
	//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidJ;

					pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshBottom;
					pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshBottom;
					pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsBottom;
					pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsBottom;

				   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jStartReal;
				   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jStartReal+((pnt_config->int_SpaceOrder+1)/2-1);
				   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iStartReal;
				   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iEndReal;
				   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kStartReal;
				   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kEndReal;

				   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jStartGhosts;
				   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jStartReal-1;
				   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iStartReal;
				   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iEndReal;
				   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kStartReal;
				   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kEndReal;


				   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iStartGhosts;
				   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iEndGhosts;
				   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jStartReal;
				   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jStartReal+((pnt_config->int_SpaceOrder+1)/2-1);
				   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kStartGhosts;
				   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kEndGhosts;

				   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iStartGhosts;
				   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iEndGhosts;
				   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jStartGhosts;
				   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jStartReal-1;
				   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kStartGhosts;
				   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kEndGhosts;
				   if(MESHDIMENSIONS==2)
				   {
					   pnt_config->MPI_intKStartSend_WithGhosts[i]=0;
					   pnt_config->MPI_intKEndSend_WithGhosts[i]=0;
					   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=0;
					   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=0;

					   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
							   (pnt_config->int_SpaceOrder+1)/2*4*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
				   }
				   else
				   {
					   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
							   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
				   }

				   pnt_config->MPI_intTransferSizeMesh[i]=
						   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
				}
					if(i==pnt_config->InterfaceNeighbourTop)
				{
	//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshJ;
	//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshJ;
	//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowJ;
	//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowJ;
	//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidJ;
	//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidJ;

					pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshTop;
					pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshTop;
					pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsTop;
					pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsTop;

				   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
				   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jEndReal;
				   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iStartReal;
				   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iEndReal;
				   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kStartReal;
				   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kEndReal;

				   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jEndReal+1;
				   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jEndGhosts;
				   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iStartReal;
				   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iEndReal;
				   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kStartReal;
				   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kEndReal;


				   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iStartGhosts;
				   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iEndGhosts;
				   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
				   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jEndReal;
				   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kStartGhosts;
				   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kEndGhosts;

				   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iStartGhosts;
				   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iEndGhosts;
				   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jEndReal+1;
				   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jEndGhosts;
				   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kStartGhosts;
				   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kEndGhosts;
				   if(MESHDIMENSIONS==2)
				   {
					   pnt_config->MPI_intKStartSend_WithGhosts[i]=0;
					   pnt_config->MPI_intKEndSend_WithGhosts[i]=0;
					   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=0;
					   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=0;

					   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
							   (pnt_config->int_SpaceOrder+1)/2*4*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
				   }
				   else
				   {
					   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
							   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;
				   }

				   pnt_config->MPI_intTransferSizeMesh[i]=
						   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;

				}
					if(i==pnt_config->InterfaceNeighbourBehind)
				{
	//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshK;
	//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshK;
	//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowK;
	//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowK;
	//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidK;
	//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidK;

					pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshBehind;
					pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshBehind;
					pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsBehind;
					pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsBehind;

				   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kStartReal;
				   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kStartReal+((pnt_config->int_SpaceOrder+1)/2-1);
				   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iStartReal;
				   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iEndReal;
				   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jStartReal;
				   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jEndReal;

				   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kStartGhosts;
				   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kStartReal-1;
				   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iStartReal;
				   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iEndReal;
				   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jStartReal;
				   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jEndReal;


				   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iStartGhosts;
				   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iEndGhosts;
				   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jStartGhosts;
				   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jEndGhosts;
				   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kStartReal;
				   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kStartReal+((pnt_config->int_SpaceOrder+1)/2-1);

				   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iStartGhosts;
				   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iEndGhosts;
				   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jStartGhosts;
				   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jEndGhosts;
				   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kStartGhosts;
				   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kStartReal-1;

				   pnt_config->MPI_intTransferSizeMesh[i]=
						   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_iMeshPointsGhostCells;
				   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
						   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_iMeshPointsGhostCells;

				}
					if(i==pnt_config->InterfaceNeighbourInFront)
				{
	//				pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshK;
	//				pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshK;
	//				pnt_config->MPI_SendBufferFlow[i]=pnt_config->bufferSendFlowK;
	//				pnt_config->MPI_RecieveBufferFlow[i]=pnt_config->bufferRecieveFlowK;
	//				pnt_config->MPI_SendBufferViscid[i]=pnt_config->bufferSendViscidK;
	//				pnt_config->MPI_RecieveBufferViscid[i]=pnt_config->bufferRecieveViscidK;

					pnt_config->MPI_SendBufferMesh[i]=pnt_config->bufferSendMeshInFront;
					pnt_config->MPI_RecieveBufferMesh[i]=pnt_config->bufferRecieveMeshInFront;
					pnt_config->MPI_SendBufferFlowWithGhosts[i]=pnt_config->bufferSendFlowWithGhostsInFront;
					pnt_config->MPI_RecieveBufferFlowWithGhosts[i]=pnt_config->bufferRecieveFlowWithGhostsInFront;

				   pnt_config->MPI_intKStartSend[i]=pnt_config->int_kEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
				   pnt_config->MPI_intKEndSend[i]=pnt_config->int_kEndReal;
				   pnt_config->MPI_intIStartSend[i]=pnt_config->int_iStartReal;
				   pnt_config->MPI_intIEndSend[i]=pnt_config->int_iEndReal;
				   pnt_config->MPI_intJStartSend[i]=pnt_config->int_jStartReal;
				   pnt_config->MPI_intJEndSend[i]=pnt_config->int_jEndReal;

				   pnt_config->MPI_intKStartRecieve[i]=pnt_config->int_kEndReal+1;
				   pnt_config->MPI_intKEndRecieve[i]=pnt_config->int_kEndGhosts;
				   pnt_config->MPI_intIStartRecieve[i]=pnt_config->int_iStartReal;
				   pnt_config->MPI_intIEndRecieve[i]=pnt_config->int_iEndReal;
				   pnt_config->MPI_intJStartRecieve[i]=pnt_config->int_jStartReal;
				   pnt_config->MPI_intJEndRecieve[i]=pnt_config->int_jEndReal;


				   pnt_config->MPI_intIStartSend_WithGhosts[i]=pnt_config->int_iStartGhosts;
				   pnt_config->MPI_intIEndSend_WithGhosts[i]=pnt_config->int_iEndGhosts;
				   pnt_config->MPI_intJStartSend_WithGhosts[i]=pnt_config->int_jStartGhosts;
				   pnt_config->MPI_intJEndSend_WithGhosts[i]=pnt_config->int_jEndGhosts;
				   pnt_config->MPI_intKStartSend_WithGhosts[i]=pnt_config->int_kEndReal-((pnt_config->int_SpaceOrder+1)/2-1);
				   pnt_config->MPI_intKEndSend_WithGhosts[i]=pnt_config->int_kEndReal;

				   pnt_config->MPI_intIStartRecieve_WithGhosts[i]=pnt_config->int_iStartGhosts;
				   pnt_config->MPI_intIEndRecieve_WithGhosts[i]=pnt_config->int_iEndGhosts;
				   pnt_config->MPI_intJStartRecieve_WithGhosts[i]=pnt_config->int_jStartGhosts;
				   pnt_config->MPI_intJEndRecieve_WithGhosts[i]=pnt_config->int_jEndGhosts;
				   pnt_config->MPI_intKStartRecieve_WithGhosts[i]=pnt_config->int_kEndReal+1;
				   pnt_config->MPI_intKEndRecieve_WithGhosts[i]=pnt_config->int_kEndGhosts;

				   pnt_config->MPI_intTransferSizeMesh[i]=
						   (pnt_config->int_SpaceOrder+1)/2*13*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_iMeshPointsGhostCells;
				   pnt_config->MPI_intTransferSizeFlow_WithGhosts[i]=
						   (pnt_config->int_SpaceOrder+1)/2*5*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_iMeshPointsGhostCells;

				}
					check_TransformationMatrix(
						i,
						pnt_config);
		}

//	printf("SHOCK: Rank: %d (%s)",pnt_config->MPI_rank,pnt_config->Zonename);
//	printf(": Left: %d, index: %d",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourLeft],pnt_config->InterfaceNeighbourLeft);
//	printf(": Right: %d, index: %d ",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourRight],pnt_config->InterfaceNeighbourRight);
//	printf(": Bottom: %d, index: %d",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourBottom],pnt_config->InterfaceNeighbourBottom);
//	printf(": Top: %d, index: %d",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourTop],pnt_config->InterfaceNeighbourTop);
//	printf(": Behind: %d, index: %d ",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourBehind],pnt_config->InterfaceNeighbourBehind);
//	printf(": InFront: %d, index: %d\n",pnt_config->MPI_rankNeighbours[pnt_config->InterfaceNeighbourInFront],pnt_config->InterfaceNeighbourInFront);
//
//	for (i=0;i<pnt_config->NumberInterfaces;i++)
//	{
//		printf("myrank: %d | rank: %d | interface:%d | transform: %d %d %d | donorname: %s\n",
//				pnt_config->MPI_rank,
//				 pnt_config->MPI_rankNeighbours[i],i,
//				pnt_config->TransformMatrixOfInterface[i][0],
//				pnt_config->TransformMatrixOfInterface[i][1],
//				pnt_config->TransformMatrixOfInterface[i][2],
//				pnt_config->Donorname[i]);
//
//	}

	AllocMemory(
			pnt_config,
			pnt_mesh,
			pnt_U_lastStep,
			pnt_U_RK,
			pnt_Flux,
			pnt_Flux_PlusHalf,
			pnt_Q,
			pnt_Q_sum,
			pnt_Film,
			pnt_U_backup1,
			pnt_U_backup2);

	if(pnt_config->MPI_rank==0){printf("SHOCK: Speicherallokierung fertig!\n");}

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				ijk=k*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells+j*pnt_config->int_iMeshPointsGhostCells+i;

				i2=i-pnt_config->int_iStartReal;
				j2=j-pnt_config->int_jStartReal;
				k2=k-pnt_config->int_kStartReal;


				ijk2=k2*pnt_config->int_jMeshPoints*pnt_config->int_iMeshPoints+j2*pnt_config->int_iMeshPoints+i2;

				if( x->dt==RealSingle )
					pnt_mesh->x[ijk]=(FLT)( (float*)x->ptr )[ijk2];
				else
					pnt_mesh->x[ijk]=(FLT)( (double*)x->ptr )[ijk2];

				if( y->dt==RealSingle )
					pnt_mesh->y[ijk]=(FLT)( (float*)y->ptr )[ijk2];
				else
					pnt_mesh->y[ijk]=(FLT)( (double*)y->ptr )[ijk2];

				if(MESHDIMENSIONS==3)
				{
					if( z->dt==RealSingle )
						pnt_mesh->z[ijk]=(FLT)( (float*)z->ptr )[ijk2];
					else
						pnt_mesh->z[ijk]=(FLT)( (double*)z->ptr )[ijk2];
				}
				else
				{
					pnt_mesh->z[ijk]=0.0;
				}
			}
		}
	}

	Initialize(
			pnt_config,
			pnt_mesh,
			pnt_U_lastStep);

	if (pnt_config->int_initializeType==0)
	{
		if(pnt_config->MPI_rank==0){printf("SHOCK: Neu-Initialisierung des Rechengebietes fertig!\n");}
	}

	if(
			(u->ptr!=NULL)&&
			(abs(pnt_config->int_initializeType)==1)
			)
	{
		for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
		{
			for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
			{
				for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
	//				ijk=k*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells+j*pnt_config->int_iMeshPointsGhostCells+i;

					i2=i-pnt_config->int_iStartReal;
					j2=j-pnt_config->int_jStartReal;
					k2=k-pnt_config->int_kStartReal;

					ijk2=k2*pnt_config->int_jMeshPoints*pnt_config->int_iMeshPoints+j2*pnt_config->int_iMeshPoints+i2;

					if( u->dt==RealSingle )
						pnt_U_lastStep->u[ijk]=(FLT)( (float*)u->ptr )[ijk2];
					else
						pnt_U_lastStep->u[ijk]=(FLT)( (double*)u->ptr )[ijk2];

					if( v->dt==RealSingle )
						pnt_U_lastStep->v[ijk]=(FLT)( (float*)v->ptr )[ijk2];
					else
						pnt_U_lastStep->v[ijk]=(FLT)( (double*)v->ptr )[ijk2];


					FLT total_z_size=0.1;
					if( w->dt==RealSingle )
					{
						//				Creating 3D-Waves in order to obtain 3D-flow
						if((strcmp(pnt_config->BC_Bottom,pnt_config->BCWallViscous)==0)&&(pnt_config->int_initializeType==-1))
						{
							//printf("Creating 3D-waves...\n");
							pnt_U_lastStep->w[ijk]=0.01*sin(pnt_mesh->z[ijk]/total_z_size*3.0*2.0*MY_PI);
						}
						else
						{
							pnt_U_lastStep->w[ijk]=( (float*)w->ptr )[ijk2];
						}

					}
					else
					{
						//				Creating 3D-Waves in order to obtain 3D-flow
						if((strcmp(pnt_config->BC_Bottom,pnt_config->BCWallViscous)==0)&&(pnt_config->int_initializeType==-1))
						{
							//printf("Creating 3D-waves...\n");
							pnt_U_lastStep->w[ijk]=0.01*sin(pnt_mesh->z[ijk]/total_z_size*3.0*2.0*MY_PI);
						}
						else
						{
							pnt_U_lastStep->w[ijk]=(FLT)( (double*)w->ptr )[ijk2];
						}
					}

					if( p->dt==RealSingle )
						pnt_U_lastStep->p[ijk]=(FLT)( (float*)p->ptr )[ijk2];
					else
						pnt_U_lastStep->p[ijk]=(FLT)( (double*)p->ptr )[ijk2];

					if( rho->dt==RealSingle )
						pnt_U_lastStep->rho[ijk]=(FLT)( (float*)rho->ptr )[ijk2];
					else
						pnt_U_lastStep->rho[ijk]=(FLT)( (double*)rho->ptr )[ijk2];


					pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+
											pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);

				}
			}
		}
		if(pnt_config->MPI_rank==0){printf("SHOCK: Import der Ergebnisse fertig! (Iteration: %d, Time: %g)\n",
				pnt_config->int_StartIteration,
				(double)pnt_config->start_Time);}
	}

	if(pnt_config->int_initializeType==1)
	{

//			IBC
		if((pnt_config->flag_IBC==1)&&(pnt_config->flag_IBC_Moving==1))
		{
			pnt_config->IBC_MovingLastPosition=IBC_getActualPosition(pnt_config);
			pnt_config->IBC_MovingActualPosition=IBC_getActualPosition(pnt_config);

			IBC_Actual2Last(
					pnt_config,
					pnt_mesh);

			IBC_prepare(
					pnt_config,
					pnt_mesh);
		}

	}
}

void AllocMemory(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf,
		struct strct_Flux * pnt_Q,
		struct strct_Flux * pnt_Q_sum,
		struct strct_Film * pnt_Film,
		struct strct_U * pnt_U_backup1,
		struct strct_U * pnt_U_backup2)
{
	AllocMemoryMesh(
			pnt_config,
			pnt_mesh);
	AllocMemoryStrctU(
			pnt_config,
			pnt_U_lastStep);
	AllocMemoryStrctU(
			pnt_config,
			pnt_U_RK);
	AllocMemoryStrctFlux(
			pnt_config,
			pnt_Flux);
	AllocMemoryStrctFlux(
			pnt_config,
			pnt_Flux_PlusHalf);
	AllocMemoryStrctFlux(
			pnt_config,
			pnt_Q);
	AllocMemoryStrctFlux(
			pnt_config,
			pnt_Q_sum);
	AllocMemoryStrctFilm(
			pnt_config,
			pnt_Film);

	AllocMemoryBackup(
			pnt_config,
			pnt_U_backup1,
			pnt_U_backup2);
}

void AllocMemoryBackup(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U_backup1,
		struct strct_U * pnt_U_backup2)
{
	//	Backup1
	pnt_U_backup1->rho = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_U_backup1->u = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_U_backup1->v = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_U_backup1->w = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_U_backup1->p = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_U_backup1->e = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));

	//	Backup2
	pnt_U_backup2->rho = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_U_backup2->u = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_U_backup2->v = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_U_backup2->w = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_U_backup2->p = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_U_backup2->e = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
}

void FreeMemory(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf,
		struct strct_Flux * pnt_Q,
		struct strct_Flux * pnt_Q_sum,
		struct strct_Film * pnt_Film)
{
	FreeMemoryMesh(
			pnt_mesh,
			pnt_config);
	FreeMemoryStrctU(
			pnt_U_lastStep);
	FreeMemoryStrctU(
			pnt_U_RK);
	FreeMemoryStrctFlux(
			pnt_Flux);
	FreeMemoryStrctFlux(
			pnt_Flux_PlusHalf);
	FreeMemoryStrctFlux(
			pnt_Q);
	FreeMemoryStrctFlux(
			pnt_Q_sum);
	FreeMemoryStrctFilm(
			pnt_Film);
	
/*
	for(i=0;i<pnt_config->NumberInterfaces;i++)
	{
		free(pnt_config->Donorname[i]);
		free(pnt_config->Interfacename[i]);
		free(pnt_config->TransformMatrixOfInterface[i]);
		free(pnt_config->RangeOfInterface[i]);
		free(pnt_config->DonorRangeOfInterface[i]);
		free(pnt_config->RotationCenter[i]);
		free(pnt_config->RotationAngle[i]);
		free(pnt_config->Translation[i]);
	}
	free(pnt_config->Donorname);
	free(pnt_config->Interfacename);
	free(pnt_config->TransformMatrixOfInterface);
	free(pnt_config->RangeOfInterface);
	free(pnt_config->DonorRangeOfInterface);
	free(pnt_config->RotationCenter);
	free(pnt_config->RotationAngle);
	free(pnt_config->Translation);
*/

	free(pnt_config->bufferSendFlowWithGhostsLeft);
	free(pnt_config->bufferSendFlowWithGhostsRight);
	free(pnt_config->bufferSendFlowWithGhostsBottom);
	free(pnt_config->bufferSendFlowWithGhostsTop);
	free(pnt_config->bufferRecieveFlowWithGhostsLeft);
	free(pnt_config->bufferRecieveFlowWithGhostsRight);
	free(pnt_config->bufferRecieveFlowWithGhostsBottom);
	free(pnt_config->bufferRecieveFlowWithGhostsTop);

	free(pnt_config->bufferSendMeshLeft);
	free(pnt_config->bufferSendMeshRight);
	free(pnt_config->bufferSendMeshBottom);
	free(pnt_config->bufferSendMeshTop);
	free(pnt_config->bufferRecieveMeshLeft);
	free(pnt_config->bufferRecieveMeshRight);
	free(pnt_config->bufferRecieveMeshBottom);
	free(pnt_config->bufferRecieveMeshTop);

	if(MESHDIMENSIONS==3)
	{
		free(pnt_config->bufferSendFlowWithGhostsBehind);
		free(pnt_config->bufferSendFlowWithGhostsInFront);
		free(pnt_config->bufferRecieveFlowWithGhostsBehind);
		free(pnt_config->bufferRecieveFlowWithGhostsInFront);

		free(pnt_config->bufferSendMeshBehind);
		free(pnt_config->bufferSendMeshInFront);
		free(pnt_config->bufferRecieveMeshBehind);
		free(pnt_config->bufferRecieveMeshInFront);
	}
	free(pnt_config->MPI_intTransferSizeMesh);
	free(pnt_config->MPI_intTransferSizeFlow_WithGhosts);

	free(pnt_config->MPI_SendBufferMesh);
	free(pnt_config->MPI_RecieveBufferMesh);
	free(pnt_config->MPI_SendBufferFlowWithGhosts);
	free(pnt_config->MPI_RecieveBufferFlowWithGhosts);

	free(pnt_config->MPI_intIStartSend);
	free(pnt_config->MPI_intIEndSend);
	free(pnt_config->MPI_intJStartSend);
	free(pnt_config->MPI_intJEndSend);
	free(pnt_config->MPI_intKStartSend);
	free(pnt_config->MPI_intKEndSend);

	free(pnt_config->MPI_intIStartRecieve);
	free(pnt_config->MPI_intIEndRecieve);
	free(pnt_config->MPI_intJStartRecieve);
	free(pnt_config->MPI_intJEndRecieve);
	free(pnt_config->MPI_intKStartRecieve);
	free(pnt_config->MPI_intKEndRecieve);

	free(pnt_config->MPI_intIStartSend_WithGhosts);
	free(pnt_config->MPI_intIEndSend_WithGhosts);
	free(pnt_config->MPI_intJStartSend_WithGhosts);
	free(pnt_config->MPI_intJEndSend_WithGhosts);
	free(pnt_config->MPI_intKStartSend_WithGhosts);
	free(pnt_config->MPI_intKEndSend_WithGhosts);

	free(pnt_config->MPI_intIStartRecieve_WithGhosts);
	free(pnt_config->MPI_intIEndRecieve_WithGhosts);
	free(pnt_config->MPI_intJStartRecieve_WithGhosts);
	free(pnt_config->MPI_intJEndRecieve_WithGhosts);
	free(pnt_config->MPI_intKStartRecieve_WithGhosts);
	free(pnt_config->MPI_intKEndRecieve_WithGhosts);

	free(pnt_config->MPI_rankNeighbours);

	free(pnt_config->MPI_dblTransformation_xi_x);
	free(pnt_config->MPI_dblTransformation_xi_y);
	free(pnt_config->MPI_dblTransformation_xi_z);
	free(pnt_config->MPI_dblTransformation_eta_x);
	free(pnt_config->MPI_dblTransformation_eta_y);
	free(pnt_config->MPI_dblTransformation_eta_z);
	free(pnt_config->MPI_dblTransformation_zeta_x);
	free(pnt_config->MPI_dblTransformation_zeta_y);
	free(pnt_config->MPI_dblTransformation_zeta_z);

	free(pnt_config->MPI_intTransformation_IMax);
	free(pnt_config->MPI_intTransformation_JMax);
	free(pnt_config->MPI_intTransformation_KMax);
	free(pnt_config->MPI_intTransformation_IMax_Mesh);
	free(pnt_config->MPI_intTransformation_JMax_Mesh);
	free(pnt_config->MPI_intTransformation_KMax_Mesh);
	free(pnt_config->MPI_intTransformation_flag_I0_I);
	free(pnt_config->MPI_intTransformation_flag_I0_J);
	free(pnt_config->MPI_intTransformation_flag_I0_K);
	free(pnt_config->MPI_intTransformation_flag_J0_I);
	free(pnt_config->MPI_intTransformation_flag_J0_J);
	free(pnt_config->MPI_intTransformation_flag_J0_K);
	free(pnt_config->MPI_intTransformation_flag_K0_I);
	free(pnt_config->MPI_intTransformation_flag_K0_J);
	free(pnt_config->MPI_intTransformation_flag_K0_K);
	free(pnt_config->MPI_intTransformation_Offset_I);
	free(pnt_config->MPI_intTransformation_Offset_J);
	free(pnt_config->MPI_intTransformation_Offset_K);
	free(pnt_config->MPI_intTransformation_Offset_I_Ghosts);
	free(pnt_config->MPI_intTransformation_Offset_J_Ghosts);
	free(pnt_config->MPI_intTransformation_Offset_K_Ghosts);


}

void AllocMemoryMesh(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{
	//Speicherallokierung für die einzelnen Strukturen
	// 1)Gitter
	// 2)Erhaltungsgrößen lastStep
	// 3)Erhaltungsgrößen nextStep
	pnt_mesh->BC_Corrector_xiMomentum = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->BC_Corrector_etaMomentum = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->BC_Corrector_zetaMomentum = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));

	//Gitter
	pnt_mesh->x = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->y = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->z = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));

	pnt_mesh->x_extrapolate = (FLT *)calloc((pnt_config->int_iMeshPointsGhostCells+2)*(pnt_config->int_jMeshPointsGhostCells+2)*(pnt_config->int_kMeshPointsGhostCells+2), sizeof(FLT ));
	pnt_mesh->y_extrapolate = (FLT *)calloc((pnt_config->int_iMeshPointsGhostCells+2)*(pnt_config->int_jMeshPointsGhostCells+2)*(pnt_config->int_kMeshPointsGhostCells+2), sizeof(FLT ));
	pnt_mesh->z_extrapolate = (FLT *)calloc((pnt_config->int_iMeshPointsGhostCells+2)*(pnt_config->int_jMeshPointsGhostCells+2)*(pnt_config->int_kMeshPointsGhostCells+2), sizeof(FLT ));

	pnt_mesh->xi_x = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->xi_y = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->xi_z = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));

	pnt_mesh->eta_x = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->eta_y = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->eta_z = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));

	pnt_mesh->zeta_x = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->zeta_y = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->zeta_z = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));

	pnt_mesh->jacobian = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));

	//	Immerged BC
	pnt_mesh->flag_IBC = (int *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(int ));
	pnt_mesh->flag_IBC_last = (int *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(int ));

	//	Pressure Waves
	if(pnt_config->flag_PressureWaves==1)
	{
		pnt_mesh->flag_PressureWaves = (int *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(int ));
		pnt_mesh->startPressure_PressureWaves = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
		pnt_mesh->startDensity_PressureWaves = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	}

	pnt_mesh->xiFluss_Faktor = (FLT *)calloc(30*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->etaFluss_Faktor = (FLT *)calloc(30*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_mesh->zetaFluss_Faktor = (FLT *)calloc(30*pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
}

void AllocMemoryStrctU(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_strctU)
{
	pnt_strctU->rho = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->u = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->v = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->w = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->p = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->e = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->theta1 = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->theta2 = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->theta3 = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->c = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->gradRho = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->T = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->mue = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->Lambda2 = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->MachNumber = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));

	pnt_strctU->u_xi = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->u_eta = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->u_zeta = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->v_xi = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->v_eta = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->v_zeta = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->w_xi = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->w_eta = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->w_zeta = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->T_xi = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->T_eta = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctU->T_zeta = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
}

void FreeMemoryStrctU(
		struct strct_U * pnt_strctU)
{
	free(pnt_strctU->rho);
	free(pnt_strctU->u);
	free(pnt_strctU->v);
	free(pnt_strctU->w);
	free(pnt_strctU->p);
	free(pnt_strctU->e);
	free(pnt_strctU->theta1);
	free(pnt_strctU->theta2);
	free(pnt_strctU->theta3);
	free(pnt_strctU->c);
	free(pnt_strctU->gradRho);
	free(pnt_strctU->T);
	free(pnt_strctU->mue);
	free(pnt_strctU->Lambda2);
	free(pnt_strctU->MachNumber);
	free(pnt_strctU->u_xi);
	free(pnt_strctU->u_eta);
	free(pnt_strctU->u_zeta);
	free(pnt_strctU->v_xi);
	free(pnt_strctU->v_eta);
	free(pnt_strctU->v_zeta);
	free(pnt_strctU->w_xi);
	free(pnt_strctU->w_eta);
	free(pnt_strctU->w_zeta);
	free(pnt_strctU->T_xi);
	free(pnt_strctU->T_eta);
	free(pnt_strctU->T_zeta);
}

void AllocMemoryStrctFilm(
		struct strct_configuration * pnt_config,
		struct strct_Film * pnt_strctFilm)
{
	pnt_strctFilm->rho = (FLT *)calloc(pnt_config->int_Samples*pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints, sizeof(FLT ));
	pnt_strctFilm->u = (FLT *)calloc(pnt_config->int_Samples*pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints, sizeof(FLT ));
	pnt_strctFilm->v = (FLT *)calloc(pnt_config->int_Samples*pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints, sizeof(FLT ));
	pnt_strctFilm->w = (FLT *)calloc(pnt_config->int_Samples*pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints, sizeof(FLT ));
	pnt_strctFilm->p = (FLT *)calloc(pnt_config->int_Samples*pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints, sizeof(FLT ));
	pnt_strctFilm->gradRho = (FLT *)calloc(pnt_config->int_Samples*pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints, sizeof(FLT ));
	pnt_strctFilm->Lambda2 = (FLT *)calloc(pnt_config->int_Samples*pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints, sizeof(FLT ));
	pnt_strctFilm->MachNumber = (FLT *)calloc(pnt_config->int_Samples*pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints, sizeof(FLT ));
	pnt_strctFilm->time_dim = (FLT *)calloc(pnt_config->int_Samples, sizeof(FLT ));
}

void FreeMemoryStrctFilm(
		struct strct_Film * pnt_strctFilm)
{
	free(pnt_strctFilm->rho);
	free(pnt_strctFilm->u);
	free(pnt_strctFilm->v);
	free(pnt_strctFilm->w);
	free(pnt_strctFilm->p);
	free(pnt_strctFilm->gradRho);
	free(pnt_strctFilm->Lambda2);
	free(pnt_strctFilm->MachNumber);

}


void AllocMemoryStrctFlux(
		struct strct_configuration * pnt_config,
		struct strct_Flux * pnt_strctFlux)
{
	pnt_strctFlux->Mass = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctFlux->xiMomentum = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctFlux->etaMomentum = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctFlux->zetaMomentum = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
	pnt_strctFlux->Energy = (FLT *)calloc(pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells, sizeof(FLT ));
}

void FreeMemoryStrctFlux(
		struct strct_Flux * pnt_strctFlux)
{
	free(pnt_strctFlux->Mass);
	free(pnt_strctFlux->xiMomentum);
	free(pnt_strctFlux->etaMomentum);
	free(pnt_strctFlux->zetaMomentum);
	free(pnt_strctFlux->Energy);
}


void FreeMemoryMesh(
		struct strct_mesh * pnt_mesh,
		struct strct_configuration * pnt_config)
{
	free(pnt_mesh->x);
	free(pnt_mesh->y);
	free(pnt_mesh->z);
	free(pnt_mesh->xi_x);
	free(pnt_mesh->xi_y);
	free(pnt_mesh->xi_z);
	free(pnt_mesh->eta_x);
	free(pnt_mesh->eta_y);
	free(pnt_mesh->eta_z);
	free(pnt_mesh->zeta_x);
	free(pnt_mesh->zeta_y);
	free(pnt_mesh->zeta_z);
	free(pnt_mesh->jacobian);
	free(pnt_mesh->flag_IBC);

	//	Pressure Waves
	if(pnt_config->flag_PressureWaves==1)
	{
		free(pnt_mesh->flag_PressureWaves);
		free(pnt_mesh->startPressure_PressureWaves);
		free(pnt_mesh->startDensity_PressureWaves);
	}


	free(pnt_mesh->xiFluss_Faktor);
	free(pnt_mesh->etaFluss_Faktor);
	free(pnt_mesh->zetaFluss_Faktor);
}


//Das Rechengebiet wird initialisiert.
//Je nachdem welche Option gewählt wurde, werden die Werte mit Standardwerten oder mit
//vorhandenen Startlösungen belegt
void Initialize(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	int i,j,k,ijk;

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				pnt_U_lastStep->rho[ijk]=pnt_config->InitializeValues_rho0;
				pnt_U_lastStep->u[ijk]=pnt_config->u_inflow;
				pnt_U_lastStep->v[ijk]=pnt_config->v_inflow;
				pnt_U_lastStep->w[ijk]=pnt_config->w_inflow;
				pnt_U_lastStep->p[ijk]=pnt_config->InitializeValues_p0;
				pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+
						pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);

			}
		}
	}
}

void InitializeLaminarBoundary(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	int i,j,k,ijk,ijPlus1k;
	FLT rho,p,u,v;
	FLT x_start;
	FLT eta,delta_y,T_unendl,T_ad;
	FLT grenzschicht_dicke_x;

	rho=1.0;
	p=1.0;
	u=1.0;
	v=1.0;


	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				x_start=pnt_config->LaminarBoundary_xStart;
				p=1.0;
				rho=1.0;
				u=1.0;
				v=0.0;
				if(pnt_mesh->x[ijk]>x_start)
				{
					ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;
					delta_y=pnt_mesh->y[ijPlus1k]-pnt_mesh->y[ijk];
					T_unendl=p/rho;
					T_ad=T_unendl*(1+(pnt_config->gammaNumber-1.0)/2.0*pow(pnt_config->machNumber,2.0));

					grenzschicht_dicke_x=4.9*sqrt((pnt_mesh->x[ijk]-x_start)/pnt_config->reynoldsNumber);
					eta=(pnt_mesh->y[ijk]+delta_y*0.5)/
							grenzschicht_dicke_x;
					if (eta<1.0)
					{
//							approximation der blasius-lösung
//						u=1.0-1.4567*exp(-1.0*eta)-1.2956*exp(-1.0*eta)*(eta-1.0)-0.8392*exp(-2.0*eta);
//						Loesung nach Schroeder
						u=1.5*eta-0.5*eta*eta*eta;
//							Temperaturverlauf nach Busemann/Crocco
//							printf(">%f %f %f< ",pnt_config->gammaNumber,pnt_config->machNumber,T_ad);
						rho=p/(T_unendl*
								((pnt_config->gammaNumber-1.0)/2.0*pow(pnt_config->machNumber,2.0)*u*(1.0-u)+u*(T_unendl-T_ad)/T_unendl+T_ad));

					}
				}


				pnt_U_lastStep->rho[ijk]=rho;
				pnt_U_lastStep->u[ijk]=u;
				pnt_U_lastStep->v[ijk]=v;
				pnt_U_lastStep->w[ijk]=0.0;
				pnt_U_lastStep->p[ijk]=p;
				pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+
						pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
			}
		}
	}
}


void InitializeSpecialConditions(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	int i,j,k,ijk,ij0k;
	int gebiet;
	FLT rho,p,u,v,w,e;
	FLT eta,delta_y,T_unendl,T_ad;
	FLT distance;
	FLT theta1,theta2;
	FLT u_tmp,v_tmp;	

	gebiet=0;
	rho=1.0;
	p=1.0;
	u=1.0;
	v=1.0;
	w=0.0;

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;



				switch (pnt_config->int_specialInitializeType)
				{
//				Validatio: 2D-Riemann Test
				case 22:
					if((pnt_mesh->x[ijk]<0.5)&&(pnt_mesh->y[ijk]<0.5)){gebiet=1;} //LU
					if((pnt_mesh->x[ijk]<0.5)&&(pnt_mesh->y[ijk]>=0.5)){gebiet=0;} //LO
					if((pnt_mesh->x[ijk]>=0.5)&&(pnt_mesh->y[ijk]<0.5)){gebiet=3;} //RU
					if((pnt_mesh->x[ijk]>=0.5)&&(pnt_mesh->y[ijk]>=0.5)){gebiet=2;} //RO

					switch (gebiet)
					{
						case 0:p=1.0;rho=2.0;u=0.0;v=-0.3;break;
						case 1:p=0.4;rho=1.0625;u=0.0;v=0.2145;break;
						case 2:p=1.0;rho=1.0;u=0.0;v=-0.4;break;
						case 3:p=0.4;rho=0.5197;u=0.0;v=-1.1259;break;
					}
					break;

//				Validatio: 1D-Riemann Test
				case 2:
					if(pnt_mesh->x[ijk]<=pnt_config->InitializeValues_xBorder){gebiet=0;} //L
					if(pnt_mesh->x[ijk]>pnt_config->InitializeValues_xBorder){gebiet=1;} //R

					switch (gebiet)
					{
						case 0:
							p=pnt_config->InitializeValues_p0;
							rho=pnt_config->InitializeValues_rho0;
							u=pnt_config->InitializeValues_u0;
							v=0.0;
							break;
						case 1:
							p=pnt_config->InitializeValues_p1;
							rho=pnt_config->InitializeValues_rho1;
							u=pnt_config->InitializeValues_u1;
							v=0.0;
							break;
					}
					break;

//				Validation: interaction of oblique shock with wall
				case 3:
					if((pnt_mesh->y[ijk]>0.0)&&(pnt_mesh->x[ijk]<=(2.664-(pnt_mesh->y[ijk]+0.0001)/tan(27.90066761/180.0*MY_PI))))
					{gebiet=0;}
					else
					{gebiet=1;}

					switch (gebiet)
					{
					case 0:
						rho=pnt_config->InitializeValues_rho1;
						u=0.0;
						v=0.0;
						p=pnt_config->InitializeValues_p1;
						break;
					case 1:
						rho=pnt_config->InitializeValues_rho0;
						u=0.0;
						v=0.0;
						p=pnt_config->InitializeValues_p0;
						break;
					}
					break;

//				SVO
				case 4:
					if((pnt_mesh->y[ijk]>pnt_config->IBC_yKolben)||(pnt_mesh->x[ijk]<100.0))
					{gebiet=0;}
					else
					{gebiet=1;}

					switch (gebiet)
					{
						case 0:
							rho=pnt_config->InitializeValues_rho1;
							u=0.0;
							v=0.0;
							p=pnt_config->InitializeValues_p1;
							break;
						case 1:
							rho=pnt_config->InitializeValues_rho0;
							u=0.0;
							v=0.0;
							p=pnt_config->InitializeValues_p0;
							break;
					}

					break;

//				Laminar Boundary Layer for airfoil
				case 5:

					p=pnt_config->InitializeValues_p0;
					rho=pnt_config->InitializeValues_rho0;
					u=pnt_config->u_inflow;
					v=pnt_config->v_inflow;
					
					e=(0.5*((u*u)+(v*v))+p/rho/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
					
					distance=sqrt(
					pow((pnt_mesh->x[ijk]-0.),2.)
					+pow((pnt_mesh->y[ijk]-0.),2.)
					);					
					u=u*tanh(4*distance);
					v=v*tanh(4*distance);					
					
					theta1=u;
					theta2=v;
					
					if(strcmp(pnt_config->BC_Bottom,pnt_config->BCWallViscous)==0)
					{
						ij0k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+0*pnt_config->int_kMeshPointsGhostCells+k;
						delta_y=sqrt(
						(pnt_mesh->y[ijk]-pnt_mesh->y[ij0k])*(pnt_mesh->y[ijk]-pnt_mesh->y[ij0k])+
						(pnt_mesh->x[ijk]-pnt_mesh->x[ij0k])*(pnt_mesh->x[ijk]-pnt_mesh->x[ij0k])
						);
						T_unendl=p/rho;
						T_ad=T_unendl*(1+(pnt_config->gammaNumber-1.0)/2.0*pow(pnt_config->machNumber,2.0));

						eta=delta_y/(4.9*pnt_mesh->x[ijk]/sqrt(pnt_config->reynoldsNumber*pnt_mesh->x[ijk]));
						if (eta<7.0)
						{
//							approximation der blasius-lösung
							theta1=1.0-1.4567*exp(-1.0*eta)-1.2956*exp(-1.0*eta)*(eta-1.0)-0.8392*exp(-2.0*eta);
//							Temperaturverlauf nach Busemann/Crocco
//							printf(">%f %f %f< ",pnt_config->gammaNumber,pnt_config->machNumber,T_ad);
							rho=p/(T_unendl*((pnt_config->gammaNumber-1.0)/2.0*pow(pnt_config->machNumber,2.0)*u*(1.0-u)+u*(T_unendl-T_ad)/T_unendl+T_ad));

						
							u_tmp=fabs(-(theta1*pnt_mesh->eta_y[ijk]-theta2*pnt_mesh->xi_y[ijk]))/
							fabs((-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]));
							if((pnt_mesh->y[ijk])>=0.0)
							{
							v_tmp=-(-theta1*(pnt_mesh->eta_x[ijk])+theta2*(pnt_mesh->xi_x[ijk]))/
							(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]);				
							}
							else
							{
							v_tmp=-((-theta1*pnt_mesh->eta_x[ijk])+(theta2*pnt_mesh->xi_x[ijk]))/
							fabs(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]);	
							}
					
							u=u_tmp/sqrt(u_tmp*u_tmp+v_tmp*v_tmp)*theta1;
							v=v_tmp/sqrt(u_tmp*u_tmp+v_tmp*v_tmp)*theta1;
						
						}
					
						
					}
					
					p=(e-0.5*((u*u)+(v*v)))*rho*(pnt_config->gammaNumber-1.0)/pnt_config->Upsilon;					
					
					
					break;

//				Stoss-STK
				case 6:
					if(pnt_mesh->x[ijk]<-0.1+0.01*sin(2.0*MY_PI*(pnt_mesh->z[ijk]-0.025)/0.1)){gebiet=0;} //L
					else					 {gebiet=1;} //R


					switch (gebiet)
					{
						case 0:
							p=pnt_config->InitializeValues_p0;
							rho=pnt_config->InitializeValues_rho0;
							u=1.0;
							v=0.0;
							break;

						case 1:
							p=pnt_config->InitializeValues_p1;
							rho=pnt_config->InitializeValues_rho1;
							u=0.0;
							v=0.0;
							break;
					}
					break;

//				gitterangepasste Startloesung (Stromlinien entlang xi-Linien)
				case 7:
					distance=sqrt(
					pow((pnt_mesh->x[ijk]-0.),2.)
					+pow((pnt_mesh->y[ijk]-0.),2.)
					);
					p=pnt_config->InitializeValues_p0;
					rho=pnt_config->InitializeValues_rho0;

					theta1=pnt_config->u_inflow;
					theta2=pnt_config->v_inflow;
					if(distance<1.)
					{
					u_tmp=fabs(-(theta1*pnt_mesh->eta_y[ijk]-theta2*pnt_mesh->xi_y[ijk]))/
					fabs((-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]));
					if((pnt_mesh->y[ijk])>=0.0)
					{
					v_tmp=-(-theta1*(pnt_mesh->eta_x[ijk])+theta2*(pnt_mesh->xi_x[ijk]))/
					(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]);				
					}
					else
					{
					v_tmp=-((-theta1*pnt_mesh->eta_x[ijk])+(theta2*pnt_mesh->xi_x[ijk]))/
					fabs(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]);	
					}
					
					u=u_tmp/sqrt(u_tmp*u_tmp+v_tmp*v_tmp)*tanh(4*distance);
					v=v_tmp/sqrt(u_tmp*u_tmp+v_tmp*v_tmp)*tanh(4*distance);
					
					}
					else
					{
					u=pnt_config->u_inflow;
					v=pnt_config->v_inflow;
					}
					
					break;

				}

				pnt_U_lastStep->rho[ijk]=rho;
				pnt_U_lastStep->u[ijk]=u;
				pnt_U_lastStep->v[ijk]=v;
				pnt_U_lastStep->w[ijk]=w;
				pnt_U_lastStep->p[ijk]=p;
				pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
			}
		}
	}
}


void InitializeMPI(
		int argc,
		char *argv[],
		struct strct_configuration * pnt_config)
{
//	MPI_Cart_create(
//			MPI_COMM_WORLD,
//			pnt_config->int_meshDimensions,
//			pnt_config->MPI_intArray_NoCPUs,
//			pnt_config->MPI_intArray_PeriodicBC,
//			pnt_config->MPI_int_ReorderCPUs,
//			&pnt_config->MPI_comm);
//
//	//Zuweisung Koordinate <-> MyID//
//	MPI_Comm_size(
//			pnt_config->MPI_comm,
//			&pnt_config->MPI_size);
//	MPI_Comm_rank(
//			pnt_config->MPI_comm,
//			&pnt_config->MPI_rank);
//	MPI_Cart_coords(
//			pnt_config->MPI_comm,
//			pnt_config->MPI_rank,
//			pnt_config->int_meshDimensions,
//			pnt_config->MPI_coordinates);
//
//
//    //ID der Nachbarprozesse//
//    MPI_Cart_shift(pnt_config->MPI_comm, 0, 1, &pnt_config->InterfaceNeighbourLeft[0], &pnt_config->InterfaceNeighbourRight[0]);
//    MPI_Cart_shift(pnt_config->MPI_comm, 1, -1, &pnt_config->InterfaceNeighbourTop[0], &pnt_config->InterfaceNeighbourBottom[0]);
//    MPI_Cart_shift(pnt_config->MPI_comm, 2, -1, &pnt_config->InterfaceNeighbourInFront[0], &pnt_config->InterfaceNeighbourBehind[0]);
//
////    printf( "CPU %d of %d is ready!\n", pnt_config->MPI_rank+1, pnt_config->MPI_size);
//
////Bei einem oder 2 freien Rändern in i-Richtung müssen die GC berücksichtigt werden
//
//    pnt_config->MPI_datatype = MPI_FLT;
//    MPI_Type_size(pnt_config->MPI_datatype, &pnt_config->MPI_intMyDatatypeSize);
//
//    if (pnt_config->int_meshDimensions>=1)
//    {
//    	pnt_config->MPI_intArray_Subsize[0]=pnt_config->int_iMeshPoints;
//    	pnt_config->MPI_intArray_Size[0]=pnt_config->int_iMeshPoints*pnt_config->MPI_intArray_NoCPUs[0];
//    	pnt_config->MPI_intArray_StartPoint[0]=pnt_config->MPI_coordinates[0]*pnt_config->MPI_intArray_Subsize[0];
//    }
//
//    if (pnt_config->int_meshDimensions>=2)
//	{
//    	pnt_config->MPI_intArray_Subsize[1]=pnt_config->int_jMeshPoints;
//    	pnt_config->MPI_intArray_Size[1]=pnt_config->int_jMeshPoints*pnt_config->MPI_intArray_NoCPUs[1];
//    	pnt_config->MPI_intArray_StartPoint[1]=pnt_config->MPI_coordinates[1]*pnt_config->MPI_intArray_Subsize[1];
//	}
//
//    if (pnt_config->int_meshDimensions>=3)
//	{
//    	pnt_config->MPI_intArray_Subsize[2]=pnt_config->int_kMeshPoints;
//        pnt_config->MPI_intArray_Size[2]=pnt_config->int_kMeshPoints*pnt_config->MPI_intArray_NoCPUs[2];
//        pnt_config->MPI_intArray_StartPoint[2]=pnt_config->MPI_coordinates[2]*pnt_config->MPI_intArray_Subsize[2];
//	}
//
//
//
//
//    MPI_Type_create_subarray(
//    		pnt_config->int_meshDimensions,
//    		pnt_config->MPI_intArray_Size,
//    		pnt_config->MPI_intArray_Subsize,
//    		pnt_config->MPI_intArray_StartPoint,
//    		MPI_ORDER_C,
//    		pnt_config->MPI_datatype,
//    		&pnt_config->MPI_subarray);
//    MPI_Type_commit(&pnt_config->MPI_subarray);
}

void DefineParameters(struct strct_configuration * pnt_config)
{
	pnt_config->c0_dim=sqrt(pnt_config->gammaNumber*pnt_config->gasConstantNumber*pnt_config->T0_dim);
	pnt_config->u0_dim=pnt_config->c0_dim*pnt_config->machNumber;

	pnt_config->AlphaNonRef=0.25*pnt_config->L0_dim/pnt_config->u0_dim;

	pnt_config->IBC_yKolben=40.0;
	pnt_config->IBC_alphaKolben=27.5;
//	pnt_config->IBC_alphaKolben=14.0;

	pnt_config->IBC_MovingLastPosition=IBC_getActualPosition(pnt_config);
	pnt_config->IBC_MovingActualPosition=IBC_getActualPosition(pnt_config);

	pnt_config->is_avrg=0.0;
	pnt_config->is_avrg_counter=0.0;
	pnt_config->is_maximum=0.0;
	pnt_config->is_minimum=999.0;

	sprintf(pnt_config->BCFarfield,"BCFarfield");
	sprintf(pnt_config->BCInflow,"BCInflow");
	sprintf(pnt_config->BCOutflow,"BCOutflow");
	sprintf(pnt_config->BCOutflowSubsonic,"BCOutflowSubsonic");
	sprintf(pnt_config->BCWallInviscid,"BCWallInviscid");
	sprintf(pnt_config->BCWallViscous,"BCWallViscous");
	sprintf(pnt_config->BCInflowSupersonic,"BCInflowSupersonic");
	sprintf(pnt_config->BCInflowSubsonic,"BCInflowSubsonic");
	sprintf(pnt_config->BCWallViscousIsothermal,"BCWallViscousIsothermal");

	sprintf(pnt_config->BC_Left,"-");
	sprintf(pnt_config->BC_Right,"-");
	sprintf(pnt_config->BC_Top,"-");
	sprintf(pnt_config->BC_Bottom,"-");
	sprintf(pnt_config->BC_Behind,"-");
	sprintf(pnt_config->BC_InFront,"-");

	pnt_config->flag_NAN=0;

	pnt_config->start_Time=0.0;

	pnt_config->time_dim=0.0;
	pnt_config->time_dim_lastAction=0.0;

	pnt_config->int_IterationsBetweenSamples=(int)((pnt_config->int_TotalIterations-pnt_config->int_StartSampling)/pnt_config->int_Samples);
	if (pnt_config->int_IterationsBetweenSamples<1)
	{
		pnt_config->int_IterationsBetweenSamples=1000000;
	}

	pnt_config->int_StartIteration=0;
	pnt_config->int_actualIteration=pnt_config->int_StartIteration;

	pnt_config->comm_time=0.0;

	pnt_config->int_conservationEquations=5;

    if(SPACEORDER==5)
    {
//    	In Anlehnung an:
//    	Efficient Implementation of Weichted ENO Schemes, Jiang, G-S,1996
    	pnt_config->wenoEpsilon=1.e-6;
    }
    if(SPACEORDER==9)
    {
//    	In Anlehnung an:
//    	Monoticity Preserving Weighted Essentially Non-oscillatory Schemes with Increasingly High Order of Accuracy, Balsara D.,2000
//    	pnt_config->wenoEpsilon=1.e-10;
    	pnt_config->wenoEpsilon=1.e-6;
    }
//    In Anlehnung an:
//    Mapped weighted essentially non-oscillatory schemes: Achieving optimal order near critical points, Henrick, A.K.,2005
	pnt_config->wenoEpsilon=10.*pow(MY_FLT_MIN,0.5);

	pnt_config->wenoP=2.;
	pnt_config->wenoOptimalerKoeffizient_W9[0]=1./126.;
	pnt_config->wenoOptimalerKoeffizient_W9[1]=10./63.;
	pnt_config->wenoOptimalerKoeffizient_W9[2]=10./21.;
	pnt_config->wenoOptimalerKoeffizient_W9[3]=20./63.;
	pnt_config->wenoOptimalerKoeffizient_W9[4]=5./126.;

	pnt_config->wenoOptimalerKoeffizient_W5[0]=1./10.;
	pnt_config->wenoOptimalerKoeffizient_W5[1]=3./5.;
	pnt_config->wenoOptimalerKoeffizient_W5[2]=3./10.;


	pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient = (FLT *)calloc((SPACEORDER+1), sizeof(FLT ));
	pnt_config->ZD_Interpolation_Koeffizient = (FLT *)calloc((SPACEORDER+1), sizeof(FLT ));
	pnt_config->ZD_Ableitung_Koeffizient = (FLT *)calloc((SPACEORDER+2), sizeof(FLT ));

	if(SPACEORDER==9)
	{
		/*
		 * Die Koeffizienten fuer die viskosen Terme werden hier gespeichert. Die Mentalitaet ist hierbei dass die Ableitung des aeusseren
		 * Flusses dF/dx approximiert wird ueber dF/dx=(F_i+1/2-F_i-1/2)/dx.
		 * In den Zwischenpunkten werden dann die Ableitungen (du/dx) und Variablen (u,v,..) mit den unten stehenden Koeffizienten berechnet.
		 * In der finalen Berechnung der Flussableitung dF/dx wird schliesslich die hohe Ordnung (6 oder 10) erreicht.
		 */


		// checked (30.03.2015 Koeffizientenvergleich.xls)
		pnt_config->ZD_Interpolation_Koeffizient[0] =       1. / 1260.;
		pnt_config->ZD_Interpolation_Koeffizient[1] =     -23. / 2520.;
		pnt_config->ZD_Interpolation_Koeffizient[2] =     127. / 2520.;
		pnt_config->ZD_Interpolation_Koeffizient[3] =    -473. / 2520.;
		pnt_config->ZD_Interpolation_Koeffizient[4] =    1627. / 2520.;
		pnt_config->ZD_Interpolation_Koeffizient[5] =    1627. / 2520.;
		pnt_config->ZD_Interpolation_Koeffizient[6] =    -473. / 2520.;
		pnt_config->ZD_Interpolation_Koeffizient[7] =     127. / 2520.;
		pnt_config->ZD_Interpolation_Koeffizient[8] =     -23. / 2520.;
		pnt_config->ZD_Interpolation_Koeffizient[9] =       1. / 1260.;

		// checked (30.03.2015 Koeffizientenvergleich.xls)
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[0] =     -3. /  9450.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[1] =    351. / 75600.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[2] =  -2649. / 75600.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[3] =  15351. / 75600.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[4] =-110649. / 75600.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[5] = 110649. / 75600.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[6] = -15351. / 75600.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[7] =   2649. / 75600.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[8] =   -351. / 75600.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[9] =      3. /  9450.;

		// checked (30.03.2015 Koeffizientenvergleich.xls)
		// Im Unterschied zu Prof. Klioutchnikovs Version wird die Ableitung u_xi=(u_i+1/2-u_i-1/2)/xi nicht
		//mittels der zweifachen Anwendung von Interpolationskoeffizienten berechnet sondern DIREKT mittels dieser Ableitungskoeffizienten
		pnt_config->ZD_Ableitung_Koeffizient[0]=   -1./1260.;
		pnt_config->ZD_Ableitung_Koeffizient[1]=   25./2520.;
		pnt_config->ZD_Ableitung_Koeffizient[2]= -150./2520.;
		pnt_config->ZD_Ableitung_Koeffizient[3]=  600./2520.;
		pnt_config->ZD_Ableitung_Koeffizient[4]=-2100./2520.;
		pnt_config->ZD_Ableitung_Koeffizient[5]=    0.;
		pnt_config->ZD_Ableitung_Koeffizient[6]= 2100./2520.;
		pnt_config->ZD_Ableitung_Koeffizient[7]= -600./2520.;
		pnt_config->ZD_Ableitung_Koeffizient[8]=  150./2520.;
		pnt_config->ZD_Ableitung_Koeffizient[9]=  -25./2520.;
		pnt_config->ZD_Ableitung_Koeffizient[10]=   1./1260.;


		/*
		// checked (30.03.2015 Koeffizientenvergleich.xls)
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[0]=      24./75600.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[1]=    -375./75600.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[2]=    3000./75600.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[3]=  -18000./75600.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[4]=  126000./75600.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[5]= -221298./75600.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[6]=  126000./75600.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[7]=  -18000./75600.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[8]=    3000./75600.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[9]=    -375./75600.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[10]=     24./75600.;
		*/
	}
	else
	{
		// checked (30.03.2015 Koeffizientenvergleich.xls)
		pnt_config->ZD_Interpolation_Koeffizient[0] =   1. / 60.;
		pnt_config->ZD_Interpolation_Koeffizient[1] =  -8. / 60.;
		pnt_config->ZD_Interpolation_Koeffizient[2] =  37. / 60.;
		pnt_config->ZD_Interpolation_Koeffizient[3] =  37. / 60.;
		pnt_config->ZD_Interpolation_Koeffizient[4] =  -8. / 60.;
		pnt_config->ZD_Interpolation_Koeffizient[5] =   1. / 60.;

		// checked (30.03.2015 Koeffizientenvergleich.xls)
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[0] =    -1. / 90.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[1] =    25. / 180.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[2] =  -245. / 180.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[3] =   245. / 180.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[4] =   -25. / 180.;
		pnt_config->ZD_AbleitungZwischenPunkt_Koeffizient[5] =     1. / 90.;


		// checked (30.03.2015 Koeffizientenvergleich.xls)
		/*
		 * Im Unterschied zu Prof. Klioutchnikovs Version wird die Ableitung u_xi=(u_i+1/2-u_i-1/2)/xi nicht
		 * mittels der zweifachen Anwendung von Interpolationskoeffizienten berechnet sondern DIREKT mittels dieser Ableitungskoeffizienten
		 */
		pnt_config->ZD_Ableitung_Koeffizient[0]=    -1./60.;
		pnt_config->ZD_Ableitung_Koeffizient[1]=     3./20.;
		pnt_config->ZD_Ableitung_Koeffizient[2]=    -3./4.;
		pnt_config->ZD_Ableitung_Koeffizient[3]=     0.;
		pnt_config->ZD_Ableitung_Koeffizient[4]=     3./4.;
		pnt_config->ZD_Ableitung_Koeffizient[5]=    -3./20.;
		pnt_config->ZD_Ableitung_Koeffizient[6]=     1./60.;

		/*
		// checked (30.03.2015 Koeffizientenvergleich.xls)
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[0]=      2./180.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[1]=    -27./180.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[2]=    270./180.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[3]=   -490./180.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[4]=    270./180;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[5]=    -27./180.;
		pnt_config->ZD_ZweiteAbleitung_Koeffizient[6]=      2./180.;
		*/
	}


	pnt_config->deltaXi=1.0;
	pnt_config->deltaEta=1.0;
	pnt_config->deltaZeta=1.0;

	if(pnt_config->int_TimeOrder==4)
	{
		pnt_config->RK_U_n_Faktor[0]=1.;
		pnt_config->RK_U_n_Faktor[1]=1.;
		pnt_config->RK_U_n_Faktor[2]=1.;
		pnt_config->RK_U_n_Faktor[3]=1.;

		pnt_config->RK_U_ABC_Faktor[0]=0.;
		pnt_config->RK_U_ABC_Faktor[1]=0.;
		pnt_config->RK_U_ABC_Faktor[2]=0.;
		pnt_config->RK_U_ABC_Faktor[3]=0.;

		pnt_config->RK_Q_Faktor[0]=0.5;
		pnt_config->RK_Q_Faktor[1]=0.5;
		pnt_config->RK_Q_Faktor[2]=1.;
		pnt_config->RK_Q_Faktor[3]=1./6.;

		pnt_config->RK_Q_Summe_Flag[0]=0.;
		pnt_config->RK_Q_Summe_Flag[1]=0.;
		pnt_config->RK_Q_Summe_Flag[2]=0.;
		pnt_config->RK_Q_Summe_Flag[3]=1./6.;

		pnt_config->RK_Q_Summe_Faktor[0]=1.;
		pnt_config->RK_Q_Summe_Faktor[1]=2.;
		pnt_config->RK_Q_Summe_Faktor[2]=2.;
		pnt_config->RK_Q_Summe_Faktor[3]=0.;
	}
	if(pnt_config->int_TimeOrder==3)
	{
		pnt_config->RK_U_n_Faktor[0]=1.;
		pnt_config->RK_U_n_Faktor[1]=3./4.;
		pnt_config->RK_U_n_Faktor[2]=1./3.;

		pnt_config->RK_U_ABC_Faktor[0]=0.;
		pnt_config->RK_U_ABC_Faktor[1]=1./4.;
		pnt_config->RK_U_ABC_Faktor[2]=2./3.;

		pnt_config->RK_Q_Faktor[0]=1.;
		pnt_config->RK_Q_Faktor[1]=1./4.;
		pnt_config->RK_Q_Faktor[2]=2./3.;

		pnt_config->RK_Q_Summe_Flag[0]=0.;
		pnt_config->RK_Q_Summe_Flag[1]=0.;
		pnt_config->RK_Q_Summe_Flag[2]=0.;

		pnt_config->RK_Q_Summe_Faktor[0]=0.;
		pnt_config->RK_Q_Summe_Faktor[1]=0.;
		pnt_config->RK_Q_Summe_Faktor[2]=0.;
	}


	pnt_config->Upsilon=1.0/(pnt_config->gammaNumber*(pnt_config->machNumber*pnt_config->machNumber));
	pnt_config->Psi=1.0/pnt_config->reynoldsNumber;
	pnt_config->SutherlandConstant=100.4/pnt_config->T0_dim;

    //Gamma wird nach dem Import des Gitters bei postprocessload definiert

	pnt_config->flag_reinitialization=0;

	strcpy(pnt_config->ManufacturedSolution_L2_Delta_name,"Density");
	pnt_config->ManufacturedSolution_L2_last=1.0L;
	pnt_config->ManufacturedSolution_L2_last_pressure=1.0L;
	pnt_config->ManufacturedSolution_L2_counter=0;
	pnt_config->ManufacturedSolution_last_Q_Mass=1.0L;
	pnt_config->ManufacturedSolution_last_Q_xiMomentum=1.0L;
	pnt_config->ManufacturedSolution_last_Q_etaMomentum=1.0L;
	pnt_config->ManufacturedSolution_last_Q_Energy=1.0L;
	pnt_config->all_L2_norm_rho=0.0L;
	pnt_config->all_L2_norm_pressure=0.0L;
	pnt_config->all_Linf_norm_rho=0.0L;
	pnt_config->all_Linf_norm_pressure=0.0L;
}

void CalcRungeKutta(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q,
		struct strct_Flux * pnt_Q_sum,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf)
{
	int int_RKSchritt;
	FLT rho_letzterRK_Schritt;
//	-----------------------------------------------------
//	Erster bis Vierter RK-Schritt
//	Folgende Variablen enthalten die Faktoren für die 4 Schritte:
//	RKSchrittFaktor[4];
//	RKQsummeFlag[4];
//	RKQsummeFaktor[4];
//	-----------------------------------------------------

	DeleteQ(
			pnt_config,
			pnt_Q_sum);

	if(pnt_config->flag_constantZValues==1)
	{
		WriteConstantZValues(
				pnt_config,
				pnt_U_RK);

		changeMeshInto2D(
				pnt_config);
	}


	if((pnt_config->flag_IBC_Moving==1)&&(pnt_config->flag_IBC==1))
	{
		pnt_config->IBC_MovingLastPosition=pnt_config->IBC_MovingActualPosition;
		pnt_config->IBC_MovingActualPosition=IBC_getActualPosition(pnt_config);

		IBC_Actual2Last(
				pnt_config,
				pnt_mesh);

		IBC_prepare(
				pnt_config,
				pnt_mesh);

		IBC_BornCells(
				pnt_config,
				pnt_mesh,
				pnt_U_lastStep,
				pnt_U_RK);

        /*IBC_set(
                pnt_config,
                pnt_mesh,
                pnt_U_lastStep);*/

	}

	if(pnt_config->flag_PressureWaves==1)
	{
		inducePressureWaves(
				pnt_config,
				pnt_mesh,
				pnt_U_lastStep);
		inducePressureWaves(
				pnt_config,
				pnt_mesh,
				pnt_U_RK);
	}

	for(int_RKSchritt=0;int_RKSchritt<pnt_config->int_TimeOrder;int_RKSchritt++)
	//for(int_RKSchritt=0;int_RKSchritt<2;int_RKSchritt++)
	{
//				In Q werden die Flüsse für jede Richtung gespeichert. Am Anfang jedes RK-Schrittes muss dieser gelöscht werden.
		DeleteQ(
				pnt_config,
				pnt_Q);


//		Strömungsparameter (u,v,w,e,rho) werden übertragen
		TransferFlowParameterWithGhosts(
				pnt_config,
				pnt_mesh,
				pnt_U_RK);

		if (pnt_config->flag_Inviscid!=1)
		{
			TransferFlowParameterWithGhosts(
					pnt_config,
					pnt_mesh,
					pnt_U_RK);

			if (MESHDIMENSIONS==3)
			{
				TransferFlowParameterWithGhosts(
								pnt_config,
								pnt_mesh,
								pnt_U_RK);
			}
		}

		CalcValues(
				pnt_config,
				pnt_mesh,
				pnt_U_RK);

		SetAllBoundaryConditions(
				pnt_config,
				pnt_mesh,
				pnt_U_RK,
				pnt_U_lastStep);


//###########################################
//		REIBUNGSFREIE FLUSSBERECHNUNG
//##########################################
		if(pnt_config->flag_IBC==1)
		{
			IBC_ApplyBC4FluxInXi(
					pnt_config,
					pnt_mesh,
					pnt_U_RK);
		}

		CalcFluxesInXiDirection(
				pnt_config,
				pnt_mesh,
				pnt_U_RK,
				pnt_Flux,
				pnt_Flux_PlusHalf,
				pnt_Q);
		//ijk=4*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+4*pnt_config->int_kMeshPointsGhostCells+4;
		//if((pnt_mesh->x[ijk]==pnt_mesh->y[ijk])&&(pnt_mesh->y[ijk]==pnt_mesh->z[ijk]))
		//printf("Nach InviscidX: mass: %.10e - x: %.10e - y: %.10e - z: %.10e - energy: %.10e\n",pnt_Q->Mass[ijk],pnt_Q->xiMomentum[ijk],pnt_Q->etaMomentum[ijk],pnt_Q->zetaMomentum[ijk],pnt_Q->Energy[ijk]);
		//DeleteQ(pnt_config,pnt_Q);
		if(pnt_config->flag_IBC==1)
		{
			IBC_ApplyBC4FluxInEta(
					pnt_config,
					pnt_mesh,
					pnt_U_RK);
		}

		CalcFluxesInEtaDirection(
				pnt_config,
				pnt_mesh,
				pnt_U_RK,
				pnt_Flux,
				pnt_Flux_PlusHalf,
				pnt_Q);
		//ijk=4*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+4*pnt_config->int_kMeshPointsGhostCells+4;
		//if((pnt_mesh->x[ijk]==pnt_mesh->y[ijk])&&(pnt_mesh->y[ijk]==pnt_mesh->z[ijk]))
		//printf("Nach InviscidY: mass: %.10e - x: %.10e - y: %.10e - z: %.10e - energy: %.10e\n",pnt_Q->Mass[ijk],pnt_Q->xiMomentum[ijk],pnt_Q->etaMomentum[ijk],pnt_Q->zetaMomentum[ijk],pnt_Q->Energy[ijk]);
		//DeleteQ(pnt_config,pnt_Q);
#if MESHDIMENSIONS==3
		if(pnt_config->flag_IBC==1)
		{
			IBC_ApplyBC4FluxInZeta(
					pnt_config,
					pnt_mesh,
					pnt_U_RK);
		}

		CalcFluxesInZetaDirection(
				pnt_config,
				pnt_mesh,
				pnt_U_RK,
				pnt_Flux,
				pnt_Flux_PlusHalf,
				pnt_Q);
#endif
		//ijk=4*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+4*pnt_config->int_kMeshPointsGhostCells+4;
		//if((pnt_mesh->x[ijk]==pnt_mesh->y[ijk])&&(pnt_mesh->y[ijk]==pnt_mesh->z[ijk]))
		//printf("Nach InviscidZ: mass: %.10e - x: %.10e - y: %.10e - z: %.10e - energy: %.10e\n",pnt_Q->Mass[ijk],pnt_Q->xiMomentum[ijk],pnt_Q->etaMomentum[ijk],pnt_Q->zetaMomentum[ijk],pnt_Q->Energy[ijk]);
		//DeleteQ(pnt_config,pnt_Q);

//###########################################
//		REIBUNGSBEHAFTETE FLUSSBERECHNUNG
//##########################################
		if (pnt_config->flag_Inviscid!=1)
		{

			CalcDeviationsForDirectViscidFluxes(
						pnt_config,
						pnt_mesh,
						pnt_U_RK);

			if(pnt_config->flag_IBC==1)
			{
				IBC_ApplyBC4FluxInXi(
						pnt_config,
						pnt_mesh,
						pnt_U_RK);
			}

			CalcViscidFluxesInXiDirectionDirectly(
					pnt_config,
					pnt_mesh,
					pnt_U_RK,
					pnt_Q);
			//ijk=4*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+4*pnt_config->int_kMeshPointsGhostCells+4;
			//if((pnt_mesh->x[ijk]==pnt_mesh->y[ijk])&&(pnt_mesh->y[ijk]==pnt_mesh->z[ijk]))
			//printf("Nach ViscousX: mass: %.10e - x: %.10e - y: %.10e - z: %.10e - energy: %.10e\n",pnt_Q->Mass[ijk],pnt_Q->xiMomentum[ijk],pnt_Q->etaMomentum[ijk],pnt_Q->zetaMomentum[ijk],pnt_Q->Energy[ijk]);
			//DeleteQ(pnt_config,pnt_Q);
			if(pnt_config->flag_IBC==1)
			{
				IBC_ApplyBC4FluxInEta(
						pnt_config,
						pnt_mesh,
						pnt_U_RK);
			}
			CalcViscidFluxesInEtaDirectionDirectly(
					pnt_config,
					pnt_mesh,
					pnt_U_RK,
					pnt_Q);
			//ijk=4*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+4*pnt_config->int_kMeshPointsGhostCells+4;
			//if((pnt_mesh->x[ijk]==pnt_mesh->y[ijk])&&(pnt_mesh->y[ijk]==pnt_mesh->z[ijk]))
			//printf("Nach ViscousY: mass: %.10e - x: %.10e - y: %.10e - z: %.10e - energy: %.10e\n",pnt_Q->Mass[ijk],pnt_Q->xiMomentum[ijk],pnt_Q->etaMomentum[ijk],pnt_Q->zetaMomentum[ijk],pnt_Q->Energy[ijk]);
			//DeleteQ(pnt_config,pnt_Q);

#if MESHDIMENSIONS==3
			if(pnt_config->flag_IBC==1)
			{
				IBC_ApplyBC4FluxInZeta(
						pnt_config,
						pnt_mesh,
						pnt_U_RK);
			}

			CalcViscidFluxesInZetaDirectionDirectly(
					pnt_config,
					pnt_mesh,
					pnt_U_RK,
					pnt_Q);
			//ijk=4*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+4*pnt_config->int_kMeshPointsGhostCells+4;
			//if((pnt_mesh->x[ijk]==pnt_mesh->y[ijk])&&(pnt_mesh->y[ijk]==pnt_mesh->z[ijk]))
			//printf("Nach ViscousZ: mass: %.10e - x: %.10e - y: %.10e - z: %.10e - energy: %.10e\n",pnt_Q->Mass[ijk],pnt_Q->xiMomentum[ijk],pnt_Q->etaMomentum[ijk],pnt_Q->zetaMomentum[ijk],pnt_Q->Energy[ijk]);
			//DeleteQ(pnt_config,pnt_Q);
#endif
		}
		if(pnt_config->flag_rotation_symmetric==1)
		{
			AddRotationSymmetricFluxes(
					pnt_config,
					pnt_mesh,
					pnt_U_RK,
					pnt_Q);
		}

		if(pnt_config->flag_ManufacturedSolution==1)
		{
			AddManufacturedSolutionSource(
					pnt_config,
					pnt_mesh,
					pnt_Q);
		}

		//ijk=4*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+4*pnt_config->int_kMeshPointsGhostCells+4;
		//if((pnt_mesh->x[ijk]==pnt_mesh->y[ijk])&&(pnt_mesh->y[ijk]==pnt_mesh->z[ijk]))
		//printf("Nach AddManu: mass: %.10e - x: %.10e - y: %.10e - z: %.10e - energy: %.10e\n",pnt_Q->Mass[ijk],pnt_Q->xiMomentum[ijk],pnt_Q->etaMomentum[ijk],pnt_Q->zetaMomentum[ijk],pnt_Q->Energy[ijk]);
		//DeleteQ(pnt_config,pnt_Q);

		for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
		{
			for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
			{
				for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

					pnt_Q_sum->Mass[ijk]=pnt_Q_sum->Mass[ijk]+pnt_config->RK_Q_Summe_Faktor[int_RKSchritt]*pnt_Q->Mass[ijk];
					pnt_Q_sum->xiMomentum[ijk]=pnt_Q_sum->xiMomentum[ijk]+pnt_config->RK_Q_Summe_Faktor[int_RKSchritt]*pnt_Q->xiMomentum[ijk];
					pnt_Q_sum->etaMomentum[ijk]=pnt_Q_sum->etaMomentum[ijk]+pnt_config->RK_Q_Summe_Faktor[int_RKSchritt]*pnt_Q->etaMomentum[ijk];
					pnt_Q_sum->zetaMomentum[ijk]=pnt_Q_sum->zetaMomentum[ijk]+pnt_config->RK_Q_Summe_Faktor[int_RKSchritt]*pnt_Q->zetaMomentum[ijk];
					pnt_Q_sum->Energy[ijk]=pnt_Q_sum->Energy[ijk]+pnt_config->RK_Q_Summe_Faktor[int_RKSchritt]*pnt_Q->Energy[ijk];

					rho_letzterRK_Schritt=pnt_U_RK->rho[ijk];

					pnt_U_RK->rho[ijk]=
											pnt_config->RK_U_n_Faktor[int_RKSchritt]*pnt_U_lastStep->rho[ijk]+
											pnt_config->RK_U_ABC_Faktor[int_RKSchritt]*pnt_U_RK->rho[ijk]+
											pnt_config->numericalTau/pnt_mesh->jacobian[ijk]*
											(
											pnt_config->RK_Q_Faktor[int_RKSchritt]*pnt_Q->Mass[ijk]+
											pnt_config->RK_Q_Summe_Flag[int_RKSchritt]*pnt_Q_sum->Mass[ijk]
											);
					pnt_U_RK->u[ijk]=
											(pnt_config->RK_U_n_Faktor[int_RKSchritt]*pnt_U_lastStep->u[ijk]*pnt_U_lastStep->rho[ijk]+
											pnt_config->RK_U_ABC_Faktor[int_RKSchritt]*pnt_U_RK->u[ijk]*rho_letzterRK_Schritt+
											pnt_config->numericalTau/pnt_mesh->jacobian[ijk]*
											(
											pnt_config->RK_Q_Faktor[int_RKSchritt]*pnt_Q->xiMomentum[ijk]+
											pnt_config->RK_Q_Summe_Flag[int_RKSchritt]*pnt_Q_sum->xiMomentum[ijk]
											))/pnt_U_RK->rho[ijk];
					pnt_U_RK->v[ijk]=
											(pnt_config->RK_U_n_Faktor[int_RKSchritt]*pnt_U_lastStep->v[ijk]*pnt_U_lastStep->rho[ijk]+
											pnt_config->RK_U_ABC_Faktor[int_RKSchritt]*pnt_U_RK->v[ijk]*rho_letzterRK_Schritt+
											pnt_config->numericalTau/pnt_mesh->jacobian[ijk]*
											(
											pnt_config->RK_Q_Faktor[int_RKSchritt]*pnt_Q->etaMomentum[ijk]+
											pnt_config->RK_Q_Summe_Flag[int_RKSchritt]*pnt_Q_sum->etaMomentum[ijk]
											))/pnt_U_RK->rho[ijk];

#if MESHDIMENSIONS==3
					pnt_U_RK->w[ijk]=
											(pnt_config->RK_U_n_Faktor[int_RKSchritt]*pnt_U_lastStep->w[ijk]*pnt_U_lastStep->rho[ijk]+
											pnt_config->RK_U_ABC_Faktor[int_RKSchritt]*pnt_U_RK->w[ijk]*rho_letzterRK_Schritt+
											pnt_config->numericalTau/pnt_mesh->jacobian[ijk]*
											(
											pnt_config->RK_Q_Faktor[int_RKSchritt]*pnt_Q->zetaMomentum[ijk]+
											pnt_config->RK_Q_Summe_Flag[int_RKSchritt]*pnt_Q_sum->zetaMomentum[ijk]
											))/pnt_U_RK->rho[ijk];
#else
					pnt_U_RK->w[ijk]=0.0;
#endif

					pnt_U_RK->e[ijk]=
											(pnt_config->RK_U_n_Faktor[int_RKSchritt]*pnt_U_lastStep->e[ijk]*pnt_U_lastStep->rho[ijk]+
											pnt_config->RK_U_ABC_Faktor[int_RKSchritt]*pnt_U_RK->e[ijk]*rho_letzterRK_Schritt+
											pnt_config->numericalTau/pnt_mesh->jacobian[ijk]*
											(
											pnt_config->RK_Q_Faktor[int_RKSchritt]*pnt_Q->Energy[ijk]+
											pnt_config->RK_Q_Summe_Flag[int_RKSchritt]*pnt_Q_sum->Energy[ijk]
											))/pnt_U_RK->rho[ijk];

				}
			}
		}

		//Um sicherstellen, dass die Berechnungen innerhalb der IBC keine NAN/INF erzeugen,
		//werden die Berechnungen dort direkt ueberschrieben
		if(pnt_config->flag_IBC==1)
		{
			IBC_set(
					pnt_config,
					pnt_mesh,
					pnt_U_RK);
		}
	}

	if(pnt_config->flag_constantZValues==1)
	{
		changeMeshInto3D(
				pnt_config);
	}
	if(pnt_config->flag_PressureWaves==1)
	{
		inducePressureWaves(
				pnt_config,
				pnt_mesh,
				pnt_U_lastStep);
	}
}

void CalcValues(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
		int i,j,k,ijk;
		for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
		{
			for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
			{
				for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

	//Im Vergleich zu Klioutchnikovs Code fehlt hier die Subtraktion durch die Jakobische.
	//Diese steckt aber bereits im xi_x []= (y_eta*z_zeta-y_zeta*z_eta)/jacobian]
#if MESHDIMENSIONS==2
						pnt_U->theta1[ijk]=
								pnt_U->u[ijk]*pnt_mesh->xi_x[ijk]+
								pnt_U->v[ijk]*pnt_mesh->xi_y[ijk];

						pnt_U->theta2[ijk]=
								pnt_U->u[ijk]*pnt_mesh->eta_x[ijk]+
								pnt_U->v[ijk]*pnt_mesh->eta_y[ijk];

						pnt_U->p[ijk]=
								fabs
								(
								pnt_U->rho[ijk]*
								(pnt_config->gammaNumber-1)/
								pnt_config->Upsilon*
								(
								pnt_U->e[ijk]
								-0.5*
								(
								(pnt_U->u[ijk]*pnt_U->u[ijk])+
								(pnt_U->v[ijk]*pnt_U->v[ijk])
								)
								));
#else
						pnt_U->theta1[ijk]=
								pnt_U->u[ijk]*pnt_mesh->xi_x[ijk]+
								pnt_U->v[ijk]*pnt_mesh->xi_y[ijk]+
								pnt_U->w[ijk]*pnt_mesh->xi_z[ijk];
						pnt_U->theta2[ijk]=
								pnt_U->u[ijk]*pnt_mesh->eta_x[ijk]+
								pnt_U->v[ijk]*pnt_mesh->eta_y[ijk]+
								pnt_U->w[ijk]*pnt_mesh->eta_z[ijk];
						pnt_U->theta3[ijk]=
								pnt_U->u[ijk]*pnt_mesh->zeta_x[ijk]+
								pnt_U->v[ijk]*pnt_mesh->zeta_y[ijk]+
								pnt_U->w[ijk]*pnt_mesh->zeta_z[ijk];

						pnt_U->p[ijk]=
								fabs
								(
								pnt_U->rho[ijk]*
								(pnt_config->gammaNumber-1)/
								pnt_config->Upsilon*
								(
								pnt_U->e[ijk]
								-0.5*
								(
								(pnt_U->u[ijk]*pnt_U->u[ijk])+
								(pnt_U->v[ijk]*pnt_U->v[ijk])+
								(pnt_U->w[ijk]*pnt_U->w[ijk])
								)
								));
#endif

					pnt_U->T[ijk]=
							fabs(pnt_U->p[ijk]/pnt_U->rho[ijk]);


					pnt_U->c[ijk]=
							sqrt(pnt_config->Upsilon*
							pnt_config->gammaNumber*pnt_U->p[ijk]/pnt_U->rho[ijk]);

					pnt_U->mue[ijk]=((1.0+pnt_config->SutherlandConstant)*pow(pnt_U->p[ijk]/pnt_U->rho[ijk],1.5)/
							(pnt_U->p[ijk]/pnt_U->rho[ijk]+pnt_config->SutherlandConstant));

					if(pnt_config->flag_ManufacturedSolution==1)
					{
						pnt_U->mue[ijk]=1.0;
					}

				}
			}
		}
}


void WriteValuesFromU1ToU2(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U1,
		struct strct_U * pnt_U2)
{
	int i,j,k,ijk;
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				pnt_U2->rho[ijk]=pnt_U1->rho[ijk];
				pnt_U2->p[ijk]=pnt_U1->p[ijk];
				pnt_U2->u[ijk]=pnt_U1->u[ijk];
				pnt_U2->v[ijk]=pnt_U1->v[ijk];
				pnt_U2->w[ijk]=pnt_U1->w[ijk];
				pnt_U2->e[ijk]=pnt_U1->e[ijk];

			}
		}
	}
}

void WriteValuesFromUToFilm(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U,
		struct strct_Film * pnt_Film,
		struct strct_mesh * pnt_mesh)
{
	int i,j,k,ijk,ijkFilm;

	pnt_Film->time_dim[pnt_config->int_actualSample]=pnt_config->time_dim;

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal_original; k <= pnt_config->int_kEndReal_original; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkFilm=(i-pnt_config->int_iStartReal)*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints+(j-pnt_config->int_jStartReal)*pnt_config->int_kMeshPoints+(k-pnt_config->int_kStartReal_original)+
						pnt_config->int_actualSample*pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints;

				pnt_Film->rho[ijkFilm]=pnt_U->rho[ijk];
				pnt_Film->u[ijkFilm]=pnt_U->u[ijk];
				pnt_Film->v[ijkFilm]=pnt_U->v[ijk];
				pnt_Film->w[ijkFilm]=pnt_U->w[ijk];
				pnt_Film->p[ijkFilm]=pnt_U->p[ijk];
				pnt_Film->gradRho[ijkFilm]=pnt_U->gradRho[ijk];
				pnt_Film->MachNumber[ijkFilm]=pnt_U->MachNumber[ijk];
				pnt_Film->Lambda2[ijkFilm]=pnt_U->Lambda2[ijk];

//				Werte innerhalb des Solids (IBC)
				if(pnt_mesh->flag_IBC[ijk]==1)
				{
					if(pnt_config->flag_IBC_Moving==1){pnt_Film->u[ijkFilm]=//pnt_config->IBC_MovingSpeed/pnt_config->u0_dim;
						(pnt_config->IBC_MovingActualPosition-pnt_config->IBC_MovingLastPosition)/pnt_config->numericalTau ;}
					else{pnt_Film->u[ijkFilm]=-10.0;}
					pnt_Film->v[ijkFilm]=-10.0;
					pnt_Film->w[ijkFilm]=-10.0;
					pnt_Film->p[ijkFilm]=-10.0;
					pnt_Film->rho[ijkFilm]=-10.0;
					pnt_Film->gradRho[ijkFilm]=-10.0;
					pnt_Film->MachNumber[ijkFilm]=-10.0;
					pnt_Film->Lambda2[ijkFilm]=-10.0;
				}
			}
		}
	}
}


//Die dreidimensionalen Matrizen fuer x,y,z... muessen in einen eindimensionalen Array überführt werden.
//Das entspricht einer Transformation eines 4D-Arrays nach 1D
//i spricht Elemente an, über j werden Scheiben angesprochen und über k ganze Flächen.
//Mit den Werten 0-12 werden die Würfel (3D-Matrix) der einzelnen Parameter (x,y,z) angesprochen
void WriteValuesFromMeshToBuffer(
		FLT * buffer,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		int iStart,
		int iEnd,
		int jStart,
		int jEnd,
		int kStart,
		int kEnd)
{
	int i,j,k,ijk,ijk13;
	int i_buffer,j_buffer,k_buffer;
	int imax,jmax,kmax;

	imax=(iEnd-iStart+1);
	jmax=(jEnd-jStart+1);
	kmax=(kEnd-kStart+1);

	for (i=iStart; i <= iEnd; i++)
	{
		for (j=jStart; j <= jEnd; j++)
		{
			for (k=kStart; k <= kEnd; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				i_buffer=(i-iStart);
				j_buffer=(j-jStart);
				k_buffer=(k-kStart);

				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*0;
				buffer[ijk13]=pnt_mesh->x[ijk];
				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*1;
				buffer[ijk13]=pnt_mesh->y[ijk];
				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*2;
				buffer[ijk13]=pnt_mesh->z[ijk];

				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*3;
				buffer[ijk13]=pnt_mesh->xi_x[ijk];
				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*4;
				buffer[ijk13]=pnt_mesh->xi_y[ijk];
				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*5;
				buffer[ijk13]=pnt_mesh->xi_z[ijk];

				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*6;
				buffer[ijk13]=pnt_mesh->eta_x[ijk];
				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*7;
				buffer[ijk13]=pnt_mesh->eta_y[ijk];
				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*8;
				buffer[ijk13]=pnt_mesh->eta_z[ijk];

				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*9;
				buffer[ijk13]=pnt_mesh->zeta_x[ijk];
				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*10;
				buffer[ijk13]=pnt_mesh->zeta_y[ijk];
				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*11;
				buffer[ijk13]=pnt_mesh->zeta_z[ijk];

				ijk13=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*12;
				buffer[ijk13]=pnt_mesh->jacobian[ijk];
			}
		}
	}
}

void WriteValuesFromBufferToMesh(
		FLT * buffer,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		int iStart,
		int iEnd,
		int jStart,
		int jEnd,
		int kStart,
		int kEnd,
		int interface)
{
	int i,j,k,ijk,ijk13;

	int i_buffer,j_buffer,k_buffer;

	for (i=iStart; i <= iEnd; i++)
	{
		for (j=jStart; j <= jEnd; j++)
		{
			for (k=kStart; k <= kEnd; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				i_buffer=(i-iStart);
				j_buffer=(j-jStart);
				k_buffer=(k-kStart);


				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,0);
				pnt_mesh->x[ijk]=buffer[ijk13]+pnt_config->Translation[interface][0];
				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,1);
				pnt_mesh->y[ijk]=buffer[ijk13]+pnt_config->Translation[interface][1];
				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,2);
				pnt_mesh->z[ijk]=buffer[ijk13]+pnt_config->Translation[interface][2];


				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,3);
				pnt_mesh->xi_x[ijk]=pnt_config->MPI_dblTransformation_xi_x[interface]*buffer[ijk13];
				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,4);
				pnt_mesh->xi_y[ijk]=pnt_config->MPI_dblTransformation_xi_y[interface]*buffer[ijk13];
				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,5);
				pnt_mesh->xi_z[ijk]=pnt_config->MPI_dblTransformation_xi_z[interface]*buffer[ijk13];

				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,6);
				pnt_mesh->eta_x[ijk]=pnt_config->MPI_dblTransformation_eta_x[interface]*buffer[ijk13];
				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,7);
				pnt_mesh->eta_y[ijk]=pnt_config->MPI_dblTransformation_eta_y[interface]*buffer[ijk13];
				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,8);
				pnt_mesh->eta_z[ijk]=pnt_config->MPI_dblTransformation_eta_z[interface]*buffer[ijk13];

				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,9);
				pnt_mesh->zeta_x[ijk]=pnt_config->MPI_dblTransformation_zeta_x[interface]*buffer[ijk13];
				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,10);
				pnt_mesh->zeta_y[ijk]=pnt_config->MPI_dblTransformation_zeta_y[interface]*buffer[ijk13];
				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,11);
				pnt_mesh->zeta_z[ijk]=pnt_config->MPI_dblTransformation_zeta_z[interface]*buffer[ijk13];

				ijk13=get_ijkTransformMesh(pnt_config,interface,i_buffer,j_buffer,k_buffer,12);
				pnt_mesh->jacobian[ijk]=buffer[ijk13];
			}
		}
	}
}


//Die dreidimensionalen Matrizen für x,y,z... müssen in einen eindimensionalen Array überführt werden.
//Das entspricht einer Transformation eines 4D-Arrays nach 1D
//i spricht Elemente an, über j werden Scheiben angesprochen und über k ganze Flächen.
//Mit den Werten 0-12 werden die Würfel (3D-Matrix) der einzelnen Parameter (p,rho,u,v,w) angesprochen
void WriteValuesFromUToBuffer(
		FLT * buffer,
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U,
		int iStart,
		int iEnd,
		int jStart,
		int jEnd,
		int kStart,
		int kEnd)
{
	int i,j,k,ijk,ijk5,ijk4;
	int i_buffer,j_buffer,k_buffer;
	int imax,jmax,kmax;

	imax=(iEnd-iStart+1);
	jmax=(jEnd-jStart+1);
	kmax=(kEnd-kStart+1);

	for (i=iStart; i <= iEnd; i++)
	{
		for (j=jStart; j <= jEnd; j++)
		{
			for (k=kStart; k <= kEnd; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				i_buffer=(i-iStart);
				j_buffer=(j-jStart);
				k_buffer=(k-kStart);

				if(MESHDIMENSIONS==2)
				{
					ijk4=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*0;
					buffer[ijk4]=pnt_U->u[ijk];
					ijk4=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*1;
					buffer[ijk4]=pnt_U->v[ijk];
					ijk4=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*2;
					buffer[ijk4]=pnt_U->e[ijk];
					ijk4=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*3;
					buffer[ijk4]=pnt_U->rho[ijk];
				}
				else
				{
					ijk5=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*0;
					buffer[ijk5]=pnt_U->u[ijk];
					ijk5=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*1;
					buffer[ijk5]=pnt_U->v[ijk];
					ijk5=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*2;
					buffer[ijk5]=pnt_U->w[ijk];
					ijk5=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*3;
					buffer[ijk5]=pnt_U->e[ijk];
					ijk5=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*4;
					buffer[ijk5]=pnt_U->rho[ijk];
				}
			}
		}
	}
}


void WriteValuesFromUAndMeshToBuffer(
		FLT * buffer,
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U,
		struct strct_mesh * pnt_mesh,
		int iStart,
		int iEnd,
		int jStart,
		int jEnd,
		int kStart,
		int kEnd)
{
	int i,j,k,ijk,ijk8;
	int i_buffer,j_buffer,k_buffer;
	int imax,jmax,kmax;

	imax=(iEnd-iStart+1);
	jmax=(jEnd-jStart+1);
	kmax=(kEnd-kStart+1);

	for (i=iStart; i <= iEnd; i++)
	{
		for (j=jStart; j <= jEnd; j++)
		{
			for (k=kStart; k <= kEnd; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				i_buffer=(i-iStart);
				j_buffer=(j-jStart);
				k_buffer=(k-kStart);

				ijk8=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*0;
				buffer[ijk8]=pnt_mesh->x[ijk];
				ijk8=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*1;
				buffer[ijk8]=pnt_mesh->y[ijk];
				ijk8=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*2;
				buffer[ijk8]=pnt_mesh->z[ijk];
				ijk8=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*3;
				buffer[ijk8]=pnt_U->u[ijk];
				ijk8=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*4;
				buffer[ijk8]=pnt_U->v[ijk];
				ijk8=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*5;
				buffer[ijk8]=pnt_U->w[ijk];
				ijk8=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*6;
				buffer[ijk8]=pnt_U->p[ijk];
				ijk8=i_buffer+j_buffer*imax+k_buffer*imax*jmax+imax*jmax*kmax*7;
				buffer[ijk8]=pnt_U->rho[ijk];
			}
		}
	}
}

void WriteValuesFromBufferToU(
		FLT * buffer,
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U,
		int iStart,
		int iEnd,
		int jStart,
		int jEnd,
		int kStart,
		int kEnd,
		int interface)
{
	int i,j,k,ijk,ijk5,ijk4;
	int i_buffer,j_buffer,k_buffer;


	for (i=iStart; i <= iEnd; i++)
	{
		for (j=jStart; j <= jEnd; j++)
		{
			for (k=kStart; k <= kEnd; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				i_buffer=(i-iStart);
				j_buffer=(j-jStart);
				k_buffer=(k-kStart);

				if(MESHDIMENSIONS==2)
				{
					ijk4=get_ijkTransformWithGhosts(pnt_config,interface,i_buffer,j_buffer,k_buffer,0);
					pnt_U->u[ijk]=buffer[ijk4];
					ijk4=get_ijkTransformWithGhosts(pnt_config,interface,i_buffer,j_buffer,k_buffer,1);
					pnt_U->v[ijk]=buffer[ijk4];
					ijk4=get_ijkTransformWithGhosts(pnt_config,interface,i_buffer,j_buffer,k_buffer,2);
					pnt_U->e[ijk]=buffer[ijk4];
					ijk4=get_ijkTransformWithGhosts(pnt_config,interface,i_buffer,j_buffer,k_buffer,3);
					pnt_U->rho[ijk]=buffer[ijk4];
				}
				else
				{
					ijk5=get_ijkTransformWithGhosts(pnt_config,interface,i_buffer,j_buffer,k_buffer,0);
					pnt_U->u[ijk]=buffer[ijk5];
					ijk5=get_ijkTransformWithGhosts(pnt_config,interface,i_buffer,j_buffer,k_buffer,1);
					pnt_U->v[ijk]=buffer[ijk5];
					ijk5=get_ijkTransformWithGhosts(pnt_config,interface,i_buffer,j_buffer,k_buffer,2);
					pnt_U->w[ijk]=buffer[ijk5];
					ijk5=get_ijkTransformWithGhosts(pnt_config,interface,i_buffer,j_buffer,k_buffer,3);
					pnt_U->e[ijk]=buffer[ijk5];
					ijk5=get_ijkTransformWithGhosts(pnt_config,interface,i_buffer,j_buffer,k_buffer,4);
					pnt_U->rho[ijk]=buffer[ijk5];
				}
			}
		}
	}
}


void DeleteQ(
		struct strct_configuration * pnt_config,
		struct strct_Flux * pnt_Q)
{
	int i,j,k,ijk;
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				pnt_Q->Mass[ijk]=0.0;
				pnt_Q->xiMomentum[ijk]=0.0;
				pnt_Q->etaMomentum[ijk]=0.0;
				pnt_Q->zetaMomentum[ijk]=0.0;
				pnt_Q->Energy[ijk]=0.0;
			}
		}
	}
}

//Bei den Randbedingungen ist auf die Reihenfolge zu achten.
//1)Die i-Randbedingungen werden gesetzt, dann der x-Transfer: reale j,k
//2)Die j-Randbedingungen werden gesetzt, dann der y-Transfer: reale k i mit GC
//3)Die k-Randbedingungen werden gesetzt, dann der z-Transfer: i,j mit GC
//void TransferFlowParameter(
//		struct strct_configuration * pnt_config,
//		struct strct_mesh * pnt_mesh,
//		struct strct_U * pnt_U,
//		FLT * bufferSendLeft,
//		FLT * bufferSendRight,
//		FLT * bufferSendBottom,
//		FLT * bufferSendTop,
//		FLT * bufferSendBehind,
//		FLT * bufferSendInFront,
//		FLT * bufferRecieveLeft,
//		FLT * bufferRecieveRight,
//		FLT * bufferRecieveBottom,
//		FLT * bufferRecieveTop,
//		FLT * bufferRecieveBehind,
//		FLT * bufferRecieveInFront)
//{
////#############################################################
////Die Buffer für den Transport in i-Richtung werden geschrieben
////#############################################################
//	WriteValuesFromUToBuffer(
//		bufferSendLeft,
//		pnt_config,
//		pnt_U,
//		pnt_config->int_iStartReal,
//		(pnt_config->int_iStartReal+((pnt_config->int_order+1)/2-1)),
//		pnt_config->int_jStartReal,
//		pnt_config->int_jEndReal,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	WriteValuesFromUToBuffer(
//		bufferSendRight,
//		pnt_config,
//		pnt_U,
//		(pnt_config->int_iEndReal-((pnt_config->int_order+1)/2-1)),
//		pnt_config->int_iEndReal,
//		pnt_config->int_jStartReal,
//		pnt_config->int_jEndReal,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
////	Der Transfer wird durchgeführt. Erst geschieht die Kommunikation nach links, dann nach rechts
//	MPI_Sendrecv(
//		bufferSendRight,
//		pnt_config->MPI_intSizeITransfer,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourRight[0],
//		10,
//		bufferRecieveLeft,
//		pnt_config->MPI_intSizeITransfer,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourLeft[0],
//		10,
//		pnt_config->MPI_comm,
//		&pnt_config->MPI_status);
//	MPI_Sendrecv(
//		bufferSendLeft,
//		pnt_config->MPI_intSizeITransfer,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourLeft[0],
//		11,
//		bufferRecieveRight,
//		pnt_config->MPI_intSizeITransfer,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourRight[0],
//		11,
//		pnt_config->MPI_comm,
//		&pnt_config->MPI_status);
//	//Die empfangene Buffer für wird in die GhostCells geschrieben
//	if (pnt_config->InterfaceNeighbourLeft[0]!=MPI_PROC_NULL)
//	{
//	WriteValuesFromBufferToU(
//		bufferRecieveLeft,
//		pnt_config,
//		pnt_U,
//		pnt_mesh,
//		pnt_config->int_iStartGhosts,
//		pnt_config->int_iStartReal-1,
//		pnt_config->int_jStartReal,
//		pnt_config->int_jEndReal,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	}
//
//	if (pnt_config->InterfaceNeighbourRight[0]!=MPI_PROC_NULL)
//	{
//	WriteValuesFromBufferToU(
//		bufferRecieveRight,
//		pnt_config,
//		pnt_U,
//		pnt_mesh,
//		pnt_config->int_iEndReal+1,
//		pnt_config->int_iEndGhosts,
//		pnt_config->int_jStartReal,
//		pnt_config->int_jEndReal,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	}
//
//
////#############################################################
////Die Buffer für den Transport in j-Richtung werden geschrieben
////#############################################################
//	WriteValuesFromUToBuffer(
//		bufferSendBottom,
//		pnt_config,
//		pnt_U,
//		pnt_config->int_iStartReal,
//		pnt_config->int_iEndReal,
//		pnt_config->int_jStartReal,
//		(pnt_config->int_jStartReal+((pnt_config->int_order+1)/2-1)),
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	WriteValuesFromUToBuffer(
//		bufferSendTop,
//		pnt_config,
//		pnt_U,
//		pnt_config->int_iStartReal,
//		pnt_config->int_iEndReal,
//		(pnt_config->int_jEndReal-((pnt_config->int_order+1)/2-1)),
//		pnt_config->int_jEndReal,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
////	Der Transfer wird durchgeführt. Erst geschieht die Kommunikation nach unten, dann nach oben
//	MPI_Sendrecv(
//		bufferSendTop,
//		pnt_config->MPI_intSizeJTransfer,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourTop[0],
//		12,
//		bufferRecieveBottom,
//		pnt_config->MPI_intSizeJTransfer,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourBottom[0],
//		12,
//		pnt_config->MPI_comm,
//		&pnt_config->MPI_status);
//	MPI_Sendrecv(
//		bufferSendBottom,
//		pnt_config->MPI_intSizeJTransfer,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourBottom[0],
//		13,
//		bufferRecieveTop,
//		pnt_config->MPI_intSizeJTransfer,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourTop[0],
//		13,
//		pnt_config->MPI_comm,
//		&pnt_config->MPI_status);
//	//Die empfangene Buffer wird in die GhostCells geschrieben
//	if (pnt_config->InterfaceNeighbourBottom[0]!=MPI_PROC_NULL)
//	{
//	WriteValuesFromBufferToU(
//		bufferRecieveBottom,
//		pnt_config,
//		pnt_U,
//		pnt_mesh,
//		pnt_config->int_iStartReal,
//		pnt_config->int_iEndReal,
//		pnt_config->int_jStartGhosts,
//		pnt_config->int_jStartReal-1,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	}
//
//
//	if (pnt_config->InterfaceNeighbourTop[0]!=MPI_PROC_NULL)
//	{
//	WriteValuesFromBufferToU(
//		bufferRecieveTop,
//		pnt_config,
//		pnt_U,
//		pnt_mesh,
//		pnt_config->int_iStartReal,
//		pnt_config->int_iEndReal,
//		pnt_config->int_jEndReal+1,
//		pnt_config->int_jEndGhosts,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	}
//
//
////#############################################################
////Die Buffer für den Transport in k-Richtung werden geschrieben
////#############################################################
//	if(pnt_config->int_kMeshPoints>1)
//	{
//		WriteValuesFromUToBuffer(
//			bufferSendBehind,
//			pnt_config,
//			pnt_U,
//			pnt_config->int_iStartReal,
//			pnt_config->int_iEndReal,
//			pnt_config->int_jStartReal,
//			pnt_config->int_jEndReal,
//			pnt_config->int_kStartReal,
//			(pnt_config->int_kStartReal+((pnt_config->int_order+1)/2-1)));
//		WriteValuesFromUToBuffer(
//			bufferSendInFront,
//			pnt_config,
//			pnt_U,
//			pnt_config->int_iStartReal,
//			pnt_config->int_iEndReal,
//			pnt_config->int_jStartReal,
//			pnt_config->int_jEndReal,
//			(pnt_config->int_kEndReal-((pnt_config->int_order+1)/2-1)),
//			pnt_config->int_kEndReal);
//	//	Der Transfer wird durchgeführt. Erst geschieht die Kommunikation nach unten, dann nach oben
//		MPI_Sendrecv(
//			bufferSendInFront,
//			pnt_config->MPI_intSizeKTransfer,
//			MPI_FLT,
//			pnt_config->InterfaceNeighbourInFront[0],
//			14,
//			bufferRecieveBehind,
//			pnt_config->MPI_intSizeKTransfer,
//			MPI_FLT,
//			pnt_config->InterfaceNeighbourBehind[0],
//			14,
//			pnt_config->MPI_comm,
//			&pnt_config->MPI_status);
//		MPI_Sendrecv(
//			bufferSendBehind,
//			pnt_config->MPI_intSizeKTransfer,
//			MPI_FLT,
//			pnt_config->InterfaceNeighbourBehind[0],
//			15,
//			bufferRecieveInFront,
//			pnt_config->MPI_intSizeKTransfer,
//			MPI_FLT,
//			pnt_config->InterfaceNeighbourInFront[0],
//			15,
//			pnt_config->MPI_comm,
//			&pnt_config->MPI_status);
//		//Die empfangene Buffer wird in die GhostCells geschrieben
//		if (pnt_config->InterfaceNeighbourBehind[0]!=MPI_PROC_NULL)
//		{
//		WriteValuesFromBufferToU(
//			bufferRecieveBehind,
//			pnt_config,
//			pnt_U,
//			pnt_mesh,
//			pnt_config->int_iStartReal,
//			pnt_config->int_iEndReal,
//			pnt_config->int_jStartReal,
//			pnt_config->int_jEndReal,
//			pnt_config->int_kStartGhosts,
//			pnt_config->int_kStartReal-1);
//		}
//
//		if (pnt_config->InterfaceNeighbourInFront[0]!=MPI_PROC_NULL)
//		{
//		WriteValuesFromBufferToU(
//			bufferRecieveInFront,
//			pnt_config,
//			pnt_U,
//			pnt_mesh,
//			pnt_config->int_iStartReal,
//			pnt_config->int_iEndReal,
//			pnt_config->int_jStartReal,
//			pnt_config->int_jEndReal,
//			pnt_config->int_kEndReal+1,
//			pnt_config->int_kEndGhosts);
//		}
//	}
//}


//void TransferViscidFlowParameter(
//		struct strct_configuration * pnt_config,
//		struct strct_mesh * pnt_mesh,
//		struct strct_ZD * pnt_ZD,
//		FLT * bufferSendLeftViscid,
//		FLT * bufferSendRightViscid,
//		FLT * bufferSendBottomViscid,
//		FLT * bufferSendTopViscid,
//		FLT * bufferSendBehindViscid,
//		FLT * bufferSendInFrontViscid,
//		FLT * bufferRecieveLeftViscid,
//		FLT * bufferRecieveRightViscid,
//		FLT * bufferRecieveBottomViscid,
//		FLT * bufferRecieveTopViscid,
//		FLT * bufferRecieveBehindViscid,
//		FLT * bufferRecieveInFrontViscid)
//{
////#############################################################
////Die Buffer für den Transport in i-Richtung werden geschrieben
////#############################################################
//	WriteValuesFromZDToBuffer(
//		bufferSendLeftViscid,
//		pnt_config,
//		pnt_ZD,
//		pnt_config->int_iStartReal,
//		(pnt_config->int_iStartReal+((pnt_config->int_order+1)/2-1)),
//		pnt_config->int_jStartReal,
//		pnt_config->int_jEndReal,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	WriteValuesFromZDToBuffer(
//		bufferSendRightViscid,
//		pnt_config,
//		pnt_ZD,
//		(pnt_config->int_iEndReal-((pnt_config->int_order+1)/2-1)),
//		pnt_config->int_iEndReal,
//		pnt_config->int_jStartReal,
//		pnt_config->int_jEndReal,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
////	Der Transfer wird durchgeführt. Erst geschieht die Kommunikation nach links, dann nach rechts
//	MPI_Sendrecv(
//		bufferSendRightViscid,
//		pnt_config->MPI_intSizeITransferViscid,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourRight[0],
//		16,
//		bufferRecieveLeftViscid,
//		pnt_config->MPI_intSizeITransferViscid,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourLeft[0],
//		16,
//		pnt_config->MPI_comm,
//		&pnt_config->MPI_status);
//	MPI_Sendrecv(
//		bufferSendLeftViscid,
//		pnt_config->MPI_intSizeITransferViscid,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourLeft[0],
//		17,
//		bufferRecieveRightViscid,
//		pnt_config->MPI_intSizeITransferViscid,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourRight[0],
//		17,
//		pnt_config->MPI_comm,
//		&pnt_config->MPI_status);
//	//Die empfangene Buffer für wird in die GhostCells geschrieben
//	if (pnt_config->InterfaceNeighbourLeft[0]!=MPI_PROC_NULL)
//	{
//	WriteValuesFromBufferToZD(
//		bufferRecieveLeftViscid,
//		pnt_config,
//		pnt_ZD,
//		pnt_mesh,
//		pnt_config->int_iStartGhosts,
//		pnt_config->int_iStartReal-1,
//		pnt_config->int_jStartReal,
//		pnt_config->int_jEndReal,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	}
//
//	if (pnt_config->InterfaceNeighbourRight[0]!=MPI_PROC_NULL)
//	{
//	WriteValuesFromBufferToZD(
//		bufferRecieveRightViscid,
//		pnt_config,
//		pnt_ZD,
//		pnt_mesh,
//		pnt_config->int_iEndReal+1,
//		pnt_config->int_iEndGhosts,
//		pnt_config->int_jStartReal,
//		pnt_config->int_jEndReal,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	}
//
//
////#############################################################
////Die Buffer für den Transport in j-Richtung werden geschrieben
////#############################################################
//	WriteValuesFromZDToBuffer(
//		bufferSendBottomViscid,
//		pnt_config,
//		pnt_ZD,
//		pnt_config->int_iStartReal,
//		pnt_config->int_iEndReal,
//		pnt_config->int_jStartReal,
//		(pnt_config->int_jStartReal+((pnt_config->int_order+1)/2-1)),
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	WriteValuesFromZDToBuffer(
//		bufferSendTopViscid,
//		pnt_config,
//		pnt_ZD,
//		pnt_config->int_iStartReal,
//		pnt_config->int_iEndReal,
//		(pnt_config->int_jEndReal-((pnt_config->int_order+1)/2-1)),
//		pnt_config->int_jEndReal,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
////	Der Transfer wird durchgeführt. Erst geschieht die Kommunikation nach unten, dann nach oben
//	MPI_Sendrecv(
//		bufferSendTopViscid,
//		pnt_config->MPI_intSizeJTransferViscid,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourTop[0],
//		18,
//		bufferRecieveBottomViscid,
//		pnt_config->MPI_intSizeJTransferViscid,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourBottom[0],
//		18,
//		pnt_config->MPI_comm,
//		&pnt_config->MPI_status);
//	MPI_Sendrecv(
//		bufferSendBottomViscid,
//		pnt_config->MPI_intSizeJTransferViscid,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourBottom[0],
//		19,
//		bufferRecieveTopViscid,
//		pnt_config->MPI_intSizeJTransferViscid,
//		MPI_FLT,
//		pnt_config->InterfaceNeighbourTop[0],
//		19,
//		pnt_config->MPI_comm,
//		&pnt_config->MPI_status);
//	//Die empfangene Buffer wird in die GhostCells geschrieben
//	if (pnt_config->InterfaceNeighbourBottom[0]!=MPI_PROC_NULL)
//	{
//	WriteValuesFromBufferToZD(
//		bufferRecieveBottomViscid,
//		pnt_config,
//		pnt_ZD,
//		pnt_mesh,
//		pnt_config->int_iStartReal,
//		pnt_config->int_iEndReal,
//		pnt_config->int_jStartGhosts,
//		pnt_config->int_jStartReal-1,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	}
//
//
//	if (pnt_config->InterfaceNeighbourTop[0]!=MPI_PROC_NULL)
//	{
//	WriteValuesFromBufferToZD(
//		bufferRecieveTopViscid,
//		pnt_config,
//		pnt_ZD,
//		pnt_mesh,
//		pnt_config->int_iStartReal,
//		pnt_config->int_iEndReal,
//		pnt_config->int_jEndReal+1,
//		pnt_config->int_jEndGhosts,
//		pnt_config->int_kStartReal,
//		pnt_config->int_kEndReal);
//	}
//
//
////#############################################################
////Die Buffer für den Transport in k-Richtung werden geschrieben
////#############################################################
//	if(pnt_config->int_kMeshPoints>1)
//	{
//		WriteValuesFromZDToBuffer(
//			bufferSendBehindViscid,
//			pnt_config,
//			pnt_ZD,
//			pnt_config->int_iStartReal,
//			pnt_config->int_iEndReal,
//			pnt_config->int_jStartReal,
//			pnt_config->int_jEndReal,
//			pnt_config->int_kStartReal,
//			(pnt_config->int_kStartReal+((pnt_config->int_order+1)/2-1)));
//		WriteValuesFromZDToBuffer(
//			bufferSendInFrontViscid,
//			pnt_config,
//			pnt_ZD,
//			pnt_config->int_iStartReal,
//			pnt_config->int_iEndReal,
//			pnt_config->int_jStartReal,
//			pnt_config->int_jEndReal,
//			(pnt_config->int_kEndReal-((pnt_config->int_order+1)/2-1)),
//			pnt_config->int_kEndReal);
//	//	Der Transfer wird durchgeführt. Erst geschieht die Kommunikation nach unten, dann nach oben
//		MPI_Sendrecv(
//			bufferSendInFrontViscid,
//			pnt_config->MPI_intSizeKTransferViscid,
//			MPI_FLT,
//			pnt_config->InterfaceNeighbourInFront[0],
//			20,
//			bufferRecieveBehindViscid,
//			pnt_config->MPI_intSizeKTransferViscid,
//			MPI_FLT,
//			pnt_config->InterfaceNeighbourBehind[0],
//			20,
//			pnt_config->MPI_comm,
//			&pnt_config->MPI_status);
//		MPI_Sendrecv(
//			bufferSendBehindViscid,
//			pnt_config->MPI_intSizeKTransferViscid,
//			MPI_FLT,
//			pnt_config->InterfaceNeighbourBehind[0],
//			21,
//			bufferRecieveInFrontViscid,
//			pnt_config->MPI_intSizeKTransferViscid,
//			MPI_FLT,
//			pnt_config->InterfaceNeighbourInFront[0],
//			21,
//			pnt_config->MPI_comm,
//			&pnt_config->MPI_status);
//		//Die empfangene Buffer wird in die GhostCells geschrieben
//		if (pnt_config->InterfaceNeighbourBehind[0]!=MPI_PROC_NULL)
//		{
//		WriteValuesFromBufferToZD(
//			bufferRecieveBehindViscid,
//			pnt_config,
//			pnt_ZD,
//			pnt_mesh,
//			pnt_config->int_iStartReal,
//			pnt_config->int_iEndReal,
//			pnt_config->int_jStartReal,
//			pnt_config->int_jEndReal,
//			pnt_config->int_kStartGhosts,
//			pnt_config->int_kStartReal-1);
//		}
//
//		if (pnt_config->InterfaceNeighbourInFront[0]!=MPI_PROC_NULL)
//		{
//		WriteValuesFromBufferToZD(
//			bufferRecieveInFrontViscid,
//			pnt_config,
//			pnt_ZD,
//			pnt_mesh,
//			pnt_config->int_iStartReal,
//			pnt_config->int_iEndReal,
//			pnt_config->int_jStartReal,
//			pnt_config->int_jEndReal,
//			pnt_config->int_kEndReal+1,
//			pnt_config->int_kEndGhosts);
//		}
//	}
//}


//void TransferMeshParameter(
//		struct strct_configuration * pnt_config,
//		struct strct_mesh * pnt_mesh)
//{
//	int interface;
//
//	for(interface=0;interface<pnt_config->NumberInterfaces;interface++)
//	{
//		if((pnt_config->MemberOfInterface[interface][0]==pnt_config->MPI_rank)||(pnt_config->MemberOfInterface[interface][1]==pnt_config->MPI_rank))
//		{
//			WriteValuesFromMeshToBuffer(
//				pnt_config->MPI_SendBufferMesh[interface],
//				pnt_config,
//				pnt_mesh,
//				pnt_config->MPI_intIStartSend_WithGhosts[interface],
//				pnt_config->MPI_intIEndSend_WithGhosts[interface],
//				pnt_config->MPI_intJStartSend_WithGhosts[interface],
//				pnt_config->MPI_intJEndSend_WithGhosts[interface],
//				pnt_config->MPI_intKStartSend_WithGhosts[interface],
//				pnt_config->MPI_intKEndSend_WithGhosts[interface]);
//
//			MPI_Sendrecv(
//				pnt_config->MPI_SendBufferMesh[interface],
//				pnt_config->MPI_intTransferSizeMesh[interface],
//				MPI_FLT,
//				pnt_config->MPI_rankNeighbours[interface],
//				interface,
//				pnt_config->MPI_RecieveBufferMesh[interface],
//				pnt_config->MPI_intTransferSizeMesh[interface],
//				MPI_FLT,
//				pnt_config->MPI_rankNeighbours[interface],
//				interface,
//				pnt_config->MPI_comm,
//				&pnt_config->MPI_status);
//
//			WriteValuesFromBufferToMesh(
//				pnt_config->MPI_RecieveBufferMesh[interface],
//				pnt_config,
//				pnt_mesh,
//				pnt_config->MPI_intIStartRecieve_WithGhosts[interface],
//				pnt_config->MPI_intIEndRecieve_WithGhosts[interface],
//				pnt_config->MPI_intJStartRecieve_WithGhosts[interface],
//				pnt_config->MPI_intJEndRecieve_WithGhosts[interface],
//				pnt_config->MPI_intKStartRecieve_WithGhosts[interface],
//				pnt_config->MPI_intKEndRecieve_WithGhosts[interface],
//				interface);
//
//		}
//
//	}
//}

void TransferMeshParameter(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{
	int interface, flag;
	int no_comm=0;
	int counter_comm=0;
	int periodic[2];
	int i_p;
	i_p=0;

	int buffer_copied[pnt_config->NumberInterfaces];
	MPI_Request send_request[pnt_config->NumberInterfaces];
	MPI_Request recieve_request[pnt_config->NumberInterfaces];
	for(interface=0;interface<pnt_config->NumberInterfaces;interface++)
	{
		buffer_copied[interface]=0;
		send_request[interface]=MPI_REQUEST_NULL;
		recieve_request[interface]=MPI_REQUEST_NULL;
	}

	for(interface=0;interface<pnt_config->NumberInterfaces;interface++)
	{
		if (pnt_config->MPI_rank!=pnt_config->MPI_rankNeighbours[interface])
		{
			MPI_Irecv(pnt_config->MPI_RecieveBufferMesh[interface],
					 pnt_config->MPI_intTransferSizeMesh[interface],
				 MPI_FLT,
				 pnt_config->MPI_rankNeighbours[interface],
				 pnt_config->MPI_tag[interface],
				 pnt_config->MPI_comm,
				 &recieve_request[interface]);
		}
		else
		{
			periodic[i_p]=interface;
			i_p++;
		}
	}

	for(interface=0;interface<pnt_config->NumberInterfaces;interface++)
	{
		no_comm++;
		WriteValuesFromMeshToBuffer(
			pnt_config->MPI_SendBufferMesh[interface],
			pnt_config,
			pnt_mesh,
			pnt_config->MPI_intIStartSend_WithGhosts[interface],
			pnt_config->MPI_intIEndSend_WithGhosts[interface],
			pnt_config->MPI_intJStartSend_WithGhosts[interface],
			pnt_config->MPI_intJEndSend_WithGhosts[interface],
			pnt_config->MPI_intKStartSend_WithGhosts[interface],
			pnt_config->MPI_intKEndSend_WithGhosts[interface]);

		if (pnt_config->MPI_rank!=pnt_config->MPI_rankNeighbours[interface])
		{
			MPI_Isend(pnt_config->MPI_SendBufferMesh[interface],
				 pnt_config->MPI_intTransferSizeMesh[interface],
				 MPI_FLT,
				 pnt_config->MPI_rankNeighbours[interface],
				 pnt_config->MPI_tag[interface],
				 pnt_config->MPI_comm,
				 &send_request[interface]);
		}
	}

	if(i_p>0)
	{
		memcpy(pnt_config->MPI_RecieveBufferMesh[periodic[0]],pnt_config->MPI_SendBufferMesh[periodic[1]],pnt_config->MPI_intTransferSizeMesh[periodic[1]]*sizeof(FLT));
		memcpy(pnt_config->MPI_RecieveBufferMesh[periodic[1]],pnt_config->MPI_SendBufferMesh[periodic[0]],pnt_config->MPI_intTransferSizeMesh[periodic[0]]*sizeof(FLT));
	}

	do
	{
		for(interface=0;interface<pnt_config->NumberInterfaces;interface++)
		{
			MPI_Test(&recieve_request[interface],&flag,MPI_STATUSES_IGNORE);
			if ((flag==true)&&(buffer_copied[interface]==0))
			{ //If neighbour in direction i exists
				counter_comm++;
				WriteValuesFromBufferToMesh(
					pnt_config->MPI_RecieveBufferMesh[interface],
					pnt_config,
					pnt_mesh,
					pnt_config->MPI_intIStartRecieve_WithGhosts[interface],
					pnt_config->MPI_intIEndRecieve_WithGhosts[interface],
					pnt_config->MPI_intJStartRecieve_WithGhosts[interface],
					pnt_config->MPI_intJEndRecieve_WithGhosts[interface],
					pnt_config->MPI_intKStartRecieve_WithGhosts[interface],
					pnt_config->MPI_intKEndRecieve_WithGhosts[interface],
					interface);
				buffer_copied[interface]=1;
			}
		}
	}while(counter_comm<no_comm);

	MPI_Waitall( pnt_config->NumberInterfaces, send_request, MPI_STATUSES_IGNORE);
}

void TransferFlowParameterWithGhosts(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	double comm_t0,comm_t1;
	comm_t0 = MPI_Wtime();

	int interface, flag;
	int no_comm=0;
	int counter_comm=0;
	int periodic[2];
	int i_p;
	i_p=0;

	int buffer_copied[pnt_config->NumberInterfaces];
	MPI_Request send_request[pnt_config->NumberInterfaces];
	MPI_Request recieve_request[pnt_config->NumberInterfaces];
	for(interface=0;interface<pnt_config->NumberInterfaces;interface++)
	{
		buffer_copied[interface]=0;
		send_request[interface]=MPI_REQUEST_NULL;
		recieve_request[interface]=MPI_REQUEST_NULL;
	}

	for(interface=0;interface<pnt_config->NumberInterfaces;interface++)
	{
		if (pnt_config->MPI_rank!=pnt_config->MPI_rankNeighbours[interface])
		{
			MPI_Irecv(pnt_config->MPI_RecieveBufferFlowWithGhosts[interface],
					 pnt_config->MPI_intTransferSizeFlow_WithGhosts[interface],
				 MPI_FLT,
				 pnt_config->MPI_rankNeighbours[interface],
				 pnt_config->MPI_tag[interface],
				 pnt_config->MPI_comm,
				 &recieve_request[interface]);
		}
		else
		{
			periodic[i_p]=interface;
			i_p++;
		}
	}

	for(interface=0;interface<pnt_config->NumberInterfaces;interface++)
	{
		no_comm++;
		WriteValuesFromUToBuffer(
			pnt_config->MPI_SendBufferFlowWithGhosts[interface],
			pnt_config,
			pnt_U,
			pnt_config->MPI_intIStartSend_WithGhosts[interface],
			pnt_config->MPI_intIEndSend_WithGhosts[interface],
			pnt_config->MPI_intJStartSend_WithGhosts[interface],
			pnt_config->MPI_intJEndSend_WithGhosts[interface],
			pnt_config->MPI_intKStartSend_WithGhosts[interface],
			pnt_config->MPI_intKEndSend_WithGhosts[interface]);

		if (pnt_config->MPI_rank!=pnt_config->MPI_rankNeighbours[interface])
		{
			MPI_Isend(pnt_config->MPI_SendBufferFlowWithGhosts[interface],
				 pnt_config->MPI_intTransferSizeFlow_WithGhosts[interface],
				 MPI_FLT,
				 pnt_config->MPI_rankNeighbours[interface],
				 pnt_config->MPI_tag[interface],
				 pnt_config->MPI_comm,
				 &send_request[interface]);
		}
	}

	if(i_p>0)
	{
		memcpy(pnt_config->MPI_RecieveBufferFlowWithGhosts[periodic[0]],pnt_config->MPI_SendBufferFlowWithGhosts[periodic[1]],pnt_config->MPI_intTransferSizeFlow_WithGhosts[periodic[1]]*sizeof(FLT));
		memcpy(pnt_config->MPI_RecieveBufferFlowWithGhosts[periodic[1]],pnt_config->MPI_SendBufferFlowWithGhosts[periodic[0]],pnt_config->MPI_intTransferSizeFlow_WithGhosts[periodic[0]]*sizeof(FLT));
	}

	do
	{
		for(interface=0;interface<pnt_config->NumberInterfaces;interface++)
		{
			MPI_Test(&recieve_request[interface],&flag,MPI_STATUSES_IGNORE);
			if ((flag==true)&&(buffer_copied[interface]==0))
			{ //If neighbour in direction i exists
				counter_comm++;
				WriteValuesFromBufferToU(
					pnt_config->MPI_RecieveBufferFlowWithGhosts[interface],
					pnt_config,
					pnt_U,
					pnt_config->MPI_intIStartRecieve_WithGhosts[interface],
					pnt_config->MPI_intIEndRecieve_WithGhosts[interface],
					pnt_config->MPI_intJStartRecieve_WithGhosts[interface],
					pnt_config->MPI_intJEndRecieve_WithGhosts[interface],
					pnt_config->MPI_intKStartRecieve_WithGhosts[interface],
					pnt_config->MPI_intKEndRecieve_WithGhosts[interface],
					interface);
				buffer_copied[interface]=1;
			}
		}
	}while(counter_comm<no_comm);

	MPI_Waitall( pnt_config->NumberInterfaces, send_request, MPI_STATUSES_IGNORE);

	comm_t1 = MPI_Wtime();
	pnt_config->comm_time=pnt_config->comm_time+(comm_t1-comm_t0);
}


void CalcValuesForPost(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
		int i,j,k,ijk,iMinus1jk,iPlus1jk,ijMinus1k,ijPlus1k,ijkMinus1,ijkPlus1;

		for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
		{
			for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
			{
				for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

					iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					iPlus1jk=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

					ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;
					ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;

					ijkMinus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+(k-1);
					ijkPlus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+(k+1);


					pnt_U_lastStep->MachNumber[ijk]= //correct since 27.9.2012
							(
									sqrt(pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk]+
										pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk]+
										pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))/
									sqrt(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk])*pnt_config->machNumber;

					if(MESHDIMENSIONS==2)
					{
						pnt_U_lastStep->gradRho[ijk]=
						sqrt(
						pow((0.5*pnt_U_lastStep->rho[iPlus1jk]-0.5*pnt_U_lastStep->rho[iMinus1jk])/(pnt_config->deltaXi)*pnt_mesh->xi_x[ijk],2.0)+
						pow((0.5*pnt_U_lastStep->rho[ijPlus1k]-0.5*pnt_U_lastStep->rho[ijMinus1k])/(pnt_config->deltaEta)*pnt_mesh->eta_x[ijk],2.0)+
						pow((0.5*pnt_U_lastStep->rho[iPlus1jk]-0.5*pnt_U_lastStep->rho[iMinus1jk])/(pnt_config->deltaXi)*pnt_mesh->xi_y[ijk],2.0)+
						pow((0.5*pnt_U_lastStep->rho[ijPlus1k]-0.5*pnt_U_lastStep->rho[ijMinus1k])/(pnt_config->deltaEta)*pnt_mesh->eta_y[ijk],2.0)
						);
					}
					else
					{
						pnt_U_lastStep->gradRho[ijk]=
						sqrt(
						pow((0.5*pnt_U_lastStep->rho[iPlus1jk]-0.5*pnt_U_lastStep->rho[iMinus1jk])/(pnt_config->deltaXi)*pnt_mesh->xi_x[ijk],2.0)+
						pow((0.5*pnt_U_lastStep->rho[ijPlus1k]-0.5*pnt_U_lastStep->rho[ijMinus1k])/(pnt_config->deltaEta)*pnt_mesh->eta_x[ijk],2.0)+
						pow((0.5*pnt_U_lastStep->rho[ijkPlus1]-0.5*pnt_U_lastStep->rho[ijkMinus1])/(pnt_config->deltaZeta)*pnt_mesh->zeta_x[ijk],2.0)+
						pow((0.5*pnt_U_lastStep->rho[iPlus1jk]-0.5*pnt_U_lastStep->rho[iMinus1jk])/(pnt_config->deltaXi)*pnt_mesh->xi_y[ijk],2.0)+
						pow((0.5*pnt_U_lastStep->rho[ijPlus1k]-0.5*pnt_U_lastStep->rho[ijMinus1k])/(pnt_config->deltaEta)*pnt_mesh->eta_y[ijk],2.0)+
						pow((0.5*pnt_U_lastStep->rho[ijkPlus1]-0.5*pnt_U_lastStep->rho[ijkMinus1])/(pnt_config->deltaZeta)*pnt_mesh->zeta_y[ijk],2.0)+
						pow((0.5*pnt_U_lastStep->rho[iPlus1jk]-0.5*pnt_U_lastStep->rho[iMinus1jk])/(pnt_config->deltaXi)*pnt_mesh->xi_z[ijk],2.0)+
						pow((0.5*pnt_U_lastStep->rho[ijPlus1k]-0.5*pnt_U_lastStep->rho[ijMinus1k])/(pnt_config->deltaEta)*pnt_mesh->eta_z[ijk],2.0)+
						pow((0.5*pnt_U_lastStep->rho[ijkPlus1]-0.5*pnt_U_lastStep->rho[ijkMinus1])/(pnt_config->deltaZeta)*pnt_mesh->zeta_z[ijk],2.0)
						);

						pnt_U_lastStep->Lambda2[ijk]=CalcLambda2(i,j,k,pnt_config,pnt_mesh,pnt_U_lastStep);
					}

					if(pnt_mesh->flag_IBC[ijk]==1)
					{
						if(pnt_config->flag_IBC_Moving==1)
						{pnt_U_lastStep->u[ijk]=(pnt_config->IBC_MovingActualPosition-pnt_config->IBC_MovingLastPosition)/pnt_config->numericalTau;}
						else{pnt_U_lastStep->u[ijk]=0.;}
						pnt_U_lastStep->v[ijk]=0.;
						pnt_U_lastStep->w[ijk]=0.;
						pnt_U_lastStep->gradRho[ijk]=100000.;
						pnt_U_lastStep->MachNumber[ijk]=100000.;
						pnt_U_lastStep->Lambda2[ijk]=100000.;
					}
				}
			}
		}
}

void filecopy(char *fileIn, char *fileOut)
{
	FILE *in, *out;
	char ch;

	if((in=fopen(fileIn, "rb")) == NULL) {
	printf("Cannot open input file.\n");
	exit(1);
	}
	if((out=fopen(fileOut, "wb")) == NULL) {
	printf("Cannot open output file.\n");
	exit(1);
	}

	while(!feof(in)) {
	ch = getc(in);
	if(ferror(in)) {
	  printf("Read Error");
	  clearerr(in);
	  break;
	} else {
	  if(!feof(in)) putc(ch, out);
	  if(ferror(out)) {
		printf("Write Error");
		clearerr(out);
		break;
	  }
	}
	}
	fclose(in);
	fclose(out);
}


void getNeighbour(
		struct strct_configuration * pnt_config)
{
	int i,j,k,l;
	int MPI_tag_factor[6];
	MPI_tag_factor[0]=1;
	MPI_tag_factor[1]=3;
	MPI_tag_factor[2]=5;
	MPI_tag_factor[3]=7;
	MPI_tag_factor[4]=11;
	MPI_tag_factor[5]=13;

	pnt_config->InterfaceNeighbourLeft=NO_NEIGHBOUR;
	pnt_config->InterfaceNeighbourRight=NO_NEIGHBOUR;
	pnt_config->InterfaceNeighbourBottom=NO_NEIGHBOUR;
	pnt_config->InterfaceNeighbourTop=NO_NEIGHBOUR;
	pnt_config->InterfaceNeighbourBehind=NO_NEIGHBOUR;
	pnt_config->InterfaceNeighbourInFront=NO_NEIGHBOUR;


	int dim=pnt_config->int_meshDimensions;
	for(i=0;i<pnt_config->NumberInterfaces;i++)
	{
		j=-1;
		do
		{
			j++;
		}while(strcmp(pnt_config->Donorname[i],pnt_config->ZonenameAll[j])!=0);
//		printf("my rank %d my name %s",pnt_config->MPI_rank,pnt_config->Zonename);
//		printf("donor %s zonenameall %s neighbour %d\n",pnt_config->Donorname[i],pnt_config->ZonenameAll[j],j);

		pnt_config->MPI_rankNeighbours[i]=j;

		if((pnt_config->RangeOfInterface[i][0]==1)&&(pnt_config->RangeOfInterface[i][0+dim]==1))
		{
			pnt_config->InterfaceNeighbourLeft=i;
		}
		if((pnt_config->RangeOfInterface[i][0]!=1)&&(pnt_config->RangeOfInterface[i][0+dim]!=1))
		{
			pnt_config->InterfaceNeighbourRight=i;
		}
		if((pnt_config->RangeOfInterface[i][1]==1)&&(pnt_config->RangeOfInterface[i][1+dim]==1))
		{
			pnt_config->InterfaceNeighbourBottom=i;
		}
		if((pnt_config->RangeOfInterface[i][1]!=1)&&(pnt_config->RangeOfInterface[i][1+dim]!=1))
		{
			pnt_config->InterfaceNeighbourTop=i;
		}
		if(dim==3)
		{
			if((pnt_config->RangeOfInterface[i][2]==1)&&(pnt_config->RangeOfInterface[i][2+dim]==1))
			{
				pnt_config->InterfaceNeighbourBehind=i;
			}
			if((pnt_config->RangeOfInterface[i][2]!=1)&&(pnt_config->RangeOfInterface[i][2+dim]!=1))
			{
				pnt_config->InterfaceNeighbourInFront=i;
			}
		}

		int flag;
		void *v;
		int value_MPI_TAG_UB;
		MPI_Comm_get_attr( MPI_COMM_WORLD, MPI_TAG_UB, &v, &flag );
		value_MPI_TAG_UB = *(int*)v;
//		printf("MAX TAG: %d\n",value_MPI_TAG_UB);

		if(pnt_config->MPI_rank<=pnt_config->MPI_rankNeighbours[i])
		{
			pnt_config->MPI_tag[i]=0;
			for(k=0;k<2*pnt_config->int_meshDimensions;k++)
	        {
				pnt_config->MPI_tag[i]+=(pnt_config->RangeOfInterface[i][k]*MPI_tag_factor[k]);
	        }
			while(pnt_config->MPI_tag[i]>value_MPI_TAG_UB)
			{
				pnt_config->MPI_tag[i]=pnt_config->MPI_tag[i]-value_MPI_TAG_UB;
			}
		}
		else
		{
			pnt_config->MPI_tag[i]=0;
			for(k=0;k<2*pnt_config->int_meshDimensions;k++)
	        {
				pnt_config->MPI_tag[i]+=(pnt_config->DonorRangeOfInterface[i][k]*MPI_tag_factor[k]);
	        }
			while(pnt_config->MPI_tag[i]>value_MPI_TAG_UB)
			{
				pnt_config->MPI_tag[i]=pnt_config->MPI_tag[i]-value_MPI_TAG_UB;
			}
		}
		for(l=0;l<i;l++)
		{
			if ((pnt_config->MPI_tag[i]==pnt_config->MPI_tag[l])&&(pnt_config->MPI_rankNeighbours[i]==pnt_config->MPI_rankNeighbours[l]))
			{
				printf("ERROR: Zwei Interfaces mit identischer tag entdeckt:");
				printf("tag1/2: rank %d - partner %d - TAG: %ld\n",pnt_config->MPI_rank,pnt_config->MPI_rankNeighbours[l],pnt_config->MPI_tag[l]);
			}
		}


	}
}

void check_TransformationMatrix(
		int interface,
		struct strct_configuration * pnt_config)
{

//	Anhand der Pointrange und der DonorPointrange wird das Vorzeichen der Laufkoordinate (i,j,k) bestimmt.
//	Dieses ist durch Pointwise nicht physikalisch richtig in der Transformmatrix.
//	Zunächst werden die Punkte gesucht, die am Interface konstant sind. Sind diese beim Geber 1 und beim mir ungleich 1
//	ist die laufrichtung konstant, so dass sich ein Vorzeichen von insgesamt plus ergibt.
//	Wenn Sie bei mir 1 sind müssen Sie beim geber ungleich 1 sein damit es insgesamt positiv wird
//	Die anderen Fälle ergeben ein negatives Vorzeichen



	int vorzeichen_me,vorzeichen_donor,vorzeichen_gesamt;
	vorzeichen_donor=1;
	vorzeichen_me=1;
//	Ueberpruefung des Vorzeichen fuer die eigene Seite
	if (pnt_config->int_meshDimensions==2)
	{
		if(pnt_config->RangeOfInterface[interface][0]==pnt_config->RangeOfInterface[interface][2])
		{
			if(pnt_config->RangeOfInterface[interface][0]!=1){vorzeichen_me=1;}
			else{vorzeichen_me=-1;}
		}
		else
		{
			if(pnt_config->RangeOfInterface[interface][1]!=1){vorzeichen_me=1;}
			else{vorzeichen_me=-1;}
		}
	}

	if (pnt_config->int_meshDimensions==3)
	{
		if(pnt_config->RangeOfInterface[interface][0]==pnt_config->RangeOfInterface[interface][3])
		{
			if(pnt_config->RangeOfInterface[interface][0]!=1){vorzeichen_me=1;}
			else{vorzeichen_me=-1;}
		}
		else if(pnt_config->RangeOfInterface[interface][1]==pnt_config->RangeOfInterface[interface][4])
		{
			if(pnt_config->RangeOfInterface[interface][1]!=1){vorzeichen_me=1;}
			else{vorzeichen_me=-1;}
		}
		else
		{
			if(pnt_config->RangeOfInterface[interface][2]!=1){vorzeichen_me=1;}
			else{vorzeichen_me=-1;}
		}
	}

//	Ueberpruefung des Vorzeichen fuer die Donorseite
	if (pnt_config->int_meshDimensions==2)
	{
		if(pnt_config->DonorRangeOfInterface[interface][0]==pnt_config->DonorRangeOfInterface[interface][2])
		{
			if(pnt_config->DonorRangeOfInterface[interface][0]==1){vorzeichen_donor=1;}
			else{vorzeichen_donor=-1;}
		}
		else
		{
			if(pnt_config->DonorRangeOfInterface[interface][1]==1){vorzeichen_donor=1;}
			else{vorzeichen_donor=-1;}
		}
	}

	if (pnt_config->int_meshDimensions==3)
	{
		if(pnt_config->DonorRangeOfInterface[interface][0]==pnt_config->DonorRangeOfInterface[interface][3])
		{
			if(pnt_config->DonorRangeOfInterface[interface][0]==1){vorzeichen_donor=1;}
			else{vorzeichen_donor=-1;}
		}
		else if(pnt_config->DonorRangeOfInterface[interface][1]==pnt_config->DonorRangeOfInterface[interface][4])
		{
			if(pnt_config->DonorRangeOfInterface[interface][1]==1){vorzeichen_donor=1;}
			else{vorzeichen_donor=-1;}
		}
		else
		{
			if(pnt_config->DonorRangeOfInterface[interface][2]==1){vorzeichen_donor=1;}
			else{vorzeichen_donor=-1;}
		}
	}

	vorzeichen_gesamt=vorzeichen_me*vorzeichen_donor;

	if((interface==pnt_config->InterfaceNeighbourLeft)||(interface==pnt_config->InterfaceNeighbourRight))
	{
		pnt_config->TransformMatrixOfInterface[interface][0]=
				vorzeichen_gesamt*abs(pnt_config->TransformMatrixOfInterface[interface][0]);
	}

	if((interface==pnt_config->InterfaceNeighbourBottom)||(interface==pnt_config->InterfaceNeighbourTop))
	{
		pnt_config->TransformMatrixOfInterface[interface][1]=
				vorzeichen_gesamt*abs(pnt_config->TransformMatrixOfInterface[interface][1]);
	}

	if((interface==pnt_config->InterfaceNeighbourBehind)||(interface==pnt_config->InterfaceNeighbourInFront))
	{
		pnt_config->TransformMatrixOfInterface[interface][2]=
				vorzeichen_gesamt*abs(pnt_config->TransformMatrixOfInterface[interface][2]);
	}

	pnt_config->MPI_intTransformation_flag_I0_I[interface]= 0;
	pnt_config->MPI_intTransformation_flag_I0_J[interface]= 0;
	pnt_config->MPI_intTransformation_flag_I0_K[interface]= 0;
	pnt_config->MPI_intTransformation_flag_J0_I[interface]= 0;
	pnt_config->MPI_intTransformation_flag_J0_J[interface]= 0;
	pnt_config->MPI_intTransformation_flag_J0_K[interface]= 0;
	pnt_config->MPI_intTransformation_flag_K0_I[interface]= 0;
	pnt_config->MPI_intTransformation_flag_K0_J[interface]= 0;
	pnt_config->MPI_intTransformation_flag_K0_K[interface]= 0;
	pnt_config->MPI_intTransformation_Offset_I[interface]= 0;
	pnt_config->MPI_intTransformation_Offset_J[interface]= 0;
	pnt_config->MPI_intTransformation_Offset_K[interface]= 0;

	switch(pnt_config->TransformMatrixOfInterface[interface][0])
	{
		case 1:
			pnt_config->MPI_intTransformation_IMax[interface]=
					pnt_config->MPI_intIEndRecieve[interface]-pnt_config->MPI_intIStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_IMax_Mesh[interface]=
					pnt_config->MPI_intIEndRecieve_WithGhosts[interface]-pnt_config->MPI_intIStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_flag_I0_I[interface]= 1;

			pnt_config->MPI_dblTransformation_xi_x[interface]=1.0;
			pnt_config->MPI_dblTransformation_xi_y[interface]=1.0;
			pnt_config->MPI_dblTransformation_xi_z[interface]=1.0;
			break;
		case 2:
			pnt_config->MPI_intTransformation_JMax[interface]=
					pnt_config->MPI_intIEndRecieve[interface]-pnt_config->MPI_intIStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_JMax_Mesh[interface]=
					pnt_config->MPI_intIEndRecieve_WithGhosts[interface]-pnt_config->MPI_intIStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_flag_I0_J[interface]= 1;

			pnt_config->MPI_dblTransformation_xi_x[interface]=1.0;
			pnt_config->MPI_dblTransformation_xi_y[interface]=1.0;
			pnt_config->MPI_dblTransformation_xi_z[interface]=1.0;
			break;
		case 3:
			pnt_config->MPI_intTransformation_KMax[interface]=
					pnt_config->MPI_intIEndRecieve[interface]-pnt_config->MPI_intIStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_KMax_Mesh[interface]=
					pnt_config->MPI_intIEndRecieve_WithGhosts[interface]-pnt_config->MPI_intIStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_flag_I0_K[interface]= 1;

			pnt_config->MPI_dblTransformation_xi_x[interface]=1.0;
			pnt_config->MPI_dblTransformation_xi_y[interface]=1.0;
			pnt_config->MPI_dblTransformation_xi_z[interface]=1.0;
			break;
		case -1:
			pnt_config->MPI_intTransformation_IMax[interface]=
					pnt_config->MPI_intIEndRecieve[interface]-pnt_config->MPI_intIStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_flag_I0_I[interface]= 1;
			pnt_config->MPI_intTransformation_Offset_I[interface]=
					pnt_config->MPI_intTransformation_IMax[interface]-1;

			pnt_config->MPI_intTransformation_IMax_Mesh[interface]=
					pnt_config->MPI_intIEndRecieve_WithGhosts[interface]-pnt_config->MPI_intIStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_Offset_I_Ghosts[interface]=
					pnt_config->MPI_intTransformation_IMax_Mesh[interface]-1;

			pnt_config->MPI_dblTransformation_xi_x[interface]=-1.0;
			pnt_config->MPI_dblTransformation_xi_y[interface]=-1.0;
			pnt_config->MPI_dblTransformation_xi_z[interface]=-1.0;
			break;
		case -2:
			pnt_config->MPI_intTransformation_JMax[interface]=
					pnt_config->MPI_intIEndRecieve[interface]-pnt_config->MPI_intIStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_flag_I0_J[interface]=1;
			pnt_config->MPI_intTransformation_Offset_J[interface]=
					pnt_config->MPI_intTransformation_JMax[interface]-1;

			pnt_config->MPI_intTransformation_JMax_Mesh[interface]=
					pnt_config->MPI_intIEndRecieve_WithGhosts[interface]-pnt_config->MPI_intIStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_Offset_J_Ghosts[interface]=
					pnt_config->MPI_intTransformation_JMax_Mesh[interface]-1;

			pnt_config->MPI_dblTransformation_xi_x[interface]=-1.0;
			pnt_config->MPI_dblTransformation_xi_y[interface]=-1.0;
			pnt_config->MPI_dblTransformation_xi_z[interface]=-1.0;
			break;
		case -3:
			pnt_config->MPI_intTransformation_KMax[interface]=
					pnt_config->MPI_intIEndRecieve[interface]-pnt_config->MPI_intIStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_flag_I0_K[interface]= 1;
			pnt_config->MPI_intTransformation_Offset_K[interface]=
					pnt_config->MPI_intTransformation_KMax[interface]-1;

			pnt_config->MPI_intTransformation_KMax_Mesh[interface]=
					pnt_config->MPI_intIEndRecieve_WithGhosts[interface]-pnt_config->MPI_intIStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_Offset_K_Ghosts[interface]=
					pnt_config->MPI_intTransformation_KMax_Mesh[interface]-1;

			pnt_config->MPI_dblTransformation_xi_x[interface]=-1.0;
			pnt_config->MPI_dblTransformation_xi_y[interface]=-1.0;
			pnt_config->MPI_dblTransformation_xi_z[interface]=-1.0;
			break;
	}

	switch(pnt_config->TransformMatrixOfInterface[interface][1])
	{
		case 1:
			pnt_config->MPI_intTransformation_IMax[interface]=
					pnt_config->MPI_intJEndRecieve[interface]-pnt_config->MPI_intJStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_IMax_Mesh[interface]=
					pnt_config->MPI_intJEndRecieve_WithGhosts[interface]-pnt_config->MPI_intJStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_flag_J0_I[interface]= 1;

			pnt_config->MPI_dblTransformation_eta_x[interface]=1.0;
			pnt_config->MPI_dblTransformation_eta_y[interface]=1.0;
			pnt_config->MPI_dblTransformation_eta_z[interface]=1.0;
			break;
		case 2:
			pnt_config->MPI_intTransformation_JMax[interface]=
					pnt_config->MPI_intJEndRecieve[interface]-pnt_config->MPI_intJStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_JMax_Mesh[interface]=
					pnt_config->MPI_intJEndRecieve_WithGhosts[interface]-pnt_config->MPI_intJStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_flag_J0_J[interface]= 1;

			pnt_config->MPI_dblTransformation_eta_x[interface]=1.0;
			pnt_config->MPI_dblTransformation_eta_y[interface]=1.0;
			pnt_config->MPI_dblTransformation_eta_z[interface]=1.0;
			break;
		case 3:
			pnt_config->MPI_intTransformation_KMax[interface]=
					pnt_config->MPI_intJEndRecieve[interface]-pnt_config->MPI_intJStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_KMax_Mesh[interface]=
					pnt_config->MPI_intJEndRecieve_WithGhosts[interface]-pnt_config->MPI_intJStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_flag_J0_K[interface]= 1;

			pnt_config->MPI_dblTransformation_eta_x[interface]=1.0;
			pnt_config->MPI_dblTransformation_eta_y[interface]=1.0;
			pnt_config->MPI_dblTransformation_eta_z[interface]=1.0;
			break;
		case -1:
			pnt_config->MPI_intTransformation_IMax[interface]=
					pnt_config->MPI_intJEndRecieve[interface]-pnt_config->MPI_intJStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_flag_J0_I[interface]= 1;
			pnt_config->MPI_intTransformation_Offset_I[interface]=
					pnt_config->MPI_intTransformation_IMax[interface]-1;

			pnt_config->MPI_intTransformation_IMax_Mesh[interface]=
					pnt_config->MPI_intJEndRecieve_WithGhosts[interface]-pnt_config->MPI_intJStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_Offset_I_Ghosts[interface]=
					pnt_config->MPI_intTransformation_IMax_Mesh[interface]-1;

			pnt_config->MPI_dblTransformation_eta_x[interface]=-1.0;
			pnt_config->MPI_dblTransformation_eta_y[interface]=-1.0;
			pnt_config->MPI_dblTransformation_eta_z[interface]=-1.0;
			break;
		case -2:
			pnt_config->MPI_intTransformation_JMax[interface]=
					pnt_config->MPI_intJEndRecieve[interface]-pnt_config->MPI_intJStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_flag_J0_J[interface]= 1;
			pnt_config->MPI_intTransformation_Offset_J[interface]=
					pnt_config->MPI_intTransformation_JMax[interface]-1;

			pnt_config->MPI_intTransformation_JMax_Mesh[interface]=
					pnt_config->MPI_intJEndRecieve_WithGhosts[interface]-pnt_config->MPI_intJStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_Offset_J_Ghosts[interface]=
					pnt_config->MPI_intTransformation_JMax_Mesh[interface]-1;

			pnt_config->MPI_dblTransformation_eta_x[interface]=-1.0;
			pnt_config->MPI_dblTransformation_eta_y[interface]=-1.0;
			pnt_config->MPI_dblTransformation_eta_z[interface]=-1.0;
			break;
		case -3:
			pnt_config->MPI_intTransformation_KMax[interface]=
					pnt_config->MPI_intJEndRecieve[interface]-pnt_config->MPI_intJStartRecieve[interface]+1;
			pnt_config->MPI_intTransformation_flag_J0_K[interface]= 1;
			pnt_config->MPI_intTransformation_Offset_K[interface]=
					pnt_config->MPI_intTransformation_KMax[interface]-1;

			pnt_config->MPI_intTransformation_KMax_Mesh[interface]=
					pnt_config->MPI_intJEndRecieve_WithGhosts[interface]-pnt_config->MPI_intJStartRecieve_WithGhosts[interface]+1;
			pnt_config->MPI_intTransformation_Offset_K_Ghosts[interface]=
					pnt_config->MPI_intTransformation_KMax_Mesh[interface]-1;

			pnt_config->MPI_dblTransformation_eta_x[interface]=-1.0;
			pnt_config->MPI_dblTransformation_eta_y[interface]=-1.0;
			pnt_config->MPI_dblTransformation_eta_z[interface]=-1.0;
			break;
	}

	if (pnt_config->int_meshDimensions==3)
	{
		switch(pnt_config->TransformMatrixOfInterface[interface][2])
		{
			case 1:
				pnt_config->MPI_intTransformation_IMax[interface]=
						pnt_config->MPI_intKEndRecieve[interface]-pnt_config->MPI_intKStartRecieve[interface]+1;
				pnt_config->MPI_intTransformation_IMax_Mesh[interface]=
						pnt_config->MPI_intKEndRecieve_WithGhosts[interface]-pnt_config->MPI_intKStartRecieve_WithGhosts[interface]+1;
				pnt_config->MPI_intTransformation_flag_K0_I[interface]= 1;

				pnt_config->MPI_dblTransformation_zeta_x[interface]=1.0;
				pnt_config->MPI_dblTransformation_zeta_y[interface]=1.0;
				pnt_config->MPI_dblTransformation_zeta_z[interface]=1.0;
				break;
			case 2:
				pnt_config->MPI_intTransformation_JMax[interface]=
						pnt_config->MPI_intKEndRecieve[interface]-pnt_config->MPI_intKStartRecieve[interface]+1;
				pnt_config->MPI_intTransformation_JMax_Mesh[interface]=
						pnt_config->MPI_intKEndRecieve_WithGhosts[interface]-pnt_config->MPI_intKStartRecieve_WithGhosts[interface]+1;
				pnt_config->MPI_intTransformation_flag_K0_J[interface]= 1;

				pnt_config->MPI_dblTransformation_zeta_x[interface]=1.0;
				pnt_config->MPI_dblTransformation_zeta_y[interface]=1.0;
				pnt_config->MPI_dblTransformation_zeta_z[interface]=1.0;
				break;
			case 3:
				pnt_config->MPI_intTransformation_KMax[interface]=
						pnt_config->MPI_intKEndRecieve[interface]-pnt_config->MPI_intKStartRecieve[interface]+1;
				pnt_config->MPI_intTransformation_KMax_Mesh[interface]=
						pnt_config->MPI_intKEndRecieve_WithGhosts[interface]-pnt_config->MPI_intKStartRecieve_WithGhosts[interface]+1;
				pnt_config->MPI_intTransformation_flag_K0_K[interface]= 1;

				pnt_config->MPI_dblTransformation_zeta_x[interface]=1.0;
				pnt_config->MPI_dblTransformation_zeta_y[interface]=1.0;
				pnt_config->MPI_dblTransformation_zeta_z[interface]=1.0;
				break;
			case -1:
				pnt_config->MPI_intTransformation_IMax[interface]=
						pnt_config->MPI_intKEndRecieve[interface]-pnt_config->MPI_intKStartRecieve[interface]+1;
				pnt_config->MPI_intTransformation_flag_K0_I[interface]= 1;
				pnt_config->MPI_intTransformation_Offset_I[interface]=
						pnt_config->MPI_intTransformation_IMax[interface]-1;

				pnt_config->MPI_intTransformation_IMax_Mesh[interface]=
						pnt_config->MPI_intKEndRecieve_WithGhosts[interface]-pnt_config->MPI_intKStartRecieve_WithGhosts[interface]+1;
				pnt_config->MPI_intTransformation_Offset_I_Ghosts[interface]=
						pnt_config->MPI_intTransformation_IMax_Mesh[interface]-1;

				pnt_config->MPI_dblTransformation_zeta_x[interface]=-1.0;
				pnt_config->MPI_dblTransformation_zeta_y[interface]=-1.0;
				pnt_config->MPI_dblTransformation_zeta_z[interface]=-1.0;
				break;
			case -2:
				pnt_config->MPI_intTransformation_JMax[interface]=
						pnt_config->MPI_intKEndRecieve[interface]-pnt_config->MPI_intKStartRecieve[interface]+1;
				pnt_config->MPI_intTransformation_flag_K0_J[interface]= 1;
				pnt_config->MPI_intTransformation_Offset_J[interface]=
						pnt_config->MPI_intTransformation_JMax[interface]-1;

				pnt_config->MPI_intTransformation_JMax_Mesh[interface]=
						pnt_config->MPI_intKEndRecieve_WithGhosts[interface]-pnt_config->MPI_intKStartRecieve_WithGhosts[interface]+1;
				pnt_config->MPI_intTransformation_Offset_J_Ghosts[interface]=
						pnt_config->MPI_intTransformation_JMax_Mesh[interface]-1;

				pnt_config->MPI_dblTransformation_zeta_x[interface]=-1.0;
				pnt_config->MPI_dblTransformation_zeta_y[interface]=-1.0;
				pnt_config->MPI_dblTransformation_zeta_z[interface]=-1.0;
				break;
			case -3:
				pnt_config->MPI_intTransformation_KMax[interface]=
						pnt_config->MPI_intKEndRecieve[interface]-pnt_config->MPI_intKStartRecieve[interface]+1;
				pnt_config->MPI_intTransformation_flag_K0_K[interface]= 1;
				pnt_config->MPI_intTransformation_Offset_K[interface]=
						pnt_config->MPI_intTransformation_KMax[interface]-1;

				pnt_config->MPI_intTransformation_KMax_Mesh[interface]=
						pnt_config->MPI_intKEndRecieve_WithGhosts[interface]-pnt_config->MPI_intKStartRecieve_WithGhosts[interface]+1;
				pnt_config->MPI_intTransformation_Offset_K_Ghosts[interface]=
						pnt_config->MPI_intTransformation_KMax_Mesh[interface]-1;

				pnt_config->MPI_dblTransformation_zeta_x[interface]=-1.0;
				pnt_config->MPI_dblTransformation_zeta_y[interface]=-1.0;
				pnt_config->MPI_dblTransformation_zeta_z[interface]=-1.0;
				break;
		}
	}
	else
	{
		pnt_config->MPI_intTransformation_KMax[interface]=1;
		pnt_config->MPI_intTransformation_KMax_Mesh[interface]=1;
		pnt_config->MPI_dblTransformation_zeta_x[interface]=1.0;
		pnt_config->MPI_dblTransformation_zeta_y[interface]=1.0;
		pnt_config->MPI_dblTransformation_zeta_z[interface]=1.0;
	}
}

int get_ijkTransformWithGhosts(
		struct strct_configuration * pnt_config,
		int interface,
		int i,
		int j,
		int k,
		int c)
{
	int i_buffer,j_buffer,k_buffer;
	i_buffer=pnt_config->MPI_intTransformation_flag_I0_I[interface]*i+pnt_config->MPI_intTransformation_flag_J0_I[interface]*j+pnt_config->MPI_intTransformation_flag_K0_I[interface]*k;
	j_buffer=pnt_config->MPI_intTransformation_flag_I0_J[interface]*i+pnt_config->MPI_intTransformation_flag_J0_J[interface]*j+pnt_config->MPI_intTransformation_flag_K0_J[interface]*k;
	k_buffer=pnt_config->MPI_intTransformation_flag_I0_K[interface]*i+pnt_config->MPI_intTransformation_flag_J0_K[interface]*j+pnt_config->MPI_intTransformation_flag_K0_K[interface]*k;

	int result;

//	result=abs(i_buffer-pnt_config->MPI_intTransformation_Offset_I[interface])+
//			abs(j_buffer-pnt_config->MPI_intTransformation_Offset_J[interface])*pnt_config->MPI_intTransformation_IMax[interface]+
//			abs(k_buffer-pnt_config->MPI_intTransformation_Offset_K[interface])*pnt_config->MPI_intTransformation_IMax[interface]*pnt_config->MPI_intTransformation_JMax[interface]+
//			c*pnt_config->MPI_intTransformation_IMax[interface]*pnt_config->MPI_intTransformation_JMax[interface]*pnt_config->MPI_intTransformation_KMax[interface];

	result=abs(i_buffer-pnt_config->MPI_intTransformation_Offset_I_Ghosts[interface])+
			abs(j_buffer-pnt_config->MPI_intTransformation_Offset_J_Ghosts[interface])*pnt_config->MPI_intTransformation_IMax_Mesh[interface]+
			abs(k_buffer-pnt_config->MPI_intTransformation_Offset_K_Ghosts[interface])*pnt_config->MPI_intTransformation_IMax_Mesh[interface]*pnt_config->MPI_intTransformation_JMax_Mesh[interface]+
			c*pnt_config->MPI_intTransformation_IMax_Mesh[interface]*pnt_config->MPI_intTransformation_JMax_Mesh[interface]*pnt_config->MPI_intTransformation_KMax_Mesh[interface];

	return result;
}

int get_ijkTransformMesh(
		struct strct_configuration * pnt_config,
		int interface,
		int i,
		int j,
		int k,
		int c)
{
	int i_buffer,j_buffer,k_buffer;
	i_buffer=pnt_config->MPI_intTransformation_flag_I0_I[interface]*i+pnt_config->MPI_intTransformation_flag_J0_I[interface]*j+pnt_config->MPI_intTransformation_flag_K0_I[interface]*k;
	j_buffer=pnt_config->MPI_intTransformation_flag_I0_J[interface]*i+pnt_config->MPI_intTransformation_flag_J0_J[interface]*j+pnt_config->MPI_intTransformation_flag_K0_J[interface]*k;
	k_buffer=pnt_config->MPI_intTransformation_flag_I0_K[interface]*i+pnt_config->MPI_intTransformation_flag_J0_K[interface]*j+pnt_config->MPI_intTransformation_flag_K0_K[interface]*k;

	if ((c>=3)&&(c<=5))
	{
		switch(pnt_config->TransformMatrixOfInterface[interface][0])
		{
			case 1:
				break;
			case 2:
				c=c+3;
				break;
			case 3:
				c=c+6;
				break;
			case -1:
				break;
			case -2:
				c=c+3;
				break;
			case -3:
				c=c+6;
				break;
		}
	}
	else if ((c>=6)&&(c<=8))
	{
		switch(pnt_config->TransformMatrixOfInterface[interface][1])
		{
			case 1:
				c=c-3;
				break;
			case 2:
				break;
			case 3:
				c=c+3;
				break;
			case -1:
				c=c-3;
				break;
			case -2:
				break;
			case -3:
				c=c+3;
				break;
		}
	}
	else if ((c>=9)&&(c<=11))
	{
		if (pnt_config->int_meshDimensions==3)
		{
			switch(pnt_config->TransformMatrixOfInterface[interface][2])
			{
				case 1:
					c=c-6;
					break;
				case 2:
					c=c-3;
					break;
				case 3:
					break;
				case -1:
					c=c-6;
					break;
				case -2:
					c=c-3;
					break;
				case -3:
					break;
			}
		}
	}


	int result;

	result=abs(i_buffer-pnt_config->MPI_intTransformation_Offset_I_Ghosts[interface])+
			abs(j_buffer-pnt_config->MPI_intTransformation_Offset_J_Ghosts[interface])*pnt_config->MPI_intTransformation_IMax_Mesh[interface]+
			abs(k_buffer-pnt_config->MPI_intTransformation_Offset_K_Ghosts[interface])*pnt_config->MPI_intTransformation_IMax_Mesh[interface]*pnt_config->MPI_intTransformation_JMax_Mesh[interface]+
			c*pnt_config->MPI_intTransformation_IMax_Mesh[interface]*pnt_config->MPI_intTransformation_JMax_Mesh[interface]*pnt_config->MPI_intTransformation_KMax_Mesh[interface];

	return result;
}


//Sofern die Randpunkte eines Interfaces übereinanderliegen (Abstand kleiner als 1e-6) wird für einen Prozessor das Gebiet verkleinert.
//Dafür wird der eigene Rank mit dem Rank des Partners verlglichen und der eigene rank größer ist
//werden die zu sendenden Gitterpunkte und das Rechengebiet selber angepasst
void check_Connectivity(
		struct strct_configuration * pnt_config,
		struct strct_mesh* pnt_mesh)
{
	int interface;
	int ijk_real,ijk_ghost;
	int ghost;
	int i,j,k;
	FLT difference;
	for(interface=0;interface<pnt_config->NumberInterfaces;interface++)
	{
		if(interface==pnt_config->InterfaceNeighbourLeft)
		{
			i=pnt_config->int_iStartReal;
			j=pnt_config->int_jMid;
			k=pnt_config->int_kMid;
			ghost=pnt_config->int_iStartReal-1;
			ijk_real=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
			ijk_ghost=ghost*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
			difference=
					fabs(pnt_mesh->x[ijk_real]-pnt_mesh->x[ijk_ghost])+
					fabs(pnt_mesh->y[ijk_real]-pnt_mesh->y[ijk_ghost])+
					fabs(pnt_mesh->z[ijk_real]-pnt_mesh->z[ijk_ghost]);


			if(difference<1e-6)
			{
				pnt_config->MPI_intIStartSend[interface]=pnt_config->MPI_intIStartSend[interface]+1;
				pnt_config->MPI_intIEndSend[interface]=pnt_config->MPI_intIEndSend[interface]+1;
			}
		}
		if(interface==pnt_config->InterfaceNeighbourRight)
		{
			i=pnt_config->int_iEndReal;
			j=pnt_config->int_jMid;
			k=pnt_config->int_kMid;
			ghost=pnt_config->int_iEndReal+1;
			ijk_real=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
			ijk_ghost=ghost*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
			difference=
					fabs(pnt_mesh->x[ijk_real]-pnt_mesh->x[ijk_ghost])+
					fabs(pnt_mesh->y[ijk_real]-pnt_mesh->y[ijk_ghost])+
					fabs(pnt_mesh->z[ijk_real]-pnt_mesh->z[ijk_ghost]);

			if(difference<1e-6)
			{
				pnt_config->MPI_intIStartSend[interface]=pnt_config->MPI_intIStartSend[interface]-1;
				pnt_config->MPI_intIEndSend[interface]=pnt_config->MPI_intIEndSend[interface]-1;
			}
		}
		if(interface==pnt_config->InterfaceNeighbourBottom)
		{
			j=pnt_config->int_jStartReal;
			i=pnt_config->int_iMid;
			k=pnt_config->int_kMid;
			ghost=pnt_config->int_jStartReal-1;
			ijk_real=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
			ijk_ghost=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+ghost*pnt_config->int_kMeshPointsGhostCells+k;
			difference=
					fabs(pnt_mesh->x[ijk_real]-pnt_mesh->x[ijk_ghost])+
					fabs(pnt_mesh->y[ijk_real]-pnt_mesh->y[ijk_ghost])+
					fabs(pnt_mesh->z[ijk_real]-pnt_mesh->z[ijk_ghost]);

			if(difference<1e-6)
			{
				pnt_config->MPI_intJStartSend[interface]=pnt_config->MPI_intJStartSend[interface]+1;
				pnt_config->MPI_intJEndSend[interface]=pnt_config->MPI_intJEndSend[interface]+1;
			}

		}
		if(interface==pnt_config->InterfaceNeighbourTop)
		{
			j=pnt_config->int_jEndReal;
			i=pnt_config->int_iMid;
			k=pnt_config->int_kMid;
			ghost=pnt_config->int_jEndReal+1;
			ijk_real=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
			ijk_ghost=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+ghost*pnt_config->int_kMeshPointsGhostCells+k;
			difference=
					fabs(pnt_mesh->x[ijk_real]-pnt_mesh->x[ijk_ghost])+
					fabs(pnt_mesh->y[ijk_real]-pnt_mesh->y[ijk_ghost])+
					fabs(pnt_mesh->z[ijk_real]-pnt_mesh->z[ijk_ghost]);

			if(difference<1e-6)
			{
				pnt_config->MPI_intJStartSend[interface]=pnt_config->MPI_intJStartSend[interface]-1;
				pnt_config->MPI_intJEndSend[interface]=pnt_config->MPI_intJEndSend[interface]-1;
			}

		}
		if(MESHDIMENSIONS==3)
		{
			if(interface==pnt_config->InterfaceNeighbourBehind)
			{
				k=pnt_config->int_kStartReal;
				i=pnt_config->int_iMid;
				j=pnt_config->int_jMid;
				ghost=pnt_config->int_kStartReal-1;
				ijk_real=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijk_ghost=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+ghost;
				difference=
						fabs(pnt_mesh->x[ijk_real]-pnt_mesh->x[ijk_ghost])+
						fabs(pnt_mesh->y[ijk_real]-pnt_mesh->y[ijk_ghost])+
						fabs(pnt_mesh->z[ijk_real]-pnt_mesh->z[ijk_ghost]);

				if(difference<1e-6)
				{
					pnt_config->MPI_intKStartSend[interface]=pnt_config->MPI_intKStartSend[interface]+1;
					pnt_config->MPI_intKEndSend[interface]=pnt_config->MPI_intKEndSend[interface]+1;
				}

			}
			if(interface==pnt_config->InterfaceNeighbourInFront)
			{
				k=pnt_config->int_kEndReal;
				i=pnt_config->int_iMid;
				j=pnt_config->int_jMid;
				ghost=pnt_config->int_kEndReal+1;
				ijk_real=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijk_ghost=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+ghost;
				difference=
						fabs(pnt_mesh->x[ijk_real]-pnt_mesh->x[ijk_ghost])+
						fabs(pnt_mesh->y[ijk_real]-pnt_mesh->y[ijk_ghost])+
						fabs(pnt_mesh->z[ijk_real]-pnt_mesh->z[ijk_ghost]);

				if(difference<1e-6)
				{
					pnt_config->MPI_intKStartSend[interface]=pnt_config->MPI_intKStartSend[interface]-1;
					pnt_config->MPI_intKEndSend[interface]=pnt_config->MPI_intKEndSend[interface]-1;
				}

			}
		}

	}
}


int check_CGNSFile(
		struct strct_configuration * pnt_config)
{

	FILE *datei;

	datei=fopen(pnt_config->chr_MeshPath,"r");
	if(datei == NULL)
	{

		if(pnt_config->MPI_rank==0){printf("SHOCK: CGNS-Datei '%s' wurde nicht gefunden.\n",pnt_config->chr_MeshPath);}
		return 0;
	}
	else
	{
		return 1;
	}
}

int check_ConfigFile(
		struct strct_configuration * pnt_config)
{

	FILE *datei;

	datei=fopen(pnt_config->chr_configPath,"r");
	if(datei == NULL)
	{

		if(pnt_config->MPI_rank==0){printf("SHOCK: Config-Datei '%s' wurde nicht gefunden.\n",pnt_config->chr_configPath);}
		return 0;
	}
	else
	{
		return 1;
	}
}

void IBC_Actual2Last(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{
	for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				pnt_mesh->flag_IBC_last[ijk]=pnt_mesh->flag_IBC[ijk];
			}
		}
	}
}

void IBC_prepare(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{
	//####################
	//	Custom
	//####################
	FLT x_min,x_max;
	FLT y_min,y_max;
	FLT z_min,z_max;

	//####################
	//	VG
	//####################
	FLT VG_height,VG_length;
	int iMinus1jk,ijk_tmp;
	FLT distance;
	int i_max,j_max,i_min;
	FLT VG_start_x,VG_start_x_rank;
	FLT VG_start_y,VG_start_y_rank;
	//####################
	//	piston
	//####################
	FLT actual_x;
	FLT y_kolben;
	FLT alpha_kolben;

	//####################
	//	small piston
	//####################
	FLT y_kolben_max;
	//####################
	
	switch (pnt_config->IBC_Type)
	{
		case 0:
		//####################
		//	Custom
		//####################
		actual_x=pnt_config->IBC_MovingActualPosition;

		x_min=actual_x-pnt_config->IBC_SizeX/2.;
		x_max=actual_x+pnt_config->IBC_SizeX/2.;
		y_min=pnt_config->IBC_StartpositionY-pnt_config->IBC_SizeY/2.;
		y_max=pnt_config->IBC_StartpositionY+pnt_config->IBC_SizeY/2.;
		z_min=pnt_config->IBC_StartpositionZ-pnt_config->IBC_SizeZ/2.;
		z_max=pnt_config->IBC_StartpositionZ+pnt_config->IBC_SizeZ/2.;

		for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
		{
			for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
			{
				for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

					if
					(
						(pnt_mesh->x[ijk]>=x_min)&&(pnt_mesh->x[ijk]<=x_max)
						&&
						(pnt_mesh->y[ijk]>=y_min)&&(pnt_mesh->y[ijk]<=y_max)
						&&
						(pnt_mesh->z[ijk]>=z_min)&&(pnt_mesh->z[ijk]<=z_max)
					)
					{
						pnt_mesh->flag_IBC[ijk]=1;
					}
					else
					{
						pnt_mesh->flag_IBC[ijk]=0;
					}
				}
			}
		}
		break;		


		case 1:
		//####################
		//	VG
		//####################
		if(pnt_config->MPI_rank==0){printf("SHOCK: IBC for VG is set.\n");}
		VG_height= 0.4/80.; //hoehe: 0.4mm -> entdimensionieren mit c=80mm
		VG_length= 1.0/80.; //laenge: 1.0mm -> entdimensionieren mit c=80mm
		//Dies sind die genauen Koordinaten des ersten Gitterpunktes des VG (Ecke links unten)
		VG_start_x=0.650723;
		VG_start_y=0.0554763;

		for (i=pnt_config->int_iStartGhosts+1; i <= pnt_config->int_iEndGhosts; i++)
		{
			j=pnt_config->int_jStartReal;
			k=pnt_config->int_kStartReal;
			ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
			iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

			if(
					(pnt_mesh->x[iMinus1jk]<0.6505)&&
					(pnt_mesh->x[ijk]>=0.6505)
					)

			{
				//Bestimmung der linken unteren Ecke des VG (j=1,i=i_min)
				i_min=i;
				//Speichern der rank abhaengigen Position der linken unteren Ecke fuer Berechnung der
				//korrektren Laenge
				VG_start_x_rank=pnt_mesh->x[ijk];
				VG_start_y_rank=pnt_mesh->y[ijk];

				//Bestimmung der linken oberen Ecke des VG innerhalb des Rechengebietes des CPU (j=j_max,i=i_min)
				i=i_min;
				j=pnt_config->int_jStartGhosts-1;
				do{
					j++;
					ijk_tmp=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					distance=sqrt(
							pow((pnt_mesh->x[ijk_tmp]-VG_start_x),2)+
							pow((pnt_mesh->y[ijk_tmp]-VG_start_y),2));
				}while((distance<VG_height)&&(j<=pnt_config->int_jEndGhosts));
				j_max=j;
				//Wenn diese Bedingung erfuellt ist, ist auch der unterste Gitterpunkt zu weit entfernt
				//so dass keine IBC gesetzt werden duerfen
				if(j_max==pnt_config->int_jStartGhosts){break;}

				j=pnt_config->int_jStartReal;
				i=i_min;
				do{
					i++;
					ijk_tmp=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					distance=sqrt(
							pow((pnt_mesh->x[ijk_tmp]-VG_start_x_rank),2)+
							pow((pnt_mesh->y[ijk_tmp]-VG_start_y_rank),2));
				}while((distance<VG_length)&&(i<=pnt_config->int_iEndGhosts));
				i_max=i-1;

				ijk=i_min*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_jStartReal*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_kStartReal;
				ijk_tmp=i_max*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j_max*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_kStartReal;
				printf("SHOCK: VG defined between i=%d->%d, j=0->%d\n",
						i_min,i_max,j_max);

				FLT delta_z;
				delta_z=0.1/127.;

				for (i=i_min; i <= i_max; i++)
				{
					for (j=pnt_config->int_jStartGhosts; j < j_max; j++)
					{
						for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
						{
							ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

							if((pnt_mesh->z[ijk]>=14.*delta_z)&&(pnt_mesh->z[ijk]<=27.*delta_z))
							{
								ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
								pnt_mesh->flag_IBC[ijk]=1;
							}

							if((pnt_mesh->z[ijk]>=57.*delta_z)&&(pnt_mesh->z[ijk]<=70.*delta_z))
							{
								ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
								pnt_mesh->flag_IBC[ijk]=1;
							}

							if((pnt_mesh->z[ijk]>=100.*delta_z)&&(pnt_mesh->z[ijk]<=113.*delta_z))
							{
								ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
								pnt_mesh->flag_IBC[ijk]=1;
							}
						}
					}
				}
			}
		}
		break;


		case 2:
		//####################
		//	piston
		//####################
		actual_x=pnt_config->IBC_MovingActualPosition;
		y_kolben=pnt_config->IBC_yKolben;
		alpha_kolben=pnt_config->IBC_alphaKolben;
		k=0;

		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				if(
				(pnt_mesh->x[ijk]>0.0)
				&&
				(pnt_mesh->y[ijk]>y_kolben)
				&&
				(pnt_mesh->x[ijk]>=(actual_x-(pnt_mesh->y[ijk]-y_kolben)/tan(alpha_kolben/180.0*MY_PI)))
					)
				{
					pnt_mesh->flag_IBC[ijk]=1;
				}
				else
				{
					pnt_mesh->flag_IBC[ijk]=0;
				}
			}
		}
		break;


		case 3:
		//####################
		//	small piston
		//####################
		actual_x=pnt_config->IBC_MovingActualPosition;
		y_kolben=pnt_config->IBC_yKolben;
		y_kolben_max=y_kolben+3.0;
		alpha_kolben=pnt_config->IBC_alphaKolben;
		k=0;

		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				if(
				(pnt_mesh->x[ijk]>0.0)
				&&
				(pnt_mesh->y[ijk]>y_kolben)
				&&
				(pnt_mesh->y[ijk]<y_kolben_max)
				&&
				(pnt_mesh->x[ijk]>=(actual_x-(pnt_mesh->y[ijk]-y_kolben)/tan(alpha_kolben/180.0*MY_PI)))
					)
				{
					pnt_mesh->flag_IBC[ijk]=1;
				}
				else
				{
					pnt_mesh->flag_IBC[ijk]=0;
				}
			}
		}
		break;
	}
}

void IBC_set(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				if(pnt_mesh->flag_IBC[ijk]==1)
				{
					if(pnt_config->flag_IBC_Moving==1)
					{
						pnt_U->u[ijk]=//pnt_config->IBC_MovingSpeed/pnt_config->u0_dim;}
							(pnt_config->IBC_MovingActualPosition-pnt_config->IBC_MovingLastPosition)/pnt_config->numericalTau;
					}
					else
					{
						pnt_U->u[ijk]=0.0;
					}
					pnt_U->v[ijk]=0.0;
					pnt_U->w[ijk]=0.0;
					pnt_U->theta1[ijk]=
							pnt_U->u[ijk]*pnt_mesh->xi_x[ijk]+
							pnt_U->v[ijk]*pnt_mesh->xi_y[ijk]+
							pnt_U->w[ijk]*pnt_mesh->xi_z[ijk];
					pnt_U->theta2[ijk]=
							pnt_U->u[ijk]*pnt_mesh->eta_x[ijk]+
							pnt_U->v[ijk]*pnt_mesh->eta_y[ijk]+
							pnt_U->w[ijk]*pnt_mesh->eta_z[ijk];
					pnt_U->theta3[ijk]=
							pnt_U->u[ijk]*pnt_mesh->zeta_x[ijk]+
							pnt_U->v[ijk]*pnt_mesh->zeta_y[ijk]+
							pnt_U->w[ijk]*pnt_mesh->zeta_z[ijk];
					pnt_U->e[ijk]=(0.5*((pnt_U->u[ijk]*pnt_U->u[ijk])+(pnt_U->v[ijk]*pnt_U->v[ijk])+(pnt_U->w[ijk]*pnt_U->w[ijk]))+
												1.0/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
                    
					pnt_U->p[ijk]=1.0;
					pnt_U->rho[ijk]=1.0;

				    pnt_U->T[ijk]=1.0;

					pnt_U->c[ijk]=
							sqrt(pnt_config->Upsilon*
							pnt_config->gammaNumber*pnt_U->p[ijk]/pnt_U->rho[ijk]);

					pnt_U->mue[ijk]=((1.0+pnt_config->SutherlandConstant)*pow(pnt_U->p[ijk]/pnt_U->rho[ijk],1.5)/
							(pnt_U->p[ijk]/pnt_U->rho[ijk]+pnt_config->SutherlandConstant));

					pnt_mesh->BC_Corrector_xiMomentum[ijk]=1.0;
					pnt_mesh->BC_Corrector_etaMomentum[ijk]=1.0;
					pnt_mesh->BC_Corrector_zetaMomentum[ijk]=1.0;
				}
			}
		}
	}
}

void inducePressureWavesPlateau(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	FLT r, Ma_r, c, delta_t;
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				if(pnt_mesh->flag_PressureWaves[ijk]==1)
				{
					r=sqrt(pow((pnt_mesh->x[ijk]-pnt_config->pw_x1),2)+
							pow((pnt_mesh->y[ijk]-pnt_config->pw_y1),2));

					Ma_r=pnt_config->machNumber*(pnt_mesh->x[ijk]-pnt_config->pw_x1)/r
							+sqrt(1.-pnt_config->machNumber*pow((pnt_mesh->y[ijk]-pnt_config->pw_y1),2)/r);

					c=pnt_U_lastStep->c[ijk];

					delta_t=(r-pnt_config->pw_r0*Ma_r)/(Ma_r*c)*pnt_config->L0_dim/pnt_config->u0_dim;

					pnt_U_lastStep->p[ijk]=pnt_mesh->startPressure_PressureWaves[ijk]+
							pnt_config->pw_amplitude*sin(2.0*MY_PI*pnt_config->pw_frequency*(pnt_config->time_dim-delta_t))*sqrt(Ma_r*pnt_config->pw_r0/r);

//					pnt_U_lastStep->p[ijk]=pnt_mesh->startPressure_PressureWaves[ijk]+pnt_config->pw_amplitude*sin(2.0*MY_PI*pnt_config->pw_frequency*(pnt_config->time_dim));

				    pnt_U_lastStep->T[ijk]=pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk];					
					pnt_U_lastStep->rho[ijk]=pow(pnt_U_lastStep->p[ijk],1.0/pnt_config->gammaNumber);


					pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+
							pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);

					pnt_U_lastStep->c[ijk]=sqrt(pnt_config->Upsilon*pnt_config->gammaNumber
							*pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]);

					pnt_U_lastStep->mue[ijk]=((1.0+pnt_config->SutherlandConstant)*pow(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk],1.5)/
							(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]+pnt_config->SutherlandConstant));

				}

			}
		}
	}
}

void inducePressureWaves(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	FLT alpha,r,r_quelle,u_alpha,lambda_alpha,f,lambda,x_Q,y_Q;
	f=pnt_config->pw_frequency/pnt_config->u0_dim*pnt_config->L0_dim;
	lambda=(1.0/pnt_config->machNumber)/f;
	r_quelle=pnt_config->pw_r0;
	x_Q=pnt_config->pw_x0-r_quelle*pnt_config->machNumber;
	y_Q=pnt_config->pw_y0;

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				if(pnt_mesh->flag_PressureWaves[ijk]==1)
				{
					alpha=atan((y_Q-pnt_mesh->y[ijk])/(x_Q-pnt_mesh->x[ijk]));
					if(pnt_mesh->x[ijk]<x_Q){alpha=MY_PI-alpha;}

					u_alpha=1.0*cos(alpha)+sqrt(pow(1.0/pnt_config->machNumber,2.)-pow(1.0*sin(alpha),2.));

					r=sqrt(pow((pnt_mesh->x[ijk]-x_Q),2)+
							pow((pnt_mesh->y[ijk]-y_Q),2));
					lambda_alpha=u_alpha/f;

					pnt_U_lastStep->p[ijk]=pnt_mesh->startPressure_PressureWaves[ijk]+
							pnt_config->pw_amplitude
							*sin(2.0*MY_PI*(pnt_config->pw_frequency*(pnt_config->time_dim-pnt_config->start_Time)-r/lambda_alpha+r_quelle/lambda))
					//*sqrt(lambda_alpha/r*r_quelle/lambda);
					*sqrt(1.0/r);

//					pnt_U_lastStep->p[ijk]=pnt_mesh->startPressure_PressureWaves[ijk]+pnt_config->pw_amplitude*sin(2.0*MY_PI*pnt_config->pw_frequency*(pnt_config->time_dim));

				    pnt_U_lastStep->T[ijk]=pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk];					
					pnt_U_lastStep->rho[ijk]=pow(pnt_U_lastStep->p[ijk],1.0/pnt_config->gammaNumber);

					pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+
							pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);

					pnt_U_lastStep->c[ijk]=sqrt(pnt_config->Upsilon*pnt_config->gammaNumber
							*pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]);

					pnt_U_lastStep->mue[ijk]=((1.0+pnt_config->SutherlandConstant)*pow(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk],1.5)/
							(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]+pnt_config->SutherlandConstant));

				}

			}
		}
	}
}

void preparePressureWaves(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	//Bestimmung des Radius, der genau eine Wellenlaenge enthaelt
	FLT c_tmp,f_tmp;
	f_tmp=pnt_config->pw_frequency/pnt_config->u0_dim*pnt_config->L0_dim;
	c_tmp=(1.0/pnt_config->machNumber);
	pnt_config->pw_r0=c_tmp/f_tmp;

	FLT distance;
	int BC_check;
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				pnt_mesh->flag_PressureWaves[ijk]=0;

				if(pnt_config->pw_UseBC)
				{
					if(
					(i<pnt_config->int_iStartReal)||(j<pnt_config->int_jStartReal)||(k<pnt_config->int_kStartReal)||
					(i>pnt_config->int_iEndReal)||(j>pnt_config->int_jEndReal)||(k>pnt_config->int_kEndReal))
					{
						BC_check=1;
					}
					else
					{
						BC_check=0;
					}

				}
				else
				{
					BC_check=1;
				}

//				distance=sqrt(
//						pow((pnt_mesh->x[ijk]-pnt_config->pw_x0),2)+
//						pow((pnt_mesh->y[ijk]-pnt_config->pw_y0),2)+
//						pow((pnt_mesh->z[ijk]-pnt_config->pw_z0),2));

//				Variante von Volker
				distance=sqrt(
						pow((pnt_mesh->x[ijk]-pnt_config->pw_x0),2)
						+pow((pnt_mesh->y[ijk]-pnt_config->pw_y0),2));

				if((distance<=pnt_config->pw_r0)&&(BC_check))
				{
					pnt_mesh->flag_PressureWaves[ijk]=1;
					if(pnt_config->pw_UseFlowAverage)
					{
						pnt_mesh->startPressure_PressureWaves[ijk]=pnt_U_lastStep->p[ijk];
						pnt_mesh->startDensity_PressureWaves[ijk]=pnt_U_lastStep->rho[ijk];
					}
					else
					{
						pnt_mesh->startPressure_PressureWaves[ijk]=1.0;
						pnt_mesh->startDensity_PressureWaves[ijk]=1.0;
					}
				}

				distance=sqrt(
						pow((pnt_mesh->x[ijk]-pnt_config->pw_x1),2)+
						pow((pnt_mesh->y[ijk]-pnt_config->pw_y1),2)+
						pow((pnt_mesh->z[ijk]-pnt_config->pw_z1),2));

				if((distance<=pnt_config->pw_r0)&&(pnt_config->pw_numberSources==2)&&(BC_check))
				{
					pnt_mesh->flag_PressureWaves[ijk]=1;
					if(pnt_config->pw_UseFlowAverage)
					{
						pnt_mesh->startPressure_PressureWaves[ijk]=pnt_U_lastStep->p[ijk];
						pnt_mesh->startDensity_PressureWaves[ijk]=pnt_U_lastStep->rho[ijk];
					}
					else
					{
						pnt_mesh->startPressure_PressureWaves[ijk]=1.0;
						pnt_mesh->startDensity_PressureWaves[ijk]=1.0;
					}
				}


			}
		}
	}
}

void preparePressureHistoryValues(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{
	FLT *distance,*distance_tmp;
	distance = (FLT *)calloc(pnt_config->PressureHistory_No, sizeof(FLT ));
	distance_tmp = (FLT *)calloc(pnt_config->PressureHistory_No, sizeof(FLT ));

	int p;

	for (p=0;p<pnt_config->PressureHistory_No;p++)
	{
		distance[p]=999999999.9;
	}

	for (p=0;p<pnt_config->PressureHistory_No;p++)
	{
		for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
		{
			for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
			{
				for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

					if (MESHDIMENSIONS==2)
					{
						distance_tmp[p]=sqrt(
									pow((pnt_mesh->x[ijk]-pnt_config->PressureHistory_x_P[p]),2.)+
									pow((pnt_mesh->y[ijk]-pnt_config->PressureHistory_y_P[p]),2.));
					}
					if (MESHDIMENSIONS==3)
					{
						distance_tmp[p]=sqrt(
									pow((pnt_mesh->x[ijk]-pnt_config->PressureHistory_x_P[p]),2.)+
									pow((pnt_mesh->y[ijk]-pnt_config->PressureHistory_y_P[p]),2.)+
									pow((pnt_mesh->z[ijk]-pnt_config->PressureHistory_z_P[p]),2.));
					}

					if(distance_tmp[p]<distance[p])
					{
						if((i>=pnt_config->int_iStartReal && i <= pnt_config->int_iEndReal)&&(j>=pnt_config->int_jStartReal && j <= pnt_config->int_jEndReal)&&(k>=pnt_config->int_kStartReal && k <= pnt_config->int_kEndReal))
						{
							pnt_config->ijk_PressureHistory_P[p]=ijk;
							pnt_config->flag_PressureHistory_P[p]=1;
							pnt_config->PressureHistory_x_P_real[p]=pnt_mesh->x[ijk];
							pnt_config->PressureHistory_y_P_real[p]=pnt_mesh->y[ijk];
							pnt_config->PressureHistory_z_P_real[p]=pnt_mesh->z[ijk];
						}
						else
						{
							pnt_config->flag_PressureHistory_P[p]=0;
						}
						distance[p]=distance_tmp[p];
					}
				}
			}
		}
	}

	for (p=0;p<pnt_config->PressureHistory_No;p++)
	{
		if(pnt_config->flag_PressureHistory_P[p]==1)
		{
			printf("Druckaufzeichnung bei x=%g, y=%g, z=%g (Vorgabe: x=%g, y=%g, z=%g) von Rank %d\n",
					(double)pnt_config->PressureHistory_x_P_real[p],
					(double)pnt_config->PressureHistory_y_P_real[p],
					(double)pnt_config->PressureHistory_z_P_real[p],
					(double)pnt_config->PressureHistory_x_P[p],
					(double)pnt_config->PressureHistory_y_P[p],
					(double)pnt_config->PressureHistory_z_P[p],
					pnt_config->MPI_rank);
		}
	}
}

void prepareVelocityHistoryValues(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{
	FLT *distance,*distance_tmp;
	distance = (FLT *)calloc(pnt_config->VelocityHistory_No, sizeof(FLT ));
	distance_tmp = (FLT *)calloc(pnt_config->VelocityHistory_No, sizeof(FLT ));

	int p;

	for (p=0;p<pnt_config->VelocityHistory_No;p++)
	{
		distance[p]=999999999.9;
	}

	for (p=0;p<pnt_config->VelocityHistory_No;p++)
	{
		for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
		{
			for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
			{
				for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

					if (MESHDIMENSIONS==2)
					{
						distance_tmp[p]=sqrt(
									pow((pnt_mesh->x[ijk]-pnt_config->VelocityHistory_x_P[p]),2.)+
									pow((pnt_mesh->y[ijk]-pnt_config->VelocityHistory_y_P[p]),2.));
					}
					if (MESHDIMENSIONS==3)
					{
						distance_tmp[p]=sqrt(
									pow((pnt_mesh->x[ijk]-pnt_config->VelocityHistory_x_P[p]),2.)+
									pow((pnt_mesh->y[ijk]-pnt_config->VelocityHistory_y_P[p]),2.)+
									pow((pnt_mesh->z[ijk]-pnt_config->VelocityHistory_z_P[p]),2.));
					}

					if(distance_tmp[p]<distance[p])
					{
						if((i>=pnt_config->int_iStartReal && i <= pnt_config->int_iEndReal)&&(j>=pnt_config->int_jStartReal && j <= pnt_config->int_jEndReal)&&(k>=pnt_config->int_kStartReal && k <= pnt_config->int_kEndReal))
						{
							pnt_config->ijk_VelocityHistory_P[p]=ijk;
							pnt_config->flag_VelocityHistory_P[p]=1;
							pnt_config->VelocityHistory_x_P_real[p]=pnt_mesh->x[ijk];
							pnt_config->VelocityHistory_y_P_real[p]=pnt_mesh->y[ijk];
							pnt_config->VelocityHistory_z_P_real[p]=pnt_mesh->z[ijk];
						}
						else
						{
							pnt_config->flag_VelocityHistory_P[p]=0;
						}
						distance[p]=distance_tmp[p];
					}
				}
			}
		}
	}

	for (p=0;p<pnt_config->VelocityHistory_No;p++)
	{
		if(pnt_config->flag_VelocityHistory_P[p]==1)
		{
			printf("Geschwindigkeitsaufzeichnung bei x=%g, y=%g, z=%g (Vorgabe: x=%g, y=%g, z=%g) von Rank %d\n",
					(double)pnt_config->VelocityHistory_x_P_real[p],
					(double)pnt_config->VelocityHistory_y_P_real[p],
					(double)pnt_config->VelocityHistory_z_P_real[p],
					(double)pnt_config->VelocityHistory_x_P[p],
					(double)pnt_config->VelocityHistory_y_P[p],
					(double)pnt_config->VelocityHistory_z_P[p],
					pnt_config->MPI_rank);
		}
	}
}

void InitializeVortex(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	FLT x_quer,y_quer;
	FLT x_wirb_zentr,y_wirb_zentr;
	FLT radius;
	FLT beta;
	FLT faktor_quer;
	FLT r_wirb_max;

	x_wirb_zentr=pnt_config->Vortex_x_wirb_zentr;
	y_wirb_zentr=pnt_config->Vortex_y_wirb_zentr;

	faktor_quer=pnt_config->Vortex_faktor_quer;
	r_wirb_max=pnt_config->Vortex_r_wirb_max;


	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				x_quer = pnt_mesh->x[ijk]-x_wirb_zentr;
				y_quer = pnt_mesh->y[ijk]-y_wirb_zentr;

				x_quer *=  faktor_quer;   y_quer *=  faktor_quer;

				radius = x_quer*x_quer + y_quer*y_quer;

				if(radius <= r_wirb_max*r_wirb_max)
				{
					beta = pnt_config->Vortex_beta;

					pnt_U_lastStep->u[ijk] += y_quer*beta/(2.*MY_PI)*exp((1.- radius)/2.);

					pnt_U_lastStep->v[ijk] -= x_quer*beta/(2.*MY_PI)*exp((1.- radius)/2.);

					beta = 1.;

					pnt_U_lastStep->rho[ijk] = pow((pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]
					- (pnt_config->gammaNumber - 1.)*beta*beta/
					(8.*pnt_config->gammaNumber*MY_PI*MY_PI)/pnt_config->Upsilon*exp(1. - radius)), (1./(pnt_config->gammaNumber - 1.)));

					pnt_U_lastStep->p[ijk] = pow(pnt_U_lastStep->rho[ijk], pnt_config->gammaNumber);


					pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+
							pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
				}

			}
		}
	}
}

void DefineFilesPath(
        struct strct_configuration * pnt_config)
{
	int i=strlen(pnt_config->chr_MeshPath);
	int j,k,k2;
	while(pnt_config->chr_MeshPath[i] != '/' && i>-1)
	{
		i--;
	}


	i++;
	k=0;
	k2=0;
	for(j=i;j<=strlen(pnt_config->chr_MeshPath);j++)
	{
		if(pnt_config->chr_MeshPath[j] == '.')
		{
		  pnt_config->chr_DivisionFile[k2]='.';
		  pnt_config->chr_DivisionFile[k2+1]='s';
		  pnt_config->chr_DivisionFile[k2+2]='d';
		  pnt_config->chr_DivisionFile[k2+3]='d';
		  pnt_config->chr_DivisionFile[k2+4]='c';
		  pnt_config->chr_DivisionFile[k2+5]='\0';
		  k2+=6;
		}

		pnt_config->chr_MeshFile[k]=pnt_config->chr_MeshPath[j];
		pnt_config->chr_DivisionFile[k2]=pnt_config->chr_MeshPath[j];


		k++;
		k2++;
	}
	k=0;
	for(j=0;j<i;j++)
	{
		pnt_config->chr_folder[k]=pnt_config->chr_MeshPath[j];
		k++;
	}

	pnt_config->chr_folder[k]=pnt_config->chr_MeshPath[strlen(pnt_config->chr_MeshPath)];
	if (i==0){strcpy(pnt_config->chr_folder,"./");}

	sprintf(pnt_config->chr_PressureHistoryFile,"PressureHistory_%s",pnt_config->chr_MeshFile);
	sprintf(pnt_config->chr_VelocityHistoryFile,"VelocityHistory_%s",pnt_config->chr_MeshFile);

	sprintf(pnt_config->chr_DivisionPath,"%s%s",pnt_config->chr_folder,pnt_config->chr_DivisionFile);
	sprintf(pnt_config->chr_PressureHistoryPath,"%s%s",pnt_config->chr_folder,pnt_config->chr_PressureHistoryFile);
	sprintf(pnt_config->chr_VelocityHistoryPath,"%s%s",pnt_config->chr_folder,pnt_config->chr_VelocityHistoryFile);


	strcpy(pnt_config->chr_MeshPathOriginal,pnt_config->chr_MeshPath);
}

void WriteValuesForPressureHistory(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U_lastStep)
{
	int p;

	pnt_config->PressureHistory_time[pnt_config->int_actualIteration-(pnt_config->int_StartIteration+1)]=pnt_config->time_dim;
	for (p=0;p<pnt_config->PressureHistory_No;p++)
	{
		if(pnt_config->flag_PressureHistory_P[p]==1)
		{
			pnt_config->PressureHistory_pressure[p][pnt_config->int_actualIteration-(pnt_config->int_StartIteration+1)]=pnt_U_lastStep->p[pnt_config->ijk_PressureHistory_P[p]];
		}
	}


}

void WriteValuesForVelocityHistory(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U_lastStep)
{
	int p;

	pnt_config->VelocityHistory_time[pnt_config->int_actualIteration-(pnt_config->int_StartIteration+1)]=pnt_config->time_dim;
	for (p=0;p<pnt_config->VelocityHistory_No;p++)
	{
		if(pnt_config->flag_VelocityHistory_P[p]==1)
		{
			pnt_config->VelocityHistory_VelocityX[p][pnt_config->int_actualIteration-(pnt_config->int_StartIteration+1)]=pnt_U_lastStep->u[pnt_config->ijk_VelocityHistory_P[p]];
			pnt_config->VelocityHistory_VelocityY[p][pnt_config->int_actualIteration-(pnt_config->int_StartIteration+1)]=pnt_U_lastStep->v[pnt_config->ijk_VelocityHistory_P[p]];
			pnt_config->VelocityHistory_VelocityZ[p][pnt_config->int_actualIteration-(pnt_config->int_StartIteration+1)]=pnt_U_lastStep->w[pnt_config->ijk_VelocityHistory_P[p]];
		}
	}


}

int checkNAN(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U,
		struct strct_mesh * pnt_mesh)
{
	int *NANCounterArray = (int *)calloc(pnt_config->MPI_size,sizeof(int));
	int NANCounter=0;
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				if( (isinf(pnt_U->rho[ijk])||isnan(pnt_U->rho[ijk])) && (pnt_mesh->flag_IBC[ijk]==0))
				{
					NANCounter++;
					/*
					printf("SHOCK: NAN bei x=%g, y=%g, z=%g\n",
							(double)pnt_mesh->x[ijk],
							(double)pnt_mesh->y[ijk],
							(double)pnt_mesh->z[ijk]);
							*/
//					MPI_Abort(pnt_config->MPI_comm,13370);
					pnt_config->flag_NAN=1;

					pnt_U->u[ijk]=pnt_U->v[ijk]=pnt_U->w[ijk]=pnt_U->p[ijk]=pnt_U->rho[ijk]=1000000;
				}
			}
		}
	}
	MPI_Allgather(&NANCounter, 1, MPI_INT, NANCounterArray, 1, MPI_INT,pnt_config->MPI_comm);
	int i,sum;
	sum=0;
	for(i=0;i<pnt_config->MPI_size;i++)
	{
		sum+=NANCounterArray[i];
	}
	free(NANCounterArray);
	if (sum>0)
	{
//		if(pnt_config->MPI_rank==0){printf("\nSHOCK: ---->  NAN bei %d\n",pnt_config->int_actualIteration);}
		return 1;
	}
	else
	{
		return 0;
	}
}

void IBC_BornCells(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK)
{
	int iMinus1jk,iPlus1jk,ijMinus1k,ijPlus1k;

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				if((pnt_mesh->flag_IBC_last[ijk]==1)&&(pnt_mesh->flag_IBC[ijk]==0))
				{

					iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					iPlus1jk=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;
					ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;

					pnt_U_RK->u[ijk]=
							(((1-pnt_mesh->flag_IBC[iMinus1jk])*pnt_U_RK->u[iMinus1jk]+
							(1-pnt_mesh->flag_IBC[iPlus1jk])*pnt_U_RK->u[iPlus1jk]+
							(1-pnt_mesh->flag_IBC[ijMinus1k])*pnt_U_RK->u[ijMinus1k]+
							(1-pnt_mesh->flag_IBC[ijPlus1k])*pnt_U_RK->u[ijPlus1k])/
							(4-pnt_mesh->flag_IBC[iMinus1jk]-pnt_mesh->flag_IBC[iPlus1jk]-
							pnt_mesh->flag_IBC[ijMinus1k]-pnt_mesh->flag_IBC[ijPlus1k]));

					pnt_U_RK->v[ijk]=
							(((1-pnt_mesh->flag_IBC[iMinus1jk])*pnt_U_RK->v[iMinus1jk]+
							(1-pnt_mesh->flag_IBC[iPlus1jk])*pnt_U_RK->v[iPlus1jk]+
							(1-pnt_mesh->flag_IBC[ijMinus1k])*pnt_U_RK->v[ijMinus1k]+
							(1-pnt_mesh->flag_IBC[ijPlus1k])*pnt_U_RK->v[ijPlus1k])/
							(4-pnt_mesh->flag_IBC[iMinus1jk]-pnt_mesh->flag_IBC[iPlus1jk]-
							pnt_mesh->flag_IBC[ijMinus1k]-pnt_mesh->flag_IBC[ijPlus1k]));

					pnt_U_RK->p[ijk]=
							(((1-pnt_mesh->flag_IBC[iMinus1jk])*pnt_U_RK->p[iMinus1jk]+
							(1-pnt_mesh->flag_IBC[iPlus1jk])*pnt_U_RK->p[iPlus1jk]+
							(1-pnt_mesh->flag_IBC[ijMinus1k])*pnt_U_RK->p[ijMinus1k]+
							(1-pnt_mesh->flag_IBC[ijPlus1k])*pnt_U_RK->p[ijPlus1k])/
							(4-pnt_mesh->flag_IBC[iMinus1jk]-pnt_mesh->flag_IBC[iPlus1jk]-
							pnt_mesh->flag_IBC[ijMinus1k]-pnt_mesh->flag_IBC[ijPlus1k]));

					pnt_U_RK->rho[ijk]=
							(((1-pnt_mesh->flag_IBC[iMinus1jk])*pnt_U_RK->rho[iMinus1jk]+
							(1-pnt_mesh->flag_IBC[iPlus1jk])*pnt_U_RK->rho[iPlus1jk]+
							(1-pnt_mesh->flag_IBC[ijMinus1k])*pnt_U_RK->rho[ijMinus1k]+
							(1-pnt_mesh->flag_IBC[ijPlus1k])*pnt_U_RK->rho[ijPlus1k])/
							(4-pnt_mesh->flag_IBC[iMinus1jk]-pnt_mesh->flag_IBC[iPlus1jk]-
							pnt_mesh->flag_IBC[ijMinus1k]-pnt_mesh->flag_IBC[ijPlus1k]));

					pnt_U_RK->theta1[ijk]=
							pnt_U_RK->u[ijk]*pnt_mesh->xi_x[ijk]+
							pnt_U_RK->v[ijk]*pnt_mesh->xi_y[ijk];

					pnt_U_RK->theta2[ijk]=
							pnt_U_RK->u[ijk]*pnt_mesh->eta_x[ijk]+
							pnt_U_RK->v[ijk]*pnt_mesh->eta_y[ijk];

//					pnt_U_RK->p[ijk]=pnt_config->InitializeValues_p1/20.;
//					pnt_U_RK->rho[ijk]=pnt_config->InitializeValues_rho1/20.;
//					pnt_U_RK->u[ijk]=0.0;
//					pnt_U_RK->v[ijk]=0.0;
//					pnt_U_RK->w[ijk]=0.0;
					pnt_U_RK->e[ijk]=(0.5*((pnt_U_RK->u[ijk]*pnt_U_RK->u[ijk])+(pnt_U_RK->v[ijk]*pnt_U_RK->v[ijk])+(pnt_U_RK->w[ijk]*pnt_U_RK->w[ijk]))+
							pnt_U_RK->p[ijk]/pnt_U_RK->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);

//					pnt_U_RK->theta1[ijk]=0.0;
//					pnt_U_RK->theta2[ijk]=0.0;
//					pnt_U_RK->theta3[ijk]=0.0;

					pnt_U_RK->T[ijk]=
						fabs(pnt_U_RK->p[ijk]/pnt_U_RK->rho[ijk]);

					pnt_U_RK->c[ijk]=
						sqrt(pnt_config->Upsilon*
						pnt_config->gammaNumber*pnt_U_RK->p[ijk]/pnt_U_RK->rho[ijk]);

					pnt_U_RK->mue[ijk]=((1.0+pnt_config->SutherlandConstant)*pow(pnt_U_RK->p[ijk]/pnt_U_RK->rho[ijk],1.5)/
						(pnt_U_RK->p[ijk]/pnt_U_RK->rho[ijk]+pnt_config->SutherlandConstant));

					//###################################

					pnt_U_lastStep->p[ijk]=pnt_U_RK->p[ijk];
					pnt_U_lastStep->rho[ijk]=pnt_U_RK->rho[ijk];
					pnt_U_lastStep->u[ijk]=pnt_U_RK->u[ijk];
					pnt_U_lastStep->v[ijk]=pnt_U_RK->v[ijk];
					pnt_U_lastStep->w[ijk]=pnt_U_RK->w[ijk];
					pnt_U_lastStep->e[ijk]=pnt_U_RK->e[ijk];
					pnt_U_lastStep->theta1[ijk]=pnt_U_RK->theta1[ijk];
					pnt_U_lastStep->theta2[ijk]=pnt_U_RK->theta2[ijk];
					pnt_U_lastStep->theta3[ijk]=pnt_U_RK->theta3[ijk];
					pnt_U_lastStep->T[ijk]=pnt_U_RK->T[ijk];
					pnt_U_lastStep->c[ijk]=pnt_U_RK->c[ijk];
					pnt_U_lastStep->mue[ijk]=pnt_U_RK->mue[ijk];

					pnt_mesh->BC_Corrector_xiMomentum[ijk]=1.0;
					pnt_mesh->BC_Corrector_etaMomentum[ijk]=1.0;
					pnt_mesh->BC_Corrector_zetaMomentum[ijk]=1.0;
				}
			}

		}
	}
}

void IBC_ApplyBC4FluxInXi(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int ijk2,ijk3,int_direction_corrector,int_symmetryIndex,ijkSymmetry;

	FLT corrector[3]={-1.,1.,1.};

	for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
	{
	for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
	{
		for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
		{
			ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
			if(i==0){
					ijk2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					ijk3=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
            }else if (i==pnt_config->int_iEndGhosts){
					ijk2=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					ijk3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
            }else{
					ijk2=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					ijk3=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
			}
			if(
					((pnt_mesh->flag_IBC[ijk]==1)&&(pnt_mesh->flag_IBC[ijk2]==0))
					||
					((pnt_mesh->flag_IBC[ijk]==1)&&(pnt_mesh->flag_IBC[ijk3]==0))
			)
			{
				if((pnt_mesh->flag_IBC[ijk]==1)&&(pnt_mesh->flag_IBC[ijk3]==0))
				{int_direction_corrector=-1;}
				else{int_direction_corrector=1;}

				for (int_symmetryIndex=1; int_symmetryIndex <=(pnt_config->int_SpaceOrder+1)/2 ; int_symmetryIndex++)
				{
					ijk=(i+(int_symmetryIndex-1)*int_direction_corrector)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					ijkSymmetry=(i-int_symmetryIndex*int_direction_corrector)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

					if(
							((i+(int_symmetryIndex-1)*int_direction_corrector)<=pnt_config->int_iEndGhosts)&&
							((i+(int_symmetryIndex-1)*int_direction_corrector)>=pnt_config->int_iStartGhosts)&&
							((i-int_symmetryIndex*int_direction_corrector)<=pnt_config->int_iEndGhosts)&&
							((i-int_symmetryIndex*int_direction_corrector)>=pnt_config->int_iStartGhosts)
							)
					{
						if((pnt_config->flag_IBC_Moving==1) && (pnt_mesh->y[ijk] > 0.5))
						{
							WriteMovingWallNoSlipIsothermalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}
						else if((pnt_config->flag_IBC_Moving==1) && (pnt_mesh->y[ijk] <= 0.5))
						{
							WriteMovingWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}
						/*else if(pnt_config->flag_Inviscid==1)
						{
							WriteWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}*/
						else
						{
							//Bei festen IBC kann die Metrik kopiert werden.
							//Bei bewegten IBC erscheint das Problem der "born" Cells, die ihre urspruengliche Metric brauchen.
							copyMetric(pnt_config,pnt_mesh,ijk,ijkSymmetry);
							WriteWallNoSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}
					}
				}
			}
		}
	}
	}
}

void IBC_ApplyBC4FluxInEta(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int ijk2,ijk3,int_direction_corrector,int_symmetryIndex,ijkSymmetry;

	FLT corrector[3]={1.,-1.,1.};

	for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
	{
	for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
	{
		for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
		{
			ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
			if(j==0){
					ijk2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					ijk3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;
            }else if (j==pnt_config->int_jEndGhosts){
					ijk2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;
					ijk3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
            }else{
					ijk2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;
					ijk3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;
			}

			if(
					((pnt_mesh->flag_IBC[ijk]==1)&&(pnt_mesh->flag_IBC[ijk2]==0))
					||
					((pnt_mesh->flag_IBC[ijk]==1)&&(pnt_mesh->flag_IBC[ijk3]==0))
			)
			{
				if((pnt_mesh->flag_IBC[ijk]==1)&&(pnt_mesh->flag_IBC[ijk3]==0))
				{int_direction_corrector=-1;}
				else{int_direction_corrector=1;}

				for (int_symmetryIndex=1; int_symmetryIndex <=(pnt_config->int_SpaceOrder+1)/2 ; int_symmetryIndex++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+(int_symmetryIndex-1)*int_direction_corrector)*pnt_config->int_kMeshPointsGhostCells+k;
					ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-int_symmetryIndex*int_direction_corrector)*pnt_config->int_kMeshPointsGhostCells+k;

					if(
							((j+(int_symmetryIndex-1)*int_direction_corrector)<=pnt_config->int_jEndGhosts)&&
							((j+(int_symmetryIndex-1)*int_direction_corrector)>=pnt_config->int_jStartGhosts)&&
							((j-int_symmetryIndex*int_direction_corrector)<=pnt_config->int_jEndGhosts)&&
							((j-int_symmetryIndex*int_direction_corrector)>=pnt_config->int_jStartGhosts)
							)
					{
						if((pnt_config->flag_IBC_Moving==1) && (pnt_mesh->y[ijk] > 0.5))
						{
							WriteMovingWallNoSlipIsothermalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}
						else if((pnt_config->flag_IBC_Moving==1) && (pnt_mesh->y[ijk] <= 0.5))
						{
							WriteMovingWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}
						/*else if(pnt_config->flag_Inviscid==1)
						{
							WriteWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}*/
						else
						{
							WriteWallNoSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}
					}
				}
			}
		}
	}
	}
}

void IBC_ApplyBC4FluxInZeta(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int ijk2,ijk3,int_direction_corrector,int_symmetryIndex,ijkSymmetry;

	FLT corrector[3]={1.,1.,-1.};

	for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
	{
	for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
	{
		for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
		{
			ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
			if(k==0){
					ijk2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					ijk3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+1;
            }else if(k==pnt_config->int_kEndGhosts){
					ijk2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-1;
					ijk3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
            }else{
					ijk2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-1;
					ijk3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+1;
			}
			if(
					((pnt_mesh->flag_IBC[ijk]==1)&&(pnt_mesh->flag_IBC[ijk2]==0))
					||
					((pnt_mesh->flag_IBC[ijk]==1)&&(pnt_mesh->flag_IBC[ijk3]==0))
			)
			{
				if((pnt_mesh->flag_IBC[ijk]==1)&&(pnt_mesh->flag_IBC[ijk3]==0))
				{int_direction_corrector=-1;}
				else{int_direction_corrector=1;}


				for (int_symmetryIndex=1; int_symmetryIndex <= (pnt_config->int_SpaceOrder+1)/2; int_symmetryIndex++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+(k+(int_symmetryIndex-1)*int_direction_corrector);
					ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+(k-int_symmetryIndex*int_direction_corrector);

					if(
							((k+(int_symmetryIndex-1)*int_direction_corrector)<=pnt_config->int_kEndGhosts)&&
							((k+(int_symmetryIndex-1)*int_direction_corrector)>=pnt_config->int_kStartGhosts)&&
							((k-int_symmetryIndex*int_direction_corrector)<=pnt_config->int_kEndGhosts)&&
							((k-int_symmetryIndex*int_direction_corrector)>=pnt_config->int_kStartGhosts)
							)
					{
						if((pnt_config->flag_IBC_Moving==1) && (pnt_mesh->y[ijk] > 0.5))
						{
							WriteMovingWallNoSlipIsothermalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}
						else if((pnt_config->flag_IBC_Moving==1) && (pnt_mesh->y[ijk] <= 0.5))
						{
							WriteMovingWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}
						/*else if(pnt_config->flag_Inviscid==1)
						{
							WriteWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}*/
						else
						{
							WriteWallNoSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
						}
					}
				}
			}
		}
	}
	}
}

void changeMeshInto2D(
		struct strct_configuration * pnt_config)
{
	pnt_config->int_kStartReal=pnt_config->int_kMid;
	pnt_config->int_kEndReal=pnt_config->int_kMid;

	pnt_config->int_kStartGhosts=pnt_config->int_kMid;
	pnt_config->int_kEndGhosts=pnt_config->int_kMid;
}

void changeMeshInto3D(
		struct strct_configuration * pnt_config)
{
	pnt_config->int_kStartReal=pnt_config->int_kStartReal_original;
	pnt_config->int_kEndReal=pnt_config->int_kEndReal_original;

	pnt_config->int_kStartGhosts=pnt_config->int_kStartGhosts_original;
	pnt_config->int_kEndGhosts=pnt_config->int_kEndGhosts_original;
}

void WriteConstantZValues(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U_lastStep)
{
	int ijk,ijk_0;
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
 				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
 				ijk_0=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+pnt_config->int_kMid;

				pnt_U_lastStep->e[ijk]=pnt_U_lastStep->e[ijk_0];
				pnt_U_lastStep->rho[ijk]=pnt_U_lastStep->rho[ijk_0];
				pnt_U_lastStep->p[ijk]=pnt_U_lastStep->p[ijk_0];
				pnt_U_lastStep->u[ijk]=pnt_U_lastStep->u[ijk_0];
				pnt_U_lastStep->v[ijk]=pnt_U_lastStep->v[ijk_0];
				pnt_U_lastStep->w[ijk]=0.0;
			}
		}
	}
}

FLT CalcLambda2(
		int i,
		int j,
		int k,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	int index0_xi,index1_xi;
	int index0_eta,index1_eta;
	int index0_zeta,index1_zeta;
	FLT u_iMinusHalf, u_iPlusHalf, u_jMinusHalf, u_jPlusHalf, u_kMinusHalf, u_kPlusHalf;
	FLT v_iMinusHalf, v_iPlusHalf, v_jMinusHalf, v_jPlusHalf, v_kMinusHalf, v_kPlusHalf;
	FLT w_iMinusHalf, w_iPlusHalf, w_jMinusHalf, w_jPlusHalf, w_kMinusHalf, w_kPlusHalf;

	int ijk,m;
	FLT p,q,rrr,al,yy1,yy2,yy3,lam1,lam2,lam3,lam,lamax,lamin;
	FLT xi_x,eta_x,zeta_x;
	FLT xi_y,eta_y,zeta_y;
	FLT xi_z,eta_z,zeta_z;

	FLT udx,udy,udz;
	FLT vdx,vdy,vdz;
	FLT wdx,wdy,wdz;

	FLT u_xi,u_eta,u_zeta;
	FLT v_xi,v_eta,v_zeta;
	FLT w_xi,w_eta,w_zeta;

	FLT om12,om13,om23,oo11,oo12,oo13,oo22,oo23,oo33;
	FLT s11,s12,s13,s22,s23,s33,ss11,ss12,ss13,ss22,ss23,ss33;
	FLT a11,a12,a13,a22,a23,a33,f1,f2,f3;

	lam =0.;

	/****************************************************************************/
	ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

	xi_x=pnt_mesh->xi_x[ijk];
	xi_y=pnt_mesh->xi_y[ijk];
	xi_z=pnt_mesh->xi_z[ijk];
	eta_x=pnt_mesh->eta_x[ijk];
	eta_y=pnt_mesh->eta_y[ijk];
	eta_z=pnt_mesh->eta_z[ijk];
	zeta_x=pnt_mesh->zeta_x[ijk];
	zeta_y=pnt_mesh->zeta_y[ijk];
	zeta_z=pnt_mesh->zeta_z[ijk];

	u_iMinusHalf=0.;		v_iMinusHalf=0.;			w_iMinusHalf=0.;
	u_jMinusHalf=0.;		v_jMinusHalf=0.;			w_jMinusHalf=0.;
	u_kMinusHalf=0.;		v_kMinusHalf=0.;			w_kMinusHalf=0.;
	u_iPlusHalf=0.;			v_iPlusHalf=0.;				w_iPlusHalf=0.;
	u_jPlusHalf=0.;			v_jPlusHalf=0.;				w_jPlusHalf=0.;
	u_kPlusHalf=0.;			v_kPlusHalf=0.;				w_kPlusHalf=0.;


	for(m=0;m<=SPACEORDER;m++)
	{
		index0_xi=(i+(m-(SPACEORDER+1)/2))*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
		index1_xi=(i+(m-(SPACEORDER+1)/2)+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
		index0_eta=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+(m-(SPACEORDER+1)/2))*pnt_config->int_kMeshPointsGhostCells+k;
		index1_eta=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+(m-(SPACEORDER+1)/2+1))*pnt_config->int_kMeshPointsGhostCells+k;
		index0_zeta=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+(k+(m-(SPACEORDER+1)/2));
		index1_zeta=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+(k+(m-(SPACEORDER+1)/2+1));

		u_iMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->u[index0_xi];		u_iPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->u[index1_xi];
		v_iMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->v[index0_xi];		v_iPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->v[index1_xi];
		w_iMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->w[index0_xi];		w_iPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->w[index1_xi];
		u_jMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->u[index0_eta];		u_jPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->u[index1_eta];
		v_jMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->v[index0_eta];		v_jPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->v[index1_eta];
		w_jMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->w[index0_eta];		w_jPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->w[index1_eta];
		u_kMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->u[index0_zeta];	u_kPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->u[index1_zeta];
		v_kMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->v[index0_zeta];	v_kPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->v[index1_zeta];
		w_kMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->w[index0_zeta];	w_kPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[m]*pnt_U_lastStep->w[index1_zeta];
	}

	u_xi=(u_iPlusHalf-u_iMinusHalf);
	v_xi=(v_iPlusHalf-v_iMinusHalf);
	w_xi=(w_iPlusHalf-w_iMinusHalf);

	u_eta=(u_jPlusHalf-u_jMinusHalf);
	v_eta=(v_jPlusHalf-v_jMinusHalf);
	w_eta=(w_jPlusHalf-w_jMinusHalf);

	u_zeta=(u_kPlusHalf-u_kMinusHalf);
	v_zeta=(v_kPlusHalf-v_kMinusHalf);
	w_zeta=(w_kPlusHalf-w_kMinusHalf);

	udx=u_xi*xi_x+u_eta*eta_x+u_zeta*zeta_x;
	udy=u_xi*xi_y+u_eta*eta_y+u_zeta*zeta_y;
	udz=u_xi*xi_z+u_eta*eta_z+u_zeta*zeta_z;
	vdx=v_xi*xi_x+v_eta*eta_x+v_zeta*zeta_x;
	vdy=v_xi*xi_y+v_eta*eta_y+v_zeta*zeta_y;
	vdz=v_xi*xi_z+v_eta*eta_z+v_zeta*zeta_z;
	wdx=w_xi*xi_x+w_eta*eta_x+w_zeta*zeta_x;
	wdy=w_xi*xi_y+w_eta*eta_y+w_zeta*zeta_y;
	wdz=w_xi*xi_z+w_eta*eta_z+w_zeta*zeta_z;

	/****************************************************************************/

	om12 = 0.5 * (udy-vdx);
	om13 = 0.5 * (udz-wdx);
	om23 = 0.5 * (vdz-wdy);

	/****************************************************************************/

	s11 = 0.5 * (udx);
	s12 = 0.5 * (udy+vdx);
	s13 = 0.5 * (udz+wdx);
	s22 = 0.5 * (vdy);
	s23 = 0.5 * (vdz+wdy);
	s33 = 0.5 * (wdz);

	/*****************************************************************************/

	oo11 =-om12*om12-om13*om13;
	oo12 =-om13*om23;
	oo13 =om12*om23;
	oo22 =-om12*om12-om23*om23;
	oo23 =-om12*om13;
	oo33 =-om13*om13-om23*om23;

	/*****************************************************************************/

	ss11 =s11*s11+s12*s12+s13*s13;
	ss12 =s11*s12+s12*s22+s13*s23;
	ss13 =s11*s13+s12*s23+s13*s33;
	ss22 =s12*s12+s22*s22+s23*s23;
	ss23 =s12*s13+s22*s23+s23*s33;
	ss33 =s13*s13+s23*s23+s33*s33;

	/******************************************************************************/

	a11 = ss11+oo11;
	a12 = ss12+oo12;
	a13 = ss13+oo13;
	a22 = ss22+oo22;
	a23 = ss23+oo23;
	a33 = ss33+oo33;

	/******************************************************************************/

	f1 = -a11-a22-a33;
	f2 = -a23*a23-a12*a12-a13*a13+a11*a33+a11*a22+a22*a33;
	f3 = -a11*a22*a33+a11*a23*a23+a22*a13*a13+a33*a12*a12-2.*a12*a13*a23;

	/******************************************************************************/

	p = -f1*f1/3.+f2;    q = 2.*f1*f1*f1/27.-f1*f2/3.+f3;

	if(p < 0.){rrr=-p;}
	else      {rrr=p; }
	rrr = sqrt(rrr/3.);

	al = acos(q/(2.*rrr*rrr*rrr));

	yy1 =  -2.*rrr*cos(al/3.);
	yy2 = 2.*rrr*cos(MY_PI/3.-al/3.);
	yy3 = 2.*rrr*cos(MY_PI/3.+al/3.);

	/**********************************************************************/
	lam1 = yy1 - f1/3;    lam2 = yy2 - f1/3;    lam3 = yy3 - f1/3;
	/*****************************************************/
	/*****************************************************/
	/****************************************************/

	if(lam2>=lam1)
	{
		lamax=lam2;
	}
	else
	{
		lamax=lam1;
	}
	if(lam3>=lamax)
	{
		lamax=lam3;
	}

	if(lam2<=lam1)
	{
		lamin=lam2;
	}
	else
	{
		lamin=lam1;
	}
	if(lam3<=lamin)
	{
		lamin=lam3;
	}

	if(lam1==lam2 || lam2==lam3)
	{
		lam = lam2;
	}

	if(lam1>lamin && lam1<lamax)
	{
		lam = lam1;
	}
	else
	{
		if(lam2>lamin && lam2<lamax)
		{
			lam = lam2;
		}
		else
		{
			if(lam3>lamin && lam3<lamax)
			{
				lam = lam3;
			}
		}
	}

	if(lam >= 0.)
	{
		lam =0.;
	}

	return lam;

}


void AddRotationSymmetricFluxes(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q)
{
//Hier werden die konvektiven Anteile fuer ein rotationssymetrisches Gitter hinzugefuegt.
//Die entsprechenden viskosen Terme sind in CalcViscidFluxesXi/Eta enthalten

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				pnt_Q->Mass[ijk]-=pnt_mesh->jacobian[ijk]*pnt_U_RK->rho[ijk]*pnt_U_RK->v[ijk]/pnt_mesh->y[ijk];
				pnt_Q->xiMomentum[ijk]-=pnt_mesh->jacobian[ijk]*pnt_U_RK->rho[ijk]*pnt_U_RK->u[ijk]*pnt_U_RK->v[ijk]/pnt_mesh->y[ijk];
				pnt_Q->etaMomentum[ijk]-=pnt_mesh->jacobian[ijk]*pnt_U_RK->rho[ijk]*pnt_U_RK->v[ijk]*pnt_U_RK->v[ijk]/pnt_mesh->y[ijk];
				pnt_Q->Energy[ijk]-=pnt_mesh->jacobian[ijk]*pnt_U_RK->rho[ijk]*pnt_U_RK->v[ijk]*
						(pnt_U_RK->e[ijk]+pnt_config->Upsilon*pnt_U_RK->p[ijk]/pnt_U_RK->rho[ijk])
						/pnt_mesh->y[ijk];
			}
		}
	}
}

FLT IBC_getActualPosition(
		struct strct_configuration * pnt_config)
{
	FLT A,B,C,D,E,F,G,H,I,J,K,Y,X;
	X=pnt_config->time_dim;
	if((pnt_config->int_actualIteration%pnt_config->IBC_MovingStepsize==0)||
			(pnt_config->int_actualIteration==pnt_config->int_StartIteration))
	{
		if(pnt_config->IBC_MovingType==1)
		{
			X=X*1000.*pnt_config->IBC_SpeedFactor;

			A =  3.84021840481025777958E-01;
			B =  2.68347216478741421031E-02;
			C =  3.68769715145374510357E+00;
			D = -1.06913919412601088332E+00;
			E =  1.78847617010102349910E-01;
			F = -1.53841857081949747593E-02;
			G =  7.04907734706520833386E-04;
			H = -1.71741718123494413884E-05;
			I =  1.93590692040886198262E-07;
			J = -3.10882096431328766368E-10;
			K = -7.68466394126276780998E-12;

			Y = A + B*X + C*pow(X,2) + D*pow(X,3) + E*pow(X,4) + F*pow(X,5) + G*pow(X,6) + H*pow(X,7) + I*pow(X,8) + J*pow(X,9) + K*pow(X,10);
			return (pnt_config->IBC_StartpositionX+Y);
		}
		else if(pnt_config->IBC_MovingType==2)
		{
			X=X*1000.*pnt_config->IBC_SpeedFactor;

//			stopps within grid (leads to NAN)
//			if (X>14.0140142441)
//			{
//				Y=161.;
//			}

//			stopps outside of grid
			if (X>16.5)
			{
				Y=200.;
			}

			else
			{
				//new piston movement bases on EWM manuscript p.24
				Y = 100.*(1.-cos(X/16.5*MY_PI));
			}

			return (pnt_config->IBC_StartpositionX+Y);
		}
		else
		{
			Y=pnt_config->time_dim*pnt_config->IBC_MovingSpeed/pnt_config->L0_dim;

//			Y = 0.5*sin(X*1000);
//			if(X>((MY_PI/2)/1000.))
//			{
//			Y=0.5;
//			}

//			SVO_Validation
//			if(Y>0.5)
//			{
//				Y=0.5;
//			}
			return (pnt_config->IBC_StartpositionX+Y);
		}
	}
	else
	{
		return (pnt_config->IBC_MovingActualPosition);
	}
}

void freeVTA
(
		vta* x,vta* y,vta* z,vta* u,vta* v,vta* w,vta* rho,vta* p)
{
	if( x->ptr )
		free( x->ptr );
	if( y->ptr )
		free( y->ptr );
	if( z->ptr )
		free( z->ptr );
	if( u->ptr )
		free( u->ptr );
	if( v->ptr )
		free( v->ptr );
	if( w->ptr )
		free( w->ptr );
	if( p->ptr )
		free( p->ptr );
	if( rho->ptr )
		free( rho->ptr );
}


void check_Metric(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{
	int ijk,i,j,k;
	int flag_negjacobian;
	flag_negjacobian=0;
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				pnt_mesh->BC_Corrector_xiMomentum[ijk]=1.0;
				pnt_mesh->BC_Corrector_etaMomentum[ijk]=1.0;
				pnt_mesh->BC_Corrector_zetaMomentum[ijk]=1.0;
				if(pnt_mesh->jacobian[ijk]<0.0)
				{
					flag_negjacobian++;
//					pnt_mesh->jacobian[ijk]=-pnt_mesh->jacobian[ijk];
				}
			}
		}
	}
	if(flag_negjacobian>0)
	{
		j=pnt_config->int_jMid;
		j=pnt_config->int_jMid;
		k=pnt_config->int_kMid;
/*		printf("mid %d %d %d\n",pnt_config->int_iMid,pnt_config->int_jMid,pnt_config->int_kMid);
		printf("start real %d %d %d\n",pnt_config->int_iStartReal,pnt_config->int_jStartReal,pnt_config->int_kStartReal);
		printf("end real %d %d %d\n",pnt_config->int_iEndReal,pnt_config->int_jEndReal,pnt_config->int_kEndReal);
		printf("start ghosts %d %d %d\n",pnt_config->int_iStartGhosts,pnt_config->int_jStartGhosts,pnt_config->int_kStartGhosts);
		printf("end ghosts %d %d %d\n",pnt_config->int_iEndGhosts,pnt_config->int_jEndGhosts,pnt_config->int_kEndGhosts);
		ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k; */

		printf("ERROR: Rank %d (x_mid=%g y_mid=%g z_mid=%g) hat %d Gitterpunkte mit negativer Jacobian.\n",
				pnt_config->MPI_rank,
				(double)pnt_mesh->x[ijk],
				(double)pnt_mesh->y[ijk],
				(double)pnt_mesh->z[ijk],
				flag_negjacobian);
		printf("ERROR: Die Gruende sind meistens linksdrehendes Koordinatensystem oder sich ueberkreuzende GhostCells (keine Orthogonalitaet am Rand).\n");
	}
}

void CreateViscidMetric(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{
	int i,j,k,ijk,ijkMAX;

	ijkMAX=pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				pnt_mesh->xiFluss_Faktor[0*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(4.0/3.0*pow(pnt_mesh->xi_x[ijk],2.0)+pow(pnt_mesh->xi_y[ijk],2.0)+pow(pnt_mesh->xi_z[ijk],2.0));
				pnt_mesh->xiFluss_Faktor[1*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(4.0/3.0*pnt_mesh->eta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->eta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->eta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[2*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(4.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[3*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->xi_x[ijk]*pnt_mesh->xi_y[ijk];
				pnt_mesh->xiFluss_Faktor[4*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->eta_y[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->eta_x[ijk]*pnt_mesh->xi_y[ijk]);
				pnt_mesh->xiFluss_Faktor[5*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_y[ijk]);
				pnt_mesh->xiFluss_Faktor[6*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->xi_x[ijk]*pnt_mesh->xi_z[ijk];
				pnt_mesh->xiFluss_Faktor[7*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->eta_z[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->eta_x[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[8*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_z[ijk]);

				pnt_mesh->xiFluss_Faktor[9*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->xi_x[ijk]*pnt_mesh->xi_y[ijk];
				pnt_mesh->xiFluss_Faktor[10*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_y[ijk]*pnt_mesh->xi_x[ijk]-2.0/3.0* pnt_mesh->eta_x[ijk]*pnt_mesh->xi_y[ijk]);
				pnt_mesh->xiFluss_Faktor[11*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_x[ijk]-2.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_y[ijk]);
				pnt_mesh->xiFluss_Faktor[12*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pow(pnt_mesh->xi_x[ijk],2.0)+4.0/3.0*pow(pnt_mesh->xi_y[ijk],2.0)+pow(pnt_mesh->xi_z[ijk],2.0));
				pnt_mesh->xiFluss_Faktor[13*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_x[ijk]*pnt_mesh->xi_x[ijk]+4.0/3.0*pnt_mesh->eta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->eta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[14*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_x[ijk]+4.0/3.0*pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[15*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->xi_y[ijk]*pnt_mesh->xi_z[ijk];
				pnt_mesh->xiFluss_Faktor[16*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->eta_z[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->eta_y[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[17*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_z[ijk]);

				pnt_mesh->xiFluss_Faktor[18*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->xi_x[ijk]*pnt_mesh->xi_z[ijk];
				pnt_mesh->xiFluss_Faktor[19*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_z[ijk]*pnt_mesh->xi_x[ijk]-2.0/3.0*pnt_mesh->eta_x[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[20*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_x[ijk]-2.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[21*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->xi_y[ijk]*pnt_mesh->xi_z[ijk];
				pnt_mesh->xiFluss_Faktor[22*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_z[ijk]*pnt_mesh->xi_y[ijk]-2.0/3.0*pnt_mesh->eta_y[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[23*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_y[ijk]-2.0/3.0*pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[24*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pow(pnt_mesh->xi_x[ijk],2.0)+pow(pnt_mesh->xi_y[ijk],2.0)+4.0/3.0* pow(pnt_mesh->xi_z[ijk],2.0));
				pnt_mesh->xiFluss_Faktor[25*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->eta_y[ijk]*pnt_mesh->xi_y[ijk]+4.0/3.0*pnt_mesh->eta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[26*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_y[ijk]+4.0/3.0*pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_z[ijk]);

				pnt_mesh->xiFluss_Faktor[27*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pow(pnt_mesh->xi_x[ijk],2.0)+pow(pnt_mesh->xi_y[ijk],2.0)+pow(pnt_mesh->xi_z[ijk],2.0));
				pnt_mesh->xiFluss_Faktor[28*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->eta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->eta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->xiFluss_Faktor[29*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_z[ijk]);


				pnt_mesh->etaFluss_Faktor[0*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(4.0/3.0*pnt_mesh->eta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->eta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->eta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->etaFluss_Faktor[1*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(4.0/3.0* pow(pnt_mesh->eta_x[ijk],2.0)+pow(pnt_mesh->eta_y[ijk],2.0)+pow(pnt_mesh->eta_z[ijk],2.0));
				pnt_mesh->etaFluss_Faktor[2*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(4.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->etaFluss_Faktor[3*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_y[ijk]*pnt_mesh->xi_x[ijk]-2.0/3.0*pnt_mesh->eta_x[ijk]*pnt_mesh->xi_y[ijk]);
				pnt_mesh->etaFluss_Faktor[4*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->eta_x[ijk]*pnt_mesh->eta_y[ijk];
				pnt_mesh->etaFluss_Faktor[5*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_x[ijk]+pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_y[ijk]);
				pnt_mesh->etaFluss_Faktor[6*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_z[ijk]*pnt_mesh->xi_x[ijk]-2.0/3.0*pnt_mesh->eta_x[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->etaFluss_Faktor[7*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->eta_x[ijk]*pnt_mesh->eta_z[ijk];
				pnt_mesh->etaFluss_Faktor[8*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_x[ijk]+pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_z[ijk]);

				pnt_mesh->etaFluss_Faktor[9*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->eta_y[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->eta_x[ijk]*pnt_mesh->xi_y[ijk]);
				pnt_mesh->etaFluss_Faktor[10*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->eta_x[ijk]*pnt_mesh->eta_y[ijk];
				pnt_mesh->etaFluss_Faktor[11*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_x[ijk]-2.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_y[ijk]);
				pnt_mesh->etaFluss_Faktor[12*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_x[ijk]*pnt_mesh->xi_x[ijk]+4.0/3.0*pnt_mesh->eta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->eta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->etaFluss_Faktor[13*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pow(pnt_mesh->eta_x[ijk],2.0)+4.0/3.0*pow(pnt_mesh->eta_y[ijk],2.0)+pow(pnt_mesh->eta_z[ijk],2.0));
				pnt_mesh->etaFluss_Faktor[14*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_x[ijk]+4.0/3.0*pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->etaFluss_Faktor[15*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_z[ijk]*pnt_mesh->xi_y[ijk]-2.0/3.0*pnt_mesh->eta_y[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->etaFluss_Faktor[16*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->eta_y[ijk]*pnt_mesh->eta_z[ijk];
				pnt_mesh->etaFluss_Faktor[17*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_z[ijk]);

				pnt_mesh->etaFluss_Faktor[18*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->eta_z[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->eta_x[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->etaFluss_Faktor[19*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->eta_x[ijk]*pnt_mesh->eta_z[ijk];
				pnt_mesh->etaFluss_Faktor[20*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_x[ijk]-2.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->etaFluss_Faktor[21*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->eta_z[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->eta_y[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->etaFluss_Faktor[22*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->eta_y[ijk]*pnt_mesh->eta_z[ijk];
				pnt_mesh->etaFluss_Faktor[23*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_y[ijk]-2.0/3.0*pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->etaFluss_Faktor[24*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->eta_y[ijk]*pnt_mesh->xi_y[ijk]+4.0/3.0*pnt_mesh->eta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->etaFluss_Faktor[25*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pow(pnt_mesh->eta_x[ijk],2.0)+pow(pnt_mesh->eta_y[ijk],2.0)+4.0/3.0*pow(pnt_mesh->eta_z[ijk],2.0));
				pnt_mesh->etaFluss_Faktor[26*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_y[ijk]+4.0/3.0*pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_z[ijk]);

				pnt_mesh->etaFluss_Faktor[27*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->eta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->eta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->eta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->etaFluss_Faktor[28*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pow(pnt_mesh->eta_x[ijk],2.0)+pow(pnt_mesh->eta_y[ijk],2.0)+pow(pnt_mesh->eta_z[ijk],2.0));
				pnt_mesh->etaFluss_Faktor[29*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_z[ijk]);

				pnt_mesh->zetaFluss_Faktor[0*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(4.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[1*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(4.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[2*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(4.0/3.0*pow(pnt_mesh->zeta_x[ijk],2.0)+pow(pnt_mesh->zeta_y[ijk],2.0)+pow(pnt_mesh->zeta_z[ijk],2.0));
				pnt_mesh->zetaFluss_Faktor[3*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_x[ijk]-2.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_y[ijk]);
				pnt_mesh->zetaFluss_Faktor[4*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_x[ijk]-2.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_y[ijk]);
				pnt_mesh->zetaFluss_Faktor[5*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->zeta_x[ijk]*pnt_mesh->zeta_y[ijk];
				pnt_mesh->zetaFluss_Faktor[6*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_x[ijk]-2.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[7*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_x[ijk]-2.0/3.0*pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[8*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->zeta_x[ijk]*pnt_mesh->zeta_z[ijk];

				pnt_mesh->zetaFluss_Faktor[9*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_y[ijk]);
				pnt_mesh->zetaFluss_Faktor[10*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_x[ijk]+pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_y[ijk]);
				pnt_mesh->zetaFluss_Faktor[11*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->zeta_x[ijk]*pnt_mesh->zeta_y[ijk];
				pnt_mesh->zetaFluss_Faktor[12*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_x[ijk]+4.0/3.0*pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[13*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_x[ijk]+4.0/3.0*pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[14*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pow(pnt_mesh->zeta_x[ijk],2.0)+4.0/3.0*pow(pnt_mesh->zeta_y[ijk],2.0)+pow(pnt_mesh->zeta_z[ijk],2.0));
				pnt_mesh->zetaFluss_Faktor[15*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_y[ijk]-2.0/3.0*pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[16*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_y[ijk]-2.0/3.0*pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[17*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->zeta_y[ijk]*pnt_mesh->zeta_z[ijk];

				pnt_mesh->zetaFluss_Faktor[18*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[19*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_x[ijk]+pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[20*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->zeta_x[ijk]*pnt_mesh->zeta_z[ijk];
				pnt_mesh->zetaFluss_Faktor[21*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[22*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(-2.0/3.0* pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[23*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*1.0/3.0* pnt_mesh->zeta_y[ijk]*pnt_mesh->zeta_z[ijk];
				pnt_mesh->zetaFluss_Faktor[24*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_y[ijk]+4.0/3.0*pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[25*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_y[ijk]+4.0/3.0*pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[26*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pow(pnt_mesh->zeta_x[ijk],2.0)+pow(pnt_mesh->zeta_y[ijk],2.0)+4.0/3.0*pow(pnt_mesh->zeta_z[ijk],2.0));

				pnt_mesh->zetaFluss_Faktor[27*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->xi_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->xi_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->xi_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[28*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pnt_mesh->zeta_x[ijk]*pnt_mesh->eta_x[ijk]+pnt_mesh->zeta_y[ijk]*pnt_mesh->eta_y[ijk]+pnt_mesh->zeta_z[ijk]*pnt_mesh->eta_z[ijk]);
				pnt_mesh->zetaFluss_Faktor[29*ijkMAX+ijk]=pnt_mesh->jacobian[ijk]*(pow(pnt_mesh->zeta_x[ijk],2.0)+pow(pnt_mesh->zeta_y[ijk],2.0)+pow(pnt_mesh->zeta_z[ijk],2.0));
			}
		}
	}
	/*
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				int l,ijkl;
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				for(l=0;l<30;l++)
				{
					ijkl=l*ijkMAX+ijk;
					pnt_mesh->xiFluss_Faktor[ijkl]=1.0;
					pnt_mesh->etaFluss_Faktor[ijkl]=1.0;
					pnt_mesh->zetaFluss_Faktor[ijkl]=1.0;
				}
			}
		}
	}
	*/
}

void CreateMetric(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{
	FLT x_xi,x_eta,x_zeta;
	FLT y_xi,y_eta,y_zeta;
	FLT z_xi,z_eta,z_zeta;
	int i,j,k,ijk;
	int iPlus1jk,iMinus1jk;
	int ijPlus1k,ijMinus1k;
	int ijkPlus1,ijkMinus1;

	//	Erneute Übertragung des Gitters. Hier sind nun die Gitterinformationen der Nachbar-CPUs enthalten
	int ijk_extrapolate,i_extrapolate,j_extrapolate,k_extrapolate;
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		i_extrapolate=i+1;
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			j_extrapolate=j+1;
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				if(MESHDIMENSIONS==3)
				{
					k_extrapolate=k+1;
				}
				else
				{
					k_extrapolate=0;
				}
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijk_extrapolate=i_extrapolate*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+k_extrapolate;

				pnt_mesh->x_extrapolate[ijk_extrapolate]=pnt_mesh->x[ijk];
				pnt_mesh->y_extrapolate[ijk_extrapolate]=pnt_mesh->y[ijk];
				pnt_mesh->z_extrapolate[ijk_extrapolate]=pnt_mesh->z[ijk];

			}
		}
	}


	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		i_extrapolate=i+1;
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			j_extrapolate=j+1;
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				if(MESHDIMENSIONS==3)
				{
					k_extrapolate=k+1;
				}
				else
				{
					k_extrapolate=0;
				}
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				iPlus1jk=(i_extrapolate+1)*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+k_extrapolate;
				iMinus1jk=(i_extrapolate-1)*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+k_extrapolate;
				ijPlus1k=i_extrapolate*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+(j_extrapolate+1)*pnt_config->int_kMeshPointsGhostCells_extrapolate+k_extrapolate;
				ijMinus1k=i_extrapolate*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+(j_extrapolate-1)*pnt_config->int_kMeshPointsGhostCells_extrapolate+k_extrapolate;
				ijkPlus1=i_extrapolate*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+k_extrapolate+1;
				ijkMinus1=i_extrapolate*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+k_extrapolate-1;

				if(MESHDIMENSIONS==2)
				{
					x_xi=(-0.5*pnt_mesh->x_extrapolate[iMinus1jk]+0.5*pnt_mesh->x_extrapolate[iPlus1jk]);// /pnt_config->deltaXi;
					y_xi=(-0.5*pnt_mesh->y_extrapolate[iMinus1jk]+0.5*pnt_mesh->y_extrapolate[iPlus1jk]);// /pnt_config->deltaXi;
					z_xi=0.0;
					x_eta=(-0.5*pnt_mesh->x_extrapolate[ijMinus1k]+0.5*pnt_mesh->x_extrapolate[ijPlus1k]);// /pnt_config->deltaEta;
					y_eta=(-0.5*pnt_mesh->y_extrapolate[ijMinus1k]+0.5*pnt_mesh->y_extrapolate[ijPlus1k]);// /pnt_config->deltaEta;
					z_eta=0.0;
					x_zeta=0.0;
					y_zeta=0.0;
					z_zeta=1.0;
				}
				else
				{
					x_xi=(-0.5*pnt_mesh->x_extrapolate[iMinus1jk]+0.5*pnt_mesh->x_extrapolate[iPlus1jk]);// /pnt_config->deltaXi;
					y_xi=(-0.5*pnt_mesh->y_extrapolate[iMinus1jk]+0.5*pnt_mesh->y_extrapolate[iPlus1jk]);// /pnt_config->deltaXi;
					z_xi=(-0.5*pnt_mesh->z_extrapolate[iMinus1jk]+0.5*pnt_mesh->z_extrapolate[iPlus1jk]);// /pnt_config->deltaXi;
					x_eta=(-0.5*pnt_mesh->x_extrapolate[ijMinus1k]+0.5*pnt_mesh->x_extrapolate[ijPlus1k]);// /pnt_config->deltaEta;
					y_eta=(-0.5*pnt_mesh->y_extrapolate[ijMinus1k]+0.5*pnt_mesh->y_extrapolate[ijPlus1k]);// /pnt_config->deltaEta;
					z_eta=(-0.5*pnt_mesh->z_extrapolate[ijMinus1k]+0.5*pnt_mesh->z_extrapolate[ijPlus1k]);// /pnt_config->deltaEta;
					x_zeta=(-0.5*pnt_mesh->x_extrapolate[ijkMinus1]+0.5*pnt_mesh->x_extrapolate[ijkPlus1]);// /pnt_config->deltaZeta;
					y_zeta=(-0.5*pnt_mesh->y_extrapolate[ijkMinus1]+0.5*pnt_mesh->y_extrapolate[ijkPlus1]);// /pnt_config->deltaZeta;
					z_zeta=(-0.5*pnt_mesh->z_extrapolate[ijkMinus1]+0.5*pnt_mesh->z_extrapolate[ijkPlus1]);// /pnt_config->deltaZeta;

					y_xi=0.0;
					z_xi=0.0;
					x_eta=0.0;
					z_eta=0.0;
					x_zeta=0.0;
					y_zeta=0.0;

				}

				pnt_mesh->jacobian[ijk]=
						x_xi*y_eta*z_zeta+x_eta*y_zeta*z_xi+x_zeta*y_xi*z_eta
						-x_zeta*y_eta*z_xi-x_eta*y_xi*z_zeta-x_xi*y_zeta*z_eta;


				pnt_mesh->xi_x[ijk]=(y_eta*z_zeta-y_zeta*z_eta)/pnt_mesh->jacobian[ijk];
				pnt_mesh->xi_y[ijk]=(x_zeta*z_eta-x_eta*z_zeta)/pnt_mesh->jacobian[ijk];
				pnt_mesh->xi_z[ijk]=(x_eta*y_zeta-x_zeta*y_eta)/pnt_mesh->jacobian[ijk];
				pnt_mesh->eta_x[ijk]=(y_zeta*z_xi-y_xi*z_zeta)/pnt_mesh->jacobian[ijk];
				pnt_mesh->eta_y[ijk]=(x_xi*z_zeta-x_zeta*z_xi)/pnt_mesh->jacobian[ijk];
				pnt_mesh->eta_z[ijk]=(x_zeta*y_xi-x_xi*y_zeta)/pnt_mesh->jacobian[ijk];
				pnt_mesh->zeta_x[ijk]=(y_xi*z_eta-y_eta*z_xi)/pnt_mesh->jacobian[ijk];
				pnt_mesh->zeta_y[ijk]=(x_eta*z_xi-x_xi*z_eta)/pnt_mesh->jacobian[ijk];
				pnt_mesh->zeta_z[ijk]=(x_xi*y_eta-x_eta*y_xi)/pnt_mesh->jacobian[ijk];
			}
		}
	}

	free(pnt_mesh->x_extrapolate);
	free(pnt_mesh->y_extrapolate);
	free(pnt_mesh->z_extrapolate);
}

void ExtrapolateGhostCells(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{

//	######################
//	EXTRAPOLATION DER GHOSTCELLS
//	######################
	int i,j,k,ijk,ijkSymmetry,ijkSymmetry2;
	int int_symmetryIndex;
	int iMinus1jk,iPlus1jk,ijMinus1k,ijPlus1k;
#if MESHDIMENSIONS==3
	int ijkMinus1,ijkPlus1;
#endif

//Left
	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iPlus1jk=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry2=(int_symmetryIndex+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				pnt_mesh->x[ijk]=pnt_mesh->x[iPlus1jk]+pnt_mesh->x[ijkSymmetry]-pnt_mesh->x[ijkSymmetry2];
				pnt_mesh->y[ijk]=pnt_mesh->y[iPlus1jk]+pnt_mesh->y[ijkSymmetry]-pnt_mesh->y[ijkSymmetry2];
#if MESHDIMENSIONS==3
				pnt_mesh->z[ijk]=pnt_mesh->z[iPlus1jk]+pnt_mesh->z[ijkSymmetry]-pnt_mesh->z[ijkSymmetry2];
#endif
			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}

//Right
	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry2=(int_symmetryIndex-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;


				pnt_mesh->x[ijk]=pnt_mesh->x[iMinus1jk]+pnt_mesh->x[ijkSymmetry]-pnt_mesh->x[ijkSymmetry2];
				pnt_mesh->y[ijk]=pnt_mesh->y[iMinus1jk]+pnt_mesh->y[ijkSymmetry]-pnt_mesh->y[ijkSymmetry2];
#if MESHDIMENSIONS==3
				pnt_mesh->z[ijk]=pnt_mesh->z[iMinus1jk]+pnt_mesh->z[ijkSymmetry]-pnt_mesh->z[ijkSymmetry2];
#endif


			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}

//Bottom
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(int_symmetryIndex+1)*pnt_config->int_kMeshPointsGhostCells+k;

				pnt_mesh->x[ijk]=pnt_mesh->x[ijPlus1k]+pnt_mesh->x[ijkSymmetry]-pnt_mesh->x[ijkSymmetry2];
				pnt_mesh->y[ijk]=pnt_mesh->y[ijPlus1k]+pnt_mesh->y[ijkSymmetry]-pnt_mesh->y[ijkSymmetry2];
#if MESHDIMENSIONS==3
				pnt_mesh->z[ijk]=pnt_mesh->z[ijPlus1k]+pnt_mesh->z[ijkSymmetry]-pnt_mesh->z[ijkSymmetry2];
#endif
			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}

//Top
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(int_symmetryIndex-1)*pnt_config->int_kMeshPointsGhostCells+k;

				pnt_mesh->x[ijk]=pnt_mesh->x[ijMinus1k]+pnt_mesh->x[ijkSymmetry]-pnt_mesh->x[ijkSymmetry2];
				pnt_mesh->y[ijk]=pnt_mesh->y[ijMinus1k]+pnt_mesh->y[ijkSymmetry]-pnt_mesh->y[ijkSymmetry2];
#if MESHDIMENSIONS==3
				pnt_mesh->z[ijk]=pnt_mesh->z[ijMinus1k]+pnt_mesh->z[ijkSymmetry]-pnt_mesh->z[ijkSymmetry2];
#endif

			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}

#if MESHDIMENSIONS==3
//Behind
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkPlus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+1;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;
				ijkSymmetry2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex+1;

				pnt_mesh->x[ijk]=pnt_mesh->x[ijkPlus1]+pnt_mesh->x[ijkSymmetry]-pnt_mesh->x[ijkSymmetry2];
				pnt_mesh->y[ijk]=pnt_mesh->y[ijkPlus1]+pnt_mesh->y[ijkSymmetry]-pnt_mesh->y[ijkSymmetry2];
				pnt_mesh->z[ijk]=pnt_mesh->z[ijkPlus1]+pnt_mesh->z[ijkSymmetry]-pnt_mesh->z[ijkSymmetry2];

				int_symmetryIndex=int_symmetryIndex+1;
			}
		}
	}

//InFront
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkMinus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-1;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;
				ijkSymmetry2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex-1;

				pnt_mesh->x[ijk]=pnt_mesh->x[ijkMinus1]+pnt_mesh->x[ijkSymmetry]-pnt_mesh->x[ijkSymmetry2];
				pnt_mesh->y[ijk]=pnt_mesh->y[ijkMinus1]+pnt_mesh->y[ijkSymmetry]-pnt_mesh->y[ijkSymmetry2];
				pnt_mesh->z[ijk]=pnt_mesh->z[ijkMinus1]+pnt_mesh->z[ijkSymmetry]-pnt_mesh->z[ijkSymmetry2];


				int_symmetryIndex=int_symmetryIndex-1;
			}
		}
	}
#endif

	//	######################
	//	EXTRAPOLATION DER GHOSTCELLS FÜR DAS ERWEITERTE GITTER
	//	######################

	//	Zunächst Übertragung des Gitters. Hier ist erstmal nur das reale Gitter enthalten
	int ijk_extrapolate,i_extrapolate,j_extrapolate,k_extrapolate;
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		i_extrapolate=i+1;
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			j_extrapolate=j+1;
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				if(MESHDIMENSIONS==3)
				{
					k_extrapolate=k+1;
				}
				else
				{
					k_extrapolate=0;
				}
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijk_extrapolate=i_extrapolate*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+k_extrapolate;

				pnt_mesh->x_extrapolate[ijk_extrapolate]=pnt_mesh->x[ijk];
				pnt_mesh->y_extrapolate[ijk_extrapolate]=pnt_mesh->y[ijk];
#if MESHDIMENSIONS==3
				pnt_mesh->z_extrapolate[ijk_extrapolate]=pnt_mesh->z[ijk];
#endif
			}
		}
	}

//Left
	int_symmetryIndex=pnt_config->int_iStartReal_extrapolate;
	for (i=pnt_config->int_iStartReal_extrapolate-1; i >= pnt_config->int_iStartGhosts_extrapolate; i--)
	{
		for (j=pnt_config->int_jStartReal_extrapolate; j <= pnt_config->int_jEndReal_extrapolate; j++)
		{
			for (k=pnt_config->int_kStartReal_extrapolate; k <= pnt_config->int_kEndReal_extrapolate; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				iPlus1jk=(i+1)*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijkSymmetry2=(int_symmetryIndex+1)*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;

				pnt_mesh->x_extrapolate[ijk]=pnt_mesh->x_extrapolate[iPlus1jk]+pnt_mesh->x_extrapolate[ijkSymmetry]-pnt_mesh->x_extrapolate[ijkSymmetry2];
				pnt_mesh->y_extrapolate[ijk]=pnt_mesh->y_extrapolate[iPlus1jk]+pnt_mesh->y_extrapolate[ijkSymmetry]-pnt_mesh->y_extrapolate[ijkSymmetry2];
#if MESHDIMENSIONS==3
				pnt_mesh->z_extrapolate[ijk]=pnt_mesh->z_extrapolate[iPlus1jk]+pnt_mesh->z_extrapolate[ijkSymmetry]-pnt_mesh->z_extrapolate[ijkSymmetry2];
#endif
			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}

//Right
	int_symmetryIndex=pnt_config->int_iEndReal_extrapolate;
	for (i=pnt_config->int_iEndReal_extrapolate+1; i <= pnt_config->int_iEndGhosts_extrapolate; i++)
	{
		for (j=pnt_config->int_jStartReal_extrapolate; j <= pnt_config->int_jEndReal_extrapolate; j++)
		{
			for (k=pnt_config->int_kStartReal_extrapolate; k <= pnt_config->int_kEndReal_extrapolate; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijkSymmetry2=(int_symmetryIndex-1)*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;


				pnt_mesh->x_extrapolate[ijk]=pnt_mesh->x_extrapolate[iMinus1jk]+pnt_mesh->x_extrapolate[ijkSymmetry]-pnt_mesh->x_extrapolate[ijkSymmetry2];
				pnt_mesh->y_extrapolate[ijk]=pnt_mesh->y_extrapolate[iMinus1jk]+pnt_mesh->y_extrapolate[ijkSymmetry]-pnt_mesh->y_extrapolate[ijkSymmetry2];
#if MESHDIMENSIONS==3
				pnt_mesh->z_extrapolate[ijk]=pnt_mesh->z_extrapolate[iMinus1jk]+pnt_mesh->z_extrapolate[ijkSymmetry]-pnt_mesh->z_extrapolate[ijkSymmetry2];
#endif

			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}

//Bottom
	for (i=pnt_config->int_iStartGhosts_extrapolate; i <= pnt_config->int_iEndGhosts_extrapolate; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal_extrapolate;
		for (j=pnt_config->int_jStartReal_extrapolate-1; j >= pnt_config->int_jStartGhosts_extrapolate; j--)
		{
			for (k=pnt_config->int_kStartReal_extrapolate; k <= pnt_config->int_kEndReal_extrapolate; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+(j+1)*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijkSymmetry2=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+(int_symmetryIndex+1)*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;

				pnt_mesh->x_extrapolate[ijk]=pnt_mesh->x_extrapolate[ijPlus1k]+pnt_mesh->x_extrapolate[ijkSymmetry]-pnt_mesh->x_extrapolate[ijkSymmetry2];
				pnt_mesh->y_extrapolate[ijk]=pnt_mesh->y_extrapolate[ijPlus1k]+pnt_mesh->y_extrapolate[ijkSymmetry]-pnt_mesh->y_extrapolate[ijkSymmetry2];
#if MESHDIMENSIONS==3
				pnt_mesh->z_extrapolate[ijk]=pnt_mesh->z_extrapolate[ijPlus1k]+pnt_mesh->z_extrapolate[ijkSymmetry]-pnt_mesh->z_extrapolate[ijkSymmetry2];
#endif

			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}

//Top
	for (i=pnt_config->int_iStartGhosts_extrapolate; i <= pnt_config->int_iEndGhosts_extrapolate; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal_extrapolate;
		for (j=pnt_config->int_jEndReal_extrapolate+1; j <= pnt_config->int_jEndGhosts_extrapolate; j++)
		{
			for (k=pnt_config->int_kStartReal_extrapolate; k <= pnt_config->int_kEndReal_extrapolate; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+(j-1)*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijkSymmetry2=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+(int_symmetryIndex-1)*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;

				pnt_mesh->x_extrapolate[ijk]=pnt_mesh->x_extrapolate[ijMinus1k]+pnt_mesh->x_extrapolate[ijkSymmetry]-pnt_mesh->x_extrapolate[ijkSymmetry2];
				pnt_mesh->y_extrapolate[ijk]=pnt_mesh->y_extrapolate[ijMinus1k]+pnt_mesh->y_extrapolate[ijkSymmetry]-pnt_mesh->y_extrapolate[ijkSymmetry2];
#if MESHDIMENSIONS==3
				pnt_mesh->z_extrapolate[ijk]=pnt_mesh->z_extrapolate[ijMinus1k]+pnt_mesh->z_extrapolate[ijkSymmetry]-pnt_mesh->z_extrapolate[ijkSymmetry2];
#endif

			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}

#if MESHDIMENSIONS==3
//Behind
	for (i=pnt_config->int_iStartGhosts_extrapolate; i <= pnt_config->int_iEndGhosts_extrapolate; i++)
	{
		for (j=pnt_config->int_jStartGhosts_extrapolate; j <= pnt_config->int_jEndGhosts_extrapolate; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal_extrapolate;
			for (k=pnt_config->int_kStartReal_extrapolate-1; k >= pnt_config->int_kStartGhosts_extrapolate; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijkPlus1=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k+1;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+int_symmetryIndex;
				ijkSymmetry2=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+int_symmetryIndex+1;

				pnt_mesh->x_extrapolate[ijk]=pnt_mesh->x_extrapolate[ijkPlus1]+pnt_mesh->x_extrapolate[ijkSymmetry]-pnt_mesh->x_extrapolate[ijkSymmetry2];
				pnt_mesh->y_extrapolate[ijk]=pnt_mesh->y_extrapolate[ijkPlus1]+pnt_mesh->y_extrapolate[ijkSymmetry]-pnt_mesh->y_extrapolate[ijkSymmetry2];
				pnt_mesh->z_extrapolate[ijk]=pnt_mesh->z_extrapolate[ijkPlus1]+pnt_mesh->z_extrapolate[ijkSymmetry]-pnt_mesh->z_extrapolate[ijkSymmetry2];

				int_symmetryIndex=int_symmetryIndex+1;
			}
		}
	}

//InFront
	for (i=pnt_config->int_iStartGhosts_extrapolate; i <= pnt_config->int_iEndGhosts_extrapolate; i++)
	{
		for (j=pnt_config->int_jStartGhosts_extrapolate; j <= pnt_config->int_jEndGhosts_extrapolate; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal_extrapolate;
			for (k=pnt_config->int_kEndReal_extrapolate+1; k <= pnt_config->int_kEndGhosts_extrapolate; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k;
				ijkMinus1=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+k-1;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+int_symmetryIndex;
				ijkSymmetry2=i*pnt_config->int_jMeshPointsGhostCells_extrapolate*pnt_config->int_kMeshPointsGhostCells_extrapolate+j*pnt_config->int_kMeshPointsGhostCells_extrapolate+int_symmetryIndex-1;

				pnt_mesh->x_extrapolate[ijk]=pnt_mesh->x_extrapolate[ijkMinus1]+pnt_mesh->x_extrapolate[ijkSymmetry]-pnt_mesh->x_extrapolate[ijkSymmetry2];
				pnt_mesh->y_extrapolate[ijk]=pnt_mesh->y_extrapolate[ijkMinus1]+pnt_mesh->y_extrapolate[ijkSymmetry]-pnt_mesh->y_extrapolate[ijkSymmetry2];
				pnt_mesh->z_extrapolate[ijk]=pnt_mesh->z_extrapolate[ijkMinus1]+pnt_mesh->z_extrapolate[ijkSymmetry]-pnt_mesh->z_extrapolate[ijkSymmetry2];


				int_symmetryIndex=int_symmetryIndex-1;
			}
		}
	}
#endif
}

void mirrorGhostCellsMetric(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;

	if((strcmp(pnt_config->BC_Left,pnt_config->BCWallInviscid)==0)
		||
		(strcmp(pnt_config->BC_Left,pnt_config->BCWallViscous)==0)
		||
		(strcmp(pnt_config->BC_Left,pnt_config->BCFarfield)==0)
		||
		(strcmp(pnt_config->BC_Left,pnt_config->BCWallViscousIsothermal)==0))
	{
	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				copyMetric(pnt_config,pnt_mesh,ijk,ijkSymmetry);
			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
	}

	if((strcmp(pnt_config->BC_Right,pnt_config->BCWallInviscid)==0)
		||
		(strcmp(pnt_config->BC_Right,pnt_config->BCWallViscous)==0)
		||
		(strcmp(pnt_config->BC_Right,pnt_config->BCFarfield)==0)
		||
		(strcmp(pnt_config->BC_Right,pnt_config->BCWallViscousIsothermal)==0))
	{
	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				copyMetric(pnt_config,pnt_mesh,ijk,ijkSymmetry);
			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
	}

	if((strcmp(pnt_config->BC_Bottom,pnt_config->BCWallInviscid)==0)
		||
		(strcmp(pnt_config->BC_Bottom,pnt_config->BCWallViscous)==0)
		||
		(strcmp(pnt_config->BC_Bottom,pnt_config->BCFarfield)==0)
		||
		(strcmp(pnt_config->BC_Bottom,pnt_config->BCWallViscousIsothermal)==0))
	{
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				copyMetric(pnt_config,pnt_mesh,ijk,ijkSymmetry);
			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
	}



	if((strcmp(pnt_config->BC_Top,pnt_config->BCWallInviscid)==0)
		||
		(strcmp(pnt_config->BC_Top,pnt_config->BCWallViscous)==0)
		||
		(strcmp(pnt_config->BC_Top,pnt_config->BCFarfield)==0)
		||
		(strcmp(pnt_config->BC_Top,pnt_config->BCWallViscousIsothermal)==0))		
	{
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				copyMetric(pnt_config,pnt_mesh,ijk,ijkSymmetry);
			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
	}



	if((strcmp(pnt_config->BC_Behind,pnt_config->BCWallInviscid)==0)
		||
		(strcmp(pnt_config->BC_Behind,pnt_config->BCWallViscous)==0)
		||
		(strcmp(pnt_config->BC_Behind,pnt_config->BCFarfield)==0)
		||
		(strcmp(pnt_config->BC_Behind,pnt_config->BCWallViscousIsothermal)==0))		
	{
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				copyMetric(pnt_config,pnt_mesh,ijk,ijkSymmetry);

				int_symmetryIndex=int_symmetryIndex+1;
			}
		}
	}
	}


	if((strcmp(pnt_config->BC_InFront,pnt_config->BCWallInviscid)==0)
		||
		(strcmp(pnt_config->BC_InFront,pnt_config->BCWallViscous)==0)
		||
		(strcmp(pnt_config->BC_InFront,pnt_config->BCFarfield)==0)
		||
		(strcmp(pnt_config->BC_InFront,pnt_config->BCWallViscousIsothermal)==0))		
	{
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				copyMetric(pnt_config,pnt_mesh,ijk,ijkSymmetry);

				int_symmetryIndex=int_symmetryIndex-1;
			}
		}
	}
	}

}

void copyMetric(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		int ijkDestination,
		int ijkSource)
{
	pnt_mesh->jacobian[ijkDestination]=pnt_mesh->jacobian[ijkSource];
	pnt_mesh->xi_x[ijkDestination]=pnt_mesh->xi_x[ijkSource];
	pnt_mesh->xi_y[ijkDestination]=pnt_mesh->xi_y[ijkSource];
	pnt_mesh->xi_z[ijkDestination]=pnt_mesh->xi_z[ijkSource];
	pnt_mesh->eta_x[ijkDestination]=pnt_mesh->eta_x[ijkSource];
	pnt_mesh->eta_y[ijkDestination]=pnt_mesh->eta_y[ijkSource];
	pnt_mesh->eta_z[ijkDestination]=pnt_mesh->eta_z[ijkSource];
	pnt_mesh->zeta_x[ijkDestination]=pnt_mesh->zeta_x[ijkSource];
	pnt_mesh->zeta_y[ijkDestination]=pnt_mesh->zeta_y[ijkSource];
	pnt_mesh->zeta_z[ijkDestination]=pnt_mesh->zeta_z[ijkSource];
}


/*----------------------------------------------------------------------------
! This file is part of psOpen
!
! Version 0.9
!
! Copyright (C) 2013 Jens Henrik Goebbert <jens.henrik.goebbert()rwth-aachen.de>
!
!    psOpen is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!    Please consider mentioning the copyright owner in any publication
!    if results of psOpen are used.
!----------------------------------------------------------------------------*/

#ifdef IBM

#include <spi/include/kernel/memory.h>

extern void print_memusage_c()
{
   printf("SHOCK: IBM memory usage per rank\n");
  uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;

  Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);

  printf("SHOCK: current MEMSIZE heap  : %.2f/%.2f stack: %.2f/%.2f mmap: %.2f mbyte\n", (FLT.heap/(1024*1024), (FLT.heapavail/(1024*1024),
                                                                              (FLT)stack/(1024*1024), (FLT)stackavail/(1024*1024),
                                                                              (FLT)mmap/(1024*1024));
  printf("SHOCK: current MEMSIZE shared: %.2f persist: %.2f guard: %.2f mbyte\n", (FLT)shared/(1024*1024),
                                                                       (FLT)persist/(1024*1024),
                                                                        (FLT)guard/(1024*1024));

}
#elif INTEL

#include <unistd.h> // for sysconf(_SC_PAGESIZE)
#include <linux/limits.h> // for PATH_MAX

void readone(FILE *f, long long int *x)            { if(fscanf(f, "%lld ", x)==EOF){printf("Error in fscanf\n");} }
void readunsigned(FILE *f, unsigned long long *x) { if(fscanf(f, "%llu ", x)==EOF){printf("Error in fscanf\n");} }
void readstr(FILE *f, char *x)                      { if(fscanf(f, "%s ", x)==EOF){printf("Error in fscanf\n");} }
void readchar(FILE *f, char *x)                     { if(fscanf(f, "%c ", x)==EOF){printf("Error in fscanf\n");} }

extern void print_memusage_c()
{

  printf("SHOCK: Intel memory usage per rank\n");
  long page_size = sysconf(_SC_PAGESIZE);
  //long s = -1;
  FILE *f = fopen("/proc/self/stat", "r");
  if (!f) return ;

  // example: 24773 (cat) R
  long long int pid;    readone(f,&pid);    // process id
  char tcomm[PATH_MAX];  readstr(f,tcomm);   // process name
  char state;            readchar(f,&state); // process status (R==running, ...)

  // example: 7627 24773 7627 34827 24773
  long long int ppid;     readone(f,&ppid);
  long long int pgid;     readone(f,&pgid);
  long long int sid;      readone(f,&sid);
  long long int tty_nr;   readone(f,&tty_nr);
  long long int tty_pgrp; readone(f,&tty_pgrp);

  // example: 4202496 220 0 0 0 0 0
  long long int flags;    readone(f,&flags);
  long long int min_dbl;  readone(f,&min_dbl);
  long long int cmin_dbl; readone(f,&cmin_dbl);
  long long int maj_dbl;  readone(f,&maj_dbl);
  long long int cmaj_dbl; readone(f,&cmaj_dbl);
  long long int utime;    readone(f,&utime);
  long long int stimev;   readone(f,&stimev);

  // example: 0 0 20 0 1 0
  long long int cutime;        readone(f,&cutime);
  long long int cstime;        readone(f,&cstime);
  long long int priority;      readone(f,&priority);    // process priority
  long long int nicev;         readone(f,&nicev);       // process nice value
  long long int num_threads;   readone(f,&num_threads); // no. process threads
  long long int it_real_value; readone(f,&it_real_value);

  // example: 29350452
  unsigned long long start_time; readunsigned(f,&start_time); // start time as UNIX time

  // example: 4210688 123 18446744073709551615 4194304 4235780 140734041535920 140734041533112 215179179504
  long long int vsize;      readone(f,&vsize);                                // virtual memory size (in bytes)
  long long int rss;        readone(f,&rss);        rss   = rss   *page_size; // resident set size (rss) is the portion of a process's memory that is held in RAM.
  long long int rsslim;     readone(f,&rsslim);     rsslim= rsslim*page_size; // limit for resident set size (rss) (in pages)
  long long int start_code; readone(f,&start_code);
  long long int end_code;   readone(f,&end_code);
  long long int start_stack;readone(f,&start_stack);
  long long int esp;        readone(f,&esp);
  long long int eip;        readone(f,&eip);

  // example: 0 0 0 0 0 0 0 17 3 0 0
  long long int pending;     readone(f,&pending);
  long long int blocked;     readone(f,&blocked);
  long long int sigign;      readone(f,&sigign);
  long long int sigcatch;    readone(f,&sigcatch);
  long long int wchan;       readone(f,&wchan);
  long long int zero1;       readone(f,&zero1);
  long long int zero2;       readone(f,&zero2);
  long long int exit_signal; readone(f,&exit_signal);
  long long int cpu;         readone(f,&cpu);
  long long int rt_priority; readone(f,&rt_priority);
  long long int policy;      readone(f,&policy);

  // example: 0 0 0
  // unknown

  fclose (f);

  printf("SHOCK: current MEMSIZE RSS  : %.2f mbyte\n", (double)rss/(1024*1024));
  printf("SHOCK: current MEMSIZE VSIZE: %.2f mbyte\n", (double)vsize/(1024*1024));
}
#else

extern void print_memusage_c()
{
printf("SHOCK: current MEMSIZE RSS  : n/a\n");
printf("SHOCK: current MEMSIZE VSIZE: n/a\n");
}
#endif
