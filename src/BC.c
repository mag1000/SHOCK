#include "BC.h"
#include "SHOCK.h"
#include "WENO.h"
#include "ZD.h"
#include "Functions.h"
#include "ManufacturedSolution.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "string.h"

void SetAllBoundaryConditions(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	if (pnt_config->InterfaceNeighbourLeft==NO_NEIGHBOUR)
	{
		ControlBoundaryLowerI(
				pnt_config,
				pnt_mesh,
				pnt_U,
				pnt_U_lastStep);
	}
	if (pnt_config->InterfaceNeighbourRight==NO_NEIGHBOUR)
	{
		ControlBoundaryUpperI(
				pnt_config,
				pnt_mesh,
				pnt_U,
				pnt_U_lastStep);
	}

	if (pnt_config->InterfaceNeighbourBottom==NO_NEIGHBOUR)
	{
		ControlBoundaryLowerJ(
				pnt_config,
				pnt_mesh,
				pnt_U,
				pnt_U_lastStep);
	}
	if (pnt_config->InterfaceNeighbourTop==NO_NEIGHBOUR)
	{
		ControlBoundaryUpperJ(
				pnt_config,
				pnt_mesh,
				pnt_U,
				pnt_U_lastStep);
	}


	if(MESHDIMENSIONS==3)
	{

		if (pnt_config->InterfaceNeighbourBehind==NO_NEIGHBOUR)
		{
			ControlBoundaryLowerK(
					pnt_config,
					pnt_mesh,
					pnt_U,
					pnt_U_lastStep);
		}
		if (pnt_config->InterfaceNeighbourInFront==NO_NEIGHBOUR)
		{
			ControlBoundaryUpperK(
					pnt_config,
					pnt_mesh,
					pnt_U,
					pnt_U_lastStep);
		}
	}
}

void ControlBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{

	if(strcmp(pnt_config->BC_Left,pnt_config->BCWallViscous)==0)
	{
		WriteWallNoSlipBoundaryLowerI(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Left,pnt_config->BCWallInviscid)==0)
	{
		WriteWallSlipBoundaryLowerI(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if((strcmp(pnt_config->BC_Left,pnt_config->BCInflow)==0)||
			(strcmp(pnt_config->BC_Left,pnt_config->BCInflowSubsonic)==0)||
			(strcmp(pnt_config->BC_Left,pnt_config->BCInflowSupersonic)==0))
	{
		if(pnt_config->flag_BC_option_inflow_normal_sub)
			WriteInflowSubsonicNormalBoundaryLowerI(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_riemann_sub)
			WriteInflowSubsonicRiemannBoundaryLowerI(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_isentrop_sub)
			WriteInflowSubsonicIsentropBoundaryLowerI(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_normal_super)
			WriteInflowSupersonicNormalBoundaryLowerI(
					pnt_config,
					pnt_mesh,
					pnt_U);
	}
	else if((strcmp(pnt_config->BC_Left,pnt_config->BCOutflow)==0)||
			(strcmp(pnt_config->BC_Left,pnt_config->BCOutflowSubsonic)==0))
	{
		if(pnt_config->flag_BC_option_outflow_normal_sub)
			WriteOutflowSubsonicNormalBoundaryLowerI(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_riemann_sub)
			WriteOutflowSubsonicRiemannBoundaryLowerI(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_rudy_sub)
			WriteOutflowSubsonicRudyBoundaryLowerI(
					pnt_config,
					pnt_mesh,
					pnt_U,
					pnt_U_lastStep);
	}

	else if(strcmp(pnt_config->BC_Left,pnt_config->BCFarfield)==0)
	{
		WriteFarfieldBoundaryLowerI(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Left,pnt_config->BCWallViscousIsothermal)==0)
	{
		WriteWallNoSlipIsothermalBoundaryLowerI(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Left,pnt_config->BCManufacturedSolution)==0)
	{
		WriteBCManufacturedSolutionLowerI(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else
	{
		printf("SHOCK: Error: Rank %d hat links weder Randbedingung noch Nachbar!\n",pnt_config->MPI_rank);
	}
}

void ControlBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	if(strcmp(pnt_config->BC_Bottom,pnt_config->BCWallViscous)==0)
	{
		WriteWallNoSlipBoundaryLowerJ(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Bottom,pnt_config->BCWallInviscid)==0)
	{
		WriteWallSlipBoundaryLowerJ(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if((strcmp(pnt_config->BC_Bottom,pnt_config->BCInflow)==0)||
			(strcmp(pnt_config->BC_Bottom,pnt_config->BCInflowSubsonic)==0)||
			(strcmp(pnt_config->BC_Bottom,pnt_config->BCInflowSupersonic)==0))
	{
		if(pnt_config->flag_BC_option_inflow_normal_sub)
			WriteInflowSubsonicNormalBoundaryLowerJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_riemann_sub)
			WriteInflowSubsonicRiemannBoundaryLowerJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_isentrop_sub)
			WriteInflowSubsonicIsentropBoundaryLowerJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_normal_super)
			WriteInflowSupersonicNormalBoundaryLowerJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
	}
	else if((strcmp(pnt_config->BC_Bottom,pnt_config->BCOutflow)==0)||
			(strcmp(pnt_config->BC_Bottom,pnt_config->BCOutflowSubsonic)==0))
	{
		if(pnt_config->flag_BC_option_outflow_normal_sub)
			WriteOutflowSubsonicNormalBoundaryLowerJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_riemann_sub)
			WriteOutflowSubsonicRiemannBoundaryLowerJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_rudy_sub)
			WriteOutflowSubsonicRudyBoundaryLowerJ(
					pnt_config,
					pnt_mesh,
					pnt_U,
					pnt_U_lastStep);
	}
	else if(strcmp(pnt_config->BC_Bottom,pnt_config->BCFarfield)==0)
	{
		WriteFarfieldBoundaryLowerJ(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Bottom,pnt_config->BCWallViscousIsothermal)==0)
	{
		WriteWallNoSlipIsothermalBoundaryLowerJ(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Bottom,pnt_config->BCManufacturedSolution)==0)
	{
		WriteBCManufacturedSolutionLowerJ(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else
	{
		printf("SHOCK: Error: Rank %d hat unten weder Randbedingung noch Nachbar!\n",pnt_config->MPI_rank);
	}
//	else if(strcmp(pnt_config->BC_Bottom,pnt_config->BCInflowSupersonic)==0)
//	{
//		WriteInflowSupersonicNormalBoundaryLowerJ(
//				pnt_config,
//				pnt_mesh,
//				pnt_U,
//				pnt_ZD);
//	}
}

void ControlBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	if(strcmp(pnt_config->BC_Behind,pnt_config->BCWallViscous)==0)
	{
		WriteWallNoSlipBoundaryLowerK(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Behind,pnt_config->BCWallInviscid)==0)
	{
		WriteWallSlipBoundaryLowerK(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if((strcmp(pnt_config->BC_Behind,pnt_config->BCInflow)==0)||
			(strcmp(pnt_config->BC_Behind,pnt_config->BCInflowSubsonic)==0)||
			(strcmp(pnt_config->BC_Behind,pnt_config->BCInflowSupersonic)==0))
	{
		if(pnt_config->flag_BC_option_inflow_normal_sub)
			WriteInflowSubsonicNormalBoundaryLowerK(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_riemann_sub)
			WriteInflowSubsonicRiemannBoundaryLowerK(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_isentrop_sub)
			WriteInflowSubsonicIsentropBoundaryLowerK(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_normal_super)
			WriteInflowSupersonicNormalBoundaryLowerK(
					pnt_config,
					pnt_mesh,
					pnt_U);
	}
	else if((strcmp(pnt_config->BC_Behind,pnt_config->BCOutflow)==0)||
			(strcmp(pnt_config->BC_Behind,pnt_config->BCOutflowSubsonic)==0))
	{
		if(pnt_config->flag_BC_option_outflow_normal_sub)
			WriteOutflowSubsonicNormalBoundaryLowerK(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_riemann_sub)
			WriteOutflowSubsonicRiemannBoundaryLowerK(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_rudy_sub)
			WriteOutflowSubsonicRudyBoundaryLowerK(
					pnt_config,
					pnt_mesh,
					pnt_U,
					pnt_U_lastStep);
	}
	else if(strcmp(pnt_config->BC_Behind,pnt_config->BCFarfield)==0)
	{
		WriteFarfieldBoundaryLowerK(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Behind,pnt_config->BCWallViscousIsothermal)==0)
	{
		WriteWallNoSlipIsothermalBoundaryLowerK(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Behind,pnt_config->BCManufacturedSolution)==0)
	{
		WriteBCManufacturedSolutionLowerK(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else
	{
		printf("SHOCK: Error: Rank %d hat hinten weder Randbedingung noch Nachbar!\n",pnt_config->MPI_rank);
	}
}

void ControlBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	if(strcmp(pnt_config->BC_Right,pnt_config->BCWallViscous)==0)
	{
		WriteWallNoSlipBoundaryUpperI(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Right,pnt_config->BCWallInviscid)==0)
	{
		WriteWallSlipBoundaryUpperI(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if((strcmp(pnt_config->BC_Right,pnt_config->BCInflow)==0)||
			(strcmp(pnt_config->BC_Right,pnt_config->BCInflowSubsonic)==0)||
			(strcmp(pnt_config->BC_Right,pnt_config->BCInflowSupersonic)==0))
	{
		if(pnt_config->flag_BC_option_inflow_normal_sub)
			WriteInflowSubsonicNormalBoundaryUpperI(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_riemann_sub)
			WriteInflowSubsonicRiemannBoundaryUpperI(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_isentrop_sub)
			WriteInflowSubsonicIsentropBoundaryUpperI(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_normal_super)
			WriteInflowSupersonicNormalBoundaryUpperI(
					pnt_config,
					pnt_mesh,
					pnt_U);
	}
	else if((strcmp(pnt_config->BC_Right,pnt_config->BCOutflow)==0)||
			(strcmp(pnt_config->BC_Right,pnt_config->BCOutflowSubsonic)==0))
	{
		if(pnt_config->flag_BC_option_outflow_normal_sub)
			WriteOutflowSubsonicNormalBoundaryUpperI(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_riemann_sub)
			WriteOutflowSubsonicRiemannBoundaryUpperI(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_rudy_sub)
			WriteOutflowSubsonicRudyBoundaryUpperI(
					pnt_config,
					pnt_mesh,
					pnt_U,
					pnt_U_lastStep);
	}
	else if(strcmp(pnt_config->BC_Right,pnt_config->BCFarfield)==0)
	{
		WriteFarfieldBoundaryUpperI(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Right,pnt_config->BCWallViscousIsothermal)==0)
	{
		WriteWallNoSlipIsothermalBoundaryUpperI(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Right,pnt_config->BCManufacturedSolution)==0)
	{
		WriteBCManufacturedSolutionUpperI(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else
	{
		printf("SHOCK: Error: Rank %d hat rechts weder Randbedingung noch Nachbar!\n",pnt_config->MPI_rank);
	}
}

void ControlBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	if(strcmp(pnt_config->BC_Top,pnt_config->BCWallViscous)==0)
	{
		WriteWallNoSlipBoundaryUpperJ(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Top,pnt_config->BCWallInviscid)==0)
	{
		WriteWallSlipBoundaryUpperJ(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if((strcmp(pnt_config->BC_Top,pnt_config->BCInflow)==0)||
			(strcmp(pnt_config->BC_Top,pnt_config->BCInflowSubsonic)==0)||
			(strcmp(pnt_config->BC_Top,pnt_config->BCInflowSupersonic)==0))
	{
		if(pnt_config->flag_BC_option_inflow_normal_sub)
			WriteInflowSubsonicNormalBoundaryUpperJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_riemann_sub)
			WriteInflowSubsonicRiemannBoundaryUpperJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_isentrop_sub)
			WriteInflowSubsonicIsentropBoundaryUpperJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_normal_super)
			WriteInflowSupersonicNormalBoundaryUpperJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
	}
	else if((strcmp(pnt_config->BC_Top,pnt_config->BCOutflow)==0)||
			(strcmp(pnt_config->BC_Top,pnt_config->BCOutflowSubsonic)==0))
	{
		if(pnt_config->flag_BC_option_outflow_normal_sub)
			WriteOutflowSubsonicNormalBoundaryUpperJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_riemann_sub)
			WriteOutflowSubsonicRiemannBoundaryUpperJ(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_rudy_sub)
			WriteOutflowSubsonicRudyBoundaryUpperJ(
					pnt_config,
					pnt_mesh,
					pnt_U,
					pnt_U_lastStep);
	}
	else if(strcmp(pnt_config->BC_Top,pnt_config->BCFarfield)==0)
	{
		WriteFarfieldBoundaryUpperJ(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Top,pnt_config->BCWallViscousIsothermal)==0)
	{
		WriteWallNoSlipIsothermalBoundaryUpperJ(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_Top,pnt_config->BCManufacturedSolution)==0)
	{
		WriteBCManufacturedSolutionUpperJ(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else
	{
		printf("SHOCK: Error: Rank %d hat oben weder Randbedingung noch Nachbar!\n",pnt_config->MPI_rank);
	}
}

void ControlBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	if(strcmp(pnt_config->BC_InFront,pnt_config->BCWallViscous)==0)
	{
		WriteWallNoSlipBoundaryUpperK(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_InFront,pnt_config->BCWallInviscid)==0)
	{
		WriteWallSlipBoundaryUpperK(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if((strcmp(pnt_config->BC_InFront,pnt_config->BCInflow)==0)||
			(strcmp(pnt_config->BC_InFront,pnt_config->BCInflowSubsonic)==0)||
			(strcmp(pnt_config->BC_InFront,pnt_config->BCInflowSupersonic)==0))
	{
		if(pnt_config->flag_BC_option_inflow_normal_sub)
			WriteInflowSubsonicNormalBoundaryUpperK(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_riemann_sub)
			WriteInflowSubsonicRiemannBoundaryUpperK(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_isentrop_sub)
			WriteInflowSubsonicIsentropBoundaryUpperK(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_inflow_normal_super)
			WriteInflowSupersonicNormalBoundaryUpperK(
					pnt_config,
					pnt_mesh,
					pnt_U);
	}
	else if((strcmp(pnt_config->BC_InFront,pnt_config->BCOutflow)==0)||
			(strcmp(pnt_config->BC_InFront,pnt_config->BCOutflowSubsonic)==0))
	{
		if(pnt_config->flag_BC_option_outflow_normal_sub)
			WriteOutflowSubsonicNormalBoundaryUpperK(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_riemann_sub)
			WriteOutflowSubsonicRiemannBoundaryUpperK(
					pnt_config,
					pnt_mesh,
					pnt_U);
		if(pnt_config->flag_BC_option_outflow_rudy_sub)
			WriteOutflowSubsonicRudyBoundaryUpperK(
					pnt_config,
					pnt_mesh,
					pnt_U,
					pnt_U_lastStep);
	}
	else if(strcmp(pnt_config->BC_InFront,pnt_config->BCFarfield)==0)
	{
		WriteFarfieldBoundaryUpperK(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_InFront,pnt_config->BCWallViscousIsothermal)==0)
	{
		WriteWallNoSlipIsothermalBoundaryUpperK(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else if(strcmp(pnt_config->BC_InFront,pnt_config->BCManufacturedSolution)==0)
	{
		WriteBCManufacturedSolutionUpperK(
				pnt_config,
				pnt_mesh,
				pnt_U);
	}
	else
	{
		printf("SHOCK: Error: Rank %d hat vorne weder Randbedingung noch Nachbar!\n",pnt_config->MPI_rank);
	}
}

void WriteFarfieldBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteFarfieldBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);

			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
}

void WriteFarfieldBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteFarfieldBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
}

void WriteFarfieldBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteFarfieldBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);

				int_symmetryIndex=int_symmetryIndex+1;


			}
		}
	}
}

void WriteFarfieldBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteFarfieldBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
}

void WriteFarfieldBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteFarfieldBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);

				int_symmetryIndex=int_symmetryIndex-1;

			}
		}
	}
}

void WriteFarfieldBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteFarfieldBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
}

void WriteWallNoSlipBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{

				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallNoSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
}

void WriteWallNoSlipBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};


	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{

				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteWallNoSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



				int_symmetryIndex=int_symmetryIndex+1;

			}
		}

	}
}

void WriteWallNoSlipBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};


	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{

				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallNoSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
}

void WriteWallNoSlipBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};


	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{

				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteWallNoSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



				int_symmetryIndex=int_symmetryIndex-1;

			}
		}
	}
}

void WriteWallNoSlipBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallNoSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
}

void WriteWallNoSlipBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallNoSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
}

void WriteWallSlipBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};


	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);




			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
}

void WriteWallSlipBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};


	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



				int_symmetryIndex=int_symmetryIndex+1;

			}
		}

	}
}

void WriteWallSlipBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};


	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
}

void WriteWallSlipBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};


	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex-1;
		}

	}
}

void WriteWallSlipBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};


	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



				int_symmetryIndex=int_symmetryIndex-1;

			}
		}

	}
}

void WriteWallSlipBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};


	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallSlipBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
}

void WriteInflowSubsonicRiemannBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;

	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
}

void WriteInflowSubsonicRiemannBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteInflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex-1;
			}
		}
	}
}

void WriteInflowSubsonicRiemannBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
}

void WriteInflowSubsonicRiemannBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
}

void WriteInflowSubsonicRiemannBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
}

void WriteInflowSubsonicRiemannBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteInflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex+1;
			}
		}
	}
}

void WriteInflowSupersonicNormalBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSupersonicNormalBoundary(corrector,ijk,pnt_config,pnt_mesh,pnt_U);


			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
}

void WriteInflowSupersonicNormalBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSupersonicNormalBoundary(corrector,ijk,pnt_config,pnt_mesh,pnt_U);



			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
}

void WriteInflowSupersonicNormalBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSupersonicNormalBoundary(corrector,ijk,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex+1;

			}
		}

	}
}

void WriteInflowSupersonicNormalBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSupersonicNormalBoundary(corrector,ijk,pnt_config,pnt_mesh,pnt_U);


			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
}

void WriteInflowSupersonicNormalBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSupersonicNormalBoundary(corrector,ijk,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
}

void WriteInflowSupersonicNormalBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSupersonicNormalBoundary(corrector,ijk,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex-1;
			}
		}
	}
}

void WriteInflowSubsonicIsentropBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicIsentropBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);

			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
}

void WriteInflowSubsonicIsentropBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicIsentropBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
}

void WriteInflowSubsonicIsentropBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteInflowSubsonicIsentropBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex+1;
			}
		}
	}
}

void WriteInflowSubsonicIsentropBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicIsentropBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
}

void WriteInflowSubsonicIsentropBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicIsentropBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
}

void WriteInflowSubsonicIsentropBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteInflowSubsonicIsentropBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex-1;
			}
		}
	}
}

void WriteInflowSubsonicNormalBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
}

void WriteInflowSubsonicNormalBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
}

void WriteInflowSubsonicNormalBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteInflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex+1;
			}
		}
	}
}

void WriteInflowSubsonicNormalBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
}

void WriteInflowSubsonicNormalBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteInflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
}

void WriteInflowSubsonicNormalBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteInflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex-1;
			}
		}
	}
}

void WriteOutflowSubsonicRudyBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicRudyBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U,pnt_U_lastStep);

			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
}

void WriteOutflowSubsonicRudyBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicRudyBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U,pnt_U_lastStep);
			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
}

void WriteOutflowSubsonicRudyBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{

				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteOutflowSubsonicRudyBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U,pnt_U_lastStep);

				int_symmetryIndex=int_symmetryIndex+1;
			}
		}

	}
}

void WriteOutflowSubsonicRudyBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicRudyBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U,pnt_U_lastStep);


			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
}

void WriteOutflowSubsonicRudyBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicRudyBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U,pnt_U_lastStep);


			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
}

void WriteOutflowSubsonicRudyBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteOutflowSubsonicRudyBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U,pnt_U_lastStep);


				int_symmetryIndex=int_symmetryIndex-1;
			}
		}
	}
}

void WriteOutflowSubsonicRiemannBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
}

void WriteOutflowSubsonicRiemannBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteOutflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex-1;
			}
		}
	}
}

void WriteOutflowSubsonicRiemannBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
}

void WriteOutflowSubsonicRiemannBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
}

void WriteOutflowSubsonicRiemannBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
}

void WriteOutflowSubsonicRiemannBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteOutflowSubsonicRiemannBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex+1;
			}
		}
	}
}

void WriteOutflowSubsonicNormalBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
}

void WriteOutflowSubsonicNormalBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteOutflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex-1;
			}
		}
	}
}

void WriteOutflowSubsonicNormalBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);



			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
}

void WriteOutflowSubsonicNormalBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
}

void WriteOutflowSubsonicNormalBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteOutflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);
			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
}

void WriteOutflowSubsonicNormalBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteOutflowSubsonicNormalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


				int_symmetryIndex=int_symmetryIndex+1;
			}
		}
	}
}

void WriteInflowSubsonicIsentropBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	pnt_U->theta1[ijk]=pnt_U->theta1[ijkSymmetry];
	pnt_U->theta2[ijk]=pnt_U->theta2[ijkSymmetry];
	pnt_U->theta3[ijk]=pnt_U->theta3[ijkSymmetry];

	pnt_U->u[ijk]=-(pnt_U->theta1[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]-pnt_U->theta1[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk])/(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);
	pnt_U->v[ijk]=-(-pnt_U->theta1[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]+pnt_U->theta1[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk])/(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);
	pnt_U->w[ijk]=-(-pnt_U->theta1[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_U->theta1[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->zeta_y[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->zeta_x[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk])/(pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);


	pnt_U->rho[ijk]=
			pow(
			pow(pnt_config->rho_inflow,(pnt_config->gammaNumber-1.0))-
			pow(pnt_config->rho_inflow,pnt_config->gammaNumber)/pnt_config->p_inflow*
			(pnt_config->gammaNumber-1.0)/pnt_config->gammaNumber*
			pnt_config->gammaNumber*pnt_config->machNumber*pnt_config->machNumber*
			(pnt_U->u[ijk]*pnt_U->u[ijk]/2.0)
			,1.0/(pnt_config->gammaNumber-1.0));
	pnt_U->p[ijk]=
			pnt_config->p_inflow*pow(pnt_U->rho[ijk]/pnt_config->rho_inflow,pnt_config->gammaNumber);

	pnt_U->e[ijk]=(0.5*((pnt_U->u[ijk]*pnt_U->u[ijk])+(pnt_U->v[ijk]*pnt_U->v[ijk])+(pnt_U->w[ijk]*pnt_U->w[ijk]))+
			pnt_U->p[ijk]/pnt_U->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
	pnt_U->c[ijk]=sqrt(pnt_config->Upsilon*
			pnt_config->gammaNumber*pnt_U->p[ijk]/pnt_U->rho[ijk]);




//				Viskose Randbedingungen
	pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
	pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];
}

void WriteInflowSubsonicRiemannBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	FLT cinf,d__1,d__2,rrr,cii;
	cinf=1.0/pnt_config->machNumber;
	cii=pnt_U->c[ijkSymmetry];


	pnt_U->u[ijk] = 0.5*(pnt_U->u[ijkSymmetry]+pnt_config->u_inflow)+(cinf-cii)/(pnt_config->gammaNumber-1.);
	pnt_U->v[ijk] = pnt_config->v_inflow;
	pnt_U->w[ijk] = pnt_config->w_inflow;

	d__1 = (0.5*(cii+cinf)+0.25*(pnt_config->gammaNumber-1.)*(pnt_config->u_inflow-pnt_U->u[ijkSymmetry]))/cinf;
	d__2 = 2./(pnt_config->gammaNumber-1.);
	rrr = pow(d__1, d__2);
	pnt_U->rho[ijk]= pnt_config->rho_inflow*rrr;

	d__2 = 2.*pnt_config->gammaNumber/(pnt_config->gammaNumber-1.);
	rrr = pow(d__1, d__2);
	pnt_U->p[ijk]= pnt_config->p_inflow*rrr;

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


	pnt_U->c[ijk]=sqrt(pnt_config->Upsilon*
			pnt_config->gammaNumber*pnt_U->p[ijk]/pnt_U->rho[ijk]);


	pnt_U->e[ijk]=(0.5*((pnt_U->u[ijk]*pnt_U->u[ijk])+(pnt_U->v[ijk]*pnt_U->v[ijk])+(pnt_U->w[ijk]*pnt_U->w[ijk]))+
			pnt_U->p[ijk]/pnt_U->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);


//				Viskose Randbedingungen
	pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
	pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];
}

void WriteInflowSubsonicNormalBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	pnt_U->u[ijk]=pnt_config->u_inflow;
	pnt_U->v[ijk]=pnt_config->v_inflow;
	pnt_U->w[ijk]=pnt_config->w_inflow;

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

	pnt_U->p[ijk]=pnt_U->p[ijkSymmetry];
//				T_unendlich=p_unendlich/rho_unendlich
//				1.RB: T_unendlich=1.0  ---> Ma_unendlich=const. ueber Einströmrand
//				T_unendlich=1.0 ---> p_unendlich=rho_unendlich
	pnt_U->rho[ijk]=pnt_U->p[ijk];
	//Veraendert am 24.10.2014:
	pnt_U->p[ijk]=pnt_config->p_inflow;
	pnt_U->rho[ijk]=pnt_config->rho_inflow;

	pnt_U->e[ijk]=(0.5*((pnt_U->u[ijk]*pnt_U->u[ijk])+(pnt_U->v[ijk]*pnt_U->v[ijk])+(pnt_U->w[ijk]*pnt_U->w[ijk]))+
			pnt_U->p[ijk]/pnt_U->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
	pnt_U->c[ijk]=sqrt(pnt_config->Upsilon*
			pnt_config->gammaNumber*pnt_U->p[ijk]/pnt_U->rho[ijk]);




//				Viskose Randbedingungen
	pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
	pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];
	//Veraendert am 24.10.2014:
}

void WriteInflowSupersonicNormalBoundary(
		FLT corrector[3],
		int ijk,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	pnt_U->p[ijk]=pnt_config->p_inflow;
	pnt_U->rho[ijk]=pnt_config->rho_inflow;

	pnt_U->u[ijk]=pnt_config->u_inflow;
	pnt_U->v[ijk]=pnt_config->v_inflow;
	pnt_U->w[ijk]=pnt_config->w_inflow;

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
			pnt_U->p[ijk]/pnt_U->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
	pnt_U->c[ijk]=sqrt(pnt_config->Upsilon*
			pnt_config->gammaNumber*pnt_U->p[ijk]/pnt_U->rho[ijk]);


//				Viskose Randbedingungen
	pnt_U->T[ijk]=pnt_config->p_inflow/pnt_config->rho_inflow;
	pnt_U->mue[ijk]=((1.0+pnt_config->SutherlandConstant)*pow(pnt_U->p[ijk]/pnt_U->rho[ijk],1.5)/
			(pnt_U->p[ijk]/pnt_U->rho[ijk]+pnt_config->SutherlandConstant));
}

void WriteOutflowSubsonicNormalBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	pnt_U->theta1[ijk]=pnt_U->theta1[ijkSymmetry];
	pnt_U->theta2[ijk]=pnt_U->theta2[ijkSymmetry];
	pnt_U->theta3[ijk]=pnt_U->theta3[ijkSymmetry];

	pnt_U->u[ijk]=-(pnt_U->theta1[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]-pnt_U->theta1[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk])/(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);
	pnt_U->v[ijk]=-(-pnt_U->theta1[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]+pnt_U->theta1[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk])/(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);
	pnt_U->w[ijk]=-(-pnt_U->theta1[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_U->theta1[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->zeta_y[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->zeta_x[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk])/(pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);


	pnt_U->p[ijk]=pnt_config->p_out;
	pnt_U->rho[ijk]=pnt_U->rho[ijkSymmetry];
	pnt_U->e[ijk]=pnt_U->e[ijkSymmetry];
	pnt_U->c[ijk]=pnt_U->c[ijkSymmetry];


//				Viskose Randbedingungen
	pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
	pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];
}

void WriteOutflowSubsonicRiemannBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	FLT cinf,d__1,d__2,rrr,cii;
	cinf=1.0/pnt_config->machNumber;
	cii=pnt_U->c[ijkSymmetry];

	pnt_U->u[ijk] = 0.5*(pnt_U->u[ijkSymmetry]+pnt_config->u_inflow)-(cinf-cii)/(pnt_config->gammaNumber-1.);
	pnt_U->v[ijk] = pnt_U->v[ijkSymmetry];
	pnt_U->w[ijk] = pnt_U->w[ijkSymmetry];

	d__1 = (0.5*(cii+cinf)-0.25*(pnt_config->gammaNumber-1.)*(pnt_config->u_inflow-pnt_U->u[ijkSymmetry]))/cii;
	d__2 = 2./(pnt_config->gammaNumber-1.);
	rrr = pow(d__1, d__2);
	pnt_U->rho[ijk]= pnt_U->rho[ijkSymmetry]*rrr;

	d__2 = 2.*pnt_config->gammaNumber/(pnt_config->gammaNumber-1.);
	rrr = pow(d__1, d__2);
	pnt_U->p[ijk]= pnt_U->p[ijkSymmetry]*rrr;

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


	pnt_U->c[ijk]=sqrt(pnt_config->Upsilon*
			pnt_config->gammaNumber*pnt_U->p[ijk]/pnt_U->rho[ijk]);


	pnt_U->e[ijk]=(0.5*((pnt_U->u[ijk]*pnt_U->u[ijk])+(pnt_U->v[ijk]*pnt_U->v[ijk])+(pnt_U->w[ijk]*pnt_U->w[ijk]))+
			pnt_U->p[ijk]/pnt_U->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);


//				Viskose Randbedingungen
	pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
	pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];
}

void WriteOutflowSubsonicRudyBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep)

{
	pnt_U->theta1[ijk]=pnt_U->theta1[ijkSymmetry];
	pnt_U->theta2[ijk]=pnt_U->theta2[ijkSymmetry];
	pnt_U->theta3[ijk]=pnt_U->theta3[ijkSymmetry];

	pnt_U->u[ijk]=-(pnt_U->theta1[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]-pnt_U->theta1[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk])/(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);
	pnt_U->v[ijk]=-(-pnt_U->theta1[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]+pnt_U->theta1[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk])/(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);
	pnt_U->w[ijk]=-(-pnt_U->theta1[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_U->theta1[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->zeta_y[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->zeta_x[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk])/(pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);


	pnt_U->p[ijk]=
			(pnt_U_lastStep->p[ijk]+pnt_config->AlphaNonRef*pnt_config->numericalTau*pnt_config->p_out+
			pnt_U_lastStep->rho[ijk]*pnt_U_lastStep->c[ijkSymmetry]*(pnt_U->u[ijk]-pnt_U_lastStep->u[ijk])/pnt_config->Upsilon)
			/(1.0+pnt_config->AlphaNonRef*pnt_config->numericalTau);

	pnt_U->rho[ijk]=pnt_U->p[ijk]/pnt_U->p[ijkSymmetry]*pnt_U->rho[ijkSymmetry];
	pnt_U->e[ijk]=pnt_U->e[ijkSymmetry];
	pnt_U->c[ijk]=pnt_U->c[ijkSymmetry];

//				Viskose Randbedingungen
	pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
	pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];
}

void WriteWallNoSlipBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	pnt_U->theta1[ijk]=-pnt_U->theta1[ijkSymmetry];
	pnt_U->theta2[ijk]=-pnt_U->theta2[ijkSymmetry];
	pnt_U->theta3[ijk]=-pnt_U->theta3[ijkSymmetry];

	pnt_U->u[ijk]=-pnt_U->u[ijkSymmetry];
	pnt_U->v[ijk]=-pnt_U->v[ijkSymmetry];
	pnt_U->w[ijk]=-pnt_U->w[ijkSymmetry];

	pnt_mesh->BC_Corrector_xiMomentum[ijk]=-1.0;
	pnt_mesh->BC_Corrector_etaMomentum[ijk]=-1.0;
	pnt_mesh->BC_Corrector_zetaMomentum[ijk]=-1.0;

	pnt_U->rho[ijk]=pnt_U->rho[ijkSymmetry];
	pnt_U->p[ijk]=pnt_U->p[ijkSymmetry];
	pnt_U->e[ijk]=pnt_U->e[ijkSymmetry];
	pnt_U->c[ijk]=pnt_U->c[ijkSymmetry];


//				Viskose Randbedingungen
	pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
	pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];
}

void WriteWallSlipBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	pnt_U->theta1[ijk]=corrector[0]*pnt_U->theta1[ijkSymmetry];
	pnt_U->theta2[ijk]=corrector[1]*pnt_U->theta2[ijkSymmetry];
	pnt_U->theta3[ijk]=corrector[2]*pnt_U->theta3[ijkSymmetry];

	pnt_U->u[ijk]=-(pnt_U->theta1[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]-pnt_U->theta1[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk])/(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);
	pnt_U->v[ijk]=-(-pnt_U->theta1[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]+pnt_U->theta1[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk])/(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);
	pnt_U->w[ijk]=-(-pnt_U->theta1[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_U->theta1[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->zeta_y[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->zeta_x[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk])/(pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);

	pnt_mesh->BC_Corrector_xiMomentum[ijk]=corrector[0];
	pnt_mesh->BC_Corrector_etaMomentum[ijk]=corrector[1];
	pnt_mesh->BC_Corrector_zetaMomentum[ijk]=corrector[2];


//	printf("u:%g v:%g || u':%g v':%g\n",pnt_U->u[ijk],pnt_U->v[ijk],pnt_U->u[ijkSymmetry],pnt_U->v[ijkSymmetry]);

	pnt_U->rho[ijk]=pnt_U->rho[ijkSymmetry];
	pnt_U->p[ijk]=pnt_U->p[ijkSymmetry];
	pnt_U->e[ijk]=pnt_U->e[ijkSymmetry];
	pnt_U->c[ijk]=pnt_U->c[ijkSymmetry];

//				Viskose Randbedingungen
	pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
	pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];
}

void WriteMovingWallSlipBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{   
	//pnt_U->u[ijk]=-pnt_U->u[ijkSymmetry]+2.0*(pnt_config->IBC_MovingActualPosition-pnt_config->IBC_MovingLastPosition)/pnt_config->numericalTau;
	pnt_U->u[ijk]=(pnt_config->IBC_MovingActualPosition-pnt_config->IBC_MovingLastPosition)/pnt_config->numericalTau;
	pnt_U->v[ijk]=-pnt_U->v[ijkSymmetry];
	pnt_U->w[ijk]=-pnt_U->w[ijkSymmetry];

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
			
	pnt_mesh->BC_Corrector_xiMomentum[ijk]=1.0;
	pnt_mesh->BC_Corrector_etaMomentum[ijk]=1.0;
	pnt_mesh->BC_Corrector_zetaMomentum[ijk]=1.0;

    pnt_U->p[ijk]=pnt_U->p[ijkSymmetry];
	pnt_U->rho[ijk]=pnt_U->rho[ijkSymmetry];
	pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];
	
	pnt_U->c[ijk]=pnt_U->c[ijkSymmetry];
    pnt_U->e[ijk]=(0.5*((pnt_U->u[ijk]*pnt_U->u[ijk])+(pnt_U->v[ijk]*pnt_U->v[ijk])+(pnt_U->w[ijk]*pnt_U->w[ijk]))+
				pnt_U->p[ijk]/pnt_U->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
				
//				Viskose Randbedingungen
	pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
		
}
void WriteMovingWallNoSlipIsothermalBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	FLT prandtl_local;
	//pnt_U->u[ijk]=-pnt_U->u[ijkSymmetry]+2.0*(pnt_config->IBC_MovingActualPosition-pnt_config->IBC_MovingLastPosition)/pnt_config->numericalTau;
	pnt_U->u[ijk]=(pnt_config->IBC_MovingActualPosition-pnt_config->IBC_MovingLastPosition)/pnt_config->numericalTau;
	pnt_U->v[ijk]=-pnt_U->v[ijkSymmetry];
	pnt_U->w[ijk]=-pnt_U->w[ijkSymmetry];

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

	pnt_mesh->BC_Corrector_xiMomentum[ijk]=1.0;
	pnt_mesh->BC_Corrector_etaMomentum[ijk]=1.0;
	pnt_mesh->BC_Corrector_zetaMomentum[ijk]=1.0;

	pnt_U->p[ijk]=pnt_U->p[ijkSymmetry];
	pnt_U->rho[ijk]=pnt_U->rho[ijkSymmetry];
	pnt_U->e[ijk]=(0.5*((pnt_U->u[ijk]*pnt_U->u[ijk])+(pnt_U->v[ijk]*pnt_U->v[ijk])+(pnt_U->w[ijk]*pnt_U->w[ijk]))+
				pnt_U->p[ijk]/pnt_U->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
	pnt_U->c[ijk]=sqrt(pnt_config->Upsilon*pnt_config->gammaNumber*pnt_U->p[ijk]/pnt_U->rho[ijk]);

					
//				Viskose Randbedingungen
	pnt_U->mue[ijk]=((1.0+pnt_config->SutherlandConstant)*pow(pnt_U->p[ijk]/pnt_U->rho[ijk],1.5)/
	                (pnt_U->p[ijk]/pnt_U->rho[ijk]+pnt_config->SutherlandConstant));

//  adiabat moving wand
//	pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];

//  isothermal moving wand
	pnt_U->T[ijk]=1.0;
	
// um die Wärmeflüsse am Wand zu kontrollieren.			
    prandtl_local=0.72;
    pnt_config->Gamma[ijk]=1.0/((pnt_config->gammaNumber-1.0)*pow(pnt_config->machNumber,2.0)*pnt_config->reynoldsNumber*prandtl_local);			
		
}

void WriteFarfieldBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	pnt_U->u[ijk]=pnt_U->u[ijkSymmetry];
	pnt_U->v[ijk]=pnt_U->v[ijkSymmetry];
	pnt_U->w[ijk]=pnt_U->w[ijkSymmetry];

//	pnt_U->theta1[ijk]=
//			pnt_U->u[ijk]*pnt_mesh->xi_x[ijk]+
//			pnt_U->v[ijk]*pnt_mesh->xi_y[ijk]+
//			pnt_U->w[ijk]*pnt_mesh->xi_z[ijk];
//	pnt_U->theta2[ijk]=
//			pnt_U->u[ijk]*pnt_mesh->eta_x[ijk]+
//			pnt_U->v[ijk]*pnt_mesh->eta_y[ijk]+
//			pnt_U->w[ijk]*pnt_mesh->eta_z[ijk];
//	pnt_U->theta3[ijk]=
//			pnt_U->u[ijk]*pnt_mesh->zeta_x[ijk]+
//			pnt_U->v[ijk]*pnt_mesh->zeta_y[ijk]+
//			pnt_U->w[ijk]*pnt_mesh->zeta_z[ijk];

	pnt_U->theta1[ijk]=pnt_U->theta1[ijkSymmetry];
	pnt_U->theta2[ijk]=pnt_U->theta2[ijkSymmetry];
	pnt_U->theta3[ijk]=pnt_U->theta3[ijkSymmetry];
//
//	pnt_U->u[ijk]=-(pnt_U->theta1[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]-pnt_U->theta1[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk])/(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);
//	pnt_U->v[ijk]=-(-pnt_U->theta1[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]+pnt_U->theta1[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk])/(-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);
//	pnt_U->w[ijk]=-(-pnt_U->theta1[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]+pnt_U->theta1[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]+pnt_U->theta2[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->zeta_y[ijk]-pnt_U->theta2[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->zeta_x[ijk]-pnt_U->theta3[ijk]*pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]+pnt_U->theta3[ijk]*pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk])/(pnt_mesh->xi_x[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_z[ijk]-pnt_mesh->xi_x[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_y[ijk]-pnt_mesh->xi_y[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_z[ijk]+pnt_mesh->xi_y[ijk]*pnt_mesh->eta_z[ijk]*pnt_mesh->zeta_x[ijk]+pnt_mesh->xi_z[ijk]*pnt_mesh->eta_x[ijk]*pnt_mesh->zeta_y[ijk]-pnt_mesh->xi_z[ijk]*pnt_mesh->eta_y[ijk]*pnt_mesh->zeta_x[ijk]);


	pnt_U->rho[ijk]=pnt_U->rho[ijkSymmetry];
	pnt_U->p[ijk]=pnt_U->p[ijkSymmetry];
	pnt_U->e[ijk]=pnt_U->e[ijkSymmetry];
	pnt_U->c[ijk]=pnt_U->c[ijkSymmetry];

//				Viskose Randbedingungen
	pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
	pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];
}

void WriteMovingBCUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT x_shock,u_shock;
	u_shock=10.0;

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				x_shock = 1./6. + pnt_mesh->y[ijk]/tan(60.*MY_PI/180.) + pnt_config->time_dim * u_shock / sin(60.*MY_PI/180.);
				if (pnt_mesh->x[ijk] > x_shock)
				{

				pnt_U->rho[ijk]=pnt_U->rho[ijkSymmetry];
				pnt_U->p[ijk]=pnt_U->p[ijkSymmetry];
				pnt_U->e[ijk]=pnt_U->e[ijkSymmetry];
				pnt_U->c[ijk]=pnt_U->c[ijkSymmetry];

				pnt_U->theta1[ijk]=pnt_U->theta1[ijkSymmetry];
				pnt_U->theta2[ijk]=pnt_U->theta2[ijkSymmetry];
				pnt_U->theta3[ijk]=pnt_U->theta3[ijkSymmetry];

				pnt_U->u[ijk]=pnt_U->u[ijkSymmetry];
				pnt_U->v[ijk]=pnt_U->v[ijkSymmetry];
				pnt_U->w[ijk]=pnt_U->w[ijkSymmetry];
//				Viskose Randbedingungen

				pnt_U->T[ijk]=pnt_U->T[ijkSymmetry];
				pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
				}
				else
				{
				pnt_U->p[ijk]=pnt_config->p_inflow;
				pnt_U->rho[ijk]=pnt_config->rho_inflow;
				pnt_U->u[ijk]=pnt_config->u_inflow;
				pnt_U->v[ijk]=pnt_config->v_inflow;
				pnt_U->w[ijk]=pnt_config->w_inflow;

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
						pnt_U->p[ijk]/pnt_U->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
				pnt_U->c[ijk]=sqrt(pnt_config->Upsilon*
						pnt_config->gammaNumber*pnt_U->p[ijk]/pnt_U->rho[ijk]);




//				Viskose Randbedingungen

				pnt_U->T[ijk]=1.0;
				pnt_U->mue[ijk]=((1.0+pnt_config->SutherlandConstant)*pow(pnt_U->p[ijk]/pnt_U->rho[ijk],1.5)/
						(pnt_U->p[ijk]/pnt_U->rho[ijk]+pnt_config->SutherlandConstant));
				}
			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
}

void WriteWallNoSlipIsothermalBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{   
    FLT prandtl_local;
	pnt_U->theta1[ijk]=-pnt_U->theta1[ijkSymmetry];
	pnt_U->theta2[ijk]=-pnt_U->theta2[ijkSymmetry];
	pnt_U->theta3[ijk]=-pnt_U->theta3[ijkSymmetry];

	pnt_U->u[ijk]=-pnt_U->u[ijkSymmetry];
	pnt_U->v[ijk]=-pnt_U->v[ijkSymmetry];
	pnt_U->w[ijk]=-pnt_U->w[ijkSymmetry];

	pnt_mesh->BC_Corrector_xiMomentum[ijk]=-1.0;
	pnt_mesh->BC_Corrector_etaMomentum[ijk]=-1.0;
	pnt_mesh->BC_Corrector_zetaMomentum[ijk]=-1.0;

    pnt_U->p[ijk]=pnt_U->p[ijkSymmetry];
	pnt_U->rho[ijk]=pnt_U->rho[ijkSymmetry];
	
	pnt_U->e[ijk]=pnt_U->e[ijkSymmetry];
	pnt_U->c[ijk]=pnt_U->c[ijkSymmetry];

//				Viskose Randbedingungen
	pnt_U->mue[ijk]=pnt_U->mue[ijkSymmetry];
	pnt_U->T[ijk]=pnt_config->T_wall;
	
	// um die Wärmeflüsse an Wand zu kontrollieren.			
    prandtl_local=0.72;
    pnt_config->Gamma[ijk]=1.0/((pnt_config->gammaNumber-1.0)*pow(pnt_config->machNumber,2.0)*pnt_config->reynoldsNumber*prandtl_local);			
}


void WriteWallNoSlipIsothermalBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{

				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallNoSlipIsothermalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);

			}
			int_symmetryIndex=int_symmetryIndex+1;
		}

	}
}

void WriteWallNoSlipIsothermalBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{

				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteWallNoSlipIsothermalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);

				int_symmetryIndex=int_symmetryIndex+1;

			}
		}

	}
}

void WriteWallNoSlipIsothermalBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,-1.,1.};

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{

				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallNoSlipIsothermalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
			int_symmetryIndex=int_symmetryIndex-1;
		}
	}
}

void WriteWallNoSlipIsothermalBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={1.,1.,-1.};

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{

				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

				WriteWallNoSlipIsothermalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);

				int_symmetryIndex=int_symmetryIndex-1;

			}
		}
	}
}

void WriteWallNoSlipIsothermalBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallNoSlipIsothermalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);

			}
		}
		int_symmetryIndex=int_symmetryIndex-1;
	}
}

void WriteWallNoSlipIsothermalBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U)
{
	int i,j,k,ijk,ijkSymmetry;
	int int_symmetryIndex;
	FLT corrector[3]={-1.,1.,1.};

	int_symmetryIndex=pnt_config->int_iStartReal;
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteWallNoSlipIsothermalBoundary(corrector,ijk,ijkSymmetry,pnt_config,pnt_mesh,pnt_U);


			}
		}
		int_symmetryIndex=int_symmetryIndex+1;
	}
}

