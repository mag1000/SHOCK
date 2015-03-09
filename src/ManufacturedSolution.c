// Manufactured Solutions
// These functions are needed for the implementation of the Manufactured Solution Approach
// in which an analytical method is used to verify the code. It can be used to ensure that the implementation is bug-free
// This is done by comparing the numerical approximated solution to the exact solution of the governing equations.
// The exact solution must be derived a priori and then, this is added as a source term. The resulting exact analytical solution
// is known and the difference to the numerical solution is the numerical error.
// This error must decrease with the order of the code when the code is bug-free

//#define _ISOC99_SOURCE  1

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "float.h"
#include "string.h"
#include "cgnslib.h"

#include "SHOCK.h"
#include "Functions.h"
#include "ManufacturedSolution.h"

int i,j,k,ijk;

void AddManufacturedSolutionSource(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_Flux * pnt_Q)
{
	long double L=(long double)1.0;
	long double rho_ref=(long double)1.0;
	long double u_ref=(long double)pnt_config->u0_dim;
	long double p_ref=(long double)rho_ref*u_ref*u_ref*pnt_config->Upsilon;


	long double rho_0,rho_x,rho_y,a_rho_x,a_rho_y;
	long double u_0,u_x,u_y,a_u_x,a_u_y;
	long double v_0,v_x,v_y,a_v_x,a_v_y;
	long double p_0,p_x,p_y,a_p_x,a_p_y;

	if (pnt_config->ManufacturedSolution_case==1)
	{
	//supersonic
	rho_0=1.0/rho_ref;	rho_x=0.15/rho_ref;		rho_y=-0.1/rho_ref;		a_rho_x=1.0;	a_rho_y=0.5;
	u_0=800.0/u_ref;	u_x=50.0/u_ref;			u_y=-30.0/u_ref;		a_u_x=1.5;		a_u_y=0.6;
	v_0=800.0/u_ref;	v_x=-75.0/u_ref;		v_y=40.0/u_ref;			a_v_x=0.5;		a_v_y=2./3.;
	p_0=100000.0/p_ref;	p_x=0.2*100000.0/p_ref;	p_y=0.5*100000.0/p_ref;	a_p_x=2.0;		a_p_y=1.0;
	}
	else // (pnt_config->ManufacturedSolution_case==0)
	{
	//subsonic
	rho_0=1.0/rho_ref;	rho_x=0.15/rho_ref;		rho_y=-0.1/rho_ref;		a_rho_x=1.0;	a_rho_y=0.5;
	u_0=70.0/u_ref;		u_x=5.0/u_ref;			u_y=-7.0/u_ref;			a_u_x=1.5;		a_u_y=0.6;
	v_0=90.0/u_ref;		v_x=-15.0/u_ref;		v_y=8.5/u_ref;			a_v_x=0.5;		a_v_y=2./3.;
	p_0=100000.0/p_ref;	p_x=0.2*100000.0/p_ref;	p_y=0.5*100000.0/p_ref;	a_p_x=2.0;		a_p_y=1.0;
	}

	long double mass,x_momentum,y_momentum,energy;
	long double CoordX,CoordY,Upsilon,gammaNumber;

	Upsilon=(long double)pnt_config->Upsilon;
	gammaNumber=(long double)pnt_config->gammaNumber;

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				CoordX=(long double)pnt_mesh->x[ijk];
				CoordY=(long double)pnt_mesh->y[ijk];

				mass = (M_PIl*a_u_x*u_x*cosl((CoordX*M_PIl*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L)))/L+(M_PIl*a_v_y*v_y*cosl((CoordY*M_PIl*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L)))/L+(M_PIl*a_rho_x*rho_x*cosl((CoordX*M_PIl*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L)))/L-(M_PIl*a_rho_y*rho_y*sinl((CoordY*M_PIl*a_rho_y)/L)*(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L)))/L;
				x_momentum = (M_PIl*a_rho_x*rho_x*cosl((CoordX*M_PIl*a_rho_x)/L)*powl(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L),2.0))/L-(M_PIl*Upsilon*a_p_x*p_x*sinl((CoordX*M_PIl*a_p_x)/L))/L+(M_PIl*a_u_x*u_x*cosl((CoordX*M_PIl*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))*(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L))*2.0)/L+(M_PIl*a_v_y*v_y*cosl((CoordY*M_PIl*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))*(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L)))/L-(M_PIl*a_u_y*u_y*sinl((CoordY*M_PIl*a_u_y)/L)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))*(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L)))/L-(M_PIl*a_rho_y*rho_y*sinl((CoordY*M_PIl*a_rho_y)/L)*(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L))*(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L)))/L;
				y_momentum = -(M_PIl*a_rho_y*rho_y*sinl((CoordY*M_PIl*a_rho_y)/L)*powl(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L),2.0))/L+(M_PIl*Upsilon*a_p_y*p_y*cosl((CoordY*M_PIl*a_p_y)/L))/L+(M_PIl*a_u_x*u_x*cosl((CoordX*M_PIl*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))*(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L)))/L+(M_PIl*a_v_y*v_y*cosl((CoordY*M_PIl*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))*(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L))*2.0)/L+(M_PIl*a_rho_x*rho_x*cosl((CoordX*M_PIl*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L))*(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L)))/L-(M_PIl*a_v_x*v_x*sinl((CoordX*M_PIl*a_v_x)/L)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))*(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L)))/L;
				energy = -(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L))*((rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))*(-(M_PIl*a_u_x*u_x*cosl((CoordX*M_PIl*a_u_x)/L)*(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L)))/L+(M_PIl*a_v_x*v_x*sinl((CoordX*M_PIl*a_v_x)/L)*(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L)))/L+(M_PIl*Upsilon*a_p_x*p_x*sinl((CoordX*M_PIl*a_p_x)/L))/(L*(gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L)))+(M_PIl*Upsilon*a_rho_x*rho_x*cosl((CoordX*M_PIl*a_rho_x)/L)*(p_0+p_x*cosl((CoordX*M_PIl*a_p_x)/L)+p_y*sinl((CoordY*M_PIl*a_p_y)/L))*1.0/powl(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L),2.0))/(L*(gammaNumber-1.0)))-(M_PIl*a_rho_x*rho_x*cosl((CoordX*M_PIl*a_rho_x)/L)*(powl(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L),2.0)*(1.0/2.0)+powl(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L),2.0)*(1.0/2.0)+(Upsilon*(p_0+p_x*cosl((CoordX*M_PIl*a_p_x)/L)+p_y*sinl((CoordY*M_PIl*a_p_y)/L)))/((gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L)))))/L+(M_PIl*Upsilon*a_p_x*p_x*sinl((CoordX*M_PIl*a_p_x)/L))/L)+(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L))*((rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))*((M_PIl*a_v_y*v_y*cosl((CoordY*M_PIl*a_v_y)/L)*(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L)))/L-(M_PIl*a_u_y*u_y*sinl((CoordY*M_PIl*a_u_y)/L)*(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L)))/L+(M_PIl*Upsilon*a_p_y*p_y*cosl((CoordY*M_PIl*a_p_y)/L))/(L*(gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L)))+(M_PIl*Upsilon*a_rho_y*rho_y*sinl((CoordY*M_PIl*a_rho_y)/L)*(p_0+p_x*cosl((CoordX*M_PIl*a_p_x)/L)+p_y*sinl((CoordY*M_PIl*a_p_y)/L))*1.0/powl(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L),2.0))/(L*(gammaNumber-1.0)))+(M_PIl*Upsilon*a_p_y*p_y*cosl((CoordY*M_PIl*a_p_y)/L))/L-(M_PIl*a_rho_y*rho_y*sinl((CoordY*M_PIl*a_rho_y)/L)*(powl(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L),2.0)*(1.0/2.0)+powl(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L),2.0)*(1.0/2.0)+(Upsilon*(p_0+p_x*cosl((CoordX*M_PIl*a_p_x)/L)+p_y*sinl((CoordY*M_PIl*a_p_y)/L)))/((gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L)))))/L)+(M_PIl*a_u_x*u_x*cosl((CoordX*M_PIl*a_u_x)/L)*((rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))*(powl(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L),2.0)*(1.0/2.0)+powl(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L),2.0)*(1.0/2.0)+(Upsilon*(p_0+p_x*cosl((CoordX*M_PIl*a_p_x)/L)+p_y*sinl((CoordY*M_PIl*a_p_y)/L)))/((gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))))+Upsilon*(p_0+p_x*cosl((CoordX*M_PIl*a_p_x)/L)+p_y*sinl((CoordY*M_PIl*a_p_y)/L))))/L+(M_PIl*a_v_y*v_y*cosl((CoordY*M_PIl*a_v_y)/L)*((rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))*(powl(u_0+u_y*cosl((CoordY*M_PIl*a_u_y)/L)+u_x*sinl((CoordX*M_PIl*a_u_x)/L),2.0)*(1.0/2.0)+powl(v_0+v_x*cosl((CoordX*M_PIl*a_v_x)/L)+v_y*sinl((CoordY*M_PIl*a_v_y)/L),2.0)*(1.0/2.0)+(Upsilon*(p_0+p_x*cosl((CoordX*M_PIl*a_p_x)/L)+p_y*sinl((CoordY*M_PIl*a_p_y)/L)))/((gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*M_PIl*a_rho_y)/L)+rho_x*sinl((CoordX*M_PIl*a_rho_x)/L))))+Upsilon*(p_0+p_x*cosl((CoordX*M_PIl*a_p_x)/L)+p_y*sinl((CoordY*M_PIl*a_p_y)/L))))/L;



				pnt_Q->Mass[ijk]+=pnt_mesh->jacobian[ijk]*mass;

				pnt_Q->xiMomentum[ijk]+=pnt_mesh->jacobian[ijk]*x_momentum;

				pnt_Q->etaMomentum[ijk]+=pnt_mesh->jacobian[ijk]*y_momentum;

				pnt_Q->Energy[ijk]+=pnt_mesh->jacobian[ijk]*energy;
			}
		}
	}
}

void InitializeManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
		{
			for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
			{
				for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					pnt_mesh->BC_Corrector[ijk]=1.0;
					WriteManufacturedSolution(
							pnt_config,
							pnt_mesh,
							pnt_U_lastStep,
							ijk);
				}
			}
		}
}

void WriteBCManufacturedSolutionLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteManufacturedSolution(
						pnt_config,
						pnt_mesh,
						pnt_U_lastStep,
						ijk);
			}
		}
	}
}

void WriteBCManufacturedSolutionLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteManufacturedSolution(
						pnt_config,
						pnt_mesh,
						pnt_U_lastStep,
						ijk);
			}
		}
	}
}

void WriteBCManufacturedSolutionLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteManufacturedSolution(
						pnt_config,
						pnt_mesh,
						pnt_U_lastStep,
						ijk);
			}
		}
	}
}

void WriteBCManufacturedSolutionUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteManufacturedSolution(
						pnt_config,
						pnt_mesh,
						pnt_U_lastStep,
						ijk);
			}
		}
	}
}

void WriteBCManufacturedSolutionUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteManufacturedSolution(
						pnt_config,
						pnt_mesh,
						pnt_U_lastStep,
						ijk);
			}
		}
	}
}

void WriteBCManufacturedSolutionUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				WriteManufacturedSolution(
						pnt_config,
						pnt_mesh,
						pnt_U_lastStep,
						ijk);
			}
		}
	}
}

void ConfigureManufacturedSolution(
		struct strct_configuration * pnt_config)
{
	sprintf(pnt_config->BCManufacturedSolution,"BCManufacturedSolution");

	if (pnt_config->ManufacturedSolution_case==0)
	{
		if (pnt_config->InterfaceNeighbourLeft==NO_NEIGHBOUR)
			{strcpy( pnt_config->BC_Left,pnt_config->BCManufacturedSolution );}

		if (pnt_config->InterfaceNeighbourBottom==NO_NEIGHBOUR)
			{strcpy( pnt_config->BC_Bottom,pnt_config->BCManufacturedSolution );}

		if (pnt_config->InterfaceNeighbourTop==NO_NEIGHBOUR)
			{strcpy( pnt_config->BC_Top,pnt_config->BCManufacturedSolution );}

		if (pnt_config->InterfaceNeighbourRight==NO_NEIGHBOUR)
			{strcpy( pnt_config->BC_Right,pnt_config->BCManufacturedSolution );}
	}

	if (pnt_config->ManufacturedSolution_case==1)
	{
		if (pnt_config->InterfaceNeighbourLeft==NO_NEIGHBOUR)
			{strcpy( pnt_config->BC_Left,pnt_config->BCManufacturedSolution );}

		if (pnt_config->InterfaceNeighbourBottom==NO_NEIGHBOUR)
			{strcpy( pnt_config->BC_Bottom,pnt_config->BCManufacturedSolution );}

		if (pnt_config->InterfaceNeighbourTop==NO_NEIGHBOUR)
			{strcpy( pnt_config->BC_Top,pnt_config->BCFarfield );}

		if (pnt_config->InterfaceNeighbourRight==NO_NEIGHBOUR)
			{strcpy( pnt_config->BC_Right,pnt_config->BCFarfield );}
	}
}


void WriteManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		int ijk)
{
	long double L=(long double)1.0;
	long double rho_ref=(long double)1.0;
	long double u_ref=(long double)pnt_config->u0_dim;
	long double p_ref=(long double)rho_ref*u_ref*u_ref*pnt_config->Upsilon;


	long double rho_0,rho_x,rho_y,a_rho_x,a_rho_y;
	long double u_0,u_x,u_y,a_u_x,a_u_y;
	long double v_0,v_x,v_y,a_v_x,a_v_y;
	long double p_0,p_x,p_y,a_p_x,a_p_y;

	if (pnt_config->ManufacturedSolution_case==1)
	{
	//supersonic
	rho_0=1.0/rho_ref;	rho_x=0.15/rho_ref;		rho_y=-0.1/rho_ref;		a_rho_x=1.0;	a_rho_y=0.5;
	u_0=800.0/u_ref;	u_x=50.0/u_ref;			u_y=-30.0/u_ref;		a_u_x=1.5;		a_u_y=0.6;
	v_0=800.0/u_ref;	v_x=-75.0/u_ref;		v_y=40.0/u_ref;			a_v_x=0.5;		a_v_y=2./3.;
	p_0=100000.0/p_ref;	p_x=0.2*100000.0/p_ref;	p_y=0.5*100000.0/p_ref;	a_p_x=2.0;		a_p_y=1.0;
	}
	else //if (pnt_config->ManufacturedSolution_case==0)
	{
	//subsonic
	rho_0=1.0/rho_ref;	rho_x=0.15/rho_ref;		rho_y=-0.1/rho_ref;		a_rho_x=1.0;	a_rho_y=0.5;
	u_0=70.0/u_ref;		u_x=5.0/u_ref;			u_y=-7.0/u_ref;			a_u_x=1.5;		a_u_y=0.6;
	v_0=90.0/u_ref;		v_x=-15.0/u_ref;		v_y=8.5/u_ref;			a_v_x=0.5;		a_v_y=2./3.;
	p_0=100000.0/p_ref;	p_x=0.2*100000.0/p_ref;	p_y=0.5*100000.0/p_ref;	a_p_x=2.0;		a_p_y=1.0;
	}

	pnt_U_lastStep->rho[ijk]=(rho_0+rho_x*sinl(a_rho_x*M_PIl*pnt_mesh->x[ijk]/L)+rho_y*cosl(a_rho_y*M_PIl*pnt_mesh->y[ijk]/L));
	pnt_U_lastStep->u[ijk]=(u_0+u_x*sinl(a_u_x*M_PIl*pnt_mesh->x[ijk]/L)+u_y*cosl(a_u_y*M_PIl*pnt_mesh->y[ijk]/L));
	pnt_U_lastStep->v[ijk]=(v_0+v_x*cosl(a_v_x*M_PIl*pnt_mesh->x[ijk]/L)+v_y*sinl(a_v_y*M_PIl*pnt_mesh->y[ijk]/L));
	pnt_U_lastStep->p[ijk]=(p_0+p_x*cosl(a_p_x*M_PIl*pnt_mesh->x[ijk]/L)+p_y*sinl(a_p_y*M_PIl*pnt_mesh->y[ijk]/L));

	pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+
							pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);

	pnt_U_lastStep->T[ijk]=
			fabsl(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]);

	pnt_U_lastStep->c[ijk]=
			sqrt(pnt_config->Upsilon*
			pnt_config->gammaNumber*pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]);

	pnt_U_lastStep->mue[ijk]=((1.0+pnt_config->SutherlandConstant)*powl(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk],1.5)/
			(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]+pnt_config->SutherlandConstant));

}

long double GetRhoManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		int ijk)
{
	long double rho_ref=1.0;
	long double L=1.0;
	long double rho_0,rho_x,rho_y,a_rho_x,a_rho_y;


	rho_0=1.0/rho_ref;	rho_x=0.15/rho_ref;	rho_y=-0.1/rho_ref;	a_rho_x=1.0;	a_rho_y=0.5;

	long double rho=(rho_0+rho_x*sinl(a_rho_x*M_PIl*pnt_mesh->x[ijk]/L)+rho_y*cosl(a_rho_y*M_PIl*pnt_mesh->y[ijk]/L));
	return rho;
}

long double GetPressureManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		int ijk)
{
	long double rho_ref=(long double)1.0;
	long double u_ref=(long double)pnt_config->u0_dim;
	long double p_ref=(long double)rho_ref*u_ref*u_ref*pnt_config->Upsilon;
	long double L=(long double)1.0;

	long double p_0,p_x,p_y,a_p_x,a_p_y;

	p_0=100000.0/p_ref;	p_x=0.2*100000.0/p_ref;	p_y=0.5*100000.0/p_ref;	a_p_x=2.0;		a_p_y=1.0;


	long double pressure;
	pressure=(p_0+p_x*cosl(a_p_x*M_PIl*pnt_mesh->x[ijk]/L)+p_y*sinl(a_p_y*M_PIl*pnt_mesh->y[ijk]/L));
	return pressure;
}

void ErrorManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	long double Linf_norm_rho=0.0;
	long double Linf_norm_pressure=0.0;
	long double L2_norm_rho=0.0;
	long double L2_norm_pressure=0.0;
	long double rho_exact,pressure_exact;
	long double N=pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints;
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kEndReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				rho_exact=GetRhoManufacturedSolution(
						pnt_config,
						pnt_mesh,
						pnt_U_lastStep,
						ijk);
				L2_norm_rho+=powl((rho_exact-pnt_U_lastStep->rho[ijk]),2.0)/(N);
				if(Linf_norm_rho<fabsl(rho_exact-pnt_U_lastStep->rho[ijk])){Linf_norm_rho=fabsl(rho_exact-pnt_U_lastStep->rho[ijk]);}

				pressure_exact=GetPressureManufacturedSolution(
										pnt_config,
										pnt_mesh,
										pnt_U_lastStep,
										ijk);
				L2_norm_pressure+=powl((pressure_exact-pnt_U_lastStep->p[ijk]),2.0)/(N);
				if(Linf_norm_pressure<fabsl(pressure_exact-pnt_U_lastStep->p[ijk])){Linf_norm_pressure=fabsl(pressure_exact-pnt_U_lastStep->p[ijk]);}

			}
		}
	}
	long double* all_L2_norm_rho;
	all_L2_norm_rho = (long double *)malloc(1*sizeof(long double));
	MPI_Reduce( &L2_norm_rho, all_L2_norm_rho,1,MPI_LONG_DOUBLE,MPI_SUM,0,pnt_config->MPI_comm);
	all_L2_norm_rho[0]=sqrt(all_L2_norm_rho[0]);
	if(pnt_config->MPI_rank==0){printf("SHOCK: ManufacturedSolution: L2_norm(rho):%.8le\n",
			(double)all_L2_norm_rho[0]);}

	long double* all_L2_norm_pressure;
	all_L2_norm_pressure = (long double *)malloc(1*sizeof(long double));
	MPI_Reduce( &L2_norm_pressure, all_L2_norm_pressure,1,MPI_LONG_DOUBLE,MPI_SUM,0,pnt_config->MPI_comm);
	all_L2_norm_pressure[0]=sqrt(all_L2_norm_pressure[0]);
	if(pnt_config->MPI_rank==0){printf("SHOCK: ManufacturedSolution: L2_norm(pressure):%.8le\n",
			(double)all_L2_norm_pressure[0]);}

	long double* all_Linf_norm_rho;
	all_Linf_norm_rho = (long double *)malloc(1*sizeof(long double));
	MPI_Reduce( &Linf_norm_rho, all_Linf_norm_rho,1,MPI_LONG_DOUBLE,MPI_MAX,0,pnt_config->MPI_comm);
	if(pnt_config->MPI_rank==0){printf("SHOCK: ManufacturedSolution: Linf_norm(rho):%.8le\n",
			(double)all_Linf_norm_rho[0]);}

	long double* all_Linf_norm_pressure;
	all_Linf_norm_pressure = (long double *)malloc(1*sizeof(long double));
	MPI_Reduce( &Linf_norm_pressure, all_Linf_norm_pressure,1,MPI_LONG_DOUBLE,MPI_MAX,0,pnt_config->MPI_comm);
	if(pnt_config->MPI_rank==0){printf("SHOCK: ManufacturedSolution: Linf_norm(pressure):%.8le\n",
			(double)all_Linf_norm_pressure[0]);}


	if(pnt_config->MPI_rank==0)
	{
		FILE * file0;
		char filename[200];
		sprintf(filename,"ManufacturedSolutions_W%d.dat",SPACEORDER);
		file0=fopen(filename,"a");
		fprintf(file0," %d %le %le %le %le %le\n",
				PRECISION,
				(double)(pnt_mesh->y[0]-pnt_mesh->y[1]),
				(double)all_L2_norm_rho[0],
				(double)all_L2_norm_pressure[0],
				(double)all_Linf_norm_rho[0],
				(double)all_Linf_norm_pressure[0]);
		fclose(file0);


	}

}
