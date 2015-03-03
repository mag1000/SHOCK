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
	float L=pnt_config->flt_L0_dim;
	float p_ref=100000.0;
	float u_ref=pnt_config->flt_u0_dim;
	float rho_ref=p_ref/287./pnt_config->flt_T0_dim;


	float rho_0,rho_x,rho_y,a_rho_x,a_rho_y;
	float u_0,u_x,u_y,a_u_x,a_u_y;
	float v_0,v_x,v_y,a_v_x,a_v_y;
	float p_0,p_x,p_y,a_p_x,a_p_y;

	rho_0=1.0;	rho_x=0.15;	rho_y=-0.1;	a_rho_x=1.0;	a_rho_y=0.5;
	u_0=800.0;	u_x=50.0;	u_y=-30.0;	a_u_x=1.5;		a_u_y=0.6;
	v_0=800.0;	v_x=-75.0;	v_y=40.0;	a_v_x=0.5;		a_v_y=2./3.;
	p_0=100000.0;	p_x=0.2*100000.0;	p_y=0.5*100000.0;	a_p_x=2.0;		a_p_y=1.0;



	float mass,x_momentum,y_momentum,energy;

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				mass = (M_PI*a_u_x*u_x*cos((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L)))/L + (M_PI*a_v_y*v_y*cos((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L)))/L + (M_PI*a_rho_x*rho_x*cos((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L)*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)))/L - (M_PI*a_rho_y*rho_y*sin((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L)*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)))/L;

				x_momentum = (M_PI*a_rho_x*rho_x*cos((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L)*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)))/L - (M_PI*a_p_x*pnt_config->flt_Upsilon*p_x*sin((M_PI*pnt_mesh->x[ijk]*a_p_x)/L))/L + (2*M_PI*a_u_x*u_x*cos((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)))/L + (M_PI*a_v_y*v_y*cos((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)))/L - (M_PI*a_u_y*u_y*sin((M_PI*pnt_mesh->y[ijk]*a_u_y)/L)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)))/L - (M_PI*a_rho_y*rho_y*sin((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L)*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)))/L;

				y_momentum =(M_PI*a_p_y*pnt_config->flt_Upsilon*p_y*cos((M_PI*pnt_mesh->y[ijk]*a_p_y)/L))/L - (M_PI*a_rho_y*rho_y*sin((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L)*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L))*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)))/L + (M_PI*a_u_x*u_x*cos((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)))/L + (2*M_PI*a_v_y*v_y*cos((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)))/L + (M_PI*a_rho_x*rho_x*cos((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L)*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)))/L - (M_PI*a_v_x*v_x*sin((M_PI*pnt_mesh->x[ijk]*a_v_x)/L)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)))/L;

				energy = (v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L))*((rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))*((M_PI*a_v_y*v_y*cos((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)))/L - (M_PI*a_u_y*u_y*sin((M_PI*pnt_mesh->y[ijk]*a_u_y)/L)*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)))/L + (M_PI*a_p_y*p_y*cos((M_PI*pnt_mesh->y[ijk]*a_p_y)/L))/(L*pnt_config->flt_machNumber*pnt_config->flt_machNumber*pnt_config->flt_gammaNumber*(pnt_config->flt_gammaNumber - 1)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))) + (M_PI*a_rho_y*rho_y*sin((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L)*(p_0 + p_x*cos((M_PI*pnt_mesh->x[ijk]*a_p_x)/L) + p_y*sin((M_PI*pnt_mesh->y[ijk]*a_p_y)/L)))/(L*pnt_config->flt_machNumber*pnt_config->flt_machNumber*pnt_config->flt_gammaNumber*(pnt_config->flt_gammaNumber - 1)*(rho_0 + rho_y*cos((pnt_mesh->y[ijk]*M_PI*a_rho_y)/L) + rho_x*sin((pnt_mesh->x[ijk]*M_PI*a_rho_x)/L))*(rho_0 + rho_y*cos((pnt_mesh->y[ijk]*M_PI*a_rho_y)/L) + rho_x*sin((pnt_mesh->x[ijk]*M_PI*a_rho_x)/L)))) - (M_PI*a_rho_y*rho_y*sin((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L)*((u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))/2 + (v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L))*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L))/2 + (p_0 + p_x*cos((M_PI*pnt_mesh->x[ijk]*a_p_x)/L) + p_y*sin((M_PI*pnt_mesh->y[ijk]*a_p_y)/L))/(pnt_config->flt_machNumber*pnt_config->flt_machNumber*pnt_config->flt_gammaNumber*(pnt_config->flt_gammaNumber - 1)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L)))))/L + (M_PI*a_p_y*pnt_config->flt_Upsilon*p_y*cos((M_PI*pnt_mesh->y[ijk]*a_p_y)/L))/L) - (u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))*((rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))*((M_PI*a_v_x*v_x*sin((M_PI*pnt_mesh->x[ijk]*a_v_x)/L)*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)))/L - (M_PI*a_u_x*u_x*cos((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)))/L + (M_PI*a_p_x*p_x*sin((M_PI*pnt_mesh->x[ijk]*a_p_x)/L))/(L*pnt_config->flt_machNumber*pnt_config->flt_machNumber*pnt_config->flt_gammaNumber*(pnt_config->flt_gammaNumber - 1)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))) + (M_PI*a_rho_x*rho_x*cos((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L)*(p_0 + p_x*cos((M_PI*pnt_mesh->x[ijk]*a_p_x)/L) + p_y*sin((M_PI*pnt_mesh->y[ijk]*a_p_y)/L)))/(L*pnt_config->flt_machNumber*pnt_config->flt_machNumber*pnt_config->flt_gammaNumber*(pnt_config->flt_gammaNumber - 1)*(rho_0 + rho_y*cos((pnt_mesh->y[ijk]*M_PI*a_rho_y)/L) + rho_x*sin((pnt_mesh->x[ijk]*M_PI*a_rho_x)/L))*(rho_0 + rho_y*cos((pnt_mesh->y[ijk]*M_PI*a_rho_y)/L) + rho_x*sin((pnt_mesh->x[ijk]*M_PI*a_rho_x)/L)))) - (M_PI*a_rho_x*rho_x*cos((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L)*((u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))/2 + (v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L))*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L))/2 + (p_0 + p_x*cos((M_PI*pnt_mesh->x[ijk]*a_p_x)/L) + p_y*sin((M_PI*pnt_mesh->y[ijk]*a_p_y)/L))/(pnt_config->flt_machNumber*pnt_config->flt_machNumber*pnt_config->flt_gammaNumber*(pnt_config->flt_gammaNumber - 1)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L)))))/L + (M_PI*a_p_x*pnt_config->flt_Upsilon*p_x*sin((M_PI*pnt_mesh->x[ijk]*a_p_x)/L))/L) + (M_PI*a_u_x*u_x*cos((M_PI*pnt_mesh->x[ijk]*a_u_x)/L)*(pnt_config->flt_Upsilon*(p_0 + p_x*cos((M_PI*pnt_mesh->x[ijk]*a_p_x)/L) + p_y*sin((M_PI*pnt_mesh->y[ijk]*a_p_y)/L)) + (rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))*((u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))/2 + (v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L))*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L))/2 + (p_0 + p_x*cos((M_PI*pnt_mesh->x[ijk]*a_p_x)/L) + p_y*sin((M_PI*pnt_mesh->y[ijk]*a_p_y)/L))/(pnt_config->flt_machNumber*pnt_config->flt_machNumber*pnt_config->flt_gammaNumber*(pnt_config->flt_gammaNumber - 1)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))))))/L + (M_PI*a_v_y*v_y*cos((M_PI*pnt_mesh->y[ijk]*a_v_y)/L)*(pnt_config->flt_Upsilon*(p_0 + p_x*cos((M_PI*pnt_mesh->x[ijk]*a_p_x)/L) + p_y*sin((M_PI*pnt_mesh->y[ijk]*a_p_y)/L)) + (rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))*((u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))*(u_0 + u_y*cos((M_PI*pnt_mesh->y[ijk]*a_u_y)/L) + u_x*sin((M_PI*pnt_mesh->x[ijk]*a_u_x)/L))/2 + (v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L))*(v_0 + v_x*cos((M_PI*pnt_mesh->x[ijk]*a_v_x)/L) + v_y*sin((M_PI*pnt_mesh->y[ijk]*a_v_y)/L))/2 + (p_0 + p_x*cos((M_PI*pnt_mesh->x[ijk]*a_p_x)/L) + p_y*sin((M_PI*pnt_mesh->y[ijk]*a_p_y)/L))/(pnt_config->flt_machNumber*pnt_config->flt_machNumber*pnt_config->flt_gammaNumber*(pnt_config->flt_gammaNumber - 1)*(rho_0 + rho_y*cos((M_PI*pnt_mesh->y[ijk]*a_rho_y)/L) + rho_x*sin((M_PI*pnt_mesh->x[ijk]*a_rho_x)/L))))))/L;



				pnt_Q->Mass[ijk]+=pnt_mesh->jacobian[ijk]*mass/rho_ref;

				pnt_Q->xiMomentum[ijk]+=pnt_mesh->jacobian[ijk]*x_momentum/(rho_ref*u_ref);

				pnt_Q->etaMomentum[ijk]+=pnt_mesh->jacobian[ijk]*y_momentum/(rho_ref*u_ref);

				pnt_Q->Energy[ijk]+=pnt_mesh->jacobian[ijk]*energy/(rho_ref*u_ref*u_ref);
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
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
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
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
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
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
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
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
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
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
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
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
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
	float L=pnt_config->flt_L0_dim;
	float p_ref=100000.0;
	float u_ref=pnt_config->flt_u0_dim;
	float rho_ref=p_ref/287./pnt_config->flt_T0_dim;

	float rho_0,rho_x,rho_y,a_rho_x,a_rho_y;
	float u_0,u_x,u_y,a_u_x,a_u_y;
	float v_0,v_x,v_y,a_v_x,a_v_y;
	float p_0,p_x,p_y,a_p_x,a_p_y;

	rho_0=1.0;	rho_x=0.15;	rho_y=-0.1;	a_rho_x=1.0;	a_rho_y=0.5;
	u_0=800.0;	u_x=50.0;	u_y=-30.0;	a_u_x=1.5;		a_u_y=0.6;
	v_0=800.0;	v_x=-75.0;	v_y=40.0;	a_v_x=0.5;		a_v_y=2./3.;
	p_0=100000.0;	p_x=0.2*100000.0;	p_y=0.5*100000.0;	a_p_x=2.0;		a_p_y=1.0;

	pnt_U_lastStep->rho[ijk]=1./rho_ref*(rho_0+rho_x*sin(a_rho_x*M_PI*pnt_mesh->x[ijk]/L)+rho_y*cos(a_rho_y*M_PI*pnt_mesh->y[ijk]/L));
	pnt_U_lastStep->u[ijk]=1./u_ref*(u_0+u_x*sin(a_u_x*M_PI*pnt_mesh->x[ijk]/L)+u_y*cos(a_u_y*M_PI*pnt_mesh->y[ijk]/L));
	pnt_U_lastStep->v[ijk]=1./u_ref*(v_0+v_x*cos(a_v_x*M_PI*pnt_mesh->x[ijk]/L)+v_y*sin(a_v_y*M_PI*pnt_mesh->y[ijk]/L));
	pnt_U_lastStep->p[ijk]=1./p_ref*(p_0+p_x*cos(a_p_x*M_PI*pnt_mesh->x[ijk]/L)+p_y*sin(a_p_y*M_PI*pnt_mesh->y[ijk]/L));
	pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+
							pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->flt_gammaNumber-1.0)*pnt_config->flt_Upsilon);
}
