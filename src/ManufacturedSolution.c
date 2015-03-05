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
	double L=1.0;
	double p_ref=100000.0;
	double u_ref=pnt_config->dbl_u0_dim;
	double rho_ref=1.0;


	double rho_0,rho_x,rho_y,a_rho_x,a_rho_y;
	double u_0,u_x,u_y,a_u_x,a_u_y;
	double v_0,v_x,v_y,a_v_x,a_v_y;
	double p_0,p_x,p_y,a_p_x,a_p_y;

	if (pnt_config->ManufacturedSolution_case==1)
	{
	//supersonic
	rho_0=1.0;	rho_x=0.15;	rho_y=-0.1;	a_rho_x=1.0;	a_rho_y=0.5;
	u_0=800.0;	u_x=50.0;	u_y=-30.0;	a_u_x=1.5;		a_u_y=0.6;
	v_0=800.0;	v_x=-75.0;	v_y=40.0;	a_v_x=0.5;		a_v_y=2./3.;
	p_0=100000.0;	p_x=0.2*100000.0;	p_y=0.5*100000.0;	a_p_x=2.0;		a_p_y=1.0;
	}
	if (pnt_config->ManufacturedSolution_case==0)
	{
	//subsonic
	rho_0=1.0;	rho_x=0.15;	rho_y=-0.1;	a_rho_x=1.0;	a_rho_y=0.5;
	u_0=70.0;	u_x=5.0;	u_y=-7.0;	a_u_x=1.5;		a_u_y=0.6;
	v_0=90.0;	v_x=-15.0;	v_y=8.5;	a_v_x=0.5;		a_v_y=2./3.;
	p_0=100000.0;	p_x=0.2*100000.0;	p_y=0.5*100000.0;	a_p_x=2.0;		a_p_y=1.0;
	}

	double mass,x_momentum,y_momentum,energy;

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				mass = (3.141592653589793*a_u_x*u_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L)))/L+(3.141592653589793*a_v_y*v_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L)))/L+(3.141592653589793*a_rho_x*rho_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L)*(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L)))/L-(3.141592653589793*a_rho_y*rho_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)*(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)))/L;

				x_momentum = -(3.141592653589793*a_p_x*p_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_p_x)/L))/L+(3.141592653589793*a_rho_x*rho_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L)*pow(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L),2.0))/L+(3.141592653589793*a_u_x*u_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))*(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L))*2.0)/L+(3.141592653589793*a_v_y*v_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))*(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L)))/L-(3.141592653589793*a_u_y*u_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))*(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)))/L-(3.141592653589793*a_rho_y*rho_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)*(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L))*(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)))/L;

				y_momentum = (3.141592653589793*a_p_y*p_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_p_y)/L))/L-(3.141592653589793*a_rho_y*rho_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)*pow(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L),2.0))/L+(3.141592653589793*a_u_x*u_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))*(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)))/L+(3.141592653589793*a_v_y*v_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))*(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L))*2.0)/L+(3.141592653589793*a_rho_x*rho_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L)*(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L))*(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)))/L-(3.141592653589793*a_v_x*v_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))*(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L)))/L;

				energy = -(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L))*((rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))*(-(3.141592653589793*a_u_x*u_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L)*(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L)))/L+(3.141592653589793*a_v_x*v_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)*(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)))/L+(3.141592653589793*a_p_x*p_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_p_x)/L))/(L*(pnt_config->dbl_gammaNumber-1.0)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L)))+(3.141592653589793*a_rho_x*rho_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L)*(p_0+p_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_p_x)/L)+p_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_p_y)/L))*1.0/pow(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L),2.0))/(L*(pnt_config->dbl_gammaNumber-1.0)))+(3.141592653589793*a_p_x*p_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_p_x)/L))/L-(3.141592653589793*a_rho_x*rho_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L)*(pow(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L),2.0)*(1.0/2.0)+pow(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L),2.0)*(1.0/2.0)+(p_0+p_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_p_x)/L)+p_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_p_y)/L))/((pnt_config->dbl_gammaNumber-1.0)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L)))))/L)+(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L))*((rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))*((3.141592653589793*a_v_y*v_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)*(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)))/L-(3.141592653589793*a_u_y*u_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)*(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L)))/L+(3.141592653589793*a_p_y*p_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_p_y)/L))/(L*(pnt_config->dbl_gammaNumber-1.0)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L)))+(3.141592653589793*a_rho_y*rho_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)*(p_0+p_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_p_x)/L)+p_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_p_y)/L))*1.0/pow(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L),2.0))/(L*(pnt_config->dbl_gammaNumber-1.0)))+(3.141592653589793*a_p_y*p_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_p_y)/L))/L-(3.141592653589793*a_rho_y*rho_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)*(pow(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L),2.0)*(1.0/2.0)+pow(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L),2.0)*(1.0/2.0)+(p_0+p_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_p_x)/L)+p_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_p_y)/L))/((pnt_config->dbl_gammaNumber-1.0)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L)))))/L)+(3.141592653589793*a_u_x*u_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L)*(p_0+p_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_p_x)/L)+p_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_p_y)/L)+(pow(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L),2.0)*(1.0/2.0)+pow(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L),2.0)*(1.0/2.0)+(p_0+p_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_p_x)/L)+p_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_p_y)/L))/((pnt_config->dbl_gammaNumber-1.0)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))))*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))))/L+(3.141592653589793*a_v_y*v_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L)*(p_0+p_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_p_x)/L)+p_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_p_y)/L)+(pow(u_0+u_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_u_y)/L)+u_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_u_x)/L),2.0)*(1.0/2.0)+pow(v_0+v_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_v_x)/L)+v_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_v_y)/L),2.0)*(1.0/2.0)+(p_0+p_x*cos((pnt_mesh->x[ijk]*3.141592653589793*a_p_x)/L)+p_y*sin((pnt_mesh->y[ijk]*3.141592653589793*a_p_y)/L))/((pnt_config->dbl_gammaNumber-1.0)*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))))*(rho_0+rho_y*cos((pnt_mesh->y[ijk]*3.141592653589793*a_rho_y)/L)+rho_x*sin((pnt_mesh->x[ijk]*3.141592653589793*a_rho_x)/L))))/L;



				pnt_Q->Mass[ijk]+=pnt_mesh->jacobian[ijk]*mass/(rho_ref*u_ref/L);

				pnt_Q->xiMomentum[ijk]+=pnt_mesh->jacobian[ijk]*x_momentum/(rho_ref*u_ref*u_ref/L);

				pnt_Q->etaMomentum[ijk]+=pnt_mesh->jacobian[ijk]*y_momentum/(rho_ref*u_ref*u_ref/L);

				pnt_Q->Energy[ijk]+=pnt_mesh->jacobian[ijk]*energy/(u_ref*p_ref/L)*pnt_config->dbl_Upsilon;
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
	int ijkSymmetry, int_symmetryIndex;
	int_symmetryIndex=pnt_config->int_iStartReal;

	for (i=pnt_config->int_iStartReal-1; i >= pnt_config->int_iStartGhosts; i--)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

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
	int ijkSymmetry, int_symmetryIndex;
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		int_symmetryIndex=pnt_config->int_jStartReal;
		for (j=pnt_config->int_jStartReal-1; j >= pnt_config->int_jStartGhosts; j--)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

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
	int ijkSymmetry, int_symmetryIndex;
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			int_symmetryIndex=pnt_config->int_kStartReal;
			for (k=pnt_config->int_kStartReal-1; k >= pnt_config->int_kStartGhosts; k--)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

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
	int ijkSymmetry, int_symmetryIndex;
	int_symmetryIndex=pnt_config->int_iEndReal;
	for (i=pnt_config->int_iEndReal+1; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=int_symmetryIndex*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

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
	int ijkSymmetry, int_symmetryIndex;
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		int_symmetryIndex=pnt_config->int_jEndReal;
		for (j=pnt_config->int_jEndReal+1; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex*pnt_config->int_kMeshPointsGhostCells+k;

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
	int ijkSymmetry, int_symmetryIndex;
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			int_symmetryIndex=pnt_config->int_kEndReal;
			for (k=pnt_config->int_kEndReal+1; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkSymmetry=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+int_symmetryIndex;

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
	double L=1.0;
	double p_ref=100000.0;
	double u_ref=pnt_config->dbl_u0_dim;
	double rho_ref=1.0;

	double rho_0,rho_x,rho_y,a_rho_x,a_rho_y;
	double u_0,u_x,u_y,a_u_x,a_u_y;
	double v_0,v_x,v_y,a_v_x,a_v_y;
	double p_0,p_x,p_y,a_p_x,a_p_y;

	if (pnt_config->ManufacturedSolution_case==1)
	{
	//supersonic
	rho_0=1.0;	rho_x=0.15;	rho_y=-0.1;	a_rho_x=1.0;	a_rho_y=0.5;
	u_0=800.0;	u_x=50.0;	u_y=-30.0;	a_u_x=1.5;		a_u_y=0.6;
	v_0=800.0;	v_x=-75.0;	v_y=40.0;	a_v_x=0.5;		a_v_y=2./3.;
	p_0=100000.0;	p_x=0.2*100000.0;	p_y=0.5*100000.0;	a_p_x=2.0;		a_p_y=1.0;
	}
	if (pnt_config->ManufacturedSolution_case==0)
	{
	//subsonic
	rho_0=1.0;	rho_x=0.15;	rho_y=-0.1;	a_rho_x=1.0;	a_rho_y=0.5;
	u_0=70.0;	u_x=5.0;	u_y=-7.0;	a_u_x=1.5;		a_u_y=0.6;
	v_0=90.0;	v_x=-15.0;	v_y=8.5;	a_v_x=0.5;		a_v_y=2./3.;
	p_0=100000.0;	p_x=0.2*100000.0;	p_y=0.5*100000.0;	a_p_x=2.0;		a_p_y=1.0;
	}

	pnt_U_lastStep->rho[ijk]=1./rho_ref*(rho_0+rho_x*sin(a_rho_x*M_PI*pnt_mesh->x[ijk]/L)+rho_y*cos(a_rho_y*M_PI*pnt_mesh->y[ijk]/L));
	pnt_U_lastStep->u[ijk]=1./u_ref*(u_0+u_x*sin(a_u_x*M_PI*pnt_mesh->x[ijk]/L)+u_y*cos(a_u_y*M_PI*pnt_mesh->y[ijk]/L));
	pnt_U_lastStep->v[ijk]=1./u_ref*(v_0+v_x*cos(a_v_x*M_PI*pnt_mesh->x[ijk]/L)+v_y*sin(a_v_y*M_PI*pnt_mesh->y[ijk]/L));
	pnt_U_lastStep->p[ijk]=1./p_ref*(p_0+p_x*cos(a_p_x*M_PI*pnt_mesh->x[ijk]/L)+p_y*sin(a_p_y*M_PI*pnt_mesh->y[ijk]/L));

	pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+
							pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->dbl_gammaNumber-1.0)*pnt_config->dbl_Upsilon);

	pnt_U_lastStep->T[ijk]=
			fabs(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]);

	pnt_U_lastStep->c[ijk]=
			sqrt(pnt_config->dbl_Upsilon*
			pnt_config->dbl_gammaNumber*pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]);

	pnt_U_lastStep->mue[ijk]=((1.0+pnt_config->dbl_SutherlandConstant)*pow(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk],1.5)/
			(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]+pnt_config->dbl_SutherlandConstant));

}

double GetRhoManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		int ijk)
{
	double L=pnt_config->dbl_L0_dim;
	double rho_ref=1.0;

	double rho_0,rho_x,rho_y,a_rho_x,a_rho_y;


	rho_0=1.0;	rho_x=0.15;	rho_y=-0.1;	a_rho_x=1.0;	a_rho_y=0.5;

	double rho=1./rho_ref*(rho_0+rho_x*sin(a_rho_x*M_PI*pnt_mesh->x[ijk]/L)+rho_y*cos(a_rho_y*M_PI*pnt_mesh->y[ijk]/L));
	return rho;
}

double GetPressureManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		int ijk)
{
	double L=pnt_config->dbl_L0_dim;
	double p_ref=100000.0;
	double p_0,p_x,p_y,a_p_x,a_p_y;

	p_0=100000.0;	p_x=0.2*100000.0;	p_y=0.5*100000.0;	a_p_x=2.0;		a_p_y=1.0;


	double pressure;
	pressure=1./p_ref*(p_0+p_x*cos(a_p_x*M_PI*pnt_mesh->x[ijk]/L)+p_y*sin(a_p_y*M_PI*pnt_mesh->y[ijk]/L));
	return pressure;
}

void ErrorManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep)
{
	double L2_norm_rho=0.0;
	double L2_norm_pressure=0.0;
	double rho_exact,pressure_exact;
	double N=pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints;
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
				L2_norm_rho+=pow((rho_exact-pnt_U_lastStep->rho[ijk]),2.0)/(N);

				pressure_exact=GetPressureManufacturedSolution(
										pnt_config,
										pnt_mesh,
										pnt_U_lastStep,
										ijk);
				L2_norm_pressure+=pow((pressure_exact-pnt_U_lastStep->p[ijk]),2.0)/(N);

			}
		}
	}
	double* all_L2_norm_rho;
	all_L2_norm_rho = (double *)malloc(1*sizeof(double));
	MPI_Reduce( &L2_norm_rho, all_L2_norm_rho,1,MPI_DOUBLE,MPI_SUM,0,pnt_config->MPI_comm);
	all_L2_norm_rho[0]=sqrt(all_L2_norm_rho[0]);
	if(pnt_config->MPI_rank==0){printf("SHOCK: ManufacturedSolution: L2_norm(rho):%.8le\n",all_L2_norm_rho[0]);}

	double* all_L2_norm_pressure;
	all_L2_norm_pressure = (double *)malloc(1*sizeof(double));
	MPI_Reduce( &L2_norm_pressure, all_L2_norm_pressure,1,MPI_DOUBLE,MPI_SUM,0,pnt_config->MPI_comm);
	all_L2_norm_pressure[0]=sqrt(all_L2_norm_pressure[0]);
	if(pnt_config->MPI_rank==0){printf("SHOCK: ManufacturedSolution: L2_norm(pressure):%.8le\n",all_L2_norm_pressure[0]);}
}