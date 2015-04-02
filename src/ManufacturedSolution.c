// Manufactured Solutions
// These functions are needed for the implementation of the Manufactured Solution Approach
// in which an analytical method is used to verify the code. It can be used to ensure that the implementation is bug-free
// This is done by comparing the numerical approximated solution to the exact solution of the governing equations.
// The exact solution must be derived a priori and then, this is added as a source term. The resulting exact analytical solution
// is known and the difference to the numerical solution is the numerical error.
// This error must decrease with the order of the code when the code is bug-free

#define MY_PI_L 3.14159265358979323846264338327950288419716939937510L

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
	long double L=1.0L;

	long double rho_0,rho_x,rho_y,rho_z,a_rho_x,a_rho_y,a_rho_z;
	long double u_0,u_x,u_y,u_z,a_u_x,a_u_y,a_u_z;
	long double v_0,v_x,v_y,v_z,a_v_x,a_v_y,a_v_z;
	long double w_0,w_x,w_y,w_z,a_w_x,a_w_y,a_w_z;
	long double p_0,p_x,p_y,p_z,a_p_x,a_p_y,a_p_z;

	rho_0=pnt_config->ManufacturedSolution_param_rho[0];
	rho_x=pnt_config->ManufacturedSolution_param_rho[1];
	rho_y=pnt_config->ManufacturedSolution_param_rho[2];
	rho_z=pnt_config->ManufacturedSolution_param_rho[3];
	a_rho_x=pnt_config->ManufacturedSolution_param_rho[4];
	a_rho_y=pnt_config->ManufacturedSolution_param_rho[5];
	a_rho_z=pnt_config->ManufacturedSolution_param_rho[6];

	u_0=pnt_config->ManufacturedSolution_param_u[0];
	u_x=pnt_config->ManufacturedSolution_param_u[1];
	u_y=pnt_config->ManufacturedSolution_param_u[2];
	u_z=pnt_config->ManufacturedSolution_param_u[3];
	a_u_x=pnt_config->ManufacturedSolution_param_u[4];
	a_u_y=pnt_config->ManufacturedSolution_param_u[5];
	a_u_z=pnt_config->ManufacturedSolution_param_u[6];

	v_0=pnt_config->ManufacturedSolution_param_v[0];
	v_x=pnt_config->ManufacturedSolution_param_v[1];
	v_y=pnt_config->ManufacturedSolution_param_v[2];
	v_z=pnt_config->ManufacturedSolution_param_v[3];
	a_v_x=pnt_config->ManufacturedSolution_param_v[4];
	a_v_y=pnt_config->ManufacturedSolution_param_v[5];
	a_v_z=pnt_config->ManufacturedSolution_param_v[6];

	w_0=pnt_config->ManufacturedSolution_param_w[0];
	w_x=pnt_config->ManufacturedSolution_param_w[1];
	w_y=pnt_config->ManufacturedSolution_param_w[2];
	w_z=pnt_config->ManufacturedSolution_param_w[3];
	a_w_x=pnt_config->ManufacturedSolution_param_w[4];
	a_w_y=pnt_config->ManufacturedSolution_param_w[5];
	a_w_z=pnt_config->ManufacturedSolution_param_w[6];

	p_0=pnt_config->ManufacturedSolution_param_p[0];
	p_x=pnt_config->ManufacturedSolution_param_p[1];
	p_y=pnt_config->ManufacturedSolution_param_p[2];
	p_z=pnt_config->ManufacturedSolution_param_p[3];
	a_p_x=pnt_config->ManufacturedSolution_param_p[4];
	a_p_y=pnt_config->ManufacturedSolution_param_p[5];
	a_p_z=pnt_config->ManufacturedSolution_param_p[6];


	long double mass,x_momentum,y_momentum,z_momentum,energy;
	long double energy_flux_1,energy_flux_2,energy_flux_3,energy_flux_4;
	long double CoordX,CoordY,CoordZ,Upsilon,Psi,Gamma,gammaNumber,mue;



	Upsilon=(long double)pnt_config->Upsilon;
	Psi=(long double)pnt_config->Psi;
	Gamma=(long double)(1.0L/(((long double)pnt_config->gammaNumber-1.0L)*powl((long double)pnt_config->machNumber,2.0L)*(long double)pnt_config->reynoldsNumber*(long double)pnt_config->prandtlNumber));
	gammaNumber=(long double)pnt_config->gammaNumber;

	mue=1.0L/pnt_config->reynoldsNumber;
	//long double value_old,value_new;


	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				CoordX=(long double)pnt_mesh->x[ijk];
				CoordY=(long double)pnt_mesh->y[ijk];


#if MESHDIMENSIONS==2
				if (pnt_config->flag_Inviscid==1)
				{
					mass = (MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))/L+(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L;
					x_momentum = (MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0))/L-(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/L+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*2.0)/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L-(MY_PI_L*a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L;
					y_momentum = -(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L),2.0))/L+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/L+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L))*2.0)/L+(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L-(MY_PI_L*a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L;
					energy = -(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(-(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L+(MY_PI_L*a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L+(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/(L*(gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))+(MY_PI_L*Upsilon*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L),2.0))/(L*(gammaNumber-1.0)))-(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0)*(1.0/2.0)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L),2.0)*(1.0/2.0)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))))/L+(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/L)+(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L))*((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*((MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L-(MY_PI_L*a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/(L*(gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))+(MY_PI_L*Upsilon*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L),2.0))/(L*(gammaNumber-1.0)))+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/L-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0)*(1.0/2.0)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L),2.0)*(1.0/2.0)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))))/L)+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0)*(1.0/2.0)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L),2.0)*(1.0/2.0)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))))+Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0)*(1.0/2.0)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L),2.0)*(1.0/2.0)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))))+Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))))/L;
				}
				else
				{
					mass = (MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))/L+(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L;
					x_momentum = (MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L))/L-(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/L+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_u_y*a_u_y)*mue*u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_u_x*a_u_x)*mue*u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)*(4.0L/3.0L)+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*2.0L)/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L-(MY_PI_L*a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L;
					y_momentum = -(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L),2.0L))/L+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/L+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_v_x*a_v_x)*mue*v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_v_y*a_v_y)*mue*v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)*(4.0L/3.0L)+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L))*2.0L)/L+(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L-(MY_PI_L*a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L;
					energy = -(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*(-(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L+(MY_PI_L*a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L+(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/(L*(gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))+(MY_PI_L*Upsilon*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L),2.0L))/(L*(gammaNumber-1.0L)))-(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))))/L+(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/L)+(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L))*((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))*((MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)))/L-(MY_PI_L*a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/(L*(gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))+(MY_PI_L*Upsilon*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L),2.0L))/(L*(gammaNumber-1.0L)))-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)))))/L+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/L)-Gamma*mue*(-(1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_p_x*a_p_x)*p_x*cosl((CoordX*MY_PI_L*a_p_x)/L))/(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_rho_x*a_rho_x)*rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L),2.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_rho_x*a_rho_x)*(rho_x*rho_x)*powl(cosl((CoordX*MY_PI_L*a_rho_x)/L),2.0L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L),3.0L)*2.0L+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*a_p_x*a_rho_x*p_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*sinl((CoordX*MY_PI_L*a_p_x)/L)*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L),2.0L)*2.0L)-Gamma*mue*(-(1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_p_y*a_p_y)*p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))/(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_rho_y*a_rho_y)*(rho_y*rho_y)*powl(sinl((CoordY*MY_PI_L*a_rho_y)/L),2.0L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L),3.0L)*2.0L+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_rho_y*a_rho_y)*rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L),2.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*a_p_y*a_rho_y*p_y*rho_y*cosl((CoordY*MY_PI_L*a_p_y)/L)*sinl((CoordY*MY_PI_L*a_rho_y)/L)*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L),2.0L)*2.0L)+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*((powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))+Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*((powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L))+Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))))/L-1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_u_x*mue*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*2.0L-a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L))*(2.0L/3.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_v_y*mue*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)-a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*2.0L)*(2.0L/3.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_u_y*a_u_y)*mue*u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_v_x*a_v_x)*mue*v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L))-1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_u_y*mue*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)*(a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)+a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L))-1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_v_x*mue*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L)*(a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)+a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_u_x*a_u_x)*mue*u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(4.0L/3.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_v_y*a_v_y)*mue*v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L))*(4.0L/3.0L);
				}
#endif

#if MESHDIMENSIONS==3
				CoordZ=(long double)pnt_mesh->z[ijk];
				mass = (MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L+(MY_PI_L*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L-(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L;

				if (pnt_config->flag_Inviscid==1)
				{
					x_momentum = (MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L))/L-(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/L+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*2.0L)/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L+(MY_PI_L*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L-(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L-(MY_PI_L*a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L-(MY_PI_L*a_u_z*u_z*sinl((CoordZ*MY_PI_L*a_u_z)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L;
					y_momentum = -(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L))/L+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/L+(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L-(MY_PI_L*a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L))*2.0L)/L+(MY_PI_L*a_v_z*v_z*cosl((CoordZ*MY_PI_L*a_v_z)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L-(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L;
					z_momentum = (MY_PI_L*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L))/L-(MY_PI_L*Upsilon*a_p_z*p_z*sinl((CoordZ*MY_PI_L*a_p_z)/L))/L+(MY_PI_L*a_w_x*w_x*cosl((CoordX*MY_PI_L*a_w_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L+(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*a_w_y*w_y*cosl((CoordY*MY_PI_L*a_w_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L-(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L))*2.0L)/L-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L;
					energy = -((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(-(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L-(MY_PI_L*a_w_x*w_x*cosl((CoordX*MY_PI_L*a_w_x)/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L+(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/(L*(gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))+(MY_PI_L*Upsilon*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L))/(L*(gammaNumber-1.0L)))-(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))))/L+(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))+((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(-(MY_PI_L*a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L+(MY_PI_L*a_w_y*w_y*cosl((CoordY*MY_PI_L*a_w_y)/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/(L*(gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))+(MY_PI_L*Upsilon*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L))/(L*(gammaNumber-1.0L)))-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))))/L+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L))-((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*((MY_PI_L*a_u_z*u_z*sinl((CoordZ*MY_PI_L*a_u_z)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L-(MY_PI_L*a_v_z*v_z*cosl((CoordZ*MY_PI_L*a_v_z)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L+(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*Upsilon*a_p_z*p_z*sinl((CoordZ*MY_PI_L*a_p_z)/L))/(L*(gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))+(MY_PI_L*Upsilon*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L))/(L*(gammaNumber-1.0L)))-(MY_PI_L*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))))/L+(MY_PI_L*Upsilon*a_p_z*p_z*sinl((CoordZ*MY_PI_L*a_p_z)/L))/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L))+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))+(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))))))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))+(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))))))/L-(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))+(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))))))/L;
				}
				else
				{
					mass = (MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L+(MY_PI_L*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L-(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L;
					x_momentum = (MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L))/L-(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/L+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_u_y*a_u_y)*mue*u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_u_z*a_u_z)*mue*u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_u_x*a_u_x)*mue*u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)*(4.0L/3.0L)+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*2.0L)/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L+(MY_PI_L*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L-(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L-(MY_PI_L*a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L-(MY_PI_L*a_u_z*u_z*sinl((CoordZ*MY_PI_L*a_u_z)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L;
					y_momentum = -(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L))/L+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/L+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_v_x*a_v_x)*mue*v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_v_y*a_v_y)*mue*v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)*(4.0L/3.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_v_z*a_v_z)*mue*v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)+(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L-(MY_PI_L*a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L))*2.0L)/L+(MY_PI_L*a_v_z*v_z*cosl((CoordZ*MY_PI_L*a_v_z)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L-(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L;
					z_momentum = (MY_PI_L*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L))/L-(MY_PI_L*Upsilon*a_p_z*p_z*sinl((CoordZ*MY_PI_L*a_p_z)/L))/L+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_w_z*a_w_z)*mue*w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)*(4.0L/3.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_w_x*a_w_x)*mue*w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_w_y*a_w_y)*mue*w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)+(MY_PI_L*a_w_x*w_x*cosl((CoordX*MY_PI_L*a_w_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))/L+(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*a_w_y*w_y*cosl((CoordY*MY_PI_L*a_w_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L-(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L))*2.0L)/L-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L))*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L;
					energy = -((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(-(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L-(MY_PI_L*a_w_x*w_x*cosl((CoordX*MY_PI_L*a_w_x)/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L+(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/(L*(gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))+(MY_PI_L*Upsilon*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L))/(L*(gammaNumber-1.0L)))-(MY_PI_L*a_rho_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))))/L+(MY_PI_L*Upsilon*a_p_x*p_x*sinl((CoordX*MY_PI_L*a_p_x)/L))/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))+((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(-(MY_PI_L*a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L+(MY_PI_L*a_w_y*w_y*cosl((CoordY*MY_PI_L*a_w_y)/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/(L*(gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))+(MY_PI_L*Upsilon*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L))/(L*(gammaNumber-1.0L)))-(MY_PI_L*a_rho_y*rho_y*sinl((CoordY*MY_PI_L*a_rho_y)/L)*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))))/L+(MY_PI_L*Upsilon*a_p_y*p_y*cosl((CoordY*MY_PI_L*a_p_y)/L))/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L))-((rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*((MY_PI_L*a_u_z*u_z*sinl((CoordZ*MY_PI_L*a_u_z)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)))/L-(MY_PI_L*a_v_z*v_z*cosl((CoordZ*MY_PI_L*a_v_z)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)))/L+(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)))/L+(MY_PI_L*Upsilon*a_p_z*p_z*sinl((CoordZ*MY_PI_L*a_p_z)/L))/(L*(gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))+(MY_PI_L*Upsilon*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L))/(L*(gammaNumber-1.0L)))-(MY_PI_L*a_rho_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)))))/L+(MY_PI_L*Upsilon*a_p_z*p_z*sinl((CoordZ*MY_PI_L*a_p_z)/L))/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L))-Gamma*mue*(-(1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_p_x*a_p_x)*p_x*cosl((CoordX*MY_PI_L*a_p_x)/L))/(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_rho_x*a_rho_x)*(rho_x*rho_x)*powl(cosl((CoordX*MY_PI_L*a_rho_x)/L),2.0L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),3.0L)*2.0L+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_rho_x*a_rho_x)*rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*a_p_x*a_rho_x*p_x*rho_x*cosl((CoordX*MY_PI_L*a_rho_x)/L)*sinl((CoordX*MY_PI_L*a_p_x)/L)*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L)*2.0L)-Gamma*mue*(-(1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_p_z*a_p_z)*p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L))/(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_rho_z*a_rho_z)*(rho_z*rho_z)*powl(cosl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),3.0L)*2.0L+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_rho_z*a_rho_z)*rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*a_p_z*a_rho_z*p_z*rho_z*cosl((CoordZ*MY_PI_L*a_rho_z)/L)*sinl((CoordZ*MY_PI_L*a_p_z)/L)*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L)*2.0L)-Gamma*mue*(-(1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_p_y*a_p_y)*p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))/(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_rho_y*a_rho_y)*(rho_y*rho_y)*powl(sinl((CoordY*MY_PI_L*a_rho_y)/L),2.0L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),3.0L)*2.0L+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*(a_rho_y*a_rho_y)*rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*a_p_y*a_rho_y*p_y*rho_y*cosl((CoordY*MY_PI_L*a_p_y)/L)*sinl((CoordY*MY_PI_L*a_rho_y)/L)*1.0L/powl(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L),2.0L)*2.0L)+(MY_PI_L*a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))+(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))))))/L+(MY_PI_L*a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))+(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))))))/L-(MY_PI_L*a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L))+(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))*(powl(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L),2.0L)*(1.0L/2.0L)+powl(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L),2.0L)*(1.0L/2.0L)+powl(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L),2.0L)*(1.0L/2.0L)+(Upsilon*(p_0+p_x*cosl((CoordX*MY_PI_L*a_p_x)/L)+p_z*cosl((CoordZ*MY_PI_L*a_p_z)/L)+p_y*sinl((CoordY*MY_PI_L*a_p_y)/L)))/((gammaNumber-1.0L)*(rho_0+rho_y*cosl((CoordY*MY_PI_L*a_rho_y)/L)+rho_x*sinl((CoordX*MY_PI_L*a_rho_x)/L)+rho_z*sinl((CoordZ*MY_PI_L*a_rho_z)/L))))))/L-1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_v_z*mue*v_z*cosl((CoordZ*MY_PI_L*a_v_z)/L)*(a_v_z*v_z*cosl((CoordZ*MY_PI_L*a_v_z)/L)+a_w_y*w_y*cosl((CoordY*MY_PI_L*a_w_y)/L))-1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_w_y*mue*w_y*cosl((CoordY*MY_PI_L*a_w_y)/L)*(a_v_z*v_z*cosl((CoordZ*MY_PI_L*a_v_z)/L)+a_w_y*w_y*cosl((CoordY*MY_PI_L*a_w_y)/L))-1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_w_x*mue*w_x*cosl((CoordX*MY_PI_L*a_w_x)/L)*(a_w_x*w_x*cosl((CoordX*MY_PI_L*a_w_x)/L)-a_u_z*u_z*sinl((CoordZ*MY_PI_L*a_u_z)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_u_z*mue*u_z*sinl((CoordZ*MY_PI_L*a_u_z)/L)*(a_w_x*w_x*cosl((CoordX*MY_PI_L*a_w_x)/L)-a_u_z*u_z*sinl((CoordZ*MY_PI_L*a_u_z)/L))-1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_u_y*mue*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)*(a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)+a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L))-1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_v_x*mue*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L)*(a_u_y*u_y*sinl((CoordY*MY_PI_L*a_u_y)/L)+a_v_x*v_x*sinl((CoordX*MY_PI_L*a_v_x)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_u_y*a_u_y)*mue*u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_u_z*a_u_z)*mue*u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_u_x*a_u_x)*mue*u_x*sinl((CoordX*MY_PI_L*a_u_x)/L)*(u_0+u_y*cosl((CoordY*MY_PI_L*a_u_y)/L)+u_z*cosl((CoordZ*MY_PI_L*a_u_z)/L)+u_x*sinl((CoordX*MY_PI_L*a_u_x)/L))*(4.0L/3.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_v_x*a_v_x)*mue*v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_w_z*a_w_z)*mue*w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L))*(4.0L/3.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_v_y*a_v_y)*mue*v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L))*(4.0L/3.0L)+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_v_z*a_v_z)*mue*v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L)*(v_0+v_x*cosl((CoordX*MY_PI_L*a_v_x)/L)+v_y*sinl((CoordY*MY_PI_L*a_v_y)/L)+v_z*sinl((CoordZ*MY_PI_L*a_v_z)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_w_x*a_w_x)*mue*w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L))+1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*(a_w_y*a_w_y)*mue*w_y*sinl((CoordY*MY_PI_L*a_w_y)/L)*(w_0+w_z*cosl((CoordZ*MY_PI_L*a_w_z)/L)+w_x*sinl((CoordX*MY_PI_L*a_w_x)/L)+w_y*sinl((CoordY*MY_PI_L*a_w_y)/L))-1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_u_x*mue*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*(a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)*2.0L-a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)+a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L))*(2.0L/3.0L)-1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_v_y*mue*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*(-a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)+a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)*2.0L+a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L))*(2.0L/3.0L)-1.0L/(L*L)*(MY_PI_L*MY_PI_L)*Psi*a_w_z*mue*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*(a_u_x*u_x*cosl((CoordX*MY_PI_L*a_u_x)/L)+a_v_y*v_y*cosl((CoordY*MY_PI_L*a_v_y)/L)+a_w_z*w_z*sinl((CoordZ*MY_PI_L*a_w_z)/L)*2.0L)*(2.0L/3.0L);
				}


				pnt_Q->zetaMomentum[ijk]+=pnt_mesh->jacobian[ijk]*z_momentum;
#endif



				pnt_Q->Mass[ijk]+=pnt_mesh->jacobian[ijk]*mass;
				pnt_Q->xiMomentum[ijk]+=pnt_mesh->jacobian[ijk]*x_momentum;
				pnt_Q->etaMomentum[ijk]+=pnt_mesh->jacobian[ijk]*y_momentum;
				pnt_Q->Energy[ijk]+=pnt_mesh->jacobian[ijk]*energy;
			}
		}
	}
	/*if (pnt_config->MPI_rank==0){printf("SHOCK: Fehler durch PRECISION:\n "
			"old: %.128Le\n "
			"new: %.128Le\n "
			"delta: %.128Le\n",
			value_old,
			value_new,
			(value_new-value_old)
			);}*/
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
					pnt_mesh->BC_Corrector_xiMomentum[ijk]=1.0;
					pnt_mesh->BC_Corrector_etaMomentum[ijk]=1.0;
					pnt_mesh->BC_Corrector_zetaMomentum[ijk]=1.0;
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
			{strcpy( pnt_config->BC_Top,pnt_config->BCManufacturedSolution );}

		if (pnt_config->InterfaceNeighbourRight==NO_NEIGHBOUR)
			{strcpy( pnt_config->BC_Right,pnt_config->BCManufacturedSolution );}
#if MESHDIMENSIONS==3
		if (pnt_config->InterfaceNeighbourBehind==NO_NEIGHBOUR)
			{strcpy( pnt_config->BC_Behind,pnt_config->BCManufacturedSolution );}

		if (pnt_config->InterfaceNeighbourInFront==NO_NEIGHBOUR)
			{strcpy( pnt_config->BC_InFront,pnt_config->BCManufacturedSolution );}
#endif
	}

	long double rho_0,rho_x,rho_y,rho_z,a_rho_x,a_rho_y,a_rho_z;
	long double u_0,u_x,u_y,u_z,a_u_x,a_u_y,a_u_z;
	long double v_0,v_x,v_y,v_z,a_v_x,a_v_y,a_v_z;
	long double w_0,w_x,w_y,w_z,a_w_x,a_w_y,a_w_z;
	long double p_0,p_x,p_y,p_z,a_p_x,a_p_y,a_p_z;

	long double rho_ref;//=1.0L;
	long double u_ref;//=(long double)(pnt_config->machNumber*sqrtl(pnt_config->gammaNumber*pnt_config->gasConstantNumber*pnt_config->T0_dim));
	long double p_ref;//=(long double)rho_ref*u_ref*u_ref/(pnt_config->gammaNumber*pnt_config->machNumber*pnt_config->machNumber);

	p_ref=100000.0L;
	u_ref=(long double)(pnt_config->machNumber*sqrtl(pnt_config->gammaNumber*pnt_config->gasConstantNumber*pnt_config->T0_dim));
	rho_ref=1.0L;


	p_z=0.0L;	a_p_z=0.0L;
	u_z=0.0L;	a_u_z=0.0L;
	v_z=0.0L;	a_v_z=0.0L;
	rho_z=0.0L;	a_rho_z=0.0L;
	w_0=0.0L;	w_x=0.0L;	w_y=0.0L;	w_z=0.0L;	a_w_x=0.0L;	a_w_y=0.0L;	a_w_z=0.0L;

#if MESHDIMENSIONS==2
	if (pnt_config->ManufacturedSolution_case==1)
	{
		if(pnt_config->flag_Inviscid==1)
		{
			//supersonic-Euler
			rho_0=1.0L/rho_ref;		rho_x=0.15L/rho_ref;		rho_y=-0.1L/rho_ref;		a_rho_x=1.0L;	a_rho_y=0.5L;
			u_0=800.0L/u_ref;		u_x=50.0L/u_ref;			u_y=-30.0L/u_ref;			a_u_x=1.5L;		a_u_y=0.6L;
			v_0=800.0L/u_ref;		v_x=-75.0L/u_ref;			v_y=40.0L/u_ref;			a_v_x=0.5L;		a_v_y=2.L/3.L;
			p_0=100000.0L/p_ref;	p_x=0.2L*100000.0L/p_ref;	p_y=0.5L*100000.0L/p_ref;	a_p_x=2.0L;		a_p_y=1.0L;
		}
		else
		{
			//supersonic-NavierStokes
			rho_0=1.0L/rho_ref;		rho_x=0.15L/rho_ref;		rho_y=-0.1L/rho_ref;		a_rho_x=1.0L;	a_rho_y=0.5L;
			u_0=800.0L/u_ref;		u_x=50.0L/u_ref;			u_y=-30.0L/u_ref;			a_u_x=1.5L;		a_u_y=0.6L;
			v_0=800.0L/u_ref;		v_x=-75.0L/u_ref;			v_y=40.0L/u_ref;			a_v_x=0.5L;		a_v_y=1.5L;
			p_0=100000.0L/p_ref;	p_x=0.2L*100000.0L/p_ref;	p_y=0.5L*100000.0L/p_ref;	a_p_x=2.0L/3.0L;a_p_y=1.0L;
		}
	}
	else // (pnt_config->ManufacturedSolution_case==0)
	{
		if(pnt_config->flag_Inviscid==1)
		{
			//subsonic-Euler
			rho_0=1.0L/rho_ref;		rho_x=0.15L/rho_ref;		rho_y=-0.1L/rho_ref;		a_rho_x=1.0L;	a_rho_y=0.5L;
			u_0=70.0L/u_ref;		u_x=5.0L/u_ref;				u_y=-7.0L/u_ref;			a_u_x=1.5L;		a_u_y=0.6L;
			v_0=90.0L/u_ref;		v_x=-15.0L/u_ref;			v_y=8.5L/u_ref;				a_v_x=0.5L;		a_v_y=2.L/3.L;
			p_0=100000.0L/p_ref;	p_x=0.2L*100000.0L/p_ref;	p_y=0.5L*100000.0L/p_ref;	a_p_x=2.0L;		a_p_y=1.0L;
		}
		else
		{
			//subsonic-NavierStokes
			rho_0=1.0L/rho_ref;		rho_x=0.1L/rho_ref;			rho_y=0.15L/rho_ref;		a_rho_x=0.75L;	a_rho_y=1.0L;
			u_0=70.0L/u_ref;		u_x=4.0L/u_ref;				u_y=-12.0L/u_ref;			a_u_x=5.0L/3.0L;a_u_y=1.5L;
			v_0=90.0L/u_ref;		v_x=-20.0L/u_ref;			v_y=4.0L/u_ref;				a_v_x=1.5L;		a_v_y=1.0L;
			p_0=100000.0L/p_ref;	p_x=-0.3L*100000.0L/p_ref;	p_y=0.2L*100000.0L/p_ref;	a_p_x=1.0L;		a_p_y=1.25L;
		}
	}
#endif
#if MESHDIMENSIONS==3
	//supersonic-Euler
	rho_0=1.0L/rho_ref;		rho_x=0.15L/rho_ref;		rho_y=-0.1L/rho_ref;		rho_z=-0.12L/rho_ref;		a_rho_x=1.0L;	a_rho_y=0.5L;	a_rho_z=1.5L;
	u_0=800.0L/u_ref;		u_x=50.0L/u_ref;			u_y=-30.0L/u_ref;			u_z=-18.0L/u_ref;			a_u_x=1.5L;		a_u_y=0.6L;		a_u_z=0.5L;
	v_0=800.0L/u_ref;		v_x=-75.0L/u_ref;			v_y=40.0L/u_ref;			v_z=-30.0L/u_ref;			a_v_x=0.5L;		a_v_y=2.L/3.L;	a_v_z=1.25L;
	w_0=800.0L/u_ref;		w_x=15.0L/u_ref;			w_y=-25.0L/u_ref;			w_z=35.0L/u_ref;			a_w_x=1.0L/3.0L;a_w_y=1.5L;		a_w_z=1.0L;
	p_0=100000.0L/p_ref;	p_x=0.2L*100000.0L/p_ref;	p_y=0.5L*100000.0L/p_ref;	p_z=-0.35L*100000.0L/p_ref;	a_p_x=2.0L;		a_p_y=1.0L;		a_p_z=1.L/3.L;
#endif

	pnt_config->ManufacturedSolution_param_rho[0]=rho_0;
	pnt_config->ManufacturedSolution_param_rho[1]=rho_x;
	pnt_config->ManufacturedSolution_param_rho[2]=rho_y;
	pnt_config->ManufacturedSolution_param_rho[3]=rho_z;
	pnt_config->ManufacturedSolution_param_rho[4]=a_rho_x;
	pnt_config->ManufacturedSolution_param_rho[5]=a_rho_y;
	pnt_config->ManufacturedSolution_param_rho[6]=a_rho_z;

	pnt_config->ManufacturedSolution_param_u[0]=u_0;
	pnt_config->ManufacturedSolution_param_u[1]=u_x;
	pnt_config->ManufacturedSolution_param_u[2]=u_y;
	pnt_config->ManufacturedSolution_param_u[3]=u_z;
	pnt_config->ManufacturedSolution_param_u[4]=a_u_x;
	pnt_config->ManufacturedSolution_param_u[5]=a_u_y;
	pnt_config->ManufacturedSolution_param_u[6]=a_u_z;

	pnt_config->ManufacturedSolution_param_v[0]=v_0;
	pnt_config->ManufacturedSolution_param_v[1]=v_x;
	pnt_config->ManufacturedSolution_param_v[2]=v_y;
	pnt_config->ManufacturedSolution_param_v[3]=v_z;
	pnt_config->ManufacturedSolution_param_v[4]=a_v_x;
	pnt_config->ManufacturedSolution_param_v[5]=a_v_y;
	pnt_config->ManufacturedSolution_param_v[6]=a_v_z;

	pnt_config->ManufacturedSolution_param_w[0]=w_0;
	pnt_config->ManufacturedSolution_param_w[1]=w_x;
	pnt_config->ManufacturedSolution_param_w[2]=w_y;
	pnt_config->ManufacturedSolution_param_w[3]=w_z;
	pnt_config->ManufacturedSolution_param_w[4]=a_w_x;
	pnt_config->ManufacturedSolution_param_w[5]=a_w_y;
	pnt_config->ManufacturedSolution_param_w[6]=a_w_z;

	pnt_config->ManufacturedSolution_param_p[0]=p_0;
	pnt_config->ManufacturedSolution_param_p[1]=p_x;
	pnt_config->ManufacturedSolution_param_p[2]=p_y;
	pnt_config->ManufacturedSolution_param_p[3]=p_z;
	pnt_config->ManufacturedSolution_param_p[4]=a_p_x;
	pnt_config->ManufacturedSolution_param_p[5]=a_p_y;
	pnt_config->ManufacturedSolution_param_p[6]=a_p_z;

}


void WriteManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		int ijk)
{
	long double rho_0,rho_x,rho_y,rho_z,a_rho_x,a_rho_y,a_rho_z;
	long double u_0,u_x,u_y,u_z,a_u_x,a_u_y,a_u_z;
	long double v_0,v_x,v_y,v_z,a_v_x,a_v_y,a_v_z;
	long double w_0,w_x,w_y,w_z,a_w_x,a_w_y,a_w_z;
	long double p_0,p_x,p_y,p_z,a_p_x,a_p_y,a_p_z;

	rho_0=pnt_config->ManufacturedSolution_param_rho[0];
	rho_x=pnt_config->ManufacturedSolution_param_rho[1];
	rho_y=pnt_config->ManufacturedSolution_param_rho[2];
	rho_z=pnt_config->ManufacturedSolution_param_rho[3];
	a_rho_x=pnt_config->ManufacturedSolution_param_rho[4];
	a_rho_y=pnt_config->ManufacturedSolution_param_rho[5];
	a_rho_z=pnt_config->ManufacturedSolution_param_rho[6];

	u_0=pnt_config->ManufacturedSolution_param_u[0];
	u_x=pnt_config->ManufacturedSolution_param_u[1];
	u_y=pnt_config->ManufacturedSolution_param_u[2];
	u_z=pnt_config->ManufacturedSolution_param_u[3];
	a_u_x=pnt_config->ManufacturedSolution_param_u[4];
	a_u_y=pnt_config->ManufacturedSolution_param_u[5];
	a_u_z=pnt_config->ManufacturedSolution_param_u[6];

	v_0=pnt_config->ManufacturedSolution_param_v[0];
	v_x=pnt_config->ManufacturedSolution_param_v[1];
	v_y=pnt_config->ManufacturedSolution_param_v[2];
	v_z=pnt_config->ManufacturedSolution_param_v[3];
	a_v_x=pnt_config->ManufacturedSolution_param_v[4];
	a_v_y=pnt_config->ManufacturedSolution_param_v[5];
	a_v_z=pnt_config->ManufacturedSolution_param_v[6];

	w_0=pnt_config->ManufacturedSolution_param_w[0];
	w_x=pnt_config->ManufacturedSolution_param_w[1];
	w_y=pnt_config->ManufacturedSolution_param_w[2];
	w_z=pnt_config->ManufacturedSolution_param_w[3];
	a_w_x=pnt_config->ManufacturedSolution_param_w[4];
	a_w_y=pnt_config->ManufacturedSolution_param_w[5];
	a_w_z=pnt_config->ManufacturedSolution_param_w[6];

	p_0=pnt_config->ManufacturedSolution_param_p[0];
	p_x=pnt_config->ManufacturedSolution_param_p[1];
	p_y=pnt_config->ManufacturedSolution_param_p[2];
	p_z=pnt_config->ManufacturedSolution_param_p[3];
	a_p_x=pnt_config->ManufacturedSolution_param_p[4];
	a_p_y=pnt_config->ManufacturedSolution_param_p[5];
	a_p_z=pnt_config->ManufacturedSolution_param_p[6];

	long double CoordX,CoordY,CoordZ;

	CoordX=(long double)pnt_mesh->x[ijk];
	CoordY=(long double)pnt_mesh->y[ijk];
	CoordZ=(long double)pnt_mesh->z[ijk];
	long double L=1.0L;

	pnt_U_lastStep->rho[ijk]=(rho_0+rho_x*sinl(a_rho_x*MY_PI_L*CoordX/L)+rho_y*cosl(a_rho_y*MY_PI_L*CoordY/L)+rho_z*sinl(a_rho_z*MY_PI_L*CoordZ/L));
	pnt_U_lastStep->u[ijk]=(u_0+u_x*sinl(a_u_x*MY_PI_L*CoordX/L)+u_y*cosl(a_u_y*MY_PI_L*CoordY/L)+u_z*cosl(a_u_z*MY_PI_L*CoordZ/L));
	pnt_U_lastStep->v[ijk]=(v_0+v_x*cosl(a_v_x*MY_PI_L*CoordX/L)+v_y*sinl(a_v_y*MY_PI_L*CoordY/L)+v_z*sinl(a_v_z*MY_PI_L*CoordZ/L));
	pnt_U_lastStep->w[ijk]=(w_0+w_x*sinl(a_w_x*MY_PI_L*CoordX/L)+w_y*sinl(a_w_y*MY_PI_L*CoordY/L)+w_z*cosl(a_w_z*MY_PI_L*CoordZ/L));
	pnt_U_lastStep->p[ijk]=(p_0+p_x*cosl(a_p_x*MY_PI_L*CoordX/L)+p_y*sinl(a_p_y*MY_PI_L*CoordY/L)+p_z*cosl(a_p_z*MY_PI_L*CoordZ/L));

	pnt_U_lastStep->e[ijk]=(0.5*((pnt_U_lastStep->u[ijk]*pnt_U_lastStep->u[ijk])+(pnt_U_lastStep->v[ijk]*pnt_U_lastStep->v[ijk])+(pnt_U_lastStep->w[ijk]*pnt_U_lastStep->w[ijk]))+
							pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]/(pnt_config->gammaNumber-1.0)*pnt_config->Upsilon);
	pnt_U_lastStep->T[ijk]=
			fabsl(pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]);

	pnt_U_lastStep->c[ijk]=
			sqrt(pnt_config->Upsilon*
			pnt_config->gammaNumber*pnt_U_lastStep->p[ijk]/pnt_U_lastStep->rho[ijk]);

	pnt_U_lastStep->mue[ijk]=1.0/pnt_config->reynoldsNumber;


}

long double GetRhoManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		int ijk)
{
	long double L=1.0L;
	long double CoordX,CoordY,CoordZ;
	long double rho_0,rho_x,rho_y,rho_z,a_rho_x,a_rho_y,a_rho_z;

	rho_0=pnt_config->ManufacturedSolution_param_rho[0];
	rho_x=pnt_config->ManufacturedSolution_param_rho[1];
	rho_y=pnt_config->ManufacturedSolution_param_rho[2];
	rho_z=pnt_config->ManufacturedSolution_param_rho[3];
	a_rho_x=pnt_config->ManufacturedSolution_param_rho[4];
	a_rho_y=pnt_config->ManufacturedSolution_param_rho[5];
	a_rho_z=pnt_config->ManufacturedSolution_param_rho[6];

	CoordX=(long double)pnt_mesh->x[ijk];
	CoordY=(long double)pnt_mesh->y[ijk];
	CoordZ=(long double)pnt_mesh->z[ijk];

	long double rho;

	rho=(rho_0+rho_x*sinl(a_rho_x*MY_PI_L*CoordX/L)+rho_y*cosl(a_rho_y*MY_PI_L*CoordY/L)+rho_z*sinl(a_rho_z*MY_PI_L*CoordZ/L));

	return rho;
}

long double GetPressureManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		int ijk)
{
	long double L=1.0L;
	long double CoordX,CoordY,CoordZ;
	long double p_0,p_x,p_y,p_z,a_p_x,a_p_y,a_p_z;

	p_0=pnt_config->ManufacturedSolution_param_p[0];
	p_x=pnt_config->ManufacturedSolution_param_p[1];
	p_y=pnt_config->ManufacturedSolution_param_p[2];
	p_z=pnt_config->ManufacturedSolution_param_p[3];
	a_p_x=pnt_config->ManufacturedSolution_param_p[4];
	a_p_y=pnt_config->ManufacturedSolution_param_p[5];
	a_p_z=pnt_config->ManufacturedSolution_param_p[6];

	CoordX=(long double)pnt_mesh->x[ijk];
	CoordY=(long double)pnt_mesh->y[ijk];
	CoordZ=(long double)pnt_mesh->z[ijk];

	long double pressure;

	pressure=(p_0+p_x*cosl(a_p_x*MY_PI_L*CoordX/L)+p_y*sinl(a_p_y*MY_PI_L*CoordY/L)+p_z*cosl(a_p_z*MY_PI_L*CoordZ/L));

	return pressure;
}

void ErrorManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_Film * pnt_Film,
		struct strct_Flux * pnt_Q_sum,
		int flag)
{
	long double Linf_norm_rho=0.0L;
	long double Linf_norm_pressure=0.0L;
	long double L2_norm_rho=0.0L;
	long double L2_norm_pressure=0.0L;
	long double rho_exact,pressure_exact;
	long double N=pnt_config->int_iMeshPoints*pnt_config->int_jMeshPoints*pnt_config->int_kMeshPoints;
	int ijk1,ijk2;

	ijk1=1*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+0*pnt_config->int_kMeshPointsGhostCells+0;
	ijk2=2*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+0*pnt_config->int_kMeshPointsGhostCells+0;
	if(flag==1)
	{
		if(pnt_config->MPI_rank==0)
		{
			FILE * file0;
			char filename[200];
			sprintf(filename,"ManufacturedSolutions_W%d_%dP.dat",SPACEORDER,PRECISION);
			file0=fopen(filename,"a");
			fprintf(file0," %d %.20Le %Le %Le %Le %Le\n",
					PRECISION,
					(long double)(pnt_mesh->x[ijk2]-pnt_mesh->x[ijk1]),
					pnt_config->all_L2_norm_rho,
					pnt_config->all_L2_norm_pressure,
					pnt_config->all_Linf_norm_rho,
					pnt_config->all_Linf_norm_pressure);
			fclose(file0);
		}
	}

	if(flag==0)
	{
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
					L2_norm_rho+=powl((rho_exact-pnt_U_lastStep->rho[ijk]),2.0L)/(N);
					if(Linf_norm_rho<fabsl(rho_exact-pnt_U_lastStep->rho[ijk])){Linf_norm_rho=fabsl(rho_exact-pnt_U_lastStep->rho[ijk]);}

					pressure_exact=GetPressureManufacturedSolution(
											pnt_config,
											pnt_mesh,
											pnt_U_lastStep,
											ijk);
					L2_norm_pressure+=powl((pressure_exact-pnt_U_lastStep->p[ijk]),2.0L)/(N);
					if(Linf_norm_pressure<fabsl(pressure_exact-pnt_U_lastStep->p[ijk])){Linf_norm_pressure=fabsl(pressure_exact-pnt_U_lastStep->p[ijk]);}

				}
			}
		}

		long double* all_L2_norm_rho;
		all_L2_norm_rho = (long double *)malloc(1*sizeof(long double));
		MPI_Reduce( &L2_norm_rho, all_L2_norm_rho,1,MPI_LONG_DOUBLE,MPI_SUM,0,pnt_config->MPI_comm);
		all_L2_norm_rho[0]=sqrtl(all_L2_norm_rho[0]);


		long double* all_L2_norm_pressure;
		all_L2_norm_pressure = (long double *)malloc(1*sizeof(long double));
		MPI_Reduce( &L2_norm_pressure, all_L2_norm_pressure,1,MPI_LONG_DOUBLE,MPI_SUM,0,pnt_config->MPI_comm);
		all_L2_norm_pressure[0]=sqrtl(all_L2_norm_pressure[0]);


		long double* all_Linf_norm_rho;
		all_Linf_norm_rho = (long double *)malloc(1*sizeof(long double));
		MPI_Reduce( &Linf_norm_rho, all_Linf_norm_rho,1,MPI_LONG_DOUBLE,MPI_MAX,0,pnt_config->MPI_comm);


		long double* all_Linf_norm_pressure;
		all_Linf_norm_pressure = (long double *)malloc(1*sizeof(long double));
		MPI_Reduce( &Linf_norm_pressure, all_Linf_norm_pressure,1,MPI_LONG_DOUBLE,MPI_MAX,0,pnt_config->MPI_comm);


		//Monitoring of Convergence by Density
		long double delta_rho,delta_pressure;
		strcpy(pnt_config->ManufacturedSolution_L2_Delta_name,"Density");
		delta_rho=(long double)fabsl(pnt_config->ManufacturedSolution_L2_last-all_L2_norm_rho[0]);
		delta_pressure=(long double)fabsl(pnt_config->ManufacturedSolution_L2_last_pressure-all_L2_norm_pressure[0]);
		pnt_config->ManufacturedSolution_L2_Delta=fabsl(pnt_config->ManufacturedSolution_L2_last-all_L2_norm_rho[0]);
		pnt_config->ManufacturedSolution_L2_last=all_L2_norm_rho[0];
		pnt_config->ManufacturedSolution_L2_last_pressure=all_L2_norm_pressure[0];



		if(pnt_config->MPI_rank==0)
		{
			FILE * file0;
			char filename[200];
			sprintf(filename,"Residual_W%d_%dP.dat",SPACEORDER,PRECISION);

			if (pnt_config->int_actualIteration==1)
			{
				file0=fopen(filename,"w");
				fprintf(file0,"VARIABLES = \"Iteration\" \"L2(rho)\" \"Residual(rho)\" \"L2(pressure)\" \"Residual(pressure)\" \"Residual(Mass)\" \"Residual(xiMomentum)\" \"Residual(etaMomentum)\" \"Residual(Energy)\"\n");
				fprintf(file0,"TITLE=\"Convergence\"\n");
				fprintf(file0,"ZONE T=\"W%d-%d\", F=POINT, I=0, DT=(DOUBLE)\n",SPACEORDER,PRECISION);
			}
			else
			{
				file0=fopen(filename,"a");
			}
			fprintf(file0," %d %Le %Le %Le %Le %Le %Le %Le %Le\n",
					pnt_config->int_actualIteration,
					(long double)all_L2_norm_rho[0],
					delta_rho,
					(long double)all_L2_norm_pressure[0],
					delta_pressure,
					fabsl(pnt_config->ManufacturedSolution_last_Q_Mass-(long double)pnt_Q_sum->Mass[ijk]),
					fabsl(pnt_config->ManufacturedSolution_last_Q_xiMomentum-(long double)pnt_Q_sum->xiMomentum[ijk]),
					fabsl(pnt_config->ManufacturedSolution_last_Q_etaMomentum-(long double)pnt_Q_sum->etaMomentum[ijk]),
					fabsl(pnt_config->ManufacturedSolution_last_Q_Energy-(long double)pnt_Q_sum->Energy[ijk]));
			pnt_config->ManufacturedSolution_last_Q_Mass=pnt_Q_sum->Mass[ijk];
			pnt_config->ManufacturedSolution_last_Q_xiMomentum=pnt_Q_sum->xiMomentum[ijk];
			pnt_config->ManufacturedSolution_last_Q_etaMomentum=pnt_Q_sum->etaMomentum[ijk];
			pnt_config->ManufacturedSolution_last_Q_Energy=pnt_Q_sum->Energy[ijk];

			fclose(file0);
		}

		/*if(pnt_config->MPI_rank==0){printf("SHOCK: L2_norm(rho):%.10Le (delta:%.8Le)\n",
				all_L2_norm_rho[0],
				delta);}*/

		MPI_Bcast(&pnt_config->ManufacturedSolution_L2_Delta,1,MPI_LONG_DOUBLE,0,pnt_config->MPI_comm);

		if((pnt_config->ManufacturedSolution_L2_Delta<CONV_ERROR)&&(pnt_config->int_actualIteration>500))
		{

			MPI_Bcast(all_L2_norm_rho,1,MPI_LONG_DOUBLE,0,pnt_config->MPI_comm);
			MPI_Bcast(all_L2_norm_pressure,1,MPI_LONG_DOUBLE,0,pnt_config->MPI_comm);
			MPI_Bcast(all_Linf_norm_rho,1,MPI_LONG_DOUBLE,0,pnt_config->MPI_comm);
			MPI_Bcast(all_Linf_norm_pressure,1,MPI_LONG_DOUBLE,0,pnt_config->MPI_comm);

			for(i=0;i<14;i++)
			{
				pnt_config->ManufacturedSolution_L2_Converged[i]=pnt_config->ManufacturedSolution_L2_Converged[i+1];
			}
			pnt_config->ManufacturedSolution_L2_Converged[14]=pnt_config->int_actualIteration;

			//gemittelte Variante
			//pnt_config->all_L2_norm_rho+=all_L2_norm_rho[0]/15.L;
			//pnt_config->all_L2_norm_pressure+=all_L2_norm_pressure[0]/15.L;
			//pnt_config->all_Linf_norm_rho+=all_Linf_norm_rho[0]/15.L;
			//pnt_config->all_Linf_norm_pressure+=all_Linf_norm_pressure[0]/15.L;

			//Variante letzter Fehlerwert
			pnt_config->all_L2_norm_rho=all_L2_norm_rho[0];
			pnt_config->all_L2_norm_pressure=all_L2_norm_pressure[0];
			pnt_config->all_Linf_norm_rho=all_Linf_norm_rho[0];
			pnt_config->all_Linf_norm_pressure=all_Linf_norm_pressure[0];

			pnt_config->ManufacturedSolution_L2_counter++;

			if(pnt_config->int_actualIteration==pnt_config->ManufacturedSolution_L2_Converged[0]+14)
			{
				if(pnt_config->MPI_rank==0){printf("SHOCK: %d. time: Convergence limit (%.8Le) of L2_norm(last Delta:%.8Le) reached. Exit!\n",
						pnt_config->ManufacturedSolution_L2_counter,
						CONV_ERROR,
						pnt_config->ManufacturedSolution_L2_Delta);}
				pnt_config->int_actualIteration=pnt_config->int_EndIteration+100;

//				for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
//				{
//					for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
//					{
//						for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
//						{
//							ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
////							rho_exact=GetRhoManufacturedSolution(
////									pnt_config,
////									pnt_mesh,
////									pnt_U_lastStep,
////									ijk);
////							pnt_U_lastStep->rho[ijk]=fabs(pnt_U_lastStep->rho[ijk]-rho_exact);
//
//							pnt_U_lastStep->rho[ijk]=pnt_U_lastStep->T[ijk];
//						}
//					}
//				}



				pnt_config->int_actualSample=pnt_config->int_Samples-1;
				WriteValuesFromUToFilm(
					pnt_config,
					pnt_U_lastStep,
					pnt_Film,
					pnt_mesh);

			}
		}
	}


}
