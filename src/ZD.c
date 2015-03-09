#include "SHOCK.h"
#include "Functions.h"
#include "ZD.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void CalcViscidFluxesInXiDirectionDirectly(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q)
{
	int ijk,i,j,k,l,iMinus1jk,iPlus1jk;
	int ijkMAX;
	int indexMinus_xi, indexPlus_xi;

	FLT u_xi_MinusHalf, u_eta_MinusHalf, u_zeta_MinusHalf;
	FLT v_xi_MinusHalf, v_eta_MinusHalf, v_zeta_MinusHalf;
 	FLT w_xi_MinusHalf, w_eta_MinusHalf, w_zeta_MinusHalf;
 	FLT T_xi_MinusHalf, T_eta_MinusHalf, T_zeta_MinusHalf;
	FLT u_xi_PlusHalf, u_eta_PlusHalf, u_zeta_PlusHalf;
	FLT v_xi_PlusHalf, v_eta_PlusHalf, v_zeta_PlusHalf;
 	FLT w_xi_PlusHalf, w_eta_PlusHalf, w_zeta_PlusHalf;
 	FLT T_xi_PlusHalf, T_eta_PlusHalf, T_zeta_PlusHalf;

	FLT u_iMinusHalf, u_iPlusHalf;
	FLT v_iMinusHalf, v_iPlusHalf;
	FLT w_iMinusHalf, w_iPlusHalf;
	FLT T_iMinusHalf, T_iPlusHalf;

	FLT XiMomRotSymm,EtaMomRotSymm,EnergyRotSymm;
	XiMomRotSymm=0.0;
	EtaMomRotSymm=0.0;
	EnergyRotSymm=0.0;

	ijkMAX=pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;


	FLT Koeffizient_PlusHalf[30];
	FLT Koeffizient_MinusHalf[30];

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iPlus1jk=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				u_xi_MinusHalf=0.;		v_xi_MinusHalf=0.;			w_xi_MinusHalf=0.;		T_xi_MinusHalf=0.;
				u_eta_MinusHalf=0.;		v_eta_MinusHalf=0.;			w_eta_MinusHalf=0.;		T_eta_MinusHalf=0.;
				u_zeta_MinusHalf=0.;	v_zeta_MinusHalf=0.;		w_zeta_MinusHalf=0.;	T_zeta_MinusHalf=0.;
				u_xi_PlusHalf=0.;		v_xi_PlusHalf=0.;			w_xi_PlusHalf=0.;		T_xi_PlusHalf=0.;
				u_eta_PlusHalf=0.;		v_eta_PlusHalf=0.;			w_eta_PlusHalf=0.;		T_eta_PlusHalf=0.;
				u_zeta_PlusHalf=0.;		v_zeta_PlusHalf=0.;			w_zeta_PlusHalf=0.;		T_zeta_PlusHalf=0.;

				u_iMinusHalf=0.;		v_iMinusHalf=0.;			w_iMinusHalf=0.;		T_iMinusHalf=0.;
				u_iPlusHalf=0.;			v_iPlusHalf=0.;				w_iPlusHalf=0.;			T_iPlusHalf=0.;


				for(l=0;l<=SPACEORDER;l++)
				{
					indexMinus_xi=(i+l-(SPACEORDER+1)/2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					indexPlus_xi=(i+l-(SPACEORDER+1)/2+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

					u_iMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u[indexMinus_xi];		u_iPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u[indexPlus_xi];
					v_iMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v[indexMinus_xi];		v_iPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v[indexPlus_xi];
					T_iMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T[indexMinus_xi];		T_iPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T[indexPlus_xi];

					u_xi_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->u[indexMinus_xi];	u_xi_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->u[indexPlus_xi];
					v_xi_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->v[indexMinus_xi];	v_xi_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->v[indexPlus_xi];
					T_xi_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->T[indexMinus_xi];	T_xi_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->T[indexPlus_xi];

					u_eta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_eta[indexMinus_xi];	u_eta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_eta[indexPlus_xi];
					v_eta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_eta[indexMinus_xi];	v_eta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_eta[indexPlus_xi];
					T_eta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_eta[indexMinus_xi];	T_eta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_eta[indexPlus_xi];

#if MESHDIMENSIONS==3

						w_iMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w[indexMinus_xi];		w_iPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w[indexPlus_xi];

						w_xi_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->w[indexMinus_xi];	w_xi_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->w[indexPlus_xi];

						w_eta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_eta[indexMinus_xi];	w_eta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_eta[indexPlus_xi];

						u_zeta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_zeta[indexMinus_xi];	u_zeta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_zeta[indexPlus_xi];
						v_zeta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_zeta[indexMinus_xi];	v_zeta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_zeta[indexPlus_xi];
						w_zeta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_zeta[indexMinus_xi];	w_zeta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_zeta[indexPlus_xi];
						T_zeta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_zeta[indexMinus_xi];	T_zeta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_zeta[indexPlus_xi];
#endif
				}


				for(l=0;l<30;l++)
				{
					Koeffizient_MinusHalf[l]=
							0.5*
							(pnt_U_RK->mue[iMinus1jk]*pnt_mesh->xiFluss_Faktor[l*ijkMAX+iMinus1jk]+
							pnt_U_RK->mue[ijk]*pnt_mesh->xiFluss_Faktor[l*ijkMAX+ijk]);
					Koeffizient_PlusHalf[l]=
							0.5*
							(pnt_U_RK->mue[ijk]*pnt_mesh->xiFluss_Faktor[l*ijkMAX+ijk]+
							pnt_U_RK->mue[iPlus1jk]*pnt_mesh->xiFluss_Faktor[l*ijkMAX+iPlus1jk]);
				}


				if(pnt_config->flag_rotation_symmetric==1)
				{
					XiMomRotSymm=pnt_mesh->jacobian[ijk]*(1./3.*pnt_mesh->xi_x[ijk]*0.5*(v_xi_MinusHalf+v_xi_PlusHalf)
									+pnt_mesh->xi_y[ijk]*0.5*(u_xi_MinusHalf+u_xi_PlusHalf))
									*pnt_config->Psi*pnt_U_RK->mue[ijk]/pnt_mesh->y[ijk];

					EtaMomRotSymm=pnt_mesh->jacobian[ijk]*2./3.*(pnt_mesh->xi_y[ijk]*0.5*(v_xi_MinusHalf+v_xi_PlusHalf)
									+pnt_mesh->xi_x[ijk]*0.5*(u_xi_MinusHalf+u_xi_PlusHalf)+pnt_U_RK->v[ijk]/pnt_mesh->y[ijk])
									*pnt_config->Psi*pnt_U_RK->mue[ijk]/pnt_mesh->y[ijk];

					EnergyRotSymm=pnt_mesh->jacobian[ijk]*((1./3.*pnt_U_RK->u[ijk]*pnt_mesh->xi_x[ijk]+2./3.*pnt_U_RK->v[ijk]*pnt_mesh->xi_y[ijk])*0.5*(v_xi_MinusHalf+v_xi_PlusHalf)
									+(pnt_U_RK->u[ijk]*pnt_mesh->xi_y[ijk]+2./3.*pnt_U_RK->v[ijk]*pnt_mesh->xi_x[ijk])*0.5*(u_xi_MinusHalf+u_xi_PlusHalf)
									-2./3.*pnt_U_RK->v[ijk]*pnt_U_RK->v[ijk]/pnt_mesh->y[ijk])
									*pnt_config->Psi*pnt_U_RK->mue[ijk]/pnt_mesh->y[ijk]
									+0.5*(pnt_config->Gamma[iPlus1jk]+pnt_config->Gamma[ijk])*pnt_mesh->xi_y[ijk]*0.5*(T_xi_MinusHalf+T_xi_PlusHalf)/pnt_mesh->y[ijk];
				}

				pnt_Q->xiMomentum[ijk]=pnt_Q->xiMomentum[ijk]+
						pnt_config->Psi*(
					u_xi_PlusHalf*Koeffizient_PlusHalf[0]-u_xi_MinusHalf*Koeffizient_MinusHalf[0]+
					u_eta_PlusHalf*Koeffizient_PlusHalf[1]-u_eta_MinusHalf*Koeffizient_MinusHalf[1]+
					u_zeta_PlusHalf*Koeffizient_PlusHalf[2]-u_zeta_MinusHalf*Koeffizient_MinusHalf[2]+
					v_xi_PlusHalf*Koeffizient_PlusHalf[3]-v_xi_MinusHalf*Koeffizient_MinusHalf[3]+
					v_eta_PlusHalf*Koeffizient_PlusHalf[4]-v_eta_MinusHalf*Koeffizient_MinusHalf[4]+
					v_zeta_PlusHalf*Koeffizient_PlusHalf[5]-v_zeta_MinusHalf*Koeffizient_MinusHalf[5]+
					w_xi_PlusHalf*Koeffizient_PlusHalf[6]-w_xi_MinusHalf*Koeffizient_MinusHalf[6]+
					w_eta_PlusHalf*Koeffizient_PlusHalf[7]-w_eta_MinusHalf*Koeffizient_MinusHalf[7]+
					w_zeta_PlusHalf*Koeffizient_PlusHalf[8]-w_zeta_MinusHalf*Koeffizient_MinusHalf[8])+

					pnt_config->flag_rotation_symmetric*XiMomRotSymm;

				pnt_Q->etaMomentum[ijk]=pnt_Q->etaMomentum[ijk]+
						pnt_config->Psi*(
					u_xi_PlusHalf*Koeffizient_PlusHalf[9]-u_xi_MinusHalf*Koeffizient_MinusHalf[9]+
					u_eta_PlusHalf*Koeffizient_PlusHalf[10]-u_eta_MinusHalf*Koeffizient_MinusHalf[10]+
					u_zeta_PlusHalf*Koeffizient_PlusHalf[11]-u_zeta_MinusHalf*Koeffizient_MinusHalf[11]+
					v_xi_PlusHalf*Koeffizient_PlusHalf[12]-v_xi_MinusHalf*Koeffizient_MinusHalf[12]+
					v_eta_PlusHalf*Koeffizient_PlusHalf[13]-v_eta_MinusHalf*Koeffizient_MinusHalf[13]+
					v_zeta_PlusHalf*Koeffizient_PlusHalf[14]-v_zeta_MinusHalf*Koeffizient_MinusHalf[14]+
					w_xi_PlusHalf*Koeffizient_PlusHalf[15]-w_xi_MinusHalf*Koeffizient_MinusHalf[15]+
					w_eta_PlusHalf*Koeffizient_PlusHalf[16]-w_eta_MinusHalf*Koeffizient_MinusHalf[16]+
					w_zeta_PlusHalf*Koeffizient_PlusHalf[17]-w_zeta_MinusHalf*Koeffizient_MinusHalf[17])+

					pnt_config->flag_rotation_symmetric*EtaMomRotSymm;

#if MESHDIMENSIONS==3
					pnt_Q->zetaMomentum[ijk]=pnt_Q->zetaMomentum[ijk]+
							pnt_config->Psi*(
						u_xi_PlusHalf*Koeffizient_PlusHalf[18]-u_xi_MinusHalf*Koeffizient_MinusHalf[18]+
						u_eta_PlusHalf*Koeffizient_PlusHalf[19]-u_eta_MinusHalf*Koeffizient_MinusHalf[19]+
						u_zeta_PlusHalf*Koeffizient_PlusHalf[20]-u_zeta_MinusHalf*Koeffizient_MinusHalf[20]+
						v_xi_PlusHalf*Koeffizient_PlusHalf[21]-v_xi_MinusHalf*Koeffizient_MinusHalf[21]+
						v_eta_PlusHalf*Koeffizient_PlusHalf[22]-v_eta_MinusHalf*Koeffizient_MinusHalf[22]+
						v_zeta_PlusHalf*Koeffizient_PlusHalf[23]-v_zeta_MinusHalf*Koeffizient_MinusHalf[23]+
						w_xi_PlusHalf*Koeffizient_PlusHalf[24]-w_xi_MinusHalf*Koeffizient_MinusHalf[24]+
						w_eta_PlusHalf*Koeffizient_PlusHalf[25]-w_eta_MinusHalf*Koeffizient_MinusHalf[25]+
						w_zeta_PlusHalf*Koeffizient_PlusHalf[26]-w_zeta_MinusHalf*Koeffizient_MinusHalf[26]);
#endif

				pnt_Q->Energy[ijk]=pnt_Q->Energy[ijk]+
						0.5*(pnt_config->Gamma[iPlus1jk]+pnt_config->Gamma[ijk])*(
					T_xi_PlusHalf*Koeffizient_PlusHalf[27]-T_xi_MinusHalf*Koeffizient_MinusHalf[27]+
					T_eta_PlusHalf*Koeffizient_PlusHalf[28]-T_eta_MinusHalf*Koeffizient_MinusHalf[28]+
					T_zeta_PlusHalf*Koeffizient_PlusHalf[29]-T_zeta_MinusHalf*Koeffizient_MinusHalf[29])+
						pnt_config->Psi*(
					u_iPlusHalf*u_xi_PlusHalf*Koeffizient_PlusHalf[0]-u_iMinusHalf*u_xi_MinusHalf*Koeffizient_MinusHalf[0]+
					u_iPlusHalf*u_eta_PlusHalf*Koeffizient_PlusHalf[1]-u_iMinusHalf*u_eta_MinusHalf*Koeffizient_MinusHalf[1]+
					u_iPlusHalf*u_zeta_PlusHalf*Koeffizient_PlusHalf[2]-u_iMinusHalf*u_zeta_MinusHalf*Koeffizient_MinusHalf[2]+
					u_iPlusHalf*v_xi_PlusHalf*Koeffizient_PlusHalf[3]-u_iMinusHalf*v_xi_MinusHalf*Koeffizient_MinusHalf[3]+
					u_iPlusHalf*v_eta_PlusHalf*Koeffizient_PlusHalf[4]-u_iMinusHalf*v_eta_MinusHalf*Koeffizient_MinusHalf[4]+
					u_iPlusHalf*v_zeta_PlusHalf*Koeffizient_PlusHalf[5]-u_iMinusHalf*v_zeta_MinusHalf*Koeffizient_MinusHalf[5]+
					u_iPlusHalf*w_xi_PlusHalf*Koeffizient_PlusHalf[6]-u_iMinusHalf*w_xi_MinusHalf*Koeffizient_MinusHalf[6]+
					u_iPlusHalf*w_eta_PlusHalf*Koeffizient_PlusHalf[7]-u_iMinusHalf*w_eta_MinusHalf*Koeffizient_MinusHalf[7]+
					u_iPlusHalf*w_zeta_PlusHalf*Koeffizient_PlusHalf[8]-u_iMinusHalf*w_zeta_MinusHalf*Koeffizient_MinusHalf[8]+

					v_iPlusHalf*u_xi_PlusHalf*Koeffizient_PlusHalf[9]-v_iMinusHalf*u_xi_MinusHalf*Koeffizient_MinusHalf[9]+
					v_iPlusHalf*u_eta_PlusHalf*Koeffizient_PlusHalf[10]-v_iMinusHalf*u_eta_MinusHalf*Koeffizient_MinusHalf[10]+
					v_iPlusHalf*u_zeta_PlusHalf*Koeffizient_PlusHalf[11]-v_iMinusHalf*u_zeta_MinusHalf*Koeffizient_MinusHalf[11]+
					v_iPlusHalf*v_xi_PlusHalf*Koeffizient_PlusHalf[12]-v_iMinusHalf*v_xi_MinusHalf*Koeffizient_MinusHalf[12]+
					v_iPlusHalf*v_eta_PlusHalf*Koeffizient_PlusHalf[13]-v_iMinusHalf*v_eta_MinusHalf*Koeffizient_MinusHalf[13]+
					v_iPlusHalf*v_zeta_PlusHalf*Koeffizient_PlusHalf[14]-v_iMinusHalf*v_zeta_MinusHalf*Koeffizient_MinusHalf[14]+
					v_iPlusHalf*w_xi_PlusHalf*Koeffizient_PlusHalf[15]-v_iMinusHalf*w_xi_MinusHalf*Koeffizient_MinusHalf[15]+
					v_iPlusHalf*w_eta_PlusHalf*Koeffizient_PlusHalf[16]-v_iMinusHalf*w_eta_MinusHalf*Koeffizient_MinusHalf[16]+
					v_iPlusHalf*w_zeta_PlusHalf*Koeffizient_PlusHalf[17]-v_iMinusHalf*w_zeta_MinusHalf*Koeffizient_MinusHalf[17]+

					w_iPlusHalf*u_xi_PlusHalf*Koeffizient_PlusHalf[18]-w_iMinusHalf*u_xi_MinusHalf*Koeffizient_MinusHalf[18]+
					w_iPlusHalf*u_eta_PlusHalf*Koeffizient_PlusHalf[19]-w_iMinusHalf*u_eta_MinusHalf*Koeffizient_MinusHalf[19]+
					w_iPlusHalf*u_zeta_PlusHalf*Koeffizient_PlusHalf[20]-w_iMinusHalf*u_zeta_MinusHalf*Koeffizient_MinusHalf[20]+
					w_iPlusHalf*v_xi_PlusHalf*Koeffizient_PlusHalf[21]-w_iMinusHalf*v_xi_MinusHalf*Koeffizient_MinusHalf[21]+
					w_iPlusHalf*v_eta_PlusHalf*Koeffizient_PlusHalf[22]-w_iMinusHalf*v_eta_MinusHalf*Koeffizient_MinusHalf[22]+
					w_iPlusHalf*v_zeta_PlusHalf*Koeffizient_PlusHalf[23]-w_iMinusHalf*v_zeta_MinusHalf*Koeffizient_MinusHalf[23]+
					w_iPlusHalf*w_xi_PlusHalf*Koeffizient_PlusHalf[24]-w_iMinusHalf*w_xi_MinusHalf*Koeffizient_MinusHalf[24]+
					w_iPlusHalf*w_eta_PlusHalf*Koeffizient_PlusHalf[25]-w_iMinusHalf*w_eta_MinusHalf*Koeffizient_MinusHalf[25]+
					w_iPlusHalf*w_zeta_PlusHalf*Koeffizient_PlusHalf[26]-w_iMinusHalf*w_zeta_MinusHalf*Koeffizient_MinusHalf[26])+

					pnt_config->flag_rotation_symmetric*EnergyRotSymm;
			}
		}
	}
}



void CalcViscidFluxesInEtaDirectionDirectly(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q)
{
	int ijk,i,j,k,l,ijMinus1k,ijPlus1k;
	int ijkMAX;
	int indexMinus_eta, indexPlus_eta;
	FLT u_xi_MinusHalf, u_eta_MinusHalf, u_zeta_MinusHalf;
	FLT v_xi_MinusHalf, v_eta_MinusHalf, v_zeta_MinusHalf;
 	FLT w_xi_MinusHalf, w_eta_MinusHalf, w_zeta_MinusHalf;
 	FLT T_xi_MinusHalf, T_eta_MinusHalf, T_zeta_MinusHalf;
	FLT u_xi_PlusHalf, u_eta_PlusHalf, u_zeta_PlusHalf;
	FLT v_xi_PlusHalf, v_eta_PlusHalf, v_zeta_PlusHalf;
 	FLT w_xi_PlusHalf, w_eta_PlusHalf, w_zeta_PlusHalf;
 	FLT T_xi_PlusHalf, T_eta_PlusHalf, T_zeta_PlusHalf;

	FLT u_jMinusHalf, u_jPlusHalf;
	FLT v_jMinusHalf, v_jPlusHalf;
	FLT w_jMinusHalf, w_jPlusHalf;
	FLT T_jMinusHalf, T_jPlusHalf;

	FLT XiMomRotSymm,EtaMomRotSymm,EnergyRotSymm;
	XiMomRotSymm=0.0;
	EtaMomRotSymm=0.0;
	EnergyRotSymm=0.0;

	ijkMAX=pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;


	FLT Koeffizient_PlusHalf[30];
	FLT Koeffizient_MinusHalf[30];

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;
				ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;

				u_xi_MinusHalf=0.;		v_xi_MinusHalf=0.;			w_xi_MinusHalf=0.;		T_xi_MinusHalf=0.;
				u_eta_MinusHalf=0.;		v_eta_MinusHalf=0.;			w_eta_MinusHalf=0.;		T_eta_MinusHalf=0.;
				u_zeta_MinusHalf=0.;	v_zeta_MinusHalf=0.;		w_zeta_MinusHalf=0.;	T_zeta_MinusHalf=0.;
				u_xi_PlusHalf=0.;		v_xi_PlusHalf=0.;			w_xi_PlusHalf=0.;		T_xi_PlusHalf=0.;
				u_eta_PlusHalf=0.;		v_eta_PlusHalf=0.;			w_eta_PlusHalf=0.;		T_eta_PlusHalf=0.;
				u_zeta_PlusHalf=0.;		v_zeta_PlusHalf=0.;			w_zeta_PlusHalf=0.;		T_zeta_PlusHalf=0.;

				u_jMinusHalf=0.;		v_jMinusHalf=0.;			w_jMinusHalf=0.;		T_jMinusHalf=0.;
				u_jPlusHalf=0.;			v_jPlusHalf=0.;				w_jPlusHalf=0.;			T_jPlusHalf=0.;

				for(l=0;l<=SPACEORDER;l++)
				{
					indexMinus_eta=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+l-(SPACEORDER+1)/2)*pnt_config->int_kMeshPointsGhostCells+k;
					indexPlus_eta=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+l-(SPACEORDER+1)/2+1)*pnt_config->int_kMeshPointsGhostCells+k;

					u_jMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u[indexMinus_eta];		u_jPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u[indexPlus_eta];
					v_jMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v[indexMinus_eta];		v_jPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v[indexPlus_eta];
					T_jMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T[indexMinus_eta];		T_jPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T[indexPlus_eta];

					u_eta_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->u[indexMinus_eta];	u_eta_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->u[indexPlus_eta];
					v_eta_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->v[indexMinus_eta];	v_eta_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->v[indexPlus_eta];
					T_eta_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->T[indexMinus_eta];	T_eta_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->T[indexPlus_eta];

					u_xi_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_xi[indexMinus_eta];	u_xi_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_xi[indexPlus_eta];
					v_xi_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_xi[indexMinus_eta];	v_xi_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_xi[indexPlus_eta];
					T_xi_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_xi[indexMinus_eta];	T_xi_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_xi[indexPlus_eta];

#if MESHDIMENSIONS==3
						w_jMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w[indexMinus_eta];		w_jPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w[indexPlus_eta];

						w_eta_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->w[indexMinus_eta];	w_eta_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->w[indexPlus_eta];

						w_xi_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_xi[indexMinus_eta];	w_xi_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_xi[indexPlus_eta];

						u_zeta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_zeta[indexMinus_eta];	u_zeta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_zeta[indexPlus_eta];
						v_zeta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_zeta[indexMinus_eta];	v_zeta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_zeta[indexPlus_eta];
						w_zeta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_zeta[indexMinus_eta];	w_zeta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_zeta[indexPlus_eta];
						T_zeta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_zeta[indexMinus_eta];	T_zeta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_zeta[indexPlus_eta];
#endif
				}


				for(l=0;l<30;l++)
				{
					Koeffizient_MinusHalf[l]=
							0.5*
							(pnt_U_RK->mue[ijMinus1k]*pnt_mesh->etaFluss_Faktor[l*ijkMAX+ijMinus1k]+
							pnt_U_RK->mue[ijk]*pnt_mesh->etaFluss_Faktor[l*ijkMAX+ijk]);
					Koeffizient_PlusHalf[l]=
							0.5*
							(pnt_U_RK->mue[ijk]*pnt_mesh->etaFluss_Faktor[l*ijkMAX+ijk]+
							pnt_U_RK->mue[ijPlus1k]*pnt_mesh->etaFluss_Faktor[l*ijkMAX+ijPlus1k]);
				}


				if(pnt_config->flag_rotation_symmetric==1)
				{
					XiMomRotSymm=pnt_mesh->jacobian[ijk]*(1./3.*pnt_mesh->eta_x[ijk]*0.5*(v_eta_MinusHalf+v_eta_PlusHalf)
							-pnt_mesh->xi_y[ijk]*0.5*(u_eta_MinusHalf+u_eta_PlusHalf))*pnt_config->Psi*pnt_U_RK->mue[ijk]/pnt_mesh->y[ijk];

					EtaMomRotSymm=pnt_mesh->jacobian[ijk]*2./3.*(pnt_mesh->xi_y[ijk]*0.5*(v_eta_MinusHalf+v_eta_PlusHalf)
							-pnt_mesh->eta_x[ijk]*0.5*(u_eta_MinusHalf+u_eta_PlusHalf))*pnt_config->Psi*pnt_U_RK->mue[ijk]/pnt_mesh->y[ijk];

					EnergyRotSymm=pnt_mesh->jacobian[ijk]*((1./3.*pnt_U_RK->u[ijk]*pnt_mesh->eta_x[ijk]+2./3.*pnt_U_RK->v[ijk]*pnt_mesh->xi_y[ijk])*0.5*(v_eta_MinusHalf+v_eta_PlusHalf)
							+(pnt_U_RK->u[ijk]*pnt_mesh->xi_y[ijk]-2./3.*pnt_U_RK->v[ijk]*pnt_mesh->eta_x[ijk])*0.5*(u_eta_MinusHalf+u_eta_PlusHalf))
							*pnt_config->Psi*pnt_U_RK->mue[ijk]/pnt_mesh->y[ijk]
							-0.5*(pnt_config->Gamma[ijPlus1k]+pnt_config->Gamma[ijk])*pnt_mesh->xi_y[ijk]*0.5*(T_eta_MinusHalf+T_eta_PlusHalf)/pnt_mesh->y[ijk];
				}


				pnt_Q->xiMomentum[ijk]=pnt_Q->xiMomentum[ijk]+
						pnt_config->Psi*(
					u_xi_PlusHalf*Koeffizient_PlusHalf[0]-u_xi_MinusHalf*Koeffizient_MinusHalf[0]+
					u_eta_PlusHalf*Koeffizient_PlusHalf[1]-u_eta_MinusHalf*Koeffizient_MinusHalf[1]+
					u_zeta_PlusHalf*Koeffizient_PlusHalf[2]-u_zeta_MinusHalf*Koeffizient_MinusHalf[2]+
					v_xi_PlusHalf*Koeffizient_PlusHalf[3]-v_xi_MinusHalf*Koeffizient_MinusHalf[3]+
					v_eta_PlusHalf*Koeffizient_PlusHalf[4]-v_eta_MinusHalf*Koeffizient_MinusHalf[4]+
					v_zeta_PlusHalf*Koeffizient_PlusHalf[5]-v_zeta_MinusHalf*Koeffizient_MinusHalf[5]+
					w_xi_PlusHalf*Koeffizient_PlusHalf[6]-w_xi_MinusHalf*Koeffizient_MinusHalf[6]+
					w_eta_PlusHalf*Koeffizient_PlusHalf[7]-w_eta_MinusHalf*Koeffizient_MinusHalf[7]+
					w_zeta_PlusHalf*Koeffizient_PlusHalf[8]-w_zeta_MinusHalf*Koeffizient_MinusHalf[8])+

					pnt_config->flag_rotation_symmetric*XiMomRotSymm;


				pnt_Q->etaMomentum[ijk]=pnt_Q->etaMomentum[ijk]+
						pnt_config->Psi*(
					u_xi_PlusHalf*Koeffizient_PlusHalf[9]-u_xi_MinusHalf*Koeffizient_MinusHalf[9]+
					u_eta_PlusHalf*Koeffizient_PlusHalf[10]-u_eta_MinusHalf*Koeffizient_MinusHalf[10]+
					u_zeta_PlusHalf*Koeffizient_PlusHalf[11]-u_zeta_MinusHalf*Koeffizient_MinusHalf[11]+
					v_xi_PlusHalf*Koeffizient_PlusHalf[12]-v_xi_MinusHalf*Koeffizient_MinusHalf[12]+
					v_eta_PlusHalf*Koeffizient_PlusHalf[13]-v_eta_MinusHalf*Koeffizient_MinusHalf[13]+
					v_zeta_PlusHalf*Koeffizient_PlusHalf[14]-v_zeta_MinusHalf*Koeffizient_MinusHalf[14]+
					w_xi_PlusHalf*Koeffizient_PlusHalf[15]-w_xi_MinusHalf*Koeffizient_MinusHalf[15]+
					w_eta_PlusHalf*Koeffizient_PlusHalf[16]-w_eta_MinusHalf*Koeffizient_MinusHalf[16]+
					w_zeta_PlusHalf*Koeffizient_PlusHalf[17]-w_zeta_MinusHalf*Koeffizient_MinusHalf[17])+

					pnt_config->flag_rotation_symmetric*EtaMomRotSymm;


#if MESHDIMENSIONS==3
					pnt_Q->zetaMomentum[ijk]=pnt_Q->zetaMomentum[ijk]+
							pnt_config->Psi*(
						u_xi_PlusHalf*Koeffizient_PlusHalf[18]-u_xi_MinusHalf*Koeffizient_MinusHalf[18]+
						u_eta_PlusHalf*Koeffizient_PlusHalf[19]-u_eta_MinusHalf*Koeffizient_MinusHalf[19]+
						u_zeta_PlusHalf*Koeffizient_PlusHalf[20]-u_zeta_MinusHalf*Koeffizient_MinusHalf[20]+
						v_xi_PlusHalf*Koeffizient_PlusHalf[21]-v_xi_MinusHalf*Koeffizient_MinusHalf[21]+
						v_eta_PlusHalf*Koeffizient_PlusHalf[22]-v_eta_MinusHalf*Koeffizient_MinusHalf[22]+
						v_zeta_PlusHalf*Koeffizient_PlusHalf[23]-v_zeta_MinusHalf*Koeffizient_MinusHalf[23]+
						w_xi_PlusHalf*Koeffizient_PlusHalf[24]-w_xi_MinusHalf*Koeffizient_MinusHalf[24]+
						w_eta_PlusHalf*Koeffizient_PlusHalf[25]-w_eta_MinusHalf*Koeffizient_MinusHalf[25]+
						w_zeta_PlusHalf*Koeffizient_PlusHalf[26]-w_zeta_MinusHalf*Koeffizient_MinusHalf[26]);
#endif

				pnt_Q->Energy[ijk]=pnt_Q->Energy[ijk]+
						0.5*(pnt_config->Gamma[ijPlus1k]+pnt_config->Gamma[ijk])*(
					T_xi_PlusHalf*Koeffizient_PlusHalf[27]-T_xi_MinusHalf*Koeffizient_MinusHalf[27]+
					T_eta_PlusHalf*Koeffizient_PlusHalf[28]-T_eta_MinusHalf*Koeffizient_MinusHalf[28]+
					T_zeta_PlusHalf*Koeffizient_PlusHalf[29]-T_zeta_MinusHalf*Koeffizient_MinusHalf[29])+
						pnt_config->Psi*(
					u_jPlusHalf*u_xi_PlusHalf*Koeffizient_PlusHalf[0]-u_jMinusHalf*u_xi_MinusHalf*Koeffizient_MinusHalf[0]+
					u_jPlusHalf*u_eta_PlusHalf*Koeffizient_PlusHalf[1]-u_jMinusHalf*u_eta_MinusHalf*Koeffizient_MinusHalf[1]+
					u_jPlusHalf*u_zeta_PlusHalf*Koeffizient_PlusHalf[2]-u_jMinusHalf*u_zeta_MinusHalf*Koeffizient_MinusHalf[2]+
					u_jPlusHalf*v_xi_PlusHalf*Koeffizient_PlusHalf[3]-u_jMinusHalf*v_xi_MinusHalf*Koeffizient_MinusHalf[3]+
					u_jPlusHalf*v_eta_PlusHalf*Koeffizient_PlusHalf[4]-u_jMinusHalf*v_eta_MinusHalf*Koeffizient_MinusHalf[4]+
					u_jPlusHalf*v_zeta_PlusHalf*Koeffizient_PlusHalf[5]-u_jMinusHalf*v_zeta_MinusHalf*Koeffizient_MinusHalf[5]+
					u_jPlusHalf*w_xi_PlusHalf*Koeffizient_PlusHalf[6]-u_jMinusHalf*w_xi_MinusHalf*Koeffizient_MinusHalf[6]+
					u_jPlusHalf*w_eta_PlusHalf*Koeffizient_PlusHalf[7]-u_jMinusHalf*w_eta_MinusHalf*Koeffizient_MinusHalf[7]+
					u_jPlusHalf*w_zeta_PlusHalf*Koeffizient_PlusHalf[8]-u_jMinusHalf*w_zeta_MinusHalf*Koeffizient_MinusHalf[8]+

					v_jPlusHalf*u_xi_PlusHalf*Koeffizient_PlusHalf[9]-v_jMinusHalf*u_xi_MinusHalf*Koeffizient_MinusHalf[9]+
					v_jPlusHalf*u_eta_PlusHalf*Koeffizient_PlusHalf[10]-v_jMinusHalf*u_eta_MinusHalf*Koeffizient_MinusHalf[10]+
					v_jPlusHalf*u_zeta_PlusHalf*Koeffizient_PlusHalf[11]-v_jMinusHalf*u_zeta_MinusHalf*Koeffizient_MinusHalf[11]+
					v_jPlusHalf*v_xi_PlusHalf*Koeffizient_PlusHalf[12]-v_jMinusHalf*v_xi_MinusHalf*Koeffizient_MinusHalf[12]+
					v_jPlusHalf*v_eta_PlusHalf*Koeffizient_PlusHalf[13]-v_jMinusHalf*v_eta_MinusHalf*Koeffizient_MinusHalf[13]+
					v_jPlusHalf*v_zeta_PlusHalf*Koeffizient_PlusHalf[14]-v_jMinusHalf*v_zeta_MinusHalf*Koeffizient_MinusHalf[14]+
					v_jPlusHalf*w_xi_PlusHalf*Koeffizient_PlusHalf[15]-v_jMinusHalf*w_xi_MinusHalf*Koeffizient_MinusHalf[15]+
					v_jPlusHalf*w_eta_PlusHalf*Koeffizient_PlusHalf[16]-v_jMinusHalf*w_eta_MinusHalf*Koeffizient_MinusHalf[16]+
					v_jPlusHalf*w_zeta_PlusHalf*Koeffizient_PlusHalf[17]-v_jMinusHalf*w_zeta_MinusHalf*Koeffizient_MinusHalf[17]+

					w_jPlusHalf*u_xi_PlusHalf*Koeffizient_PlusHalf[18]-w_jMinusHalf*u_xi_MinusHalf*Koeffizient_MinusHalf[18]+
					w_jPlusHalf*u_eta_PlusHalf*Koeffizient_PlusHalf[19]-w_jMinusHalf*u_eta_MinusHalf*Koeffizient_MinusHalf[19]+
					w_jPlusHalf*u_zeta_PlusHalf*Koeffizient_PlusHalf[20]-w_jMinusHalf*u_zeta_MinusHalf*Koeffizient_MinusHalf[20]+
					w_jPlusHalf*v_xi_PlusHalf*Koeffizient_PlusHalf[21]-w_jMinusHalf*v_xi_MinusHalf*Koeffizient_MinusHalf[21]+
					w_jPlusHalf*v_eta_PlusHalf*Koeffizient_PlusHalf[22]-w_jMinusHalf*v_eta_MinusHalf*Koeffizient_MinusHalf[22]+
					w_jPlusHalf*v_zeta_PlusHalf*Koeffizient_PlusHalf[23]-w_jMinusHalf*v_zeta_MinusHalf*Koeffizient_MinusHalf[23]+
					w_jPlusHalf*w_xi_PlusHalf*Koeffizient_PlusHalf[24]-w_jMinusHalf*w_xi_MinusHalf*Koeffizient_MinusHalf[24]+
					w_jPlusHalf*w_eta_PlusHalf*Koeffizient_PlusHalf[25]-w_jMinusHalf*w_eta_MinusHalf*Koeffizient_MinusHalf[25]+
					w_jPlusHalf*w_zeta_PlusHalf*Koeffizient_PlusHalf[26]-w_jMinusHalf*w_zeta_MinusHalf*Koeffizient_MinusHalf[26])+

					pnt_config->flag_rotation_symmetric*EnergyRotSymm;
			}
		}
	}
}

void CalcViscidFluxesInZetaDirectionDirectly(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q)
{
	int ijk,i,j,k,l,ijkMinus1,ijkPlus1;
	int ijkMAX;
	int indexMinus_zeta, indexPlus_zeta;
	FLT u_xi_MinusHalf, u_eta_MinusHalf, u_zeta_MinusHalf;
	FLT v_xi_MinusHalf, v_eta_MinusHalf, v_zeta_MinusHalf;
 	FLT w_xi_MinusHalf, w_eta_MinusHalf, w_zeta_MinusHalf;
 	FLT T_xi_MinusHalf, T_eta_MinusHalf, T_zeta_MinusHalf;
	FLT u_xi_PlusHalf, u_eta_PlusHalf, u_zeta_PlusHalf;
	FLT v_xi_PlusHalf, v_eta_PlusHalf, v_zeta_PlusHalf;
 	FLT w_xi_PlusHalf, w_eta_PlusHalf, w_zeta_PlusHalf;
 	FLT T_xi_PlusHalf, T_eta_PlusHalf, T_zeta_PlusHalf;

	FLT u_kMinusHalf, u_kPlusHalf;
	FLT v_kMinusHalf, v_kPlusHalf;
	FLT w_kMinusHalf, w_kPlusHalf;
	FLT T_kMinusHalf, T_kPlusHalf;




	ijkMAX=pnt_config->int_iMeshPointsGhostCells*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells;


	FLT Koeffizient_PlusHalf[30];
	FLT Koeffizient_MinusHalf[30];

	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijkMinus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-1;
				ijkPlus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+1;

				u_xi_MinusHalf=0.;		v_xi_MinusHalf=0.;			w_xi_MinusHalf=0.;		T_xi_MinusHalf=0.;
				u_eta_MinusHalf=0.;		v_eta_MinusHalf=0.;			w_eta_MinusHalf=0.;		T_eta_MinusHalf=0.;
				u_zeta_MinusHalf=0.;	v_zeta_MinusHalf=0.;		w_zeta_MinusHalf=0.;	T_zeta_MinusHalf=0.;
				u_xi_PlusHalf=0.;		v_xi_PlusHalf=0.;			w_xi_PlusHalf=0.;		T_xi_PlusHalf=0.;
				u_eta_PlusHalf=0.;		v_eta_PlusHalf=0.;			w_eta_PlusHalf=0.;		T_eta_PlusHalf=0.;
				u_zeta_PlusHalf=0.;		v_zeta_PlusHalf=0.;			w_zeta_PlusHalf=0.;		T_zeta_PlusHalf=0.;

				u_kMinusHalf=0.;		v_kMinusHalf=0.;			w_kMinusHalf=0.;		T_kMinusHalf=0.;
				u_kPlusHalf=0.;			v_kPlusHalf=0.;				w_kPlusHalf=0.;			T_kPlusHalf=0.;


				for(l=0;l<=SPACEORDER;l++)
				{
					indexMinus_zeta=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+(j+l-(SPACEORDER+1)/2);
					indexPlus_zeta=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+(j+l-(SPACEORDER+1)/2)+1;

					u_kMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u[indexMinus_zeta];		u_kPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u[indexPlus_zeta];
					v_kMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v[indexMinus_zeta];		v_kPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v[indexPlus_zeta];
					w_kMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w[indexMinus_zeta];		w_kPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w[indexPlus_zeta];
					T_kMinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T[indexMinus_zeta];		T_kPlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T[indexPlus_zeta];

					u_zeta_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->u[indexMinus_zeta];	u_zeta_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->u[indexPlus_zeta];
					v_zeta_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->v[indexMinus_zeta];	v_zeta_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->v[indexPlus_zeta];
					w_zeta_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->w[indexMinus_zeta];	w_zeta_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->w[indexPlus_zeta];
					T_zeta_MinusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->T[indexMinus_zeta];	T_zeta_PlusHalf+=pnt_config->ZD_ZweiteAbleitungZwischenPunkt_Koeffizient[l]*pnt_U_RK->T[indexPlus_zeta];

					u_xi_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_xi[indexMinus_zeta];	u_xi_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_xi[indexPlus_zeta];
					v_xi_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_xi[indexMinus_zeta];	v_xi_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_xi[indexPlus_zeta];
					w_xi_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_xi[indexMinus_zeta];	w_xi_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_xi[indexPlus_zeta];
					T_xi_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_xi[indexMinus_zeta];	T_xi_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_xi[indexPlus_zeta];

					u_eta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_eta[indexMinus_zeta];	u_eta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->u_eta[indexPlus_zeta];
					v_eta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_eta[indexMinus_zeta];	v_eta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->v_eta[indexPlus_zeta];
					w_eta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_eta[indexMinus_zeta];	w_eta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->w_eta[indexPlus_zeta];
					T_eta_MinusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_eta[indexMinus_zeta];	T_eta_PlusHalf+=pnt_config->ZD_Interpolation_Koeffizient[l]*pnt_U_RK->T_eta[indexPlus_zeta];
				}


				for(l=0;l<30;l++)
				{
					Koeffizient_MinusHalf[l]=
							0.5*
							(pnt_U_RK->mue[ijkMinus1]*pnt_mesh->zetaFluss_Faktor[l*ijkMAX+ijkMinus1]+
							pnt_U_RK->mue[ijk]*pnt_mesh->zetaFluss_Faktor[l*ijkMAX+ijk]);
					Koeffizient_PlusHalf[l]=
							0.5*
							(pnt_U_RK->mue[ijk]*pnt_mesh->zetaFluss_Faktor[l*ijkMAX+ijk]+
							pnt_U_RK->mue[ijkPlus1]*pnt_mesh->zetaFluss_Faktor[l*ijkMAX+ijkPlus1]);
				}



				pnt_Q->xiMomentum[ijk]=pnt_Q->xiMomentum[ijk]+
						pnt_config->Psi*(
					u_xi_PlusHalf*Koeffizient_PlusHalf[0]-u_xi_MinusHalf*Koeffizient_MinusHalf[0]+
					u_eta_PlusHalf*Koeffizient_PlusHalf[1]-u_eta_MinusHalf*Koeffizient_MinusHalf[1]+
					u_zeta_PlusHalf*Koeffizient_PlusHalf[2]-u_zeta_MinusHalf*Koeffizient_MinusHalf[2]+
					v_xi_PlusHalf*Koeffizient_PlusHalf[3]-v_xi_MinusHalf*Koeffizient_MinusHalf[3]+
					v_eta_PlusHalf*Koeffizient_PlusHalf[4]-v_eta_MinusHalf*Koeffizient_MinusHalf[4]+
					v_zeta_PlusHalf*Koeffizient_PlusHalf[5]-v_zeta_MinusHalf*Koeffizient_MinusHalf[5]+
					w_xi_PlusHalf*Koeffizient_PlusHalf[6]-w_xi_MinusHalf*Koeffizient_MinusHalf[6]+
					w_eta_PlusHalf*Koeffizient_PlusHalf[7]-w_eta_MinusHalf*Koeffizient_MinusHalf[7]+
					w_zeta_PlusHalf*Koeffizient_PlusHalf[8]-w_zeta_MinusHalf*Koeffizient_MinusHalf[8]);

				pnt_Q->etaMomentum[ijk]=pnt_Q->etaMomentum[ijk]+
						pnt_config->Psi*(
					u_xi_PlusHalf*Koeffizient_PlusHalf[9]-u_xi_MinusHalf*Koeffizient_MinusHalf[9]+
					u_eta_PlusHalf*Koeffizient_PlusHalf[10]-u_eta_MinusHalf*Koeffizient_MinusHalf[10]+
					u_zeta_PlusHalf*Koeffizient_PlusHalf[11]-u_zeta_MinusHalf*Koeffizient_MinusHalf[11]+
					v_xi_PlusHalf*Koeffizient_PlusHalf[12]-v_xi_MinusHalf*Koeffizient_MinusHalf[12]+
					v_eta_PlusHalf*Koeffizient_PlusHalf[13]-v_eta_MinusHalf*Koeffizient_MinusHalf[13]+
					v_zeta_PlusHalf*Koeffizient_PlusHalf[14]-v_zeta_MinusHalf*Koeffizient_MinusHalf[14]+
					w_xi_PlusHalf*Koeffizient_PlusHalf[15]-w_xi_MinusHalf*Koeffizient_MinusHalf[15]+
					w_eta_PlusHalf*Koeffizient_PlusHalf[16]-w_eta_MinusHalf*Koeffizient_MinusHalf[16]+
					w_zeta_PlusHalf*Koeffizient_PlusHalf[17]-w_zeta_MinusHalf*Koeffizient_MinusHalf[17]);

				pnt_Q->zetaMomentum[ijk]=pnt_Q->zetaMomentum[ijk]+
						pnt_config->Psi*(
					u_xi_PlusHalf*Koeffizient_PlusHalf[18]-u_xi_MinusHalf*Koeffizient_MinusHalf[18]+
					u_eta_PlusHalf*Koeffizient_PlusHalf[19]-u_eta_MinusHalf*Koeffizient_MinusHalf[19]+
					u_zeta_PlusHalf*Koeffizient_PlusHalf[20]-u_zeta_MinusHalf*Koeffizient_MinusHalf[20]+
					v_xi_PlusHalf*Koeffizient_PlusHalf[21]-v_xi_MinusHalf*Koeffizient_MinusHalf[21]+
					v_eta_PlusHalf*Koeffizient_PlusHalf[22]-v_eta_MinusHalf*Koeffizient_MinusHalf[22]+
					v_zeta_PlusHalf*Koeffizient_PlusHalf[23]-v_zeta_MinusHalf*Koeffizient_MinusHalf[23]+
					w_xi_PlusHalf*Koeffizient_PlusHalf[24]-w_xi_MinusHalf*Koeffizient_MinusHalf[24]+
					w_eta_PlusHalf*Koeffizient_PlusHalf[25]-w_eta_MinusHalf*Koeffizient_MinusHalf[25]+
					w_zeta_PlusHalf*Koeffizient_PlusHalf[26]-w_zeta_MinusHalf*Koeffizient_MinusHalf[26]);

				pnt_Q->Energy[ijk]=pnt_Q->Energy[ijk]+
						0.5*(pnt_config->Gamma[ijkPlus1]+pnt_config->Gamma[ijk])*(
					T_xi_PlusHalf*Koeffizient_PlusHalf[27]-T_xi_MinusHalf*Koeffizient_MinusHalf[27]+
					T_eta_PlusHalf*Koeffizient_PlusHalf[28]-T_eta_MinusHalf*Koeffizient_MinusHalf[28]+
					T_zeta_PlusHalf*Koeffizient_PlusHalf[29]-T_zeta_MinusHalf*Koeffizient_MinusHalf[29])+
						pnt_config->Psi*(
					u_kPlusHalf*u_xi_PlusHalf*Koeffizient_PlusHalf[0]-u_kMinusHalf*u_xi_MinusHalf*Koeffizient_MinusHalf[0]+
					u_kPlusHalf*u_eta_PlusHalf*Koeffizient_PlusHalf[1]-u_kMinusHalf*u_eta_MinusHalf*Koeffizient_MinusHalf[1]+
					u_kPlusHalf*u_zeta_PlusHalf*Koeffizient_PlusHalf[2]-u_kMinusHalf*u_zeta_MinusHalf*Koeffizient_MinusHalf[2]+
					u_kPlusHalf*v_xi_PlusHalf*Koeffizient_PlusHalf[3]-u_kMinusHalf*v_xi_MinusHalf*Koeffizient_MinusHalf[3]+
					u_kPlusHalf*v_eta_PlusHalf*Koeffizient_PlusHalf[4]-u_kMinusHalf*v_eta_MinusHalf*Koeffizient_MinusHalf[4]+
					u_kPlusHalf*v_zeta_PlusHalf*Koeffizient_PlusHalf[5]-u_kMinusHalf*v_zeta_MinusHalf*Koeffizient_MinusHalf[5]+
					u_kPlusHalf*w_xi_PlusHalf*Koeffizient_PlusHalf[6]-u_kMinusHalf*w_xi_MinusHalf*Koeffizient_MinusHalf[6]+
					u_kPlusHalf*w_eta_PlusHalf*Koeffizient_PlusHalf[7]-u_kMinusHalf*w_eta_MinusHalf*Koeffizient_MinusHalf[7]+
					u_kPlusHalf*w_zeta_PlusHalf*Koeffizient_PlusHalf[8]-u_kMinusHalf*w_zeta_MinusHalf*Koeffizient_MinusHalf[8]+

					v_kPlusHalf*u_xi_PlusHalf*Koeffizient_PlusHalf[9]-v_kMinusHalf*u_xi_MinusHalf*Koeffizient_MinusHalf[9]+
					v_kPlusHalf*u_eta_PlusHalf*Koeffizient_PlusHalf[10]-v_kMinusHalf*u_eta_MinusHalf*Koeffizient_MinusHalf[10]+
					v_kPlusHalf*u_zeta_PlusHalf*Koeffizient_PlusHalf[11]-v_kMinusHalf*u_zeta_MinusHalf*Koeffizient_MinusHalf[11]+
					v_kPlusHalf*v_xi_PlusHalf*Koeffizient_PlusHalf[12]-v_kMinusHalf*v_xi_MinusHalf*Koeffizient_MinusHalf[12]+
					v_kPlusHalf*v_eta_PlusHalf*Koeffizient_PlusHalf[13]-v_kMinusHalf*v_eta_MinusHalf*Koeffizient_MinusHalf[13]+
					v_kPlusHalf*v_zeta_PlusHalf*Koeffizient_PlusHalf[14]-v_kMinusHalf*v_zeta_MinusHalf*Koeffizient_MinusHalf[14]+
					v_kPlusHalf*w_xi_PlusHalf*Koeffizient_PlusHalf[15]-v_kMinusHalf*w_xi_MinusHalf*Koeffizient_MinusHalf[15]+
					v_kPlusHalf*w_eta_PlusHalf*Koeffizient_PlusHalf[16]-v_kMinusHalf*w_eta_MinusHalf*Koeffizient_MinusHalf[16]+
					v_kPlusHalf*w_zeta_PlusHalf*Koeffizient_PlusHalf[17]-v_kMinusHalf*w_zeta_MinusHalf*Koeffizient_MinusHalf[17]+

					w_kPlusHalf*u_xi_PlusHalf*Koeffizient_PlusHalf[18]-w_kMinusHalf*u_xi_MinusHalf*Koeffizient_MinusHalf[18]+
					w_kPlusHalf*u_eta_PlusHalf*Koeffizient_PlusHalf[19]-w_kMinusHalf*u_eta_MinusHalf*Koeffizient_MinusHalf[19]+
					w_kPlusHalf*u_zeta_PlusHalf*Koeffizient_PlusHalf[20]-w_kMinusHalf*u_zeta_MinusHalf*Koeffizient_MinusHalf[20]+
					w_kPlusHalf*v_xi_PlusHalf*Koeffizient_PlusHalf[21]-w_kMinusHalf*v_xi_MinusHalf*Koeffizient_MinusHalf[21]+
					w_kPlusHalf*v_eta_PlusHalf*Koeffizient_PlusHalf[22]-w_kMinusHalf*v_eta_MinusHalf*Koeffizient_MinusHalf[22]+
					w_kPlusHalf*v_zeta_PlusHalf*Koeffizient_PlusHalf[23]-w_kMinusHalf*v_zeta_MinusHalf*Koeffizient_MinusHalf[23]+
					w_kPlusHalf*w_xi_PlusHalf*Koeffizient_PlusHalf[24]-w_kMinusHalf*w_xi_MinusHalf*Koeffizient_MinusHalf[24]+
					w_kPlusHalf*w_eta_PlusHalf*Koeffizient_PlusHalf[25]-w_kMinusHalf*w_eta_MinusHalf*Koeffizient_MinusHalf[25]+
					w_kPlusHalf*w_zeta_PlusHalf*Koeffizient_PlusHalf[26]-w_kMinusHalf*w_zeta_MinusHalf*Koeffizient_MinusHalf[26]);
			}
		}
	}
}


#include "SHOCK.h"
#include "Functions.h"
#include "ZD.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int i,j,k;
int ijk;
int ijkPlus1,ijkPlus2,ijkPlus3,ijkPlus4,ijkPlus5;
int ijkMinus1,ijkMinus2,ijkMinus3,ijkMinus4,ijkMinus5;
int iPlus1jk,iPlus2jk,iPlus3jk,iPlus4jk,iPlus5jk;
int iMinus1jk,iMinus2jk,iMinus3jk,iMinus4jk,iMinus5jk;
int ijPlus1k,ijPlus2k,ijPlus3k,ijPlus4k,ijPlus5k;
int ijMinus1k,ijMinus2k,ijMinus3k,ijMinus4k,ijMinus5k;

FLT xImpuls_xi;
FLT yImpuls_xi;
FLT zImpuls_xi;
FLT Energie_xi;

FLT xImpuls_eta;
FLT yImpuls_eta;
FLT zImpuls_eta;
FLT Energie_eta;

FLT xImpuls_zeta;
FLT yImpuls_zeta;
FLT zImpuls_zeta;
FLT Energie_zeta;

FLT xImpuls_xi_help[11];
FLT yImpuls_xi_help[11];
FLT zImpuls_xi_help[11];
FLT Energie_xi_help[11];

FLT xImpuls_eta_help[11];
FLT yImpuls_eta_help[11];
FLT zImpuls_eta_help[11];
FLT Energie_eta_help[11];

FLT xImpuls_zeta_help[11];
FLT yImpuls_zeta_help[11];
FLT zImpuls_zeta_help[11];
FLT Energie_zeta_help[11];

void CalcDeviationsForViscidFluxes(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD)
{
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iPlus1jk=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iPlus2jk=(i+2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iPlus3jk=(i+3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iMinus2jk=(i-2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iMinus3jk=(i-3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;
				ijPlus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+2)*pnt_config->int_kMeshPointsGhostCells+k;
				ijPlus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+3)*pnt_config->int_kMeshPointsGhostCells+k;
				ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;
				ijMinus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-2)*pnt_config->int_kMeshPointsGhostCells+k;
				ijMinus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-3)*pnt_config->int_kMeshPointsGhostCells+k;
				ijkPlus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+1;
				ijkPlus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+2;
				ijkPlus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+3;
				ijkMinus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-1;
				ijkMinus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-2;
				ijkMinus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-3;

#if SPACEORDER==9
					iPlus4jk=(i+4)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					iPlus5jk=(i+5)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					iMinus4jk=(i-4)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					iMinus5jk=(i-5)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					ijPlus4k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+4)*pnt_config->int_kMeshPointsGhostCells+k;
					ijPlus5k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+5)*pnt_config->int_kMeshPointsGhostCells+k;
					ijMinus4k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-4)*pnt_config->int_kMeshPointsGhostCells+k;
					ijMinus5k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-5)*pnt_config->int_kMeshPointsGhostCells+k;
					ijkPlus4=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+4;
					ijkPlus5=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+5;
					ijkMinus4=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-4;
					ijkMinus5=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-5;
#endif

#if SPACEORDER==5
//					Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
					pnt_ZD->u_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->u[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[iPlus3jk]);
					pnt_ZD->v_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->v[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[iPlus3jk]);
					pnt_ZD->T_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->T[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[iPlus3jk]);

//					Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
					pnt_ZD->u_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->u[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[ijPlus3k]);
					pnt_ZD->v_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->v[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[ijPlus3k]);
					pnt_ZD->T_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->T[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[ijPlus3k]);

#if MESHDIMENSIONS==3

						pnt_ZD->w_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->w[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[iPlus3jk]);

						pnt_ZD->w_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->w[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[ijPlus3k]);

//						Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
						pnt_ZD->u_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->u[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[ijkPlus3]);
						pnt_ZD->v_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->v[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[ijkPlus3]);
						pnt_ZD->w_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->w[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[ijkPlus3]);
						pnt_ZD->T_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->T[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[ijkPlus3]);
#endif
#endif

#if SPACEORDER==9
//					Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
					pnt_ZD->u_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[iMinus5jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[iMinus4jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->u[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->u[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->u[iPlus3jk]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->u[iPlus4jk]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->u[iPlus5jk]);
					pnt_ZD->v_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[iMinus5jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[iMinus4jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->v[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->v[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->v[iPlus3jk]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->v[iPlus4jk]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->v[iPlus5jk]);
					pnt_ZD->T_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[iMinus5jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[iMinus4jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->T[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->T[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->T[iPlus3jk]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->T[iPlus4jk]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->T[iPlus5jk]);

//					Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
					pnt_ZD->u_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[ijMinus5k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[ijMinus4k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->u[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->u[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->u[ijPlus3k]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->u[ijPlus4k]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->u[ijPlus5k]);
					pnt_ZD->v_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[ijMinus5k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[ijMinus4k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->v[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->v[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->v[ijPlus3k]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->v[ijPlus4k]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->v[ijPlus5k]);
					pnt_ZD->T_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[ijMinus5k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[ijMinus4k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->T[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->T[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->T[ijPlus3k]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->T[ijPlus4k]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->T[ijPlus5k]);

#if MESHDIMENSIONS==3

						pnt_ZD->w_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[iMinus5jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[iMinus4jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->w[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->w[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->w[iPlus3jk]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->w[iPlus4jk]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->w[iPlus5jk]);

						pnt_ZD->w_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[ijMinus5k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[ijMinus4k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->w[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->w[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->w[ijPlus3k]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->w[ijPlus4k]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->w[ijPlus5k]);

//						Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
						pnt_ZD->u_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[ijkMinus5]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[ijkMinus4]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->u[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->u[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->u[ijkPlus3]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->u[ijkPlus4]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->u[ijkPlus5]);
						pnt_ZD->v_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[ijkMinus5]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[ijkMinus4]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->v[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->v[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->v[ijkPlus3]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->v[ijkPlus4]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->v[ijkPlus5]);
						pnt_ZD->w_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[ijkMinus5]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[ijkMinus4]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->w[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->w[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->w[ijkPlus3]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->w[ijkPlus4]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->w[ijkPlus5]);
						pnt_ZD->T_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[ijkMinus5]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[ijkMinus4]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->T[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->T[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->T[ijkPlus3]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->T[ijkPlus4]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->T[ijkPlus5]);
#endif
#endif

			}
		}
	}
}

void CalcDeviationsForDirectViscidFluxes(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK)
{
	if(pnt_config->flag_IBC==1)
	{
		IBC_ApplyBC4FluxInXi(
				pnt_config,
				pnt_mesh,
				pnt_U_RK);
	}
//	Ableitungen in Xi-Richtung
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iPlus1jk=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iPlus2jk=(i+2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iPlus3jk=(i+3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iMinus2jk=(i-2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iMinus3jk=(i-3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;


#if SPACEORDER==9
					iPlus4jk=(i+4)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					iPlus5jk=(i+5)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					iMinus4jk=(i-4)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					iMinus5jk=(i-5)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
#endif

#if SPACEORDER==5
//					Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
					pnt_U_RK->u_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->u[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[iPlus3jk]);
					pnt_U_RK->v_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->v[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[iPlus3jk]);
					pnt_U_RK->T_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->T[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[iPlus3jk]);

#if MESHDIMENSIONS==3
						pnt_U_RK->w_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->w[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[iPlus3jk]);
#endif
#endif

#if SPACEORDER==9
//					Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
					pnt_U_RK->u_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[iMinus5jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[iMinus4jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->u[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->u[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->u[iPlus3jk]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->u[iPlus4jk]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->u[iPlus5jk]);
					pnt_U_RK->v_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[iMinus5jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[iMinus4jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->v[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->v[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->v[iPlus3jk]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->v[iPlus4jk]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->v[iPlus5jk]);
					pnt_U_RK->T_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[iMinus5jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[iMinus4jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->T[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->T[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->T[iPlus3jk]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->T[iPlus4jk]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->T[iPlus5jk]);

#if MESHDIMENSIONS==3
						pnt_U_RK->w_xi[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[iMinus5jk]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[iMinus4jk]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[iMinus3jk]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->w[iMinus2jk]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[iMinus1jk]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[iPlus1jk]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->w[iPlus2jk]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->w[iPlus3jk]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->w[iPlus4jk]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->w[iPlus5jk]);
#endif
#endif

			}
		}
	}
	
	if(pnt_config->flag_IBC==1)
	{
		IBC_ApplyBC4FluxInEta(
				pnt_config,
				pnt_mesh,
				pnt_U_RK);
	}
//	Ableitungen in Eta-Richtung
	for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartGhosts; k <= pnt_config->int_kEndGhosts; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;
				ijPlus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+2)*pnt_config->int_kMeshPointsGhostCells+k;
				ijPlus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+3)*pnt_config->int_kMeshPointsGhostCells+k;
				ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;
				ijMinus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-2)*pnt_config->int_kMeshPointsGhostCells+k;
				ijMinus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-3)*pnt_config->int_kMeshPointsGhostCells+k;

#if SPACEORDER==9
					ijPlus4k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+4)*pnt_config->int_kMeshPointsGhostCells+k;
					ijPlus5k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+5)*pnt_config->int_kMeshPointsGhostCells+k;
					ijMinus4k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-4)*pnt_config->int_kMeshPointsGhostCells+k;
					ijMinus5k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-5)*pnt_config->int_kMeshPointsGhostCells+k;
#endif

#if SPACEORDER==5
//					Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
					pnt_U_RK->u_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->u[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[ijPlus3k]);
					pnt_U_RK->v_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->v[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[ijPlus3k]);
					pnt_U_RK->T_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->T[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[ijPlus3k]);

#if MESHDIMENSIONS==3
						pnt_U_RK->w_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->w[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[ijPlus3k]);
#endif
#endif

#if SPACEORDER==9
//					Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
					pnt_U_RK->u_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[ijMinus5k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[ijMinus4k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->u[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->u[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->u[ijPlus3k]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->u[ijPlus4k]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->u[ijPlus5k]);
					pnt_U_RK->v_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[ijMinus5k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[ijMinus4k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->v[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->v[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->v[ijPlus3k]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->v[ijPlus4k]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->v[ijPlus5k]);
					pnt_U_RK->T_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[ijMinus5k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[ijMinus4k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->T[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->T[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->T[ijPlus3k]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->T[ijPlus4k]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->T[ijPlus5k]);

#if MESHDIMENSIONS==3
						pnt_U_RK->w_eta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[ijMinus5k]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[ijMinus4k]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[ijMinus3k]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->w[ijMinus2k]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[ijMinus1k]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[ijPlus1k]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->w[ijPlus2k]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->w[ijPlus3k]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->w[ijPlus4k]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->w[ijPlus5k]);
#endif
#endif

			}
		}
	}

	if(pnt_config->flag_IBC==1)
	{
		IBC_ApplyBC4FluxInZeta(
				pnt_config,
				pnt_mesh,
				pnt_U_RK);
	}
//	Ableitungen in Zeta-Richtung
#if MESHDIMENSIONS==3
		for (i=pnt_config->int_iStartGhosts; i <= pnt_config->int_iEndGhosts; i++)
		{
			for (j=pnt_config->int_jStartGhosts; j <= pnt_config->int_jEndGhosts; j++)
			{
				for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
				{
					ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					ijkPlus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+1;
					ijkPlus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+2;
					ijkPlus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+3;
					ijkMinus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-1;
					ijkMinus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-2;
					ijkMinus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-3;

#if SPACEORDER==9
						ijkPlus4=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+4;
						ijkPlus5=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+5;
						ijkMinus4=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-4;
						ijkMinus5=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-5;
#endif

#if SPACEORDER==5
//						Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
						pnt_U_RK->u_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->u[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[ijkPlus3]);
						pnt_U_RK->v_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->v[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[ijkPlus3]);
						pnt_U_RK->w_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->w[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[ijkPlus3]);
						pnt_U_RK->T_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[5]*pnt_U_RK->T[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[ijkPlus3]);
#endif

#if SPACEORDER==9
//						Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
						pnt_U_RK->u_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->u[ijkMinus5]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->u[ijkMinus4]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->u[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->u[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->u[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->u[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->u[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->u[ijkPlus3]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->u[ijkPlus4]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->u[ijkPlus5]);
						pnt_U_RK->v_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->v[ijkMinus5]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->v[ijkMinus4]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->v[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->v[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->v[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->v[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->v[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->v[ijkPlus3]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->v[ijkPlus4]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->v[ijkPlus5]);
						pnt_U_RK->w_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->w[ijkMinus5]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->w[ijkMinus4]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->w[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->w[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->w[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->w[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->w[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->w[ijkPlus3]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->w[ijkPlus4]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->w[ijkPlus5]);
						pnt_U_RK->T_zeta[ijk]=(pnt_config->ZD_Ableitung_Koeffizient[0]*pnt_U_RK->T[ijkMinus5]+pnt_config->ZD_Ableitung_Koeffizient[1]*pnt_U_RK->T[ijkMinus4]+pnt_config->ZD_Ableitung_Koeffizient[2]*pnt_U_RK->T[ijkMinus3]+pnt_config->ZD_Ableitung_Koeffizient[3]*pnt_U_RK->T[ijkMinus2]+pnt_config->ZD_Ableitung_Koeffizient[4]*pnt_U_RK->T[ijkMinus1]+pnt_config->ZD_Ableitung_Koeffizient[6]*pnt_U_RK->T[ijkPlus1]+pnt_config->ZD_Ableitung_Koeffizient[7]*pnt_U_RK->T[ijkPlus2]+pnt_config->ZD_Ableitung_Koeffizient[8]*pnt_U_RK->T[ijkPlus3]+pnt_config->ZD_Ableitung_Koeffizient[9]*pnt_U_RK->T[ijkPlus4]+pnt_config->ZD_Ableitung_Koeffizient[10]*pnt_U_RK->T[ijkPlus5]);
#endif

				}
			}
		}
#endif
}

void CalcTauQForViscidFluxes(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD)
{
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;

				CalcTauQValues(ijk,pnt_config,pnt_mesh,pnt_U_RK,pnt_ZD);
			}
		}
	}
}


void CalcViscidFluxesInZetaDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD,
		struct strct_Flux * pnt_Q)
{
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iPlus1jk=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iPlus2jk=(i+2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iPlus3jk=(i+3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iMinus2jk=(i-2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iMinus3jk=(i-3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijPlus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+2)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijPlus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+3)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijMinus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-2)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijMinus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-3)*pnt_config->int_kMeshPointsGhostCells+k;
				ijkPlus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+1;
				ijkPlus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+2;
				ijkPlus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+3;
				ijkMinus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-1;
				ijkMinus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-2;
				ijkMinus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-3;

#if SPACEORDER==9
//					iPlus4jk=(i+4)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//					iPlus5jk=(i+5)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//					iMinus4jk=(i-4)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//					iMinus5jk=(i-5)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//					ijPlus4k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+4)*pnt_config->int_kMeshPointsGhostCells+k;
//					ijPlus5k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+5)*pnt_config->int_kMeshPointsGhostCells+k;
//					ijMinus4k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-4)*pnt_config->int_kMeshPointsGhostCells+k;
//					ijMinus5k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-5)*pnt_config->int_kMeshPointsGhostCells+k;
					ijkPlus4=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+4;
					ijkPlus5=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+5;
					ijkMinus4=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-4;
					ijkMinus5=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-5;
#endif


#if SPACEORDER==9
					xImpuls_zeta_help[0]=(pnt_mesh->jacobian[ijkMinus5]*(pnt_mesh->zeta_x[ijkMinus5]*pnt_ZD->tau_xx[ijkMinus5]+pnt_mesh->zeta_y[ijkMinus5]*pnt_ZD->tau_xy[ijkMinus5]+pnt_mesh->zeta_z[ijkMinus5]*pnt_ZD->tau_xz[ijkMinus5]));
					xImpuls_zeta_help[1]=(pnt_mesh->jacobian[ijkMinus4]*(pnt_mesh->zeta_x[ijkMinus4]*pnt_ZD->tau_xx[ijkMinus4]+pnt_mesh->zeta_y[ijkMinus4]*pnt_ZD->tau_xy[ijkMinus4]+pnt_mesh->zeta_z[ijkMinus4]*pnt_ZD->tau_xz[ijkMinus4]));
					xImpuls_zeta_help[2]=(pnt_mesh->jacobian[ijkMinus3]*(pnt_mesh->zeta_x[ijkMinus3]*pnt_ZD->tau_xx[ijkMinus3]+pnt_mesh->zeta_y[ijkMinus3]*pnt_ZD->tau_xy[ijkMinus3]+pnt_mesh->zeta_z[ijkMinus3]*pnt_ZD->tau_xz[ijkMinus3]));
					xImpuls_zeta_help[3]=(pnt_mesh->jacobian[ijkMinus2]*(pnt_mesh->zeta_x[ijkMinus2]*pnt_ZD->tau_xx[ijkMinus2]+pnt_mesh->zeta_y[ijkMinus2]*pnt_ZD->tau_xy[ijkMinus2]+pnt_mesh->zeta_z[ijkMinus2]*pnt_ZD->tau_xz[ijkMinus2]));
					xImpuls_zeta_help[4]=(pnt_mesh->jacobian[ijkMinus1]*(pnt_mesh->zeta_x[ijkMinus1]*pnt_ZD->tau_xx[ijkMinus1]+pnt_mesh->zeta_y[ijkMinus1]*pnt_ZD->tau_xy[ijkMinus1]+pnt_mesh->zeta_z[ijkMinus1]*pnt_ZD->tau_xz[ijkMinus1]));
					xImpuls_zeta_help[6]=(pnt_mesh->jacobian[ijkPlus1]*(pnt_mesh->zeta_x[ijkPlus1]*pnt_ZD->tau_xx[ijkPlus1]+pnt_mesh->zeta_y[ijkPlus1]*pnt_ZD->tau_xy[ijkPlus1]+pnt_mesh->zeta_z[ijkPlus1]*pnt_ZD->tau_xz[ijkPlus1]));
					xImpuls_zeta_help[7]=(pnt_mesh->jacobian[ijkPlus2]*(pnt_mesh->zeta_x[ijkPlus2]*pnt_ZD->tau_xx[ijkPlus2]+pnt_mesh->zeta_y[ijkPlus2]*pnt_ZD->tau_xy[ijkPlus2]+pnt_mesh->zeta_z[ijkPlus2]*pnt_ZD->tau_xz[ijkPlus2]));
					xImpuls_zeta_help[8]=(pnt_mesh->jacobian[ijkPlus3]*(pnt_mesh->zeta_x[ijkPlus3]*pnt_ZD->tau_xx[ijkPlus3]+pnt_mesh->zeta_y[ijkPlus3]*pnt_ZD->tau_xy[ijkPlus3]+pnt_mesh->zeta_z[ijkPlus3]*pnt_ZD->tau_xz[ijkPlus3]));
					xImpuls_zeta_help[9]=(pnt_mesh->jacobian[ijkPlus4]*(pnt_mesh->zeta_x[ijkPlus4]*pnt_ZD->tau_xx[ijkPlus4]+pnt_mesh->zeta_y[ijkPlus4]*pnt_ZD->tau_xy[ijkPlus4]+pnt_mesh->zeta_z[ijkPlus4]*pnt_ZD->tau_xz[ijkPlus4]));
					xImpuls_zeta_help[10]=(pnt_mesh->jacobian[ijkPlus5]*(pnt_mesh->zeta_x[ijkPlus5]*pnt_ZD->tau_xx[ijkPlus5]+pnt_mesh->zeta_y[ijkPlus5]*pnt_ZD->tau_xy[ijkPlus5]+pnt_mesh->zeta_z[ijkPlus5]*pnt_ZD->tau_xz[ijkPlus5]));

					xImpuls_zeta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*xImpuls_zeta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*xImpuls_zeta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*xImpuls_zeta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[3]*xImpuls_zeta_help[3]+pnt_config->ZD_Ableitung_Koeffizient[4]*xImpuls_zeta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[6]*xImpuls_zeta_help[6]+pnt_config->ZD_Ableitung_Koeffizient[7]*xImpuls_zeta_help[7]+pnt_config->ZD_Ableitung_Koeffizient[8]*xImpuls_zeta_help[8]+pnt_config->ZD_Ableitung_Koeffizient[9]*xImpuls_zeta_help[9]+pnt_config->ZD_Ableitung_Koeffizient[10]*xImpuls_zeta_help[10]);


					yImpuls_zeta_help[0]=(pnt_mesh->jacobian[ijkMinus5]*(pnt_mesh->zeta_x[ijkMinus5]*pnt_ZD->tau_xy[ijkMinus5]+pnt_mesh->zeta_y[ijkMinus5]*pnt_ZD->tau_yy[ijkMinus5]+pnt_mesh->zeta_z[ijkMinus5]*pnt_ZD->tau_yz[ijkMinus5]));
					yImpuls_zeta_help[1]=(pnt_mesh->jacobian[ijkMinus4]*(pnt_mesh->zeta_x[ijkMinus4]*pnt_ZD->tau_xy[ijkMinus4]+pnt_mesh->zeta_y[ijkMinus4]*pnt_ZD->tau_yy[ijkMinus4]+pnt_mesh->zeta_z[ijkMinus4]*pnt_ZD->tau_yz[ijkMinus4]));
					yImpuls_zeta_help[2]=(pnt_mesh->jacobian[ijkMinus3]*(pnt_mesh->zeta_x[ijkMinus3]*pnt_ZD->tau_xy[ijkMinus3]+pnt_mesh->zeta_y[ijkMinus3]*pnt_ZD->tau_yy[ijkMinus3]+pnt_mesh->zeta_z[ijkMinus3]*pnt_ZD->tau_yz[ijkMinus3]));
					yImpuls_zeta_help[3]=(pnt_mesh->jacobian[ijkMinus2]*(pnt_mesh->zeta_x[ijkMinus2]*pnt_ZD->tau_xy[ijkMinus2]+pnt_mesh->zeta_y[ijkMinus2]*pnt_ZD->tau_yy[ijkMinus2]+pnt_mesh->zeta_z[ijkMinus2]*pnt_ZD->tau_yz[ijkMinus2]));
					yImpuls_zeta_help[4]=(pnt_mesh->jacobian[ijkMinus1]*(pnt_mesh->zeta_x[ijkMinus1]*pnt_ZD->tau_xy[ijkMinus1]+pnt_mesh->zeta_y[ijkMinus1]*pnt_ZD->tau_yy[ijkMinus1]+pnt_mesh->zeta_z[ijkMinus1]*pnt_ZD->tau_yz[ijkMinus1]));
					yImpuls_zeta_help[6]=(pnt_mesh->jacobian[ijkPlus1]*(pnt_mesh->zeta_x[ijkPlus1]*pnt_ZD->tau_xy[ijkPlus1]+pnt_mesh->zeta_y[ijkPlus1]*pnt_ZD->tau_yy[ijkPlus1]+pnt_mesh->zeta_z[ijkPlus1]*pnt_ZD->tau_yz[ijkPlus1]));
					yImpuls_zeta_help[7]=(pnt_mesh->jacobian[ijkPlus2]*(pnt_mesh->zeta_x[ijkPlus2]*pnt_ZD->tau_xy[ijkPlus2]+pnt_mesh->zeta_y[ijkPlus2]*pnt_ZD->tau_yy[ijkPlus2]+pnt_mesh->zeta_z[ijkPlus2]*pnt_ZD->tau_yz[ijkPlus2]));
					yImpuls_zeta_help[8]=(pnt_mesh->jacobian[ijkPlus3]*(pnt_mesh->zeta_x[ijkPlus3]*pnt_ZD->tau_xy[ijkPlus3]+pnt_mesh->zeta_y[ijkPlus3]*pnt_ZD->tau_yy[ijkPlus3]+pnt_mesh->zeta_z[ijkPlus3]*pnt_ZD->tau_yz[ijkPlus3]));
					yImpuls_zeta_help[9]=(pnt_mesh->jacobian[ijkPlus4]*(pnt_mesh->zeta_x[ijkPlus4]*pnt_ZD->tau_xy[ijkPlus4]+pnt_mesh->zeta_y[ijkPlus4]*pnt_ZD->tau_yy[ijkPlus4]+pnt_mesh->zeta_z[ijkPlus4]*pnt_ZD->tau_yz[ijkPlus4]));
					yImpuls_zeta_help[10]=(pnt_mesh->jacobian[ijkPlus5]*(pnt_mesh->zeta_x[ijkPlus5]*pnt_ZD->tau_xy[ijkPlus5]+pnt_mesh->zeta_y[ijkPlus5]*pnt_ZD->tau_yy[ijkPlus5]+pnt_mesh->zeta_z[ijkPlus5]*pnt_ZD->tau_yz[ijkPlus5]));

					yImpuls_zeta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*yImpuls_zeta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*yImpuls_zeta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*yImpuls_zeta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[3]*yImpuls_zeta_help[3]+pnt_config->ZD_Ableitung_Koeffizient[4]*yImpuls_zeta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[6]*yImpuls_zeta_help[6]+pnt_config->ZD_Ableitung_Koeffizient[7]*yImpuls_zeta_help[7]+pnt_config->ZD_Ableitung_Koeffizient[8]*yImpuls_zeta_help[8]+pnt_config->ZD_Ableitung_Koeffizient[9]*yImpuls_zeta_help[9]+pnt_config->ZD_Ableitung_Koeffizient[10]*yImpuls_zeta_help[10]);


#if MESHDIMENSIONS==3
						zImpuls_zeta_help[0]=(pnt_mesh->jacobian[ijkMinus5]*(pnt_mesh->zeta_x[ijkMinus5]*pnt_ZD->tau_xz[ijkMinus5]+pnt_mesh->zeta_y[ijkMinus5]*pnt_ZD->tau_yz[ijkMinus5]+pnt_mesh->zeta_z[ijkMinus5]*pnt_ZD->tau_zz[ijkMinus5]));
						zImpuls_zeta_help[1]=(pnt_mesh->jacobian[ijkMinus4]*(pnt_mesh->zeta_x[ijkMinus4]*pnt_ZD->tau_xz[ijkMinus4]+pnt_mesh->zeta_y[ijkMinus4]*pnt_ZD->tau_yz[ijkMinus4]+pnt_mesh->zeta_z[ijkMinus4]*pnt_ZD->tau_zz[ijkMinus4]));
						zImpuls_zeta_help[2]=(pnt_mesh->jacobian[ijkMinus3]*(pnt_mesh->zeta_x[ijkMinus3]*pnt_ZD->tau_xz[ijkMinus3]+pnt_mesh->zeta_y[ijkMinus3]*pnt_ZD->tau_yz[ijkMinus3]+pnt_mesh->zeta_z[ijkMinus3]*pnt_ZD->tau_zz[ijkMinus3]));
						zImpuls_zeta_help[3]=(pnt_mesh->jacobian[ijkMinus2]*(pnt_mesh->zeta_x[ijkMinus2]*pnt_ZD->tau_xz[ijkMinus2]+pnt_mesh->zeta_y[ijkMinus2]*pnt_ZD->tau_yz[ijkMinus2]+pnt_mesh->zeta_z[ijkMinus2]*pnt_ZD->tau_zz[ijkMinus2]));
						zImpuls_zeta_help[4]=(pnt_mesh->jacobian[ijkMinus1]*(pnt_mesh->zeta_x[ijkMinus1]*pnt_ZD->tau_xz[ijkMinus1]+pnt_mesh->zeta_y[ijkMinus1]*pnt_ZD->tau_yz[ijkMinus1]+pnt_mesh->zeta_z[ijkMinus1]*pnt_ZD->tau_zz[ijkMinus1]));
						zImpuls_zeta_help[6]=(pnt_mesh->jacobian[ijkPlus1]*(pnt_mesh->zeta_x[ijkPlus1]*pnt_ZD->tau_xz[ijkPlus1]+pnt_mesh->zeta_y[ijkPlus1]*pnt_ZD->tau_yz[ijkPlus1]+pnt_mesh->zeta_z[ijkPlus1]*pnt_ZD->tau_zz[ijkPlus1]));
						zImpuls_zeta_help[7]=(pnt_mesh->jacobian[ijkPlus2]*(pnt_mesh->zeta_x[ijkPlus2]*pnt_ZD->tau_xz[ijkPlus2]+pnt_mesh->zeta_y[ijkPlus2]*pnt_ZD->tau_yz[ijkPlus2]+pnt_mesh->zeta_z[ijkPlus2]*pnt_ZD->tau_zz[ijkPlus2]));
						zImpuls_zeta_help[8]=(pnt_mesh->jacobian[ijkPlus3]*(pnt_mesh->zeta_x[ijkPlus3]*pnt_ZD->tau_xz[ijkPlus3]+pnt_mesh->zeta_y[ijkPlus3]*pnt_ZD->tau_yz[ijkPlus3]+pnt_mesh->zeta_z[ijkPlus3]*pnt_ZD->tau_zz[ijkPlus3]));
						zImpuls_zeta_help[9]=(pnt_mesh->jacobian[ijkPlus4]*(pnt_mesh->zeta_x[ijkPlus4]*pnt_ZD->tau_xz[ijkPlus4]+pnt_mesh->zeta_y[ijkPlus4]*pnt_ZD->tau_yz[ijkPlus4]+pnt_mesh->zeta_z[ijkPlus4]*pnt_ZD->tau_zz[ijkPlus4]));
						zImpuls_zeta_help[10]=(pnt_mesh->jacobian[ijkPlus5]*(pnt_mesh->zeta_x[ijkPlus5]*pnt_ZD->tau_xz[ijkPlus5]+pnt_mesh->zeta_y[ijkPlus5]*pnt_ZD->tau_yz[ijkPlus5]+pnt_mesh->zeta_z[ijkPlus5]*pnt_ZD->tau_zz[ijkPlus5]));

						zImpuls_zeta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*zImpuls_zeta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*zImpuls_zeta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*zImpuls_zeta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[3]*zImpuls_zeta_help[3]+pnt_config->ZD_Ableitung_Koeffizient[4]*zImpuls_zeta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[6]*zImpuls_zeta_help[6]+pnt_config->ZD_Ableitung_Koeffizient[7]*zImpuls_zeta_help[7]+pnt_config->ZD_Ableitung_Koeffizient[8]*zImpuls_zeta_help[8]+pnt_config->ZD_Ableitung_Koeffizient[9]*zImpuls_zeta_help[9]+pnt_config->ZD_Ableitung_Koeffizient[10]*zImpuls_zeta_help[10]);
#endif


					Energie_zeta_help[0]=(pnt_mesh->jacobian[ijkMinus5]*(pnt_mesh->zeta_x[ijkMinus5]*pnt_ZD->q_x[ijkMinus5]+pnt_mesh->zeta_y[ijkMinus5]*pnt_ZD->q_y[ijkMinus5]+pnt_mesh->zeta_z[ijkMinus5]*pnt_ZD->q_z[ijkMinus5]));
					Energie_zeta_help[1]=(pnt_mesh->jacobian[ijkMinus4]*(pnt_mesh->zeta_x[ijkMinus4]*pnt_ZD->q_x[ijkMinus4]+pnt_mesh->zeta_y[ijkMinus4]*pnt_ZD->q_y[ijkMinus4]+pnt_mesh->zeta_z[ijkMinus4]*pnt_ZD->q_z[ijkMinus4]));
					Energie_zeta_help[2]=(pnt_mesh->jacobian[ijkMinus3]*(pnt_mesh->zeta_x[ijkMinus3]*pnt_ZD->q_x[ijkMinus3]+pnt_mesh->zeta_y[ijkMinus3]*pnt_ZD->q_y[ijkMinus3]+pnt_mesh->zeta_z[ijkMinus3]*pnt_ZD->q_z[ijkMinus3]));
					Energie_zeta_help[3]=(pnt_mesh->jacobian[ijkMinus2]*(pnt_mesh->zeta_x[ijkMinus2]*pnt_ZD->q_x[ijkMinus2]+pnt_mesh->zeta_y[ijkMinus2]*pnt_ZD->q_y[ijkMinus2]+pnt_mesh->zeta_z[ijkMinus2]*pnt_ZD->q_z[ijkMinus2]));
					Energie_zeta_help[4]=(pnt_mesh->jacobian[ijkMinus1]*(pnt_mesh->zeta_x[ijkMinus1]*pnt_ZD->q_x[ijkMinus1]+pnt_mesh->zeta_y[ijkMinus1]*pnt_ZD->q_y[ijkMinus1]+pnt_mesh->zeta_z[ijkMinus1]*pnt_ZD->q_z[ijkMinus1]));
					Energie_zeta_help[6]=(pnt_mesh->jacobian[ijkPlus1]*(pnt_mesh->zeta_x[ijkPlus1]*pnt_ZD->q_x[ijkPlus1]+pnt_mesh->zeta_y[ijkPlus1]*pnt_ZD->q_y[ijkPlus1]+pnt_mesh->zeta_z[ijkPlus1]*pnt_ZD->q_z[ijkPlus1]));
					Energie_zeta_help[7]=(pnt_mesh->jacobian[ijkPlus2]*(pnt_mesh->zeta_x[ijkPlus2]*pnt_ZD->q_x[ijkPlus2]+pnt_mesh->zeta_y[ijkPlus2]*pnt_ZD->q_y[ijkPlus2]+pnt_mesh->zeta_z[ijkPlus2]*pnt_ZD->q_z[ijkPlus2]));
					Energie_zeta_help[8]=(pnt_mesh->jacobian[ijkPlus3]*(pnt_mesh->zeta_x[ijkPlus3]*pnt_ZD->q_x[ijkPlus3]+pnt_mesh->zeta_y[ijkPlus3]*pnt_ZD->q_y[ijkPlus3]+pnt_mesh->zeta_z[ijkPlus3]*pnt_ZD->q_z[ijkPlus3]));
					Energie_zeta_help[9]=(pnt_mesh->jacobian[ijkPlus4]*(pnt_mesh->zeta_x[ijkPlus4]*pnt_ZD->q_x[ijkPlus4]+pnt_mesh->zeta_y[ijkPlus4]*pnt_ZD->q_y[ijkPlus4]+pnt_mesh->zeta_z[ijkPlus4]*pnt_ZD->q_z[ijkPlus4]));
					Energie_zeta_help[10]=(pnt_mesh->jacobian[ijkPlus5]*(pnt_mesh->zeta_x[ijkPlus5]*pnt_ZD->q_x[ijkPlus5]+pnt_mesh->zeta_y[ijkPlus5]*pnt_ZD->q_y[ijkPlus5]+pnt_mesh->zeta_z[ijkPlus5]*pnt_ZD->q_z[ijkPlus5]));

					Energie_zeta=
							pnt_config->ZD_Ableitung_Koeffizient[0]*(Energie_zeta_help[0]*pnt_config->Gamma[ijkMinus5]+pnt_config->Psi*(xImpuls_zeta_help[0]*pnt_U_RK->u[ijkMinus5]+yImpuls_zeta_help[0]*pnt_U_RK->v[ijkMinus5]+zImpuls_zeta_help[0]*pnt_U_RK->w[ijkMinus5]))+
							pnt_config->ZD_Ableitung_Koeffizient[1]*(Energie_zeta_help[1]*pnt_config->Gamma[ijkMinus4]+pnt_config->Psi*(xImpuls_zeta_help[1]*pnt_U_RK->u[ijkMinus4]+yImpuls_zeta_help[1]*pnt_U_RK->v[ijkMinus4]+zImpuls_zeta_help[1]*pnt_U_RK->w[ijkMinus4]))+
							pnt_config->ZD_Ableitung_Koeffizient[2]*(Energie_zeta_help[2]*pnt_config->Gamma[ijkMinus3]+pnt_config->Psi*(xImpuls_zeta_help[2]*pnt_U_RK->u[ijkMinus3]+yImpuls_zeta_help[2]*pnt_U_RK->v[ijkMinus3]+zImpuls_zeta_help[2]*pnt_U_RK->w[ijkMinus3]))+
							pnt_config->ZD_Ableitung_Koeffizient[3]*(Energie_zeta_help[3]*pnt_config->Gamma[ijkMinus2]+pnt_config->Psi*(xImpuls_zeta_help[3]*pnt_U_RK->u[ijkMinus2]+yImpuls_zeta_help[3]*pnt_U_RK->v[ijkMinus2]+zImpuls_zeta_help[3]*pnt_U_RK->w[ijkMinus2]))+
							pnt_config->ZD_Ableitung_Koeffizient[4]*(Energie_zeta_help[4]*pnt_config->Gamma[ijkMinus1]+pnt_config->Psi*(xImpuls_zeta_help[4]*pnt_U_RK->u[ijkMinus1]+yImpuls_zeta_help[4]*pnt_U_RK->v[ijkMinus1]+zImpuls_zeta_help[4]*pnt_U_RK->w[ijkMinus1]))+
							pnt_config->ZD_Ableitung_Koeffizient[6]*(Energie_zeta_help[6]*pnt_config->Gamma[ijkPlus1]+pnt_config->Psi*(xImpuls_zeta_help[6]*pnt_U_RK->u[ijkPlus1]+yImpuls_zeta_help[6]*pnt_U_RK->v[ijkPlus1]+zImpuls_zeta_help[6]*pnt_U_RK->w[ijkPlus1]))+
							pnt_config->ZD_Ableitung_Koeffizient[7]*(Energie_zeta_help[7]*pnt_config->Gamma[ijkPlus2]+pnt_config->Psi*(xImpuls_zeta_help[7]*pnt_U_RK->u[ijkPlus2]+yImpuls_zeta_help[7]*pnt_U_RK->v[ijkPlus2]+zImpuls_zeta_help[7]*pnt_U_RK->w[ijkPlus2]))+
							pnt_config->ZD_Ableitung_Koeffizient[8]*(Energie_zeta_help[8]*pnt_config->Gamma[ijkPlus3]+pnt_config->Psi*(xImpuls_zeta_help[8]*pnt_U_RK->u[ijkPlus3]+yImpuls_zeta_help[8]*pnt_U_RK->v[ijkPlus3]+zImpuls_zeta_help[8]*pnt_U_RK->w[ijkPlus3]))+
							pnt_config->ZD_Ableitung_Koeffizient[9]*(Energie_zeta_help[9]*pnt_config->Gamma[ijkPlus4]+pnt_config->Psi*(xImpuls_zeta_help[9]*pnt_U_RK->u[ijkPlus4]+yImpuls_zeta_help[9]*pnt_U_RK->v[ijkPlus4]+zImpuls_zeta_help[9]*pnt_U_RK->w[ijkPlus4]))+
							pnt_config->ZD_Ableitung_Koeffizient[10]*(Energie_zeta_help[10]*pnt_config->Gamma[ijkPlus5]+pnt_config->Psi*(xImpuls_zeta_help[10]*pnt_U_RK->u[ijkPlus5]+yImpuls_zeta_help[10]*pnt_U_RK->v[ijkPlus5]+zImpuls_zeta_help[10]*pnt_U_RK->w[ijkPlus5]));
#endif
#if SPACEORDER==5
					xImpuls_zeta_help[0]=(pnt_mesh->jacobian[ijkMinus3]*(pnt_mesh->zeta_x[ijkMinus3]*pnt_ZD->tau_xx[ijkMinus3]+pnt_mesh->zeta_y[ijkMinus3]*pnt_ZD->tau_xy[ijkMinus3]+pnt_mesh->zeta_z[ijkMinus3]*pnt_ZD->tau_xz[ijkMinus3]));
					xImpuls_zeta_help[1]=(pnt_mesh->jacobian[ijkMinus2]*(pnt_mesh->zeta_x[ijkMinus2]*pnt_ZD->tau_xx[ijkMinus2]+pnt_mesh->zeta_y[ijkMinus2]*pnt_ZD->tau_xy[ijkMinus2]+pnt_mesh->zeta_z[ijkMinus2]*pnt_ZD->tau_xz[ijkMinus2]));
					xImpuls_zeta_help[2]=(pnt_mesh->jacobian[ijkMinus1]*(pnt_mesh->zeta_x[ijkMinus1]*pnt_ZD->tau_xx[ijkMinus1]+pnt_mesh->zeta_y[ijkMinus1]*pnt_ZD->tau_xy[ijkMinus1]+pnt_mesh->zeta_z[ijkMinus1]*pnt_ZD->tau_xz[ijkMinus1]));
					xImpuls_zeta_help[4]=(pnt_mesh->jacobian[ijkPlus1]*(pnt_mesh->zeta_x[ijkPlus1]*pnt_ZD->tau_xx[ijkPlus1]+pnt_mesh->zeta_y[ijkPlus1]*pnt_ZD->tau_xy[ijkPlus1]+pnt_mesh->zeta_z[ijkPlus1]*pnt_ZD->tau_xz[ijkPlus1]));
					xImpuls_zeta_help[5]=(pnt_mesh->jacobian[ijkPlus2]*(pnt_mesh->zeta_x[ijkPlus2]*pnt_ZD->tau_xx[ijkPlus2]+pnt_mesh->zeta_y[ijkPlus2]*pnt_ZD->tau_xy[ijkPlus2]+pnt_mesh->zeta_z[ijkPlus2]*pnt_ZD->tau_xz[ijkPlus2]));
					xImpuls_zeta_help[6]=(pnt_mesh->jacobian[ijkPlus3]*(pnt_mesh->zeta_x[ijkPlus3]*pnt_ZD->tau_xx[ijkPlus3]+pnt_mesh->zeta_y[ijkPlus3]*pnt_ZD->tau_xy[ijkPlus3]+pnt_mesh->zeta_z[ijkPlus3]*pnt_ZD->tau_xz[ijkPlus3]));

					xImpuls_zeta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*xImpuls_zeta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*xImpuls_zeta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*xImpuls_zeta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[4]*xImpuls_zeta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[5]*xImpuls_zeta_help[5]+pnt_config->ZD_Ableitung_Koeffizient[6]*xImpuls_zeta_help[6]);


					yImpuls_zeta_help[0]=(pnt_mesh->jacobian[ijkMinus3]*(pnt_mesh->zeta_x[ijkMinus3]*pnt_ZD->tau_xy[ijkMinus3]+pnt_mesh->zeta_y[ijkMinus3]*pnt_ZD->tau_yy[ijkMinus3]+pnt_mesh->zeta_z[ijkMinus3]*pnt_ZD->tau_yz[ijkMinus3]));
					yImpuls_zeta_help[1]=(pnt_mesh->jacobian[ijkMinus2]*(pnt_mesh->zeta_x[ijkMinus2]*pnt_ZD->tau_xy[ijkMinus2]+pnt_mesh->zeta_y[ijkMinus2]*pnt_ZD->tau_yy[ijkMinus2]+pnt_mesh->zeta_z[ijkMinus2]*pnt_ZD->tau_yz[ijkMinus2]));
					yImpuls_zeta_help[2]=(pnt_mesh->jacobian[ijkMinus1]*(pnt_mesh->zeta_x[ijkMinus1]*pnt_ZD->tau_xy[ijkMinus1]+pnt_mesh->zeta_y[ijkMinus1]*pnt_ZD->tau_yy[ijkMinus1]+pnt_mesh->zeta_z[ijkMinus1]*pnt_ZD->tau_yz[ijkMinus1]));
					yImpuls_zeta_help[4]=(pnt_mesh->jacobian[ijkPlus1]*(pnt_mesh->zeta_x[ijkPlus1]*pnt_ZD->tau_xy[ijkPlus1]+pnt_mesh->zeta_y[ijkPlus1]*pnt_ZD->tau_yy[ijkPlus1]+pnt_mesh->zeta_z[ijkPlus1]*pnt_ZD->tau_yz[ijkPlus1]));
					yImpuls_zeta_help[5]=(pnt_mesh->jacobian[ijkPlus2]*(pnt_mesh->zeta_x[ijkPlus2]*pnt_ZD->tau_xy[ijkPlus2]+pnt_mesh->zeta_y[ijkPlus2]*pnt_ZD->tau_yy[ijkPlus2]+pnt_mesh->zeta_z[ijkPlus2]*pnt_ZD->tau_yz[ijkPlus2]));
					yImpuls_zeta_help[6]=(pnt_mesh->jacobian[ijkPlus3]*(pnt_mesh->zeta_x[ijkPlus3]*pnt_ZD->tau_xy[ijkPlus3]+pnt_mesh->zeta_y[ijkPlus3]*pnt_ZD->tau_yy[ijkPlus3]+pnt_mesh->zeta_z[ijkPlus3]*pnt_ZD->tau_yz[ijkPlus3]));

					yImpuls_zeta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*yImpuls_zeta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*yImpuls_zeta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*yImpuls_zeta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[4]*yImpuls_zeta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[5]*yImpuls_zeta_help[5]+pnt_config->ZD_Ableitung_Koeffizient[6]*yImpuls_zeta_help[6]);


#if MESHDIMENSIONS==3
						zImpuls_zeta_help[0]=(pnt_mesh->jacobian[ijkMinus3]*(pnt_mesh->zeta_x[ijkMinus3]*pnt_ZD->tau_xz[ijkMinus3]+pnt_mesh->zeta_y[ijkMinus3]*pnt_ZD->tau_yz[ijkMinus3]+pnt_mesh->zeta_z[ijkMinus3]*pnt_ZD->tau_zz[ijkMinus3]));
						zImpuls_zeta_help[1]=(pnt_mesh->jacobian[ijkMinus2]*(pnt_mesh->zeta_x[ijkMinus2]*pnt_ZD->tau_xz[ijkMinus2]+pnt_mesh->zeta_y[ijkMinus2]*pnt_ZD->tau_yz[ijkMinus2]+pnt_mesh->zeta_z[ijkMinus2]*pnt_ZD->tau_zz[ijkMinus2]));
						zImpuls_zeta_help[2]=(pnt_mesh->jacobian[ijkMinus1]*(pnt_mesh->zeta_x[ijkMinus1]*pnt_ZD->tau_xz[ijkMinus1]+pnt_mesh->zeta_y[ijkMinus1]*pnt_ZD->tau_yz[ijkMinus1]+pnt_mesh->zeta_z[ijkMinus1]*pnt_ZD->tau_zz[ijkMinus1]));
						zImpuls_zeta_help[4]=(pnt_mesh->jacobian[ijkPlus1]*(pnt_mesh->zeta_x[ijkPlus1]*pnt_ZD->tau_xz[ijkPlus1]+pnt_mesh->zeta_y[ijkPlus1]*pnt_ZD->tau_yz[ijkPlus1]+pnt_mesh->zeta_z[ijkPlus1]*pnt_ZD->tau_zz[ijkPlus1]));
						zImpuls_zeta_help[5]=(pnt_mesh->jacobian[ijkPlus2]*(pnt_mesh->zeta_x[ijkPlus2]*pnt_ZD->tau_xz[ijkPlus2]+pnt_mesh->zeta_y[ijkPlus2]*pnt_ZD->tau_yz[ijkPlus2]+pnt_mesh->zeta_z[ijkPlus2]*pnt_ZD->tau_zz[ijkPlus2]));
						zImpuls_zeta_help[6]=(pnt_mesh->jacobian[ijkPlus3]*(pnt_mesh->zeta_x[ijkPlus3]*pnt_ZD->tau_xz[ijkPlus3]+pnt_mesh->zeta_y[ijkPlus3]*pnt_ZD->tau_yz[ijkPlus3]+pnt_mesh->zeta_z[ijkPlus3]*pnt_ZD->tau_zz[ijkPlus3]));

						zImpuls_zeta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*zImpuls_zeta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*zImpuls_zeta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*zImpuls_zeta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[4]*zImpuls_zeta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[5]*zImpuls_zeta_help[5]+pnt_config->ZD_Ableitung_Koeffizient[6]*zImpuls_zeta_help[6]);
#endif


					Energie_zeta_help[0]=(pnt_mesh->jacobian[ijkMinus3]*(pnt_mesh->zeta_x[ijkMinus3]*pnt_ZD->q_x[ijkMinus3]+pnt_mesh->zeta_y[ijkMinus3]*pnt_ZD->q_y[ijkMinus3]+pnt_mesh->zeta_z[ijkMinus3]*pnt_ZD->q_z[ijkMinus3]));
					Energie_zeta_help[1]=(pnt_mesh->jacobian[ijkMinus2]*(pnt_mesh->zeta_x[ijkMinus2]*pnt_ZD->q_x[ijkMinus2]+pnt_mesh->zeta_y[ijkMinus2]*pnt_ZD->q_y[ijkMinus2]+pnt_mesh->zeta_z[ijkMinus2]*pnt_ZD->q_z[ijkMinus2]));
					Energie_zeta_help[2]=(pnt_mesh->jacobian[ijkMinus1]*(pnt_mesh->zeta_x[ijkMinus1]*pnt_ZD->q_x[ijkMinus1]+pnt_mesh->zeta_y[ijkMinus1]*pnt_ZD->q_y[ijkMinus1]+pnt_mesh->zeta_z[ijkMinus1]*pnt_ZD->q_z[ijkMinus1]));
					Energie_zeta_help[4]=(pnt_mesh->jacobian[ijkPlus1]*(pnt_mesh->zeta_x[ijkPlus1]*pnt_ZD->q_x[ijkPlus1]+pnt_mesh->zeta_y[ijkPlus1]*pnt_ZD->q_y[ijkPlus1]+pnt_mesh->zeta_z[ijkPlus1]*pnt_ZD->q_z[ijkPlus1]));
					Energie_zeta_help[5]=(pnt_mesh->jacobian[ijkPlus2]*(pnt_mesh->zeta_x[ijkPlus2]*pnt_ZD->q_x[ijkPlus2]+pnt_mesh->zeta_y[ijkPlus2]*pnt_ZD->q_y[ijkPlus2]+pnt_mesh->zeta_z[ijkPlus2]*pnt_ZD->q_z[ijkPlus2]));
					Energie_zeta_help[6]=(pnt_mesh->jacobian[ijkPlus3]*(pnt_mesh->zeta_x[ijkPlus3]*pnt_ZD->q_x[ijkPlus3]+pnt_mesh->zeta_y[ijkPlus3]*pnt_ZD->q_y[ijkPlus3]+pnt_mesh->zeta_z[ijkPlus3]*pnt_ZD->q_z[ijkPlus3]));

					Energie_zeta=
							pnt_config->ZD_Ableitung_Koeffizient[0]*(Energie_zeta_help[0]*pnt_config->Gamma[ijkMinus3]+pnt_config->Psi*(xImpuls_zeta_help[0]*pnt_U_RK->u[ijkMinus3]+yImpuls_zeta_help[0]*pnt_U_RK->v[ijkMinus3]+zImpuls_zeta_help[0]*pnt_U_RK->w[ijkMinus3]))+
							pnt_config->ZD_Ableitung_Koeffizient[1]*(Energie_zeta_help[1]*pnt_config->Gamma[ijkMinus2]+pnt_config->Psi*(xImpuls_zeta_help[1]*pnt_U_RK->u[ijkMinus2]+yImpuls_zeta_help[1]*pnt_U_RK->v[ijkMinus2]+zImpuls_zeta_help[1]*pnt_U_RK->w[ijkMinus2]))+
							pnt_config->ZD_Ableitung_Koeffizient[2]*(Energie_zeta_help[2]*pnt_config->Gamma[ijkMinus1]+pnt_config->Psi*(xImpuls_zeta_help[2]*pnt_U_RK->u[ijkMinus1]+yImpuls_zeta_help[2]*pnt_U_RK->v[ijkMinus1]+zImpuls_zeta_help[2]*pnt_U_RK->w[ijkMinus1]))+
							pnt_config->ZD_Ableitung_Koeffizient[4]*(Energie_zeta_help[4]*pnt_config->Gamma[ijkPlus1]+pnt_config->Psi*(xImpuls_zeta_help[4]*pnt_U_RK->u[ijkPlus1]+yImpuls_zeta_help[4]*pnt_U_RK->v[ijkPlus1]+zImpuls_zeta_help[4]*pnt_U_RK->w[ijkPlus1]))+
							pnt_config->ZD_Ableitung_Koeffizient[5]*(Energie_zeta_help[5]*pnt_config->Gamma[ijkPlus2]+pnt_config->Psi*(xImpuls_zeta_help[5]*pnt_U_RK->u[ijkPlus2]+yImpuls_zeta_help[5]*pnt_U_RK->v[ijkPlus2]+zImpuls_zeta_help[5]*pnt_U_RK->w[ijkPlus2]))+
							pnt_config->ZD_Ableitung_Koeffizient[6]*(Energie_zeta_help[6]*pnt_config->Gamma[ijkPlus3]+pnt_config->Psi*(xImpuls_zeta_help[6]*pnt_U_RK->u[ijkPlus3]+yImpuls_zeta_help[6]*pnt_U_RK->v[ijkPlus3]+zImpuls_zeta_help[6]*pnt_U_RK->w[ijkPlus3]));
#endif

//				Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
				pnt_Q->xiMomentum[ijk]=pnt_Q->xiMomentum[ijk]+xImpuls_zeta;
				pnt_Q->etaMomentum[ijk]=pnt_Q->etaMomentum[ijk]+yImpuls_zeta;
#if MESHDIMENSIONS==3
					pnt_Q->zetaMomentum[ijk]=pnt_Q->zetaMomentum[ijk]+zImpuls_zeta;
#endif
#if MESHDIMENSIONS==2
					pnt_Q->zetaMomentum[ijk]=0.0;
#endif
				pnt_Q->Energy[ijk]=pnt_Q->Energy[ijk]+Energie_zeta;
			}
		}
	}

}


void CalcViscidFluxesInEtaDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD,
		struct strct_Flux * pnt_Q)
{
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iPlus1jk=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iPlus2jk=(i+2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iPlus3jk=(i+3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iMinus2jk=(i-2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				iMinus3jk=(i-3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;
				ijPlus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+2)*pnt_config->int_kMeshPointsGhostCells+k;
				ijPlus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+3)*pnt_config->int_kMeshPointsGhostCells+k;
				ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;
				ijMinus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-2)*pnt_config->int_kMeshPointsGhostCells+k;
				ijMinus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-3)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijkPlus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+1;
//				ijkPlus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+2;
//				ijkPlus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+3;
//				ijkMinus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-1;
//				ijkMinus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-2;
//				ijkMinus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-3;

#if SPACEORDER==9
//					iPlus4jk=(i+4)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//					iPlus5jk=(i+5)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//					iMinus4jk=(i-4)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//					iMinus5jk=(i-5)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					ijPlus4k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+4)*pnt_config->int_kMeshPointsGhostCells+k;
					ijPlus5k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+5)*pnt_config->int_kMeshPointsGhostCells+k;
					ijMinus4k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-4)*pnt_config->int_kMeshPointsGhostCells+k;
					ijMinus5k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-5)*pnt_config->int_kMeshPointsGhostCells+k;
//					ijkPlus4=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+4;
//					ijkPlus5=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+5;
//					ijkMinus4=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-4;
//					ijkMinus5=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-5;
#endif


#if SPACEORDER==9
					xImpuls_eta_help[0]=(pnt_mesh->jacobian[ijMinus5k]*(pnt_mesh->eta_x[ijMinus5k]*pnt_ZD->tau_xx[ijMinus5k]+pnt_mesh->eta_y[ijMinus5k]*pnt_ZD->tau_xy[ijMinus5k]+pnt_mesh->eta_z[ijMinus5k]*pnt_ZD->tau_xz[ijMinus5k]));
					xImpuls_eta_help[1]=(pnt_mesh->jacobian[ijMinus4k]*(pnt_mesh->eta_x[ijMinus4k]*pnt_ZD->tau_xx[ijMinus4k]+pnt_mesh->eta_y[ijMinus4k]*pnt_ZD->tau_xy[ijMinus4k]+pnt_mesh->eta_z[ijMinus4k]*pnt_ZD->tau_xz[ijMinus4k]));
					xImpuls_eta_help[2]=(pnt_mesh->jacobian[ijMinus3k]*(pnt_mesh->eta_x[ijMinus3k]*pnt_ZD->tau_xx[ijMinus3k]+pnt_mesh->eta_y[ijMinus3k]*pnt_ZD->tau_xy[ijMinus3k]+pnt_mesh->eta_z[ijMinus3k]*pnt_ZD->tau_xz[ijMinus3k]));
					xImpuls_eta_help[3]=(pnt_mesh->jacobian[ijMinus2k]*(pnt_mesh->eta_x[ijMinus2k]*pnt_ZD->tau_xx[ijMinus2k]+pnt_mesh->eta_y[ijMinus2k]*pnt_ZD->tau_xy[ijMinus2k]+pnt_mesh->eta_z[ijMinus2k]*pnt_ZD->tau_xz[ijMinus2k]));
					xImpuls_eta_help[4]=(pnt_mesh->jacobian[ijMinus1k]*(pnt_mesh->eta_x[ijMinus1k]*pnt_ZD->tau_xx[ijMinus1k]+pnt_mesh->eta_y[ijMinus1k]*pnt_ZD->tau_xy[ijMinus1k]+pnt_mesh->eta_z[ijMinus1k]*pnt_ZD->tau_xz[ijMinus1k]));
					xImpuls_eta_help[6]=(pnt_mesh->jacobian[ijPlus1k]*(pnt_mesh->eta_x[ijPlus1k]*pnt_ZD->tau_xx[ijPlus1k]+pnt_mesh->eta_y[ijPlus1k]*pnt_ZD->tau_xy[ijPlus1k]+pnt_mesh->eta_z[ijPlus1k]*pnt_ZD->tau_xz[ijPlus1k]));
					xImpuls_eta_help[7]=(pnt_mesh->jacobian[ijPlus2k]*(pnt_mesh->eta_x[ijPlus2k]*pnt_ZD->tau_xx[ijPlus2k]+pnt_mesh->eta_y[ijPlus2k]*pnt_ZD->tau_xy[ijPlus2k]+pnt_mesh->eta_z[ijPlus2k]*pnt_ZD->tau_xz[ijPlus2k]));
					xImpuls_eta_help[8]=(pnt_mesh->jacobian[ijPlus3k]*(pnt_mesh->eta_x[ijPlus3k]*pnt_ZD->tau_xx[ijPlus3k]+pnt_mesh->eta_y[ijPlus3k]*pnt_ZD->tau_xy[ijPlus3k]+pnt_mesh->eta_z[ijPlus3k]*pnt_ZD->tau_xz[ijPlus3k]));
					xImpuls_eta_help[9]=(pnt_mesh->jacobian[ijPlus4k]*(pnt_mesh->eta_x[ijPlus4k]*pnt_ZD->tau_xx[ijPlus4k]+pnt_mesh->eta_y[ijPlus4k]*pnt_ZD->tau_xy[ijPlus4k]+pnt_mesh->eta_z[ijPlus4k]*pnt_ZD->tau_xz[ijPlus4k]));
					xImpuls_eta_help[10]=(pnt_mesh->jacobian[ijPlus5k]*(pnt_mesh->eta_x[ijPlus5k]*pnt_ZD->tau_xx[ijPlus5k]+pnt_mesh->eta_y[ijPlus5k]*pnt_ZD->tau_xy[ijPlus5k]+pnt_mesh->eta_z[ijPlus5k]*pnt_ZD->tau_xz[ijPlus5k]));

					xImpuls_eta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*xImpuls_eta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*xImpuls_eta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*xImpuls_eta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[3]*xImpuls_eta_help[3]+pnt_config->ZD_Ableitung_Koeffizient[4]*xImpuls_eta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[6]*xImpuls_eta_help[6]+pnt_config->ZD_Ableitung_Koeffizient[7]*xImpuls_eta_help[7]+pnt_config->ZD_Ableitung_Koeffizient[8]*xImpuls_eta_help[8]+pnt_config->ZD_Ableitung_Koeffizient[9]*xImpuls_eta_help[9]+pnt_config->ZD_Ableitung_Koeffizient[10]*xImpuls_eta_help[10]);


					yImpuls_eta_help[0]=(pnt_mesh->jacobian[ijMinus5k]*(pnt_mesh->eta_x[ijMinus5k]*pnt_ZD->tau_xy[ijMinus5k]+pnt_mesh->eta_y[ijMinus5k]*pnt_ZD->tau_yy[ijMinus5k]+pnt_mesh->eta_z[ijMinus5k]*pnt_ZD->tau_yz[ijMinus5k]));
					yImpuls_eta_help[1]=(pnt_mesh->jacobian[ijMinus4k]*(pnt_mesh->eta_x[ijMinus4k]*pnt_ZD->tau_xy[ijMinus4k]+pnt_mesh->eta_y[ijMinus4k]*pnt_ZD->tau_yy[ijMinus4k]+pnt_mesh->eta_z[ijMinus4k]*pnt_ZD->tau_yz[ijMinus4k]));
					yImpuls_eta_help[2]=(pnt_mesh->jacobian[ijMinus3k]*(pnt_mesh->eta_x[ijMinus3k]*pnt_ZD->tau_xy[ijMinus3k]+pnt_mesh->eta_y[ijMinus3k]*pnt_ZD->tau_yy[ijMinus3k]+pnt_mesh->eta_z[ijMinus3k]*pnt_ZD->tau_yz[ijMinus3k]));
					yImpuls_eta_help[3]=(pnt_mesh->jacobian[ijMinus2k]*(pnt_mesh->eta_x[ijMinus2k]*pnt_ZD->tau_xy[ijMinus2k]+pnt_mesh->eta_y[ijMinus2k]*pnt_ZD->tau_yy[ijMinus2k]+pnt_mesh->eta_z[ijMinus2k]*pnt_ZD->tau_yz[ijMinus2k]));
					yImpuls_eta_help[4]=(pnt_mesh->jacobian[ijMinus1k]*(pnt_mesh->eta_x[ijMinus1k]*pnt_ZD->tau_xy[ijMinus1k]+pnt_mesh->eta_y[ijMinus1k]*pnt_ZD->tau_yy[ijMinus1k]+pnt_mesh->eta_z[ijMinus1k]*pnt_ZD->tau_yz[ijMinus1k]));
					yImpuls_eta_help[6]=(pnt_mesh->jacobian[ijPlus1k]*(pnt_mesh->eta_x[ijPlus1k]*pnt_ZD->tau_xy[ijPlus1k]+pnt_mesh->eta_y[ijPlus1k]*pnt_ZD->tau_yy[ijPlus1k]+pnt_mesh->eta_z[ijPlus1k]*pnt_ZD->tau_yz[ijPlus1k]));
					yImpuls_eta_help[7]=(pnt_mesh->jacobian[ijPlus2k]*(pnt_mesh->eta_x[ijPlus2k]*pnt_ZD->tau_xy[ijPlus2k]+pnt_mesh->eta_y[ijPlus2k]*pnt_ZD->tau_yy[ijPlus2k]+pnt_mesh->eta_z[ijPlus2k]*pnt_ZD->tau_yz[ijPlus2k]));
					yImpuls_eta_help[8]=(pnt_mesh->jacobian[ijPlus3k]*(pnt_mesh->eta_x[ijPlus3k]*pnt_ZD->tau_xy[ijPlus3k]+pnt_mesh->eta_y[ijPlus3k]*pnt_ZD->tau_yy[ijPlus3k]+pnt_mesh->eta_z[ijPlus3k]*pnt_ZD->tau_yz[ijPlus3k]));
					yImpuls_eta_help[9]=(pnt_mesh->jacobian[ijPlus4k]*(pnt_mesh->eta_x[ijPlus4k]*pnt_ZD->tau_xy[ijPlus4k]+pnt_mesh->eta_y[ijPlus4k]*pnt_ZD->tau_yy[ijPlus4k]+pnt_mesh->eta_z[ijPlus4k]*pnt_ZD->tau_yz[ijPlus4k]));
					yImpuls_eta_help[10]=(pnt_mesh->jacobian[ijPlus5k]*(pnt_mesh->eta_x[ijPlus5k]*pnt_ZD->tau_xy[ijPlus5k]+pnt_mesh->eta_y[ijPlus5k]*pnt_ZD->tau_yy[ijPlus5k]+pnt_mesh->eta_z[ijPlus5k]*pnt_ZD->tau_yz[ijPlus5k]));

					yImpuls_eta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*yImpuls_eta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*yImpuls_eta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*yImpuls_eta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[3]*yImpuls_eta_help[3]+pnt_config->ZD_Ableitung_Koeffizient[4]*yImpuls_eta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[6]*yImpuls_eta_help[6]+pnt_config->ZD_Ableitung_Koeffizient[7]*yImpuls_eta_help[7]+pnt_config->ZD_Ableitung_Koeffizient[8]*yImpuls_eta_help[8]+pnt_config->ZD_Ableitung_Koeffizient[9]*yImpuls_eta_help[9]+pnt_config->ZD_Ableitung_Koeffizient[10]*yImpuls_eta_help[10]);


#if MESHDIMENSIONS==3
						zImpuls_eta_help[0]=(pnt_mesh->jacobian[ijMinus5k]*(pnt_mesh->eta_x[ijMinus5k]*pnt_ZD->tau_xz[ijMinus5k]+pnt_mesh->eta_y[ijMinus5k]*pnt_ZD->tau_yz[ijMinus5k]+pnt_mesh->eta_z[ijMinus5k]*pnt_ZD->tau_zz[ijMinus5k]));
						zImpuls_eta_help[1]=(pnt_mesh->jacobian[ijMinus4k]*(pnt_mesh->eta_x[ijMinus4k]*pnt_ZD->tau_xz[ijMinus4k]+pnt_mesh->eta_y[ijMinus4k]*pnt_ZD->tau_yz[ijMinus4k]+pnt_mesh->eta_z[ijMinus4k]*pnt_ZD->tau_zz[ijMinus4k]));
						zImpuls_eta_help[2]=(pnt_mesh->jacobian[ijMinus3k]*(pnt_mesh->eta_x[ijMinus3k]*pnt_ZD->tau_xz[ijMinus3k]+pnt_mesh->eta_y[ijMinus3k]*pnt_ZD->tau_yz[ijMinus3k]+pnt_mesh->eta_z[ijMinus3k]*pnt_ZD->tau_zz[ijMinus3k]));
						zImpuls_eta_help[3]=(pnt_mesh->jacobian[ijMinus2k]*(pnt_mesh->eta_x[ijMinus2k]*pnt_ZD->tau_xz[ijMinus2k]+pnt_mesh->eta_y[ijMinus2k]*pnt_ZD->tau_yz[ijMinus2k]+pnt_mesh->eta_z[ijMinus2k]*pnt_ZD->tau_zz[ijMinus2k]));
						zImpuls_eta_help[4]=(pnt_mesh->jacobian[ijMinus1k]*(pnt_mesh->eta_x[ijMinus1k]*pnt_ZD->tau_xz[ijMinus1k]+pnt_mesh->eta_y[ijMinus1k]*pnt_ZD->tau_yz[ijMinus1k]+pnt_mesh->eta_z[ijMinus1k]*pnt_ZD->tau_zz[ijMinus1k]));
						zImpuls_eta_help[6]=(pnt_mesh->jacobian[ijPlus1k]*(pnt_mesh->eta_x[ijPlus1k]*pnt_ZD->tau_xz[ijPlus1k]+pnt_mesh->eta_y[ijPlus1k]*pnt_ZD->tau_yz[ijPlus1k]+pnt_mesh->eta_z[ijPlus1k]*pnt_ZD->tau_zz[ijPlus1k]));
						zImpuls_eta_help[7]=(pnt_mesh->jacobian[ijPlus2k]*(pnt_mesh->eta_x[ijPlus2k]*pnt_ZD->tau_xz[ijPlus2k]+pnt_mesh->eta_y[ijPlus2k]*pnt_ZD->tau_yz[ijPlus2k]+pnt_mesh->eta_z[ijPlus2k]*pnt_ZD->tau_zz[ijPlus2k]));
						zImpuls_eta_help[8]=(pnt_mesh->jacobian[ijPlus3k]*(pnt_mesh->eta_x[ijPlus3k]*pnt_ZD->tau_xz[ijPlus3k]+pnt_mesh->eta_y[ijPlus3k]*pnt_ZD->tau_yz[ijPlus3k]+pnt_mesh->eta_z[ijPlus3k]*pnt_ZD->tau_zz[ijPlus3k]));
						zImpuls_eta_help[9]=(pnt_mesh->jacobian[ijPlus4k]*(pnt_mesh->eta_x[ijPlus4k]*pnt_ZD->tau_xz[ijPlus4k]+pnt_mesh->eta_y[ijPlus4k]*pnt_ZD->tau_yz[ijPlus4k]+pnt_mesh->eta_z[ijPlus4k]*pnt_ZD->tau_zz[ijPlus4k]));
						zImpuls_eta_help[10]=(pnt_mesh->jacobian[ijPlus5k]*(pnt_mesh->eta_x[ijPlus5k]*pnt_ZD->tau_xz[ijPlus5k]+pnt_mesh->eta_y[ijPlus5k]*pnt_ZD->tau_yz[ijPlus5k]+pnt_mesh->eta_z[ijPlus5k]*pnt_ZD->tau_zz[ijPlus5k]));

						zImpuls_eta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*zImpuls_eta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*zImpuls_eta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*zImpuls_eta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[3]*zImpuls_eta_help[3]+pnt_config->ZD_Ableitung_Koeffizient[4]*zImpuls_eta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[6]*zImpuls_eta_help[6]+pnt_config->ZD_Ableitung_Koeffizient[7]*zImpuls_eta_help[7]+pnt_config->ZD_Ableitung_Koeffizient[8]*zImpuls_eta_help[8]+pnt_config->ZD_Ableitung_Koeffizient[9]*zImpuls_eta_help[9]+pnt_config->ZD_Ableitung_Koeffizient[10]*zImpuls_eta_help[10]);
#endif


					Energie_eta_help[0]=(pnt_mesh->jacobian[ijMinus5k]*(pnt_mesh->eta_x[ijMinus5k]*pnt_ZD->q_x[ijMinus5k]+pnt_mesh->eta_y[ijMinus5k]*pnt_ZD->q_y[ijMinus5k]+pnt_mesh->eta_z[ijMinus5k]*pnt_ZD->q_z[ijMinus5k]));
					Energie_eta_help[1]=(pnt_mesh->jacobian[ijMinus4k]*(pnt_mesh->eta_x[ijMinus4k]*pnt_ZD->q_x[ijMinus4k]+pnt_mesh->eta_y[ijMinus4k]*pnt_ZD->q_y[ijMinus4k]+pnt_mesh->eta_z[ijMinus4k]*pnt_ZD->q_z[ijMinus4k]));
					Energie_eta_help[2]=(pnt_mesh->jacobian[ijMinus3k]*(pnt_mesh->eta_x[ijMinus3k]*pnt_ZD->q_x[ijMinus3k]+pnt_mesh->eta_y[ijMinus3k]*pnt_ZD->q_y[ijMinus3k]+pnt_mesh->eta_z[ijMinus3k]*pnt_ZD->q_z[ijMinus3k]));
					Energie_eta_help[3]=(pnt_mesh->jacobian[ijMinus2k]*(pnt_mesh->eta_x[ijMinus2k]*pnt_ZD->q_x[ijMinus2k]+pnt_mesh->eta_y[ijMinus2k]*pnt_ZD->q_y[ijMinus2k]+pnt_mesh->eta_z[ijMinus2k]*pnt_ZD->q_z[ijMinus2k]));
					Energie_eta_help[4]=(pnt_mesh->jacobian[ijMinus1k]*(pnt_mesh->eta_x[ijMinus1k]*pnt_ZD->q_x[ijMinus1k]+pnt_mesh->eta_y[ijMinus1k]*pnt_ZD->q_y[ijMinus1k]+pnt_mesh->eta_z[ijMinus1k]*pnt_ZD->q_z[ijMinus1k]));
					Energie_eta_help[6]=(pnt_mesh->jacobian[ijPlus1k]*(pnt_mesh->eta_x[ijPlus1k]*pnt_ZD->q_x[ijPlus1k]+pnt_mesh->eta_y[ijPlus1k]*pnt_ZD->q_y[ijPlus1k]+pnt_mesh->eta_z[ijPlus1k]*pnt_ZD->q_z[ijPlus1k]));
					Energie_eta_help[7]=(pnt_mesh->jacobian[ijPlus2k]*(pnt_mesh->eta_x[ijPlus2k]*pnt_ZD->q_x[ijPlus2k]+pnt_mesh->eta_y[ijPlus2k]*pnt_ZD->q_y[ijPlus2k]+pnt_mesh->eta_z[ijPlus2k]*pnt_ZD->q_z[ijPlus2k]));
					Energie_eta_help[8]=(pnt_mesh->jacobian[ijPlus3k]*(pnt_mesh->eta_x[ijPlus3k]*pnt_ZD->q_x[ijPlus3k]+pnt_mesh->eta_y[ijPlus3k]*pnt_ZD->q_y[ijPlus3k]+pnt_mesh->eta_z[ijPlus3k]*pnt_ZD->q_z[ijPlus3k]));
					Energie_eta_help[9]=(pnt_mesh->jacobian[ijPlus4k]*(pnt_mesh->eta_x[ijPlus4k]*pnt_ZD->q_x[ijPlus4k]+pnt_mesh->eta_y[ijPlus4k]*pnt_ZD->q_y[ijPlus4k]+pnt_mesh->eta_z[ijPlus4k]*pnt_ZD->q_z[ijPlus4k]));
					Energie_eta_help[10]=(pnt_mesh->jacobian[ijPlus5k]*(pnt_mesh->eta_x[ijPlus5k]*pnt_ZD->q_x[ijPlus5k]+pnt_mesh->eta_y[ijPlus5k]*pnt_ZD->q_y[ijPlus5k]+pnt_mesh->eta_z[ijPlus5k]*pnt_ZD->q_z[ijPlus5k]));

					Energie_eta=
							pnt_config->ZD_Ableitung_Koeffizient[0]*(Energie_eta_help[0]*pnt_config->Gamma[ijMinus5k]+pnt_config->Psi*(xImpuls_eta_help[0]*pnt_U_RK->u[ijMinus5k]+yImpuls_eta_help[0]*pnt_U_RK->v[ijMinus5k]+zImpuls_eta_help[0]*pnt_U_RK->w[ijMinus5k]))+
							pnt_config->ZD_Ableitung_Koeffizient[1]*(Energie_eta_help[1]*pnt_config->Gamma[ijMinus4k]+pnt_config->Psi*(xImpuls_eta_help[1]*pnt_U_RK->u[ijMinus4k]+yImpuls_eta_help[1]*pnt_U_RK->v[ijMinus4k]+zImpuls_eta_help[1]*pnt_U_RK->w[ijMinus4k]))+
							pnt_config->ZD_Ableitung_Koeffizient[2]*(Energie_eta_help[2]*pnt_config->Gamma[ijMinus3k]+pnt_config->Psi*(xImpuls_eta_help[2]*pnt_U_RK->u[ijMinus3k]+yImpuls_eta_help[2]*pnt_U_RK->v[ijMinus3k]+zImpuls_eta_help[2]*pnt_U_RK->w[ijMinus3k]))+
							pnt_config->ZD_Ableitung_Koeffizient[3]*(Energie_eta_help[3]*pnt_config->Gamma[ijMinus2k]+pnt_config->Psi*(xImpuls_eta_help[3]*pnt_U_RK->u[ijMinus2k]+yImpuls_eta_help[3]*pnt_U_RK->v[ijMinus2k]+zImpuls_eta_help[3]*pnt_U_RK->w[ijMinus2k]))+
							pnt_config->ZD_Ableitung_Koeffizient[4]*(Energie_eta_help[4]*pnt_config->Gamma[ijMinus1k]+pnt_config->Psi*(xImpuls_eta_help[4]*pnt_U_RK->u[ijMinus1k]+yImpuls_eta_help[4]*pnt_U_RK->v[ijMinus1k]+zImpuls_eta_help[4]*pnt_U_RK->w[ijMinus1k]))+
							pnt_config->ZD_Ableitung_Koeffizient[6]*(Energie_eta_help[6]*pnt_config->Gamma[ijPlus1k]+pnt_config->Psi*(xImpuls_eta_help[6]*pnt_U_RK->u[ijPlus1k]+yImpuls_eta_help[6]*pnt_U_RK->v[ijPlus1k]+zImpuls_eta_help[6]*pnt_U_RK->w[ijPlus1k]))+
							pnt_config->ZD_Ableitung_Koeffizient[7]*(Energie_eta_help[7]*pnt_config->Gamma[ijPlus2k]+pnt_config->Psi*(xImpuls_eta_help[7]*pnt_U_RK->u[ijPlus2k]+yImpuls_eta_help[7]*pnt_U_RK->v[ijPlus2k]+zImpuls_eta_help[7]*pnt_U_RK->w[ijPlus2k]))+
							pnt_config->ZD_Ableitung_Koeffizient[8]*(Energie_eta_help[8]*pnt_config->Gamma[ijPlus3k]+pnt_config->Psi*(xImpuls_eta_help[8]*pnt_U_RK->u[ijPlus3k]+yImpuls_eta_help[8]*pnt_U_RK->v[ijPlus3k]+zImpuls_eta_help[8]*pnt_U_RK->w[ijPlus3k]))+
							pnt_config->ZD_Ableitung_Koeffizient[9]*(Energie_eta_help[9]*pnt_config->Gamma[ijPlus4k]+pnt_config->Psi*(xImpuls_eta_help[9]*pnt_U_RK->u[ijPlus4k]+yImpuls_eta_help[9]*pnt_U_RK->v[ijPlus4k]+zImpuls_eta_help[9]*pnt_U_RK->w[ijPlus4k]))+
							pnt_config->ZD_Ableitung_Koeffizient[10]*(Energie_eta_help[10]*pnt_config->Gamma[ijPlus5k]+pnt_config->Psi*(xImpuls_eta_help[10]*pnt_U_RK->u[ijPlus5k]+yImpuls_eta_help[10]*pnt_U_RK->v[ijPlus5k]+zImpuls_eta_help[10]*pnt_U_RK->w[ijPlus5k]));
#endif
#if SPACEORDER==5
					xImpuls_eta_help[0]=(pnt_mesh->jacobian[ijMinus3k]*(pnt_mesh->eta_x[ijMinus3k]*pnt_ZD->tau_xx[ijMinus3k]+pnt_mesh->eta_y[ijMinus3k]*pnt_ZD->tau_xy[ijMinus3k]+pnt_mesh->eta_z[ijMinus3k]*pnt_ZD->tau_xz[ijMinus3k]));
					xImpuls_eta_help[1]=(pnt_mesh->jacobian[ijMinus2k]*(pnt_mesh->eta_x[ijMinus2k]*pnt_ZD->tau_xx[ijMinus2k]+pnt_mesh->eta_y[ijMinus2k]*pnt_ZD->tau_xy[ijMinus2k]+pnt_mesh->eta_z[ijMinus2k]*pnt_ZD->tau_xz[ijMinus2k]));
					xImpuls_eta_help[2]=(pnt_mesh->jacobian[ijMinus1k]*(pnt_mesh->eta_x[ijMinus1k]*pnt_ZD->tau_xx[ijMinus1k]+pnt_mesh->eta_y[ijMinus1k]*pnt_ZD->tau_xy[ijMinus1k]+pnt_mesh->eta_z[ijMinus1k]*pnt_ZD->tau_xz[ijMinus1k]));
					xImpuls_eta_help[4]=(pnt_mesh->jacobian[ijPlus1k]*(pnt_mesh->eta_x[ijPlus1k]*pnt_ZD->tau_xx[ijPlus1k]+pnt_mesh->eta_y[ijPlus1k]*pnt_ZD->tau_xy[ijPlus1k]+pnt_mesh->eta_z[ijPlus1k]*pnt_ZD->tau_xz[ijPlus1k]));
					xImpuls_eta_help[5]=(pnt_mesh->jacobian[ijPlus2k]*(pnt_mesh->eta_x[ijPlus2k]*pnt_ZD->tau_xx[ijPlus2k]+pnt_mesh->eta_y[ijPlus2k]*pnt_ZD->tau_xy[ijPlus2k]+pnt_mesh->eta_z[ijPlus2k]*pnt_ZD->tau_xz[ijPlus2k]));
					xImpuls_eta_help[6]=(pnt_mesh->jacobian[ijPlus3k]*(pnt_mesh->eta_x[ijPlus3k]*pnt_ZD->tau_xx[ijPlus3k]+pnt_mesh->eta_y[ijPlus3k]*pnt_ZD->tau_xy[ijPlus3k]+pnt_mesh->eta_z[ijPlus3k]*pnt_ZD->tau_xz[ijPlus3k]));

					xImpuls_eta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*xImpuls_eta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*xImpuls_eta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*xImpuls_eta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[4]*xImpuls_eta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[5]*xImpuls_eta_help[5]+pnt_config->ZD_Ableitung_Koeffizient[6]*xImpuls_eta_help[6]);


					yImpuls_eta_help[0]=(pnt_mesh->jacobian[ijMinus3k]*(pnt_mesh->eta_x[ijMinus3k]*pnt_ZD->tau_xy[ijMinus3k]+pnt_mesh->eta_y[ijMinus3k]*pnt_ZD->tau_yy[ijMinus3k]+pnt_mesh->eta_z[ijMinus3k]*pnt_ZD->tau_yz[ijMinus3k]));
					yImpuls_eta_help[1]=(pnt_mesh->jacobian[ijMinus2k]*(pnt_mesh->eta_x[ijMinus2k]*pnt_ZD->tau_xy[ijMinus2k]+pnt_mesh->eta_y[ijMinus2k]*pnt_ZD->tau_yy[ijMinus2k]+pnt_mesh->eta_z[ijMinus2k]*pnt_ZD->tau_yz[ijMinus2k]));
					yImpuls_eta_help[2]=(pnt_mesh->jacobian[ijMinus1k]*(pnt_mesh->eta_x[ijMinus1k]*pnt_ZD->tau_xy[ijMinus1k]+pnt_mesh->eta_y[ijMinus1k]*pnt_ZD->tau_yy[ijMinus1k]+pnt_mesh->eta_z[ijMinus1k]*pnt_ZD->tau_yz[ijMinus1k]));
					yImpuls_eta_help[4]=(pnt_mesh->jacobian[ijPlus1k]*(pnt_mesh->eta_x[ijPlus1k]*pnt_ZD->tau_xy[ijPlus1k]+pnt_mesh->eta_y[ijPlus1k]*pnt_ZD->tau_yy[ijPlus1k]+pnt_mesh->eta_z[ijPlus1k]*pnt_ZD->tau_yz[ijPlus1k]));
					yImpuls_eta_help[5]=(pnt_mesh->jacobian[ijPlus2k]*(pnt_mesh->eta_x[ijPlus2k]*pnt_ZD->tau_xy[ijPlus2k]+pnt_mesh->eta_y[ijPlus2k]*pnt_ZD->tau_yy[ijPlus2k]+pnt_mesh->eta_z[ijPlus2k]*pnt_ZD->tau_yz[ijPlus2k]));
					yImpuls_eta_help[6]=(pnt_mesh->jacobian[ijPlus3k]*(pnt_mesh->eta_x[ijPlus3k]*pnt_ZD->tau_xy[ijPlus3k]+pnt_mesh->eta_y[ijPlus3k]*pnt_ZD->tau_yy[ijPlus3k]+pnt_mesh->eta_z[ijPlus3k]*pnt_ZD->tau_yz[ijPlus3k]));

					yImpuls_eta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*yImpuls_eta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*yImpuls_eta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*yImpuls_eta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[4]*yImpuls_eta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[5]*yImpuls_eta_help[5]+pnt_config->ZD_Ableitung_Koeffizient[6]*yImpuls_eta_help[6]);


#if MESHDIMENSIONS==3
						zImpuls_eta_help[0]=(pnt_mesh->jacobian[ijMinus3k]*(pnt_mesh->eta_x[ijMinus3k]*pnt_ZD->tau_xz[ijMinus3k]+pnt_mesh->eta_y[ijMinus3k]*pnt_ZD->tau_yz[ijMinus3k]+pnt_mesh->eta_z[ijMinus3k]*pnt_ZD->tau_zz[ijMinus3k]));
						zImpuls_eta_help[1]=(pnt_mesh->jacobian[ijMinus2k]*(pnt_mesh->eta_x[ijMinus2k]*pnt_ZD->tau_xz[ijMinus2k]+pnt_mesh->eta_y[ijMinus2k]*pnt_ZD->tau_yz[ijMinus2k]+pnt_mesh->eta_z[ijMinus2k]*pnt_ZD->tau_zz[ijMinus2k]));
						zImpuls_eta_help[2]=(pnt_mesh->jacobian[ijMinus1k]*(pnt_mesh->eta_x[ijMinus1k]*pnt_ZD->tau_xz[ijMinus1k]+pnt_mesh->eta_y[ijMinus1k]*pnt_ZD->tau_yz[ijMinus1k]+pnt_mesh->eta_z[ijMinus1k]*pnt_ZD->tau_zz[ijMinus1k]));
						zImpuls_eta_help[4]=(pnt_mesh->jacobian[ijPlus1k]*(pnt_mesh->eta_x[ijPlus1k]*pnt_ZD->tau_xz[ijPlus1k]+pnt_mesh->eta_y[ijPlus1k]*pnt_ZD->tau_yz[ijPlus1k]+pnt_mesh->eta_z[ijPlus1k]*pnt_ZD->tau_zz[ijPlus1k]));
						zImpuls_eta_help[5]=(pnt_mesh->jacobian[ijPlus2k]*(pnt_mesh->eta_x[ijPlus2k]*pnt_ZD->tau_xz[ijPlus2k]+pnt_mesh->eta_y[ijPlus2k]*pnt_ZD->tau_yz[ijPlus2k]+pnt_mesh->eta_z[ijPlus2k]*pnt_ZD->tau_zz[ijPlus2k]));
						zImpuls_eta_help[6]=(pnt_mesh->jacobian[ijPlus3k]*(pnt_mesh->eta_x[ijPlus3k]*pnt_ZD->tau_xz[ijPlus3k]+pnt_mesh->eta_y[ijPlus3k]*pnt_ZD->tau_yz[ijPlus3k]+pnt_mesh->eta_z[ijPlus3k]*pnt_ZD->tau_zz[ijPlus3k]));

						zImpuls_eta=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*zImpuls_eta_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*zImpuls_eta_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*zImpuls_eta_help[2]+pnt_config->ZD_Ableitung_Koeffizient[4]*zImpuls_eta_help[4]+pnt_config->ZD_Ableitung_Koeffizient[5]*zImpuls_eta_help[5]+pnt_config->ZD_Ableitung_Koeffizient[6]*zImpuls_eta_help[6]);
#endif


					Energie_eta_help[0]=(pnt_mesh->jacobian[ijMinus3k]*(pnt_mesh->eta_x[ijMinus3k]*pnt_ZD->q_x[ijMinus3k]+pnt_mesh->eta_y[ijMinus3k]*pnt_ZD->q_y[ijMinus3k]+pnt_mesh->eta_z[ijMinus3k]*pnt_ZD->q_z[ijMinus3k]));
					Energie_eta_help[1]=(pnt_mesh->jacobian[ijMinus2k]*(pnt_mesh->eta_x[ijMinus2k]*pnt_ZD->q_x[ijMinus2k]+pnt_mesh->eta_y[ijMinus2k]*pnt_ZD->q_y[ijMinus2k]+pnt_mesh->eta_z[ijMinus2k]*pnt_ZD->q_z[ijMinus2k]));
					Energie_eta_help[2]=(pnt_mesh->jacobian[ijMinus1k]*(pnt_mesh->eta_x[ijMinus1k]*pnt_ZD->q_x[ijMinus1k]+pnt_mesh->eta_y[ijMinus1k]*pnt_ZD->q_y[ijMinus1k]+pnt_mesh->eta_z[ijMinus1k]*pnt_ZD->q_z[ijMinus1k]));
					Energie_eta_help[4]=(pnt_mesh->jacobian[ijPlus1k]*(pnt_mesh->eta_x[ijPlus1k]*pnt_ZD->q_x[ijPlus1k]+pnt_mesh->eta_y[ijPlus1k]*pnt_ZD->q_y[ijPlus1k]+pnt_mesh->eta_z[ijPlus1k]*pnt_ZD->q_z[ijPlus1k]));
					Energie_eta_help[5]=(pnt_mesh->jacobian[ijPlus2k]*(pnt_mesh->eta_x[ijPlus2k]*pnt_ZD->q_x[ijPlus2k]+pnt_mesh->eta_y[ijPlus2k]*pnt_ZD->q_y[ijPlus2k]+pnt_mesh->eta_z[ijPlus2k]*pnt_ZD->q_z[ijPlus2k]));
					Energie_eta_help[6]=(pnt_mesh->jacobian[ijPlus3k]*(pnt_mesh->eta_x[ijPlus3k]*pnt_ZD->q_x[ijPlus3k]+pnt_mesh->eta_y[ijPlus3k]*pnt_ZD->q_y[ijPlus3k]+pnt_mesh->eta_z[ijPlus3k]*pnt_ZD->q_z[ijPlus3k]));

					Energie_eta=
							pnt_config->ZD_Ableitung_Koeffizient[0]*(Energie_eta_help[0]*pnt_config->Gamma[ijMinus3k]+pnt_config->Psi*(xImpuls_eta_help[0]*pnt_U_RK->u[ijMinus3k]+yImpuls_eta_help[0]*pnt_U_RK->v[ijMinus3k]+zImpuls_eta_help[0]*pnt_U_RK->w[ijMinus3k]))+
							pnt_config->ZD_Ableitung_Koeffizient[1]*(Energie_eta_help[1]*pnt_config->Gamma[ijMinus2k]+pnt_config->Psi*(xImpuls_eta_help[1]*pnt_U_RK->u[ijMinus2k]+yImpuls_eta_help[1]*pnt_U_RK->v[ijMinus2k]+zImpuls_eta_help[1]*pnt_U_RK->w[ijMinus2k]))+
							pnt_config->ZD_Ableitung_Koeffizient[2]*(Energie_eta_help[2]*pnt_config->Gamma[ijMinus1k]+pnt_config->Psi*(xImpuls_eta_help[2]*pnt_U_RK->u[ijMinus1k]+yImpuls_eta_help[2]*pnt_U_RK->v[ijMinus1k]+zImpuls_eta_help[2]*pnt_U_RK->w[ijMinus1k]))+
							pnt_config->ZD_Ableitung_Koeffizient[4]*(Energie_eta_help[4]*pnt_config->Gamma[ijPlus1k]+pnt_config->Psi*(xImpuls_eta_help[4]*pnt_U_RK->u[ijPlus1k]+yImpuls_eta_help[4]*pnt_U_RK->v[ijPlus1k]+zImpuls_eta_help[4]*pnt_U_RK->w[ijPlus1k]))+
							pnt_config->ZD_Ableitung_Koeffizient[5]*(Energie_eta_help[5]*pnt_config->Gamma[ijPlus2k]+pnt_config->Psi*(xImpuls_eta_help[5]*pnt_U_RK->u[ijPlus2k]+yImpuls_eta_help[5]*pnt_U_RK->v[ijPlus2k]+zImpuls_eta_help[5]*pnt_U_RK->w[ijPlus2k]))+
							pnt_config->ZD_Ableitung_Koeffizient[6]*(Energie_eta_help[6]*pnt_config->Gamma[ijPlus3k]+pnt_config->Psi*(xImpuls_eta_help[6]*pnt_U_RK->u[ijPlus3k]+yImpuls_eta_help[6]*pnt_U_RK->v[ijPlus3k]+zImpuls_eta_help[6]*pnt_U_RK->w[ijPlus3k]));
#endif


//					Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
				pnt_Q->xiMomentum[ijk]=pnt_Q->xiMomentum[ijk]+xImpuls_eta;
				pnt_Q->etaMomentum[ijk]=pnt_Q->etaMomentum[ijk]+yImpuls_eta;
#if MESHDIMENSIONS==3
					pnt_Q->zetaMomentum[ijk]=pnt_Q->zetaMomentum[ijk]+zImpuls_eta;
#endif
#if MESHDIMENSIONS==2
					pnt_Q->zetaMomentum[ijk]=0.0;
#endif
				pnt_Q->Energy[ijk]=pnt_Q->Energy[ijk]+Energie_eta;
			}
		}
	}

}

void CalcViscidFluxesInXiDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD,
		struct strct_Flux * pnt_Q)
{
	for (i=pnt_config->int_iStartReal; i <= pnt_config->int_iEndReal; i++)
	{
		for (j=pnt_config->int_jStartReal; j <= pnt_config->int_jEndReal; j++)
		{
			for (k=pnt_config->int_kStartReal; k <= pnt_config->int_kEndReal; k++)
			{
				ijk=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iPlus1jk=(i+1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iPlus2jk=(i+2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iPlus3jk=(i+3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iMinus1jk=(i-1)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iMinus2jk=(i-2)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
				iMinus3jk=(i-3)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//				ijPlus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+1)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijPlus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+2)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijPlus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+3)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijMinus1k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-1)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijMinus2k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-2)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijMinus3k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-3)*pnt_config->int_kMeshPointsGhostCells+k;
//				ijkPlus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+1;
//				ijkPlus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+2;
//				ijkPlus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+3;
//				ijkMinus1=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-1;
//				ijkMinus2=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-2;
//				ijkMinus3=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-3;

#if SPACEORDER==9
					iPlus4jk=(i+4)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					iPlus5jk=(i+5)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					iMinus4jk=(i-4)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
					iMinus5jk=(i-5)*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k;
//					ijPlus4k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+4)*pnt_config->int_kMeshPointsGhostCells+k;
//					ijPlus5k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j+5)*pnt_config->int_kMeshPointsGhostCells+k;
//					ijMinus4k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-4)*pnt_config->int_kMeshPointsGhostCells+k;
//					ijMinus5k=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+(j-5)*pnt_config->int_kMeshPointsGhostCells+k;
//					ijkPlus4=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+4;
//					ijkPlus5=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k+5;
//					ijkMinus4=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-4;
//					ijkMinus5=i*pnt_config->int_jMeshPointsGhostCells*pnt_config->int_kMeshPointsGhostCells+j*pnt_config->int_kMeshPointsGhostCells+k-5;
#endif


#if SPACEORDER==9
					xImpuls_xi_help[0]=(pnt_mesh->jacobian[iMinus5jk]*(pnt_mesh->xi_x[iMinus5jk]*pnt_ZD->tau_xx[iMinus5jk]+pnt_mesh->xi_y[iMinus5jk]*pnt_ZD->tau_xy[iMinus5jk]+pnt_mesh->xi_z[iMinus5jk]*pnt_ZD->tau_xz[iMinus5jk]));
					xImpuls_xi_help[1]=(pnt_mesh->jacobian[iMinus4jk]*(pnt_mesh->xi_x[iMinus4jk]*pnt_ZD->tau_xx[iMinus4jk]+pnt_mesh->xi_y[iMinus4jk]*pnt_ZD->tau_xy[iMinus4jk]+pnt_mesh->xi_z[iMinus4jk]*pnt_ZD->tau_xz[iMinus4jk]));
					xImpuls_xi_help[2]=(pnt_mesh->jacobian[iMinus3jk]*(pnt_mesh->xi_x[iMinus3jk]*pnt_ZD->tau_xx[iMinus3jk]+pnt_mesh->xi_y[iMinus3jk]*pnt_ZD->tau_xy[iMinus3jk]+pnt_mesh->xi_z[iMinus3jk]*pnt_ZD->tau_xz[iMinus3jk]));
					xImpuls_xi_help[3]=(pnt_mesh->jacobian[iMinus2jk]*(pnt_mesh->xi_x[iMinus2jk]*pnt_ZD->tau_xx[iMinus2jk]+pnt_mesh->xi_y[iMinus2jk]*pnt_ZD->tau_xy[iMinus2jk]+pnt_mesh->xi_z[iMinus2jk]*pnt_ZD->tau_xz[iMinus2jk]));
					xImpuls_xi_help[4]=(pnt_mesh->jacobian[iMinus1jk]*(pnt_mesh->xi_x[iMinus1jk]*pnt_ZD->tau_xx[iMinus1jk]+pnt_mesh->xi_y[iMinus1jk]*pnt_ZD->tau_xy[iMinus1jk]+pnt_mesh->xi_z[iMinus1jk]*pnt_ZD->tau_xz[iMinus1jk]));
					xImpuls_xi_help[6]=(pnt_mesh->jacobian[iPlus1jk]*(pnt_mesh->xi_x[iPlus1jk]*pnt_ZD->tau_xx[iPlus1jk]+pnt_mesh->xi_y[iPlus1jk]*pnt_ZD->tau_xy[iPlus1jk]+pnt_mesh->xi_z[iPlus1jk]*pnt_ZD->tau_xz[iPlus1jk]));
					xImpuls_xi_help[7]=(pnt_mesh->jacobian[iPlus2jk]*(pnt_mesh->xi_x[iPlus2jk]*pnt_ZD->tau_xx[iPlus2jk]+pnt_mesh->xi_y[iPlus2jk]*pnt_ZD->tau_xy[iPlus2jk]+pnt_mesh->xi_z[iPlus2jk]*pnt_ZD->tau_xz[iPlus2jk]));
					xImpuls_xi_help[8]=(pnt_mesh->jacobian[iPlus3jk]*(pnt_mesh->xi_x[iPlus3jk]*pnt_ZD->tau_xx[iPlus3jk]+pnt_mesh->xi_y[iPlus3jk]*pnt_ZD->tau_xy[iPlus3jk]+pnt_mesh->xi_z[iPlus3jk]*pnt_ZD->tau_xz[iPlus3jk]));
					xImpuls_xi_help[9]=(pnt_mesh->jacobian[iPlus4jk]*(pnt_mesh->xi_x[iPlus4jk]*pnt_ZD->tau_xx[iPlus4jk]+pnt_mesh->xi_y[iPlus4jk]*pnt_ZD->tau_xy[iPlus4jk]+pnt_mesh->xi_z[iPlus4jk]*pnt_ZD->tau_xz[iPlus4jk]));
					xImpuls_xi_help[10]=(pnt_mesh->jacobian[iPlus5jk]*(pnt_mesh->xi_x[iPlus5jk]*pnt_ZD->tau_xx[iPlus5jk]+pnt_mesh->xi_y[iPlus5jk]*pnt_ZD->tau_xy[iPlus5jk]+pnt_mesh->xi_z[iPlus5jk]*pnt_ZD->tau_xz[iPlus5jk]));

					xImpuls_xi=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*xImpuls_xi_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*xImpuls_xi_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*xImpuls_xi_help[2]+pnt_config->ZD_Ableitung_Koeffizient[3]*xImpuls_xi_help[3]+pnt_config->ZD_Ableitung_Koeffizient[4]*xImpuls_xi_help[4]+pnt_config->ZD_Ableitung_Koeffizient[6]*xImpuls_xi_help[6]+pnt_config->ZD_Ableitung_Koeffizient[7]*xImpuls_xi_help[7]+pnt_config->ZD_Ableitung_Koeffizient[8]*xImpuls_xi_help[8]+pnt_config->ZD_Ableitung_Koeffizient[9]*xImpuls_xi_help[9]+pnt_config->ZD_Ableitung_Koeffizient[10]*xImpuls_xi_help[10]);


					yImpuls_xi_help[0]=(pnt_mesh->jacobian[iMinus5jk]*(pnt_mesh->xi_x[iMinus5jk]*pnt_ZD->tau_xy[iMinus5jk]+pnt_mesh->xi_y[iMinus5jk]*pnt_ZD->tau_yy[iMinus5jk]+pnt_mesh->xi_z[iMinus5jk]*pnt_ZD->tau_yz[iMinus5jk]));
					yImpuls_xi_help[1]=(pnt_mesh->jacobian[iMinus4jk]*(pnt_mesh->xi_x[iMinus4jk]*pnt_ZD->tau_xy[iMinus4jk]+pnt_mesh->xi_y[iMinus4jk]*pnt_ZD->tau_yy[iMinus4jk]+pnt_mesh->xi_z[iMinus4jk]*pnt_ZD->tau_yz[iMinus4jk]));
					yImpuls_xi_help[2]=(pnt_mesh->jacobian[iMinus3jk]*(pnt_mesh->xi_x[iMinus3jk]*pnt_ZD->tau_xy[iMinus3jk]+pnt_mesh->xi_y[iMinus3jk]*pnt_ZD->tau_yy[iMinus3jk]+pnt_mesh->xi_z[iMinus3jk]*pnt_ZD->tau_yz[iMinus3jk]));
					yImpuls_xi_help[3]=(pnt_mesh->jacobian[iMinus2jk]*(pnt_mesh->xi_x[iMinus2jk]*pnt_ZD->tau_xy[iMinus2jk]+pnt_mesh->xi_y[iMinus2jk]*pnt_ZD->tau_yy[iMinus2jk]+pnt_mesh->xi_z[iMinus2jk]*pnt_ZD->tau_yz[iMinus2jk]));
					yImpuls_xi_help[4]=(pnt_mesh->jacobian[iMinus1jk]*(pnt_mesh->xi_x[iMinus1jk]*pnt_ZD->tau_xy[iMinus1jk]+pnt_mesh->xi_y[iMinus1jk]*pnt_ZD->tau_yy[iMinus1jk]+pnt_mesh->xi_z[iMinus1jk]*pnt_ZD->tau_yz[iMinus1jk]));
					yImpuls_xi_help[6]=(pnt_mesh->jacobian[iPlus1jk]*(pnt_mesh->xi_x[iPlus1jk]*pnt_ZD->tau_xy[iPlus1jk]+pnt_mesh->xi_y[iPlus1jk]*pnt_ZD->tau_yy[iPlus1jk]+pnt_mesh->xi_z[iPlus1jk]*pnt_ZD->tau_yz[iPlus1jk]));
					yImpuls_xi_help[7]=(pnt_mesh->jacobian[iPlus2jk]*(pnt_mesh->xi_x[iPlus2jk]*pnt_ZD->tau_xy[iPlus2jk]+pnt_mesh->xi_y[iPlus2jk]*pnt_ZD->tau_yy[iPlus2jk]+pnt_mesh->xi_z[iPlus2jk]*pnt_ZD->tau_yz[iPlus2jk]));
					yImpuls_xi_help[8]=(pnt_mesh->jacobian[iPlus3jk]*(pnt_mesh->xi_x[iPlus3jk]*pnt_ZD->tau_xy[iPlus3jk]+pnt_mesh->xi_y[iPlus3jk]*pnt_ZD->tau_yy[iPlus3jk]+pnt_mesh->xi_z[iPlus3jk]*pnt_ZD->tau_yz[iPlus3jk]));
					yImpuls_xi_help[9]=(pnt_mesh->jacobian[iPlus4jk]*(pnt_mesh->xi_x[iPlus4jk]*pnt_ZD->tau_xy[iPlus4jk]+pnt_mesh->xi_y[iPlus4jk]*pnt_ZD->tau_yy[iPlus4jk]+pnt_mesh->xi_z[iPlus4jk]*pnt_ZD->tau_yz[iPlus4jk]));
					yImpuls_xi_help[10]=(pnt_mesh->jacobian[iPlus5jk]*(pnt_mesh->xi_x[iPlus5jk]*pnt_ZD->tau_xy[iPlus5jk]+pnt_mesh->xi_y[iPlus5jk]*pnt_ZD->tau_yy[iPlus5jk]+pnt_mesh->xi_z[iPlus5jk]*pnt_ZD->tau_yz[iPlus5jk]));

					yImpuls_xi=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*yImpuls_xi_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*yImpuls_xi_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*yImpuls_xi_help[2]+pnt_config->ZD_Ableitung_Koeffizient[3]*yImpuls_xi_help[3]+pnt_config->ZD_Ableitung_Koeffizient[4]*yImpuls_xi_help[4]+pnt_config->ZD_Ableitung_Koeffizient[6]*yImpuls_xi_help[6]+pnt_config->ZD_Ableitung_Koeffizient[7]*yImpuls_xi_help[7]+pnt_config->ZD_Ableitung_Koeffizient[8]*yImpuls_xi_help[8]+pnt_config->ZD_Ableitung_Koeffizient[9]*yImpuls_xi_help[9]+pnt_config->ZD_Ableitung_Koeffizient[10]*yImpuls_xi_help[10]);


#if MESHDIMENSIONS==3
						zImpuls_xi_help[0]=(pnt_mesh->jacobian[iMinus5jk]*(pnt_mesh->xi_x[iMinus5jk]*pnt_ZD->tau_xz[iMinus5jk]+pnt_mesh->xi_y[iMinus5jk]*pnt_ZD->tau_yz[iMinus5jk]+pnt_mesh->xi_z[iMinus5jk]*pnt_ZD->tau_zz[iMinus5jk]));
						zImpuls_xi_help[1]=(pnt_mesh->jacobian[iMinus4jk]*(pnt_mesh->xi_x[iMinus4jk]*pnt_ZD->tau_xz[iMinus4jk]+pnt_mesh->xi_y[iMinus4jk]*pnt_ZD->tau_yz[iMinus4jk]+pnt_mesh->xi_z[iMinus4jk]*pnt_ZD->tau_zz[iMinus4jk]));
						zImpuls_xi_help[2]=(pnt_mesh->jacobian[iMinus3jk]*(pnt_mesh->xi_x[iMinus3jk]*pnt_ZD->tau_xz[iMinus3jk]+pnt_mesh->xi_y[iMinus3jk]*pnt_ZD->tau_yz[iMinus3jk]+pnt_mesh->xi_z[iMinus3jk]*pnt_ZD->tau_zz[iMinus3jk]));
						zImpuls_xi_help[3]=(pnt_mesh->jacobian[iMinus2jk]*(pnt_mesh->xi_x[iMinus2jk]*pnt_ZD->tau_xz[iMinus2jk]+pnt_mesh->xi_y[iMinus2jk]*pnt_ZD->tau_yz[iMinus2jk]+pnt_mesh->xi_z[iMinus2jk]*pnt_ZD->tau_zz[iMinus2jk]));
						zImpuls_xi_help[4]=(pnt_mesh->jacobian[iMinus1jk]*(pnt_mesh->xi_x[iMinus1jk]*pnt_ZD->tau_xz[iMinus1jk]+pnt_mesh->xi_y[iMinus1jk]*pnt_ZD->tau_yz[iMinus1jk]+pnt_mesh->xi_z[iMinus1jk]*pnt_ZD->tau_zz[iMinus1jk]));
						zImpuls_xi_help[6]=(pnt_mesh->jacobian[iPlus1jk]*(pnt_mesh->xi_x[iPlus1jk]*pnt_ZD->tau_xz[iPlus1jk]+pnt_mesh->xi_y[iPlus1jk]*pnt_ZD->tau_yz[iPlus1jk]+pnt_mesh->xi_z[iPlus1jk]*pnt_ZD->tau_zz[iPlus1jk]));
						zImpuls_xi_help[7]=(pnt_mesh->jacobian[iPlus2jk]*(pnt_mesh->xi_x[iPlus2jk]*pnt_ZD->tau_xz[iPlus2jk]+pnt_mesh->xi_y[iPlus2jk]*pnt_ZD->tau_yz[iPlus2jk]+pnt_mesh->xi_z[iPlus2jk]*pnt_ZD->tau_zz[iPlus2jk]));
						zImpuls_xi_help[8]=(pnt_mesh->jacobian[iPlus3jk]*(pnt_mesh->xi_x[iPlus3jk]*pnt_ZD->tau_xz[iPlus3jk]+pnt_mesh->xi_y[iPlus3jk]*pnt_ZD->tau_yz[iPlus3jk]+pnt_mesh->xi_z[iPlus3jk]*pnt_ZD->tau_zz[iPlus3jk]));
						zImpuls_xi_help[9]=(pnt_mesh->jacobian[iPlus4jk]*(pnt_mesh->xi_x[iPlus4jk]*pnt_ZD->tau_xz[iPlus4jk]+pnt_mesh->xi_y[iPlus4jk]*pnt_ZD->tau_yz[iPlus4jk]+pnt_mesh->xi_z[iPlus4jk]*pnt_ZD->tau_zz[iPlus4jk]));
						zImpuls_xi_help[10]=(pnt_mesh->jacobian[iPlus5jk]*(pnt_mesh->xi_x[iPlus5jk]*pnt_ZD->tau_xz[iPlus5jk]+pnt_mesh->xi_y[iPlus5jk]*pnt_ZD->tau_yz[iPlus5jk]+pnt_mesh->xi_z[iPlus5jk]*pnt_ZD->tau_zz[iPlus5jk]));

						zImpuls_xi=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*zImpuls_xi_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*zImpuls_xi_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*zImpuls_xi_help[2]+pnt_config->ZD_Ableitung_Koeffizient[3]*zImpuls_xi_help[3]+pnt_config->ZD_Ableitung_Koeffizient[4]*zImpuls_xi_help[4]+pnt_config->ZD_Ableitung_Koeffizient[6]*zImpuls_xi_help[6]+pnt_config->ZD_Ableitung_Koeffizient[7]*zImpuls_xi_help[7]+pnt_config->ZD_Ableitung_Koeffizient[8]*zImpuls_xi_help[8]+pnt_config->ZD_Ableitung_Koeffizient[9]*zImpuls_xi_help[9]+pnt_config->ZD_Ableitung_Koeffizient[10]*zImpuls_xi_help[10]);
#endif


					Energie_xi_help[0]=(pnt_mesh->jacobian[iMinus5jk]*(pnt_mesh->xi_x[iMinus5jk]*pnt_ZD->q_x[iMinus5jk]+pnt_mesh->xi_y[iMinus5jk]*pnt_ZD->q_y[iMinus5jk]+pnt_mesh->xi_z[iMinus5jk]*pnt_ZD->q_z[iMinus5jk]));
					Energie_xi_help[1]=(pnt_mesh->jacobian[iMinus4jk]*(pnt_mesh->xi_x[iMinus4jk]*pnt_ZD->q_x[iMinus4jk]+pnt_mesh->xi_y[iMinus4jk]*pnt_ZD->q_y[iMinus4jk]+pnt_mesh->xi_z[iMinus4jk]*pnt_ZD->q_z[iMinus4jk]));
					Energie_xi_help[2]=(pnt_mesh->jacobian[iMinus3jk]*(pnt_mesh->xi_x[iMinus3jk]*pnt_ZD->q_x[iMinus3jk]+pnt_mesh->xi_y[iMinus3jk]*pnt_ZD->q_y[iMinus3jk]+pnt_mesh->xi_z[iMinus3jk]*pnt_ZD->q_z[iMinus3jk]));
					Energie_xi_help[3]=(pnt_mesh->jacobian[iMinus2jk]*(pnt_mesh->xi_x[iMinus2jk]*pnt_ZD->q_x[iMinus2jk]+pnt_mesh->xi_y[iMinus2jk]*pnt_ZD->q_y[iMinus2jk]+pnt_mesh->xi_z[iMinus2jk]*pnt_ZD->q_z[iMinus2jk]));
					Energie_xi_help[4]=(pnt_mesh->jacobian[iMinus1jk]*(pnt_mesh->xi_x[iMinus1jk]*pnt_ZD->q_x[iMinus1jk]+pnt_mesh->xi_y[iMinus1jk]*pnt_ZD->q_y[iMinus1jk]+pnt_mesh->xi_z[iMinus1jk]*pnt_ZD->q_z[iMinus1jk]));
					Energie_xi_help[6]=(pnt_mesh->jacobian[iPlus1jk]*(pnt_mesh->xi_x[iPlus1jk]*pnt_ZD->q_x[iPlus1jk]+pnt_mesh->xi_y[iPlus1jk]*pnt_ZD->q_y[iPlus1jk]+pnt_mesh->xi_z[iPlus1jk]*pnt_ZD->q_z[iPlus1jk]));
					Energie_xi_help[7]=(pnt_mesh->jacobian[iPlus2jk]*(pnt_mesh->xi_x[iPlus2jk]*pnt_ZD->q_x[iPlus2jk]+pnt_mesh->xi_y[iPlus2jk]*pnt_ZD->q_y[iPlus2jk]+pnt_mesh->xi_z[iPlus2jk]*pnt_ZD->q_z[iPlus2jk]));
					Energie_xi_help[8]=(pnt_mesh->jacobian[iPlus3jk]*(pnt_mesh->xi_x[iPlus3jk]*pnt_ZD->q_x[iPlus3jk]+pnt_mesh->xi_y[iPlus3jk]*pnt_ZD->q_y[iPlus3jk]+pnt_mesh->xi_z[iPlus3jk]*pnt_ZD->q_z[iPlus3jk]));
					Energie_xi_help[9]=(pnt_mesh->jacobian[iPlus4jk]*(pnt_mesh->xi_x[iPlus4jk]*pnt_ZD->q_x[iPlus4jk]+pnt_mesh->xi_y[iPlus4jk]*pnt_ZD->q_y[iPlus4jk]+pnt_mesh->xi_z[iPlus4jk]*pnt_ZD->q_z[iPlus4jk]));
					Energie_xi_help[10]=(pnt_mesh->jacobian[iPlus5jk]*(pnt_mesh->xi_x[iPlus5jk]*pnt_ZD->q_x[iPlus5jk]+pnt_mesh->xi_y[iPlus5jk]*pnt_ZD->q_y[iPlus5jk]+pnt_mesh->xi_z[iPlus5jk]*pnt_ZD->q_z[iPlus5jk]));

					Energie_xi=
							pnt_config->ZD_Ableitung_Koeffizient[0]*(Energie_xi_help[0]*pnt_config->Gamma[iMinus5jk]+pnt_config->Psi*(xImpuls_xi_help[0]*pnt_U_RK->u[iMinus5jk]+yImpuls_xi_help[0]*pnt_U_RK->v[iMinus5jk]+zImpuls_xi_help[0]*pnt_U_RK->w[iMinus5jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[1]*(Energie_xi_help[1]*pnt_config->Gamma[iMinus4jk]+pnt_config->Psi*(xImpuls_xi_help[1]*pnt_U_RK->u[iMinus4jk]+yImpuls_xi_help[1]*pnt_U_RK->v[iMinus4jk]+zImpuls_xi_help[1]*pnt_U_RK->w[iMinus4jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[2]*(Energie_xi_help[2]*pnt_config->Gamma[iMinus3jk]+pnt_config->Psi*(xImpuls_xi_help[2]*pnt_U_RK->u[iMinus3jk]+yImpuls_xi_help[2]*pnt_U_RK->v[iMinus3jk]+zImpuls_xi_help[2]*pnt_U_RK->w[iMinus3jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[3]*(Energie_xi_help[3]*pnt_config->Gamma[iMinus2jk]+pnt_config->Psi*(xImpuls_xi_help[3]*pnt_U_RK->u[iMinus2jk]+yImpuls_xi_help[3]*pnt_U_RK->v[iMinus2jk]+zImpuls_xi_help[3]*pnt_U_RK->w[iMinus2jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[4]*(Energie_xi_help[4]*pnt_config->Gamma[iMinus1jk]+pnt_config->Psi*(xImpuls_xi_help[4]*pnt_U_RK->u[iMinus1jk]+yImpuls_xi_help[4]*pnt_U_RK->v[iMinus1jk]+zImpuls_xi_help[4]*pnt_U_RK->w[iMinus1jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[6]*(Energie_xi_help[6]*pnt_config->Gamma[iPlus1jk]+pnt_config->Psi*(xImpuls_xi_help[6]*pnt_U_RK->u[iPlus1jk]+yImpuls_xi_help[6]*pnt_U_RK->v[iPlus1jk]+zImpuls_xi_help[6]*pnt_U_RK->w[iPlus1jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[7]*(Energie_xi_help[7]*pnt_config->Gamma[iPlus2jk]+pnt_config->Psi*(xImpuls_xi_help[7]*pnt_U_RK->u[iPlus2jk]+yImpuls_xi_help[7]*pnt_U_RK->v[iPlus2jk]+zImpuls_xi_help[7]*pnt_U_RK->w[iPlus2jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[8]*(Energie_xi_help[8]*pnt_config->Gamma[iPlus3jk]+pnt_config->Psi*(xImpuls_xi_help[8]*pnt_U_RK->u[iPlus3jk]+yImpuls_xi_help[8]*pnt_U_RK->v[iPlus3jk]+zImpuls_xi_help[8]*pnt_U_RK->w[iPlus3jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[9]*(Energie_xi_help[9]*pnt_config->Gamma[iPlus4jk]+pnt_config->Psi*(xImpuls_xi_help[9]*pnt_U_RK->u[iPlus4jk]+yImpuls_xi_help[9]*pnt_U_RK->v[iPlus4jk]+zImpuls_xi_help[9]*pnt_U_RK->w[iPlus4jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[10]*(Energie_xi_help[10]*pnt_config->Gamma[iPlus5jk]+pnt_config->Psi*(xImpuls_xi_help[10]*pnt_U_RK->u[iPlus5jk]+yImpuls_xi_help[10]*pnt_U_RK->v[iPlus5jk]+zImpuls_xi_help[10]*pnt_U_RK->w[iPlus5jk]));
#endif
#if SPACEORDER==5
					xImpuls_xi_help[0]=(pnt_mesh->jacobian[iMinus3jk]*(pnt_mesh->xi_x[iMinus3jk]*pnt_ZD->tau_xx[iMinus3jk]+pnt_mesh->xi_y[iMinus3jk]*pnt_ZD->tau_xy[iMinus3jk]+pnt_mesh->xi_z[iMinus3jk]*pnt_ZD->tau_xz[iMinus3jk]));
					xImpuls_xi_help[1]=(pnt_mesh->jacobian[iMinus2jk]*(pnt_mesh->xi_x[iMinus2jk]*pnt_ZD->tau_xx[iMinus2jk]+pnt_mesh->xi_y[iMinus2jk]*pnt_ZD->tau_xy[iMinus2jk]+pnt_mesh->xi_z[iMinus2jk]*pnt_ZD->tau_xz[iMinus2jk]));
					xImpuls_xi_help[2]=(pnt_mesh->jacobian[iMinus1jk]*(pnt_mesh->xi_x[iMinus1jk]*pnt_ZD->tau_xx[iMinus1jk]+pnt_mesh->xi_y[iMinus1jk]*pnt_ZD->tau_xy[iMinus1jk]+pnt_mesh->xi_z[iMinus1jk]*pnt_ZD->tau_xz[iMinus1jk]));
					xImpuls_xi_help[4]=(pnt_mesh->jacobian[iPlus1jk]*(pnt_mesh->xi_x[iPlus1jk]*pnt_ZD->tau_xx[iPlus1jk]+pnt_mesh->xi_y[iPlus1jk]*pnt_ZD->tau_xy[iPlus1jk]+pnt_mesh->xi_z[iPlus1jk]*pnt_ZD->tau_xz[iPlus1jk]));
					xImpuls_xi_help[5]=(pnt_mesh->jacobian[iPlus2jk]*(pnt_mesh->xi_x[iPlus2jk]*pnt_ZD->tau_xx[iPlus2jk]+pnt_mesh->xi_y[iPlus2jk]*pnt_ZD->tau_xy[iPlus2jk]+pnt_mesh->xi_z[iPlus2jk]*pnt_ZD->tau_xz[iPlus2jk]));
					xImpuls_xi_help[6]=(pnt_mesh->jacobian[iPlus3jk]*(pnt_mesh->xi_x[iPlus3jk]*pnt_ZD->tau_xx[iPlus3jk]+pnt_mesh->xi_y[iPlus3jk]*pnt_ZD->tau_xy[iPlus3jk]+pnt_mesh->xi_z[iPlus3jk]*pnt_ZD->tau_xz[iPlus3jk]));

					xImpuls_xi=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*xImpuls_xi_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*xImpuls_xi_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*xImpuls_xi_help[2]+pnt_config->ZD_Ableitung_Koeffizient[4]*xImpuls_xi_help[4]+pnt_config->ZD_Ableitung_Koeffizient[5]*xImpuls_xi_help[5]+pnt_config->ZD_Ableitung_Koeffizient[6]*xImpuls_xi_help[6]);


					yImpuls_xi_help[0]=(pnt_mesh->jacobian[iMinus3jk]*(pnt_mesh->xi_x[iMinus3jk]*pnt_ZD->tau_xy[iMinus3jk]+pnt_mesh->xi_y[iMinus3jk]*pnt_ZD->tau_yy[iMinus3jk]+pnt_mesh->xi_z[iMinus3jk]*pnt_ZD->tau_yz[iMinus3jk]));
					yImpuls_xi_help[1]=(pnt_mesh->jacobian[iMinus2jk]*(pnt_mesh->xi_x[iMinus2jk]*pnt_ZD->tau_xy[iMinus2jk]+pnt_mesh->xi_y[iMinus2jk]*pnt_ZD->tau_yy[iMinus2jk]+pnt_mesh->xi_z[iMinus2jk]*pnt_ZD->tau_yz[iMinus2jk]));
					yImpuls_xi_help[2]=(pnt_mesh->jacobian[iMinus1jk]*(pnt_mesh->xi_x[iMinus1jk]*pnt_ZD->tau_xy[iMinus1jk]+pnt_mesh->xi_y[iMinus1jk]*pnt_ZD->tau_yy[iMinus1jk]+pnt_mesh->xi_z[iMinus1jk]*pnt_ZD->tau_yz[iMinus1jk]));
					yImpuls_xi_help[4]=(pnt_mesh->jacobian[iPlus1jk]*(pnt_mesh->xi_x[iPlus1jk]*pnt_ZD->tau_xy[iPlus1jk]+pnt_mesh->xi_y[iPlus1jk]*pnt_ZD->tau_yy[iPlus1jk]+pnt_mesh->xi_z[iPlus1jk]*pnt_ZD->tau_yz[iPlus1jk]));
					yImpuls_xi_help[5]=(pnt_mesh->jacobian[iPlus2jk]*(pnt_mesh->xi_x[iPlus2jk]*pnt_ZD->tau_xy[iPlus2jk]+pnt_mesh->xi_y[iPlus2jk]*pnt_ZD->tau_yy[iPlus2jk]+pnt_mesh->xi_z[iPlus2jk]*pnt_ZD->tau_yz[iPlus2jk]));
					yImpuls_xi_help[6]=(pnt_mesh->jacobian[iPlus3jk]*(pnt_mesh->xi_x[iPlus3jk]*pnt_ZD->tau_xy[iPlus3jk]+pnt_mesh->xi_y[iPlus3jk]*pnt_ZD->tau_yy[iPlus3jk]+pnt_mesh->xi_z[iPlus3jk]*pnt_ZD->tau_yz[iPlus3jk]));

					yImpuls_xi=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*yImpuls_xi_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*yImpuls_xi_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*yImpuls_xi_help[2]+pnt_config->ZD_Ableitung_Koeffizient[4]*yImpuls_xi_help[4]+pnt_config->ZD_Ableitung_Koeffizient[5]*yImpuls_xi_help[5]+pnt_config->ZD_Ableitung_Koeffizient[6]*yImpuls_xi_help[6]);


#if MESHDIMENSIONS==3
						zImpuls_xi_help[0]=(pnt_mesh->jacobian[iMinus3jk]*(pnt_mesh->xi_x[iMinus3jk]*pnt_ZD->tau_xz[iMinus3jk]+pnt_mesh->xi_y[iMinus3jk]*pnt_ZD->tau_yz[iMinus3jk]+pnt_mesh->xi_z[iMinus3jk]*pnt_ZD->tau_zz[iMinus3jk]));
						zImpuls_xi_help[1]=(pnt_mesh->jacobian[iMinus2jk]*(pnt_mesh->xi_x[iMinus2jk]*pnt_ZD->tau_xz[iMinus2jk]+pnt_mesh->xi_y[iMinus2jk]*pnt_ZD->tau_yz[iMinus2jk]+pnt_mesh->xi_z[iMinus2jk]*pnt_ZD->tau_zz[iMinus2jk]));
						zImpuls_xi_help[2]=(pnt_mesh->jacobian[iMinus1jk]*(pnt_mesh->xi_x[iMinus1jk]*pnt_ZD->tau_xz[iMinus1jk]+pnt_mesh->xi_y[iMinus1jk]*pnt_ZD->tau_yz[iMinus1jk]+pnt_mesh->xi_z[iMinus1jk]*pnt_ZD->tau_zz[iMinus1jk]));
						zImpuls_xi_help[4]=(pnt_mesh->jacobian[iPlus1jk]*(pnt_mesh->xi_x[iPlus1jk]*pnt_ZD->tau_xz[iPlus1jk]+pnt_mesh->xi_y[iPlus1jk]*pnt_ZD->tau_yz[iPlus1jk]+pnt_mesh->xi_z[iPlus1jk]*pnt_ZD->tau_zz[iPlus1jk]));
						zImpuls_xi_help[5]=(pnt_mesh->jacobian[iPlus2jk]*(pnt_mesh->xi_x[iPlus2jk]*pnt_ZD->tau_xz[iPlus2jk]+pnt_mesh->xi_y[iPlus2jk]*pnt_ZD->tau_yz[iPlus2jk]+pnt_mesh->xi_z[iPlus2jk]*pnt_ZD->tau_zz[iPlus2jk]));
						zImpuls_xi_help[6]=(pnt_mesh->jacobian[iPlus3jk]*(pnt_mesh->xi_x[iPlus3jk]*pnt_ZD->tau_xz[iPlus3jk]+pnt_mesh->xi_y[iPlus3jk]*pnt_ZD->tau_yz[iPlus3jk]+pnt_mesh->xi_z[iPlus3jk]*pnt_ZD->tau_zz[iPlus3jk]));

						zImpuls_xi=pnt_config->Psi*(pnt_config->ZD_Ableitung_Koeffizient[0]*zImpuls_xi_help[0]+pnt_config->ZD_Ableitung_Koeffizient[1]*zImpuls_xi_help[1]+pnt_config->ZD_Ableitung_Koeffizient[2]*zImpuls_xi_help[2]+pnt_config->ZD_Ableitung_Koeffizient[4]*zImpuls_xi_help[4]+pnt_config->ZD_Ableitung_Koeffizient[5]*zImpuls_xi_help[5]+pnt_config->ZD_Ableitung_Koeffizient[6]*zImpuls_xi_help[6]);
#endif


					Energie_xi_help[0]=(pnt_mesh->jacobian[iMinus3jk]*(pnt_mesh->xi_x[iMinus3jk]*pnt_ZD->q_x[iMinus3jk]+pnt_mesh->xi_y[iMinus3jk]*pnt_ZD->q_y[iMinus3jk]+pnt_mesh->xi_z[iMinus3jk]*pnt_ZD->q_z[iMinus3jk]));
					Energie_xi_help[1]=(pnt_mesh->jacobian[iMinus2jk]*(pnt_mesh->xi_x[iMinus2jk]*pnt_ZD->q_x[iMinus2jk]+pnt_mesh->xi_y[iMinus2jk]*pnt_ZD->q_y[iMinus2jk]+pnt_mesh->xi_z[iMinus2jk]*pnt_ZD->q_z[iMinus2jk]));
					Energie_xi_help[2]=(pnt_mesh->jacobian[iMinus1jk]*(pnt_mesh->xi_x[iMinus1jk]*pnt_ZD->q_x[iMinus1jk]+pnt_mesh->xi_y[iMinus1jk]*pnt_ZD->q_y[iMinus1jk]+pnt_mesh->xi_z[iMinus1jk]*pnt_ZD->q_z[iMinus1jk]));
					Energie_xi_help[4]=(pnt_mesh->jacobian[iPlus1jk]*(pnt_mesh->xi_x[iPlus1jk]*pnt_ZD->q_x[iPlus1jk]+pnt_mesh->xi_y[iPlus1jk]*pnt_ZD->q_y[iPlus1jk]+pnt_mesh->xi_z[iPlus1jk]*pnt_ZD->q_z[iPlus1jk]));
					Energie_xi_help[5]=(pnt_mesh->jacobian[iPlus2jk]*(pnt_mesh->xi_x[iPlus2jk]*pnt_ZD->q_x[iPlus2jk]+pnt_mesh->xi_y[iPlus2jk]*pnt_ZD->q_y[iPlus2jk]+pnt_mesh->xi_z[iPlus2jk]*pnt_ZD->q_z[iPlus2jk]));
					Energie_xi_help[6]=(pnt_mesh->jacobian[iPlus3jk]*(pnt_mesh->xi_x[iPlus3jk]*pnt_ZD->q_x[iPlus3jk]+pnt_mesh->xi_y[iPlus3jk]*pnt_ZD->q_y[iPlus3jk]+pnt_mesh->xi_z[iPlus3jk]*pnt_ZD->q_z[iPlus3jk]));

					Energie_xi=
							pnt_config->ZD_Ableitung_Koeffizient[0]*(Energie_xi_help[0]*pnt_config->Gamma[iMinus3jk]+pnt_config->Psi*(xImpuls_xi_help[0]*pnt_U_RK->u[iMinus3jk]+yImpuls_xi_help[0]*pnt_U_RK->v[iMinus3jk]+zImpuls_xi_help[0]*pnt_U_RK->w[iMinus3jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[1]*(Energie_xi_help[1]*pnt_config->Gamma[iMinus2jk]+pnt_config->Psi*(xImpuls_xi_help[1]*pnt_U_RK->u[iMinus2jk]+yImpuls_xi_help[1]*pnt_U_RK->v[iMinus2jk]+zImpuls_xi_help[1]*pnt_U_RK->w[iMinus2jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[2]*(Energie_xi_help[2]*pnt_config->Gamma[iMinus1jk]+pnt_config->Psi*(xImpuls_xi_help[2]*pnt_U_RK->u[iMinus1jk]+yImpuls_xi_help[2]*pnt_U_RK->v[iMinus1jk]+zImpuls_xi_help[2]*pnt_U_RK->w[iMinus1jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[4]*(Energie_xi_help[4]*pnt_config->Gamma[iPlus1jk]+pnt_config->Psi*(xImpuls_xi_help[4]*pnt_U_RK->u[iPlus1jk]+yImpuls_xi_help[4]*pnt_U_RK->v[iPlus1jk]+zImpuls_xi_help[4]*pnt_U_RK->w[iPlus1jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[5]*(Energie_xi_help[5]*pnt_config->Gamma[iPlus2jk]+pnt_config->Psi*(xImpuls_xi_help[5]*pnt_U_RK->u[iPlus2jk]+yImpuls_xi_help[5]*pnt_U_RK->v[iPlus2jk]+zImpuls_xi_help[5]*pnt_U_RK->w[iPlus2jk]))+
							pnt_config->ZD_Ableitung_Koeffizient[6]*(Energie_xi_help[6]*pnt_config->Gamma[iPlus3jk]+pnt_config->Psi*(xImpuls_xi_help[6]*pnt_U_RK->u[iPlus3jk]+yImpuls_xi_help[6]*pnt_U_RK->v[iPlus3jk]+zImpuls_xi_help[6]*pnt_U_RK->w[iPlus3jk]));
#endif


//					Auf den Divisor deltaXi, deltaEta und deltaZeta wird verzichtet, da diese zu 1 gesetzt wurden
				pnt_Q->xiMomentum[ijk]=pnt_Q->xiMomentum[ijk]+xImpuls_xi;
				pnt_Q->etaMomentum[ijk]=pnt_Q->etaMomentum[ijk]+yImpuls_xi;
#if MESHDIMENSIONS==3
					pnt_Q->zetaMomentum[ijk]=pnt_Q->zetaMomentum[ijk]+zImpuls_xi;
#endif
#if MESHDIMENSIONS==2
					pnt_Q->zetaMomentum[ijk]=0.0;
#endif
				pnt_Q->Energy[ijk]=pnt_Q->Energy[ijk]+Energie_xi;
			}
		}
	}

}

void CalcTauQValues(
		int ijk,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD)
{
//	Das wenoEpsilon verhindert hier nan-Werte sobald y genau auf 0 liegt.
#if MESHDIMENSIONS==2
		pnt_ZD->tau_xx[ijk]=pnt_U_RK->mue[ijk]*(
				4./3.*(pnt_mesh->xi_x[ijk]*pnt_ZD->u_xi[ijk]+pnt_mesh->eta_x[ijk]*pnt_ZD->u_eta[ijk])
				-2./3.*
				(pnt_mesh->xi_y[ijk]*pnt_ZD->v_xi[ijk]+pnt_mesh->eta_y[ijk]*pnt_ZD->v_eta[ijk]
				+pnt_config->flag_rotation_symmetric*pnt_U_RK->v[ijk]/(pnt_mesh->y[ijk]+pnt_config->wenoEpsilon))); //Zusatzterm, der nur fuer rotationssymmetrische Simulationen beruecksichtigt wird
		pnt_ZD->tau_yy[ijk]=pnt_U_RK->mue[ijk]*(
				4./3.*(pnt_mesh->xi_y[ijk]*pnt_ZD->v_xi[ijk]+pnt_mesh->eta_y[ijk]*pnt_ZD->v_eta[ijk])
				-2./3.*
				(pnt_mesh->xi_x[ijk]*pnt_ZD->u_xi[ijk]+pnt_mesh->eta_x[ijk]*pnt_ZD->u_eta[ijk]
				+pnt_config->flag_rotation_symmetric*pnt_U_RK->v[ijk]/(pnt_mesh->y[ijk]+pnt_config->wenoEpsilon))); //Zusatzterm, der nur fuer rotationssymmetrische Simulationen beruecksichtigt wird

		pnt_ZD->tau_xy[ijk]=pnt_U_RK->mue[ijk]*(
				(pnt_mesh->xi_y[ijk]*pnt_ZD->u_xi[ijk]+pnt_mesh->eta_y[ijk]*pnt_ZD->u_eta[ijk])
				+(pnt_mesh->xi_x[ijk]*pnt_ZD->v_xi[ijk]+pnt_mesh->eta_x[ijk]*pnt_ZD->v_eta[ijk]));

		pnt_ZD->q_x[ijk]=pnt_U_RK->mue[ijk]*
				(pnt_mesh->xi_x[ijk]*pnt_ZD->T_xi[ijk]+pnt_mesh->eta_x[ijk]*pnt_ZD->T_eta[ijk]);
		pnt_ZD->q_y[ijk]=pnt_U_RK->mue[ijk]*
				(pnt_mesh->xi_y[ijk]*pnt_ZD->T_xi[ijk]+pnt_mesh->eta_y[ijk]*pnt_ZD->T_eta[ijk]);
#endif

#if MESHDIMENSIONS==3
		pnt_ZD->tau_xx[ijk]=pnt_U_RK->mue[ijk]*(
				4./3.*(pnt_mesh->xi_x[ijk]*pnt_ZD->u_xi[ijk]+pnt_mesh->eta_x[ijk]*pnt_ZD->u_eta[ijk]+pnt_mesh->zeta_x[ijk]*pnt_ZD->u_zeta[ijk])
				-2./3.*(pnt_mesh->xi_y[ijk]*pnt_ZD->v_xi[ijk]+pnt_mesh->eta_y[ijk]*pnt_ZD->v_eta[ijk]+pnt_mesh->zeta_y[ijk]*pnt_ZD->v_zeta[ijk])
				-2./3.*(pnt_mesh->xi_z[ijk]*pnt_ZD->w_xi[ijk]+pnt_mesh->eta_z[ijk]*pnt_ZD->w_eta[ijk]+pnt_mesh->zeta_z[ijk]*pnt_ZD->w_zeta[ijk]));
		pnt_ZD->tau_yy[ijk]=pnt_U_RK->mue[ijk]*(
				4./3.*(pnt_mesh->xi_y[ijk]*pnt_ZD->v_xi[ijk]+pnt_mesh->eta_y[ijk]*pnt_ZD->v_eta[ijk]+pnt_mesh->zeta_y[ijk]*pnt_ZD->v_zeta[ijk])
				-2./3.*(pnt_mesh->xi_x[ijk]*pnt_ZD->u_xi[ijk]+pnt_mesh->eta_x[ijk]*pnt_ZD->u_eta[ijk]+pnt_mesh->zeta_x[ijk]*pnt_ZD->u_zeta[ijk])
				-2./3.*(pnt_mesh->xi_z[ijk]*pnt_ZD->w_xi[ijk]+pnt_mesh->eta_z[ijk]*pnt_ZD->w_eta[ijk]+pnt_mesh->zeta_z[ijk]*pnt_ZD->w_zeta[ijk]));
		pnt_ZD->tau_zz[ijk]=pnt_U_RK->mue[ijk]*(
				4./3.*(pnt_mesh->xi_z[ijk]*pnt_ZD->w_xi[ijk]+pnt_mesh->eta_z[ijk]*pnt_ZD->w_eta[ijk]+pnt_mesh->zeta_z[ijk]*pnt_ZD->w_zeta[ijk])
				-2./3.*(pnt_mesh->xi_x[ijk]*pnt_ZD->u_xi[ijk]+pnt_mesh->eta_x[ijk]*pnt_ZD->u_eta[ijk]+pnt_mesh->zeta_x[ijk]*pnt_ZD->u_zeta[ijk])
				-2./3.*(pnt_mesh->xi_y[ijk]*pnt_ZD->v_xi[ijk]+pnt_mesh->eta_y[ijk]*pnt_ZD->v_eta[ijk]+pnt_mesh->zeta_y[ijk]*pnt_ZD->v_zeta[ijk]));

		pnt_ZD->tau_xy[ijk]=pnt_U_RK->mue[ijk]*(
				(pnt_mesh->xi_y[ijk]*pnt_ZD->u_xi[ijk]+pnt_mesh->eta_y[ijk]*pnt_ZD->u_eta[ijk]+pnt_mesh->zeta_y[ijk]*pnt_ZD->u_zeta[ijk])
				+(pnt_mesh->xi_x[ijk]*pnt_ZD->v_xi[ijk]+pnt_mesh->eta_x[ijk]*pnt_ZD->v_eta[ijk]+pnt_mesh->zeta_x[ijk]*pnt_ZD->v_zeta[ijk]));
		pnt_ZD->tau_xz[ijk]=pnt_U_RK->mue[ijk]*(
				(pnt_mesh->xi_z[ijk]*pnt_ZD->u_xi[ijk]+pnt_mesh->eta_z[ijk]*pnt_ZD->u_eta[ijk]+pnt_mesh->zeta_z[ijk]*pnt_ZD->u_zeta[ijk])
				+(pnt_mesh->xi_x[ijk]*pnt_ZD->w_xi[ijk]+pnt_mesh->eta_x[ijk]*pnt_ZD->w_eta[ijk]+pnt_mesh->zeta_x[ijk]*pnt_ZD->w_zeta[ijk]));
		pnt_ZD->tau_yz[ijk]=pnt_U_RK->mue[ijk]*(
				(pnt_mesh->xi_z[ijk]*pnt_ZD->v_xi[ijk]+pnt_mesh->eta_z[ijk]*pnt_ZD->v_eta[ijk]+pnt_mesh->zeta_z[ijk]*pnt_ZD->v_zeta[ijk])
				+(pnt_mesh->xi_y[ijk]*pnt_ZD->w_xi[ijk]+pnt_mesh->eta_y[ijk]*pnt_ZD->w_eta[ijk]+pnt_mesh->zeta_y[ijk]*pnt_ZD->w_zeta[ijk]));

		pnt_ZD->q_x[ijk]=pnt_U_RK->mue[ijk]*
				(pnt_mesh->xi_x[ijk]*pnt_ZD->T_xi[ijk]+pnt_mesh->eta_x[ijk]*pnt_ZD->T_eta[ijk]+pnt_mesh->zeta_x[ijk]*pnt_ZD->T_zeta[ijk]);
		pnt_ZD->q_y[ijk]=pnt_U_RK->mue[ijk]*
				(pnt_mesh->xi_y[ijk]*pnt_ZD->T_xi[ijk]+pnt_mesh->eta_y[ijk]*pnt_ZD->T_eta[ijk]+pnt_mesh->zeta_y[ijk]*pnt_ZD->T_zeta[ijk]);
		pnt_ZD->q_z[ijk]=pnt_U_RK->mue[ijk]*
				(pnt_mesh->xi_z[ijk]*pnt_ZD->T_xi[ijk]+pnt_mesh->eta_z[ijk]*pnt_ZD->T_eta[ijk]+pnt_mesh->zeta_z[ijk]*pnt_ZD->T_zeta[ijk]);
#endif
}
