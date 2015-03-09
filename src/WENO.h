#include "SHOCK.h"

#ifndef WENO_H
#define WENO_H

extern FLT Phi_Function_W9(
		struct strct_configuration * pnt_config,
		FLT * pnt_flux,
		FLT * pnt_deltaFlux);

extern FLT Phi_Function_W5(
		struct strct_configuration * pnt_config,
		FLT * pnt_flux,
		FLT * pnt_deltaFlux);

extern void CalcFluxesInXiDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf,
		struct strct_Flux * pnt_Q);

extern void CalcFluxesInEtaDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf,
		struct strct_Flux * pnt_Q);

extern void CalcFluxesInZetaDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf,
		struct strct_Flux * pnt_Q);

extern void CalcEigenVectorsInXiDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		int ijk,
		int iPlus1jk,
		int iMinus1jk);

extern void CalcEigenVectorsInEtaDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		int ijk,
		int ijPlus1k,
		int ijMinus1k);

extern void CalcEigenVectorsInZetaDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		int ijk,
		int ijkPlus1,
		int ijkMinus1);

extern FLT GetLambdaMaxInXiDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		int int_acutalI,
		int int_acutalJ,
		int int_acutalK);

extern FLT GetLambdaMaxInEtaDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		int int_acutalI,
		int int_acutalJ,
		int int_acutalK);

extern FLT GetLambdaMaxInZetaDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		int int_acutalI,
		int int_acutalJ,
		int int_acutalK);

extern void Theta_Function_W9(
	struct strct_Flux * pnt_Flux,
	FLT * pnt_theta,
	int ijk,
	int Plus1ijk,
	int Plus2ijk,
	int Plus3ijk,
	int Plus4ijk,
	int Minus1ijk,
	int Minus2ijk,
	int Minus3ijk);

void Theta_Function_W5(
    struct strct_Flux * pnt_Flux,
    FLT * pnt_theta,
    int ijk,
    int Plus1ijk,
    int Plus2ijk,
    int Minus1ijk);

#endif
