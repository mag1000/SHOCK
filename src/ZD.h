#include "saiwenos.h"

#ifndef ZD_H
#define ZD_H

void CalcViscidFluxesInXiDirectionDirectly(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q);

void CalcViscidFluxesInEtaDirectionDirectly(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q);

void CalcViscidFluxesInZetaDirectionDirectly(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q);

void CalcViscidFluxesInXiDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD,
		struct strct_Flux * pnt_Q);

void CalcViscidFluxesInEtaDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD,
		struct strct_Flux * pnt_Q);

void CalcViscidFluxesInZetaDirection(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD,
		struct strct_Flux * pnt_Q);

void CalcTauQForViscidFluxes(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD);

void CalcTauQValues(
		int ijk,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD);

void CalcDeviationsForViscidFluxes(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_ZD * pnt_ZD);

void CalcDeviationsForDirectViscidFluxes(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK);

#endif
