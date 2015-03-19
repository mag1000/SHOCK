#include "SHOCK.h"

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

void CalcDeviationsForDirectViscidFluxes(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK);

#endif
