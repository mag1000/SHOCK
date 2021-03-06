#include "SHOCK.h"

#ifndef MANUFACTUREDSOLUTION_H
#define MANUFACTUREDSOLUTION_H

void AddManufacturedSolutionSource(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_Flux * pnt_Q);

void ConfigureManufacturedSolution(
		struct strct_configuration * pnt_config);

void InitializeManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void WriteBCManufacturedSolutionLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void WriteBCManufacturedSolutionLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void WriteBCManufacturedSolutionLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void WriteBCManufacturedSolutionUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void WriteBCManufacturedSolutionUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void WriteBCManufacturedSolutionUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void WriteManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		int ijk);

void ErrorManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_Film * pnt_Film,
		struct strct_Flux * pnt_Q_sum,
		int flag);

long double GetRhoManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		int ijk);

long double GetPressureManufacturedSolution(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		int ijk);

#endif

