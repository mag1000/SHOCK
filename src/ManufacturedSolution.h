#include "SHOCK.h"

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L  /* pi */
#endif

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
		struct strct_U * pnt_U_lastStep);

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

