#include "SHOCK.h"

#ifndef BC_H
#define BC_H

extern void SetAllBoundaryConditions(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void ControlBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void ControlBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void ControlBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void ControlBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void ControlBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void ControlBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void WriteFarfieldBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteFarfieldBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteFarfieldBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteFarfieldBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteFarfieldBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteFarfieldBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);


extern void WriteInflowSubsonicRiemannBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicRiemannBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicRiemannBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicRiemannBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicRudyBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void WriteOutflowSubsonicRudyBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void WriteOutflowSubsonicRudyBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void WriteOutflowSubsonicRudyBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void WriteOutflowSubsonicRudyBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void WriteOutflowSubsonicRudyBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

extern void WriteInflowSubsonicRiemannBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicRiemannBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallSlipBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallSlipBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallSlipBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallSlipBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallSlipBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallSlipBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallNoSlipBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallNoSlipBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallNoSlipBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallNoSlipBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallNoSlipBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallNoSlipBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSupersonicNormalBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSupersonicNormalBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSupersonicNormalBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSupersonicNormalBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSupersonicNormalBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSupersonicNormalBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicIsentropBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicIsentropBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicIsentropBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicIsentropBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicIsentropBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicIsentropBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicNormalBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicNormalBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicNormalBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicNormalBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicNormalBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteInflowSubsonicNormalBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicRiemannBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicRiemannBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicRiemannBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicRiemannBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicRiemannBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicRiemannBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicNormalBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicNormalBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicNormalBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicNormalBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicNormalBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteOutflowSubsonicNormalBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteFarfieldBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteInflowSubsonicIsentropBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteInflowSubsonicRiemannBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteInflowSubsonicNormalBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteInflowSupersonicNormalBoundary(
		FLT corrector[3],
		int ijk,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteOutflowSubsonicNormalBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteOutflowSubsonicRiemannBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteOutflowSubsonicRudyBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_U * pnt_U_lastStep);

void WriteWallNoSlipBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteMovingWallNoSlipBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteWallSlipBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteWallNoSlipIsothermalBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);
		
void WriteMovingWallNoSlipIsothermalBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);
		
void WriteMovingWallSlipBoundary(
		FLT corrector[3],
		int ijk,
		int ijkSymmetry,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);		

void WriteWallNoSlipIsothermalBoundaryLowerJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteWallNoSlipIsothermalBoundaryLowerI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteWallNoSlipIsothermalBoundaryUpperI(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteWallNoSlipIsothermalBoundaryUpperJ(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void WriteWallNoSlipIsothermalBoundaryUpperK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void WriteWallNoSlipIsothermalBoundaryLowerK(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);
		

#endif
