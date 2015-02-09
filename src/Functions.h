#include "SHOCK.h"
#include "stdbool.h"
#include "cgnslib.h"
#include "Load.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void postprocessLoad(
		vta* x,vta* y,vta* z,vta* u,vta* v,vta* w,vta* rho,vta* p,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf,
		struct strct_Flux * pnt_Q,
		struct strct_Flux * pnt_Q_sum,
		struct strct_Film * pnt_Film,
		struct strct_U * pnt_U_backup1,
		struct strct_U * pnt_U_backup2);

void copyMetric(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		int ijkDestination,
		int ijkSource);

void mirrorGhostCellsMetric(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

void check_Metric(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

extern void CreateMetric(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

extern void CreateViscidMetric(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

extern void ExtrapolateGhostCells(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

float IBC_getActualPosition(
		struct strct_configuration * pnt_config);

extern void AllocMemory(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf,
		struct strct_Flux * pnt_Q,
		struct strct_Flux * pnt_Q_sum,
		struct strct_Film * pnt_Film,
		struct strct_U * pnt_U_backup1,
		struct strct_U * pnt_U_backup2);

void AllocMemoryBackup(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U_backup1,
		struct strct_U * pnt_U_backup2);

void IBC_Actual2Last(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

void WriteValuesForPressureHistory(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U_lastStep);

void WriteValuesForVelocityHistory(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U_lastStep);

void InitializeLaminarBoundary(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void preparePressureHistoryValues(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

void prepareVelocityHistoryValues(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

void preparePressureWaves(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

float CalcLambda2(
		int i,
		int j,
		int k,
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void inducePressureWavesPlateau(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void inducePressureWaves(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

extern void AllocMemoryMesh(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

extern void AllocMemoryStrctU(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_strctU);

extern void AllocMemoryStrctFlux(
		struct strct_configuration * pnt_config,
		struct strct_Flux * pnt_strctFlux);

extern void AllocMemoryStrctBoundaryMesh(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_BoundaryMesh,
		int iMeshPoints,
		int jMeshPoints,
		int kMeshPoints);

void FreeMemory(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf,
		struct strct_Flux * pnt_Q,
		struct strct_Flux * pnt_Q_sum,
		struct strct_Film * pnt_Film);

extern void FreeMemoryMesh(
		struct strct_mesh * pnt_mesh,
		struct strct_configuration * pnt_config);

extern void FreeMemoryStrctFilm(
		struct strct_Film * pnt_strctFilm);

extern void FreeMemoryStrctU(
		struct strct_U * pnt_structU);

extern void FreeMemoryStrctFlux(
		struct strct_Flux * pnt_strctFlux);

extern void Initialize(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

extern void InitializeSpecialConditions(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

extern void InitializeMPI(
		int argc,
		char *argv[],
		struct strct_configuration * pnt_config);

extern void InitializeMPI_cgns(
		int argc,
		char *argv[],
		struct strct_configuration * pnt_config);

extern void DefineParameters(
		struct strct_configuration * pnt_config);

void IBC_BornCells(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK);

void IBC_ApplyBC4FluxInXi(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void IBC_ApplyBC4FluxInEta(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void IBC_ApplyBC4FluxInZeta(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);


extern void CalcRungeKutta(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q,
		struct strct_Flux * pnt_Q_sum,
		struct strct_Flux * pnt_Flux,
		struct strct_Flux * pnt_Flux_PlusHalf);

extern void WriteValuesFromU1ToU2(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U1,
		struct strct_U * pnt_U2);

extern void DeleteQ(
		struct strct_configuration * pnt_config,
		struct strct_Flux * pnt_Q);

extern void CalcValuesForPost(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

extern void WriteValuesFromMeshToB(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_mesh * pnt_BoundaryMeshRecieve,
		int iStart,
		int iEnd,
		int jStart,
		int jEnd,
		int kStart,
		int kEnd);

void WriteValuesFromUAndMeshToBuffer(
		float * buffer,
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U,
		struct strct_mesh * pnt_mesh,
		int iStart,
		int iEnd,
		int jStart,
		int jEnd,
		int kStart,
		int kEnd);

int get_ijkTransformFilm(
		struct strct_configuration * pnt_config,
		int interface,
		int i,
		int j,
		int k,
		int c);

extern void TransferFlowParameterWithGhosts(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void TransferForFilm(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U,
		struct strct_Film * pnt_Film);

extern void TransferMeshParameter(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

extern void CalcValues(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

extern void TransferViscidParameter(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

extern void WriteValuesFromBufferToU(
		float * buffer,
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U,
		int iStart,
		int iEnd,
		int jStart,
		int jEnd,
		int kStart,
		int kEnd,
		int interface);

extern void filecopy(
		char *fileIn,
		char *fileOut);

void getNeighbour(
		struct strct_configuration * pnt_config);

int get_ijkTransformWithGhosts(
		struct strct_configuration * pnt_config,
		int interface,
		int i,
		int j,
		int k,
		int c);

int get_ijkTransformMesh(
		struct strct_configuration * pnt_config,
		int interface,
		int i,
		int j,
		int k,
		int c);

void getInterfaceInformations(
		struct strct_configuration * pnt_config);

extern void WriteValuesFromUToFilm(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U,
		struct strct_Film * pnt_Film,
		struct strct_mesh * pnt_mesh);

extern void AllocMemoryStrctFilm(
		struct strct_configuration * pnt_config,
		struct strct_Film * pnt_strctFilm);

extern void check_TransformationMatrix(
		int interface,
		struct strct_configuration * pnt_config);

extern void check_Connectivity(
		struct strct_configuration * pnt_config,
		struct strct_mesh* pnt_mesh);

extern int check_CGNSFile(
		struct strct_configuration * pnt_config);

extern int check_ConfigFile(
		struct strct_configuration * pnt_config);

void IBC_set(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void IBC_prepare(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

void InitializeVortex(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void DefineFilesPath(
		struct strct_configuration * pnt_config);

void SplitMeshFile(
		struct strct_configuration * pnt_config);

int checkNAN(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U_lastStep,
		struct strct_mesh * pnt_mesh);

void changeMeshInto2D(
		struct strct_configuration * pnt_config);

void changeMeshInto3D(
		struct strct_configuration * pnt_config);

void WriteConstantZValues(
		struct strct_configuration * pnt_config,
		struct strct_U * pnt_U_lastStep);

void AddRotationSymmetricFluxes(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_RK,
		struct strct_Flux * pnt_Q);

void freeVTA
(
		vta* x,vta* y,vta* z,vta* u,vta* v,vta* w,vta* rho,vta* p);

#endif
