#include "saiwenos.h"

#ifndef EXPORT_H
#define EXPORT_H

void SnapshotExport(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void SnapshotExportEachCPU(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void MeshMetricExport(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

void CGNS_SnapshotExportParallel(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

void CGNS_FilmExportParallel(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_Film * pnt_Film);

void CGNS_PressureHistoryValuesExportParallel(
		struct strct_configuration * pnt_config,
		int p_index);

void ASCii_PressureHistoryValuesExportParallel(
		struct strct_configuration * pnt_config,
		int p_index);

void ASCii_VelocityHistoryValuesExportParallel(
		struct strct_configuration * pnt_config,
		int p_index);

void NANExport(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U);

#endif
