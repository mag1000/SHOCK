#include "saiwenos.h"

#ifndef IMPORT_H
#define IMPORT_H

extern void ConfigImport(
		struct strct_configuration * pnt_config);

extern void MeshImport_CGNS(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh);

extern void MeshConfig_CGNS(
		struct strct_configuration * pnt_config);

extern void ResultImport_CGNS(
		struct strct_configuration * pnt_config,
		struct strct_mesh * pnt_mesh,
		struct strct_U * pnt_U_lastStep);

extern void BCImport_CGNS(
		struct strct_configuration * pnt_config);

#endif
