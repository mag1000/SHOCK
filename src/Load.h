#ifndef _SHOCK_LOAD_
#define _SHOCK_LOAD_

#include <stdbool.h>

#include "SHOCK.h"

/** Struct for passing variable type arrays
 *
 * Contains an array and type thereof for passing variable type data
 * such as coordinates and solution data.
 */
typedef struct {
	void* ptr; ///< Pointer to the first element of the array
	DataType_t dt; ///< Datatype, usually RealSingle or RealDouble
} vta;

/** Perform loading of data from disk
 *
 * Sets up all data which is obtained from files on disk. The source of
 * the data (i.e. both, the division- and the mesh-file) is described in
 * the configuration structure. The data is then loaded into respective
 * structures as described below.
 *
 * @param[in,out] Pointer to the configuration struct which contains
 * the path to the division file as member chr_DivisionFile and the path
 * to the mesh file as member chr_MeshFile. The following members are
 * filled with data from the files when the loading is performed:
 *
 * - int_meshDimensions from the cells' dimensions
 * - int_zonesize from the processor's current zone's size
 * - MPI_rank from the zone's ordinal in the division file
 * - MPI_rankNeighbours from the processor's neighbours in the grid
 * - TransformMatrixOfInterface from the transform matrices to the
 *   processor's neighbours
 * - BC_Left through BC_Bottom from the type of BCs found on the
 *   interfaces
 * - Translation from the translation for periodic BCs
 * - DonorRangeOfInterfaces from the donor ranges of the interfaces
 *
 *
 * @param[out,caller-allocates] VTA for X data
 *
 * @param[out,caller-allocates] VTA for Y data
 *
 * @param[out,caller-allocates] VTA for Z data
 *
 * @param[out,caller-allocates] VTA for U data
 *
 * @param[out,caller-allocates] VTA for V data
 *
 * @param[out,caller-allocates] VTA for W data
 *
 * @param[out,caller-allocates] VTA for Rho data
 *
 * @param[out,caller-allocates] VTA for P data
 *
 * @return Returns 0 if the operation suceeded, the CG error code
 * otherwise.
 */
bool loadFile( struct strct_configuration*,vta*,vta*,vta*,vta*,vta*,vta*,vta*,vta* );

#endif
