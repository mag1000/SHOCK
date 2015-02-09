#ifndef _SHOCK_SAVE
#define _SHOCK_SAVE

#include <stdbool.h>
#include <pcgnslib.h>

#include "SHOCK.h"

/** Store data to disk
 *
 * Stores the current state of the solution in a meshfile. If the target
 * file is not the same file which was used for loading, the input file
 * will be duplicated first. After that point, only solution data from
 * struct strct_U is written, because all other data are assumed to not
 * have changed.
 *
 * @param The configuration structure, pointing to the target as
 * chr_FilmFile. If chr_FilmFile is not chr_MeshFile (the input), it
 * will be replaced by a copy of chr_MeshFile. The number of time steps
 * to be found in the solution structure is to be given in int_Samples
 *
 * @param The solution structure. The member flt_time_dim designates the
 * according times. All other data will be taken from all its members
 * and written into the solution file.
 *
 * @return Returns true on success. If the saving failed, for example
 * because the target file did not exist and the input file was not
 * available, will return false.
 */
void saveFile( struct strct_configuration*,struct strct_Film* );

#endif
