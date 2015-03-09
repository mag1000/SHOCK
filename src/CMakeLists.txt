

cmake_minimum_required(VERSION 2.6)

project(SHOCKProject)

mark_as_advanced(FORCE MPI_EXTRA_LIBRARY)
mark_as_advanced(CLEAR CMAKE_C_COMPILER)

#Custom build types
SET(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING "Choose the type of build, options are: Local, Debug, Cluster, Juqueen." FORCE)

#Build type local
SET(CMAKE_C_FLAGS_LOCAL "-O3 -Wall" CACHE STRING "Flags used by C compiler for local builds." FORCE)
IF (CMAKE_BUILD_TYPE MATCHES Local)
    add_definitions("-DINTEL")
    #find_path(HDF5_INCLUDE_PATH hdf5.h CACHE FILEPATH "HDF5 include path")
    find_library(HDF5_LIBRARY hdf5 CACHE FILEPATH "HDF5 library path")
    find_path(CGNS_INCLUDE_PATH cgnslib.h CACHE FILEPATH "CGNS include path")
    include_directories(${CGNS_INCLUDE_PATH})
    find_library(CGNS_LIBRARY cgns CACHE FILEPATH "CGNS library path")
    #find_path(SZIP_INCLUDE_PATH szlib.h CACHE FILEPATH "SZIP include path")
    find_library(SZIP_LIBRARY szip CACHE FILEPATH "SZIP library path")
    #find_path(ZLIB_INCLUDE_PATH zlib.h CACHE FILEPATH "ZLIB include path")
    find_library(ZLIB_LIBRARY z CACHE FILEPATH "ZLIB library path")
    SET (SHOCK_LIBS ${CGNS_LIBRARY} ${HDF5_LIBRARY} ${SZIP_LIBRARY} ${ZLIB_LIBRARY} -lm)
ENDIF (CMAKE_BUILD_TYPE MATCHES Local)

#Build type cluster
SET(CMAKE_C_FLAGS_CLUSTER "$ENV{FLAGS_FAST}" CACHE STRING "Flags used by C compiler for builds on RWTH Cluster." FORCE)
IF (CMAKE_BUILD_TYPE MATCHES Cluster)
	add_definitions("-DINTEL")
    SET (HDF5_INCLUDE_PATH "/home/swlcfd/HDF5/include" CACHE FILEPATH "HDF5 include path")
    SET (HDF5_LIBRARY "/home/swlcfd/HDF5/lib/libhdf5.so" CACHE FILEPATH "HDF5 library path")
    SET (CGNS_INCLUDE_PATH "/home/swlcfd/CGNS/include" CACHE FILEPATH "CGNS include path")
    SET (CGNS_LIBRARY "/home/swlcfd/CGNS/lib/libcgns.so" CACHE FILEPATH "CGNS library path")
    SET (SHOCK_LIBS ${CGNS_LIBRARY} ${HDF5_LIBRARY} ${SZIP_LIBRARY} ${ZLIB_LIBRARY})
ENDIF (CMAKE_BUILD_TYPE MATCHES Cluster)

#Build type cluster_debug
SET(CMAKE_C_FLAGS_CLUSTER_DEBUG "$ENV{FLAGS_DEBUG} -Wall" CACHE STRING "Flags used by C compiler for debug-builds on RWTH Cluster." FORCE)
IF (CMAKE_BUILD_TYPE MATCHES Cluster_Debug)
    add_definitions("-DINTEL")
    SET (HDF5_INCLUDE_PATH "/home/swlcfd/HDF5/include" CACHE FILEPATH "HDF5 include path")
    SET (HDF5_LIBRARY "/home/swlcfd/HDF5/lib/libhdf5.so" CACHE FILEPATH "HDF5 library path")
    SET (CGNS_INCLUDE_PATH "/home/swlcfd/CGNS/include" CACHE FILEPATH "CGNS include path")
    SET (CGNS_LIBRARY "/home/swlcfd/CGNS/lib/libcgns.so" CACHE FILEPATH "CGNS library path")
    SET (SHOCK_LIBS ${CGNS_LIBRARY} ${HDF5_LIBRARY} ${SZIP_LIBRARY} ${ZLIB_LIBRARY})
ENDIF (CMAKE_BUILD_TYPE MATCHES Cluster_Debug)

#Build type juqueen
SET(CMAKE_C_FLAGS_JUQUEEN "-O5 -qarch=qp -qtune=qp" CACHE STRING "Flags used by C compiler for builds on Juqueen." FORCE)
IF (CMAKE_BUILD_TYPE MATCHES Juqueen)
    add_definitions("-DIBM")
    #SET(HDF5_INCLUDE_PATH "/homec/jhpc13/jhpc1301/packages/phdf5-trunk/include" CACHE FILEPATH "HDF5 include path")
    #SET(HDF5_LIBRARY "/homec/jhpc13/jhpc1301/packages/phdf5-trunk/lib/libhdf5.a" CACHE FILEPATH "HDF5 library path")
    #SET(CGNS_INCLUDE_PATH "/homec/jhpc13/jhpc1301/packages/cgns-trunk/include" CACHE FILEPATH "CGNS include path")
    #SET(CGNS_LIBRARY "/homec/jhpc13/jhpc1301/packages/cgns-trunk/lib/libcgns.a" CACHE FILEPATH "CGNS library path")
    #SET (ZLIB_INCLUDE_PATH "/homec/jhpc13/jhpc1301/packages/zlib-1.2.8/lib/include" CACHE FILEPATH "Zlib include path")
    #SET (ZLIB_LIBRARY "/homec/jhpc13/jhpc1301/packages/zlib-1.2.8/lib/lib/libz.a" CACHE FILEPATH "Zlib library path")
	SET(HDF5_INCLUDE_PATH "/bgsys/local/hdf5/include" CACHE FILEPATH "HDF5 include path")
    SET(HDF5_LIBRARY "/bgsys/local/hdf5/lib/libhdf5.a" CACHE FILEPATH "HDF5 library path")
    SET(CGNS_INCLUDE_PATH "/homec/jhpc13/jhpc1301/CGNS/include" CACHE FILEPATH "CGNS include path")
    SET(CGNS_LIBRARY "/homec/jhpc13/jhpc1301/CGNS/lib/libcgns.a" CACHE FILEPATH "CGNS library path")
    SET (ZLIB_INCLUDE_PATH "/bgsys/local/zlib/include" CACHE FILEPATH "Zlib include path")
    SET (ZLIB_LIBRARY "/bgsys/local/zlib/lib/libz.a" CACHE FILEPATH "Zlib library path")    
    include_directories(${ZLIB_INCLUDE_PATH})
    SET (SZIP_INCLUDE_PATH "/bgsys/local/szip/include" CACHE FILEPATH "Szip include path")
    SET (SZIP_LIBRARY "/bgsys/local/szip/lib/libsz.a" CACHE FILEPATH "Szip library path")
    include_directories(${SZIP_INCLUDE_PATH})
    SET (SHOCK_LIBS ${CGNS_LIBRARY} ${HDF5_LIBRARY} ${SZIP_LIBRARY} ${ZLIB_LIBRARY} )
ENDIF (CMAKE_BUILD_TYPE MATCHES Juqueen)

set(PRECISION 1 CACHE STRING "Precision: 1(float) or 2(double) or 3 (long double) or 4 (quad)")
if(NOT (PRECISION STREQUAL 1 OR PRECISION STREQUAL 2 OR PRECISION STREQUAL 3 OR PRECISION STREQUAL 4))
	message(WARNING "PRECISION must be either 1(float) or 2(double) or 3 (long double) or 4 (quad)")
	set (PRECISION 1 CACHE STRING "Precision: 1(float) or 2(double) or 3 (long double) or 4 (quad)" FORCE)
endif()

set(SPACEORDER 5 CACHE STRING "Space Order: 9 or 5")
if(NOT (SPACEORDER STREQUAL 5 OR SPACEORDER STREQUAL 9))
	message(WARNING "SPACEORDER must be either 9 or 5")
	set (SPACEORDER 9 CACHE STRING "Space Order: 9 or 5" FORCE)
endif()

set(MESHDIMENSIONS 3 CACHE STRING "Mesh Dimensions: 2 or 3")
if(NOT (MESHDIMENSIONS STREQUAL 2 OR MESHDIMENSIONS STREQUAL 3))
	message(WARNING "MESHDIMENSIONS must be either 2 or 3")
	set (MESHDIMENSIONS 2 CACHE STRING "Mesh Dimensions: 2 or 3" FORCE)
endif()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../.git)
	FIND_PACKAGE(Git)
	IF(GIT_FOUND)
		EXECUTE_PROCESS(
			COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
			WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/.."
			OUTPUT_VARIABLE "GITID"
			ERROR_QUIET
			OUTPUT_STRIP_TRAILING_WHITESPACE)
	ELSE(GIT_FOUND)
		SET(GITID 0)
	ENDIF(GIT_FOUND)
ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../.git)

add_definitions("-DGITID=\"${GITID}\"")
add_definitions("-DPRECISION=${PRECISION}")
add_definitions("-DSPACEORDER=${SPACEORDER}")
add_definitions("-DMESHDIMENSIONS=${MESHDIMENSIONS}")

mark_as_advanced(CMAKE_C_FLAGS_LOCAL CMAKE_C_FLAGS_CLUSTER CMAKE_C_FLAGS_CLUSTER_DEBUG CMAKE_C_FLAGS_JUQUEEN)
include_directories(${HDF5_INCLUDE_PATH})
include_directories(${CGNS_INCLUDE_PATH})
include_directories(${INIPARSER_INCLUDE_PATH})
include(FindMPI)

SET (SHOCK_SRCS dictionary.c iniparser.c BC.c Export.c Functions.c Import.c SHOCK.c WENO.c ZD.c Load.c Save.c ManufacturedSolution.c)
add_executable (SHOCK_W${SPACEORDER}_${MESHDIMENSIONS}D ${SHOCK_SRCS})
target_link_libraries (SHOCK_W${SPACEORDER}_${MESHDIMENSIONS}D ${SHOCK_LIBS})
