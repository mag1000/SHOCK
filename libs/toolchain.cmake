# The compilers
set(CMAKE_C_COMPILER /bgsys/drivers/ppcfloor/comm/xl/bin/mpicc)
set(CMAKE_CXX_COMPILER /bgsys/drivers/ppcfloor/comm/xl/bin/mpicxx)
set(CMAKE_Fortran_COMPILER /bgsys/drivers/ppcfloor/comm/xl/bin/mpif90)

# The MPI wrappers for the serial compilers
set(MPI_ROOT "/bgsys/drivers/ppcfloor/comm/gcc")

set(MPI_C_COMPILER ${MPI_ROOT}/bin/mpixlc_r)
set(MPI_CXX_COMPILER ${MPI_ROOT}/bin/mpixlcxx_r)
set(MPI_Fortran_COMPILER ${MPI_ROOT}/bin/mpixlf90_r)

set(MPI_C_COMPILE_FLAGS       "")
set(MPI_CXX_COMPILE_FLAGS     "")
set(MPI_Fortran_COMPILE_FLAGS "")

set(MPI_C_INCLUDE_PATH "${MPI_ROOT}/include")
set(MPI_CXX_INCLUDE_PATH "${MPI_ROOT}/include")
set(MPI_Fortran_INCLUDE_PATH "${MPI_ROOT}/include")

set(MPI_C_LINK_FLAGS "-L${MPI_ROOT}/lib -L${PAMI_ROOT}/lib -L${SPI_ROOT}/lib")
set(MPI_CXX_LINK_FLAGS ${MPI_C_LINK_FLAGS})


set(MPI_C_LIBRARIES "/bgsys/drivers/V1R2M2/ppc64/comm/lib/libmpich-xl.a;/bgsys/drivers/V1R2M2/ppc64/comm/lib/libopa-xl.a;/bgsys/drivers/V1R2M2/ppc64/comm/lib/libmpl-xl.a;/bgsys/drivers/V1R2M2/ppc64/comm/lib/libpami-gcc.a;/bgsys/drivers/V1R2M2/ppc64/spi/lib/libSPI.a;/bgsys/drivers/V1R2M2/ppc64/spi/lib/libSPI_cnk.a;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/librt.a;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/libpthread.so;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/libstdc++.a;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/libpthread.so")
set(MPI_CXX_LIBRARIES "/bgsys/drivers/V1R2M2/ppc64/comm/lib/libmpichcxx-xl.a;/bgsys/drivers/V1R2M2/ppc64/comm/lib/libmpich-xl.a;/bgsys/drivers/V1R2M2/ppc64/comm/lib/libopa-xl.a;/bgsys/drivers/V1R2M2/ppc64/comm/lib/libmpl-xl.a;/bgsys/drivers/V1R2M2/ppc64/comm/lib/libpami-gcc.a;/bgsys/drivers/V1R2M2/ppc64/spi/lib/libSPI.a;/bgsys/drivers/V1R2M2/ppc64/spi/lib/libSPI_cnk.a;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/librt.a;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/libpthread.so;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/libstdc++.a;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/libpthread.so")
set(MPI_Fortran_LIBRARIES "/bgsys/drivers/V1R2M2/ppc64/comm/lib/libmpichf90-xl.a;/bgsys/drivers/V1R2M2/ppc64/comm/lib/libmpich-xl.a;/bgsys/drivers/V1R2M2/ppc64/comm/lib/libopa-xl.a;/bgsys/drivers/V1R2M2/ppc64/comm/lib/libmpl-xl.a;/bgsys/drivers/V1R2M2/ppc64/comm/lib/libpami-gcc.a;/bgsys/drivers/V1R2M2/ppc64/spi/lib/libSPI.a;/bgsys/drivers/V1R2M2/ppc64/spi/lib/libSPI_cnk.a;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/librt.a;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/libpthread.so;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/libstdc++.a;/bgsys/drivers/toolchain/V1R2M2_base/gnu-linux/powerpc64-bgq-linux/lib/libpthread.so")

set(CMAKE_LINKER "/bgsys/drivers/ppcfloor/gnu-linux/powerpc64-bgq-linux/bin/ld")

set(CMAKE_INSTALL_PREFIX=$HOME/CGNS)

set(HDF5_INCLUDE_PATH "/bgsys/local/hdf5/include")
set(HDF5_LIBRARY "/bgsys/local/hdf5/lib/libhdf5.a")
 
set(SZIP_LIBRARY "/bgsys/local/szip/lib/libsz.a")
set(ZLIB_LIBRARY "/bgsys/local/zlib/lib/libz.a")


##############################################################

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled.aftware
set(CMAKE_FIND_ROOT_PATH
    /bgsys/drivers/ppcfloor/
    /bgsys/drivers/ppcfloor/comm/xl
    /bgsys/drivers/ppcfloor/comm/sys/
    /bgsys/drivers/ppcfloor/spi/)

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

##############################################################


set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
set(CMAKE_EXE_LINKER_FLAGS "-static")
set(MATH_LIBS "${LAPACK_FLAGS} ${ESSL_FLAGS} ${XLF_FLAGS}")

