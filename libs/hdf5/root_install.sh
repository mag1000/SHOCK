#
#SZLIB
#
export CC=mpicc
tar -zxvf szip-2.1.tar.gz 
cd szip-2.1
./configure --prefix=/usr/local
make
#make check
make install
cd ..
#
#ZLIB
#
export CC=mpicc
tar -zxvf zlib-1.2.8.tar.gz
cd zlib-1.2.8
./configure --prefix=/usr/local
make
#make check
make install
cd ..
#
#HDF5
#
export CC=mpicc
tar -zxvf hdf5-1.8.14.tar.gz
cd hdf5-1.8.14
./configure --prefix=/usr/local --with-szlib=/usr/local --with-zlib=/usr/local/include,/usr/local/lib --enable-production --enable-shared
make
#make check
make install
cd ..

rm -fr szip-2.1
rm -fr zlib-1.2.8
rm -fr hdf5-1.8.14
