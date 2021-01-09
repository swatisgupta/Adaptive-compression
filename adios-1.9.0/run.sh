#!/bin/bash
#module=netcdf/intel/2013.1.039/intelmpi/4.3.2
#mpiicc=/cell_root/software/intel/ics_2013.1.039/impi/4.1.1.036/intel64/bin/mpiicc
MPIDIR=/cell_root/software/intel/ics_2013.1.039/impi/4.1.1.036/intel64

export MPICC=mpiicc 
export MPICXX=mpiicpc 
export MPIFC=mpiifort 
export CC=mpiicc  
export CXX=mpiicpc  
export FC=ifort 
export FCFLAGS=" -pg -I${MPIDIR}/include -I/homes/ssinghal/Builds/include/lzo -I/homes/ssinghal/Builds/include"
export CFLAGS=" -pg -fPIC -I${MPIDIR}/include -I/homes/ssinghal/Builds/include -I/homes/ssinghal/Builds/include/lzo"
export LDLAGS="-pg -L${MPIDIR}/lib -lmpi -L/homes/ssinghal/Builds/lib -llzo2 -lbz2" 
export LIBS=" -pg -L${MPIDIR}/lib -lmpi -L/homes/ssinghal/Builds/lib -llzo2 -lbz2"
export CPPFLAGS="-pg -fPIC -I${MPIDIR}/include -I/homes/ssinghal/Builds/include -I/homes/ssinghal/Builds/include/lzo"
export CXXFLAG="-pg -fPIC -I${MPIDIR}/include -I/homes/ssinghal/Builds/include -I/homes/ssinghal/Builds/include/lzo" 
export echo=echo


BDIR=/homes/ssinghal/Builds
NCDIR=/homes/ssinghal/Builds
#NCDIR=/cell_root/software/netcdf/4.3.2/intel/2013.1.039/intelmpi/hdf5/1.8.13/hdf4/4.2.10/sys/

#make distclean
make clean

#./configure --prefix=$BDIR --disable-fortran --enable-timer-events --with-comp="/lustre/ssinghal/WRF-LETKF/model/model/zlib" --with-mxml=$BDIR --with-lustre=/usr/lib464 --with-nc4par=$NCDIR --with-zlib="/lustre/ssinghal/WRF-LETKF/model/model/zlib"
#./configure --prefix=$BDIR --enable-timer-events --with-comp="/lustre/ssinghal/WRF-LETKF/model/model/zlib" --with-mxml=$BDIR --with-lustre=/usr/lib464 --with-nc4par=$NCDIR --with-zlib="/lustre/ssinghal/WRF-LETKF/model/model/zlib" --with-bzip2="/homes/ssinghal/Builds"
#./configure --prefix=/homes/ssinghal/Builds --enable-timer-events --with-acomp=/usr --with-mxml=/homes/ssinghal/Builds --with-lustre=/usr/ --with-nc4par=/homes/ssinghal/Builds --with-zlib=/usr --with-bzip2=/homes/ssinghal/Builds --no-create --no-recursion
./configure --prefix=/homes/ssinghal/Builds --enable-timer-events --with-acomp=/usr --with-mxml=/homes/ssinghal/Builds      --with-lustre=/usr/ --with-nc4par=/homes/ssinghal/Builds --with-zlib=/usr --with-bzip2=/homes/ssinghal/Builds --no-create      --no-recursion
make 
make check 
make install
