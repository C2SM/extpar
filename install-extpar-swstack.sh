#! /bin/bash
#_____________________________________________________________________
#
unset CC
unset FC
unset CXX

unset CFLAGS
unset FCFLAGS
unset FFLAGS
unset CXXFLAGS

export CC=gcc
export FC=nagfor
export CXX=g++
#_____________________________________________________________________
#
src_dir=$(pwd)
prefix=$(pwd)

export PATH=$prefix/bin:$PATH
#_____________________________________________________________________
#
if [[ ! -e libtool-2.4.6.dep ]]
then
    wget ftp://ftp.gnu.org/gnu/libtool/libtool-2.4.6.tar.gz
    tar xvf libtool-2.4.6.tar.gz
    cd libtool-2.4.6
    ./configure --prefix=$prefix
    make install    
    cd $src_dir
    touch libtool-2.4.6.dep 
fi

if [[ ! -e zlib-1.2.11.dep ]]
then
    wget http://zlib.net/zlib-1.2.11.tar.gz
    tar xvf zlib-1.2.11.tar.gz
    cd zlib-1.2.11
    export INSTALLDIR=$prefix
    ./configure --prefix=$prefix
    make install    
    unset INSTALLDIR
    cd $src_dir
    touch zlib-1.2.11.dep 
fi

if [[ ! -e libaec-1.0.2.dep ]]
then
    wget https://gitlab.dkrz.de/k202009/libaec/uploads/b30adc1abf805d7454896ab83c900eb8/libaec-1.0.2.tar.gz
    tar xvf libaec-1.0.2.tar.gz
    cd libaec-1.0.2
    ./configure --prefix=$prefix
    make install    
    cd $src_dir
    touch libaec-1.0.2.dep 
fi

if [[ ! -e hdf5-1.10.1.dep ]]
then
    wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.1.tar.gz
    tar xvf hdf5-1.10.1.tar.gz 
    cd hdf5-1.10.1
    CPPFLAGS="-I$prefix/include" LDFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib" \
    ./configure --disable-fortran --disable-cxx --enable-threadsafe  --enable-unsupported \
                --prefix=$prefix --with-zlib=$prefix --with-szlib=$prefix
    make -j 4
    make install
    cd $src_dir
    touch hdf5-1.10.1.dep 
fi

if [[ ! -e netcdf-4.5.0.dep ]]
then
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.5.0.tar.gz
    tar xvf netcdf-4.5.0.tar.gz
    cd netcdf-4.5.0
    CPPFLAGS="-I$prefix/include" \
    LDFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib" \
    ./configure --enable-netcdf4 --disable-dap --prefix=$prefix
    make -j 8 install    
    cd $src_dir
    touch netcdf-4.5.0.dep 
fi

if [[ ! -e netcdf-fortran-4.4.4.dep ]]
then
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz
    tar xvf netcdf-fortran-4.4.4.tar.gz
    cd netcdf-fortran-4.4.4
    libtoolize -c
    autoreconf -i
    FCFLAGS="-fpp -kind=byte -mismatch -L$prefix/lib -Wl,-rpath,$prefix/lib" \
    CPPFLAGS="-I$prefix/include" \
    CFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib" \
    ./configure --prefix=$prefix
    make -j 1 install    
    cd $src_dir
    touch netcdf-fortran-4.4.4.dep 
fi

if [[ ! -e eccodes-2.7.0.dep ]]
then
    wget https://software.ecmwf.int/wiki/download/attachments/45757960/eccodes-2.7.0-Source.tar.gz
    tar xvf eccodes-2.7.0-Source.tar.gz
    cd eccodes-2.7.0-Source
    mkdir build
    cd build
    cmake  .. -DCMAKE_INSTALL_PREFIX=$prefix -DENABLE_NETCDF=OFF -DENABLE_AEC=ON -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_Fortran_COMPILE_OPTIONS_PIE="-Wc,-fPIE" -DCMAKE_Fortran_FLAGS="-kind=byte -mismatch" -DENABLE_EXAMPLES=OFF
    make -j 8 install    
    cd $src_dir
    touch  eccodes-2.7.0.dep 
fi
