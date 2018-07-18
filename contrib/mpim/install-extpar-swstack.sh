#! /bin/bash
#_____________________________________________________________________
#
set -eu
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
export CXX=g++

build=gcc

source $MODULESHOME/init/bash

case $build in
    nag)
        module swap gcc gcc/6.2.0        
        module swap nag nag/6.2
        export FC=nagfor
        ;;
    gcc)
        module swap gcc gcc/6.2.0
        export FC=gfortran
        ;;
    *)
        echo "ERROR: compiler selection not prepared."
        exit 3
        ;;
esac
#_____________________________________________________________________
#
prefix=$HOME/local.$build
src_dir=$(pwd)

export PATH=$prefix/bin:$PATH
#_____________________________________________________________________
#
if [[ ! -e libtool-2.4.6.dep ]]
then
    wget ftp://ftp.gnu.org/gnu/libtool/libtool-2.4.6.tar.gz
    tar xf libtool-2.4.6.tar.gz
    cd libtool-2.4.6
    ./configure --prefix=$prefix
    make install
    cd $src_dir
    touch libtool-2.4.6.dep
fi

if [[ ! -e mpich-3.2.1.dep ]]
then
    if [[ ! -e mpich-3.2.1.tar.gz ]]
    then
        wget http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz
    else
        rm -rf mpich-3.2.1
    fi
    tar xf mpich-3.2.1.tar.gz
    cd mpich-3.2.1
    case $FC in
        nagfor)
            export FFLAGS="-fpp -kind=byte -mismatch"
            export FCFLAGS="-fpp -kind=byte -mismatch"
            ;;
        *)
            export FCFLAGS=""
            ;;
    esac
    ./configure --disable-wrapper-rpath --enable-fortran=all --prefix=$prefix
    unset FCFLAGS
    make -j 8 install
    cd $src_dir
    touch mpich-3.2.1.dep
fi

if [[ ! -e zlib-1.2.11.dep ]]
then
    wget http://zlib.net/zlib-1.2.11.tar.gz
    tar xf zlib-1.2.11.tar.gz
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
    tar xf libaec-1.0.2.tar.gz
    cd libaec-1.0.2
    ./configure --prefix=$prefix
    make install
    cd $src_dir
    touch libaec-1.0.2.dep
fi

if [[ ! -e hdf5-1.10.1.dep ]]
then
    wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.1.tar.gz
    tar xf hdf5-1.10.1.tar.gz
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
    tar xf netcdf-4.5.0.tar.gz
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
    tar xf netcdf-fortran-4.4.4.tar.gz
    cd netcdf-fortran-4.4.4
    libtoolize -c
    autoreconf -i
    case $FC in
        nagfor)
            export FCFLAGS="-fpp -kind=byte -mismatch -L$prefix/lib -Wl,-rpath,$prefix/lib"
            ;;
        *)
            export FCFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib"
            ;;
    esac
    CPPFLAGS="-I$prefix/include" \
    CFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib" \
    ./configure --prefix=$prefix
    unset FCFLAGS
    make -j 1 install
    cd $src_dir
    touch netcdf-fortran-4.4.4.dep
fi

if [[ ! -e netcdf-cxx4-4.3.0.dep ]]
then
    wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-cxx4-4.3.0.tar.gz
    tar xf netcdf-cxx4-4.3.0.tar.gz
    cd netcdf-cxx4-4.3.0
    CPPFLAGS="-I$prefix/include" \
    CFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib -lnetcdf" \
    CXXFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib -lnetcdf" \
    ./configure --prefix=$prefix
    make -j 1 install
    cd $src_dir
    touch netcdf-cxx4-4.3.0.dep
fi

if [[ ! -e eccodes-2.7.0.dep ]]
then
    wget https://software.ecmwf.int/wiki/download/attachments/45757960/eccodes-2.7.0-Source.tar.gz
    tar xf eccodes-2.7.0-Source.tar.gz
    cd eccodes-2.7.0-Source
    mkdir build
    cd build
    case $FC in
        nagfor)
            CMAKE_extra_fortran_flags="-kind=byte -mismatch"
            CMAKE_extra_fortran_options="-Wc,-fPIE"
            ;;
        *)
            CMAKE_extra_fortran_flags=""
            CMAKE_extra_fortran_options=""
            ;;
    esac
    cmake  .. -DCMAKE_INSTALL_PREFIX=$prefix -DENABLE_NETCDF=OFF -DENABLE_AEC=ON -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_Fortran_COMPILE_OPTIONS_PIE="${CMAKE_extra_fortran_options}" -DCMAKE_Fortran_FLAGS="${CMAKE_extra_fortran_flags}" -DENABLE_EXAMPLES=OFF
    make -j 8 install
    cd $src_dir
    touch  eccodes-2.7.0.dep
fi
