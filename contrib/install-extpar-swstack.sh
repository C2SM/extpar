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

build=gcc

export CC=gcc
export CXX=g++

#_____________________________________________________________________
# mistral specifics

#source $MODULESHOME/init/bash

#module unload cmake
#module load cmake
#module load CMake

case $build in
    nag)
        module unload gcc
        module unload nag
        module load gcc/6.4.0        
        module load nag/6.2
        export FC=nagfor
        ;;
    gcc)
        # module unload gcc
        # module load gcc/6.4.0
        module load gcc/7.2.0        
        # module load GCC
        export FC=gfortran
        ;;
    intel)
	module unload gcc
        module load gcc/6.4.0
        module unload intel
        module load intel/18.0.4
        export FC=ifort
        ;;
    *)
        echo "ERROR: compiler selection not prepared."
        exit 3
        ;;
esac

module list

#_____________________________________________________________________
#
prefix=$HOME/icon-preprocessing/local.$build
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
    make clean
    cd $src_dir
    touch libtool-2.4.6.dep
fi

if [[ ! -e bison-3.5.dep ]]
then
    wget https://ftp.gnu.org/gnu/bison/bison-3.5.tar.gz
    tar xf bison-3.5.tar.gz
    cd bison-3.5
    ./configure --prefix=$prefix
    make -j 4 install
    make clean
    cd $src_dir
    touch bison-3.5.dep
fi

if [[ ! -e flex-2.6.3.dep ]]
then
    wget https://github.com/westes/flex/releases/download/v2.6.3/flex-2.6.3.tar.gz
    tar xf flex-2.6.3.tar.gz
    cd flex-2.6.3
    ./configure --prefix=$prefix
    make -j 4 install
    make clean
    cd $src_dir
    touch flex-2.6.3.dep
fi

if [[ ! -e cmake-3.16.4.dep ]]
then
    wget https://github.com/Kitware/CMake/releases/download/v3.16.4/cmake-3.16.4.tar.gz
    tar xf cmake-3.16.4.tar.gz
    cd cmake-3.16.4
    ./bootstrap  --prefix=$prefix
    make -j 4
    make install
    make clean
    cd $src_dir
    touch cmake-3.16.4.dep
fi

if [[ ! -e mpich-3.3.2.dep ]]
then
    if [[ ! -e mpich-3.3.2.tar.gz ]]
    then
        wget http://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz
    else
        rm -rf mpich-3.3.2
    fi
    tar xf mpich-3.3.2.tar.gz
    cd mpich-3.3.2
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
    make -j 16 install
    make clean
    cd $src_dir
    touch mpich-3.3.2.dep
fi

if [[ ! -e zlib-1.2.11.dep ]]
then
    wget http://zlib.net/zlib-1.2.11.tar.gz
    tar xf zlib-1.2.11.tar.gz
    cd zlib-1.2.11
    export INSTALLDIR=$prefix
    ./configure --prefix=$prefix
    make install
    make clean
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
    make clean
    cd $src_dir
    touch libaec-1.0.2.dep
fi

if [[ ! -e hdf5-1.10.6.dep ]]
then
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.6/src/hdf5-1.10.6.tar.gz
    tar xf hdf5-1.10.6.tar.gz
    cd hdf5-1.10.6
    CPPFLAGS="-I$prefix/include" LDFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib" \
    ./configure --disable-fortran --disable-cxx --enable-threadsafe  --enable-unsupported \
                --prefix=$prefix --with-zlib=$prefix --with-szlib=$prefix
    make -j 8
    make install
    make clean
    cd $src_dir
    touch hdf5-1.10.6.dep
fi

if [[ ! -e netcdf-c-4.7.3.dep ]]
then
    wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-c-4.7.3.tar.gz
    tar xvf netcdf-c-4.7.3.tar.gz
    cd netcdf-c-4.7.3
    CPPFLAGS="-I$prefix/include" \
    LDFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib" \
    ./configure --enable-netcdf4 --disable-dap --prefix=$prefix
    make -j 8 install    
    make clean
    cd $src_dir
    touch netcdf-c-4.7.3.dep 
fi

if [[ ! -e netcdf-fortran-4.5.2.dep ]]
then
    wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-fortran-4.5.2.tar.gz
    tar xf netcdf-fortran-4.5.2.tar.gz
    cd netcdf-fortran-4.5.2
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
    make clean
    cd $src_dir
    touch netcdf-fortran-4.5.2.dep
fi

if [[ ! -e netcdf-cxx4-4.3.1.dep ]]
then
    wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-cxx4-4.3.1.tar.gz
    tar xf netcdf-cxx4-4.3.1.tar.gz
    cd netcdf-cxx4-4.3.1
    CPPFLAGS="-I$prefix/include" \
    CFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib -lnetcdf" \
    CXXFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib -lnetcdf" \
    ./configure --prefix=$prefix
    make -j 1 install
    make clean
    cd $src_dir
    touch netcdf-cxx4-4.3.1.dep
fi

if [[ ! -e eccodes-2.16.0.dep ]]
then
    wget https://confluence.ecmwf.int/download/attachments/45757960/eccodes-2.16.0-Source.tar.gz
    tar xf eccodes-2.16.0-Source.tar.gz
    cd eccodes-2.16.0-Source
    mkdir -p build
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
    cmake  .. -DCMAKE_INSTALL_PREFIX=$prefix -DENABLE_PYTHON=OFF -DENABLE_NETCDF=OFF -DENABLE_AEC=ON -DCMAKE_C_COMPILER=$CC -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_Fortran_COMPILE_OPTIONS_PIE="${CMAKE_extra_fortran_options}" -DCMAKE_Fortran_FLAGS="${CMAKE_extra_fortran_flags}" -DENABLE_EXAMPLES=OFF
    make -j 8 install
    make clean
    cd $src_dir
    touch  eccodes-2.16.0.dep
fi

if [[ ! -e nco-4.9.1.dep ]]
then
    wget https://github.com/nco/nco/archive/4.9.1.tar.gz
    mv 4.9.1.tar.gz nco-4.9.1.tar.gz
    tar xf  nco-4.9.1.tar.gz
    cd nco-4.9.1
    CPPFLAGS="-I$prefix/include" \
    CFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib -lnetcdf" \
    CXXFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib -lnetcdf" \
    ./configure --prefix=$prefix --disable-dap --disable-ncap2 --disable-udunits --disable--udunits2
    make -j 4 install
    make clean
    cd $src_dir
    touch nco-4.9.1.dep
fi

if [[ ! -e sqlite-3310100.dep ]]
then
    wget https://www.sqlite.org/2020/sqlite-autoconf-3310100.tar.gz
    tar xvf sqlite-autoconf-3310100.tar.gz
    cd sqlite-autoconf-3310100
    ./configure --prefix=$prefix
    make -j 1 install
    make clean
    cd $src_dir
    touch sqlite-3310100.dep
fi

if [[ ! -e proj-6.3.0.dep ]]
then
    wget https://download.osgeo.org/proj/proj-6.3.0.tar.gz
    tar xvf proj-6.3.0.tar.gz
    cd proj-6.3.0
    CXXFLAGS="-std=c++11" \
    PKG_CONFIG_PATH="$prefix/lib/pkgconfig" \
    ./configure --prefix=$prefix
    make -j 4 install
    make clean
    cd $src_dir
    touch proj-6.3.0.dep
fi

if [[ ! -e cdi-1.9.7.dep ]]
then
    wget https://code.mpimet.mpg.de/attachments/download/20126/cdi-1.9.7.tar.gz
    tar xf cdi-1.9.7.tar.gz
    cd cdi-1.9.7
    CPPFLAGS="-I$prefix/include" \
    CFLAGS="-g -O3 -march=native" \
    LDFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib -lnetcdf" \
    ./configure --enable-iso-c-interface --with-netcdf=$prefix --with-eccodes=$prefix --prefix=$prefix
    make -j 4 install
    make clean
    cd $src_dir
    touch cdi-1.9.7.dep
fi

if [[ ! -e cdo-1.9.8.dep ]]
then
    rm -rf cdo*
    wget https://code.mpimet.mpg.de/attachments/download/20826/cdo-1.9.8.tar.gz
    tar xf  cdo-1.9.8.tar.gz
    cd cdo-1.9.8
    CPPFLAGS="-I$prefix/include" \
    CFLAGS="-g -O3 -march=native" \
    CXXFLAGS="-g -O3 -march=native" \
    LDFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib -lnetcdf" \
    ./configure --prefix=$prefix \
    --with-eccodes=$prefix \
    --with-hdf5=$prefix \
    --with-netcdf=$prefix \
    --with-proj=$prefix \
    --with-szlib=$prefix \
    --with-ossp-uuid \
    --with-curl
    make -j 4 install
    make clean
    cd $src_dir
    touch cdo-1.9.8.dep
fi
