#!/bin/bash

ROOT=$PWD

module load intel/oneapi/2023.2
module load compiler/2023.2.0
module load intelmpi/2018.5.274
module load gcc/9.2.0

export FC=ifort
export CC=icc
export CXX=icpc

if [ ! -d surfex-cmake ]; then
	git clone https://github.com/dareg/surfex-cmake
fi
if [ ! -d FALFILFA ]; then
	git clone https://github.com/ACCORD-NWP/FALFILFA
fi
if [ ! -d eccodes ]; then
	git clone https://github.com/ecmwf/eccodes --branch 2.38.3
fi
if [ ! -d fiat ]; then
	git clone https://github.com/ecmwf-ifs/fiat/
fi
if [ ! -d pinuts ]; then
	git clone https://github.com/ACCORD-NWP/pinuts
fi

cd eccodes || exit 1
mkdir -p build
cd build || exit 1
cmake .. -DCMAKE_INSTALL_PREFIX=$ROOT/install -DENABLE_AEC=OFF
make install -j 16 || exit 1
cd $ROOT || exit 1

cd fiat
mkdir -p build
cd build || exit 1
cmake .. -DCMAKE_INSTALL_PREFIX=$ROOT/install
make install -j || exit 1
cd $ROOT || exit 1

cd FALFILFA
mkdir -p build
cd build || exit 1
cmake .. -DCMAKE_INSTALL_PREFIX=$ROOT/install -Deccodes_ROOT=$ROOT/install \
	-DENABLE_SINGLE_PRECISION=ON \
	-DENABLE_DOUBLE_PRECISION=ON
make install -j 16 || exit 1
cd $ROOT || exit 1

#cd pinuts
#mkdir -p build
#cd build || exit 1
#cmake .. -DCMAKE_INSTALL_PREFIX=$ROOT/install \
#	-DENABLE_SINGLE_PRECISION=ON \
#	-DENABLE_DOUBLE_PRECISION=ON
#make install -j 16 || exit 1
#cd $ROOT || exit 1

cd surfex-cmake
mkdir -p build
cd build || exit 1
cmake .. -DCMAKE_INSTALL_PREFIX=$ROOT/install \
	-Dfalfilfa_ROOT=$ROOT/install \
	-DNetCDF_DIR=/home/gmap/mrpm/khatib/opt/i-2018.5.274/netcdf-4.7.1/ \
	-Deccodes_ROOT=/home/gmap/mrpm/khatib/opt/i-2018.5.274/eccodes-2.27.0/ \
	-DENABLE_SINGLE_PRECISION=ON \
	-DENABLE_DOUBLE_PRECISION=ON
make install -j 16 || exit 1
cd $ROOT || exit 1

