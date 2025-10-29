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
export INSTALL_DIR_ECCODES=$ROOT/install/eccodes
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR_ECCODES -DENABLE_AEC=OFF
make install -j 16 || exit 1
cd $ROOT || exit 1
export eccodes_ROOT=$INSTALL_DIR_ECCODES

cd fiat
mkdir -p build
cd build || exit 1
export INSTALL_DIR_FIAT=$ROOT/install/fiat
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR_FIAT
make install -j || exit 1
cd $ROOT || exit 1
export fiat_ROOT=$INSTALL_DIR_FIAT

cd FALFILFA
mkdir -p build
cd build || exit 1
export INSTALL_DIR_FALFILFA=$ROOT/install/falfilfa
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR_FALFILFA \
	-DENABLE_SINGLE_PRECISION=OFF \
	-DENABLE_DOUBLE_PRECISION=ON
make install -j 16 || exit 1
cd $ROOT || exit 1
export falfilfa_ROOT=$INSTALL_DIR_FALFILFA

cd pinuts
mkdir -p build
cd build || exit 1
export INSTALL_DIR_PINUTS=$ROOT/install/pinuts
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR_PINUTS \
	-DENABLE_SINGLE_PRECISION=OFF \
	-DENABLE_DOUBLE_PRECISION=ON
make install || exit 1 #Oct. 2025 - pinuts doesnt compile in parallel
cd $ROOT || exit 1
export pinuts_ROOT=$INSTALL_DIR_PINUTS

cd surfex-cmake
mkdir -p build
cd build || exit 1
export INSTALL_DIR_SURFEX=$ROOT/install/surfex
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR_SURFEX \
	-DNetCDF_DIR=/home/gmap/mrpm/khatib/opt/i-2018.5.274/netcdf-4.7.1/ \
	-DENABLE_SINGLE_PRECISION=OFF \
	-DENABLE_DOUBLE_PRECISION=ON
make install -j 16 || exit 1
cd $ROOT || exit 1

