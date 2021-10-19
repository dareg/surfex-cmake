#!/bin/bash

# Script to perform the entire SURFEX compilation procedure

curdir=$PWD

export ARCH=${ARCH=LXgfortran}

export VER_MPI=${VER_MPI=NOMPI}

#export OPTLEVEL=${OPTLEVEL=DEBUG}
export OPTLEVEL=${OPTLEVEL=O2}

export VER_USER=${VER_USER=""}
#export VER_USER=${VER_USER="ASSIM_Sodankyla"}
#export VER_USER=${VER_USER="MYSRC"}


export XYZ=""
source configure
export SRC_SURFEX=""

cd $curdir
source ../conf/profile_surfex-${XYZ}

cd $curdir

if [[ $VER_USER = "" ]]; then
  make
  make installmaster
else
  make user
  make installuser
fi

exit
