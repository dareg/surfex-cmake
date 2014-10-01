#!/bin/bash

export ARCH=${ARCH=LXgfortran}
export VER_MPI=${VER_MPI=NOMPI}
export OPTLEVEL=${OPTLEVEL=O2}
export VER_USER=${VER_USER=""}

export XYZ=""
source configure
export SRC_SURFEX=""
cd ..
source ../conf/profile_surfex-${XYZ}

if [[ $VER_USER = "" ]]; then
  make
  make installmaster
else
  make user
  make installuser
fi

exit
