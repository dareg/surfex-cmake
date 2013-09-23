#!/bin/bash

# Some sanity checks
if [ "$SRC_SURFEX" == "" -o "$XYZ" == "" ]; then
  echo "You need to source the config file before calling this script!"
  exit 1
else
  [ ! -d $SRC_SURFEX/MY_RUN/SODA ] && echo "The directory $SRC_SURFEX/MY_RUN/SODA was not found" && exit 1
  [ ! -f $SRC_SURFEX/exe/SODA$XYZ ] && echo "The needed binary $SRC_SURFEX/exe/SODA$XYZ was not found" && exit 1
fi

LPRT=".FALSE."
IVAR="1"
if [ $# -ne 0 -a $# -ne 1 ]; then
  echo "Usage: $0 [pert]"
else
  mbr=""
  if [ $# -eq 1 ]; then
    pert=$1
    mbr="_PERT_"$(printf "%.2d" "$pert")
    echo "Running offline for perturbation $pert"
    if [ "$pert" -gt 0 ]; then
      LPRT=".TRUE."
      IVAR="$pert"
    fi
  else
    echo "Running offline"
  fi
fi

work=$SRC_SURFEX/MY_RUN/SODA/RUN/RUN_OFFLINE$mbr
[ -d $work ] || mkdir $work
cd $work || exit 1

# Copy namelists
sed -e "s/LPRT=LPRT/LPRT=$LPRT/" \
    -e "s/IVAR=IVAR/IVAR=$IVAR/" \
    $SRC_SURFEX/MY_RUN/SODA/NAMELISTS/OPTIONS.nam_offline > OPTIONS.nam

# Ecoclimap
ln -sf $SRC_SURFEX/MY_RUN/ECOCLIMAP/*.bin .
ln -sf $SRC_SURFEX/MY_RUN/ECOCLIMAP/*.dat .

# Copy Forcing
cp $SRC_SURFEX/MY_RUN/SODA/INPUT/*.txt .

# Copy PGD file
cp  $SRC_SURFEX/MY_RUN/SODA/OUTPUT/PGD_SODA.lfi .

# Copy OFFLINE PREP file
cp  $SRC_SURFEX/MY_RUN/SODA/OUTPUT/PREP_SODA$mbr.lfi PREP_SODA.lfi

# Run SODA
[ ! -f $SRC_SURFEX/exe/OFFLINE$XYZ ] && echo "$SRC_SURFEX/exe/OFFLINE$XYZ not found" && exit 1
$SRC_SURFEX/exe/OFFLINE$XYZ

if [ -f SURFOUT.20081001_00h00.lfi ]; then
  mv SURFOUT.20081001_00h00.lfi $SRC_SURFEX/MY_RUN/SODA/OUTPUT/PREP_OFFLINE$mbr.lfi
else
  echo "Output file from SURFEX not found"
  exit 1
fi
