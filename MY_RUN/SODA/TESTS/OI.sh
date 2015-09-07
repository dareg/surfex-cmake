#!/bin/bash

# Some sanity checks
if [ "$SRC_SURFEX" == "" -o "$XYZ" == "" ]; then
  echo "You need to source the config file before calling this script!"
  exit 1
else
  dtg=$1
  suffix=$2
  nproc=$3
  extra=$4
  [ ! -d $SRC_SURFEX/MY_RUN/SODA ] && echo "The directory $SRC_SURFEX/MY_RUN/SODA was not found" && exit 1
  [ ! -f $SRC_SURFEX/exe/SODA$XYZ ] && echo "The needed binary $SRC_SURFEX/exe/SODA$XYZ was not found" && exit 1
fi

fcint=06
if [ "$dtg" == "2008093018" ]; then
  forc_format="ASCII"
else
  forc_format="NETCDF"
fi

dtgtmp=`echo $dtg | cut -c1-4`"-"`echo $dtg | cut -c5-6`"-"`echo $dtg | cut -c7-8`" "`echo $dtg | cut -c9-10`":00"
dtgtmp=`date -u --date "$dtgtmp $fcint hours" +%Y%m%d%H`
YYYY_RES=`echo $dtgtmp | cut -c1-4`
YY_RES=`echo $dtgtmp | cut -c3-4`
MM_RES=`echo $dtgtmp | cut -c5-6`
DD_RES=`echo $dtgtmp | cut -c7-8`
HH_RES=`echo $dtgtmp | cut -c9-10`

# Run offline
export OMP_NUM_THREADS=1
offline="$SRC_SURFEX/MY_RUN/SODA/TESTS/OFFLINE.sh"
if [ -f $offline ]; then
  $offline $dtg $fcint -1 $suffix $forc_format $nproc $extra || exit 1
else
  echo "No script $offline found"
  exit 1
fi

work=$SRC_SURFEX/MY_RUN/SODA/RUN/RUN_OI
[ -d $work ] || mkdir $work
cd $work || exit 1

# Copy namelists
cat $SRC_SURFEX/MY_RUN/SODA/NAMELISTS/OPTIONS.nam_$suffix $SRC_SURFEX/MY_RUN/SODA/NAMELISTS/OPTIONS.nam_soda_OI$extra > OPTIONS.nam.tmp
sed -e "s/CFORCING_FILETYPE=CFORCING_FILETYPE/CFORCING_FILETYPE=\"$forc_format\"/" OPTIONS.nam.tmp > OPTIONS.nam

# Ecoclimap
ln -sf $SRC_SURFEX/MY_RUN/ECOCLIMAP/*.bin .
ln -sf $SRC_SURFEX/MY_RUN/ECOCLIMAP/*.dat .

# Copy PGD file
cp  $SRC_SURFEX/MY_RUN/SODA/OUTPUT/PGD_SODA.$suffix .

# Copy OFFLINE PREP file
cp  $SRC_SURFEX/MY_RUN/SODA/OUTPUT/PREP_OFFLINE.$suffix PREP_SODA.$suffix

# Copy input grib file
cp $SRC_SURFEX/MY_RUN/SODA/INPUT/FG_OI_MAIN .
cp $SRC_SURFEX/MY_RUN/SODA/INPUT/CANARI* .
cp $SRC_SURFEX/MY_RUN/SODA/INPUT/OBSERVATIONS_* .
cp $SRC_SURFEX/MY_RUN/SODA/INPUT/FIRST_GUESS_* .
cp $SRC_SURFEX/MY_RUN/SODA/INPUT/SST_SIC* .
cp $SRC_SURFEX/MY_RUN/SODA/INPUT/fort.61 .
cp $SRC_SURFEX/MY_RUN/SODA/INPUT/clim_isba .
cp $SRC_SURFEX/MY_RUN/SODA/INPUT/FIRST_GUESS* .
ln -sf $SRC_SURFEX/MY_RUN/SODA/INPUT/CLIMATE.DAT_$dtg CLIMATE.DAT
ln -sf $SRC_SURFEX/MY_RUN/SODA/INPUT/LSM.DAT_$dtg LSM.DAT

# Run SODA
[ ! -f $SRC_SURFEX/exe/SODA$XYZ ] && echo "$SRC_SURFEX/exe/SODA$XYZ not found" && exit 1
mpirun -np $nproc $SRC_SURFEX/exe/SODA$XYZ || exit 1

# Save output
mv SURFOUT.${YYYY_RES}${MM_RES}${DD_RES}_${HH_RES}h00.$suffix $SRC_SURFEX/MY_RUN/SODA/OUTPUT/PREP_SODA_OI.$suffix || exit 1
