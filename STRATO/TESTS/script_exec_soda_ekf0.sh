#!/bin/ksh
#set -x
npert=2

exec1=$1
repres=$3


./pgd.exe $exec1 >> $2
./prep.exe $exec1 >> $2

AAAAMMJJRR=2007071006
aa=`echo $AAAAMMJJRR | cut -c3-4`
mm=`echo $AAAAMMJJRR | cut -c5-6`
jj=`echo $AAAAMMJJRR | cut -c7-8`
RR=`echo $AAAAMMJJRR | cut -c9-10`

AAAAMMJJRR_end=2007071106

while [ $AAAAMMJJRR  -le $AAAAMMJJRR_end ]; do

  echo $AAAAMMJJRR >> $2

  cp -f TESTS/CAS_PART/SODA/FORC/'20'$aa$mm$jj$RR'_FORCING_7_3.nc' FORCING.nc
  cp -f TESTS/FORCAGES/HIVER/Params_config.txt_soda_ekf$jj Params_config.txt

  AAAAMMJJRRobs=`TESTS/CAS_PART/SODA/SCRIPTS/smsdate $AAAAMMJJRR 24`
  AAAAMMJJRR=$AAAAMMJJRRobs
  aa=`echo $AAAAMMJJRR | cut -c3-4`
  mm=`echo $AAAAMMJJRR | cut -c5-6`
  jj=`echo $AAAAMMJJRR | cut -c7-8`
  RR=`echo $AAAAMMJJRR | cut -c9-10`  

# Run control + perturbed runs
  pert=0
  pert2=1
  while [ $pert -le $npert ]; do
 
    echo "####### SODA_EKF_PERT"$pert >> $2

    if [ "$pert" -gt 0 ]; then
    	cp -f OPTIONS.nam OPTIONS.nam_save
    	sed -e "s/NIVAR = $pert2/NIVAR = $pert/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/LPRT = F/LPRT = T/g" OPTIONS.nam_save > OPTIONS.nam
    fi

    ./offline.exe $exec1 >> $2
    mv -f SURFOUT.nc PREP_${aa}${mm}${jj}H${RR}_EKF_PERT${pert}.nc
    if [ "$pert" -eq 0 ]; then
      mv -f *.OUT.nc $repres
    fi

    pert2=$pert
    pert=$(( $pert + 1 ))

  done

  mv -f PREP.nc PREP_INIT.nc

  if [ -f TESTS/CAS_PART/SODA/OBSDATA/CANARI_NATURE_$aa$mm${jj}H${RR}.DAT ]
  then
    ln -sf TESTS/CAS_PART/SODA/OBSDATA/CANARI_NATURE_$aa$mm${jj}H${RR}.DAT OBSERVATIONS_$aa$mm${jj}H${RR}.DAT
  else
    ln -sf TESTS/CAS_PART/SODA/OBSDATA/CANARI_NATURE_NULL.DAT OBSERVATIONS_$aa$mm${jj}H${RR}.DAT
  fi
  ln -s OBSERVATIONS_$aa$mm${jj}H${RR}.DAT CANARI_NATURE_$aa$mm${jj}H${RR}.DAT

  cp -f OPTIONS.nam OPTIONS.nam_save
  sed -e "s/LPRT = T/LPRT = F/g" OPTIONS.nam_save > OPTIONS.nam
  cp -f OPTIONS.nam OPTIONS.nam_save
  sed -e "s/NIVAR = $npert/NIVAR = 1/g" OPTIONS.nam_save > OPTIONS.nam
  cp -f OPTIONS.nam OPTIONS.nam_save
  sed -e "s/LASSIM = F/LASSIM = T/g" OPTIONS.nam_save > OPTIONS.nam
  cp -f OPTIONS.nam OPTIONS.nam_save
  sed -e "s/CPREPFILE = \"PREP\"/CPREPFILE = \"PREP_${aa}${mm}${jj}H${RR}_EKF_PERT0\"/g" OPTIONS.nam_save > OPTIONS.nam

  echo "####### SODA_ASSIM " >> $2
 
  ./soda.exe $exec1 >> $2
  
  cp -f OPTIONS.nam OPTIONS.nam_save
  sed -e "s/LASSIM = T/LASSIM = F/g" OPTIONS.nam_save > OPTIONS.nam
  cp -f OPTIONS.nam OPTIONS.nam_save
  sed -e "s/CPREPFILE = \"PREP_${aa}${mm}${jj}H${RR}_EKF_PERT0\"/CPREPFILE = \"PREP\"/g" OPTIONS.nam_save > OPTIONS.nam
  rm -f CANARI_NATURE_$aa$mm${jj}H${RR}.DAT
  rm -f OBSERVATIONS_$aa$mm${jj}H${RR}.DAT
  
  cp -f SURFOUT.nc $repres/SURFOUT_$aa$mm${jj}.nc

  mv -f SURFOUT.nc PREP.nc
  mv -f BGROUNDout_ASSIM.0000001 BGROUNDin.0000001

done

cp -f PREP.nc $repres/SURFOUT.nc

