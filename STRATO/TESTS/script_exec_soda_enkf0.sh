#!/bin/ksh
#set -x

exec1=$1
repres=$3
nens=$4


./pgd.exe $exec1 >> $2
./prep.exe $exec1 >> $2

AAAAMMJJRR_deb=2007071006

AAAAMMJJRR_end=2007071206

AAAAMMJJRR=$AAAAMMJJRR_deb
aa=`echo $AAAAMMJJRR | cut -c3-4`
mm=`echo $AAAAMMJJRR | cut -c5-6`
jj=`echo $AAAAMMJJRR | cut -c7-8`
RR=`echo $AAAAMMJJRR | cut -c9-10`

while [ $AAAAMMJJRR  -le $AAAAMMJJRR_end ]; do

  echo $AAAAMMJJRR >> $2

  cp -f TESTS/CAS_PART/SODA/FORC/'20'$aa$mm$jj$RR'_FORCING_7_3.nc' FORCING.nc
  cp -f TESTS/FORCAGES/HIVER/Params_config.txt_soda_ekf$jj Params_config.txt

  AAAAMMJJRR_obs=$AAAAMMJJRR

  AAAAMMJJRR=`TESTS/CAS_PART/SODA/SCRIPTS/smsdate $AAAAMMJJRR_obs 24`
  aa=`echo $AAAAMMJJRR | cut -c3-4`
  mm=`echo $AAAAMMJJRR | cut -c5-6`
  jj=`echo $AAAAMMJJRR | cut -c7-8`
  RR=`echo $AAAAMMJJRR | cut -c9-10`  

# Run control + perturbed runs
  ens=1
  ens2=0
  while [ $ens -le $nens ]; do
 
    echo "####### SODA_EKF_ENS"$ens >> $2

    	cp -f OPTIONS.nam OPTIONS.nam_save
    	sed -e "s/NIE = $ens2/NIE = $ens/g" OPTIONS.nam_save > OPTIONS.nam

	if [ $AAAAMMJJRR_obs  -gt $AAAAMMJJRR_deb ]
	then
		cp -f OPTIONS.nam OPTIONS.nam_save
    	  	sed -e "s/LENS_GEN = T/LENS_GEN = F/g" OPTIONS.nam_save > OPTIONS.nam
	  	cp -f SURFOUT${ens}.nc PREP.nc
	fi

    ./offline.exe $exec1 >> $2
    mv -f SURFOUT.nc PREP_${aa}${mm}${jj}H${RR}_EKF_ENS${ens}.nc

    ens2=$ens
    ens=$(( $ens + 1 ))

  done

  if [ -f TESTS/CAS_PART/SODA/OBSDATA/CANARI_NATURE_$aa$mm${jj}H${RR}.DAT ]
  then
    ln -s TESTS/CAS_PART/SODA/OBSDATA/CANARI_NATURE_$aa$mm${jj}H${RR}.DAT OBSERVATIONS_$aa$mm${jj}H${RR}.DAT
  else
    ln -s TESTS/CAS_PART/SODA/OBSDATA/CANARI_NATURE_NULL.DAT OBSERVATIONS_$aa$mm${jj}H${RR}.DAT
  fi
  ln -s OBSERVATIONS_$aa$mm${jj}H${RR}.DAT CANARI_NATURE_$aa$mm${jj}H${RR}.DAT

  cp -f PREP.nc $repres
  cp -f PREP_${aa}${mm}${jj}H${RR}_EKF_ENS1.nc PREP.nc

  cp -f OPTIONS.nam OPTIONS.nam_save
  sed -e "s/NIE = $nens/NIE = 0/g" OPTIONS.nam_save > OPTIONS.nam
  cp -f OPTIONS.nam OPTIONS.nam_save
  sed -e "s/LASSIM = F/LASSIM = T/g" OPTIONS.nam_save > OPTIONS.nam

  echo "####### SODA_ASSIM " >> $2
 
  ./soda.exe $exec1 >> $2
  
  cp -f OPTIONS.nam OPTIONS.nam_save
  sed -e "s/LASSIM = T/LASSIM = F/g" OPTIONS.nam_save > OPTIONS.nam
  rm -f CANARI_NATURE_$aa$mm${jj}H${RR}.DAT
  rm -f OBSERVATIONS_$aa$mm${jj}H${RR}.DAT
  
  for file in SURFOUT*.nc
  do
	  file2=$(basename $file .nc)
	  cp -f $file $repres/${file2}_${aa}${mm}${jj}.nc
  done
	

done

cp -f PREP.nc $repres/SURFOUT.nc
#
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LENS_GEN = F/LENS_GEN = T/g" OPTIONS.nam_save > OPTIONS.nam
