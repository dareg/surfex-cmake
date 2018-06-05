#set -x

exec_new=$3
exec_old=$4
nens=$5

echo $1 >> $2

cp -f OPTIONS.nam OPTIONS.nam_new

./script_base.sh "NEW" "$2"

echo "NEW" >> $2

rm -f BGROUND*

./script_exec_soda_enkf0.sh $exec_new $2 "RES_NEW" $nens

./script_test_end.sh "NEW" "$2" "$1"

mv -f *.OUT.nc RES_NEW
mv -f PGD.nc RES_NEW
mv -f PREP*.nc RES_NEW
#mv -f SURFOUT.nc RES_NEW

rm -f LISTING*

rm -f ANAL_INCR*
rm -f OBS*
rm -f INNOV*
rm -f K1*
rm -f F_AGG*
rm -f A_AGG*

./script_to_old.sh 

./script_base.sh "OLD" "$2"

echo "OLD" >> $2

rm -f BGROUND*

./script_exec_soda_enkf0.sh $exec_old $2 "RES_OLD" $nens

./script_test_end.sh "OLD" "$2" "$1"

mv -f *.OUT.nc RES_OLD
mv -f PGD.nc RES_OLD
mv -f PREP*.nc RES_OLD
#mv -f SURFOUT.nc RES_OLD

rm -f LISTING*

rm -f ANAL_INCR*
rm -f OBS*
rm -f INNOV*
rm -f K1*
rm -f F_AGG*
rm -f A_AGG*

cp -f OPTIONS.nam_new OPTIONS.nam

cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" $1

cd ../..

#rm -f RES_NEW/*.OUT.nc
#rm -f RES_NEW/PGD.*
#rm -f RES_NEW/SURFOUT.*
#rm -f RES_OLD/*.OUT.nc
#rm -f RES_OLD/PGD.*
#rm -f RES_OLD/SURFOUT.*

AAAAMMJJRR=2007071006
AAAAMMJJRR_end=2007071106
while [ $AAAAMMJJRR  -le $AAAAMMJJRR_end ]; do

  AAAAMMJJRRobs=`TESTS/CAS_PART/SODA/SCRIPTS/smsdate $AAAAMMJJRR 24`
  AAAAMMJJRR=$AAAAMMJJRRobs
  aa=`echo $AAAAMMJJRR | cut -c3-4`
  mm=`echo $AAAAMMJJRR | cut -c5-6`
  jj=`echo $AAAAMMJJRR | cut -c7-8`
  RR=`echo $AAAAMMJJRR | cut -c9-10` 

  mv -f RES_NEW/SURFOUT1_$aa$mm${jj}.nc RES_NEW/SURFOUT.nc
  mv -f RES_OLD/SURFOUT1_$aa$mm${jj}.nc RES_OLD/SURFOUT.nc

  cd TESTS/PYTHON
  ./python.exe "RES_NEW/" "RES_OLD/" $1_$AAAAMMJJRR
  cd ../..

done

rm -f ANAL_INCR*
rm -f OBS*
rm -f INNOV*
rm -f K1*
rm -f F_AGG*
rm -f A_AGG*



rm -f RES_NEW/*
rm -f RES_OLD/*

