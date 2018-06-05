
exec_new=$3
exec_old=$4

echo $1 >> $2

cp -f OPTIONS.nam OPTIONS.nam_new

./script_base.sh "NEW" "$2"

echo "NEW" >> $2
./pgd.exe $exec_new >> $2
./prep.exe $exec_new >> $2
./offline.exe $exec_new >> $2

mv -f PREP.nc PREP_SODA.nc
mv -f SURFOUT.nc PREP.nc

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LASSIM = F/LASSIM = T/g" OPTIONS.nam_save > OPTIONS.nam

./soda.exe $exec_new >> $2

./script_test_end.sh "NEW" "$2" "$1"

mv -f *.OUT.nc RES_NEW
mv -f PGD.nc RES_NEW
mv -f PREP.nc RES_NEW
mv -f SURFOUT*.nc RES_NEW

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LASSIM = T/LASSIM = F/g" OPTIONS.nam_save > OPTIONS.nam

./script_to_old.sh 

./script_base.sh "OLD" "$2"

echo "OLD" >> $2
./pgd.exe $exec_old >> $2
./prep.exe $exec_old >> $2
./offline.exe $exec_old >> $2

mv -f PREP.nc PREP_SODA.nc
mv -f SURFOUT.nc PREP.nc

rm -f LISTING*


cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LASSIM = F/LASSIM = T/g" OPTIONS.nam_save > OPTIONS.nam

./soda.exe $exec_old >> $2


./script_test_end.sh "OLD" "$2" "$1"

mv -f *.OUT.nc RES_OLD
mv -f PGD.* RES_OLD
mv -f PREP.* RES_OLD
mv -f SURFOUT.* RES_OLD


cp -f OPTIONS.nam_new OPTIONS.nam

cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" $1

cd ../..

rm -f RES_NEW/*
rm -f RES_OLD/*

