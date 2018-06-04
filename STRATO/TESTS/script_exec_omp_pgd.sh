
exec_new=$3
exec_old=$4

cp -f OPTIONS.nam OPTIONS.nam_new

echo $1 >> $2

################################

./script_to_old.sh 

./script_base.sh "OLD" "$2"

echo "OLD" >> $2
./pgd.exe $exec_old >> $2
./prep.exe $exec_old >> $2
./offline.exe $exec_old >> $2

./script_test_end.sh "OLD" "$2" "$1"

mv -f *.OUT.nc RES_OLD
mv -f PGD.* RES_OLD
mv -f PREP.* RES_OLD
mv -f SURFOUT.* RES_OLD

################################

cp -f OPTIONS.nam_new OPTIONS.nam

./script_base.sh "NEW" "$2"

echo "NEW" >> $2

export OMP_NUM_THREADS=1
./pgd.exe $exec_new >> $2
mv -f PGD.nc PGD_1.nc

export OMP_NUM_THREADS=4
./pgd_omp.exe $exec_new >> $2
mv -f PGD.nc PGD_4.nc

export OMP_NUM_THREADS=1


#1 thread
mv -f PGD_1.nc PGD.nc
./prep.exe $exec_new >> $2
./offline.exe $exec_new >> $2

./script_test_end.sh "NEW" "$2" "$1"

mv -f *.OUT.nc RES_NEW
mv -f PGD.nc RES_NEW
mv -f PREP*.nc RES_NEW
mv -f SURFOUT*.nc RES_NEW

cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" "$1_1THREAD"

cd ../..
rm -f RES_NEW/*


#4 thread
mv -f PGD_4.nc PGD.nc
./prep.exe $exec_new >> $2
./offline.exe $exec_new >> $2

./script_test_end.sh "NEW" "$2" "$1"

mv -f *.OUT.nc RES_NEW
mv -f PGD.nc RES_NEW
mv -f PREP*.nc RES_NEW
mv -f SURFOUT*.nc RES_NEW

cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" "$1_4THREAD"

cd ../..
rm -f RES_NEW/*
rm -f RES_OLD/*

