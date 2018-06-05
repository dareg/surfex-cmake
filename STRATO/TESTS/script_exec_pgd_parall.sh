#set -x 

exec_new=$3
exec_old=$4

cp -f OPTIONS.nam OPTIONS.nam_new

echo $1 >> $2

################################

./script_to_old.sh 
cp -f OPTIONS.nam OPTIONS.nam_old
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


#1 thread
echo "SIMPLE" >> $2
export OMP_NUM_THREADS=1
./pgd.exe $exec_new >> $2
./prep.exe $exec_new >> $2
#add test
#mv -f PGD.* RES_NEW
#mv -f PREP.* RES_NEW
#./script_to_old.sh
#./pgd.exe $exec_old >> $2
#./prep.exe $exec_old >> $2
#add test
./offline.exe $exec_new >> $2

./script_test_end.sh "NEW_SIMPLE" "$2" "$1"

mv -f *.OUT.nc RES_NEW
mv -f PGD.nc RES_NEW
mv -f PREP*.nc RES_NEW
mv -f SURFOUT*.nc RES_NEW

cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" "$1_SIMPLE"

cd ../..
#rm -f RES_NEW/*


#4 tasks
echo "MPI4" >> $2
./pgd_mpi.exe $exec_new >> $2
./prep_mpi.exe $exec_new >> $2
./offline.exe $exec_new >> $2

./script_test_end.sh "NEW_MPI4" "$2" "$1"

mv -f *.OUT.nc RES_NEW
mv -f PGD.nc RES_NEW
mv -f PREP*.nc RES_NEW
mv -f SURFOUT*.nc RES_NEW

cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" "$1_MPI4"

cd ../..
#rm -f RES_NEW/*

#2 tasks 2 threads
echo "MPI2OMP2" >> $2
export OMP_NUM_THREADS=2
./pgd_mpiomp.exe $exec_new >> $2
./prep_mpiomp.exe $exec_new >> $2
./offline.exe $exec_new >> $2

./script_test_end.sh "NEW_MPI2OMP2" "$2" "$1"

mv -f *.OUT.nc RES_NEW
mv -f PGD.nc RES_NEW
mv -f PREP*.nc RES_NEW
mv -f SURFOUT*.nc RES_NEW

cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" "$1_MPI2OMP2"

cd ../..
#rm -f RES_NEW/*


#4 thread
#echo "OMP4" >> $2
#export OMP_NUM_THREADS=4
#./pgd_omp.exe $exec_new >> $2
#./prep_omp.exe $exec_new >> $2
#./offline.exe $exec_new >> $2
#
#./script_test_end.sh "NEW_OMP4" "$2" "$1"
#
#mv -f *.OUT.nc RES_NEW
#mv -f PGD.nc RES_NEW
#mv -f PREP*.nc RES_NEW
#mv -f SURFOUT*.nc RES_NEW
#
#
#export OMP_NUM_THREADS=1


#cd TESTS/PYTHON
#
#./python.exe "RES_NEW/" "RES_OLD/" "$1_OMP4"

#cd ../..
rm -f RES_NEW/*
rm -f RES_OLD/*

