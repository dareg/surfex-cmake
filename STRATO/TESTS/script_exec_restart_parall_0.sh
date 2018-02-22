
exec_new=$3
exec_old=$4

cp -f OPTIONS.nam OPTIONS.nam_new

echo $1 >> $2


#################old en premier

./script_to_old.sh 

./script_restart.sh "OLD" "$2"

echo "OLD" >> $2
./offline.exe $exec_old >> $2

./script_test_end.sh "OLD" "$2" "$1"

cp -f PGD.* RES_OLD
cp -f PREP.* RES_OLD
mv -f *.OUT.nc RES_OLD
mv -f SURFOUT.* RES_OLD

#################new en second

cp -f OPTIONS.nam_new OPTIONS.nam

./script_restart.sh "NEW" "$2"

echo "NEW" >> $2

##########1: new simple
echo "SIMPLE" >> $2
./offline.exe $exec_new >> $2

./script_test_end.sh "NEW_SIMPLE" "$2" "$1"

mv -f *.OUT.nc RES_NEW
mv -f SURFOUT.* RES_NEW


cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" "$1_SIMPLE"

cd ../..

rm -f RES_NEW/*

##########1: new MPI4
echo "MPI4" >> $2
./offline_mpi.exe $exec_new >> $2

./script_test_end.sh "NEW_MPI4_LIN" "$2" "$1"

mv -f *.OUT.nc RES_NEW
mv -f SURFOUT.* RES_NEW


cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" "$1_MPI4"

cd ../..

rm -f RES_NEW/*

##########1: new OMP4
echo "OMP4" >> $2
export OMP_NUM_THREADS=4
./offline_omp.exe $exec_new >> $2

./script_test_end.sh "NEW_OMP4" "$2" "$1"

mv -f *.OUT.nc RES_NEW
mv -f SURFOUT.* RES_NEW


cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" "$1_OMP4"

cd ../..

rm -f RES_NEW/*

##########1: new MPI2OMP2
echo "MPI2OMP2" >> $2
export OMP_NUM_THREADS=2
./offline_mpiomp.exe $exec_new >> $2

./script_test_end.sh "NEW_MPI2OMP2" "$2" "$1"

mv -f PGD.* RES_NEW
mv -f PREP.* RES_NEW
mv -f *.OUT.nc RES_NEW
mv -f SURFOUT.* RES_NEW


cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" "$1_MPI2OMP2"

cd ../..

rm -f RES_NEW/*
rm -f RES_OLD/*

export OMP_NUM_THREADS=1

