
exec_new=$3
exec_old=$4

echo $1 >> $2

cp -f OPTIONS.nam OPTIONS.nam_new

./script_restart.sh "NEW" "$2"

echo "NEW" >> $2
./offline.exe $exec_new >> $2

./script_test_end.sh "NEW" "$2" "$1"

cp -f PGD.* RES_NEW
cp -f PREP.* RES_NEW
mv -f *.OUT.nc RES_NEW
mv -f SURFOUT.* RES_NEW

./script_to_old.sh 

./script_restart.sh "OLD" "$2"

echo "OLD" >> $2
./offline.exe $exec_old >> $2

./script_test_end.sh "OLD" "$2" "$1"

mv -f *.OUT.nc RES_OLD
mv -f SURFOUT.* RES_OLD
mv -f PGD.* RES_OLD
mv -f PREP.* RES_OLD



cp -f OPTIONS.nam_new OPTIONS.nam

cd TESTS/PYTHON

./python.exe "RES_NEW/" "RES_OLD/" $1

cd ../..

rm -f RES_NEW/*
rm -f RES_OLD/*

