set -x 


dir_run="/home/faroux/STRATO"
exec_new="_new"
exec_old="_old"

cd $dir_run

cd TESTS
#./get_exec.sh
#./make_exec.sh 0
cd ..

mkdir RES_NEW
mkdir RES_OLD

cp -f TESTS/_script*.sh .
cp -f TESTS/script*.sh .

cp -f TESTS/PC/*.exe .
mv -f python.exe TESTS/PYTHON

ln -s TESTS/PGD/FILES/* .

./script_pgd.sh $dir_run $exec_new $exec_old

./script_prep.sh $dir_run $exec_new $exec_old

./script_sea.sh $dir_run $exec_new $exec_old
./script_water.sh $dir_run $exec_new $exec_old
./script_flake.sh $dir_run $exec_new $exec_old
./script_isba.sh $dir_run $exec_new $exec_old
./script_teb.sh $dir_run $exec_new $exec_old

./script_csts.sh $dir_run $exec_new $exec_old
./script_offline.sh $dir_run $exec_new $exec_old

./script_cas_parts.sh $dir_run $exec_new $exec_old

mkdir NAMELIST

mv -f OPTIONS.nam_* NAMELIST


for file in TESTS/script*.sh
do
	file2=$(basename $file)
	rm -f $file2
done

for file in TESTS/PC/*.exe
do
	file2=$(basename $file)
	rm -f $file2
done

for file in TESTS/PGD/FILES/*
do
	file2=$(basename $file)
	rm -f $file2
done

for file in TESTS/PREP/FILES/*
do
	file2=$(basename $file)
	rm -f $file2
done

for file in TESTS/PREP/make*.sh
do
	file2=$(basename $file)
	rm -f $file2
done

for file in TESTS/PGD/make*.sh
do
	file2=$(basename $file)
	rm -f $file2
done


rm -f PGD*
rm -f PREP*
rm -f LISTING*


rm -f log*

rm -f Forc*.txt

rm -f Params_config* 

rm -f LISTING*
