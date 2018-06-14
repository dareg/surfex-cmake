set -x 

rm -f $tmpdir/* 

dir_run="/home/cnrm_other/cen/mrns/nheilir/SURFEX/Surfex_V81/STRATO/"
exec_new="_new"
exec_old="_old"

cd $dir_run

mkdir RES_NEW
mkdir RES_OLD

cp -f TESTS/script*.sh .

cp -f TESTS/BEAUFIX/*.exe .
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
	rm $file2
done

for file in TESTS/PC/*.exe
do
	file2=$(basename $file)
	rm $file2
done

for file in TESTS/PGD/FILES/*
do
	file2=$(basename $file)
	rm $file2
done

