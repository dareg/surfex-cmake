dir_run=$1
exec_new=$2
exec_old=$3

fname="INFO_ISBA.txt"

rm -f $fname
touch $fname

cd $1

cp -f TESTS/script_exec.sh .
cp -f TESTS/script_to_old.sh .
cp -f TESTS/script_exec_parall.sh .

rm -f Forc_*.txt
rm -f Params_config.txt

cp -f TESTS/PREP/make*.sh .
cp -f TESTS/PGD/make*.sh .

ln -s TESTS/FORCAGES/HIVER/Params_config.txt .
ln -s TESTS/FORCAGES/HIVER/Forc*.txt .

####################NAMELISTS INITIALES##########################################

#namelists avec PREP uniformes
cp -f TESTS/OPTIONS.nam_000 OPTIONS.nam_PHYSIO_UNIF
./make_prep_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_isba_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_isba_snow_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_teb_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_gd_gr_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_teb_snow_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_flake_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1.0E+20/XSST_UNIF = 283.15/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/XTS_WATER_UNIF = 1.0E+20/XTS_WATER_UNIF = 282./g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF


cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YCLAY = \"\"/YCLAY = \"clay_fao\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YSAND = \"\"/YSAND = \"sand_fao\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YSOC_TOP = \"\"/YSOC_TOP = \"soc_top\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YSOC_SUB = \"\"/YSOC_SUB = \"soc_sub\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YCTI = \"\"/YCTI = \"topo_index\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YPERM = \"\"/YPERM = \"perm_glo_10km\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF


#recup d'un cas ECOCLIMAP 
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam
sed -e "s/LECOCLIMAP\ \=\ F/LECOCLIMAP\ \=\ T/g" OPTIONS.nam > OPTIONS.nam_ECOCLIMAP

###recup d'un cas sans ECOCLIMAP (pour la grille CARTESIAN)
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_NO_ECOCLIMAP
./make_eco_tiles.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_isba_unif.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_teb_unif1.sh "OPTIONS.nam_NO_ECOCLIMAP"


cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_base

mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/LSET\_FORC\_ZS\ \=\ F/LSET\_FORC\_ZS\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_base

####################CISBA##########################################

cp -f OPTIONS.nam_base OPTIONS.nam
#./script_isba_photo.sh "OPTIONS.nam" "ISBA_2L_PHOTO" $fname $2 $3
cp -f OPTIONS.nam_base OPTIONS.nam
#./script_isba_phys.sh "OPTIONS.nam" "ISBA_2L_PHYS" $fname $2 $3

cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CISBA\ =\ \"2-L\"/CISBA\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND\_LAYER\ =\ 2/NGROUND\_LAYER\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_isba
./script_isba_photo.sh "OPTIONS.nam" "ISBA_3L_PHOTO" $fname $2 $3
mv -f OPTIONS.nam_isba OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NPATCH\ =\ 1/NPATCH\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_isba
./script_isba_phys.sh "OPTIONS.nam" "ISBA_3L_SN3L_PHYS" $fname $2 $3

echo "################### ISBA_3L_SN3L_MEB "
cp -f OPTIONS.nam_isba OPTIONS.nam_save
sed -e "s/LMEB_PATCH = F\,F\,F\,F\,F\,F\,F\,F\,F\,F\,F\,F/LMEB_PATCH = F\,F\,F\,T\,T\,T\,F\,F\,F\,F\,F\,F/g" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ISBA_3L_SN3L_MEB
./script_exec_parall.sh "ISBA_3L_SN3L_MEB" $fname $2 $3
#
##test ARRANGE_COVER
echo "############### ISBA_3L_SN_3L_TOWN_TO_ROCK"
cp -f OPTIONS.nam_isba OPTIONS.nam_save
sed -e "s/LTOWN\_TO\_ROCK\ \=\ F/LTOWN\_TO\_ROCK\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ISBA_3L_SN_3L_TOWN_TO_ROCK
./script_exec.sh "ISBA_3L_SN_3L_TOWN_TO_ROCK" $fname $2 $3
echo "############### ISBA_3L_SN_3L_WATER_TO_NATURE"
cp -f OPTIONS.nam_isba OPTIONS.nam_save
sed -e "s/LWATER\_TO\_NATURE\ \=\ F/LWATER\_TO\_NATURE\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ISBA_3L_SN_3L_WATER_TO_NATURE
./script_exec.sh "ISBA_3L_SN_3L_WATER_TO_NATURE" $fname $2 $3
echo "############### ISBA_3L_SN_3L_TOWN_TO_ROCK_WATER_TO_NATURE"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LTOWN\_TO\_ROCK\ \=\ F/LTOWN\_TO\_ROCK\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ISBA_3L_SN_3L_TOWN_TO_ROCK_WATER_TO_NATURE
./script_exec_parall.sh "ISBA_3L_SN_3L_TOWN_TO_ROCK_WATER_TO_NATURE" $fname $2 $3


cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CISBA\ =\ \"2-L\"/CISBA\ =\ \"DIF\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND\_LAYER\ =\ 2/NGROUND\_LAYER\ =\ 8/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_dif
./script_isba_photo_dif.sh "OPTIONS.nam" "ISBA_DIF8" $fname $2 $3
cp -f OPTIONS.nam OPTIONS.nam_dif
./script_isba_phys.sh "OPTIONS.nam" "ISBA_DIF8" $fname $2 $3

echo "################### ISBA_DIF8_SN3L_MEB"
cp -f OPTIONS.nam_isba OPTIONS.nam_save
sed -e "s/LMEB_PATCH = F\,F\,F\,F\,F\,F\,F\,F\,F\,F\,F\,F/LMEB_PATCH = F\,F\,F\,T\,T\,T\,F\,F\,F\,F\,F\,F/g" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ISBA_DIF8_SN3L_MEB
./script_exec_parall.sh "ISBA_DIF8_SN3L_MEB" $fname $2 $3

echo "########### ISBA_DIF8_SN3L_CO84"
cp -f OPTIONS.nam_dif OPTIONS.nam_save
sed -e "s/CPEDO\_FUNCTION\ =\ \"CH78\"/CPEDO\_FUNCTION\ =\ \"CO84\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ISBA_DIF8_SN3L_CO84
./script_exec_parall.sh "ISBA_DIF8_SN3L_CO84" $fname $2 $3


#test avec DIF grille optimale
echo "########### ISBA_DIF8_SN3L_OPT_GRID"
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CISBA\ =\ \"2-L\"/CISBA\ =\ \"DIF\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND\_LAYER\ =\ 2/NGROUND\_LAYER\ =\ 14/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSOILGRID\ =\ 0.01,0.1,1.1,2.4,3.9,5.2,6.1,7.9//g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ISBA_DIF8_SN3L_OPT_GRID
./script_exec_parall.sh "ISBA_DIF8_SN3L_OPT_GRID" $fname $2 $3


for file in TESTS/PREP/make*.sh
do
	file2=$(basename $file)
	rm $file2
done

for file in TESTS/PGD/make*.sh
do
	file2=$(basename $file)
	rm $file2
done

rm -f Params_config.txt .
rm -f Forc*.txt .

