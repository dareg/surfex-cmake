dir_run=$1
exec_new=$2
exec_old=$3

fname="INFO_WATER.txt"

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

cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
ln -s TESTS/FORCAGES/HIVER/Forc*.txt .

####################NAMELISTS INITIALES##########################################

#namelists avec PREP uniformes
cp TESTS/OPTIONS.nam_000 OPTIONS.nam_PHYSIO_UNIF
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

#ZS uniforme
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/XUNIF_ZS = 1.0E+20/XUNIF_ZS = 250./g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF

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

### il faudrait faire des tests avec une température = 0, pour la glace

####################TESTS POUR L'INTERPOLATION SST##########################################

#on vérifie juste que ça tourne, car la SST interpolée est toujours constante
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CINTERPOL\_TS\ =\ \"NONE\"/CINTERPOL\_TS\ =\ \"LINEAR\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_water_canopy.sh "OPTIONS.nam" "WATFLUX_INTERPOL_TS_LINEAR" $fname $2 $3

cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CINTERPOL\_TS\ =\ \"NONE\"/CINTERPOL\_TS\ =\ \"QUADRA\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_water_canopy.sh "OPTIONS.nam" "WATFLUX_INTERPOL_TS_QUADRA" $fname $2 $3

cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CINTERPOL\_TS\ =\ \"NONE\"/CINTERPOL\_TS\ =\ \"UNIF\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_water_canopy.sh "OPTIONS.nam" "WATFLUX_INTERPOL_TS_UNIF" $fname $2 $3

#pour appeler ice_sea_flux (sst + faible)
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/XTS\_WATER\_UNIF\ =\ 282\./XTS\_WATER\_UNIF\ =\ 271\./g" OPTIONS.nam_save > OPTIONS.nam
./script_water_canopy.sh "OPTIONS.nam" "WATFLUX_TS_WATER_UNIF271" $fname $2 $3


#tests albedo
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CWAT\_ALB\ =\ \"UNIF\"/CWAT\_ALB\ =\ \"TA96\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_water_canopy.sh "OPTIONS.nam" "WATFLUX_CWAT_ALBTA96" $fname $2 $3

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

