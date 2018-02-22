dir_run=$1
exec_new=$2
exec_old=$3

fname="INFO_SEA.txt"

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
echo "############# SEAFLUX_INTERPOL_SST_LINEAR "
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CINTERPOL\_SST\ =\ \"NONE\"/CINTERPOL\_SST\ =\ \"LINEAR\"/g" OPTIONS.nam_save > OPTIONS.nam
#CINTERPOL_SST ne fonctionne qu'en mode couplé => en offline on ne teste que le passage dans les routines, pas le résultat
#=> 1 test est suffisant
cp -f OPTIONS.nam OPTIONS.nam_SEAFLUX_INTERPOL_SST_LINEAR
#./script_exec_parall.sh "SEAFLUX_INTERPOL_SST_LINEAR" $fname $2 $3
./script_sea_canopy.sh "OPTIONS.nam" "SEAFLUX_INTERPOL_SST_LINEAR" $fname $2 $3

echo "############# SEAFLUX_INTERPOL_SST_QUADRA "
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CINTERPOL\_SST\ =\ \"NONE\"/CINTERPOL\_SST\ =\ \"QUADRA\"/g" OPTIONS.nam_save > OPTIONS.nam
#CINTERPOL_SST ne fonctionne qu'en mode couplé => en offline on ne teste que le passage dans les routines, pas le résultat
#=> 1 test est suffisant
#./script_exec_parall.sh "SEAFLUX_INTERPOL_SST_QUADRA" $fname $2 $3
./script_sea_canopy.sh "OPTIONS.nam" "SEAFLUX_INTERPOL_SST_QUADRA" $fname $2 $3

#pour appeler ice_sea_flux (sst + faible)
echo "############# SEAFLUX_SST_UNIF_271.15 "
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/XSST\_UNIF\ =\ 283\.15/XSST\_UNIF\ =\ 271\.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SEAFLUX_SST_UNIF_271.15
./script_exec.sh "SEAFLUX_SST_UNIF_271.15" $fname $2 $3


#avec des SST_DATA (évolution dans le temps du forçage SST)
echo "####################SEAFLUX_SST_DATA_ECUME"
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/LSST_DATA = F/LSST_DATA = T/g" OPTIONS.nam_save > OPTIONS.nam_sea
cp -f OPTIONS.nam_sea OPTIONS.nam
./script_sea_canopy.sh "OPTIONS.nam" "SEAFLUX_SST_DATA_ECUME" $fname $2 $3

#tests des différentes paramétrisations pour les flux
echo "###############SEAFLUX_SST_DATA_ECUME_PRECIP_PWEBB_ICHCE0.75"
cp -f OPTIONS.nam_sea OPTIONS.nam_save
sed -e "s/LPRECIP\ =\ F/LPRECIP\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LPWEBB\ =\ F/LPWEBB\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XICHCE\ =\ 0\./XICHCE\ =\ 0\.75/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_web
./script_sea_canopy.sh "OPTIONS.nam" "SEAFLUX_SST_DATA_ECUME_PRECIP_PWEBB_ICHCE0.75" $fname $2 $3

echo "#################SEAFLUX_SST_DATA_ECUME6_PRECIP_PWEBB_ICHCE0.75"
cp -f OPTIONS.nam_web OPTIONS.nam_save
sed -e "s/CSEA\_FLUX\ =\ \"ECUME\"/CSEA\_FLUX\ =\ \"ECUME6\"/g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_SEAFLUX_SST_DATA_ECUME6_PRECIP_PWEBB_ICHCE0.75
./script_exec.sh "SEAFLUX_SST_DATA_ECUME6_PRECIP_PWEBB_ICHCE0.75" $fname $2 $3

echo "###################SEAFLUX_SST_DATA_ECUME6_INTERPOL_SSS_UNIF"
cp -f OPTIONS.nam_sea OPTIONS.nam_save
sed -e "s/CSEA\_FLUX\ =\ \"ECUME\"/CSEA\_FLUX\ =\ \"ECUME6\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ecume6
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CINTERPOL_SSS =\ \"NONE\"/CINTERPOL_SSS\ =\ \"UNIF\"/g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_SEAFLUX_SST_DATA_ECUME6_INTERPOL_SSS_UNIF
./script_exec.sh "SEAFLUX_SST_DATA_ECUME6_INTERPOL_SSS_UNIF" $fname $2 $3

echo "#########SEAFLUX_SST_DATA_ECUME6_Z01"
cp -f OPTIONS.nam_ecume6 OPTIONS.nam_save
sed -e "s/NZ0\ =\ 0/NZ0\ =\ 1/g" OPTIONS.nam_save > OPTIONS.nam
./script_sea_canopy.sh "OPTIONS.nam" "SEAFLUX_SST_DATA_ECUME6_Z01" $fname $2 $3

echo "###############SEAFLUX_SST_DATA_ECUME6_Z02"
cp -f OPTIONS.nam_ecume6 OPTIONS.nam_save
sed -e "s/NZ0\ =\ 0/NZ0\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
./script_sea_canopy.sh "OPTIONS.nam" "SEAFLUX_SST_DATA_ECUME6_Z02" $fname $2 $3

echo "###############SEAFLUX_SST_DATA_DIRECT"
cp -f OPTIONS.nam_sea OPTIONS.nam_save
sed -e "s/CSEA\_FLUX\ =\ \"ECUME\"/CSEA\_FLUX\ =\ \"DIRECT\"/g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_SEAFLUX_SST_DATA_DIRECT
./script_exec.sh "SEAFLUX_SST_DATA_DIRECT" $fname $2 $3

echo "##########SEAFLUX_SST_DATA_ITERAT"
cp -f OPTIONS.nam_sea OPTIONS.nam_save
sed -e "s/CSEA\_FLUX\ =\ \"ECUME\"/CSEA\_FLUX\ =\ \"ITERAT\"/g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_SEAFLUX_SST_DATA_ITERAT
./script_exec.sh "SEAFLUX_SST_DATA_ITERAT" $fname $2 $3

echo "############SEAFLUX_SST_DATA_COARE30"
cp -f OPTIONS.nam_sea OPTIONS.nam_save
sed -e "s/CSEA\_FLUX\ =\ \"ECUME\"/CSEA\_FLUX\ =\ \"COARE3\"/g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_SEAFLUX_SST_DATA_COARE30
./script_exec_parall.sh "SEAFLUX_SST_DATA_COARE30" $fname $2 $3

echo "##############SEAFLUX_SST_DATA_COARE30_PWG"
cp -f OPTIONS.nam_sea OPTIONS.nam_save
sed -e "s/CSEA\_FLUX\ =\ \"ECUME\"/CSEA\_FLUX\ =\ \"COARE3\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LPWG\ =\ F/LPWG\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
./script_sea_canopy.sh "OPTIONS.nam" "SEAFLUX_SST_DATA_COARE30_PWG" $fname $2 $3

echo "################SEAFLUX_SST_DATA_COARE30_NGRVWAVES1"
cp -f OPTIONS.nam_sea OPTIONS.nam_save
sed -e "s/CSEA\_FLUX\ =\ \"ECUME\"/CSEA\_FLUX\ =\ \"COARE3\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGRVWAVES\ =\ 0/NGRVWAVES\ =\ 1/g" OPTIONS.nam_save > OPTIONS.nam
./script_sea_canopy.sh "OPTIONS.nam" "SEAFLUX_SST_DATA_COARE30_NGRVWAVES1" $fname $2 $3

echo "#############SEAFLUX_SST_DATA_COARE30_NGRVWAVES2"
cp -f OPTIONS.nam_sea OPTIONS.nam_save
sed -e "s/CSEA\_FLUX\ =\ \"ECUME\"/CSEA\_FLUX\ =\ \"COARE3\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGRVWAVES\ =\ 0/NGRVWAVES\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
./script_sea_canopy.sh "OPTIONS.nam" "SEAFLUX_SST_DATA_COARE30_NGRVWAVES2" $fname $2 $3

echo "################SEAFLUX_SST_DATA_SEA_ALBTA96"
#tests albedo
cp -f OPTIONS.nam_sea OPTIONS.nam_save
sed -e "s/CSEA\_ALB\ =\ \"UNIF\"/CSEA\_ALB\ =\ \"TA96\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_sea_canopy.sh "OPTIONS.nam" "SEAFLUX_SST_DATA_SEA_ALBTA96" $fname $2 $3

echo "################SEAFLUX_SST_DATA_SEA_ALBMK10"
cp -f OPTIONS.nam_sea OPTIONS.nam_save
sed -e "s/CSEA\_ALB\ =\ \"UNIF\"/CSEA\_ALB\ =\ \"MK10\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_sea_canopy.sh "OPTIONS.nam" "SEAFLUX_SST_DATA_SEA_ALBMK10" $fname $2 $3


#test gelato
echo "################ SEAFLUX_GELATO "
cp -f TESTS/OPTIONS.nam_000 OPTIONS.nam_save
sed -e "s/XUNIF_NATURE = 1\.0E+20/XUNIF_NATURE = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_SEA = 1\.0E+20/XUNIF_SEA = 1./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_WATER = 1\.0E+20/XUNIF_WATER = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TOWN = 1\.0E+20/XUNIF_TOWN = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"LONLATVAL\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NPOINTS = 127/NPOINTS = 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XX = \-88\.23/XX = \-5\.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XY = 47\.72/XY = \-65\.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDX = 0\.5/XDX = 1\.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDY = 0\.5/XDY = 1\.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CNATURE = \"ISBA\"/CNATURE = \"NONE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTOWN = \"TEB\"/CTOWN = \"NONE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CWATER = \"WATFLX\"/CWATER = \"NONE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_COVER(1) = 0./XUNIF_COVER(1) = 1./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ZS = 1.0E+20/XUNIF_ZS = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NYEAR = 1000000000/NYEAR = 1985/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH = 1000000000/NMONTH = 07/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY = 1000000000/NDAY = 01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTIME = 1\.0E+20/XTIME = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1\.0E+20/XSST_UNIF = 271\.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSSS_UNIF = 1\.0E+20/XSSS_UNIF = 33\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSIC_UNIF = 1\.0E+20/XSIC_UNIF = 1\.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSEAICE_SCHEME = \"NONE\"/CSEAICE_SCHEME = \"GELATO\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XVZ0CM = 0./XVZ0CM = 0.00001/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LNOSOF = F/LNOSOF = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LVERTSHIFT = T/LVERTSHIFT = F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XZ0SN = 0.001/XZ0SN = 0.03/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XZ0HSN = 0.0001/XZ0HSN = 0.003/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XWSNV = 5.0/XWSNV = 2.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSEA_ALB = \"UNIF\"/CSEA_ALB = \"TA96\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CINTERPOL_SST = \"NONE\"/CINTERPOL_SST = \"QUADRA\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XICHCE = 0./XICHCE = 1./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_SEABATHY = -100./XUNIF_SEABATHY = -300./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFORCING_FILETYPE = \"ASCII\"/CFORCING_FILETYPE = \"NETCDF\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTSTEP_SURF = 300./XTSTEP_SURF = 10800./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LLIMIT_QAIR = F/LLIMIT_QAIR = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDELTA_OROG = 200./XDELTA_OROG = 2200./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LDIAG_SEAICE = F/LDIAG_SEAICE = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSIC_EFOLDING_TIME = 0./XSIC_EFOLDING_TIME = 10./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSEAICE_TSTEP = 1.0E+20/XSEAICE_TSTEP = 10800./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XCD_ICE_CST = 0./XCD_ICE_CST = 0.0015/g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_SEAFLUX_GELATO
cp -f TESTS/FORCAGES/GELATO/FORCING.nc .
./script_exec.sh "SEAFLUX_GELATO" $fname $2 $3

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

