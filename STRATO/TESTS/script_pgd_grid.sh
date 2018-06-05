#set -x 

dir_run=$1
exec_new=$2
exec_old=$3

#fichier pour voir si les simuls ont tourné jusqu'à la fin 
fname="INFO_PGD_GRID.txt"

#on crée un nouveau fichier fname
rm -f $fname
touch $fname

#répertoire de run (à adapter le cas échéant), ou passer en arguemnt
cd $1

#les scripts de tests
cp -f TESTS/script_pgd_exte.sh .
cp -f TESTS/script_exec_pgd_parall.sh .
cp -f TESTS/script_exec_all_parall.sh .
#cp -f TESTS/script_exec_omp_pgd.sh .
cp -f TESTS/script_to_old.sh .

#les forçages
rm -f Forc_*.txt
rm -f Params_config.txt

cp -f TESTS/PREP/make*.sh .
cp -f TESTS/PGD/make*.sh .

#on prend les forçages d'hiver pour avoir de la neige (forçages constants)
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
ln -s TESTS/FORCAGES/HIVER/Forc*.txt .

####################TESTS DE GRILLES###########################################
#il faut qu'interpol_npts soit appelé au moins une fois par grille

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

#------------------------------------------------------------------------------

### 1.1. test CARTESIAN
echo "########## CARTESIAN "
cp -f OPTIONS.nam_NO_ECOCLIMAP OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"CARTESIAN\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_CARTESIAN
./script_exec_all_parall.sh "CARTESIAN" $fname $2 $3

# 1.2. test inifile CARTESIAN
echo "########## CARTESIAN_INIFILE_LFI"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/CSURF_FILETYPE = \"NC\"/CSURF_FILETYPE = \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_pgd_exte.sh "CARTESIAN_INIFILE_LFI" $fname $exec_new $exec_old "lfi" 
if [ ! -f PGD_BASE_OLD.lfi ] 
then
  ./script_pgd_exte.sh "CARTESIAN_INIFILE_LFI_OLD" $fname $exec_new $exec_old "lfi" 
fi
cp -f OPTIONS.nam_CARTESIAN OPTIONS.nam_save
sed -e "s/YINIFILE = \"\"/YINIFILE = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YINIFILETYPE = \"\"/YINIFILE = \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_CARTESIAN_INIFILE_LFI
./script_exec_pgd_parall.sh "CARTESIAN_INIFILE_LFI" $fname $2 $3

# 1.3. test inifile CARTESIAN pgd de la précédente version (le calculer à part)
echo "########### CARTESIAN_INIFILE_LFI_OLD"
rm -f PGD_BASE_NEW.lfi
mv -f PGD_BASE_OLD.lfi PGD_BASE.lfi
./script_exec_pgd_parall.sh "CARTESIAN_INIFILE_LFI_OLD" $fname $2 $3
rm -f PGD_BASE.lfi


# 2.1. test CONF PROJ simple
echo "########## CONF_PROJ"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_CONF_PROJ
./script_exec_all_parall.sh "CONF_PROJ" $fname $2 $3

# test CONF PROJ high res
echo "########## CONF_PROJ_INTERPOL"
cp -f OPTIONS.nam_CONF_PROJ OPTIONS.nam_save
sed -e "s/XDX = 60000./XDX = 250./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_CONF_PROJ_INTERPOL
./script_exec_pgd_parall.sh "CONF_PROJ_INTERPOL" $fname $2 $3


# 2.2. test inifile CONF PROJ
echo "########## CONF_PROJ_INIFILE_ASCII"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/CSURF_FILETYPE = \"NC\"/CSURF_FILETYPE = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_pgd_exte.sh "CONF_PROJ_INIFILE_ASCII" $fname $exec_new $exec_old "txt"
if [ ! -f PGD_BASE_OLD.txt ] 
then
  ./script_pgd_exte.sh "CONF_PROJ_INIFILE_ASCII_OLD" $fname $exec_new $exec_old "txt"
fi

cp -f OPTIONS.nam_CONF_PROJ OPTIONS.nam_save
sed -e "s/YINIFILE = \"\"/YINIFILE = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YINIFILETYPE = \"\"/YINIFILETYPE = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_CONF_PROJ_INIFILE_ASCII
./script_exec_pgd_parall.sh "CONF_PROJ_INIFILE_ASCII" $fname $2 $3

# 2.3. test inifile CONF PROJ pgd de la précédente version (le calculer à part)
echo "########## CONF_PROJ_INIFILE_ASCII_OLD"
rm -f PGD_BASE_NEW.txt
mv -f PGD_BASE_OLD.txt PGD_BASE.txt
./script_exec_pgd_parall.sh "CONF_PROJ_INIFILE_ASCII_OLD" $fname $2 $3
rm -f PGD_BASE.txt

# 3. test LONLAT REG
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"LONLAT REG\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LONLAT_REG
echo "########## LONLAT_REG"
./script_exec_all_parall.sh "LONLAT_REG" $fname $2 $3

#LONLAT REG high res
echo "########## LONLAT_REG_INTERPOL"
cp -f OPTIONS.nam_LONLAT_REG OPTIONS.nam_save
sed -e "s/XLONMIN = -10./XLONMIN = 3./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLONMAX = 80./XLONMAX = 3.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATMIN = 10./XLATMIN = 43./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATMAX = 55./XLATMAX = 43.125/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO = 2/NHALO = 4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LONLAT_REG_INTERPOL
./script_exec_pgd_parall.sh "LONLAT_REG_INTERPOL" $fname $2 $3

# 4.1. test GAUSS simple
echo "########## GAUSS"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"GAUSS\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_COVER(4) = 0./XUNIF_COVER(4) = 1./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YCLAY = \"\"/YCLAY = \"clay_fao\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO = 2/NHALO = 8/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_GAUSS_SIMPLE
./script_exec_pgd_parall.sh "GAUSS" $fname $2 $3

# 4.2. test GAUSS + LORCA_GRID 
echo "########## GAUSS_ORCA_GRID"
cp -f OPTIONS.nam_GAUSS_SIMPLE OPTIONS.nam_save 
sed -e "s/LORCA_GRID = F/LORCA_GRID = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLAT_ANT = -77./XLAT_ANT = -75./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_GAUSS_SIMPLE_ORCA
./script_exec_all_parall.sh "GAUSS_ORCA_GRID" $fname $2 $3

# 4.3. test GAUSS renversée
echo "########## GAUSS_BASCULE"
cp -f OPTIONS.nam_GAUSS_SIMPLE OPTIONS.nam_save
sed -e "s/RMUCEN = 1./RMUCEN = 0.725/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/RLOCEN = 0./RLOCEN = 0.045/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/RSTRET = 1./RSTRET = 2.4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_GAUSS_BASCULE
./script_exec_pgd_parall.sh "GAUSS_BASCULE" $fname $2 $3

# 5.1. test IGN tous points
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"IGN\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_IGN_L2E_IRREG
echo "########## IGN_L2E_IRREG"
./script_exec_all_parall.sh "IGN_L2E_IRREG" $fname $2 $3


# 5.2. test IGN coins
echo "########## IGN_L2E_REG"
cp -f OPTIONS.nam_IGN_L2E_IRREG OPTIONS.nam_save
sed -e "s/XX_LLCORNER = 1.0E+20/XX_LLCORNER = 870000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XY_LLCORNER = 1.0E+20/XY_LLCORNER = 600000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XCELLSIZE = 1.0E+20/XCELLSIZE = 8000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NCOLS = 0/NCOLS = 10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NROWS = 0/NROWS = 10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_IGN_L2E_REG
./script_exec_all_parall.sh "IGN_L2E_REG" $fname $2 $3

#IGN high res
echo "########## IGN_L2E_REG_INTERPOL"
cp -f OPTIONS.nam_IGN_L2E_REG OPTIONS.nam_save
sed -e "s/XCELLSIZE = 8000./XCELLSIZE = 250./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO = 2/NHALO = 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_IGN_L2E_REG_INTERPOL
./script_exec_pgd_parall.sh "IGN_L2E_REG_INTERPOL" $fname $2 $3


# 5.3. test IGN autres types de LAMBERT
echo "########## IGN_L1_REG"
sed -e "s/L2E/L1/g" OPTIONS.nam_IGN_L2E_REG > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_IGN_L1_REG
./script_exec_all_parall.sh "IGN_L1_REG" $fname $2 $3
echo "########## IGN_L2_REG"
sed -e "s/L2E/L2/g" OPTIONS.nam_IGN_L2E_REG > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_IGN_L2_REG
./script_exec_pgd_parall.sh "IGN_L2_REG" $fname $2 $3
echo "########## IGN_L3_REG"
sed -e "s/L2E/L3/g" OPTIONS.nam_IGN_L2E_REG > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_IGN_L3_REG
./script_exec_pgd_parall.sh "IGN_L3_REG" $fname $2 $3
echo "########## IGN_L4_REG"
sed -e "s/L2E/L4/g" OPTIONS.nam_IGN_L2E_REG > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_IGN_L4_REG
./script_exec_pgd_parall.sh "IGN_L4_REG" $fname $2 $3
echo "########## IGN_L93_REG"
sed -e "s/L2E/L93/g" OPTIONS.nam_IGN_L2E_REG > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_IGN_L93_REG
./script_exec_pgd_parall.sh "IGN_L93_REG" $fname $2 $3

# 6.1. test LONLATVAL 
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"LONLATVAL\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LONLATVAL_NOVMX1
echo "########## LONLATVAL_NOVMX1"
./script_exec_all_parall.sh "LONLATVAL_NOVMX1" $fname $2 $3

# 6.2. test LONLATVAL avec NOVMX = 2
cp -f OPTIONS.nam_LONLATVAL_NOVMX1 OPTIONS.nam_save
sed -e "s/NOVMX = 1/NOVMX = 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LONLATVAL_NOVMX2
echo "########## LONLATVAL_NOVMX2"
./script_exec_all_parall.sh "LONLATVAL_NOVMX2" $fname $2 $3


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
