dir_run=$1
exec_new=$2
exec_old=$3

fname="INFO_PGD_PHYSIO.txt"

rm -f $fname
touch $fname

cd $1

cp -f TESTS/script_pgd_exte.sh .
cp -f TESTS/script_exec_pgd_parall.sh .
cp -f TESTS/script_to_old.sh .

rm -f Forc_*.txt
rm -f Params_config.txt

cp -f TESTS/PREP/make*.sh .
cp -f TESTS/PGD/make*.sh .

cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
ln -s TESTS/FORCAGES/HIVER/Forc*.txt .

####################TESTS PHYSIO SIMPLES###########################################

#namelists avec PREP uniformes
cp TESTS/OPTIONS.nam_000 OPTIONS.nam_NO_ECOCLIMAP
./make_prep_unif.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_prep_isba_unif.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_prep_isba_snow_unif.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_prep_teb_unif.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_prep_gd_gr_unif.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_prep_teb_snow_unif.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_prep_flake_unif.sh "OPTIONS.nam_NO_ECOCLIMAP"
cp -f OPTIONS.nam_NO_ECOCLIMAP OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1.0E+20/XSST_UNIF = 283.15/g" OPTIONS.nam_save > OPTIONS.nam_NO_ECOCLIMAP
cp -f OPTIONS.nam_NO_ECOCLIMAP OPTIONS.nam_save
sed -e "s/XTS_WATER_UNIF = 1.0E+20/XTS_WATER_UNIF = 282./g" OPTIONS.nam_save > OPTIONS.nam_NO_ECOCLIMAP

#ca sans ecoclimap
./make_eco_tiles.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_isba_unif.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_teb_unif1.sh "OPTIONS.nam_NO_ECOCLIMAP"

#---------------------------------------------------------------------------------

#tests NAM_ZS


cp -f OPTIONS.nam_NO_ECOCLIMAP OPTIONS.nam_save
sed -e "s/LSET_FORC_ZS = F/LSET_FORC_ZS = T/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF

# 1.1. Test GTOPO30 OROGTYPE = "MAX"
echo "########## GTOPO30_OROG_MAX"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save 
sed -e "s/COROGTYPE = \"ENV\"/COROGTYPE = \"MAX\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_GTOPO30_OROG_MAX
./script_exec_pgd_parall.sh "GTOPO30_OROG_MAX" $fname $2 $3

# 1.2. OROGTYPE = "AVG"
echo "########## GTOPO30_OROG_AVG"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save 
sed -e "s/COROGTYPE = \"ENV\"/COROGTYPE = \"AVG\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_GTOPO30_OROG_AVG
./script_exec_pgd_parall.sh "GTOPO30_OROG_AVG" $fname $2 $3

# 1.3. OROGTYPE = "SIL"
echo "########## GTOPO30_OROG_SIL"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save 
sed -e "s/COROGTYPE = \"ENV\"/COROGTYPE = \"SIL\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_GTOPO30_OROG_SIL
./script_exec_pgd_parall.sh "GTOPO30_OROG_SIL" $fname $2 $3

# 1.4. OROGTYPE = "ENV" + XENV = 0.2
echo "########## GTOPO30_OROG_ENV_ENV0.2"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save 
sed -e "s/XENV = 0./XENV = 0.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_GTOPO30_OROG_ENV0.2
./script_exec_pgd_parall.sh "GTOPO30_OROG_ENV_ENV0.2" $fname $2 $3

# 1.5. ZSFILTER = 0 => en fait zsfilter 
echo "########## GTOPO30_OROGMAX_ZSFILTER0"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save 
sed -e "s/NZSFILTER = 1/NZSFILTER = 0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_GTOPO30_OROGMAX_ZSFILTER0
./script_exec_pgd_parall.sh "GTOPO30_OROGMAX_ZSFILTER0" $fname $2 $3

# 1.6. ZSFILTER = 3
echo "########## GTOPO30_OROGMAX_ZSFILTER3"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save 
sed -e "s/NZSFILTER = 1/NZSFILTER = 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_GTOPO30_OROGMAX_ZSFILTER3
./script_exec_pgd_parall.sh "GTOPO30_OROGMAX_ZSFILTER3" $fname $2 $3

# 1.7. LIMP_ZS = T + PGD_BASE
echo "########## LIMP_ZS"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/CSURF_FILETYPE = \"NC\"/CSURF_FILETYPE = \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_pgd_exte.sh "LIMP_ZS" $fname $exec_new $exec_old "lfi"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YZS = \"gtopo30\"/YZS = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YZSFILETYPE = \"DIRECT\"/YZSFILETYPE = \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LIMP\_ZS\ =\ F/LIMP\_ZS\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LIMP_ZS
./script_exec_pgd_parall.sh "LIMP_ZS" $fname $2 $3

# 1.8. YSLOPE
echo "########## YSLOPE"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/NIMAX = 12/NIMAX = 17/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LIMP\_ZS\ =\ F/LIMP\_ZS\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YZS = \"gtopo30\"/YZS = \"FORCING.nc\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YZSFILETYPE = \"DIRECT\"/YZSFILETYPE = \"NETCDF\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YSLOPE = \"\"/YSLOPE = \"FORCING.nc\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YSLOPEFILETYPE = \"\"/YSLOPEFILETYPE = \"NETCDF\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NYEAR = 1986/NYEAR = 2013/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH = 1/NMONTH = 8/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY = 29/NDAY = 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTIME = 0./XTIME = 21600./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFORCING_FILETYPE = \"ASCII\"/CFORCING_FILETYPE = \"NETCDF\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ZS_SLOPE_NETCDF
cp -f TESTS/FORCAGES/SLOPE/FORCING.nc .
./script_exec_pgd_parall.sh "ZS_SLOPE_NETCDF" $fname $2 $3
rm -f FORCING.nc

# 1.9. LEXPLICIT_SLOPE
echo "########## EXPLICIT_SLOPE"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/LEXPLICIT_SLOPE = F/LEXPLICIT_SLOPE = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_EXPLICIT_SLOPE
./script_exec_pgd_parall.sh "EXPLICIT_SLOPE" $fname $2 $3


# 2. Test NAM_SEABATHY

# 2.1. Test YSEABATHY = "etopo2"
echo "########## SEABATHY_ETOPO2"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YSEABATHY = \"\"/YSEABATHY = \"etopo2.nc\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SEABATHY_ETOPO2
./script_exec_pgd_parall.sh "SEABATHY_ETOPO2" $fname $2 $3


# 3. Test NAM_ISBA

echo "########## ISBA_PHYSIO_MAPS"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YCLAY = \"\"/YCLAY = \"clay_fao\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YSAND = \"\"/YSAND = \"sand_fao\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YSOC_TOP = \"\"/YSOC_TOP = \"soc_top\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YSOC_SUB = \"\"/YSOC_SUB = \"soc_sub\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YCTI = \"\"/YCTI = \"topo_index\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YPERM = \"\"/YPERM = \"perm_glo_10km\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ISBA_PHYSIO_MAPS
./script_exec_pgd_parall.sh "ISBA_PHYSIO_MAPS" $fname $2 $3


#4. Test NAM_FLAKE

sed -e "s/CWATER\ =\ \"WATFLX\"/CWATER\ =\ \"FLAKE\"/g" OPTIONS.nam_PHYSIO_UNIF > OPTIONS.nam_FLAKE

# 4.1. cas uniforme
echo "########## FLAKE_UNIF_DATA"
cp -f OPTIONS.nam_FLAKE OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FLAKE_UNIF_DATA
./script_exec_pgd_parall.sh "FLAKE_UNIF_DATA" $fname $2 $3

# 4.2. flake database
echo "########## FLAKE_DEPTH_DATA"
cp -f OPTIONS.nam_FLAKE OPTIONS.nam_save
sed -e "s/YWATER_DEPTH = \"\"/YWATER_DEPTH = \"GlobalLakeDepth\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FLAKE_DEPTH_DATA
./script_exec_pgd_parall.sh "FLAKE_DEPTH_DATA" $fname $2 $3


# 5. Test NAM_DATA_SEAFLUX

# 5.1. avec SST en entrée
echo "########### SST_DATA"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/LSST_DATA = F/LSST_DATA = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SST_DATA
./script_exec_pgd_parall.sh "SST_DATA" $fname $2 $3

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

rm -f Params_config.txt
rm -f Forc*.txt

