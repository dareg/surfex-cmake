dir_run=$1
exec_new=$2
exec_old=$3

fname="INFO_PGD_ECOCLIMAP.txt"

rm -f $fname
touch $fname

cd $1

cp -f TESTS/script_exec_pgd_parall.sh .
cp -f TESTS/script_to_old.sh .

rm -f Forc_*.txt
rm -f Params_config.txt

cp -f TESTS/PREP/make*.sh .
cp -f TESTS/PGD/make*.sh .

cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
ln -s TESTS/FORCAGES/HIVER/Forc*.txt .

####################TESTS ECOCLIMAP###########################################

### Tests sans ECOCLIMAP, différents types d'entrées

#namelist avec OPTIONS TEB activées
sed -e "s/LGARDEN = F/LGARDEN = T/g" TESTS/OPTIONS.nam_000 > OPTIONS.nam_OPT
mv -f OPTIONS.nam_OPT OPTIONS.nam_save
sed -e "s/CBEM = \"DEF\"/CBEM = \"BEM\"/g" OPTIONS.nam_save > OPTIONS.nam_OPT
mv -f OPTIONS.nam_OPT OPTIONS.nam_save
sed -e "s/LGREENROOF = F/LGREENROOF = T/g" OPTIONS.nam_save > OPTIONS.nam_OPT
mv -f OPTIONS.nam_OPT OPTIONS.nam_save
sed -e "s/XUNIF_GREENROOF = 1.0E+20/XUNIF_GREENROOF = 0.5/g" OPTIONS.nam_save > OPTIONS.nam_OPT

#namelists avec PREP uniformes
./make_prep_unif.sh "OPTIONS.nam_OPT"
./make_prep_isba_unif.sh "OPTIONS.nam_OPT"
./make_prep_isba_snow_unif.sh "OPTIONS.nam_OPT"
./make_prep_teb_unif.sh "OPTIONS.nam_OPT"
./make_prep_gd_gr_unif.sh "OPTIONS.nam_OPT"
./make_prep_teb_snow_unif.sh "OPTIONS.nam_OPT"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam_OPT"
./make_prep_flake_unif.sh "OPTIONS.nam_OPT"
cp -f OPTIONS.nam_OPT OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1.0E+20/XSST_UNIF = 283.15/g" OPTIONS.nam_save > OPTIONS.nam_OPT
cp -f OPTIONS.nam_OPT OPTIONS.nam_save
sed -e "s/XTS_WATER_UNIF = 1.0E+20/XTS_WATER_UNIF = 282./g" OPTIONS.nam_save > OPTIONS.nam_OPT

cp -f OPTIONS.nam_OPT OPTIONS.nam_save
sed -e "s/XUNIF_ZS = 1.0E+20/XUNIF_ZS = 250./g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF

#recup d'un cas ECOCLIMAP 
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam
sed -e "s/LECOCLIMAP = F/LECOCLIMAP = T/g" OPTIONS.nam > OPTIONS.nam_ECOCLIMAP
./make_greenroof_unif1.sh "OPTIONS.nam_ECOCLIMAP"

###recup d'un cas sans ECOCLIMAP 
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_NO_ECOCLIMAP
./make_bem_unif1.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_garden_unif1.sh "OPTIONS.nam_NO_ECOCLIMAP"
./make_greenroof_unif1.sh "OPTIONS.nam_NO_ECOCLIMAP"
cp -f OPTIONS.nam_NO_ECOCLIMAP OPTIONS.nam_save
sed -e "s/XUNIF_GARDEN = 1.0E+20/XUNIF_GARDEN = 0.4/g" OPTIONS.nam_save > OPTIONS.nam_NO_ECOCLIMAP

#---------------------------------------------------------------------------------

# 1.1. test sans ECOCLIMAP avec toutes les options de TEB: grille conf_proj1
echo "########## NO_ECO_ALL_TEB_OPT"
cp -f OPTIONS.nam_NO_ECOCLIMAP OPTIONS.nam
./make_eco_tiles.sh "OPTIONS.nam"
./make_isba_unif.sh "OPTIONS.nam"
./make_teb_unif2.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_NO_ECO_ALL_TEB_OPT
./script_exec_pgd_parall.sh "NO_ECO_ALL_TEB_OPT" $fname $2 $3

# 1.2. La même, mais on enlève quelques données TEB
echo "########## NO_ECO_SOME_TEB_OPT"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_isba_unif.sh "OPTIONS.nam"
./make_teb_unif3.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_NO_ECO_SOME_TEB_OPT
./script_exec_pgd_parall.sh "NO_ECO_SOME_TEB_OPT" $fname $2 $3

# 1.3. sans ECOCLIMAP avec quelques fichiers (cas GREG)
#Attention: on ne peut pas utiliser NAM_DATA_ISBA et vouloir utiliser ECOCLIMAP pour GARDEN
ln -s TESTS/PGD/DATA/data_greg/* .
echo "########## NO_ECO_TEB_OPT_FILES"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_isba_greg.sh "OPTIONS.nam"
./make_teb_greg.sh "OPTIONS.nam"
./make_garden_unif2.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"LONLAT REG\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLONMIN = -10./XLONMIN = 2.17/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLONMAX = 80./XLONMAX = 2.56/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATMIN = 10./XLATMIN = 48.76/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATMAX = 55./XLATMAX = 48.95/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NLON = 90/NLON = 10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NLAT = 45/NLAT = 10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO = 2/NHALO = 4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NO_ECO_TEB_OPT_FILES
./script_exec_pgd_parall.sh "NO_ECO_TEB_OPT_FILES" $fname $2 $3
#
for file in TESTS/PGD/DATA/data_greg/*
do
	file2=$(basename $file)
	rm -f $file2
done
#
# 1.4. sans ECOCLIMAP avec l'initialisation particulière pour TEB (cas CECILE)
ln -s TESTS/PGD/DATA/data_muscade/files .
echo "########## NO_ECO_TEB_OPT_CSV"
cp -f OPTIONS.nam_NO_ECOCLIMAP OPTIONS.nam
./make_eco_tiles.sh "OPTIONS.nam"
./make_isba_unif.sh "OPTIONS.nam"
./make_teb_cecile.sh "OPTIONS.nam"
./make_bem_unif2.sh "OPTIONS.nam"
./make_garden_unif2.sh "OPTIONS.nam"
./make_greenroof_unif2.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CGRID = \"CONF\ PROJ\"/CGRID = \"IGN\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CLAMBERT = \"L2E\"/CLAMBERT = \"L93\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XX_LLCORNER = 1.0E+20/XX_LLCORNER = 644000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XY_LLCORNER = 1.0E+20/XY_LLCORNER = 6862499./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XCELLSIZE = 1.0E+20/XCELLSIZE = 10000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NCOLS = 0/NCOLS = 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NROWS = 0/NROWS = 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NO_ECO_TEB_OPT_CSV
./script_exec_pgd_parall.sh "NO_ECO_TEB_OPT_CSV" $fname $2 $3
#
rm -f files
#
# 1.5 sans ecoclimap avec des fichiers pour ISBA
ln -s TESTS/PGD/DATA/data_mireille/* .
echo "########## NO_ECO_ISBA_FILES"
cp -f OPTIONS.nam_NO_ECOCLIMAP OPTIONS.nam
./make_eco_tiles.sh "OPTIONS.nam"
./make_isba_mireille.sh "OPTIONS.nam"
./make_teb_unif1.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLAT0 = 45./XLAT0 = 0./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLON0 = 3./XLON0 = 0./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRPK = 0./XRPK = 1.0/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XBETA = 0./XBETA = 90./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATCEN = 45./XLATCEN = 7.0/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLONCEN = 3./XLONCEN = 0./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NIMAX = 12/NIMAX = 120/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NJMAX = 8/NJMAX = 2/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDX = 60000./XDX = 6.95425E4/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDY = 60000./XDY = 1500.E3/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO = 2/NHALO = 60/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NO_ECO_ISBA_FILES
./script_exec_pgd_parall.sh "NO_ECO_ISBA_FILES" $fname $2 $3
#
for file in TESTS/PGD/DATA/data_mireille/*
do
	file2=$(basename $file)
	rm -f $file2
done
#
# 2.1. rien qu'avec ECOCLIMAP
echo "########## ALL_ECO"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ALL_ECO
./script_exec_pgd_parall.sh "ALL_ECO" $fname $2 $3


# 2.2. OPTIONS d'ECOCLIMAP 
echo "########## OPTIONS_ECO"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_save
sed -e "s/XRM_COVER = 1.E-6/XRM_COVER = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRM_COAST = 1./XRM_COAST = 0.50/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRM_LAKE = 0./XRM_LAKE = 0.02/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRM_SEA = 0./XRM_SEA = 0.02/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_OPTIONS_ECO
./script_exec_pgd_parall.sh "OPTIONS_ECO" $fname $2 $3

# 2.3. test avec XUNIF_COVER
echo "########## UNIF_COVER "
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_eco_unif.sh OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_UNIF_COVER
./script_exec_pgd_parall.sh "UNIF_COVER" $fname $2 $3

# 2.4. mixte ECOCLIMAP / fractions de TILES
echo "########## FRAC_TILES "
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO = 2/NHALO = 5/g" OPTIONS.nam_save > OPTIONS.nam
./make_eco_tiles.sh OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FRAC_TILES
./script_exec_pgd_parall.sh "FRAC_TILES" $fname $2 $3


# 3. ISBA: mixte ECOCLIMAP / NAM_DATA: temps (cas ADRIEN)
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
#si on ne met pas garden data, il essaie de lire les data d'isba et ça plante
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_all.sh "OPTIONS.nam"

# 3.1. VEGTYPES, LAI, VEG, Z0, DG, RSMIN, GAMMA, CV, ROOTFRAC - NTIME = 36,12,2,1
echo "########## ECO_DATA_ISBA_NTIME36"
cp -f OPTIONS.nam OPTIONS.nam_ECO_DATA_ISBA_NTIME36
./script_exec_pgd_parall.sh "ECO_DATA_ISBA_NTIME36" $fname $2 $3
echo "########## ECO_DATA_ISBA_NTIME12"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 36/NTIME = 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_DATA_ISBA_NTIME12
./script_exec_pgd_parall.sh "ECO_DATA_ISBA_NTIME12" $fname $2 $3
echo "########## ECO_DATA_ISBA_NTIME2"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 12/NTIME = 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_DATA_ISBA_NTIME2
./script_exec_pgd_parall.sh "ECO_DATA_ISBA_NTIME2" $fname $2 $3
echo "########## ECO_DATA_ISBA_NTIME1"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 2/NTIME = 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_DATA_ISBA_NTIME1
./script_exec_pgd_parall.sh "ECO_DATA_ISBA_NTIME1" $fname $2 $3

# 3.2. juste les vegtypes
echo "########## ECO_VGT_ISBA_NTIME36"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_vgt.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_ISBA_NTIME36
./script_exec_pgd_parall.sh "ECO_VGT_ISBA_NTIME36" $fname $2 $3
echo "########## ECO_VGT_ISBA_NTIME12"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 36/NTIME = 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_ISBA_NTIME12
./script_exec_pgd_parall.sh "ECO_VGT_ISBA_NTIME12" $fname $2 $3
echo "########## ECO_VGT_ISBA_NTIME2"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 12/NTIME = 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_ISBA_NTIME2
./script_exec_pgd_parall.sh "ECO_VGT_ISBA_NTIME2" $fname $2 $3
echo "########## ECO_VGT_ISBA_NTIME1"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 2/NTIME = 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_ISBA_NTIME1
./script_exec_pgd_parall.sh "ECO_VGT_ISBA_NTIME1" $fname $2 $3

# 3.3. juste le LAI
echo "########## ECO_LAI_ISBA_NTIME36"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_lai.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_LAI_ISBA_NTIME36
./script_exec_pgd_parall.sh "ECO_LAI_ISBA_NTIME36" $fname $2 $3
echo "########## ECO_LAI_ISBA_NTIME12"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 36/NTIME = 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_LAI_ISBA_NTIME12
./script_exec_pgd_parall.sh "ECO_LAI_ISBA_NTIME12" $fname $2 $3
echo "########## ECO_LAI_ISBA_NTIME2"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 12/NTIME = 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_LAI_ISBA_NTIME2
./script_exec_pgd_parall.sh "ECO_LAI_ISBA_NTIME2" $fname $2 $3
echo "########## ECO_LAI_ISBA_NTIME1"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 2/NTIME = 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_LAI_ISBA_NTIME1
./script_exec_pgd_parall.sh "ECO_LAI_ISBA_NTIME1" $fname $2 $3

# 3.4. juste les DG
echo "########## ECO_DG_ISBA_NTIME36"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_dg.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_DG_ISBA_NTIME36
./script_exec_pgd_parall.sh "ECO_DG_ISBA_NTIME36" $fname $2 $3
echo "########## ECO_DG_ISBA_NTIME12"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 36/NTIME = 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_DG_ISBA_NTIME12
./script_exec_pgd_parall.sh "ECO_DG_ISBA_NTIME12" $fname $2 $3
echo "########## ECO_DG_ISBA_NTIME2"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 12/NTIME = 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_DG_ISBA_NTIME2
./script_exec_pgd_parall.sh "ECO_DG_ISBA_NTIME2" $fname $2 $3
echo "########## ECO_DG_ISBA_NTIME1"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 2/NTIME = 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_DG_ISBA_NTIME1
./script_exec_pgd_parall.sh "ECO_DG_ISBA_NTIME1" $fname $2 $3

# 3.5. juste les H_TREE
echo "########## ECO_H_TREE_ISBA_NTIME36"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_h_tree.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_H_TREE_ISBA_NTIME36
./script_exec_pgd_parall.sh "ECO_H_TREE_ISBA_NTIME36" $fname $2 $3
echo "########## ECO_H_TREE_ISBA_NTIME12"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 36/NTIME = 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_H_TREE_ISBA_NTIME12
./script_exec_pgd_parall.sh "ECO_H_TREE_ISBA_NTIME12" $fname $2 $3
echo "########## ECO_H_TREE_ISBA_NTIME2"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 12/NTIME = 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_H_TREE_ISBA_NTIME2
./script_exec_pgd_parall.sh "ECO_H_TREE_ISBA_NTIME2" $fname $2 $3
echo "########## ECO_H_TREE_ISBA_NTIME1"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME = 2/NTIME = 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_H_TREE_ISBA_NTIME1
./script_exec_pgd_parall.sh "ECO_H_TREE_ISBA_NTIME1" $fname $2 $3

# 3.6. juste les RSMIN
echo "########## ECO_RSMIN_ISBA_NTIME36"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_rsmin.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_RSMIN_ISBA_NTIME36
./script_exec_pgd_parall.sh "ECO_RSMIN_ISBA_NTIME36" $fname $2 $3
echo "########## ECO_RSMIN_ISBA_NTIME12"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 36/NTIME\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_RSMIN_ISBA_NTIME12
./script_exec_pgd_parall.sh "ECO_RSMIN_ISBA_NTIME12" $fname $2 $3
echo "########## ECO_RSMIN_ISBA_NTIME2"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 12/NTIME\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_RSMIN_ISBA_NTIME2
./script_exec_pgd_parall.sh "ECO_RSMIN_ISBA_NTIME2" $fname $2 $3
echo "########## ECO_RSMIN_ISBA_NTIME1"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 2/NTIME\ =\ 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_RSMIN_ISBA_NTIME1
./script_exec_pgd_parall.sh "ECO_RSMIN_ISBA_NTIME1" $fname $2 $3

# 3.7. juste le VEG
echo "########## ECO_VEG_ISBA_NTIME36"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_veg.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_VEG_ISBA_NTIME36
./script_exec_pgd_parall.sh "ECO_VEG_ISBA_NTIME36" $fname $2 $3
echo "########## ECO_VEG_ISBA_NTIME12"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 36/NTIME\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VEG_ISBA_NTIME12
./script_exec_pgd_parall.sh "ECO_VEG_ISBA_NTIME12" $fname $2 $3
echo "########## ECO_VEG_ISBA_NTIME2"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 12/NTIME\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VEG_ISBA_NTIME2
./script_exec_pgd_parall.sh "ECO_VEG_ISBA_NTIME2" $fname $2 $3
echo "########## ECO_VEG_ISBA_NTIME1"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 2/NTIME\ =\ 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VEG_ISBA_NTIME1
./script_exec_pgd_parall.sh "ECO_VEG_ISBA_NTIME1" $fname $2 $3

# 3.8. juste le Z0
echo "########## ECO_Z0_ISBA_NTIME36"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_z0.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_Z0_ISBA_NTIME36
./script_exec_pgd_parall.sh "ECO_Z0_ISBA_NTIME36" $fname $2 $3
echo "########## ECO_Z0_ISBA_NTIME12"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 36/NTIME\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_Z0_ISBA_NTIME12
./script_exec_pgd_parall.sh "ECO_Z0_ISBA_NTIME12" $fname $2 $3
echo "########## ECO_Z0_ISBA_NTIME2"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 12/NTIME\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_Z0_ISBA_NTIME2
./script_exec_pgd_parall.sh "ECO_Z0_ISBA_NTIME2" $fname $2 $3
echo "########## ECO_Z0_ISBA_NTIME1"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 2/NTIME\ =\ 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_Z0_ISBA_NTIME1
./script_exec_pgd_parall.sh "ECO_Z0_ISBA_NTIME1" $fname $2 $3

# 3.9. juste VGT + LAI
echo "########## ECO_VGT_LAI_ISBA_NTIME36"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_vgt.sh "OPTIONS.nam"
./make_isba_adrien_lai.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_LAI_ISBA_NTIME36
./script_exec_pgd_parall.sh "ECO_VGT_LAI_ISBA_NTIME36" $fname $2 $3
echo "########## ECO_VGT_LAI_ISBA_NTIME12"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 36/NTIME\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_LAI_ISBA_NTIME12
./script_exec_pgd_parall.sh "ECO_VGT_LAI_ISBA_NTIME12" $fname $2 $3
echo "########## ECO_VGT_LAI_ISBA_NTIME2"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 12/NTIME\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_LAI_ISBA_NTIME2
./script_exec_pgd_parall.sh "ECO_VGT_LAI_ISBA_NTIME2" $fname $2 $3
echo "########## ECO_VGT_LAI_ISBA_NTIME1"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 2/NTIME\ =\ 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_LAI_ISBA_NTIME1
./script_exec_pgd_parall.sh "ECO_VGT_LAI_ISBA_NTIME1" $fname $2 $3

# 3.10. juste VGT + RSMIN
echo "########## ECO_VGT_RSMIN_ISBA_NTIME36"
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_vgt.sh "OPTIONS.nam"
./make_isba_adrien_rsmin.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_RSMIN_ISBA_NTIME36
./script_exec_pgd_parall.sh "ECO_VGT_RSMIN_ISBA_NTIME36" $fname $2 $3
echo "########## ECO_VGT_RSMIN_ISBA_NTIME12"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 36/NTIME\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_RSMIN_ISBA_NTIME12
./script_exec_pgd_parall.sh "ECO_VGT_RSMIN_ISBA_NTIME12" $fname $2 $3
echo "########## ECO_VGT_RSMIN_ISBA_NTIME2"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 12/NTIME\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_RSMIN_ISBA_NTIME2
./script_exec_pgd_parall.sh "ECO_VGT_RSMIN_ISBA_NTIME2" $fname $2 $3
echo "########## ECO_VGT_RSMIN_ISBA_NTIME1"
mv -f OPTIONS.nam OPTIONS.nam_save 
sed -e "s/NTIME\ =\ 2/NTIME\ =\ 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ECO_VGT_RSMIN_ISBA_NTIME1
./script_exec_pgd_parall.sh "ECO_VGT_RSMIN_ISBA_NTIME1" $fname $2 $3


# 4. CISBA = "DIF" NGROUND_LAYER=8 - 3 cas
mv -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_save 
sed -e "s/CISBA\ =\ \"2-L\"/CISBA\ =\ \"DIF\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW_GD\ =\ \"D95\"/CSNOW_GD\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW_LAYER\ =\ 1/NSNOW_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW_LAYER_GD\ =\ 1/NSNOW_LAYER_GD\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YPERM = \"\"/YPERM = \"perm_glo_10km\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSOILGRID = 0.01,0.1,1.1,2.4,3.9,5.2,6.1,7.9//g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND_LAYER\ =\ 2/NGROUND_LAYER\ =\ 14/g" OPTIONS.nam_save > OPTIONS.nam_DIF

# 4.1. 
echo "########## ECO_ISBA_DIF_DG_ROOTFRAC_DATA"
cp -f OPTIONS.nam_DIF OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_dif_rootfrac.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_ISBA_DIF_DG_ROOTFRAC_DATA
./script_exec_pgd_parall.sh "ECO_ISBA_DIF_DG_ROOTFRAC_DATA" $fname $2 $3

# 4.2. 
echo "########## ECO_ISBA_DIF_ROOTDEPTH_GROUNDDEPTH_EXT_LIN_DATA"
cp -f OPTIONS.nam_DIF OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_dif_rootdepth.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_ISBA_DIF_ROOTDEPTH_GROUNDDEPTH_EXT_LIN_DATA
./script_exec_pgd_parall.sh "ECO_ISBA_DIF_ROOTDEPTH_GROUNDDEPTH_EXT_LIN_DATA" $fname $2 $3

# 4.3.
echo "########## ECO_ISBA_DIF_DG_ROOTDEPTH_EXT_LIN_DATA"
cp -f OPTIONS.nam_DIF OPTIONS.nam
./make_garden_unif1.sh "OPTIONS.nam"
./make_isba_adrien_dif_root_ext.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_ECO_ISBA_DIF_DG_ROOTDEPTH_EXT_LIN_DATA
./script_exec_pgd_parall.sh "ECO_ISBA_DIF_DG_ROOTDEPTH_EXT_LIN_DATA" $fname $2 $3

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

