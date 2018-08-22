dir_run=$1
exec_new=$2
exec_old=$3

fname="INFO_PREP.txt"

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

####################TESTS POUR L'INITIALISATION PREP###########################################
## il vaut mieux activer le max d'options complexes; tester watflx et flake à chaque fois...

#namelist avec OPTIONS TEB activées
sed -e "s/LGARDEN\ =\ F/LGARDEN\ =\ T/g" TESTS/OPTIONS.nam_000 > OPTIONS.nam_OPT
mv -f OPTIONS.nam_OPT OPTIONS.nam_save
sed -e "s/CBEM\ =\ \"DEF\"/CBEM\ =\ \"BEM\"/g" OPTIONS.nam_save > OPTIONS.nam_OPT
mv -f OPTIONS.nam_OPT OPTIONS.nam_save
sed -e "s/LGREENROOF\ =\ F/LGREENROOF\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam_OPT
mv -f OPTIONS.nam_OPT OPTIONS.nam_save
sed -e "s/XUNIF_GREENROOF\ =\ 1.0E+20/XUNIF_GREENROOF\ =\ 0.4/g" OPTIONS.nam_save > OPTIONS.nam_OPT

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
sed -e "s/XUNIF_GARDEN = 0./XUNIF_GARDEN = 0.5/g" OPTIONS.nam_save > OPTIONS.nam_NO_ECOCLIMAP

#------------------------------------------------------------------------------

ln -s TESTS/PREP/FILES/* .

##########tests prep généraux 

############ TESTS AVEC UN FICHIER GRIB

#test grib aladin
echo "########## PREP_GENE_ALADIN"
cp -f OPTIONS.nam_NO_ECOCLIMAP OPTIONS.nam_save
sed -e "s/CFILE = \"\"/CFILE = \"aladin.AN.20100410.00.00\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILETYPE = \"\"/CFILETYPE = \"GRIB\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO_PREP = 2/NHALO_PREP = 35/g" OPTIONS.nam_save > OPTIONS.nam
cp -f TESTS/PREP/FILES/Params_config.txt_aladin Params_config.txt
cp -f OPTIONS.nam OPTIONS.nam_PREP_GENE_ALADIN
./script_exec_pgd_parall.sh "PREP_GENE_ALADIN" $fname $2 $3

#test grib arome
echo "########## PREP_GENE_AROME"
cp -f OPTIONS.nam_PREP_GENE_ALADIN OPTIONS.nam_save 
sed -e "s/aladin\.AN.20100410\.00\.00/arome\.AN\.20100310\.18\.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO_PREP = 35/NHALO_PREP = 76/g" OPTIONS.nam_save > OPTIONS.nam
cp -f TESTS/PREP/FILES/Params_config.txt_arome Params_config.txt
cp -f OPTIONS.nam OPTIONS.nam_PREP_GENE_AROME
./script_exec_pgd_parall.sh "PREP_GENE_AROME" $fname $2 $3

#test grib arpege
echo "########## PREP_GENE_ARPEGE"
cp -f OPTIONS.nam_PREP_GENE_ALADIN OPTIONS.nam_save 
sed -e "s/aladin\.AN.20100410\.00\.00/arpifs\.FC\.20100310\.18\.06/g" OPTIONS.nam_save > OPTIONS.nam
cp -f TESTS/PREP/FILES/Params_config.txt_arpifs Params_config.txt
cp -f OPTIONS.nam OPTIONS.nam_PREP_GENE_ARPEGE
./script_exec_pgd_parall.sh "PREP_GENE_ARPEGE" $fname $2 $3

#test grib ecmwf
echo "########## PREP_GENE_ECMWF"
cp -f OPTIONS.nam_PREP_GENE_ALADIN OPTIONS.nam_save 
sed -e "s/aladin\.AN.20100410\.00\.00/ecmwf\.OD\.20030814\.00\.06/g" OPTIONS.nam_save > OPTIONS.nam
cp -f TESTS/PREP/FILES/Params_config.txt_ecmwf Params_config.txt
cp -f OPTIONS.nam OPTIONS.nam_PREP_GENE_ECMWF
./script_exec_pgd_parall.sh "PREP_GENE_ECMWF" $fname $2 $3


#-------------------------------------------------------------------------

#test fichier PGD PREP CONF_PROJ ASCII version courante
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_EXTERN_IN
./make_prep_unif.sh "OPTIONS.nam_EXTERN_IN"
./make_prep_isba_unif.sh "OPTIONS.nam_EXTERN_IN"
./make_prep_isba_snow_unif.sh "OPTIONS.nam_EXTERN_IN"
./make_prep_teb_unif.sh "OPTIONS.nam_EXTERN_IN"
./make_prep_gd_gr_unif.sh "OPTIONS.nam_EXTERN_IN"
./make_prep_teb_snow_unif.sh "OPTIONS.nam_EXTERN_IN"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam_EXTERN_IN"
./make_prep_flake_unif.sh "OPTIONS.nam_EXTERN_IN"
cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1.0E+20/XSST_UNIF = 283.15/g" OPTIONS.nam_save > OPTIONS.nam_EXTERN_IN
cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/XTS_WATER_UNIF = 1.0E+20/XTS_WATER_UNIF = 282./g" OPTIONS.nam_save > OPTIONS.nam_EXTERN_IN


### tests avec des fichier EXTERN 
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .

echo "########## PREP_EXTE_FILE_ASCII_GENE_NEW"
cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_GENE_NEW" $fname $exec_new $exec_old "txt"
if [ ! -f PREP_BASE_OLD.txt ] 
then
  ./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_GENE_OLD" $fname $exec_new $exec_old "txt"
fi
if [ ! -f PREP_BASE_OLD.txt ] 
then
  ./script_prep_exte.sh "PREP_EXTE_FILE_EXTRAPOL_BILIN" $fname $exec_new $exec_old "txt"
fi


cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_save
sed -e "s/CFILE = \"\"/CFILE = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILETYPE = \"\"/CFILETYPE = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILEPGD = \"\"/CFILEPGD = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILEPGDTYPE = \"\"/CFILEPGDTYPE = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LGARDEN\ =\ T/LGARDEN\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LGREENROOF\ =\ T/LGREENROOF\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF\_GREENROOF\ =\ 0.4/XUNIF\_GREENROOF\ =\ 1.0E+20/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF\_GARDEN\ =\ 0.5/XUNIF\_GARDEN\ =\ 0./g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSURF_DIAG_ALBEDO = F/LSURF_DIAG_ALBEDO = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTERN

cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_GENE_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_GENE_NEW" $fname $2 $3

#test fichier PGD PREP CONF_PROJ ASCII version courante, avec extrapolation BILIN (Laure)
echo "########## PREP_EXTE_FILE_EXTRAPOL_BILIN"
cp OPTIONS.nam_PREP_EXTERN OPTIONS.nam_save
sed -e "s/NJMAX = 12/NJMAX = 20/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_EXTRAPOL_BILIN
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_EXTRAPOL_BILIN" $fname $2 $3

#test fichier PGD PREP CONF_PROJ ASCII version précédente
echo "########## PREP_EXTE_FILE_ASCII_GENE_OLD"
rm -f PGD_BASE_NEW.txt
rm -f PREP_BASE_NEW.txt
mv -f PGD_BASE_OLD.txt PGD_BASE.txt
mv -f PREP_BASE_OLD.txt PREP_BASE.txt
cp -f OPTIONS.nam_PREP_EXTE_FILE_ASCII_GENE_NEW OPTIONS.nam
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_GENE_OLD" $fname $2 $3
rm -f *_BASE*.*



#test fichier PGD PREP LONLAT_REG LFI version courante
cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"LONLAT REG\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LWRITE_EXTERN = F/LWRITE_EXTERN = T/g" OPTIONS.nam_save > OPTIONS.nam

echo "########## PREP_EXTE_FILE_LFI_GENE_NEW"
./script_prep_exte.sh "PREP_EXTE_FILE_LFI_GENE_NEW" $fname $exec_new $exec_old "lfi" 
if [ ! -f PREP_BASE_OLD.lfi ] 
then
  ./script_prep_exte.sh "PREP_EXTE_FILE_LFI_GENE_OLD" $fname $exec_new $exec_old "lfi" 
fi

cp -f OPTIONS.nam_PREP_EXTERN OPTIONS.nam_save
sed -e "s/CFILETYPE\ =\ \"ASCII\"/CFILETYPE\ =\ \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILEPGDTYPE\ =\ \"ASCII\"/CFILEPGDTYPE\ =\ \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_LFI_GENE_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_LFI_GENE_NEW" $fname $2 $3

#test fichier PGD PREP LONLAT_REG LFI version précédente
echo "########## PREP_EXTE_FILE_LFI_GENE_OLD"
rm -f PGD_BASE_NEW.lfi
rm -f PREP_BASE_NEW.lfi
mv -f PGD_BASE_OLD.lfi PGD_BASE.lfi
mv -f PREP_BASE_OLD.lfi PREP_BASE.lfi
cp -f OPTIONS.nam_PREP_EXTE_FILE_LFI_GENE_NEW OPTIONS.nam
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_LFI_GENE_OLD" $fname $2 $3
rm -f *_BASE*.*

#rajouter les cas avec NC et FA quand ce sera possible dans les deux version

#test fichier PGD PREP LONLAT_REG FA version courante
cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"FA\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"LONLAT REG\"/g" OPTIONS.nam_save > OPTIONS.nam

echo "########## PREP_EXTE_FILE_FA_GENE_NEW"
./script_prep_exte.sh "PREP_EXTE_FILE_FA_GENE_NEW" $fname $exec_new $exec_old "fa" 
if [ ! -f PREP_BASE_OLD.fa ] 
then
  ./script_prep_exte.sh "PREP_EXTE_FILE_FA_GENE_OLD" $fname $exec_new $exec_old "fa" 
fi
#
cp -f OPTIONS.nam_PREP_EXTERN OPTIONS.nam_save
sed -e "s/CFILETYPE\ =\ \"ASCII\"/CFILETYPE\ =\ \"FA\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILEPGDTYPE\ =\ \"ASCII\"/CFILEPGDTYPE\ =\ \"FA\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_FA_GENE_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_FA_GENE_NEW" $fname $2 $3

#test fichier PGD PREP LONLAT_REG LFI version précédente
echo "########## PREP_EXTE_FILE_FA_GENE_OLD"
rm -f PGD_BASE_NEW.fa
rm -f PREP_BASE_NEW.fa
mv -f PGD_BASE_OLD.fa PGD_BASE.fa
mv -f PREP_BASE_OLD.fa PREP_BASE.fa
cp -f OPTIONS.nam_PREP_EXTE_FILE_LFI_GENE_NEW OPTIONS.nam
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_FA_GENE_OLD" $fname $2 $3
rm -f *_BASE*.*

#------------------------------------------------------------

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_PREP_EXTERN_ISBA
./make_prep_unif.sh "OPTIONS.nam_PREP_EXTERN_ISBA"
./make_prep_teb_unif.sh "OPTIONS.nam_PREP_EXTERN_ISBA"
./make_prep_gd_gr_unif.sh "OPTIONS.nam_PREP_EXTERN_ISBA"
./make_prep_teb_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_ISBA"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_ISBA"
./make_prep_flake_unif.sh "OPTIONS.nam_PREP_EXTERN_ISBA"
cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1.0E+20/XSST_UNIF = 283.15/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_ISBA
cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam_save
sed -e "s/XTS_WATER_UNIF = 1.0E+20/XTS_WATER_UNIF = 282./g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_ISBA
cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam_save
sed -e "s/CFILEPGD_ISBA = \"\"/CFILEPGD_ISBA = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_ISBA
cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam_save
sed -e "s/CTYPEPGD = \"\"/CTYPEPGD = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_ISBA

cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam_save
sed -e "s/CFILE_ISBA = \"\"/CFILE_ISBA = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_ISBA_TOT
cp -f OPTIONS.nam_PREP_EXTERN_ISBA_TOT OPTIONS.nam_save
sed -e "s/CTYPE = \"\"/CTYPE = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_ISBA_TOT


##########tests PREP spécifiques ISBA
echo "########## PREP_EXTE_FILE_ASCII_ISBA_12P1P_NEW"
#test avec un fichier PREP 12 patchs en entrée, 1 patch en sortie
cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/NPATCH\ =\ 1/NPATCH\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CPHOTO\ =\ \"NON\"/CPHOTO\ =\ \"AST\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_ISBA_12P1P_NEW" $fname $exec_new $exec_old "txt"

cp -f OPTIONS.nam_PREP_EXTERN_ISBA_TOT OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO_PREP = 2/NHALO_PREP = 9/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_12P1P_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_12P1P_NEW" $fname $2 $3

echo "########## PREP_EXTE_FILE_ASCII_ISBA_12P19P_NEW"
cp -f OPTIONS.nam_PREP_EXTERN_ISBA_TOT OPTIONS.nam
cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/NPATCH\ =\ 1/NPATCH\ =\ 19/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO_PREP = 2/NHALO_PREP = 9/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_12P19P_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_12P19P_NEW" $fname $2 $3
rm -f *_BASE*.*


#test avec un fichier PREP 1 patch en entrée, 12 patchs en sortie
echo "########## PREP_EXTE_FILE_ASCII_ISBA_1P12P_NEW"

cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_ISBA_1P12P_NEW" $fname $exec_new $exec_old "txt" 
 
cp -f OPTIONS.nam_PREP_EXTERN_ISBA_TOT OPTIONS.nam_save
sed -e "s/NPATCH\ =\ 1/NPATCH\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CPHOTO\ =\ \"NON\"/CPHOTO\ =\ \"AST\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_1P12P_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_1P12P_NEW" $fname $2 $3

echo "########## PREP_EXTE_FILE_ASCII_ISBA_1P19P_NEW"
cp -f OPTIONS.nam_PREP_EXTERN_ISBA_TOT OPTIONS.nam_save
sed -e "s/NPATCH\ =\ 1/NPATCH\ =\ 19/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CPHOTO\ =\ \"NON\"/CPHOTO\ =\ \"NIT\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_1P19P_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_1P19P_NEW" $fname $2 $3


#test avec un fichier PREP 19patch en entrée, 12 patchs en sortie
echo "########## PREP_EXTE_FILE_ASCII_ISBA_19P12P_NEW"

cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/NPATCH\ =\ 1/NPATCH\ =\ 19/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CPHOTO\ =\ \"NON\"/CPHOTO\ =\ \"NIT\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_ISBA_19P12P_NEW" $fname $exec_new $exec_old "txt" 
 
cp -f OPTIONS.nam_PREP_EXTERN_ISBA_TOT OPTIONS.nam_save
sed -e "s/NPATCH\ =\ 1/NPATCH\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CPHOTO\ =\ \"NON\"/CPHOTO\ =\ \"AST\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO_PREP = 2/NHALO_PREP = 9/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_19P12P_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_19P12P_NEW" $fname $2 $3

echo "########## PREP_EXTE_FILE_ASCII_ISBA_19P1P_NEW"
cp -f OPTIONS.nam_PREP_EXTERN_ISBA_TOT OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO_PREP = 2/NHALO_PREP = 9/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_19P1P_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_19P1P_NEW" $fname $2 $3



echo "########## PREP_EXTE_FILE_ASCII_ISBA_TG_NEW "
cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam
./make_prep_isba_tg_file.sh "OPTIONS.nam"
./make_prep_isba_snow_unif.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_TG_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_TG_NEW" $fname $2 $3


echo "########## PREP_EXTE_FILE_ASCII_ISBA_HUG_NEW "
cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam
./make_prep_isba_wg_file.sh "OPTIONS.nam"
./make_prep_isba_snow_unif.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO_PREP\ \=\ 2/NHALO_PREP\ \=\ 9/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_HUG_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_HUG_NEW" $fname $2 $3
rm -f *_BASE*.*

#test avec fichiers ASCLLV pour TG et HUG
echo "########## PREP_ASCLLV_FILE_ISBA_TG_HUG "
cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam
./make_prep_isba_data.sh "OPTIONS.nam"
./make_prep_isba_snow_unif.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO\ \=\ 2/NHALO\ \=\ 15/g" OPTIONS.nam_save > OPTIONS.nam
ln -s TESTS/PREP/ISBA/*.dat .
cp -f OPTIONS.nam OPTIONS.nam_PREP_ASCLLV_FILE_ISBA_TG_HUG
./script_exec_pgd_parall.sh "PREP_ASCLLV_FILE_ISBA_TG_HUG" $fname $2 $3

#test avec fichier NETCDF pour TG
echo "########## PREP_NETCDF_FILE_ISBA_TG_NEW "
cp -f OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_TG_NEW OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG = \"PREP_BASE\"/CFILE_TG = \"TG_nc.nc\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_TG = \"ASCII\"/CTYPE_TG = \"NETCDF\"/g" OPTIONS.nam_save > OPTIONS.nam
./make_prep_isba_snow_unif.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO\ \=\ 2/NHALO\ \=\ 15/g" OPTIONS.nam_save > OPTIONS.nam
ln -s TESTS/PREP/ISBA/*.nc .
cp -f OPTIONS.nam OPTIONS.nam_PREP_NETCDF_FILE_ISBA_TG_NEW
./script_exec_pgd_parall.sh "PREP_NETCDF_FILE_ISBA_TG_NEW" $fname $2 $3

#test en changeant CISBA1 : entrée: 2-L sortie: DIF
echo "########## PREP_EXTE_FILE_ASCII_ISBA_2LDIF8_NEW "
cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_ISBA_2LDIF8_NEW" $fname $exec_new $exec_old "txt"

cp -f OPTIONS.nam_PREP_EXTERN_ISBA_TOT OPTIONS.nam_save
sed -e "s/CISBA\ =\ \"2-L\"/CISBA\ =\ \"DIF\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW_GD\ =\ \"D95\"/CSNOW_GD\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER_GD\ =\ 1/NSNOW\_LAYER_GD\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND\_LAYER\ =\ 2/NGROUND\_LAYER\ =\ 8/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_2LDIF8_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_2LDIF8_NEW" $fname $2 $3
rm -f *_BASE*.*

#test en changeant CISBA2 : entrée: DIF sortie: 3-L 
echo "########## PREP_EXTE_FILE_ASCII_ISBA_DIF83L_NEW "
cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CISBA\ =\ \"2-L\"/CISBA\ =\ \"DIF\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW_GD\ =\ \"D95\"/CSNOW_GD\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER_GD\ =\ 1/NSNOW\_LAYER_GD\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND\_LAYER\ =\ 2/NGROUND\_LAYER\ =\ 8/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_PERM\ =\ 10\./XUNIF_PERM\ =\ 1E+20/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_ISBA_DIF83L_NEW" $fname $exec_new $exec_old "txt"

cp -f OPTIONS.nam_PREP_EXTERN_ISBA_TOT OPTIONS.nam_save
sed -e "s/CISBA\ =\ \"2-L\"/CISBA\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND\_LAYER\ =\ 2/NGROUND\_LAYER\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_DIF83L_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_DIF83L_NEW" $fname $2 $3
rm -f *_BASE*.*

#-----------------------------------------------------------------------

##########tests PREP spécifiques ISBA_SNOW

#un fichier PREP exprès pour la neige: CROCUS->3-L
echo "########## PREP_EXTE_FILE_ASCII_ISBA_SNOW_CRO83L_NEW "

cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"CRO\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 8/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_ISBA_SNOW_CRO83L_NEW" $fname $exec_new $exec_old "txt"

cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam
./make_prep_isba_unif.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_SNOW = \"\"/CFILE_SNOW = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_SNOW = \"\"/CTYPE_SNOW = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILEPGD_SNOW = \"\"/CFILEPGD_SNOW = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPEPGD_SNOW = \"\"/CTYPEPGD_SNOW = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTERN_SNOW_ISBA

cp -f OPTIONS.nam_PREP_EXTERN_SNOW_ISBA OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_SNOW_CRO83L12_NEW
#ne fonctionne pas dans le trunk (bug corrigé dans NEW_PREP), à décommenter après
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_SNOW_CRO83L12_NEW" $fname $2 $3

cp -f OPTIONS.nam_PREP_EXTERN_SNOW_ISBA OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_SNOW_CRO83L3_NEW
#ne fonctionne pas dans le trunk (bug corrigé dans NEW_PREP), à décommenter après
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_SNOW_CRO83L3_NEW" $fname $2 $3
rm -f *_BASE*.*



#un fichier PREP exprès pour la neige: 3-L->D95
echo "########## PREP_EXTE_FILE_ASCII_ISBA_SNOW_3LD95_NEW "

cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_ISBA_SNOW_3LD95_NEW" $fname $exec_new $exec_old "txt"

cp -f OPTIONS.nam_PREP_EXTERN_SNOW_ISBA OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_SNOW_3LD95_NEW
#ne fonctionne pas dans le trunk (bug corrigé dans NEW_PREP), à décommenter après
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_SNOW_3LD95_NEW" $fname $2 $3
rm -f *_BASE*.*


#un fichier PREP exprès pour la neige: 3-L->EBA
echo "########## PREP_EXTE_FILE_ASCII_ISBA_SNOW_3LEBA_NEW "
cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_ISBA_SNOW_3LEBA_NEW" $fname $exec_new $exec_old "txt"

cp -f OPTIONS.nam_PREP_EXTERN_SNOW_ISBA OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"EBA\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_ISBA_SNOW_3LEBA_NEW
#ne fonctionne pas dans le trunk (bug corrigé dans NEW_PREP), à décommenter après
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_ISBA_SNOW_3LEBA_NEW" $fname $2 $3
rm -f *_BASE*.*

#LSNOW_IDEAL=T CSNOW = CRO
echo "########## PREP_ISBA_SNOW_CRO_IDEAL "
cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam
./make_prep_isba_unif.sh "OPTIONS.nam"
./make_prep_isba_snow_ideal.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSNOW_IDEAL = F/LSNOW_IDEAL = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW = \"D95\"/CSNOW = \"CRO\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW_LAYER = 1/NSNOW_LAYER = 6/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_ISBA_SNOW_CRO_IDEAL
./script_exec_pgd_parall.sh "PREP_ISBA_SNOW_CRO_IDEAL" $fname $2 $3

#LSWEMAX = T
echo "########## PREP_ISBA_LSWEMAX "
cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam
./make_prep_isba_unif.sh "OPTIONS.nam"
./make_prep_isba_snow_ideal.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSWEMAX = F/LSWEMAX = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_ISBA_LSWEMAX
./script_exec_pgd_parall.sh "PREP_ISBA_LSWEMAX" $fname $2 $3

#LSWEMAX = T XSWEMAX = 10.
echo "########## PREP_ISBA_LSWEMAX_XSWEMAX_10 "
cp -f OPTIONS.nam_PREP_ISBA_LSWEMAX OPTIONS.nam_save
sed -e "s/XSWEMAX = 500./XSWEMAX = 10./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_ISBA_LSWEMAX_XSWEMAX_10
./script_exec_pgd_parall.sh "PREP_ISBA_LSWEMAX_XSWEMAX_10" $fname $2 $3

#ZSNOW
echo "########## PREP_ISBA_SNOW_CRO_ZSNOW"
cp -f OPTIONS.nam_PREP_EXTERN_ISBA OPTIONS.nam
./make_prep_isba_unif.sh "OPTIONS.nam"
./make_prep_isba_snow_ideal.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW = \"D95\"/CSNOW = \"CRO\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW_LAYER = 1/NSNOW_LAYER = 6/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSNOW_IDEAL = F/LSNOW_IDEAL = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSNOW_PREP_PERM = T/LSNOW_PREP_PERM = F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XWSNOW = 10.,50.,100.,840.,9000.,40000./XWSNOW = 1.0E+20/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XZSNOW = 1.0E+20/XZSNOW = 0.01,0.05,0.1,0.2,0.2,0.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_ISBA_SNOW_CRO_ZSNOW
./script_exec_pgd_parall.sh "PREP_ISBA_SNOW_CRO_ZSNOW" $fname $2 $3

#LWCSNOW
echo "########## PREP_ISBA_SNOW_CRO_LWCSNOW"
cp -f OPTIONS.nam_PREP_ISBA_SNOW_CRO_ZSNOW OPTIONS.nam_save
sed -e "s/XLWCSNOW = 0./XLWCSNOW = 5,0,0,0,15,0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_ISBA_SNOW_CRO_LWCSNOW
./script_exec_pgd_parall.sh "PREP_ISBA_SNOW_CRO_LWCSNOW" $fname $2 $3

#------------------------------------------------------------------------

##########tests PREP spécifiques TEB

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_PREP_EXTERN_TEB
./make_prep_unif.sh "OPTIONS.nam_PREP_EXTERN_TEB"
./make_prep_isba_unif.sh "OPTIONS.nam_PREP_EXTERN_TEB"
./make_prep_isba_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_TEB"
./make_prep_flake_unif.sh "OPTIONS.nam_PREP_EXTERN_TEB"
./make_prep_gd_gr_unif.sh "OPTIONS.nam_PREP_EXTERN_TEB"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_TEB"
cp -f OPTIONS.nam_PREP_EXTERN_TEB OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1.0E+20/XSST_UNIF = 283.15/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_TEB
cp -f OPTIONS.nam_PREP_EXTERN_TEB OPTIONS.nam_save
sed -e "s/XTS_WATER_UNIF = 1.0E+20/XTS_WATER_UNIF = 282./g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_TEB
cp -f OPTIONS.nam_PREP_EXTERN_TEB OPTIONS.nam_save
sed -e "s/CFILEPGD_TEB = \"\"/CFILEPGD_TEB = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_TEB
cp -f OPTIONS.nam_PREP_EXTERN_TEB OPTIONS.nam_save
sed -e "s/CTYPEPGD = \"\"/CTYPEPGD = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_TEB

cp -f OPTIONS.nam_PREP_EXTERN_TEB OPTIONS.nam_save
sed -e "s/CFILE_TEB = \"\"/CFILE_TEB = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_TEB_TOT
cp -f OPTIONS.nam_PREP_EXTERN_TEB_TOT OPTIONS.nam_save
sed -e "s/CTYPE = \"\"/CTYPE = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_TEB_TOT


#test avec un fichier PREP sans TEB
echo "########## PREP_EXTE_FILE_ASCII_SANS_TEB_AVEC"

cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTOWN\ =\ \"TEB\"/CTOWN\ =\ \"NONE\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_SANS_TEB_AVEC" $fname $exec_new $exec_old "txt" 

cp -f OPTIONS.nam_PREP_EXTERN_TEB_TOT OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_SANS_TEB_AVEC
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_SANS_TEB_AVEC" $fname $2 $3
rm -f *_BASE*.*

#test avec un fichier PREP avec des nbs de couches >
echo "########## PREP_EXTE_FILE_ASCII_TEB_NBL-_NEW "

cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_TEB_NBL-_NEW" $fname $exec_new $exec_old "txt"
if [ ! -f PREP_BASE_NEW.txt ]
then
  ./script_prep_exte.sh "PREP_ASCII_FILE_TEB_TS" $fname $exec_new $exec_old "txt"
fi
if [ ! -f PREP_BASE_NEW.txt ]
then
  ./script_prep_exte.sh "PREP_ASCII_FILE_TEB_WS" $fname $exec_new $exec_old "txt"
fi

cp -f OPTIONS.nam_PREP_EXTERN_TEB_TOT OPTIONS.nam_save
sed -e "s/NROAD_LAYER\ =\ 5/NROAD_LAYER\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NROOF_LAYER\ =\ 5/NROOF_LAYER\ =\ 4/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NWALL_LAYER\ =\ 5/NWALL_LAYER\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NFLOOR_LAYER\ =\ 5/NFLOOR_LAYER\ =\ 4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_TEB_NBL-_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_TEB_NBL-_NEW" $fname $2 $3

#fichier TS ascii pour TEB
echo "########## PREP_ASCII_FILE_TEB_TS "
cp -f OPTIONS.nam_PREP_EXTERN_TEB OPTIONS.nam
./make_prep_teb_ts_file.sh "OPTIONS.nam"
./make_prep_teb_snow_unif.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_PREP_ASCII_FILE_TEB_TS
./script_exec_pgd_parall.sh "PREP_ASCII_FILE_TEB_TS" $fname $2 $3

#fichier WS ascii pour TEB
echo "########## PREP_ASCII_FILE_TEB_WS "
cp -f OPTIONS.nam_PREP_EXTERN_TEB OPTIONS.nam
./make_prep_teb_ws_file.sh "OPTIONS.nam"
./make_prep_teb_snow_unif.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_PREP_ASCII_FILE_TEB_WS
./script_exec_pgd_parall.sh "PREP_ASCII_FILE_TEB_WS" $fname $2 $3
rm -f *_BASE*.*

#test avec un fichier PREP avec des nbs de couches <
echo "########## PREP_EXTE_FILE_ASCII_TEB_NBL+_NEW "

cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NROAD_LAYER\ =\ 5/NROAD_LAYER\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NROOF_LAYER\ =\ 5/NROOF_LAYER\ =\ 4/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NWALL_LAYER\ =\ 5/NWALL_LAYER\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NFLOOR_LAYER\ =\ 5/NFLOOR_LAYER\ =\ 4/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_TEB_NBL+_NEW" $fname $exec_new $exec_old "txt"

cp -f OPTIONS.nam_PREP_EXTERN_TEB_TOT OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_TEB_NBL+_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_TEB_NBL+_NEW" $fname $2 $3
rm -f *_BASE*.*

##########tests PREP spécifiques TEB_SNOW
echo "########## PREP_TEB_SNOW_IDEAL "
cp -f OPTIONS.nam_PREP_EXTERN_TEB OPTIONS.nam
./make_prep_teb_unif.sh "OPTIONS.nam"
./make_prep_teb_snow_unif.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSNOW_IDEAL_TEB = F/LSNOW_IDEAL_TEB = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_TEB_SNOW_IDEAL
./script_exec_pgd_parall.sh "PREP_TEB_SNOW_IDEAL" $fname $2 $3
rm -f *_BASE*.*


##########tests PREP spécifiques GARDEN et GREENROOF
#### attention: actuellement la fonctionnalité fichier PREP d'entrée pour GARDEN et GREENROOF ne fonctionne pas à cause de la neige

################ce cas ne fonctionne pas dans trunk, attendre la prochaine version pour decommenter
cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_GD_GR_NEW" $fname $exec_new $exec_old "txt" 

#prep garden greenroof PREP ascii (mais initialisation de la neige)

echo "########## PREP_EXTE_FILE_ASCII_GD_GR_NEW "
cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_PREP_EXTERN_GD_GR
./make_prep_unif.sh "OPTIONS.nam_PREP_EXTERN_GD_GR"
./make_prep_isba_unif.sh "OPTIONS.nam_PREP_EXTERN_GD_GR"
./make_prep_isba_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_GD_GR"
./make_prep_teb_unif.sh "OPTIONS.nam_PREP_EXTERN_GD_GR"
./make_prep_teb_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_GD_GR"
./make_prep_flake_unif.sh "OPTIONS.nam_PREP_EXTERN_GD_GR"
cp -f OPTIONS.nam_PREP_EXTERN_GD_GR OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1.0E+20/XSST_UNIF = 283.15/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_GD_GR
cp -f OPTIONS.nam_PREP_EXTERN_GD_GR OPTIONS.nam_save
sed -e "s/XTS_WATER_UNIF = 1.0E+20/XTS_WATER_UNIF = 282./g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_GD_GR
cp -f OPTIONS.nam_PREP_EXTERN_GD_GR OPTIONS.nam_save
sed -e "s/CFILEPGD_GD = \"\"/CFILEPGD_GD = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_GD_GR
cp -f OPTIONS.nam_PREP_EXTERN_GD_GR OPTIONS.nam_save
sed -e "s/CFILEPGD_GR = \"\"/CFILEPGD_GR = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_GD_GR
cp -f OPTIONS.nam_PREP_EXTERN_GD_GR OPTIONS.nam_save
sed -e "s/CTYPEPGD = \"\"/CTYPEPGD = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_GD_GR

cp -f OPTIONS.nam_PREP_EXTERN_GD_GR OPTIONS.nam_save
sed -e "s/CFILE_GD = \"\"/CFILE_GD = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_GD_GR_TOT
cp -f OPTIONS.nam_PREP_EXTERN_GD_GR_TOT OPTIONS.nam_save
sed -e "s/CFILE_GR = \"\"/CFILE_GR = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_GD_GR_TOT
cp -f OPTIONS.nam_PREP_EXTERN_GD_GR_TOT OPTIONS.nam_save
sed -e "s/CTYPE = \"\"/CTYPE = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_GD_GR_TOT

cp -f OPTIONS.nam_PREP_EXTERN_GD_GR_TOT OPTIONS.nam_save
sed -e "s/NHALO_PREP = 2/NHALO_PREP = 3/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_GD_GR_TOT

cp -f OPTIONS.nam_PREP_EXTERN_GD_GR_TOT OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_GD_GR_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_GD_GR_NEW" $fname $2 $3

##################ce cas ne fonctionne pas dans trunk, attendre la prochaine version pour decommenter


#prep garden greenroof TG ascii
echo "########## PREP_EXTE_FILE_ASCII_GD_GR_TG_NEW "
cp -f OPTIONS.nam_PREP_EXTERN_GD_GR OPTIONS.nam
./make_prep_gd_gr_tg_file.sh "OPTIONS.nam"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_GD_GR_TG_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_GD_GR_TG_NEW" $fname $2 $3

#prep garden greenroof HUG ascii
echo "########## PREP_EXTE_FILE_ASCII_GD_GR_HUG_NEW "
cp -f OPTIONS.nam_PREP_EXTERN_GD_GR OPTIONS.nam
./make_prep_gd_gr_wg_file.sh "OPTIONS.nam"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO_PREP\ \=\ 2/NHALO_PREP\ \=\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_GD_GR_HUG_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_GD_GR_HUG_NEW" $fname $2 $3

rm -f *_BASE*.*

#test avec fichiers ASCLLV pour TG et HUG
echo "########## PREP_ASCLLV_FILE_GD_GR_TG_HUG "
cp -f OPTIONS.nam_PREP_EXTERN_GD_GR OPTIONS.nam
./make_prep_gd_gr_data.sh "OPTIONS.nam"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO\ \=\ 2/NHALO\ \=\ 15/g" OPTIONS.nam_save > OPTIONS.nam
ln -s TESTS/PREP/ISBA/*.dat .
cp -f OPTIONS.nam OPTIONS.nam_PREP_ASCLLV_FILE_GD_GR_TG_HUG
./script_exec_pgd_parall.sh "PREP_ASCLLV_FILE_GD_GR_TG_HUG" $fname $2 $3

#----------------------------------------------------------------------------------------

##########tests PREP spécifiques SEAFLUX

#test avec un fichier PREP 
echo "########## PREP_EXTE_FILE_ASCII_SEAFLUX_NEW "

cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_SEAFLUX_NEW" $fname $exec_new $exec_old "txt"

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_PREP_EXTERN_SEA
./make_prep_unif.sh "OPTIONS.nam_PREP_EXTERN_SEA"
./make_prep_isba_unif.sh "OPTIONS.nam_PREP_EXTERN_SEA"
./make_prep_isba_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_SEA"
./make_prep_teb_unif.sh "OPTIONS.nam_PREP_EXTERN_SEA"
./make_prep_teb_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_SEA"
./make_prep_gd_gr_unif.sh "OPTIONS.nam_PREP_EXTERN_SEA"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_SEA"
./make_prep_flake_unif.sh "OPTIONS.nam_PREP_EXTERN_SEA"
cp -f OPTIONS.nam_PREP_EXTERN_SEA OPTIONS.nam_save
sed -e "s/XTS_WATER_UNIF = 1.0E+20/XTS_WATER_UNIF = 282./g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_SEA
cp -f OPTIONS.nam_PREP_EXTERN_SEA OPTIONS.nam_save
sed -e "s/CFILEPGD_SEAFLX = \"\"/CFILEPGD_SEAFLX = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_SEA
cp -f OPTIONS.nam_PREP_EXTERN_SEA OPTIONS.nam_save
sed -e "s/CTYPEPGD = \"\"/CTYPEPGD = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_SEA

cp -f OPTIONS.nam_PREP_EXTERN_SEA OPTIONS.nam_save
sed -e "s/CFILE_SEAFLX = \"\"/CFILE_SEAFLX = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_SEA_TOT
cp -f OPTIONS.nam_PREP_EXTERN_SEA_TOT OPTIONS.nam_save
sed -e "s/CTYPE_SEAFLX = \"\"/CTYPE_SEAFLX = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_SEA_TOT

cp -f OPTIONS.nam_PREP_EXTERN_SEA_TOT OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_SEAFLUX_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_SEAFLUX_NEW" $fname $2 $3

rm -f *_BASE*.*

#avec un fichier NETCDF
echo "########## PREP_NETCDF_FILE_SEAFLUX_NEW "
ln -sf TESTS/PREP/FILES/mercator_20031203.nc .

cp -f OPTIONS.nam_PREP_EXTERN_SEA OPTIONS.nam_save
sed -e "s/CFILE_SEAFLX = \"\"/CFILE_SEAFLX = \"mercator_20031203.nc\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_SEAFLX = \"\"/CTYPE_SEAFLX = \"NETCDF\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLAT0\ =\ 45\./XLAT0\ =\ 43\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATCEN\ =\ 45\./XLATCEN\ =\ 43\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO_PREP = 2/NHALO_PREP = 34/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_NETCDF_FILE_SEAFLUX_NEW
./script_exec_pgd_parall.sh "PREP_NETCDF_FILE_SEAFLUX_NEW" $fname $2 $3

rm -f mercator_20031203.nc

##########tests PREP spécifiques WATFLUX
echo "########## PREP_EXTE_FILE_ASCII_WATFLUX_NEW "

cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_WATFLUX_NEW" $fname $exec_new $exec_old "txt" 

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_PREP_EXTERN_WATER
./make_prep_unif.sh "OPTIONS.nam_PREP_EXTERN_WATER"
./make_prep_isba_unif.sh "OPTIONS.nam_PREP_EXTERN_WATER"
./make_prep_isba_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_WATER"
./make_prep_teb_unif.sh "OPTIONS.nam_PREP_EXTERN_WATER"
./make_prep_teb_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_WATER"
./make_prep_gd_gr_unif.sh "OPTIONS.nam_PREP_EXTERN_WATER"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_WATER"
./make_prep_flake_unif.sh "OPTIONS.nam_PREP_EXTERN_WATER"
cp -f OPTIONS.nam_PREP_EXTERN_WATER OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1.0E+20/XSST_UNIF = 283.15/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_WATER
cp -f OPTIONS.nam_PREP_EXTERN_WATER OPTIONS.nam_save
sed -e "s/CFILEPGD_WATFLX = \"\"/CFILEPGD_WATFLX = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_WATER
cp -f OPTIONS.nam_PREP_EXTERN_WATER OPTIONS.nam_save
sed -e "s/CTYPEPGD = \"\"/CTYPEPGD = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_WATER

cp -f OPTIONS.nam_PREP_EXTERN_WATER OPTIONS.nam_save
sed -e "s/CFILE_WATFLX = \"\"/CFILE_WATFLX = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_WATER_TOT
cp -f OPTIONS.nam_PREP_EXTERN_WATER_TOT OPTIONS.nam_save
sed -e "s/CTYPE = \"\"/CTYPE = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_WATER_TOT

cp -f OPTIONS.nam_PREP_EXTERN_WATER_TOT OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_WATFLUX_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_WATFLUX_NEW" $fname $2 $3
rm -f *_BASE*.*

##########tests PREP spécifiques FLAKE

##avec un fichier PGD
echo "########## PREP_EXTE_FILE_ASCII_FLAKE_NEW "

cp -f OPTIONS.nam_EXTERN_IN OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CWATER\ =\ \"WATFLX\"/CWATER\ =\ \"FLAKE\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_prep_exte.sh "PREP_EXTE_FILE_ASCII_FLAKE_NEW" $fname $exec_new $exec_old "txt" 

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_PREP_EXTERN_FLAKE
./make_prep_unif.sh "OPTIONS.nam_PREP_EXTERN_FLAKE"
./make_prep_isba_unif.sh "OPTIONS.nam_PREP_EXTERN_FLAKE"
./make_prep_isba_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_FLAKE"
./make_prep_teb_unif.sh "OPTIONS.nam_PREP_EXTERN_FLAKE"
./make_prep_teb_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_FLAKE"
./make_prep_gd_gr_unif.sh "OPTIONS.nam_PREP_EXTERN_FLAKE"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam_PREP_EXTERN_FLAKE"
cp -f OPTIONS.nam_PREP_EXTERN_FLAKE OPTIONS.nam_save
sed -e "s/XTS_WATER_UNIF = 1.0E+20/XTS_WATER_UNIF = 282./g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_FLAKE
cp -f OPTIONS.nam_PREP_EXTERN_FLAKE OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1.0E+20/XSST_UNIF = 283.15/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_FLAKE
cp -f OPTIONS.nam_PREP_EXTERN_FLAKE OPTIONS.nam_save
sed -e "s/CFILEPGD_FLAKE = \"\"/CFILEPGD_FLAKE = \"PGD_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_FLAKE
cp -f OPTIONS.nam_PREP_EXTERN_FLAKE OPTIONS.nam_save
sed -e "s/CTYPEPGD = \"\"/CTYPEPGD = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_FLAKE

cp -f OPTIONS.nam_PREP_EXTERN_FLAKE OPTIONS.nam_save
sed -e "s/CFILE_FLAKE = \"\"/CFILE_FLAKE = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_FLAKE_TOT
cp -f OPTIONS.nam_PREP_EXTERN_FLAKE_TOT OPTIONS.nam_save
sed -e "s/CTYPE = \"\"/CTYPE = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam_PREP_EXTERN_FLAKE_TOT

cp -f OPTIONS.nam_PREP_EXTERN_FLAKE_TOT OPTIONS.nam_save
sed -e "s/CWATER\ =\ \"WATFLX\"/CWATER\ =\ \"FLAKE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_EXTE_FILE_ASCII_FLAKE_NEW
./script_exec_pgd_parall.sh "PREP_EXTE_FILE_ASCII_FLAKE_NEW" $fname $2 $3


##avec la base de données de lacs
echo "########## PREP_FLAKE_CLIM_LAKE "

cp -f OPTIONS.nam_PREP_EXTERN_FLAKE OPTIONS.nam_save
sed -e "s/CWATER\ =\ \"WATFLX\"/CWATER\ =\ \"FLAKE\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LCLIM\_LAKE\ \=\ F/LCLIM\_LAKE\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PREP_FLAKE_CLIM_LAKE
./script_exec_pgd_parall.sh "PREP_FLAKE_CLIM_LAKE" $fname $2 $3

for file in TESTS/PREP/FILES/*
do
	file2=$(basename $file)
	rm -f $file2
done

for file in TESTS/PREP/ISBA/*
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

rm -f Params_config.txt .
rm -f Forc*.txt .

