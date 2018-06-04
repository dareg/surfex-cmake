dir_run=$1
exec_new=$2
exec_old=$3

fname="INFO_OFFLINE.txt"

rm -f $fname
touch $fname

cd $1

cp -f TESTS/script_offline_exte.sh .
cp -f TESTS/script_exec.sh .
cp -f TESTS/script_exec_restart.sh .
cp -f TESTS/script_exec_restart_parall.sh .
cp -f TESTS/script_to_old.sh .
cp -f TESTS/script_exec_all_parall.sh .
cp -f TESTS/script_exec_parall.sh .

rm -f Forc_*.txt
rm -f Params_config.txt

cp -f TESTS/PREP/make*.sh .
cp -f TESTS/PGD/make*.sh .

cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
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

#ZS uniforme
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/XUNIF_ZS = 1.0E+20/XUNIF_ZS = 250./g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF

#recup d'un cas ECOCLIMAP 
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam
sed -e "s/LECOCLIMAP\ \=\ F/LECOCLIMAP\ \=\ T/g" OPTIONS.nam > OPTIONS.nam_ECOCLIMAP

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_base1

##############################################################
#pour faire les restart, on active les schémas de TEB et FLAKE

#ISBA: 3-L + SNOWCRO
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/CISBA\ =\ \"2-L\"/CISBA\ =\ \"DIF\"/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/NGROUND\_LAYER\ =\ 2/NGROUND\_LAYER\ =\ 8/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam_base

#CPHOTO='NIT' CHORT=SGH, NPATCH==12, LGLACIER
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CPHOTO\ =\ \"NON\"/CPHOTO\ =\ \"NIT\"/g" OPTIONS.nam_save > OPTIONS.nam_base
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/NPATCH\ =\ 1/NPATCH\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam_base
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CHORT\ =\ \"DEF\"/CHORT\ =\ \"SGH\"/g" OPTIONS.nam_save > OPTIONS.nam_base
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/LGLACIER\ =\ F/LGLACIER\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam_base


#TEB: on active tout
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/NTEB_PATCH\ =\ 1/NTEB_PATCH\ =\ 4/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CROAD\_DIR =\ \"UNIF\"/CROAD\_DIR\ =\ \"ORIE\"/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CWALL\_OPT\ =\ \"UNIF\"/CWALL\_OPT\ =\ \"TWO\"/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CBEM\ =\ \"DEF\"/CBEM\ =\ \"BEM\"/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/LGARDEN\ =\ F/LGARDEN\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/LGREENROOF\ =\ F/LGREENROOF\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/XUNIF\_GREENROOF\ =\ 1.0E+20/XUNIF\_GREENROOF\ =\ 0.7/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CSNOW_GD\ =\ \"D95\"/CSNOW_GD\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER_GD\ =\ 1/NSNOW\_LAYER_GD\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam_base
./make_greenroof_unif1.sh "OPTIONS.nam_base"

#on active FLAKE 
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CWATER\ =\ \"WATFLX\"/CWATER\ =\ \"FLAKE\"/g" OPTIONS.nam_save > OPTIONS.nam_base2

##############################################################

#noms des fichiers PGD PREP SURFOUT
echo "########### NOM_FICHIERS"
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/CPGDFILE\ =\ \"PGD\"/CPGDFILE\ =\ \"PGD\_TEST\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CPREPFILE\ =\ \"PREP\"/CPREPFILE\ =\ \"PREP\_TEST\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSURFFILE\ =\ \"SURFOUT\"/CSURFFILE\ =\ \"SURFOUT\_TEST\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NOMS_FICHIERS
#./script_exec.sh "NOM_FICHIERS" $fname $2 $3
rm -f PGD_TEST*
rm -f PREP_TEST*
rm -f SURFOUT_TEST*

echo "########### PRINT_RESTART_INQUIRE"
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/LPRINT\ =\ \T/LPRINT\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LRESTART\ =\ \T/LRESTART\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LINQUIRE\ =\ \T/LINQUIRE\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_PRINT_RESTART_INQUIRE
#./script_exec.sh "PRINT_RESTART_INQUIRE" $fname $2 $3

echo "########### XSTEP_SURF_OUTPUT "
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/XTSTEP\ =\ 300\./XTSTEP\ =\ 200./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTSTEP\_OUTPUT\ =\ 10800\./XTSTEP\_OUTPUT\ =\ 50400./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_XSTEP_SURF_OUTPUT
#./script_exec.sh "XSTEP_SURF_OUTPUT" $fname $2 $3


##############################################################

echo "########### LADAPT_SW "
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/LADAPT_SW = F/LADAPT_SW = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LADAPT_SW
#./script_exec.sh "LADAPT_SW" $fname $2 $3


echo "########### LT2MMW "
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/LT2MMW = F/LT2MMW = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LT2MMW
#./script_exec.sh "LT2MMW" $fname $2 $3

#CSURF_FILETYPE + RESTART + vieux fichier PGD et PREP

#NC
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
echo "########### CSURF_FILETYPE "
cp -f OPTIONS.nam_base2 OPTIONS.nam
#./script_offline_exte.sh "SURF_FILETYPE_NC_RESTART" $fname $exec_new $exec_old "nc"
#if [ ! -f PREP_OLD.nc ] 
#then
#  ./script_offline_exte.sh "SURF_FILETYPE_NC_RESTART_OLD" $fname $exec_new $exec_old "nc"
#fi
cp -f OPTIONS.nam_new OPTIONS.nam
#
echo "########### SURF_FILETYPE_NC_RESTART "
cp -f TESTS/FORCAGES/HIVER/Params_config.txt_restart Params_config.txt 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH\ =\ 1/NMONTH\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY\ =\ 29/NDAY\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SURF_FILETYPE_NC_RESTART
#./script_exec_restart_parall.sh "SURF_FILETYPE_NC_RESTART" $fname $2 $3


echo "########### SURF_FILETYPE_NC_RESTART_OLD "
rm -f PGD_NEW.nc PREP_NEW.nc
mv -f PGD_OLD.nc PGD.nc
mv -f PREP_OLD.nc PREP.nc
#./script_exec_restart.sh "SURF_FILETYPE_NC_RESTART_OLD" $fname $2 $3


#ASCII
echo "########### SURF_FILETYPE_ASCII "
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
cp -f OPTIONS.nam_base2 OPTIONS.nam_save
sed -e "s/LRESET\_BUDGETC\ =\ T/LRESET\_BUDGETC\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam_restart
#
cp -f OPTIONS.nam_restart OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SURF_FILETYPE_ASCII
#./script_exec_all_parall.sh "SURF_FILETYPE_ASCII" $fname $2 $3
rm -f PGD.txt PREP.txt SURFOUT.txt

#
echo "########### SURF_FILETYPE_ASCII_RESTART "
./script_offline_exte.sh "SURF_FILETYPE_ASCII_RESTART" $fname $exec_new $exec_old "txt"
if [ ! -f PREP_OLD.txt ]
then
  ./script_offline_exte.sh "SURF_FILETYPE_ASCII_RESTART_OLD" $fname $exec_new $exec_old "txt"
fi
cp -f OPTIONS.nam_new OPTIONS.nam
#
cp -f TESTS/FORCAGES/HIVER/Params_config.txt_restart Params_config.txt 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH\ =\ 1/NMONTH\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY\ =\ 29/NDAY\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SURF_FILETYPE_ASCII_RESTART
./script_exec_restart_parall.sh "SURF_FILETYPE_ASCII_RESTART" $fname $2 $3

#test avec de vieux fichiers PREP PGD ascii
echo "########### SURF_FILETYPE_ASCII_RESTART_OLD "
rm -f PGD_NEW.txt PREP_NEW.txt
mv -f PGD_OLD.txt PGD.txt
mv -f PREP_OLD.txt PREP.txt
./script_exec_restart.sh "SURF_FILETYPE_ASCII_RESTART_OLD" $fname $2 $3


#LFI
echo "########## SURF_FILETYPE_LFI "
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
cp -f OPTIONS.nam_restart OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
#
cp -f OPTIONS.nam OPTIONS.nam_SURF_FILETYPE_LFI
./script_exec_all_parall.sh "SURF_FILETYPE_LFI" $fname $2 $3
rm -f PGD.lfi PREP.lfi SURFOUT.lfi
#
./script_offline_exte.sh "SURF_FILETYPE_LFI_RESTART" $fname $exec_new $exec_old "lfi"
if [ ! -f PREP_OLD.lfi ]
then
  ./script_offline_exte.sh "SURF_FILETYPE_LFI_RESTART_OLD" $fname $exec_new $exec_old "lfi"
fi
cp -f OPTIONS.nam_new OPTIONS.nam
#
echo "########### SURF_FILETYPE_LFI_RESTART "
cp -f TESTS/FORCAGES/HIVER/Params_config.txt_restart Params_config.txt 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH\ =\ 1/NMONTH\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY\ =\ 29/NDAY\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SURF_FILETYPE_LFI_RESTART
./script_exec_restart_parall.sh "SURF_FILETYPE_LFI_RESTART" $fname $2 $3
#test avec de vieux fichiers PREP PGD lfi
echo "########### SURF_FILETYPE_LFI_RESTART_OLD "
rm -f PGD_NEW.lfi PREP_NEW.lfi
mv -f PGD_OLD.lfi PGD.lfi
mv -f PREP_OLD.lfi PREP.lfi
./script_exec_restart.sh "SURF_FILETYPE_LFI_RESTART_OLD" $fname $2 $3
#
#LFI avec reset_budget
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
cp -f OPTIONS.nam_restart OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LRESET\_BUDGET\ =\ F/LRESET\_BUDGET\ = T/g" OPTIONS.nam_save > OPTIONS.nam
#
./script_offline_exte.sh "SURF_FILETYPE_LFI_RESTART_RESET" $fname $exec_new $exec_old "lfi" 
if [ ! -f PREP_OLD.lfi ]
then
  ./script_offline_exte.sh "SURF_FILETYPE_LFI_RESTART_RESET_OLD" $fname $exec_new $exec_old "lfi" 
fi
cp -f OPTIONS.nam_new OPTIONS.nam
#
echo "########### SURF_FILETYPE_LFI_RESTART_RESET "
cp -f TESTS/FORCAGES/HIVER/Params_config.txt_restart Params_config.txt 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH\ =\ 1/NMONTH\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY\ =\ 29/NDAY\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SURF_FILETYPE_LFI_RESTART_RESET
./script_exec_restart_parall.sh "SURF_FILETYPE_LFI_RESTART_RESET" $fname $2 $3
#test avec de vieux fichiers PREP PGD lfi
echo "########### SURF_FILETYPE_LFI_RESTART_RESET_OLD "
rm -f PGD_NEW.lfi PREP_NEW.lfi
mv -f PGD_OLD.lfi PGD.lfi
mv -f PREP_OLD.lfi PREP.lfi
./script_exec_restart.sh "SURF_FILETYPE_LFI_RESTART_RESET_OLD" $fname $2 $3


#LFI en désactivant les diags
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
cp -f OPTIONS.nam_restart OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_restart_lfi
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LDIAG\_GRID\ =\ T/LDIAG\_GRID\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LFRAC\ =\ T/LFRAC\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/L2M\_MIN\_ZS\ =\ T/L2M\_MIN\_ZS\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSURF\_BUDGET\ =\ T/LSURF\_BUDGET\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LRAD\_BUDGET\ =\ T/LRAD\_BUDGET\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LCOEF\ =\ T/LCOEF\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSURF\_VARS\ =\ T/LSURF\_VARS\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LPGD\ =\ T/LPGD\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSURF\_EVAP\_BUDGET\ =\ T/LSURF\_EVAP\_BUDGET\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSURF\_MISC\_BUDGET\ =\ T/LSURF\_MISC\_BUDGET\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LPATCH\_BUDGET\ =\ T/LPATCH\_BUDGET\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSURF\_MISC\_DIF\ =\ T/LSURF\_MISC\_DIF\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LWATER\_BUDGET\ =\ T/LWATER\_BUDGET\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LPGD\_FIX\ =\ T/LPGD\_FIX\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LUTCI\ =\ T/LUTCI\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LDIAG\_OCEAN\ =\ T/LDIAG\_OCEAN\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LWATER\_PROFILE\ =\ T/LWATER\_PROFILE\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
#
./script_offline_exte.sh "SURF_FILETYPE_LFI_RESTART_NO_DIAG" $fname $exec_new $exec_old "lfi"
if [ ! -f PREP_OLD.lfi ]
then
  ./script_offline_exte.sh "SURF_FILETYPE_LFI_RESTART_NO_DIAG_OLD" $fname $exec_new $exec_old "lfi"
fi
cp -f OPTIONS.nam_new OPTIONS.nam
#
echo "########### SURF_FILETYPE_LFI_RESTART_NO_DIAG "
cp -f OPTIONS.nam_restart_lfi OPTIONS.nam
cp -f TESTS/FORCAGES/HIVER/Params_config.txt_restart Params_config.txt 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH\ =\ 1/NMONTH\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY\ =\ 29/NDAY\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SURF_FILETYPE_LFI_RESTART_NO_DIAG
./script_exec_restart_parall.sh "SURF_FILETYPE_LFI_RESTART_NO_DIAG" $fname $2 $3
#test avec de vieux fichiers PREP PGD lfi
echo "########### SURF_FILETYPE_LFI_RESTART_NO_DIAG_OLD "
rm -f PGD_NEW.lfi PREP_NEW.lfi
mv -f PGD_OLD.lfi PGD.lfi
mv -f PREP_OLD.lfi PREP.lfi
./script_exec_restart.sh "SURF_FILETYPE_LFI_RESTART_NO_DIAG_OLD" $fname $2 $3

#FA
echo "############ SURF_FILETYPE_FA "
rm -f *.fa
cp -f OPTIONS.nam_base2 OPTIONS.nam_save
sed -e "s/NIMAX = 12/NIMAX = 20/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NJMAX = 8/NJMAX = 20/g" OPTIONS.nam_save > OPTIONS.nam
#not possible to keep cumulated variable in a restart in FA format
cp -f OPTIONS.nam OPTIONS.nam_restart
#cp -f OPTIONS.nam OPTIONS.nam_save
#sed -e "s/LRESET\_BUDGETC\ =\ T/LRESET\_BUDGETC\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam_restart
#
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
cp -f OPTIONS.nam_restart OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ =\ \"NC\"/CSURF\_FILETYPE\ =\ \"FA\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SURF_FILETYPE_FA
./script_exec_all_parall.sh "SURF_FILETYPE_FA" $fname $2 $3
rm -f PGD.fa PREP.fa SURFOUT.fa
#
rm -f *.fa
./script_offline_exte.sh "SURF_FILETYPE_FA_RESTART" $fname $exec_new $exec_old "fa" 
if [ ! -f PREP_OLD.fa ]
then
  ./script_offline_exte.sh "SURF_FILETYPE_FA_RESTART_OLD" $fname $exec_new $exec_old "fa" 
fi
cp -f OPTIONS.nam_new OPTIONS.nam
#
echo "########### SURF_FILETYPE_FA_RESTART "
cp -f TESTS/FORCAGES/HIVER/Params_config.txt_restart Params_config.txt 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH\ =\ 1/NMONTH\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY\ =\ 29/NDAY\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SURF_FILETYPE_FA_RESTART
./script_exec_restart_parall.sh "SURF_FILETYPE_FA_RESTART" $fname $2 $3
#test avec de vieux fichiers PREP PGD fa
echo "########### SURF_FILETYPE_FA_RESTART_OLD "
rm -f PGD_NEW.fa PREP_NEW.fa
mv -f PGD_OLD.fa PGD.fa
mv -f PREP_OLD.fa PREP.fa
./script_exec_restart.sh "SURF_FILETYPE_FA_RESTART_OLD" $fname $2 $3



#CFORCING_FILETYPE
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .

#ASCII + OPTIONS
echo "########### SET_FORC_ZS "
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/LSET\_FORC\_ZS\ =\ F/LSET\_FORC\_ZS\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SET_FORC_ZS
#./script_exec.sh "SET_FORC_ZS" $fname $2 $3
echo "########### LIMIT_QAIR "
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/LLIMIT\_QAIR\ =\ F/LLIMIT\_QAIR\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LIMIT_QAIR
#./script_exec_parall.sh "LIMIT_QAIR" $fname $2 $3
echo "########### NB_READ_FORC1 "
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/NB\_READ\_FORC\ =\ 0/NB\_READ\_FORC\ =\ 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NB_READ_FORC1
#./script_exec_parall.sh "NB_READ_FORC1" $fname $2 $3 
#là il faudrait vérifier qu'on a bien la même chose quand NB_READ_FORC=0, faire un script spécial, par exemple
echo "########### NB_READ_FORC2 "
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/NB\_READ\_FORC\ =\ 0/NB\_READ\_FORC\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NB_READ_FORC2
#./script_exec.sh "NB_READ_FORC2" $fname $2 $3 
#là il faudrait vérifier qu'on a bien la même chose quand NB_READ_FORC=0, faire un script spécial, par exemple


#BINARY + OPTIONS
echo "########### FORC_BIN "
ln -s TESTS/FORCAGES/BINARY/Forc_*.bin .
cp -f TESTS/FORCAGES/BINARY/Params_config.txt .
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/CFORCING\_FILETYPE\ \=\ \"ASCII\"/CFORCING\_FILETYPE\ \=\ \"BINARY\"/g" OPTIONS.nam_save > OPTIONS.nam_bin
cp -f OPTIONS.nam_bin OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FORC_BIN
#./script_exec_parall.sh "FORC_BIN" $fname $2 $3
echo "########### FORC_BIN_NB_READ_FORC1 "
cp -f OPTIONS.nam_bin OPTIONS.nam_save
sed -e "s/NB\_READ\_FORC\ =\ 0/NB\_READ\_FORC\ =\ 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FORC_BIN_NB_READ_FORC1
#./script_exec.sh "FORC_BIN_NB_READ_FORC1" $fname $2 $3 
#là il faudrait vérifier qu'on a bien la même chose quand NB_READ_FORC=0, faire un script spécial, par exemple
echo "########### FORC_BIN_NB_READ_FORC2 "
cp -f OPTIONS.nam_bin OPTIONS.nam_save
sed -e "s/NB\_READ\_FORC\ =\ 0/NB\_READ\_FORC\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FORC_BIN_NB_READ_FORC2
#./script_exec.sh "FORC_BIN_NB_READ_FORC2" $fname $2 $3 
#là il faudrait vérifier qu'on a bien la même chose quand NB_READ_FORC=0, faire un script spécial, par exemple
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
#
rm -f Forc_*.bin
#
#NETCDF + OPTIONS
echo "########### FORC_NETCDF1 "
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/CFORCING\_FILETYPE\ \=\ \"ASCII\"/CFORCING\_FILETYPE\ \=\ \"NETCDF\"/g" OPTIONS.nam_save > OPTIONS.nam_netcdf
cp -f OPTIONS.nam_netcdf OPTIONS.nam_save
sed -e "s/LSET\_FORC\_ZS\ \=\ F/LSET\_FORC\_ZS\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_netcdf
cp -f OPTIONS.nam_netcdf OPTIONS.nam_save
sed -e "s/NDAY\ =\ 29/NDAY\ =\ 25/g" OPTIONS.nam_save > OPTIONS.nam_netcdf
cp -f OPTIONS.nam_netcdf OPTIONS.nam
cp -f TESTS/FORCAGES/NETCDF/FORCING1.nc FORCING.nc
#./script_exec.sh "FORC_NETCDF1" $fname $2 $3
echo "########### FORC_NETCDF_NB_READ_FORC1 "
cp -f OPTIONS.nam_netcdf OPTIONS.nam_save
sed -e "s/NB\_READ\_FORC\ =\ 0/NB\_READ\_FORC\ =\ 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FORC_NETCDF_NB_READ_FORC1
#./script_exec.sh "FORC_NETCDF_NB_READ_FORC1" $fname $2 $3 
#là il faudrait vérifier qu'on a bien la même chose quand NB_READ_FORC=0, faire un script spécial, par exemple
echo "########### FORC_NETCDF_NB_READ_FORC2 "
cp -f OPTIONS.nam_netcdf OPTIONS.nam_save
sed -e "s/NB\_READ\_FORC\ =\ 0/NB\_READ\_FORC\ =\ 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FORC_NETCDF_NB_READ_FORC2
#./script_exec.sh "FORC_NETCDF_NB_READ_FORC2" $fname $2 $3 
#là il faudrait vérifier qu'on a bien la même chose quand NB_READ_FORC=0, faire un script spécial, par exemple
echo "########### FORC_NETCDF_SEC "
cp -f OPTIONS.nam_netcdf OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FORC_NETCDF_SEC
cp -f TESTS/FORCAGES/NETCDF/FORCING_SEC.nc FORCING.nc
#./script_exec.sh "FORC_NETCDF_SEC" $fname $2 $3
echo "########### FORC_NETCDF_MIN "
cp -f OPTIONS.nam_netcdf OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FORC_NETCDF_MIN
cp -f TESTS/FORCAGES/NETCDF/FORCING_MIN.nc FORCING.nc
#./script_exec_parall.sh "FORC_NETCDF_MIN" $fname $2 $3
echo "########### FORC_NETCDF_HOU "
cp -f OPTIONS.nam_netcdf OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FORC_NETCDF_HOU
cp -f TESTS/FORCAGES/NETCDF/FORCING_HOU.nc FORCING.nc
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XWSNOW\_ROOF\ =\ 330\./XWSNOW\_ROOF\ =\ 0\./g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XWSNOW\_ROAD\ =\ 150\./XWSNOW\_ROAD\ =\ 0\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_hour
#./script_exec.sh "FORC_NETCDF_HOU" $fname $2 $3
echo "########### FORC_NETCDF_DAY "
cp -f OPTIONS.nam_hour OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FORC_NETCDF_DAY
cp -f TESTS/FORCAGES/NETCDF/FORCING_DAY.nc FORCING.nc
#./script_exec.sh "FORC_NETCDF_DAY" $fname $2 $3

echo "########## FORC_NETCDF_DELAYEDSTART_DATESTOP"
cp -f OPTIONS.nam_netcdf OPTIONS.nam
cp -f TESTS/FORCAGES/NETCDF/FORCING1.nc FORCING.nc
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY\ =\ 29/NDAY\ =\ 26/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTIME\ =\ 0./XTIME\ =\ 3600./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LDELAYEDSTART_NC = F/LDELAYEDSTART_NC = T/g"  OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDATESTOP = 0,0,0,0/NDATESTOP = 1986,1,27,43200/g"  OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_FORC_NETCDF_DELAYEDSTART_DATESTOP
#./script_exec.sh "FORC_NETCDF_DELAYEDSTART_DATESTOP" $fname $2 $3
rm -f FORCING.nc


#CTIMESERIES_FILETYPE: attention: du coup on ne pourra pas tester l'égalité des sorties NETCDF

echo "########## OUTPUT_NETCDF_NO_WRITE_COORD"
#test NETCDF + LWRITE_COORD=F (sorties 2d)
sed -e "s/LWRITE\_COORD\ \=\ T/LWRITE\_COORD\ \=\ F/g" OPTIONS.nam_base1 > OPTIONS.nam 
cp -rf OPTIONS.nam OPTIONS.nam_OUTPUT_NETCDF_NO_WRITE_COORD
#./script_exec_parall.sh "OUTPUT_NETCDF_NO_WRITE_COORD" $fname $2 $3


echo "########## OUTPUT_ASCII "
sed -e "s/CTIMESERIES\_FILETYPE\ \=\ \"NETCDF\"/CTIMESERIES\_FILETYPE\ \=\ \"ASCII\"/g" OPTIONS.nam_base1 > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_OUTPUT_ASCII
#./script_exec_parall.sh "OUTPUT_ASCII" $fname $2 $3
rm -f SURFOUT*

echo "########## OUTPUT_LFI"
sed -e "s/CTIMESERIES\_FILETYPE\ \=\ \"NETCDF\"/CTIMESERIES\_FILETYPE\ \=\ \"LFI  \"/g" OPTIONS.nam_base1 > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_OUTPUT_LFI
#./script_exec_parall.sh "OUTPUT_LFI" $fname $2 $3
rm -f SURFOUT*

echo "########## OUTPUT_NC"
sed -e "s/CTIMESERIES\_FILETYPE\ \=\ \"NETCDF\"/CTIMESERIES\_FILETYPE\ \=\ \"NC   \"/g" OPTIONS.nam_base1 > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_OUTPUT_NC
#./script_exec_parall.sh "OUTPUT_NC" $fname $2 $3
rm -f SURFOUT*

echo "########## OUTPUT_FA"
rm -f *.fa
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/NIMAX = 12/NIMAX = 20/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NJMAX = 8/NJMAX = 20/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTIMESERIES\_FILETYPE\ \=\ \"NETCDF\"/CTIMESERIES\_FILETYPE\ \=\ \"FA   \"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_OUTPUT_FA
#./script_exec_parall.sh "OUTPUT_FA" $fname $2 $3
echo "########## OUTPUT_FA_DIAG_FA_NOCOMPACT"
sed -e "s/LDIAG\_FA\_NOCOMPACT\ \=\ F/LDIAG\_FA\_NOCOMPACT\ \=\ T/g" OPTIONS.nam_OUTPUT_FA > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_OUTPUT_FA_DIAG_FA_NOCOMPACT
#./script_exec.sh "OUTPUT_FA_DIAG_FA_NOCOMPACT" $fname $2 $3
#ne fonctionne pas avec la nouvelle XRD....
echo "########## OUTPUT_FA_DIAG_FA_OUT_TIMENAME"
sed -e "s/LOUT\_TIMENAME\ \=\ F/LOUT\_TIMENAME\ \=\ T/g" OPTIONS.nam_OUTPUT_FA > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_OUTPUT_FA_OUT_TIMENAME
#./script_exec.sh "OUTPUT_FA_DIAG_FA_OUT_TIMENAME" $fname $2 $3
rm -f SURFOUT*

echo "########## OUTPUT_TEXTE"
sed -e "s/CTIMESERIES\_FILETYPE\ \=\ \"NETCDF\"/CTIMESERIES\_FILETYPE\ \=\ \"TEXTE \"/g" OPTIONS.nam_base1 > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_OUTPUT_TEXTE
#./script_exec.sh "OUTPUT_TEXTE" $fname $2 $3
rm -f *.TXT

echo "########## OUTPUT_BINARY"
sed -e "s/CTIMESERIES\_FILETYPE\ \=\ \"NETCDF\"/CTIMESERIES\_FILETYPE\ \=\ \"BINARY\"/g" OPTIONS.nam_base1 > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_OUTPUT_BINARY
#./script_exec_parall.sh "OUTPUT_BINARY" $fname $2 $3
rm -f *.BIN

echo "########### SHADOWS_SLOPE"
cp -f OPTIONS.nam_base1 OPTIONS.nam_save
sed -e "s/LSHADOWS_SLOPE = F/LSHADOWS_SLOPE = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SHADOWS_SLOPE
#./script_exec.sh "SHADOWS_SLOPE" $fname $2 $3

echo "########### SHADOWS_OTHER"
cp -f OPTIONS.nam_SHADOWS_SLOPE OPTIONS.nam_save
sed -e "s/LSHADOWS_OTHER = F/LSHADOWS_OTHER = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SHADOWS_OTHER
#./script_exec.sh "SHADOWS_OTHER" $fname $2 $3

#NAM_WRITE_SURF_ATM
echo "########### NAM_WRITE_SURF_ATM"
sed -e "s/LNOWRITE\_CANOPY\ \=\ F/LNOWRITE\_CANOPY\ \=\ T/g" OPTIONS.nam_base2 > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LNOWRITE\_TEXFILE\ \=\ F/LNOWRITE\_TEXFILE\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_NAM_WRITE_SURF_ATM
#./script_exec.sh "NAM_WRITE_SURF_ATM" $fname $2 $3

#NAM_WRITE_COVER_TEX
echo "########### NAM_WRITE_COVER_TEX"
sed -e "s/CLANG\ \=\ \"FR\"/CLANG\ \=\ \"EN\"/g" OPTIONS.nam_base2 > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NAM_WRITE_COVER_TEX_EN
#./script_exec.sh "NAM_WRITE_COVER_TEX_EN" $fname $2 $3

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

