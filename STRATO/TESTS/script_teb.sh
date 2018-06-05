dir_run=$1
exec_new=$2
exec_old=$3

fname="INFO_TEB.txt"

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

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_base

cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/LSET\_FORC\_ZS\ \=\ F/LSET\_FORC\_ZS\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_teb

####################CISBA##########################################

#cas initial
echo "########## TEB_SIMPLE"
cp -f OPTIONS.nam_teb OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SIMPLE
./script_teb_pgd.sh "OPTIONS.nam" "TEB_SIMPLE" $fname $2 $3

#on active BEM pour le moment
echo "########## TEB_BEM"
cp -f OPTIONS.nam_teb OPTIONS.nam_save
sed -e "s/CBEM\ =\ \"DEF\"/CBEM\ =\ \"BEM\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_teb_pgd.sh "OPTIONS.nam" "TEB_BEM" $fname $2 $3


#on active toutes les options de TEB: BEM, GARDEN, GREENROOF
cp -f OPTIONS.nam_teb OPTIONS.nam_save
sed -e "s/LGARDEN\ =\ F/LGARDEN\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LGREENROOF\ =\ F/LGREENROOF\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF\_GREENROOF\ =\ 1.0E+20/XUNIF\_GREENROOF\ =\ 0.7/g" OPTIONS.nam_save > OPTIONS.nam
./make_greenroof_unif2.sh "OPTIONS.nam"

#on active des options compliquées de TEB
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CBEM\ =\ \"DEF\"/CBEM\ =\ \"BEM\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NTEB_PATCH\ =\ 1/NTEB_PATCH\ =\ 5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CROAD\_DIR =\ \"UNIF\"/CROAD\_DIR\ =\ \"ORIE\"/g" OPTIONS.nam_save > OPTIONS.nam_step

#on active des options compliquées de ISBA 

#CPHOTO = NIT, NPATCH = 12
cp -f OPTIONS.nam_step OPTIONS.nam_save
sed -e "s/NPATCH\ =\ 1/NPATCH\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CPHOTO\ =\ \"NON\"/CPHOTO\ =\ \"NIT\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_gdgr

#CISBA = 2-L CSNOW = D95
echo "########## TEB_GARDEN_GREENROOF_BEM_NIT"
./script_teb_prep_simple.sh "OPTIONS.nam" "TEB_GARDEN_GREENROOF_BEM_NIT" $fname $2 $3

#CISBA = 3-L CSNOW = 3-L
echo "########## TEB_GARDEN_GREENROOF_BEM_NIT_3L_LSNOWDIMNC_LRESETCUMUL_LVOLUMETRIC_SNOWLIQ"
mv -f OPTIONS.nam_gdgr OPTIONS.nam_save
sed -e "s/CISBA\ =\ \"2-L\"/CISBA\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND\_LAYER\ =\ 2/NGROUND\_LAYER\ =\ 3/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
#mv -f OPTIONS.nam OPTIONS.nam_save
#sed -e "s/LSNOWDIMNC\ =\ F/LSNOWDIMNC\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
#mv -f OPTIONS.nam OPTIONS.nam_save
#sed -e "s/LRESETCUMUL\ =\ F/LRESETCUMUL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
#mv -f OPTIONS.nam OPTIONS.nam_save
#sed -e "s/LVOLUMETRIC_SNOWLIQ\ =\ F/LVOLUMETRIC_SNOWLIQ\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
./script_teb_prep_simple.sh "OPTIONS.nam" "TEB_GARDEN_GREENROOF_BEM_NIT_3L_LSNOWDIMNC_LRESETCUMUL_LVOLUMETRIC_SNOWLIQ" $fname $2 $3

#on ajoute l'irrigation 
rm -f Forc_*.txt
rm -f Params_config.txt
cp -f TESTS/FORCAGES/ETE/Params_config.txt_teb Params_config.txt
ln -s TESTS/FORCAGES/ETE/Forc*.txt .

echo "########## TEB_GARDEN_GREENROOF_BEM_3L_IRRIG"
cp -f OPTIONS.nam_step OPTIONS.nam_save
sed -e "s/NMONTH\ =\ 1/NMONTH\ =\ 7/g" OPTIONS.nam_save > OPTIONS.nam
./make_teb_irrig_unif.sh "OPTIONS.nam"
./script_teb_prep_simple.sh "OPTIONS.nam" "TEB_GARDEN_GREENROOF_BEM_3L_IRRIG" $fname $2 $3

rm -f Forc_*.txt
rm -f Params_config.txt
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
ln -s TESTS/FORCAGES/HIVER/Forc*.txt .

#on ajoute les panneaux solaires
echo "########## TEB_GARDEN_GREENROOF_BEM_3L_SOLAR_PANEL"
cp -f OPTIONS.nam_step OPTIONS.nam_save
sed -e "s/LSOLAR_PANEL\ =\ F/LSOLAR_PANEL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
./script_teb_prep_simple.sh "OPTIONS.nam" "TEB_GARDEN_GREENROOF_BEM_3L_SOLAR_PANEL" $fname $2 $3

echo "########## TEB_GARDEN_GREENROOF_BEM_3L_SOLAR_PANEL_DATA"
cp -f OPTIONS.nam_step OPTIONS.nam_save
sed -e "s/LSOLAR_PANEL\ =\ F/LSOLAR_PANEL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RESIDENTIAL = 1.0E+20/XUNIF_RESIDENTIAL = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_FRAC_PANEL = 1.0E+20/XUNIF_FRAC_PANEL = 0.3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_EMIS_PANEL = 1.0E+20/XUNIF_EMIS_PANEL = 0.93/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALB_PANEL = 1.0E+20/XUNIF_ALB_PANEL = 0.11/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_EFF_PANEL = 1.0E+20/XUNIF_EFF_PANEL = 0.14/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_gdgr
./script_teb_prep_simple.sh "OPTIONS.nam" "TEB_GARDEN_GREENROOF_BEM_3L_SOLAR_PANEL_DATA" $fname $2 $3



#CISBA = DIF CSNOW = 3-L 
echo "########## TEB_GARDEN_GREENROOF_BEM_NIT_DIF"
mv -f OPTIONS.nam_gdgr OPTIONS.nam_save
sed -e "s/CISBA\ =\ \"3-L\"/CISBA\ =\ \"DIF\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND\_LAYER\ =\ 3/NGROUND\_LAYER\ =\ 8/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW_GD\ =\ \"D95\"/CSNOW_GD\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER_GD\ =\ 1/NSNOW\_LAYER_GD\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
./script_teb_prep_simple.sh "OPTIONS.nam" "TEB_GARDEN_GREENROOF_BEM_NIT_DIF" $fname $2 $3

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

