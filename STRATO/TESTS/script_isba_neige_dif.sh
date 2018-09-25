OPTIONS_in=$1

cp -f $OPTIONS_in OPTIONS.nam_neige

echo "############# $2_SN3L "
cp -f OPTIONS.nam_neige OPTIONS.nam
./script_isba_canopy_parall.sh "OPTIONS.nam" "$2_SN3L" $3 $4 $5
echo "############# $2_SN3L_GLACIER "
cp -f OPTIONS.nam_neige OPTIONS.nam_save
sed -e "s/LGLACIER\ =\ F/LGLACIER\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SN3L_GLACIER"
./script_exec.sh "$2_SN3L_GLACIER" $3 $4 $5
echo "############# $2_SN3L_GLACIER_RIL "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOWRES\ =\ \"DEF\"/CSNOWRES\ =\ \"RIL\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SN3L_GLACIER_RIL"
./script_exec.sh "$2_SN3L_GLACIER_RIL" $3 $4 $5

echo "############# $2_SN3L_LSNOWDIMNC_LRESETCUMUL_LVOLUMETRIC_SNOWLIQ"
cp -f OPTIONS.nam_neige OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSNOWDIMNC\ =\ F/LSNOWDIMNC\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LRESETCUMUL\ =\ F/LRESETCUMUL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LVOLUMETRIC_SNOWLIQ\ =\ F/LVOLUMETRIC_SNOWLIQ\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SN3L_LSNOWDIMNC_LRESETCUMUL_LVOLUMETRIC_SNOWLIQ"
./script_isba_canopy_parall.sh "OPTIONS.nam" "$2_SN3L_LSNOWDIMNC_LRESETCUMUL_LVOLUMETRIC_SNOWLIQ" $3 $4 $5

echo "############# $2_SNCRO8 "
cp -f OPTIONS.nam_neige OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"3-L\"/CSNOW\ =\ \"CRO\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 3/NSNOW\_LAYER\ =\ 8/g" OPTIONS.nam_save > OPTIONS.nam
#mv -f OPTIONS.nam OPTIONS.nam_save
#sed -e "s/LPROSNOW\ =\ F/LPROSNOW\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8"
./script_exec.sh "$2_SNCRO8" $3 $4 $5
echo "############# $2_SNCRO8_GLACIER "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LGLACIER\ =\ F/LGLACIER\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8_GLACIER"
./script_exec.sh "$2_SNCRO8_GLACIER" $3 $4 $5
echo "############# $2_SNCRO8_GLACIER_RIL "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOWRES\ =\ \"DEF\"/CSNOWRES\ =\ \"RIL\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ril
./script_isba_canopy_parall.sh "OPTIONS.nam" "$2_SNCRO8_GLACIER_RIL" $3 $4 $5
echo "############# $2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM "
mv -f OPTIONS.nam_ril OPTIONS.nam_save
sed -e "s/LSNOWDRIFT\_SUBLIM\ =\ F/LSNOWDRIFT\_SUBLIM\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM"
cp -f OPTIONS.nam OPTIONS.nam_sublim
./script_exec_parall.sh "$2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM" $3 $4 $5
echo "############# $2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM_GA01"
mv -f OPTIONS.nam_sublim OPTIONS.nam_save
sed -e 's/CSNOWDRIFT\ =\ "DFLT"/CSNOWDRIFT\ =\ "GA01"/g' OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM_GA01"
./script_exec.sh "$2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM_GA01" $3 $4 $5
echo "############# $2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM_VI13"
mv -f OPTIONS.nam_sublim OPTIONS.nam_save
sed -e 's/CSNOWDRIFT\ =\ "DFLT"/CSNOWDRIFT\ =\ "VI13"/g' OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM_VI13"
./script_exec.sh "$2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM_VI13" $3 $4 $5
echo "############# $2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM_NONE"
mv -f OPTIONS.nam_sublim OPTIONS.nam_save
sed -e 's/CSNOWDRIFT\ =\ "DFLT"/CSNOWDRIFT\ =\ "NONE"/g' OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM_NONE"
./script_exec.sh "$2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM_NONE" $3 $4 $5

#test crocus avec 10 couches
echo "############# $2_SNCRO10 "
cp -f OPTIONS.nam_neige OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"3-L\"/CSNOW\ =\ \"CRO\"/g" OPTIONS.nam_save > OPTIONS.nam
#mv -f OPTIONS.nam OPTIONS.nam_save
#sed -e "s/LPROSNOW\ =\ F/LPROSNOW\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 3/NSNOW\_LAYER\ =\ 10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO10"
./script_exec.sh "$2_SNCRO10" $3 $4 $5

echo "############# $2_SNCRO10_LSNOWDIMNC_LRESETCUMUL_LVOLUMETRIC_SNOWLIQ"
cp -f OPTIONS.nam_neige OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSNOWDIMNC\ =\ F/LSNOWDIMNC\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LRESETCUMUL\ =\ F/LRESETCUMUL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LVOLUMETRIC_SNOWLIQ\ =\ F/LVOLUMETRIC_SNOWLIQ\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO10_LSNOWDIMNC_LRESETCUMUL_LVOLUMETRIC_SNOWLIQ"
./script_isba_canopy_parall.sh "OPTIONS.nam" "$2_SNCRO10_LSNOWDIMNC_LRESETCUMUL_LVOLUMETRIC_SNOWLIQ" $3 $4 $5

#test crocus avec 15 couches
echo "############# $2_SNCRO15 "
cp -f OPTIONS.nam_neige OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"3-L\"/CSNOW\ =\ \"CRO\"/g" OPTIONS.nam_save > OPTIONS.nam
#mv -f OPTIONS.nam OPTIONS.nam_save
#sed -e "s/LPROSNOW\ =\ F/LPROSNOW\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 3/NSNOW\_LAYER\ =\ 15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO15"
./script_exec.sh "$2_SNCRO15" $3 $4 $5

#CSNOWMETAMO 
ln -s TESTS/ISBA/drdt_bst_fit_60.nc .
cp -f OPTIONS.nam_neige OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"3-L\"/CSNOW\ =\ \"CRO\"/g" OPTIONS.nam_save > OPTIONS.nam
#mv -f OPTIONS.nam OPTIONS.nam_save
#sed -e "s/LPROSNOW\ =\ F/LPROSNOW\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 3/NSNOW\_LAYER\ =\ 8/g" OPTIONS.nam_save > OPTIONS.nam_metamo
#C13
echo "############# $2_SNCRO8_SNOWMETAMO_C13_B92"
cp -f OPTIONS.nam_metamo OPTIONS.nam_save
sed -e "s/CSNOWMETAMO = \"B92\"/CSNOWMETAMO = \"C13\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8_SNOWMETAMO_C13_B92"
./script_exec_parall.sh "$2_SNCRO8_SNOWMETAMO_C13_B92" $3 $4 $5
#T07
echo "############# $2_SNCRO8_SNOWMETAMO_T07_B92"
cp -f OPTIONS.nam_metamo OPTIONS.nam_save
sed -e "s/CSNOWMETAMO = \"B92\"/CSNOWMETAMO = \"T07\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8_SNOWMETAMO_T07_B92"
./script_exec_parall.sh "$2_SNCRO8_SNOWMETAMO_T07_B92" $3 $4 $5
#F06
echo "############# $2_SNCRO8_SNOWMETAMO_F06_B92"
cp -f OPTIONS.nam_metamo OPTIONS.nam_save
sed -e "s/CSNOWMETAMO = \"B92\"/CSNOWMETAMO = \"F06\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8_SNOWMETAMO_F06_B92"
./script_exec_parall.sh "$2_SNCRO8_SNOWMETAMO_F06_B92" $3 $4 $5

#SNOWRAD
#TAR
echo "############# $2_SNCRO8_B92_SNOWRAD_B93"
cp -f OPTIONS.nam_metamo OPTIONS.nam_save
#T17 for V9 vs TAR for V8
sed -e "s/CSNOWRAD = \"B92\"/CSNOWRAD = \"B93\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8_B92_SNOWRAD_B93"
./script_exec_parall.sh "$2_SNCRO8_B92_SNOWRAD_B93" $3 $4 $5
#
#TAR
echo "############# $2_SNCRO8_C13_SNOWRAD_T17"
cp -f OPTIONS.nam_metamo OPTIONS.nam_save
sed -e "s/CSNOWMETAMO = \"B92\"/CSNOWMETAMO = \"C13\"/g" OPTIONS.nam_save > OPTIONS.nam_c13
cp -f OPTIONS.nam_c13 OPTIONS.nam_save
sed -e "s/CSNOWRAD = \"B92\"/CSNOWRAD = \"T17\"/g" OPTIONS.nam_save > OPTIONS.nam
#T17 for V9 vs TAR for V8
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8_C13_SNOWRAD_T17"
./script_exec_parall.sh "$2_SNCRO8_C13_SNOWRAD_T17" $3 $4 $5
#
rm -f drdt_bst_fit_60.nc
