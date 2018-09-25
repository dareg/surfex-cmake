OPTIONS_in=$1

cp -f $OPTIONS_in OPTIONS.nam_neige

echo "############# $2_D95 "
cp -f OPTIONS.nam OPTIONS.nam_"$2_D95"
./script_exec.sh "$2_D95" $3 $4 $5

echo "############# $2_EBA "
cp -f OPTIONS.nam_neige OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"EBA\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_EBA"
./script_exec.sh "$2_EBA" $3 $4 $5

echo "############# $2_SN3L "
cp -f OPTIONS.nam_neige OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_sn3l
./script_isba_canopy_parall.sh "OPTIONS.nam" "$2_SN3L" $3 $4 $5
echo "############# $2_SN3L_GLACIER "
mv -f OPTIONS.nam_sn3l OPTIONS.nam_save
sed -e "s/LGLACIER\ =\ F/LGLACIER\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SN3L_GLACIER"
./script_exec.sh "$2_SN3L_GLACIER" $3 $4 $5
echo "############# $2_SN3L_GLACIER_RIL "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOWRES\ =\ \"DEF\"/CSNOWRES\ =\ \"RIL\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SN3L_GLACIER_RIL"
./script_exec.sh "$2_SN3L_GLACIER_RIL" $3 $4 $5

echo "############# $2_SNCRO8"
cp -f OPTIONS.nam_neige OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"CRO\"/g" OPTIONS.nam_save > OPTIONS.nam
#mv -f OPTIONS.nam OPTIONS.nam_save
#sed -e "s/LPROSNOW\ =\ F/LPROSNOW\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 8/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8"
./script_exec.sh "$2_SNCRO8" $3 $4 $5
echo "############# $2_SNCRO8_GLACIER"
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LGLACIER\ =\ F/LGLACIER\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO8_GLACIER"
./script_exec.sh "$2_SNCRO8_GLACIER" $3 $4 $5
echo "############# $2_SNCRO8_GLACIER_RIL"
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOWRES\ =\ \"DEF\"/CSNOWRES\ =\ \"RIL\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ril
./script_isba_canopy_parall.sh "OPTIONS.nam" "$2_SNCRO8_GLACIER_RIL" $3 $4 $5
echo "############# $2_SNCRO8_GLACIER_RIL_SNOWDRIFT_SUBLIM"
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
echo "############# $2_SNCRO10"
cp -f OPTIONS.nam_neige OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"CRO\"/g" OPTIONS.nam_save > OPTIONS.nam
#mv -f OPTIONS.nam OPTIONS.nam_save
#sed -e "s/LPROSNOW\ =\ F/LPROSNOW\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO10"
./script_exec.sh "$2_SNCRO10" $3 $4 $5

#test crocus avec 15 couches
echo "############# $2_SNCRO15"
cp -f OPTIONS.nam_neige OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"CRO\"/g" OPTIONS.nam_save > OPTIONS.nam
#mv -f OPTIONS.nam OPTIONS.nam_save
#sed -e "s/LPROSNOW\ =\ F/LPROSNOW\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ =\ 1/NSNOW\_LAYER\ =\ 15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_SNCRO15"
./script_exec.sh "$2_SNCRO15" $3 $4 $5
