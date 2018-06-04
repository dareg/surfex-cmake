OPTIONS_in=$1

cp -f $OPTIONS_in OPTIONS.nam_ocean

echo "########### $2"
#sans mercator
cp -f OPTIONS.nam_ocean OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2
./script_exec.sh "$2" $3 $4 $5

#test LMERCATOR=T
echo "########### $2_MERCATOR"
cp -f OPTIONS.nam_ocean OPTIONS.nam_save
sed -e "s/LOCEAN\_MERCATOR\ =\ F/LOCEAN\_MERCATOR\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_mercator
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR
./script_exec.sh "$2_MERCATOR" $3 $4 $5
#test current =T
echo "########### $2_MERCATOR_CURRENT"
cp -f OPTIONS.nam_mercator OPTIONS.nam_save
sed -e "s/LOCEAN\_CURRENT\ =\ F/LOCEAN\_CURRENT\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_current
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT
./script_exec.sh "$2_MERCATOR_CURRENT" $3 $4 $5
#test avec zero_flux
echo "########### $2_MERCATOR_CURRENT_ZERO_FLUX"
cp -f OPTIONS.nam_current OPTIONS.nam_save
sed -e "s/LZERO\_FLUX\ =\ F/LZERO\_FLUX\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT_ZERO_FLUX
./script_exec.sh "$2_MERCATOR_CURRENT_ZERO_FLUX" $3 $4 $5
#test avec progsst
echo "########### $2_MERCATOR_CURRENT_PROGSST"
cp -f OPTIONS.nam_current OPTIONS.nam_save
sed -e "s/LPROGSST\ =\ F/LPROGSST\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT_PROGSST
./script_exec.sh "$2_MERCATOR_CURRENT_PROGSST" $3 $4 $5
#test avec diapyc
echo "########### $2_MERCATOR_CURRENT_DIAPYC"
cp -f OPTIONS.nam_current OPTIONS.nam_save
sed -e "s/LDIAPYC\ =\ F/LDIAPYC\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT_DIAPYC
./script_exec.sh "$2_MERCATOR_CURRENT_DIAPYC" $3 $4 $5
#test avec diapyc & progsst
echo "########### $2_MERCATOR_CURRENT_DIAPYC_PROGSST"
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LPROGSST\ =\ F/LPROGSST\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT_DIAPYC_PROGSST
./script_exec_parall.sh "$2_MERCATOR_CURRENT_DIAPYC_PROGSST" $3 $4 $5


########## test avec la relaxation


#####current
echo "########### $2_MERCATOR_CURRENT_CUR_REL_TS_REL"
cp -f OPTIONS.nam_current OPTIONS.nam_save
sed -e "s/LCUR\_REL\ =\ F/LCUR\_REL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LTS\_REL\ =\ F/LTS\_REL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT_CUR_REL_TS_REL
./script_exec.sh "$2_MERCATOR_CURRENT_CUR_REL_TS_REL" $3 $4 $5
#test avec changement du time_rel
echo "########### $2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000"
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTIME\_REL\ =\ 25920000\./XTIME\_REL\ =\ 30000000\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000
./script_exec.sh "$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000" $3 $4 $5
cp -f OPTIONS.nam OPTIONS.nam_current_relax
#test avec zero_flux
echo "########### $2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_ZERO_FLUX"
cp -f OPTIONS.nam_current_relax OPTIONS.nam_save
sed -e "s/LZERO\_FLUX\ =\ F/LZERO\_FLUX\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_ZERO_FLUX
./script_exec.sh "$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_ZERO_FLUX" $3 $4 $5
#test avec corr_flux
echo "########### $2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_CORR_FLX"
cp -f OPTIONS.nam_current_relax OPTIONS.nam_save
sed -e "s/LCORR\_FLUX\ =\ F/LCORR\_FLUX\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XCORFLX\ =\ 0\./XCORFLX\ =\ 10\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_CORR_FLX
./script_exec.sh "$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_CORR_FLX" $3 $4 $5
#test avec progsst
echo "########### $2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_PROGSST"
cp -f OPTIONS.nam_current_relax OPTIONS.nam_save
sed -e "s/LPROGSST\ =\ F/LPROGSST\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_PROGSST
./script_exec.sh "$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_PROGSST" $3 $4 $5
#test avec diapyc
echo "########### $2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_DIAPYC"
cp -f OPTIONS.nam_current_relax OPTIONS.nam_save
sed -e "s/LDIAPYC\ =\ F/LDIAPYC\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_DIAPYC
./script_exec.sh "$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_DIAPYC" $3 $4 $5
#test avec diapyc & progsst
echo "########### $2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_DIAPYC_PROGSST"
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LPROGSST\ =\ F/LPROGSST\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_DIAPYC_PROGSST
./script_exec_parall.sh "$2_MERCATOR_CURRENT_CUR_REL_TS_REL_TIME_REL30000000_DIAPYC_PROGSST" $3 $4 $5




