OPTIONS_in=$1

cp -f $OPTIONS_in OPTIONS.nam_photo

cp -f OPTIONS.nam_photo OPTIONS.nam
./script_isba_neige.sh "OPTIONS.nam" "$2" $3 $4 $5

echo "############# $2_12P "
cp -f OPTIONS.nam_photo OPTIONS.nam_save
sed -e "s/NPATCH\ =\ 1/NPATCH\ =\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_12p
./script_isba_canopy_parall.sh "OPTIONS.nam" "$2_12P" $3 $4 $5


echo "############# $2_AST "
cp -f OPTIONS.nam_12p OPTIONS.nam_save
sed -e "s/CPHOTO\ =\ \"NON\"/CPHOTO\ =\ \"AST\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_AST"
./script_exec.sh "$2_AST" $3 $4 $5
echo "############# $2_AST_TR_ML "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LTR\_ML\ =\ F/LTR\_ML\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_ast
./script_isba_canopy_parall.sh "OPTIONS.nam" "$2_AST_TR_ML" $3 $4 $5
echo "############# $2_AST_TR_ML_PRM "
mv -f OPTIONS.nam_ast OPTIONS.nam_save
sed -e "s/CRESPSL\ =\ \"DEF\"/CRESPSL\ =\ \"PRM\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_AST_TR_ML_PRM"
./script_exec.sh "$2_AST_TR_ML_PRM" $3 $4 $5


#on fait un test avec l'irrigation réellement utilisée
#on se place en été
echo "############# $2_AST_TR_ML_PRM_AGRIP_IRRIG_FILE "
rm -f Forc_*.txt
rm -f Params_config.txt
cp -f TESTS/FORCAGES/ETE/Params_config.txt .
ln -s TESTS/FORCAGES/ETE/Forc*.txt .
ln -sf TESTS/ISBA/ECOCLIMAP_II_EUROP_IRRIG.prn .
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YIRRIG = \"\"/YIRRIG = \"ECOCLIMAP_II_EUROP_IRRIG.prn\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XWSNOW\_ROOF\ =\ 330\./XWSNOW\_ROOF\ =\ 0\./g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XWSNOW\_ROAD\ =\ 150\./XWSNOW\_ROAD\ =\ 0\./g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NYEAR\ =\ 1986/NYEAR\ =\ 2003/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH\ =\ 1/NMONTH\ =\ 8/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY\ =\ 29/NDAY\ =\ 10/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTIME\ =\ 0\./XTIME\ =\ 48000\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_AST_TR_ML_PRM_AGRIP_IRRIG_FILE"
cp -f OPTIONS.nam OPTIONS.nam_file
./script_exec_parall.sh "$2_AST_TR_ML_PRM_AGRIP_IRRIG_FILE" $3 $4 $5


#on fait un test avec le LAI non climatologique 
echo "############# $2_AST_TR_ML_PRM_AGRIP_IRRIG_FILE_LCLIM_LAIF "
mv -f OPTIONS.nam_file OPTIONS.nam_save
sed -e "s/LCLIM\_LAI\ =\ T/LCLIM\_LAI\ =\ F/g" OPTIONS.nam_save > OPTIONS.nam
./script_exec.sh "$2_AST_TR_ML_PRM_AGRIP_IRRIG_FILE_LCLIM_LAIF" $3 $4 $5

rm -f ECOCLIMAP_II_EUROP_IRRIG.prn

rm -f Forc_*.txt
rm -f Params_config.txt
cp -f TESTS/FORCAGES/HIVER/Params_config.txt .
ln -s TESTS/FORCAGES/HIVER/Forc*.txt .


echo "############# $2_NIT "
cp -f OPTIONS.nam_12p OPTIONS.nam_save
sed -e "s/CPHOTO\ =\ \"NON\"/CPHOTO\ =\ \"NIT\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_nit
./script_isba_neige.sh "OPTIONS.nam" "$2_NIT" $3 $4 $5
echo "############# $2_NIT_NITRO_DILU "
cp -f OPTIONS.nam_nit OPTIONS.nam_save
sed -e "s/LNITRO_DILU\ =\ \F/LNITRO_DILU\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NIT_NITRO_DILU"
./script_exec.sh "$2_NIT_NITRO_DILU" $3 $4 $5
echo "############# $2_NIT_NITRO_DILU_PRM "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CRESPSL\ =\ \"DEF\"/CRESPSL\ =\ \"PRM\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NIT_NITRO_DILU_PRM"
./script_exec.sh "$2_NIT_NITRO_DILU_PRM" $3 $4 $5
echo "############# $2_NIT_NITRO_DILU_PRM_AGRIP "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LAGRIP\ =\ F/LAGRIP\ =\ \T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NIT_NITRO_DILU_PRM_AGRIP"
./script_exec.sh "$2_NIT_NITRO_DILU_PRM_AGRIP" $3 $4 $5
echo "############# $2_NIT_NITRO_DILU_PRM_AGRIP_TR_ML "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LTR\_ML\ =\ F/LTR\_ML\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NIT_NITRO_DILU_PRM_AGRIP_TR_ML"
./script_exec.sh "$2_NIT_NITRO_DILU_PRM_AGRIP_TR_ML" $3 $4 $5



echo "############# $2_NCB "
cp -f OPTIONS.nam_12p OPTIONS.nam_save
sed -e "s/CPHOTO\ =\ \"NON\"/CPHOTO\ =\ \"NCB\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NCB"
./script_exec.sh "$2_NCB" $3 $4 $5
cp -f OPTIONS.nam OPTIONS.nam_ncb
#crespsl = PRM
echo "############# $2_NCB_PRM "
cp -f OPTIONS.nam_ncb OPTIONS.nam_save
sed -e "s/CRESPSL\ =\ \"DEF\"/CRESPSL\ =\ \"PRM\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NCB_PRM"
./script_exec.sh "$2_NCB_PRM" $3 $4 $5
echo "############# $2_NCB_PRM_TR_ML "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LTR\_ML\ =\ F/LTR\_ML\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NCB_PRM_TR_ML"
./script_exec.sh "$2_NCB_PRM_TR_ML" $3 $4 $5
#crepsl = CNT
echo "############# $2_NCB_CNT "
cp -f OPTIONS.nam_ncb OPTIONS.nam_save
sed -e "s/CRESPSL\ =\ \"DEF\"/CRESPSL\ =\ \"CNT\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NCB_CNT"
./script_exec.sh "$2_NCB_CNT" $3 $4 $5
echo "############# $2_NCB_CNT_TR_ML "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LTR\_ML\ =\ F/LTR\_ML\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NCB_CNT_TR_ML"
./script_exec.sh "$2_NCB_CNT_TR_ML" $3 $4 $5
echo "############# $2_NCB_CNT_TR_ML_NITRO_DILU "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LNITRO_DILU\ =\ F/LNITRO_DILU\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NCB_CNT_TR_ML_NITRO_DILU"
./script_exec.sh "$2_NCB_CNT_TR_ML_NITRO_DILU" $3 $4 $5
echo "############# $2_NCB_CNT_TR_ML_NITRO_DILU_SPINUPCARBS "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSPINUPCARBS\ =\ F/LSPINUPCARBS\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XCO2_START\ =\ 1\.0E+20/XCO2_START\ =\ 290./g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XCO2_END\ =\ 1\.0E+20/XCO2_END\ =\ 310./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NCB_CNT_TR_ML_NITRO_DILU_SPINUPCARBS"
./script_exec.sh "$2_NCB_CNT_TR_ML_NITRO_DILU_SPINUPCARBS" $3 $4 $5
#lagrip = t
echo "############# $2_NCB_AGRIP "
cp -f OPTIONS.nam_ncb OPTIONS.nam_save
sed -e "s/LAGRIP\ =\ F/LAGRIP\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NCB_AGRIP"
./script_exec.sh "$2_NCB_AGRIP" $3 $4 $5
echo "############# $2_NCB_AGRIP_CNT "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CRESPSL\ =\ \"DEF\"/CRESPSL\ =\ \"CNT\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NCB_AGRIP_CNT"
./script_exec.sh "$2_NCB_AGRIP_CNT" $3 $4 $5
echo "############# $2_NCB_AGRIP_CNT_SPINUPCARBW "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSPINUPCARBW\ =\ F/LSPINUPCARBW\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NNBYEARSPINS\ =\ 0/NNBYEARSPINS\ =\ 200/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NNBYEARSPINW\ =\ 0/NNBYEARSPINW\ =\ 100/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_NCB_AGRIP_CNT_SPINUPCARBW"
./script_exec.sh "$2_NCB_AGRIP_CNT_SPINUPCARBW" $3 $4 $5
