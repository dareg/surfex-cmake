OPTIONS_in=$1

cp -f $OPTIONS_in OPTIONS.nam_run

cp -f OPTIONS.nam_run OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2"
./script_exec_parall.sh "$2" $3 $4 $5

echo "############# $2_BRUT82 "
cp -f OPTIONS.nam_run OPTIONS.nam_save
sed -e "s/CZ0H\ =\ \"KAND07\"/CZ0H\ =\ \"BRUT82\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_BRUT82"
./script_exec.sh "$2_BRUT82" $3 $4 $5
echo "############# $2_MASC95 "
cp -f OPTIONS.nam_run OPTIONS.nam_save
sed -e "s/CZ0H\ =\ \"KAND07\"/CZ0H\ =\ \"MASC95\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_MASC95"
./script_exec.sh "$2_MASC95" $3 $4 $5

echo "############# $2_ROW30 "
cp -f OPTIONS.nam_run OPTIONS.nam_save
sed -e "s/CCH\_BEM\ =\ \"DOE-2\"/CCH\_BEM\ =\ \"ROW30\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_ROW30"
./script_exec.sh "$2_ROW30" $3 $4 $5

echo "############ $2_DXCOIL "
cp -f OPTIONS.nam_run OPTIONS.nam_save
sed -e "s/CCOOL\_COIL\ \=\ \"IDEAL\"/CCOOL\_COIL\ \=\ \"DXCOIL\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_DXCOIL"
./script_exec.sh "$2_DXCOIL" $3 $4 $5

echo "############ $2_FINCAP "
cp -f OPTIONS.nam_run OPTIONS.nam_save
sed -e "s/CHEAT\_COIL\ \=\ \"IDEAL\"/CHEAT\_COIL\ \=\ \"FINCAP\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_FINCAP"
./script_exec.sh "$2_FINCAP" $3 $4 $5

echo "############ $2_AUTOSIZE "
cp -f OPTIONS.nam_run OPTIONS.nam_save
sed -e "s/LAUTOSIZE\ \=\ F/LAUTOSIZE\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_AUTOSIZE"
./script_exec.sh "$2_AUTOSIZE" $3 $4 $5

echo "############ $2_DXCOIL_FINCAP_AUTOSIZE "
cp -f OPTIONS.nam_run OPTIONS.nam_save
sed -e "s/CCOOL\_COIL\ \=\ \"IDEAL\"/CCOOL\_COIL\ \=\ \"DXCOIL\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CHEAT\_COIL\ \=\ \"IDEAL\"/CHEAT\_COIL\ \=\ \"FINCAP\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LAUTOSIZE\ \=\ F/LAUTOSIZE\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_DXCOIL_FINCAP_AUTOSIZE"
./script_exec_parall.sh "$2_DXCOIL_FINCAP_AUTOSIZE" $3 $4 $5

echo "########### $2_DT_RES_OFF "
cp -f OPTIONS.nam_run OPTIONS.nam_save
sed -e "s/XDT_RES = 0./XDT_RES = 3./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDT_OFF = 0./XDT_OFF = 3./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_DT_RES_OFF"
./script_exec_parall.sh "$2_DT_RES_OFF" $3 $4 $5
