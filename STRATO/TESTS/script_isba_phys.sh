OPTIONS_in=$1

cp -f $OPTIONS_in OPTIONS.nam_phys

cp -f OPTIONS.nam_phys OPTIONS.nam
./script_exec.sh "$2" $3 $4 $5

echo "############# $2_DEEPSOIL "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/LDEEPSOIL\ =\ F/LDEEPSOIL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_DEEPSOIL"
./script_exec.sh "$2_DEEPSOIL" $3 $4 $5

echo "############# $2_TEMP_ARP "
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LTEMP_ARP\ =\ F/LTEMP_ARP\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
./script_isba_canopy_parall.sh "OPTIONS.nam" "$2_TEMP_ARP" $3 $4 $5

echo "############# $2_PHYSDOMC "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/LPHYSDOMC\ =\ F/LPHYSDOMC\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
./script_isba_canopy.sh "OPTIONS.nam" "$2_PHYSDOMC" $3 $4 $5


echo "############# $2_GB93 "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CC1DRY\ =\ \"DEF\"/CC1DRY\ =\ \"GB93\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_isba_canopy.sh "OPTIONS.nam" "$2_GB93" $3 $4 $5

echo "############# $2_PL98 "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CSCOND\ =\ \"NP89\"/CSCOND\ =\ \"PL98\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_isba_canopy.sh "OPTIONS.nam" "$2_PL98" $3 $4 $5

echo "############# $2_LWT "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CSOILFRZ\ =\ \"DEF\"/CSOILFRZ\ =\ \"LWT\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_isba_canopy.sh "OPTIONS.nam" "$2_LWT" $3 $4 $5

echo "############# $2_WET "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CALBEDO\ =\ \"MEAN\"/CALBEDO\ =\ \"WET\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_isba_canopy.sh "OPTIONS.nam" "$2_WET" $3 $4 $5
echo "############# $2_EVOL "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CALBEDO\ =\ \"MEAN\"/CALBEDO\ =\ \"EVOL\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_isba_canopy.sh "OPTIONS.nam" "$2_EVOL" $3 $4 $5
echo "############# $2_CM13 "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CALBEDO\ =\ \"MEAN\"/CALBEDO\ =\ \"CM13\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_CM13"
./script_exec.sh "$2_CM13" $3 $4 $5

echo "############# $2_HUM "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CCP\_SURF\ =\ \"DRY\"/CCP\_SURF\ =\ \"HUM\"/g" OPTIONS.nam_save > OPTIONS.nam
./script_isba_canopy.sh "OPTIONS.nam" "$2_HUM" $3 $4 $5


echo "############# $2_DT92 "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CRUNOFF\ =\ \"WSAT\"/CRUNOFF\ =\ \"DT92\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_DT92"
./script_exec.sh "$2_DT92" $3 $4 $5
#echo "############# $2_DT92_DEF "
./script_isba_canopy.sh "OPTIONS.nam" "$2_DT92_DEF" $3 $4 $5
echo "############# $2_RUNOFF_SGH "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CRUNOFF\ =\ \"WSAT\"/CRUNOFF\ =\ \"SGH\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_RUNOFF_SGH"
./script_exec.sh "$2_RUNOFF_SGH" $3 $4 $5

toto=`echo $2|cut -c -8`
if [ "$toto" != "ISBA_DIF" ] ; then
  echo "############# $2_KSAT_SGH "
  cp -f OPTIONS.nam_phys OPTIONS.nam_save
  sed -e "s/CKSAT\ =\ \"DEF\"/CKSAT\ =\ \"SGH\"/g" OPTIONS.nam_save > OPTIONS.nam
  cp -f OPTIONS.nam OPTIONS.nam_"$2_KSAT_SGH"
  ./script_exec.sh "$2_KSAT_SGH" $3 $4 $5
fi

echo "############# $2_RAIN_SGH "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CRAIN\ =\ \"DEF\"/CRAIN\ =\ \"SGH\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_RAIN_SGH"
./script_exec.sh "$2_RAIN_SGH" $3 $4 $5

echo "############# $2_HORT_SGH "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CHORT\ =\ \"DEF\"/CHORT\ =\ \"SGH\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_HORT_SGH"
./script_exec.sh "$2_HORT_SGH" $3 $4 $5

echo "############# $2_LSOC_SGH "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/LSOC\ =\ F/LSOC\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_LSOC_SGH"
./script_exec.sh "$2_LSOC_SGH" $3 $4 $5

echo "############# $2_ALL5_SGH "
if [ "$toto" != "ISBA_DIF" ] ; then
   cp -f OPTIONS.nam_phys OPTIONS.nam_save
   sed -e "s/CKSAT\ =\ \"DEF\"/CKSAT\ =\ \"SGH\"/g" OPTIONS.nam_save > OPTIONS.nam
fi
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSOC\ =\ F/LSOC\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CRUNOFF\ =\ \"WSAT\"/CRUNOFF\ =\ \"SGH\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CRAIN\ =\ \"DEF\"/CRAIN\ =\ \"SGH\"/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CHORT\ =\ \"DEF\"/CHORT\ =\ \"SGH\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_ALL5_SGH"
./script_exec_parall.sh "$2_ALL5_SGH" $3 $4 $5


echo "############# $2_MLCH "
cp -f OPTIONS.nam_phys OPTIONS.nam_save
sed -e "s/CDIFSFCOND\ =\ \"DEF\"/CDIFSFCOND\ =\ \"MLCH\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_MLCH"
./script_isba_canopy.sh "OPTIONS.nam" "$2_MLCH" $3 $4 $5

