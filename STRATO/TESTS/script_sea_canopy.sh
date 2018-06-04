OPTIONS_in=$1

cp -f $OPTIONS_in OPTIONS.nam_canopy

#sans canopy
cp -f OPTIONS.nam_canopy OPTIONS.nam
./script_1d_ocean.sh "OPTIONS.nam" "$2" $3 $4 $5

#test canopy 
cp -f OPTIONS.nam_canopy OPTIONS.nam_save
sed -e "s/LSEA\_SBL\ =\ F/LSEA\_SBL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
./script_1d_ocean.sh "OPTIONS.nam" "$2_CANOPY" $3 $4 $5

