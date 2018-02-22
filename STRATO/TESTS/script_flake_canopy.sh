OPTIONS_in=$1

cp -f $OPTIONS_in OPTIONS.nam_canopy

#sans canopy
echo "############# $2"
cp -f OPTIONS.nam_canopy OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2
./script_exec_parall.sh "$2" $3 $4 $5

#test canopy 
echo "############# $2_CANOPY  "
cp -f OPTIONS.nam_canopy OPTIONS.nam_save
sed -e "s/LWAT\_SBL\ =\ F/LWAT\_SBL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_CANOPY
./script_exec.sh "$2_CANOPY" $3 $4 $5


#sans canopy + sediments
echo "############# $2_SEDIMENTS "
cp -f OPTIONS.nam_canopy OPTIONS.nam_save
sed -e "s/LSEDIMENTS\ =\ F/LSEDIMENTS\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_SEDIMENTS
./script_exec.sh "$2_SEDIMENTS" $3 $4 $5

#test canopy + sediments
echo "############# $2_CANOPY_SEDIMENTS "
cp -f OPTIONS.nam_canopy OPTIONS.nam_save
sed -e "s/LWAT\_SBL\ =\ F/LWAT\_SBL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSEDIMENTS\ =\ F/LSEDIMENTS\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_$2_CANOPY_SEDIMENTS
./script_exec_parall.sh "$2_CANOPY_SEDIMENTS" $3 $4 $5


