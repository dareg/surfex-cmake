OPTIONS_in=$1

cp -f $OPTIONS_in OPTIONS.nam_canopy

#sans canopy
cp -f OPTIONS.nam_canopy OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2"

./script_exec.sh "$2" $3 $4 $5

#test canopy 
echo "############# $2_CANOPY "
cp -f OPTIONS.nam_canopy OPTIONS.nam_save
sed -e "s/LISBA\_CANOPY\ =\ F/LISBA\_CANOPY\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_CANOPY"
./script_exec.sh "$2_CANOPY" $3 $4 $5

#test canopy 
echo "############# $2_CANOPY_DRAG "
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LCANOPY\_DRAG\ =\ F/LCANOPY\_DRAG\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_CANOPY_DRAG"
./script_exec.sh "$2_CANOPY_DRAG" $3 $4 $5

