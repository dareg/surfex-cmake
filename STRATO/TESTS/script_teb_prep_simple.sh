OPTIONS_in=$1

cp -f $OPTIONS_in OPTIONS.nam_prep

cp -f OPTIONS.nam_prep OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2"
./script_exec_parall.sh "$2" $3 $4 $5

echo "############# $2_TWO_WALLS "
cp -f OPTIONS.nam_prep OPTIONS.nam_save
sed -e "s/CWALL\_OPT\ =\ \"UNIF\"/CWALL\_OPT\ =\ \"TWO\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_TWO_WALLS"
./script_exec.sh "$2_TWO_WALLS" $3 $4 $5

echo "############# $2_CANOPY "
cp -f OPTIONS.nam_prep OPTIONS.nam_save
sed -e "s/LTEB\_CANOPY\ =\ F/LTEB\_CANOPY\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_CANOPY"
./script_exec.sh "$2_CANOPY" $3 $4 $5

echo "############# $2_TWO_WALLS_CANOPY "
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CWALL\_OPT\ =\ \"UNIF\"/CWALL\_OPT\ =\ \"TWO\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_"$2_TWO_WALLS_CANOPY"
./script_exec_parall.sh "$2_TWO_WALLS_CANOPY" $3 $4 $5

