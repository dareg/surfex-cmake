OPTIONS_in=$1

cp -f $OPTIONS_in OPTIONS.nam_pgd

cp -f OPTIONS.nam_pgd OPTIONS.nam
./script_teb_prep.sh "OPTIONS.nam" "$2" $3 $4 $5

echo "############# $2_9TP "
cp -f OPTIONS.nam_pgd OPTIONS.nam_save
sed -e "s/NTEB_PATCH\ =\ 1/NTEB_PATCH\ =\ 9/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CROAD\_DIR =\ \"UNIF\"/CROAD\_DIR\ =\ \"ORIE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_9p
./script_teb_prep.sh "OPTIONS.nam" "$2_9TP" $3 $4 $5


echo "############ $2_9TP_OTHER_LAYERS "
cp -f OPTIONS.nam_9p OPTIONS.nam_save
sed -e "s/NROAD\_LAYER\ \=\ 5/NROAD\_LAYER\ \=\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NROOF\_LAYER\ \=\ 5/NROOF\_LAYER\ \=\ 4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NWALL\_LAYER\ \=\ 5/NWALL\_LAYER\ \=\ 6/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NFLOOR\_LAYER\ \=\ 5/NFLOOR\_LAYER\ \=\ 2/g" OPTIONS.nam_save > OPTIONS.nam
./script_teb_prep.sh "OPTIONS.nam" "$2_9TP_OTHER_LAYERS" $3 $4 $5

