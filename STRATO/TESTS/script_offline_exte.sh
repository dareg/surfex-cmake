echo "OFFLINE_EXTE NEW" >> $2
./pgd.exe $3
./prep.exe $3
./offline.exe $3
rm -f PREP.$5
mv -f PGD.$5 PGD_NEW.$5
mv -f SURFOUT.$5 PREP_NEW.$5
echo "OFFLINE_EXTE OLD" >> $2
./script_to_old.sh
./pgd.exe $4
./prep.exe $4
./offline.exe $4
rm -f PREP.$5
mv -f PGD.$5 PGD_OLD.$5
mv -f SURFOUT.$5 PREP_OLD.$5

