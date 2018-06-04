echo "PREP_EXTE NEW" >> $2
./pgd.exe $3 >> $2
./prep.exe $3 >> $2
mv -f PGD.$5 PGD_BASE_NEW.$5
mv -f PREP.$5 PREP_BASE_NEW.$5
echo "PREP_EXTE OLD" >> $2
./script_to_old.sh
./pgd.exe $4 >> $2
./prep.exe $4 >> $2
mv -f PGD.$5 PGD_BASE_OLD.$5
mv -f PREP.$5 PREP_BASE_OLD.$5

