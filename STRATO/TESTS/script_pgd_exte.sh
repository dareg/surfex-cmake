echo "PGD_EXTE NEW" >> $2
./pgd.exe $3
mv -f PGD.$5 PGD_BASE_NEW.$5
echo "PGD_EXTE OLD" >> $2
./script_to_old.sh
./pgd.exe $4
mv -f PGD.$5 PGD_BASE_OLD.$5

