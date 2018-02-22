cp -f OPTIONS.nam OPTIONS.nam_new

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "/COBS_M/d" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_old

