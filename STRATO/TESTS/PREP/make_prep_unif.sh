cp -f $1 OPTIONS.nam_save

sed -e "s/NYEAR = 1000000000/NYEAR = 1986/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH = 1000000000/NMONTH = 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY = 1000000000/NDAY = 29/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTIME = 1.0E+20/XTIME = 0./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1

