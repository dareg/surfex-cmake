cp -f $1 OPTIONS.nam_save

sed -e "s/XWSNOW = 1.0E+20/XWSNOW = 220./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTSNOW = 1.0E+20/XTSNOW = 273.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRSNOW = 1.0E+20/XRSNOW = 260./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XASNOW = 1.0E+20/XASNOW = 0.95/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSG1SNOW = 1.0E+20/XSG1SNOW = 1.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSG2SNOW = 1.0E+20/XSG2SNOW = 2.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHISTSNOW = 1.0E+20/XHISTSNOW = 3./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XAGESNOW = 1.0E+20/XAGESNOW = 10./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
