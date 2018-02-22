cp -f $1 OPTIONS.nam_save

sed -e "s/XWSNOW_ROOF = 1.0E+20/XWSNOW_ROOF = 330./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTSNOW_ROOF = 1.0E+20/XTSNOW_ROOF = 273.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRSNOW_ROOF = 1.0E+20/XRSNOW_ROOF = 600./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XASNOW_ROOF = 1.0E+20/XASNOW_ROOF = 0.95/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

sed -e "s/XWSNOW_ROAD = 1.0E+20/XWSNOW_ROAD = 150./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTSNOW_ROAD = 1.0E+20/XTSNOW_ROAD = 273.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRSNOW_ROAD = 1.0E+20/XRSNOW_ROAD = 550./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XASNOW_ROAD = 1.0E+20/XASNOW_ROAD = 0.95/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

cp -f OPTIONS.nam $1
