cp -f $1 OPTIONS.nam_save

sed -e "s/XWSNOW_GD = 1.0E+20/XWSNOW_GD = 210./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTSNOW_GD = 1.0E+20/XTSNOW_GD = 273.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRSNOW_GD = 1.0E+20/XRSNOW_GD = 700./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XASNOW_GD = 1.0E+20/XASNOW_GD = 0.95/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

sed -e "s/XWSNOW_GR = 1.0E+20/XWSNOW_GR = 200./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTSNOW_GR = 1.0E+20/XTSNOW_GR = 273.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRSNOW_GR = 1.0E+20/XRSNOW_GR = 650./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XASNOW_GR = 1.0E+20/XASNOW_GR = 0.95/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

cp -f OPTIONS.nam $1
