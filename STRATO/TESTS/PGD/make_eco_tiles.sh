cp -f $1 OPTIONS.nam_save
sed -e "s/XUNIF_WATER = 1.0E+20/XUNIF_WATER = 0.1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_SEA = 1.0E+20/XUNIF_SEA = 0.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TOWN = 1.0E+20/XUNIF_TOWN = 0.23/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_NATURE = 1.0E+20/XUNIF_NATURE = 0.52/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1

