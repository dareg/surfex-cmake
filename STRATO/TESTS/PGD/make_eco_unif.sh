cp -f $1 OPTIONS.nam_save
sed -e "s/XUNIF_COVER(1)\ = 0./XUNIF_COVER(1) = 0.1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_COVER(2)\ = 0./XUNIF_COVER(2) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_COVER(10) = 0./XUNIF_COVER(10) = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_COVER(151) = 0./XUNIF_COVER(151) = 0.15/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1

