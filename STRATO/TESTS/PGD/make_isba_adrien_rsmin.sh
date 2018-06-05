#rsmin
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RSMIN(1) = 1.0E+20/XUNIF_RSMIN(1) = 125./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
