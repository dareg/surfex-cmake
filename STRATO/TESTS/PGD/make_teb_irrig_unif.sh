cp -f $1 OPTIONS.nam_save

sed -e "s/XUNIF_GD_START_MONTH = 1.0E+20/XUNIF_GD_START_MONTH = 5./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GD_END_MONTH = 1.0E+20/XUNIF_GD_END_MONTH = 8./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GD_START_HOUR = 1.0E+20/XUNIF_GD_START_HOUR = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GD_END_HOUR = 1.0E+20/XUNIF_GD_END_HOUR = 7./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GD_24H_IRRIG = 1.0E+20/XUNIF_GD_24H_IRRIG = 4./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

sed -e "s/XUNIF_GR_START_MONTH = 1.0E+20/XUNIF_GR_START_MONTH = 6./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GR_END_MONTH = 1.0E+20/XUNIF_GR_END_MONTH = 8./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GR_START_HOUR = 1.0E+20/XUNIF_GR_START_HOUR = 3./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GR_END_HOUR = 1.0E+20/XUNIF_GR_END_HOUR = 3./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GR_24H_IRRIG = 1.0E+20/XUNIF_GR_24H_IRRIG = 2./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

sed -e "s/XUNIF_RD_START_MONTH = 1.0E+20/XUNIF_RD_START_MONTH = 8./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RD_END_MONTH = 1.0E+20/XUNIF_RD_END_MONTH = 8./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RD_START_HOUR = 1.0E+20/XUNIF_RD_START_HOUR = 5./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RD_END_HOUR = 1.0E+20/XUNIF_RD_END_HOUR = 7./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RD_24H_IRRIG = 1.0E+20/XUNIF_RD_24H_IRRIG = 1./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1

