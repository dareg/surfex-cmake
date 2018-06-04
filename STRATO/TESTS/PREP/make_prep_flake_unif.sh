cp -f $1 OPTIONS.nam_save

sed -e "s/XTS_UNIF = 1.0E+20/XTS_UNIF = 298./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_T_SNOW = 1.0E+20/XUNIF_T_SNOW = 273.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_T_ICE = 1.0E+20/XUNIF_T_ICE = 273.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_T_WML = 1.0E+20/XUNIF_T_WML = 298./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_T_BOT = 1.0E+20/XUNIF_T_BOT = 290./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_T_B1 = 1.0E+20/XUNIF_T_B1 = 290./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

sed -e "s/XUNIF_CT = 1.0E+20/XUNIF_CT = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

sed -e "s/XUNIF_H_SNOW = 1.0E+20/XUNIF_H_SNOW = 0.02/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_ICE = 1.0E+20/XUNIF_H_ICE = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_WML = 1.0E+20/XUNIF_H_WML = 8./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_B1 = 1.0E+20/XUNIF_H_B1 = 0.1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save


cp -f OPTIONS.nam $1
