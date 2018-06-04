#dg
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,1) = 1.0E+20/XUNIF_DG(1,1) = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,2) = 1.0E+20/XUNIF_DG(1,2) = 0.8/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,3) = 1.0E+20/XUNIF_DG(1,3) = 1.1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,4) = 1.0E+20/XUNIF_DG(1,4) = 1.6/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,5) = 1.0E+20/XUNIF_DG(1,5) = 2.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,6) = 1.0E+20/XUNIF_DG(1,6) = 3./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,7) = 1.0E+20/XUNIF_DG(1,7) = 4.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,8) = 1.0E+20/XUNIF_DG(1,8) = 6./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
