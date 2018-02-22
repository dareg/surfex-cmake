
#lai
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
do
	var="XUNIF_LAI(1,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.1/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,3) = 1.0E+20/XUNIF_LAI(1,3) = 0.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,4) = 1.0E+20/XUNIF_LAI(1,4) = 0.3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,5) = 1.0E+20/XUNIF_LAI(1,5) = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,6) = 1.0E+20/XUNIF_LAI(1,6) = 1.1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,7) = 1.0E+20/XUNIF_LAI(1,7) = 1.8/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,8) = 1.0E+20/XUNIF_LAI(1,8) = 2.6/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,9) = 1.0E+20/XUNIF_LAI(1,9) = 3.4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,10) = 1.0E+20/XUNIF_LAI(1,10) = 3.6/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,11) = 1.0E+20/XUNIF_LAI(1,11) = 3.4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,12) = 1.0E+20/XUNIF_LAI(1,12) = 2.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,13) = 1.0E+20/XUNIF_LAI(1,13) = 1.8/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,14) = 1.0E+20/XUNIF_LAI(1,14) = 1.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,15) = 1.0E+20/XUNIF_LAI(1,15) = 0.9/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,16) = 1.0E+20/XUNIF_LAI(1,16) = 0.6/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,17) = 1.0E+20/XUNIF_LAI(1,17) = 0.4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(1,18) = 1.0E+20/XUNIF_LAI(1,18) = 0.2/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
