#veg
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 7 8 9 10 11 12 13 14 15 16 17 18
do
	var="XUNIF_VEG(1,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.99/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
for i in 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
do
	var="XUNIF_VEG(1,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEG(1,1) = 1.0E+20/XUNIF_VEG(1,1) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEG(1,2) = 1.0E+20/XUNIF_VEG(1,2) = 0.34/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEG(1,3) = 1.0E+20/XUNIF_VEG(1,3) = 0.52/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEG(1,4) = 1.0E+20/XUNIF_VEG(1,4) = 0.62/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEG(1,5) = 1.0E+20/XUNIF_VEG(1,5) = 0.82/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEG(1,6) = 1.0E+20/XUNIF_VEG(1,6) = 0.92/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
