#vegtypes
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 19
do
	var="XUNIF_VEGTYPE($i)"
	sed -e "s/$var = 1.0E+20/$var = 0./g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
sed -e "s/XUNIF_VEGTYPE(7) = 1.0E+20/XUNIF_VEGTYPE(7) = 1./g" OPTIONS.nam_save > OPTIONS.nam


cp -f OPTIONS.nam $1
