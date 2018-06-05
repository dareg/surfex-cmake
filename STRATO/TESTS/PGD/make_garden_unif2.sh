cp -f $1 OPTIONS.nam_save
sed -e "s/NTIME = 0/NTIME = 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYP_GARDEN_HVEG = \"\"/CTYP_GARDEN_HVEG = \"CONI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYP_GARDEN_LVEG = \"\"/CTYP_GARDEN_LVEG = \"C3\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYP_GARDEN_NVEG = \"\"/CTYP_GARDEN_NVEG = \"ROCK\"/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_FRAC_HVEG = 1.0E+20/XUNIF_FRAC_HVEG = 0.3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_FRAC_LVEG = 1.0E+20/XUNIF_FRAC_LVEG = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_FRAC_NVEG = 1.0E+20/XUNIF_FRAC_NVEG = 0.2/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_HVEG = 1.0E+20/XUNIF_H_HVEG = 10./g" OPTIONS.nam_save > OPTIONS.nam

for i in "HVEG" "LVEG"
do
	for j in 5 6 7 8 9 10
	do
		cp -f OPTIONS.nam OPTIONS.nam_save
		var1="XUNIF_LAI_$i($j)"
		sed -e "s/$var1 = 1.0E+20/$var1 = 2./g" OPTIONS.nam_save > OPTIONS.nam
	done
	for j in 1 2 3
	do
		cp -f OPTIONS.nam OPTIONS.nam_save
		var1="XUNIF_LAI_$i($j)"
		sed -e "s/$var1 = 1.0E+20/$var1 = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
	done
	for j in 4 11 12
	do
		cp -f OPTIONS.nam OPTIONS.nam_save
		var1="XUNIF_LAI_$i($j)"
		sed -e "s/$var1 = 1.0E+20/$var1 = 1./g" OPTIONS.nam_save > OPTIONS.nam
	done
done

cp -f OPTIONS.nam $1

