cp -f $1 OPTIONS.nam_save
sed -e "s/NTIME_GR = 0/NTIME_GR = 1/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYP_GR = \"\"/CTYP_GR = \"SEDUM\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NLAYER_GR = 0/NLAYER_GR = 6/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI_GR(1) = 1.0E+20/XUNIF_LAI_GR(1) = 3./g" OPTIONS.nam_save > OPTIONS.nam

for i in 1 2 3 4
do
	var1="XUNIF_OM_GR($i)"
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.023/g" OPTIONS.nam_save > OPTIONS.nam

	var1="XUNIF_CLAY_GR($i)"
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.15/g" OPTIONS.nam_save > OPTIONS.nam

	var1="XUNIF_SAND_GR($i)"
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.65/g" OPTIONS.nam_save > OPTIONS.nam

done

for i in 5 6
do
	var1="XUNIF_OM_GR($i)"
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/$var1 = 1.0E+20/$var1 = 0./g" OPTIONS.nam_save > OPTIONS.nam

	var1="XUNIF_CLAY_GR($i)"
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/$var1 = 1.0E+20/$var1 = 0./g" OPTIONS.nam_save > OPTIONS.nam

	var1="XUNIF_SAND_GR($i)"
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/$var1 = 1.0E+20/$var1 = 1./g" OPTIONS.nam_save > OPTIONS.nam

done

cp -f OPTIONS.nam $1

