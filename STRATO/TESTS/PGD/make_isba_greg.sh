#time
cp -f $1 OPTIONS.nam_save
sed -e "s/NTIME\ \=\ 36/NTIME\ \=\ 12/g" OPTIONS.nam_save > OPTIONS.nam

#vegtypes
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 2 3 5 6 9 11 13 14 15 16 17 18 19
do
	var="XUNIF_VEGTYPE($i)"
	sed -e "s/$var = 1.0E+20/$var = 0./g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 4 7 8 10 12
do
	var1="CFNAM_VEGTYPE($i)"
	var2="XDATA_VEGTYPE$i.txt"
	sed -e "s/$var1 = \"\"/$var1 = \"$var2\"/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#veg
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var="XUNIF_VEG(1,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0./g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var="XUNIF_VEG(4,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.95/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var1="CFNAM_VEG(7,$i)"
	var2="XDATA_VEG${i}.txt"
	sed -e "s/$var1 = \"\"/$var1 = \"$var2\"/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#lai
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var="XUNIF_LAI(1,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0./g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 12
do
	var="XUNIF_LAI(4,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0./g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 6 7 8 9
do
	var="XUNIF_VEG(4,$i)"
	sed -e "s/$var = 1.0E+20/$var = 5.10/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(4,4) = 1.0E+20/XUNIF_LAI(4,4) = 0.82/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(4,5) = 1.0E+20/XUNIF_LAI(4,5) = 4.83/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(4,10) = 1.0E+20/XUNIF_LAI(4,10) = 4.88/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAI(4,11) = 1.0E+20/XUNIF_LAI(4,11) = 1.17/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var1="CFNAM_LAI(7,$i)"
	var2="XDATA_LAI${i}.txt"
	sed -e "s/$var1 = \"\"/$var1 = \"$var2\"/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#z0
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var="XUNIF_Z0(1,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.013/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var="XUNIF_Z0(4,$i)"
	sed -e "s/$var = 1.0E+20/$var = 1.30/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var1="CFNAM_Z0(7,$i)"
	var2="XDATA_Z0${i}.txt"
	sed -e "s/$var1 = \"\"/$var1 = \"$var2\"/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#emis
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var="XUNIF_EMIS(1,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.94/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var="XUNIF_EMIS(4,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.97/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var1="CFNAM_EMIS(7,$i)"
	var2="XDATA_EMIS${i}.txt"
	sed -e "s/$var1 = \"\"/$var1 = \"$var2\"/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#dg
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,1) = 1.0E+20/XUNIF_DG(1,1) = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,2) = 1.0E+20/XUNIF_DG(1,2) = 0.50/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,3) = 1.0E+20/XUNIF_DG(1,3) = 1.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(4,1) = 1.0E+20/XUNIF_DG(4,1) = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(4,2) = 1.0E+20/XUNIF_DG(4,2) = 2.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(4,3) = 1.0E+20/XUNIF_DG(4,3) = 3.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(7,1) = 1.0E+20/XUNIF_DG(7,1) = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_DG(7,2) = \"\"/CFNAM_DG(7,2) = \"XDATA_DG2.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_DG(7,3) = \"\"/CFNAM_DG(7,3) = \"XDATA_DG3.txt\"/g" OPTIONS.nam_save > OPTIONS.nam

#rootfrac
for i in 1 4 7
do
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_ROOTFRAC($i,1)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_ROOTFRAC($i,2)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.9/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_ROOTFRAC($i,3)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 1.0/g" OPTIONS.nam_save > OPTIONS.nam
done

#rsmin
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RSMIN(1) = 1.0E+20/XUNIF_RSMIN(1) = 150./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RSMIN(4) = 1.0E+20/XUNIF_RSMIN(4) = 150./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_RSMIN(7) = \"\"/CFNAM_RSMIN(7) = \"XDATA_RSMIN.txt\"/g" OPTIONS.nam_save > OPTIONS.nam

#gamma
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GAMMA(1) = 1.0E+20/XUNIF_GAMMA(1) = 0.04/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GAMMA(4) = 1.0E+20/XUNIF_GAMMA(4) = 0.04/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GAMMA(7) = 1.0E+20/XUNIF_GAMMA(7) = 0.00/g" OPTIONS.nam_save > OPTIONS.nam

#wrmax_cf
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_WRMAX_CF(1) = 1.0E+20/XUNIF_WRMAX_CF(1) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_WRMAX_CF(4) = 1.0E+20/XUNIF_WRMAX_CF(4) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_WRMAX_CF(7) = 1.0E+20/XUNIF_WRMAX_CF(7) = 0.20/g" OPTIONS.nam_save > OPTIONS.nam

#rgl
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RGL(1) = 1.0E+20/XUNIF_RGL(1) = 30./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RGL(4) = 1.0E+20/XUNIF_RGL(4) = 30./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RGL(7) = 1.0E+20/XUNIF_RGL(7) = 100./g" OPTIONS.nam_save > OPTIONS.nam

#cv
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CV(1) = 1.0E+20/XUNIF_CV(1) = 0.00001/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CV(4) = 1.0E+20/XUNIF_CV(4) = 0.00001/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CV(7) = 1.0E+20/XUNIF_CV(7) = 0.00002/g" OPTIONS.nam_save > OPTIONS.nam

#z0_O_Z0H
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0_O_Z0H(1) = 1.0E+20/XUNIF_Z0_O_Z0H(1) = 10./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0_O_Z0H(4) = 1.0E+20/XUNIF_Z0_O_Z0H(4) = 10./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0_O_Z0H(7) = 1.0E+20/XUNIF_Z0_O_Z0H(7) = 10./g" OPTIONS.nam_save > OPTIONS.nam

#albnir_veg
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,1) = 1.0E+20/XUNIF_ALBNIR_VEG(1,1) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,2) = 1.0E+20/XUNIF_ALBNIR_VEG(1,2) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,3) = 1.0E+20/XUNIF_ALBNIR_VEG(1,3) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,4) = 1.0E+20/XUNIF_ALBNIR_VEG(1,4) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,5) = 1.0E+20/XUNIF_ALBNIR_VEG(1,5) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,6) = 1.0E+20/XUNIF_ALBNIR_VEG(1,6) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,7) = 1.0E+20/XUNIF_ALBNIR_VEG(1,7) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,8) = 1.0E+20/XUNIF_ALBNIR_VEG(1,8) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,9) = 1.0E+20/XUNIF_ALBNIR_VEG(1,9) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,10) = 1.0E+20/XUNIF_ALBNIR_VEG(1,10) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,11) = 1.0E+20/XUNIF_ALBNIR_VEG(1,11) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(1,12) = 1.0E+20/XUNIF_ALBNIR_VEG(1,12) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
#
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,1) = 1.0E+20/XUNIF_ALBNIR_VEG(4,1) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,2) = 1.0E+20/XUNIF_ALBNIR_VEG(4,2) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,3) = 1.0E+20/XUNIF_ALBNIR_VEG(4,3) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,4) = 1.0E+20/XUNIF_ALBNIR_VEG(4,4) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,5) = 1.0E+20/XUNIF_ALBNIR_VEG(4,5) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,6) = 1.0E+20/XUNIF_ALBNIR_VEG(4,6) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,7) = 1.0E+20/XUNIF_ALBNIR_VEG(4,7) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,8) = 1.0E+20/XUNIF_ALBNIR_VEG(4,8) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,9) = 1.0E+20/XUNIF_ALBNIR_VEG(4,9) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,10) = 1.0E+20/XUNIF_ALBNIR_VEG(4,10) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,11) = 1.0E+20/XUNIF_ALBNIR_VEG(4,11) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(4,12) = 1.0E+20/XUNIF_ALBNIR_VEG(4,12) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
#
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,1) = 1.0E+20/XUNIF_ALBNIR_VEG(7,1) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,2) = 1.0E+20/XUNIF_ALBNIR_VEG(7,2) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,3) = 1.0E+20/XUNIF_ALBNIR_VEG(7,3) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,4) = 1.0E+20/XUNIF_ALBNIR_VEG(7,4) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,5) = 1.0E+20/XUNIF_ALBNIR_VEG(7,5) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,6) = 1.0E+20/XUNIF_ALBNIR_VEG(7,6) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,7) = 1.0E+20/XUNIF_ALBNIR_VEG(7,7) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,8) = 1.0E+20/XUNIF_ALBNIR_VEG(7,8) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,9) = 1.0E+20/XUNIF_ALBNIR_VEG(7,9) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,10) = 1.0E+20/XUNIF_ALBNIR_VEG(7,10) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,11) = 1.0E+20/XUNIF_ALBNIR_VEG(7,11) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBNIR_VEG(7,12) = 1.0E+20/XUNIF_ALBNIR_VEG(7,12) = 0.30/g" OPTIONS.nam_save > OPTIONS.nam

#albvis_veg
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,1) = 1.0E+20/XUNIF_ALBVIS_VEG(1,1) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,2) = 1.0E+20/XUNIF_ALBVIS_VEG(1,2) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,3) = 1.0E+20/XUNIF_ALBVIS_VEG(1,3) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,4) = 1.0E+20/XUNIF_ALBVIS_VEG(1,4) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,5) = 1.0E+20/XUNIF_ALBVIS_VEG(1,5) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,6) = 1.0E+20/XUNIF_ALBVIS_VEG(1,6) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,7) = 1.0E+20/XUNIF_ALBVIS_VEG(1,7) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,8) = 1.0E+20/XUNIF_ALBVIS_VEG(1,8) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,9) = 1.0E+20/XUNIF_ALBVIS_VEG(1,9) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,10) = 1.0E+20/XUNIF_ALBVIS_VEG(1,10) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,11) = 1.0E+20/XUNIF_ALBVIS_VEG(1,11) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(1,12) = 1.0E+20/XUNIF_ALBVIS_VEG(1,12) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
#
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,1) = 1.0E+20/XUNIF_ALBVIS_VEG(4,1) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,2) = 1.0E+20/XUNIF_ALBVIS_VEG(4,2) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,3) = 1.0E+20/XUNIF_ALBVIS_VEG(4,3) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,4) = 1.0E+20/XUNIF_ALBVIS_VEG(4,4) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,5) = 1.0E+20/XUNIF_ALBVIS_VEG(4,5) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,6) = 1.0E+20/XUNIF_ALBVIS_VEG(4,6) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,7) = 1.0E+20/XUNIF_ALBVIS_VEG(4,7) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,8) = 1.0E+20/XUNIF_ALBVIS_VEG(4,8) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,9) = 1.0E+20/XUNIF_ALBVIS_VEG(4,9) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,10) = 1.0E+20/XUNIF_ALBVIS_VEG(4,10) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,11) = 1.0E+20/XUNIF_ALBVIS_VEG(4,11) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(4,12) = 1.0E+20/XUNIF_ALBVIS_VEG(4,12) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
#
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,1) = 1.0E+20/XUNIF_ALBVIS_VEG(7,1) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,2) = 1.0E+20/XUNIF_ALBVIS_VEG(7,2) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,3) = 1.0E+20/XUNIF_ALBVIS_VEG(7,3) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,4) = 1.0E+20/XUNIF_ALBVIS_VEG(7,4) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,5) = 1.0E+20/XUNIF_ALBVIS_VEG(7,5) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,6) = 1.0E+20/XUNIF_ALBVIS_VEG(7,6) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,7) = 1.0E+20/XUNIF_ALBVIS_VEG(7,7) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,8) = 1.0E+20/XUNIF_ALBVIS_VEG(7,8) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,9) = 1.0E+20/XUNIF_ALBVIS_VEG(7,9) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,10) = 1.0E+20/XUNIF_ALBVIS_VEG(7,10) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,11) = 1.0E+20/XUNIF_ALBVIS_VEG(7,11) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBVIS_VEG(7,12) = 1.0E+20/XUNIF_ALBVIS_VEG(7,12) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam

#albuv_veg
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,1) = 1.0E+20/XUNIF_ALBUV_VEG(1,1) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,2) = 1.0E+20/XUNIF_ALBUV_VEG(1,2) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,3) = 1.0E+20/XUNIF_ALBUV_VEG(1,3) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,4) = 1.0E+20/XUNIF_ALBUV_VEG(1,4) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,5) = 1.0E+20/XUNIF_ALBUV_VEG(1,5) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,6) = 1.0E+20/XUNIF_ALBUV_VEG(1,6) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,7) = 1.0E+20/XUNIF_ALBUV_VEG(1,7) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,8) = 1.0E+20/XUNIF_ALBUV_VEG(1,8) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,9) = 1.0E+20/XUNIF_ALBUV_VEG(1,9) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,10) = 1.0E+20/XUNIF_ALBUV_VEG(1,10) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,11) = 1.0E+20/XUNIF_ALBUV_VEG(1,11) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(1,12) = 1.0E+20/XUNIF_ALBUV_VEG(1,12) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
#
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,1) = 1.0E+20/XUNIF_ALBUV_VEG(4,1) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,2) = 1.0E+20/XUNIF_ALBUV_VEG(4,2) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,3) = 1.0E+20/XUNIF_ALBUV_VEG(4,3) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,4) = 1.0E+20/XUNIF_ALBUV_VEG(4,4) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,5) = 1.0E+20/XUNIF_ALBUV_VEG(4,5) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,6) = 1.0E+20/XUNIF_ALBUV_VEG(4,6) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,7) = 1.0E+20/XUNIF_ALBUV_VEG(4,7) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,8) = 1.0E+20/XUNIF_ALBUV_VEG(4,8) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,9) = 1.0E+20/XUNIF_ALBUV_VEG(4,9) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,10) = 1.0E+20/XUNIF_ALBUV_VEG(4,10) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,11) = 1.0E+20/XUNIF_ALBUV_VEG(4,11) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALBUV_VEG(4,12) = 1.0E+20/XUNIF_ALBUV_VEG(4,12) = 0.0525/g" OPTIONS.nam_save > OPTIONS.nam
#
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,1) = \"\"/CFNAM_ALBUV_VEG(7,1) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,2) = \"\"/CFNAM_ALBUV_VEG(7,2) = \"XDATA_ALBUV_VEG.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,3) = \"\"/CFNAM_ALBUV_VEG(7,3) = \"XDATA_ALBUV_VEG.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,4) = \"\"/CFNAM_ALBUV_VEG(7,4) = \"XDATA_ALBUV_VEG.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,5) = \"\"/CFNAM_ALBUV_VEG(7,5) = \"XDATA_ALBUV_VEG.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,6) = \"\"/CFNAM_ALBUV_VEG(7,6) = \"XDATA_ALBUV_VEG.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,7) = \"\"/CFNAM_ALBUV_VEG(7,7) = \"XDATA_ALBUV_VEG.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,8) = \"\"/CFNAM_ALBUV_VEG(7,8) = \"XDATA_ALBUV_VEG.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,9) = \"\"/CFNAM_ALBUV_VEG(7,9) = \"XDATA_ALBUV_VEG.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,10) = \"\"/CFNAM_ALBUV_VEG(7,10) = \"XDATA_ALBUV_VEG.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,11) = \"\"/CFNAM_ALBUV_VEG(7,11) = \"XDATA_ALBUV_VEG.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(7,12) = \"\"/CFNAM_ALBUV_VEG(7,12) = \"XDATA_ALBUV_VEG.txt\"/g" OPTIONS.nam_save > OPTIONS.nam

for j in 1 2 3 4 5 6 7 8 9 10 11 12
do
for i in 1 4 7
do
	#albnir_soil
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_ALBNIR_SOIL($i,$j)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.3/g" OPTIONS.nam_save > OPTIONS.nam

	#albvis_soil
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_ALBVIS_SOIL($i,$j)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.1/g" OPTIONS.nam_save > OPTIONS.nam

	#albuv_soil
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_ALBUV_SOIL($i,$j)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.06/g" OPTIONS.nam_save > OPTIONS.nam

done 
done

for i in 1 4 7
do
	#gmes
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_GMES($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.001/g" OPTIONS.nam_save > OPTIONS.nam

	#re25
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_RE25($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.0000003/g" OPTIONS.nam_save > OPTIONS.nam

	#bslai
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_BSLAI($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.25/g" OPTIONS.nam_save > OPTIONS.nam

	#laimin
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_LAIMIN($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 =0.3/g" OPTIONS.nam_save > OPTIONS.nam

	#sefold
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_SEFOLD($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 31536000./g" OPTIONS.nam_save > OPTIONS.nam

	#gc
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_GC($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.00015/g" OPTIONS.nam_save > OPTIONS.nam

	#dmax
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_DMAX($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.1/g" OPTIONS.nam_save > OPTIONS.nam

	#f2i
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_F2I($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.3/g" OPTIONS.nam_save > OPTIONS.nam

	#h_tree
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_H_TREE($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 10./g" OPTIONS.nam_save > OPTIONS.nam

	#ce_nitro
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_CE_NITRO($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 4.83/g" OPTIONS.nam_save > OPTIONS.nam

	#cf_nitro
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_CF_NITRO($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 2.53/g" OPTIONS.nam_save > OPTIONS.nam

	#cna_nitro
	cp -f OPTIONS.nam OPTIONS.nam_save
	var1="XUNIF_CNA_NITRO($i)"
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.0/g" OPTIONS.nam_save > OPTIONS.nam

done

cp -f OPTIONS.nam $1
