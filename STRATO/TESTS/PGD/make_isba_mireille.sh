#time
cp -f $1 OPTIONS.nam_save
sed -e "s/NTIME\ \=\ 36/NTIME\ \=\ 12/g" OPTIONS.nam_save > OPTIONS.nam

#vegtypes
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
do
	var1="CFNAM_VEGTYPE($i)"
	var2="VEGTYPE_P$i.dat"
	sed -e "s/$var1 = \"\"/$var1 = \"$var2\"/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#veg
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var1="CFNAM_VEG(1,$i)"
	sed -e "s/$var1 = \"\"/$var1 = \"VEG_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#lai
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var1="CFNAM_LAI(4,$i)"
	sed -e "s/$var1 = \"\"/$var1 = \"LAI_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#z0
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var1="CFNAM_Z0(1,$i)"
	sed -e "s/$var1 = \"\"/$var1 = \"Z0VEG_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#emis
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var1="CFNAM_EMIS(1,$i)"
	sed -e "s/$var1 = \"\"/$var1 = \"EMIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#dg
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3
do
	var1="CFNAM_DG(1,$i)"
	var2="DG${i}_M7_P0.dat"
	sed -e "s/$var1 = \"\"/$var1 = \"$var2\"/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#rsmin
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_RSMIN(1) = \"\"/CFNAM_RSMIN(1) = \"RSMIN_M4_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam

#gamma
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_GAMMA(1) = \"\"/CFNAM_GAMMA(1) = \"GAMMA_M4_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam

#wrmax_cf
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_WRMAX_CF(1) = \"\"/CFNAM_WRMAX_CF(1) = \"WRMAX_CF_M4_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam

#rgl
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_RGL(1) = \"\"/CFNAM_RGL(1) = \"RGL_M4_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam

#cv
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_CV(1) = \"\"/CFNAM_CV(1) = \"CV_M4_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam

#z0_O_Z0H
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0_O_Z0H(1) = 1.0E+20/XUNIF_Z0_O_Z0H(1) = 10./g" OPTIONS.nam_save > OPTIONS.nam

#albnir_veg
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,1) = \"\"/CFNAM_ALBNIR_VEG(4,1) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,2) = \"\"/CFNAM_ALBNIR_VEG(4,2) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,3) = \"\"/CFNAM_ALBNIR_VEG(4,3) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,4) = \"\"/CFNAM_ALBNIR_VEG(4,4) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,5) = \"\"/CFNAM_ALBNIR_VEG(4,5) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,6) = \"\"/CFNAM_ALBNIR_VEG(4,6) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,7) = \"\"/CFNAM_ALBNIR_VEG(4,7) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,8) = \"\"/CFNAM_ALBNIR_VEG(4,8) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,9) = \"\"/CFNAM_ALBNIR_VEG(4,9) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,10) = \"\"/CFNAM_ALBNIR_VEG(4,10) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,11) = \"\"/CFNAM_ALBNIR_VEG(4,11) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_VEG(4,12) = \"\"/CFNAM_ALBNIR_VEG(4,12) = \"ALBNIR_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam

#albvis_veg
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,1) = \"\"/CFNAM_ALBVIS_VEG(4,1) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,2) = \"\"/CFNAM_ALBVIS_VEG(4,2) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,3) = \"\"/CFNAM_ALBVIS_VEG(4,3) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,4) = \"\"/CFNAM_ALBVIS_VEG(4,4) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,5) = \"\"/CFNAM_ALBVIS_VEG(4,5) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,6) = \"\"/CFNAM_ALBVIS_VEG(4,6) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,7) = \"\"/CFNAM_ALBVIS_VEG(4,7) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,8) = \"\"/CFNAM_ALBVIS_VEG(4,8) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,9) = \"\"/CFNAM_ALBVIS_VEG(4,9) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,10) = \"\"/CFNAM_ALBVIS_VEG(4,10) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,11) = \"\"/CFNAM_ALBVIS_VEG(4,11) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_VEG(4,12) = \"\"/CFNAM_ALBVIS_VEG(4,12) = \"ALBVIS_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam

#albuv_veg
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,1) = \"\"/CFNAM_ALBUV_VEG(1,1) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,2) = \"\"/CFNAM_ALBUV_VEG(1,2) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,3) = \"\"/CFNAM_ALBUV_VEG(1,3) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,4) = \"\"/CFNAM_ALBUV_VEG(1,4) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,5) = \"\"/CFNAM_ALBUV_VEG(1,5) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,6) = \"\"/CFNAM_ALBUV_VEG(1,6) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,7) = \"\"/CFNAM_ALBUV_VEG(1,7) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,8) = \"\"/CFNAM_ALBUV_VEG(1,8) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,9) = \"\"/CFNAM_ALBUV_VEG(1,9) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,10) = \"\"/CFNAM_ALBUV_VEG(1,10) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,11) = \"\"/CFNAM_ALBUV_VEG(1,11) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_VEG(1,12) = \"\"/CFNAM_ALBUV_VEG(1,12) = \"ALBUV_ISBA_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam

#albnir_soil
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,1) = \"\"/CFNAM_ALBNIR_SOIL(1,1) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,2) = \"\"/CFNAM_ALBNIR_SOIL(1,2) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,3) = \"\"/CFNAM_ALBNIR_SOIL(1,3) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,4) = \"\"/CFNAM_ALBNIR_SOIL(1,4) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,5) = \"\"/CFNAM_ALBNIR_SOIL(1,5) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,6) = \"\"/CFNAM_ALBNIR_SOIL(1,6) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,7) = \"\"/CFNAM_ALBNIR_SOIL(1,7) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,8) = \"\"/CFNAM_ALBNIR_SOIL(1,8) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,9) = \"\"/CFNAM_ALBNIR_SOIL(1,9) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,10) = \"\"/CFNAM_ALBNIR_SOIL(1,10) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,11) = \"\"/CFNAM_ALBNIR_SOIL(1,11) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBNIR_SOIL(1,12) = \"\"/CFNAM_ALBNIR_SOIL(1,12) = \"ALBNIR_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam

#albvis_soil
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,1) = \"\"/CFNAM_ALBVIS_SOIL(1,1) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,2) = \"\"/CFNAM_ALBVIS_SOIL(1,2) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,3) = \"\"/CFNAM_ALBVIS_SOIL(1,3) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,4) = \"\"/CFNAM_ALBVIS_SOIL(1,4) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,5) = \"\"/CFNAM_ALBVIS_SOIL(1,5) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,6) = \"\"/CFNAM_ALBVIS_SOIL(1,6) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,7) = \"\"/CFNAM_ALBVIS_SOIL(1,7) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,8) = \"\"/CFNAM_ALBVIS_SOIL(1,8) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,9) = \"\"/CFNAM_ALBVIS_SOIL(1,9) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,10) = \"\"/CFNAM_ALBVIS_SOIL(1,10) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,11) = \"\"/CFNAM_ALBVIS_SOIL(1,11) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBVIS_SOIL(1,12) = \"\"/CFNAM_ALBVIS_SOIL(1,12) = \"ALBVIS_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam

#albuv_soil
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,1) = \"\"/CFNAM_ALBUV_SOIL(1,1) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,2) = \"\"/CFNAM_ALBUV_SOIL(1,2) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,3) = \"\"/CFNAM_ALBUV_SOIL(1,3) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,4) = \"\"/CFNAM_ALBUV_SOIL(1,4) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,5) = \"\"/CFNAM_ALBUV_SOIL(1,5) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,6) = \"\"/CFNAM_ALBUV_SOIL(1,6) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,7) = \"\"/CFNAM_ALBUV_SOIL(1,7) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,8) = \"\"/CFNAM_ALBUV_SOIL(1,8) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,9) = \"\"/CFNAM_ALBUV_SOIL(1,9) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,10) = \"\"/CFNAM_ALBUV_SOIL(1,10) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,11) = \"\"/CFNAM_ALBUV_SOIL(1,11) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALBUV_SOIL(1,12) = \"\"/CFNAM_ALBUV_SOIL(1,12) = \"ALBUV_SOIL_M7_P0.dat\"/g" OPTIONS.nam_save > OPTIONS.nam

#gmes
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GMES(1) = 1.0E+20/XUNIF_GMES(1) = 0.000/g" OPTIONS.nam_save > OPTIONS.nam

#re25
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RE25(1) = 1.0E+20/XUNIF_RE25(1) = 0.0000001412/g" OPTIONS.nam_save > OPTIONS.nam

#bslai
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_BSLAI(1) = 1.0E+20/XUNIF_BSLAI(1) = 0./g" OPTIONS.nam_save > OPTIONS.nam

#laimin
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAIMIN(1) = 1.0E+20/XUNIF_LAIMIN(1) =0./g" OPTIONS.nam_save > OPTIONS.nam

#sefold
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_SEFOLD(1) = 1.0E+20/XUNIF_SEFOLD(1) = 0./g" OPTIONS.nam_save > OPTIONS.nam

#gc
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GC(1) = 1.0E+20/XUNIF_GC(1) = 0./g" OPTIONS.nam_save > OPTIONS.nam

#dmax
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DMAX(1) = 1.0E+20/XUNIF_DMAX(1) = 0./g" OPTIONS.nam_save > OPTIONS.nam

#f2i
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_F2I(1) = 1.0E+20/XUNIF_F2I(1) = 0./g" OPTIONS.nam_save > OPTIONS.nam

#h_tree
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_TREE(1) = 1.0E+20/XUNIF_H_TREE(1) = 0./g" OPTIONS.nam_save > OPTIONS.nam

#ce_nitro
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CE_NITRO(1) = 1.0E+20/XUNIF_CE_NITRO(1) = 0./g" OPTIONS.nam_save > OPTIONS.nam

#cf_nitro
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CF_NITRO(1) = 1.0E+20/XUNIF_CF_NITRO(1) = 0./g" OPTIONS.nam_save > OPTIONS.nam

#cna_nitro
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CNA_NITRO(1) = 1.0E+20/XUNIF_CNA_NITRO(1) = 0./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
