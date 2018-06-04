#time
cp -f $1 OPTIONS.nam_save
sed -e "s/NTIME\ \=\ 36/NTIME\ \=\ 12/g" OPTIONS.nam_save > OPTIONS.nam

#vegtypes
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 19
do
	var="XUNIF_VEGTYPE($i)"
	sed -e "s/$var = 1.0E+20/$var = 0./g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
sed -e "s/XUNIF_VEGTYPE(7) = 1.0E+20/XUNIF_VEGTYPE(7) = 1./g" OPTIONS.nam_save > OPTIONS.nam

#veg
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 10 11 12
do
	var="XUNIF_VEG(7,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0./g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
for i in 6 7 8 9
do
	var="XUNIF_VEG(7,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.9/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
sed -e "s/XUNIF_VEG(7,5) = 1.0E+20/XUNIF_VEG(7,5) = 0.5/g" OPTIONS.nam_save > OPTIONS.nam

#lai
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 10 11 12
do
	var="XUNIF_LAI(7,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0./g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
for i in 6 7 8 9
do
	var="XUNIF_LAI(7,$i)"
	sed -e "s/$var = 1.0E+20/$var = 3./g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
sed -e "s/XUNIF_LAI(7,5) = 1.0E+20/XUNIF_LAI(7,5) = 1./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

#z0
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 10 11 12
do
	var="XUNIF_Z0(7,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
for i in 6 7 8 9
do
	var="XUNIF_Z0(7,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.15/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
sed -e "s/XUNIF_Z0(7,5) = 1.0E+20/XUNIF_Z0(7,5) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

#emis
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
	var="XUNIF_EMIS(7,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.98/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

#dg
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(7,1) = 1.0E+20/XUNIF_DG(7,1) = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(7,2) = 1.0E+20/XUNIF_DG(7,2) = 1.60/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(7,3) = 1.0E+20/XUNIF_DG(7,3) = 1.60/g" OPTIONS.nam_save > OPTIONS.nam

#rsmin
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RSMIN(7) = 1.0E+20/XUNIF_RSMIN(7) = 40./g" OPTIONS.nam_save > OPTIONS.nam

#gamma
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GAMMA(7) = 1.0E+20/XUNIF_GAMMA(7) = 0./g" OPTIONS.nam_save > OPTIONS.nam

#wrmax_cf
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_WRMAX_CF(7) = 1.0E+20/XUNIF_WRMAX_CF(7) = 0.2/g" OPTIONS.nam_save > OPTIONS.nam

#rgl
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RGL(7) = 1.0E+20/XUNIF_RGL(7) = 100./g" OPTIONS.nam_save > OPTIONS.nam

#cv
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CV(7) = 1.0E+20/XUNIF_CV(7) = 0.00002/g" OPTIONS.nam_save > OPTIONS.nam

#z0_O_Z0H
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0_O_Z0H(7) = 1.0E+20/XUNIF_Z0_O_Z0H(7) = 10./g" OPTIONS.nam_save > OPTIONS.nam

#albnir_veg
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do

	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/XUNIF_ALBNIR_VEG(7,$i) = 1.0E+20/XUNIF_ALBNIR_VEG(7,$i) = 0.3/g" OPTIONS.nam_save > OPTIONS.nam

#albvis_veg
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/XUNIF_ALBVIS_VEG(7,$i) = 1.0E+20/XUNIF_ALBVIS_VEG(7,$i) = 0.1/g" OPTIONS.nam_save > OPTIONS.nam

#albuv_veg
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/XUNIF_ALBUV_VEG(7,$i) = 1.0E+20/XUNIF_ALBUV_VEG(7,$i) = 0.0425/g" OPTIONS.nam_save > OPTIONS.nam

#albnir_soil
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/XUNIF_ALBNIR_SOIL(7,$i) = 1.0E+20/XUNIF_ALBNIR_SOIL(7,$i) = 0.3/g" OPTIONS.nam_save > OPTIONS.nam

#albvis_soil
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/XUNIF_ALBVIS_SOIL(7,$i) = 1.0E+20/XUNIF_ALBVIS_SOIL(7,$i) = 0.1/g" OPTIONS.nam_save > OPTIONS.nam

#albuv_soil
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/XUNIF_ALBUV_SOIL(7,$i) = 1.0E+20/XUNIF_ALBUV_SOIL(7,$i) = 0.06/g" OPTIONS.nam_save > OPTIONS.nam

done

#gmes
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GMES(7) = 1.0E+20/XUNIF_GMES(7) = 0.003/g" OPTIONS.nam_save > OPTIONS.nam

#re25
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RE25(7) = 1.0E+20/XUNIF_RE25(7) = 0.0000003/g" OPTIONS.nam_save > OPTIONS.nam

#bslai
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_BSLAI(7) = 1.0E+20/XUNIF_BSLAI(7) = 0.06/g" OPTIONS.nam_save > OPTIONS.nam

#laimin
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LAIMIN(7) = 1.0E+20/XUNIF_LAIMIN(7) =0.3/g" OPTIONS.nam_save > OPTIONS.nam

#sefold
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_SEFOLD(7) = 1.0E+20/XUNIF_SEFOLD(7) = 5184000./g" OPTIONS.nam_save > OPTIONS.nam

#gc
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GC(7) = 1.0E+20/XUNIF_GC(7) = 0.00025/g" OPTIONS.nam_save > OPTIONS.nam

#dmax
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DMAX(7) = 1.0E+20/XUNIF_DMAX(7) = 0.1/g" OPTIONS.nam_save > OPTIONS.nam

#f2i
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_F2I(7) = 1.0E+20/XUNIF_F2I(7) = 0.3/g" OPTIONS.nam_save > OPTIONS.nam

#h_tree
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_TREE(7) = 1.0E+20/XUNIF_H_TREE(7) = 20./g" OPTIONS.nam_save > OPTIONS.nam

#ce_nitro
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CE_NITRO(7) = 1.0E+20/XUNIF_CE_NITRO(7) = 3.79/g" OPTIONS.nam_save > OPTIONS.nam

#cf_nitro
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CF_NITRO(7) = 1.0E+20/XUNIF_CF_NITRO(7) = 9.84/g" OPTIONS.nam_save > OPTIONS.nam

#cna_nitro
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CNA_NITRO(7) = 1.0E+20/XUNIF_CNA_NITRO(7) = 1.3/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
