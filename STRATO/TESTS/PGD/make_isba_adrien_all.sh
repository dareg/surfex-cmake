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


#z0
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
do
	var="XUNIF_Z0(1,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.00110514315659256/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,1) = 1.0E+20/XUNIF_Z0(1,1) = 0.0086/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,2) = 1.0E+20/XUNIF_Z0(1,2) = 0.0094/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,3) = 1.0E+20/XUNIF_Z0(1,3) = 0.0143679046080891/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,4) = 1.0E+20/XUNIF_Z0(1,4) = 0.0176833340778667/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,5) = 1.0E+20/XUNIF_Z0(1,5) = 0.022103906704237/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,6) = 1.0E+20/XUNIF_Z0(1,6) = 0.0276296224871998/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,7) = 1.0E+20/XUNIF_Z0(1,7) = 0.0353656245833476/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,8) = 1.0E+20/XUNIF_Z0(1,8) = 0.0442078134084739/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,9) = 1.0E+20/XUNIF_Z0(1,9) = 0.0530489586612144/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,10) = 1.0E+20/XUNIF_Z0(1,10) = 0.0618911474863407/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,11) = 1.0E+20/XUNIF_Z0(1,11) = 0.0707322927390811/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,12) = 1.0E+20/XUNIF_Z0(1,12) = 0.0795734379918216/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,13) = 1.0E+20/XUNIF_Z0(1,13) = 0.0740477222088588/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,14) = 1.0E+20/XUNIF_Z0(1,14) = 0.068522006425896/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,15) = 1.0E+20/XUNIF_Z0(1,15) = 0.0629962906429332/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,16) = 1.0E+20/XUNIF_Z0(1,16) = 0.0574695312875846/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,17) = 1.0E+20/XUNIF_Z0(1,17) = 0.0519438155046218/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,18) = 1.0E+20/XUNIF_Z0(1,18) = 0.0464180997216591/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0(1,19) = 1.0E+20/XUNIF_Z0(1,19) = 0.0386820976255111/g" OPTIONS.nam_save > OPTIONS.nam


#emis
cp -f OPTIONS.nam OPTIONS.nam_save
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
do
	var="XUNIF_EMIS(7,$i)"
	sed -e "s/$var = 1.0E+20/$var = 0.98/g" OPTIONS.nam_save > OPTIONS.nam
	cp -f OPTIONS.nam OPTIONS.nam_save
done

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

#rootfrac
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,1) = 1.0E+20/XUNIF_ROOTFRAC(1,1) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,2) = 1.0E+20/XUNIF_ROOTFRAC(1,2) = 0.12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,3) = 1.0E+20/XUNIF_ROOTFRAC(1,3) = 0.26/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,4) = 1.0E+20/XUNIF_ROOTFRAC(1,4) = 0.35/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,5) = 1.0E+20/XUNIF_ROOTFRAC(1,5) = 0.6/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,6) = 1.0E+20/XUNIF_ROOTFRAC(1,6) = 0.8/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,7) = 1.0E+20/XUNIF_ROOTFRAC(1,7) = 1.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,8) = 1.0E+20/XUNIF_ROOTFRAC(1,8) = 1.0/g" OPTIONS.nam_save > OPTIONS.nam

#rootdepth
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOT_DEPTH(1) = 1.0E+20/XUNIF_ROOT_DEPTH(1) = 3./g" OPTIONS.nam_save > OPTIONS.nam

#grounddepth
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GROUND_DEPTH(1) = 1.0E+20/XUNIF_GROUND_DEPTH(1) = 6./g" OPTIONS.nam_save > OPTIONS.nam

#root_extinction
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOT_EXTINCTION(1) = 1.0E+20/XUNIF_ROOT_EXTINCTION(1) = 0.7/g" OPTIONS.nam_save > OPTIONS.nam

#root_lin
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOT_LIN(1) = 1.0E+20/XUNIF_ROOT_LIN(1) = 1.1/g" OPTIONS.nam_save > OPTIONS.nam

#rsmin
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_RSMIN(1) = 1.0E+20/XUNIF_RSMIN(1) = 125./g" OPTIONS.nam_save > OPTIONS.nam

#gamma
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GAMMA(1) = 1.0E+20/XUNIF_GAMMA(1) = 0.0345/g" OPTIONS.nam_save > OPTIONS.nam

#cv
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CV(1) = 1.0E+20/XUNIF_CV(1) = 0.000011582/g" OPTIONS.nam_save > OPTIONS.nam

#h_tree
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_TREE(1) = 1.0E+20/XUNIF_H_TREE(1) = 20./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
