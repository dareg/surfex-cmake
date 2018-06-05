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

cp -f OPTIONS.nam $1
