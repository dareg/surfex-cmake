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
sed -e "s/XUNIF_DG(1,8) = 1.0E+20/XUNIF_DG(1,8) = 6.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,9) = 1.0E+20/XUNIF_DG(1,9) = 7.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,10) = 1.0E+20/XUNIF_DG(1,10) = 7.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,11) = 1.0E+20/XUNIF_DG(1,11) = 8.1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,12) = 1.0E+20/XUNIF_DG(1,12) = 9.3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,13) = 1.0E+20/XUNIF_DG(1,13) = 10.1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_DG(1,14) = 1.0E+20/XUNIF_DG(1,14) = 11.5/g" OPTIONS.nam_save > OPTIONS.nam

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
sed -e "s/XUNIF_ROOTFRAC(1,5) = 1.0E+20/XUNIF_ROOTFRAC(1,5) = 0.4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,6) = 1.0E+20/XUNIF_ROOTFRAC(1,6) = 0.47/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,7) = 1.0E+20/XUNIF_ROOTFRAC(1,7) = 0.55/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,8) = 1.0E+20/XUNIF_ROOTFRAC(1,8) = 0.61/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,9) = 1.0E+20/XUNIF_ROOTFRAC(1,9) = 0.70/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,10) = 1.0E+20/XUNIF_ROOTFRAC(1,10) = 0.83/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,11) = 1.0E+20/XUNIF_ROOTFRAC(1,11) = 0.94/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,12) = 1.0E+20/XUNIF_ROOTFRAC(1,12) = 1.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,13) = 1.0E+20/XUNIF_ROOTFRAC(1,13) = 1.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOTFRAC(1,14) = 1.0E+20/XUNIF_ROOTFRAC(1,14) = 1.0/g" OPTIONS.nam_save > OPTIONS.nam


cp -f OPTIONS.nam $1
