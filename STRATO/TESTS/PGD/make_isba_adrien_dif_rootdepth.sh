#rootdepth
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOT_DEPTH(1) = 1.0E+20/XUNIF_ROOT_DEPTH(1) = 3./g" OPTIONS.nam_save > OPTIONS.nam

#grounddepth
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GROUND_DEPTH(1) = 1.0E+20/XUNIF_GROUND_DEPTH(1) = 6./g" OPTIONS.nam_save > OPTIONS.nam

#root_extinction
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOT_EXTINCTION(1) = 1.0E+20/XUNIF_ROOT_EXTINCTION(1) = 0.95/g" OPTIONS.nam_save > OPTIONS.nam

#root_lin
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ROOT_LIN(1) = 1.0E+20/XUNIF_ROOT_LIN(1) = 0.1/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
