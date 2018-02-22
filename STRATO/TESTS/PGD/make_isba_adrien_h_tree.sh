#vegtypes
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_TREE(1) = 1.0E+20/XUNIF_H_TREE(1) = 20./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
