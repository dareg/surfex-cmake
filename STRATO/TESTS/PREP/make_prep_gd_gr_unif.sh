cp -f $1 OPTIONS.nam_save

sed -e "s/XHUG_SURF_GD = 1.0E+20/XHUG_SURF_GD = 1./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUG_ROOT_GD = 1.0E+20/XHUG_ROOT_GD = 1./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUG_DEEP_GD = 1.0E+20/XHUG_DEEP_GD = 1./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_SURF_GD = 1.0E+20/XHUGI_SURF_GD = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_ROOT_GD = 1.0E+20/XHUGI_ROOT_GD = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_DEEP_GD = 1.0E+20/XHUGI_DEEP_GD = 0.02/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_SURF_GD = 1.0E+20/XTG_SURF_GD = 289.12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_ROOT_GD = 1.0E+20/XTG_ROOT_GD = 289.12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_DEEP_GD = 1.0E+20/XTG_DEEP_GD = 273.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

sed -e "s/XHUG_SURF_GR = 1.0E+20/XHUG_SURF_GR = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUG_ROOT_GR = 1.0E+20/XHUG_ROOT_GR = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUG_DEEP_GR = 1.0E+20/XHUG_DEEP_GR = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_SURF_GR = 1.0E+20/XHUGI_SURF_GR = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_ROOT_GR = 1.0E+20/XHUGI_ROOT_GR = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_DEEP_GR = 1.0E+20/XHUGI_DEEP_GR = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_SURF_GR = 1.0E+20/XTG_SURF_GR = 291.45/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_ROOT_GR = 1.0E+20/XTG_ROOT_GR = 291.45/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_DEEP_GR = 1.0E+20/XTG_DEEP_GR = 291.45/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
