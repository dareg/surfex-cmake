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
sed -e "s/CFILE_TG_GD = \"\"/CFILE_TG_GD = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG_GR = \"\"/CFILE_TG_GR = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_TG = \"\"/CTYPE_TG = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam

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

cp -f OPTIONS.nam $1
