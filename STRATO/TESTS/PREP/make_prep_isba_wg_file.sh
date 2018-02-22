cp -f $1 OPTIONS.nam_save

sed -e "s/CFILE_HUG = \"\"/CFILE_HUG = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_HUG = \"\"/CTYPE_HUG = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_SURF = 1.0E+20/XHUGI_SURF = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_ROOT = 1.0E+20/XHUGI_ROOT = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_DEEP = 1.0E+20/XHUGI_DEEP = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_SURF = 1.0E+20/XTG_SURF = 290.16/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_ROOT = 1.0E+20/XTG_ROOT = 285.52/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_DEEP = 1.0E+20/XTG_DEEP = 273.15/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
