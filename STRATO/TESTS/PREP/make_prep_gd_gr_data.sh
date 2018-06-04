cp -f $1 OPTIONS.nam_save

sed -e "s/CFILE_HUG_SURF_GD = \"\"/CFILE_HUG_SURF_GD = \"WG1.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_HUG_ROOT_GD = \"\"/CFILE_HUG_ROOT_GD = \"WG2.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_HUG_DEEP_GD = \"\"/CFILE_HUG_DEEP_GD = \"WG2.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_SURF_GD = 1.0E+20/XHUGI_SURF_GD = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_ROOT_GD = 1.0E+20/XHUGI_ROOT_GD = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_DEEP_GD = 1.0E+20/XHUGI_DEEP_GD = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG_SURF_GD = \"\"/CFILE_TG_SURF_GD = \"TG1.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG_ROOT_GD = \"\"/CFILE_TG_ROOT_GD = \"TG2.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG_DEEP_GD = \"\"/CFILE_TG_DEEP_GD = \"TG2.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_HUG = \"\"/CTYPE_HUG = \"ASCLLV\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_TG = \"\"/CTYPE_TG = \"ASCLLV\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

sed -e "s/CFILE_HUG_SURF_GR = \"\"/CFILE_HUG_SURF_GR = \"WG1.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_HUG_ROOT_GR = \"\"/CFILE_HUG_ROOT_GR = \"WG2.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_HUG_DEEP_GR = \"\"/CFILE_HUG_DEEP_GR = \"WG2.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_SURF_GR = 1.0E+20/XHUGI_SURF_GR = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_ROOT_GR = 1.0E+20/XHUGI_ROOT_GR = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_DEEP_GR = 1.0E+20/XHUGI_DEEP_GR = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG_SURF_GR = \"\"/CFILE_TG_SURF_GR = \"TG1.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG_ROOT_GR = \"\"/CFILE_TG_ROOT_GR = \"TG2.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG_DEEP_GR = \"\"/CFILE_TG_DEEP_GR = \"TG2.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_HUG = \"\"/CTYPE_HUG = \"ASCLLV\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_TG = \"\"/CTYPE_TG = \"ASCLLV\"/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
