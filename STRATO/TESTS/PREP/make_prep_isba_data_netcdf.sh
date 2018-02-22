cp -f $1 OPTIONS.nam_save

sed -e "s/CFILE_HUG_SURF = \"\"/CFILE_HUG_SURF = \"WG1.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_HUG_ROOT = \"\"/CFILE_HUG_ROOT = \"WG2.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_HUG_DEEP = \"\"/CFILE_HUG_DEEP = \"WG3.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_SURF = 1.0E+20/XHUGI_SURF = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_ROOT = 1.0E+20/XHUGI_ROOT = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_DEEP = 1.0E+20/XHUGI_DEEP = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG = \"\"/CFILE_TG = \"TG_nc.nc\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_TG = \"\"/CTYPE_TG = \"NETCDF\"/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
