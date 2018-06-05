cp -f $1 OPTIONS.nam_save

sed -e "s/CFILE_TS = \"\"/CFILE_TS = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_TS = \"\"/CTYPE_TS = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XWS_ROOF = 1.0E+20/XWS_ROOF = 0.95/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XWS_ROAD = 1.0E+20/XWS_ROAD = 0.92/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTI_BLD = 1.0E+20/XTI_BLD = 285./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTI_ROAD = 1.0E+20/XTI_ROAD = 285./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUI_BLD = 1.0E+20/XHUI_BLD = 0.02/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save

cp -f OPTIONS.nam $1
