cp -f $1 OPTIONS.nam_save

sed -e "s/XWSNOW = 1.0E+20/XWSNOW = 10.,50.,100.,840.,9000.,40000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTSNOW = 1.0E+20/XTSNOW = 273.16,273.12,273.14,273.11,273.16,273.11/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRSNOW = 1.0E+20/XRSNOW = 900.,910.,930.,910.,940.,935./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XASNOW = 1.0E+20/XASNOW = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSG1SNOW = 1.0E+20/XSG1SNOW = 99.,97.,96.,98.,99.,97./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSG2SNOW = 1.0E+20/XSG2SNOW = 0.005,0.006,0.004,0.007,0.005,0.006/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHISTSNOW = 1.0E+20/XHISTSNOW = 2.,1.,3.,2.,2.,1./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XAGESNOW = 1.0E+20/XAGESNOW = 1000.,990.,1010.,1005.,995.,1000./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
