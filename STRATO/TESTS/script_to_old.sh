cp -f OPTIONS.nam OPTIONS.nam_new

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOWRAD = \"T17\"/CSNOWRAD = \"TAR\" /g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOWDRIFT = \" DFLT\"/LSNOWDRIFT = F /g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOWDRIFT = \"NONE\"/LSNOWDRIFT = T /g"  OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_old
