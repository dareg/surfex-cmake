cp -f $1 OPTIONS.nam_save

sed -e "s/CFNAM_BLDTYPE = \"\"/CFNAM_BLDTYPE = \"files\/bldtype.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_BLD_AGE = \"\"/CFNAM_BLD_AGE = \"files\/age_renovation.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_USETYPE = \"\"/CFNAM_USETYPE = \"files\/usetype.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CCSVDATAFILE = \"\"/CCSVDATAFILE = \"files\/BATI_HYP3.csv\"/g" OPTIONS.nam_save > OPTIONS.nam 


cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_BLD = \"\"/CFNAM_BLD = \"files\/bld_dens.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_BLD_HEIGHT = \"\"/CFNAM_BLD_HEIGHT = \"files\/height.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_WALL_O_HOR = \"\"/CFNAM_WALL_O_HOR = \"files\/wall.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_Z0_TOWN = \"\"/CFNAM_Z0_TOWN = \"files\/z0.txt\"/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ROAD_DIR = \"\"/CFNAM_ROAD_DIR = \"files\/road_dir.txt\"/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_GARDEN = \"\"/CFNAM_GARDEN = \"files\/garden.txt\"/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_TRAFFIC = 1.0E+20/XUNIF_H_TRAFFIC = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LE_TRAFFIC = 1.0E+20/XUNIF_LE_TRAFFIC = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_INDUSTRY = 1.0E+20/XUNIF_H_INDUSTRY = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LE_INDUSTRY = 1.0E+20/XUNIF_LE_INDUSTRY = 0./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
