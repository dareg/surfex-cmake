#roof
cp -f $1 OPTIONS.nam_save
sed -e "s/NPAR_ROOF_LAYER = 0/NPAR_ROOF_LAYER = 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALB_ROOF = \"\"/CFNAM_ALB_ROOF = \"XDATA_ALB_ROOF.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_EMIS_ROOF = \"\"/CFNAM_EMIS_ROOF = \"XDATA_EMIS_ROOF.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
for i in "HC" "TC" "D"
do
	for j in 1 2 3
	do
		var1="CFNAM_${i}_ROOF($j)"
		var2="XDATA_${i}_ROOF$j.txt"
		cp -f OPTIONS.nam OPTIONS.nam_save
		sed -e "s/$var1 = \"\"/$var1 = \"$var2\"/g" OPTIONS.nam_save > OPTIONS.nam
	done
done


#road
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NPAR_ROAD_LAYER = 0/NPAR_ROAD_LAYER = 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALB_ROAD = 1.0E+20/XUNIF_ALB_ROAD = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_EMIS_ROAD = 1.0E+20/XUNIF_EMIS_ROAD = 0.94/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROAD(1) = 1.0E+20/XUNIF_HC_ROAD(1) = 2000000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROAD(2) = 1.0E+20/XUNIF_HC_ROAD(2) = 2000000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROAD(3) = 1.0E+20/XUNIF_HC_ROAD(3) = 1800000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROAD(1) = 1.0E+20/XUNIF_TC_ROAD(1) = 2.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROAD(2) = 1.0E+20/XUNIF_TC_ROAD(2) = 1.50/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROAD(3) = 1.0E+20/XUNIF_TC_ROAD(3) = 0.25/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROAD(1) = 1.0E+20/XUNIF_D_ROAD(1) = 0.04/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROAD(2) = 1.0E+20/XUNIF_D_ROAD(2) = 0.37/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROAD(3) = 1.0E+20/XUNIF_D_ROAD(3) = 1.00/g" OPTIONS.nam_save > OPTIONS.nam

#wall
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NPAR_WALL_LAYER = 0/NPAR_WALL_LAYER = 4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFNAM_ALB_WALL = \"\"/CFNAM_ALB_WALL = \"XDATA_ALB_WALL.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_EMIS_WALL = 1.0E+20/XUNIF_EMIS_WALL = 0.90/g" OPTIONS.nam_save > OPTIONS.nam
for i in "HC" "TC" "D"
do
	for j in 1 2 3 4
	do
		var1="CFNAM_${i}_ROOF($j)"
		var2="XDATA_${i}_ROOF$j.txt"
		cp -f OPTIONS.nam OPTIONS.nam_save
		sed -e "s/$var1 = \"\"/$var1 = \"$var2\"/g" OPTIONS.nam_save > OPTIONS.nam
	done
done

for i in "BLD" "BLD_HEIGHT" "WALL_O_HOR" "H_TRAFFIC" "LE_TRAFFIC"
do
	var1="CFNAM_$i"
	var2="XDATA_$i.txt"	
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/$var1 = \"\"/$var1 = \"$var2\"/g" OPTIONS.nam_save > OPTIONS.nam
done

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0_TOWN = 1.0E+20/XUNIF_Z0_TOWN = 1./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_INDUSTRY = 1.0E+20/XUNIF_H_INDUSTRY = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LE_INDUSTRY = 1.0E+20/XUNIF_LE_INDUSTRY = 0./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
