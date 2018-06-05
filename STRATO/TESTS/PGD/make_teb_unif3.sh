#roof
cp -f $1 OPTIONS.nam_save
sed -e "s/NPAR_ROOF_LAYER = 0/NPAR_ROOF_LAYER = 5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALB_ROOF = 1.0E+20/XUNIF_ALB_ROOF = 0.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_EMIS_ROOF = 1.0E+20/XUNIF_EMIS_ROOF = 0.80/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROOF(1) = 1.0E+20/XUNIF_HC_ROOF(1) = 2100000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROOF(2) = 1.0E+20/XUNIF_HC_ROOF(2) = 44800./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROOF(3) = 1.0E+20/XUNIF_HC_ROOF(3) = 2100000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROOF(4) = 1.0E+20/XUNIF_HC_ROOF(4) = 75000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROOF(5) = 1.0E+20/XUNIF_HC_ROOF(5) = 2300000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROOF(1) = 1.0E+20/XUNIF_TC_ROOF(1) = 0.7/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROOF(2) = 1.0E+20/XUNIF_TC_ROOF(2) = 0.024/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROOF(3) = 1.0E+20/XUNIF_TC_ROOF(3) = 0.7/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROOF(4) = 1.0E+20/XUNIF_TC_ROOF(4) = 0.035/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROOF(5) = 1.0E+20/XUNIF_TC_ROOF(5) = 2.3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROOF(1) = 1.0E+20/XUNIF_D_ROOF(1) = 0.003/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROOF(2) = 1.0E+20/XUNIF_D_ROOF(2) = 0.060/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROOF(3) = 1.0E+20/XUNIF_D_ROOF(3) = 0.003/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROOF(4) = 1.0E+20/XUNIF_D_ROOF(4) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROOF(5) = 1.0E+20/XUNIF_D_ROOF(5) = 0.20/g" OPTIONS.nam_save > OPTIONS.nam

#road
cp -f $1 OPTIONS.nam_save
sed -e "s/NPAR_ROAD_LAYER = 0/NPAR_ROAD_LAYER = 4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALB_ROAD = 1.0E+20/XUNIF_ALB_ROAD = 0.08/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_EMIS_ROAD = 1.0E+20/XUNIF_EMIS_ROAD = 0.94/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROAD(1) = 1.0E+20/XUNIF_HC_ROAD(1) = 1825000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROAD(2) = 1.0E+20/XUNIF_HC_ROAD(2) = 2000000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROAD(3) = 1.0E+20/XUNIF_HC_ROAD(3) = 1400000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROAD(4) = 1.0E+20/XUNIF_HC_ROAD(4) = 1400000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROAD(1) = 1.0E+20/XUNIF_TC_ROAD(1) = 0.663/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROAD(2) = 1.0E+20/XUNIF_TC_ROAD(2) = 2.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROAD(3) = 1.0E+20/XUNIF_TC_ROAD(3) = 0.40/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROAD(4) = 1.0E+20/XUNIF_TC_ROAD(4) = 0.40/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROAD(1) = 1.0E+20/XUNIF_D_ROAD(1) = 0.04/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROAD(2) = 1.0E+20/XUNIF_D_ROAD(2) = 0.20/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROAD(3) = 1.0E+20/XUNIF_D_ROAD(3) = 0.50/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROAD(4) = 1.0E+20/XUNIF_D_ROAD(4) = 0.50/g" OPTIONS.nam_save > OPTIONS.nam

#wall


cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_BLD = 1.0E+20/XUNIF_BLD = 0.689/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_WALL_O_HOR = 1.0E+20/XUNIF_WALL_O_HOR = 0.75/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_TRAFFIC = 1.0E+20/XUNIF_H_TRAFFIC = 8.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LE_TRAFFIC = 1.0E+20/XUNIF_LE_TRAFFIC = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_INDUSTRY = 1.0E+20/XUNIF_H_INDUSTRY = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LE_INDUSTRY = 1.0E+20/XUNIF_LE_INDUSTRY = 0./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
