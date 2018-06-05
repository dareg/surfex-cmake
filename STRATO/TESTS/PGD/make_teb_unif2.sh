#roof
cp -f $1 OPTIONS.nam_save
sed -e "s/NPAR_ROOF_LAYER = 0/NPAR_ROOF_LAYER = 5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALB_ROOF = 1.0E+20/XUNIF_ALB_ROOF = 0.1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_EMIS_ROOF = 1.0E+20/XUNIF_EMIS_ROOF = 0.65/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROOF(1) = 1.0E+20/XUNIF_HC_ROOF(1) = 1580000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROOF(2) = 1.0E+20/XUNIF_HC_ROOF(2) = 1580000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROOF(3) = 1.0E+20/XUNIF_HC_ROOF(3) = 1580000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROOF(4) = 1.0E+20/XUNIF_HC_ROOF(4) = 1580000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROOF(5) = 1.0E+20/XUNIF_HC_ROOF(5) = 52030./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROOF(1) = 1.0E+20/XUNIF_TC_ROOF(1) = 1.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROOF(2) = 1.0E+20/XUNIF_TC_ROOF(2) = 1.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROOF(3) = 1.0E+20/XUNIF_TC_ROOF(3) = 1.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROOF(4) = 1.0E+20/XUNIF_TC_ROOF(4) = 1.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROOF(5) = 1.0E+20/XUNIF_TC_ROOF(5) = 0.03/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROOF(1) = 1.0E+20/XUNIF_D_ROOF(1) = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROOF(2) = 1.0E+20/XUNIF_D_ROOF(2) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROOF(3) = 1.0E+20/XUNIF_D_ROOF(3) = 0.18/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROOF(4) = 1.0E+20/XUNIF_D_ROOF(4) = 0.06/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROOF(5) = 1.0E+20/XUNIF_D_ROOF(5) = 0.03/g" OPTIONS.nam_save > OPTIONS.nam

#road
cp -f $1 OPTIONS.nam_save
sed -e "s/NPAR_ROAD_LAYER = 0/NPAR_ROAD_LAYER = 4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALB_ROAD = 1.0E+20/XUNIF_ALB_ROAD = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_EMIS_ROAD = 1.0E+20/XUNIF_EMIS_ROAD = 0.85/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROAD(1) = 1.0E+20/XUNIF_HC_ROAD(1) = 1740000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROAD(2) = 1.0E+20/XUNIF_HC_ROAD(2) = 1740000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROAD(3) = 1.0E+20/XUNIF_HC_ROAD(3) = 2000000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_ROAD(4) = 1.0E+20/XUNIF_HC_ROAD(4) = 1400000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROAD(1) = 1.0E+20/XUNIF_TC_ROAD(1) = 0.82/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROAD(2) = 1.0E+20/XUNIF_TC_ROAD(2) = 0.82/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROAD(3) = 1.0E+20/XUNIF_TC_ROAD(3) = 2.1000/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_ROAD(4) = 1.0E+20/XUNIF_TC_ROAD(4) = 0.4000/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROAD(1) = 1.0E+20/XUNIF_D_ROAD(1) = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROAD(2) = 1.0E+20/XUNIF_D_ROAD(2) = 0.04/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROAD(3) = 1.0E+20/XUNIF_D_ROAD(3) = 0.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_ROAD(4) = 1.0E+20/XUNIF_D_ROAD(4) = 1.0/g" OPTIONS.nam_save > OPTIONS.nam

#wall
cp -f $1 OPTIONS.nam_save
sed -e "s/NPAR_WALL_LAYER = 0/NPAR_WALL_LAYER = 4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ALB_WALL = 1.0E+20/XUNIF_ALB_WALL = 0.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_EMIS_WALL = 1.0E+20/XUNIF_EMIS_WALL = 0.65/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_WALL(1) = 1.0E+20/XUNIF_HC_WALL(1) = 1344000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_WALL(2) = 1.0E+20/XUNIF_HC_WALL(2) = 1344000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_WALL(3) = 1.0E+20/XUNIF_HC_WALL(3) = 15125./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HC_WALL(4) = 1.0E+20/XUNIF_HC_WALL(4) = 47325./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_WALL(1) = 1.0E+20/XUNIF_TC_WALL(1) =1.45/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_WALL(2) = 1.0E+20/XUNIF_TC_WALL(2) = 1.45/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_WALL(3) = 1.0E+20/XUNIF_TC_WALL(3) = 0.04/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TC_WALL(4) = 1.0E+20/XUNIF_TC_WALL(4) = 0.04/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_WALL(1) = 1.0E+20/XUNIF_D_WALL(1) = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_WALL(2) = 1.0E+20/XUNIF_D_WALL(2) = 0.23/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_WALL(3) = 1.0E+20/XUNIF_D_WALL(3) = 0.13/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_WALL(4) = 1.0E+20/XUNIF_D_WALL(4) = 0.02/g" OPTIONS.nam_save > OPTIONS.nam


cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_BLD = 1.0E+20/XUNIF_BLD = 0.35/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_BLD_HEIGHT = 1.0E+20/XUNIF_BLD_HEIGHT = 14.9/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_WALL_O_HOR = 1.0E+20/XUNIF_WALL_O_HOR = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_Z0_TOWN = 1.0E+20/XUNIF_Z0_TOWN = 1.5/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_TRAFFIC = 1.0E+20/XUNIF_H_TRAFFIC = 8./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LE_TRAFFIC = 1.0E+20/XUNIF_LE_TRAFFIC = 0.1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_H_INDUSTRY = 1.0E+20/XUNIF_H_INDUSTRY = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_LE_INDUSTRY = 1.0E+20/XUNIF_LE_INDUSTRY = 0./g" OPTIONS.nam_save > OPTIONS.nam


cp -f OPTIONS.nam $1
