#roof
cp -f $1 OPTIONS.nam_save
sed -e "s/NPAR_FLOOR_LAYER = 0/NPAR_FLOOR_LAYER = 5/g" OPTIONS.nam_save > OPTIONS.nam

for i in 1 2 3 4 5
do
	var1="XUNIF_HC_FLOOR($i)"
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/$var1 = 1.0E+20/$var1 = 2016000./g" OPTIONS.nam_save > OPTIONS.nam
done
for i in 1 2 3 4 5
do
	var1="XUNIF_TC_FLOOR($i)"
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/$var1 = 1.0E+20/$var1 = 1.95/g" OPTIONS.nam_save > OPTIONS.nam
done
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_FLOOR(1) = 1.0E+20/XUNIF_D_FLOOR(1) = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_FLOOR(2) = 1.0E+20/XUNIF_D_FLOOR(2) = 0.04/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_FLOOR(3) = 1.0E+20/XUNIF_D_FLOOR(3) = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_FLOOR(4) = 1.0E+20/XUNIF_D_FLOOR(4) = 0.04/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_FLOOR(5) = 1.0E+20/XUNIF_D_FLOOR(5) = 0.01/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TCOOL_TARGET = 1.0E+20/XUNIF_TCOOL_TARGET = 297.16/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_THEAT_TARGET = 1.0E+20/XUNIF_THEAT_TARGET = 297.16/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_F_WASTE_CAN = 1.0E+20/XUNIF_F_WASTE_CAN = 1.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_COP_RAT = 1.0E+20/XUNIF_COP_RAT = 2.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_EFF_HEAT = 1.0E+20/XUNIF_EFF_HEAT = 0.9/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_QIN = 1.0E+20/XUNIF_QIN = 5.8/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_QIN_FRAD = 1.0E+20/XUNIF_QIN_FRAD = 0.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_QIN_FLAT = 1.0E+20/XUNIF_QIN_FLAT = 0.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_SHGC = 1.0E+20/XUNIF_SHGC = 0.763/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_U_WIN = 1.0E+20/XUNIF_U_WIN = 2.716/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GR = 1.0E+20/XUNIF_GR = 0.3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_FLOOR_HEIGHT = 1.0E+20/XUNIF_FLOOR_HEIGHT = 3.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_INF = 1.0E+20/XUNIF_INF = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_V_VENT = 1.0E+20/XUNIF_V_VENT = 0.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HR_TARGET = 1.0E+20/XUNIF_HR_TARGET = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CAP_SYS_RAT = 1.0E+20/XUNIF_CAP_SYS_RAT = 90./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_M_SYS_RAT = 1.0E+20/XUNIF_M_SYS_RAT = 0.0067/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_T_ADP = 1.0E+20/XUNIF_T_ADP = 285.66/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CAP_SYS_HEAT = 1.0E+20/XUNIF_CAP_SYS_HEAT = 100./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_T_SIZE_MAX = 1.0E+20/XUNIF_T_SIZE_MAX = 301.95/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_T_SIZE_MIN = 1.0E+20/XUNIF_T_SIZE_MIN = 268.96/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_SHADE = 1.0E+20/XUNIF_SHADE = 1./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_NATVENT = 1.0E+20/XUNIF_NATVENT = 1./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1

