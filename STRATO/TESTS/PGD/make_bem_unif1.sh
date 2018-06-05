#roof
cp -f $1 OPTIONS.nam_save
sed -e "s/NPAR_FLOOR_LAYER = 0/NPAR_FLOOR_LAYER = 1/g" OPTIONS.nam_save > OPTIONS.nam

for i in 1
do
	var1="XUNIF_HC_FLOOR($i)"
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/$var1 = 1.0E+20/$var1 = 800000./g" OPTIONS.nam_save > OPTIONS.nam
done
for i in 1
do
	var1="XUNIF_TC_FLOOR($i)"
	cp -f OPTIONS.nam OPTIONS.nam_save
	sed -e "s/$var1 = 1.0E+20/$var1 = 0.29/g" OPTIONS.nam_save > OPTIONS.nam
done
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_D_FLOOR(1) = 1.0E+20/XUNIF_D_FLOOR(1) = 0.2/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_TCOOL_TARGET = 1.0E+20/XUNIF_TCOOL_TARGET = 297.16/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_THEAT_TARGET = 1.0E+20/XUNIF_THEAT_TARGET = 297.16/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_F_WASTE_CAN = 1.0E+20/XUNIF_F_WASTE_CAN = 0.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_COP_RAT = 1.0E+20/XUNIF_COP_RAT = 2.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_EFF_HEAT = 1.0E+20/XUNIF_EFF_HEAT = 0.9/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_QIN = 1.0E+20/XUNIF_QIN = 4.433/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_QIN_FRAD = 1.0E+20/XUNIF_QIN_FRAD = 0.412/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_QIN_FLAT = 1.0E+20/XUNIF_QIN_FLAT = 0.175/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_SHGC = 1.0E+20/XUNIF_SHGC = 0.425/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_SHGC_SH = 1.0E+20/XUNIF_SHGC_SH = 0.763/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_U_WIN = 1.0E+20/XUNIF_U_WIN = 4.95/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_GR = 1.0E+20/XUNIF_GR = 0.19/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_FLOOR_HEIGHT = 1.0E+20/XUNIF_FLOOR_HEIGHT = 3.07/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_INF = 1.0E+20/XUNIF_INF = 0.52164/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_V_VENT = 1.0E+20/XUNIF_V_VENT = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_HR_TARGET = 1.0E+20/XUNIF_HR_TARGET = 0.5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_F_WATER_COND = 1.0E+20/XUNIF_F_WATER_COND = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CAP_SYS_RAT = 1.0E+20/XUNIF_CAP_SYS_RAT = 100./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_M_SYS_RAT = 1.0E+20/XUNIF_M_SYS_RAT = 0.0067/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_T_ADP = 1.0E+20/XUNIF_T_ADP = 285.66/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_CAP_SYS_HEAT = 1.0E+20/XUNIF_CAP_SYS_HEAT = 90./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_T_SIZE_MAX = 1.0E+20/XUNIF_T_SIZE_MAX = 301.95/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_T_SIZE_MIN = 1.0E+20/XUNIF_T_SIZE_MIN = 268.96/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_SHADE = 1.0E+20/XUNIF_SHADE = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_NATVENT = 1.0E+20/XUNIF_NATVENT = 0./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam $1
