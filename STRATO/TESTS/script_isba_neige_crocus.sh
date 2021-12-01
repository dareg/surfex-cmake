cp -f $1 OPTIONS.nam_neige
cp -f $1 OPTIONS.nam_save

sed -i "s/XUNIF_VEGTYPE(1) = 1.0E+20/XUNIF_VEGTYPE(1) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(2) = 1.0E+20/XUNIF_VEGTYPE(2) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(3) = 1.0E+20/XUNIF_VEGTYPE(3) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(4) = 1.0E+20/XUNIF_VEGTYPE(4) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(5) = 1.0E+20/XUNIF_VEGTYPE(5) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(6) = 1.0E+20/XUNIF_VEGTYPE(6) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(7) = 1.0E+20/XUNIF_VEGTYPE(7) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(8) = 1.0E+20/XUNIF_VEGTYPE(8) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(9) = 1.0E+20/XUNIF_VEGTYPE(9) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(10) = 1.0E+20/XUNIF_VEGTYPE(10) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(11) = 1.0E+20/XUNIF_VEGTYPE(11) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(12) = 1.0E+20/XUNIF_VEGTYPE(12) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(13) = 1.0E+20/XUNIF_VEGTYPE(13) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(14) = 1.0E+20/XUNIF_VEGTYPE(14) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(15) = 1.0E+20/XUNIF_VEGTYPE(15) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(16) = 1.0E+20/XUNIF_VEGTYPE(16) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(17) = 1.0E+20/XUNIF_VEGTYPE(17) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(18) = 1.0E+20/XUNIF_VEGTYPE(18) = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEGTYPE(19) = 1.0E+20/XUNIF_VEGTYPE(19) = 0./g" OPTIONS.nam_neige

sed -i "s/XUNIF_VEG(1,1) = 1.0E+20/XUNIF_VEG(1,1) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEG(1,2) = 1.0E+20/XUNIF_VEG(1,2) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEG(1,3) = 1.0E+20/XUNIF_VEG(1,3) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEG(1,4) = 1.0E+20/XUNIF_VEG(1,4) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEG(1,5) = 1.0E+20/XUNIF_VEG(1,5) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEG(1,6) = 1.0E+20/XUNIF_VEG(1,6) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEG(1,7) = 1.0E+20/XUNIF_VEG(1,7) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEG(1,8) = 1.0E+20/XUNIF_VEG(1,8) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEG(1,9) = 1.0E+20/XUNIF_VEG(1,9) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEG(1,10) = 1.0E+20/XUNIF_VEG(1,10) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEG(1,11) = 1.0E+20/XUNIF_VEG(1,11) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_VEG(1,12) = 1.0E+20/XUNIF_VEG(1,12) = 1./g" OPTIONS.nam_neige

sed -i "s/XUNIF_LAI(1,1) = 1.0E+20/XUNIF_LAI(1,1) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAI(1,2) = 1.0E+20/XUNIF_LAI(1,2) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAI(1,3) = 1.0E+20/XUNIF_LAI(1,3) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAI(1,4) = 1.0E+20/XUNIF_LAI(1,4) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAI(1,5) = 1.0E+20/XUNIF_LAI(1,5) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAI(1,6) = 1.0E+20/XUNIF_LAI(1,6) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAI(1,7) = 1.0E+20/XUNIF_LAI(1,7) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAI(1,8) = 1.0E+20/XUNIF_LAI(1,8) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAI(1,9) = 1.0E+20/XUNIF_LAI(1,9) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAI(1,10) = 1.0E+20/XUNIF_LAI(1,10) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAI(1,11) = 1.0E+20/XUNIF_LAI(1,11) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAI(1,12) = 1.0E+20/XUNIF_LAI(1,12) = 1./g" OPTIONS.nam_neige

sed -i "s/XUNIF_Z0(1,1) = 1.0E+20/XUNIF_Z0(1,1) = 0.005/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0(1,2) = 1.0E+20/XUNIF_Z0(1,2) = 0.005/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0(1,3) = 1.0E+20/XUNIF_Z0(1,3) = 0.005/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0(1,4) = 1.0E+20/XUNIF_Z0(1,4) = 0.005/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0(1,5) = 1.0E+20/XUNIF_Z0(1,5) = 0.005/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0(1,6) = 1.0E+20/XUNIF_Z0(1,6) = 0.005/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0(1,7) = 1.0E+20/XUNIF_Z0(1,7) = 0.005/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0(1,8) = 1.0E+20/XUNIF_Z0(1,8) = 0.005/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0(1,9) = 1.0E+20/XUNIF_Z0(1,9) = 0.005/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0(1,10) = 1.0E+20/XUNIF_Z0(1,10) = 0.005/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0(1,11) = 1.0E+20/XUNIF_Z0(1,11) = 0.005/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0(1,12) = 1.0E+20/XUNIF_Z0(1,12) = 0.005/g" OPTIONS.nam_neige

sed -i "s/XUNIF_EMIS(1,1) = 1.0E+20/XUNIF_EMIS(1,1) = 0.97/g" OPTIONS.nam_neige
sed -i "s/XUNIF_EMIS(1,2) = 1.0E+20/XUNIF_EMIS(1,2) = 0.97/g" OPTIONS.nam_neige
sed -i "s/XUNIF_EMIS(1,3) = 1.0E+20/XUNIF_EMIS(1,3) = 0.97/g" OPTIONS.nam_neige
sed -i "s/XUNIF_EMIS(1,4) = 1.0E+20/XUNIF_EMIS(1,4) = 0.97/g" OPTIONS.nam_neige
sed -i "s/XUNIF_EMIS(1,5) = 1.0E+20/XUNIF_EMIS(1,5) = 0.97/g" OPTIONS.nam_neige
sed -i "s/XUNIF_EMIS(1,6) = 1.0E+20/XUNIF_EMIS(1,6) = 0.97/g" OPTIONS.nam_neige
sed -i "s/XUNIF_EMIS(1,7) = 1.0E+20/XUNIF_EMIS(1,7) = 0.97/g" OPTIONS.nam_neige
sed -i "s/XUNIF_EMIS(1,8) = 1.0E+20/XUNIF_EMIS(1,8) = 0.97/g" OPTIONS.nam_neige
sed -i "s/XUNIF_EMIS(1,9) = 1.0E+20/XUNIF_EMIS(1,9) = 0.97/g" OPTIONS.nam_neige
sed -i "s/XUNIF_EMIS(1,10) = 1.0E+20/XUNIF_EMIS(1,10) = 0.97/g" OPTIONS.nam_neige
sed -i "s/XUNIF_EMIS(1,11) = 1.0E+20/XUNIF_EMIS(1,11) = 0.97/g" OPTIONS.nam_neige
sed -i "s/XUNIF_EMIS(1,12) = 1.0E+20/XUNIF_EMIS(1,12) = 0.97/g" OPTIONS.nam_neige

sed -i "s/XUNIF_ALBNIR_VEG(1,1) = 1.0E+20/XUNIF_ALBNIR_VEG(1,1) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_VEG(1,2) = 1.0E+20/XUNIF_ALBNIR_VEG(1,2) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_VEG(1,3) = 1.0E+20/XUNIF_ALBNIR_VEG(1,3) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_VEG(1,4) = 1.0E+20/XUNIF_ALBNIR_VEG(1,4) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_VEG(1,5) = 1.0E+20/XUNIF_ALBNIR_VEG(1,5) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_VEG(1,6) = 1.0E+20/XUNIF_ALBNIR_VEG(1,6) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_VEG(1,7) = 1.0E+20/XUNIF_ALBNIR_VEG(1,7) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_VEG(1,8) = 1.0E+20/XUNIF_ALBNIR_VEG(1,8) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_VEG(1,9) = 1.0E+20/XUNIF_ALBNIR_VEG(1,9) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_VEG(1,10) = 1.0E+20/XUNIF_ALBNIR_VEG(1,10) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_VEG(1,11) = 1.0E+20/XUNIF_ALBNIR_VEG(1,11) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_VEG(1,12) = 1.0E+20/XUNIF_ALBNIR_VEG(1,12) = 0.2/g" OPTIONS.nam_neige

sed -i "s/XUNIF_ALBVIS_VEG(1,1) = 1.0E+20/XUNIF_ALBVIS_VEG(1,1) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_VEG(1,2) = 1.0E+20/XUNIF_ALBVIS_VEG(1,2) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_VEG(1,3) = 1.0E+20/XUNIF_ALBVIS_VEG(1,3) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_VEG(1,4) = 1.0E+20/XUNIF_ALBVIS_VEG(1,4) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_VEG(1,5) = 1.0E+20/XUNIF_ALBVIS_VEG(1,5) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_VEG(1,6) = 1.0E+20/XUNIF_ALBVIS_VEG(1,6) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_VEG(1,7) = 1.0E+20/XUNIF_ALBVIS_VEG(1,7) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_VEG(1,8) = 1.0E+20/XUNIF_ALBVIS_VEG(1,8) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_VEG(1,9) = 1.0E+20/XUNIF_ALBVIS_VEG(1,9) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_VEG(1,10) = 1.0E+20/XUNIF_ALBVIS_VEG(1,10) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_VEG(1,11) = 1.0E+20/XUNIF_ALBVIS_VEG(1,11) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_VEG(1,12) = 1.0E+20/XUNIF_ALBVIS_VEG(1,12) = 0.2/g" OPTIONS.nam_neige

sed -i "s/XUNIF_ALBUV_VEG(1,1) = 1.0E+20/XUNIF_ALBUV_VEG(1,1) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_VEG(1,2) = 1.0E+20/XUNIF_ALBUV_VEG(1,2) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_VEG(1,3) = 1.0E+20/XUNIF_ALBUV_VEG(1,3) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_VEG(1,4) = 1.0E+20/XUNIF_ALBUV_VEG(1,4) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_VEG(1,5) = 1.0E+20/XUNIF_ALBUV_VEG(1,5) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_VEG(1,6) = 1.0E+20/XUNIF_ALBUV_VEG(1,6) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_VEG(1,7) = 1.0E+20/XUNIF_ALBUV_VEG(1,7) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_VEG(1,8) = 1.0E+20/XUNIF_ALBUV_VEG(1,8) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_VEG(1,9) = 1.0E+20/XUNIF_ALBUV_VEG(1,9) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_VEG(1,10) = 1.0E+20/XUNIF_ALBUV_VEG(1,10) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_VEG(1,11) = 1.0E+20/XUNIF_ALBUV_VEG(1,11) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_VEG(1,12) = 1.0E+20/XUNIF_ALBUV_VEG(1,12) = 0.2/g" OPTIONS.nam_neige

sed -i "s/XUNIF_ALBNIR_SOIL(1,1) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,1) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_SOIL(1,2) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,2) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_SOIL(1,3) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,3) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_SOIL(1,4) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,4) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_SOIL(1,5) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,5) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_SOIL(1,6) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,6) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_SOIL(1,7) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,7) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_SOIL(1,8) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,8) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_SOIL(1,9) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,9) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_SOIL(1,10) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,10) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_SOIL(1,11) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,11) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBNIR_SOIL(1,12) = 1.0E+20/XUNIF_ALBNIR_SOIL(1,12) = 0.2/g" OPTIONS.nam_neige

sed -i "s/XUNIF_ALBVIS_SOIL(1,1) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,1) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_SOIL(1,2) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,2) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_SOIL(1,3) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,3) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_SOIL(1,4) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,4) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_SOIL(1,5) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,5) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_SOIL(1,6) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,6) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_SOIL(1,7) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,7) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_SOIL(1,8) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,8) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_SOIL(1,9) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,9) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_SOIL(1,10) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,10) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_SOIL(1,11) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,11) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBVIS_SOIL(1,12) = 1.0E+20/XUNIF_ALBVIS_SOIL(1,12) = 0.2/g" OPTIONS.nam_neige

sed -i "s/XUNIF_ALBUV_SOIL(1,1) = 1.0E+20/XUNIF_ALBUV_SOIL(1,1) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_SOIL(1,2) = 1.0E+20/XUNIF_ALBUV_SOIL(1,2) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_SOIL(1,3) = 1.0E+20/XUNIF_ALBUV_SOIL(1,3) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_SOIL(1,4) = 1.0E+20/XUNIF_ALBUV_SOIL(1,4) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_SOIL(1,5) = 1.0E+20/XUNIF_ALBUV_SOIL(1,5) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_SOIL(1,6) = 1.0E+20/XUNIF_ALBUV_SOIL(1,6) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_SOIL(1,7) = 1.0E+20/XUNIF_ALBUV_SOIL(1,7) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_SOIL(1,8) = 1.0E+20/XUNIF_ALBUV_SOIL(1,8) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_SOIL(1,9) = 1.0E+20/XUNIF_ALBUV_SOIL(1,9) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_SOIL(1,10) = 1.0E+20/XUNIF_ALBUV_SOIL(1,10) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_SOIL(1,11) = 1.0E+20/XUNIF_ALBUV_SOIL(1,11) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ALBUV_SOIL(1,12) = 1.0E+20/XUNIF_ALBUV_SOIL(1,12) = 0.2/g" OPTIONS.nam_neige

sed -i "s/XUNIF_GNDLITTER(1,1) = 1.0E+20/XUNIF_GNDLITTER(1,1) = 1.0/g" OPTIONS.nam_neige
sed -i "s/XUNIF_GNDLITTER(1,2) = 1.0E+20/XUNIF_GNDLITTER(1,2) = 1.0/g" OPTIONS.nam_neige
sed -i "s/XUNIF_GNDLITTER(1,3) = 1.0E+20/XUNIF_GNDLITTER(1,3) = 1.0/g" OPTIONS.nam_neige
sed -i "s/XUNIF_GNDLITTER(1,4) = 1.0E+20/XUNIF_GNDLITTER(1,4) = 1.0/g" OPTIONS.nam_neige
sed -i "s/XUNIF_GNDLITTER(1,5) = 1.0E+20/XUNIF_GNDLITTER(1,5) = 1.0/g" OPTIONS.nam_neige
sed -i "s/XUNIF_GNDLITTER(1,6) = 1.0E+20/XUNIF_GNDLITTER(1,6) = 1.0/g" OPTIONS.nam_neige
sed -i "s/XUNIF_GNDLITTER(1,7) = 1.0E+20/XUNIF_GNDLITTER(1,7) = 1.0/g" OPTIONS.nam_neige
sed -i "s/XUNIF_GNDLITTER(1,8) = 1.0E+20/XUNIF_GNDLITTER(1,8) = 1.0/g" OPTIONS.nam_neige
sed -i "s/XUNIF_GNDLITTER(1,9) = 1.0E+20/XUNIF_GNDLITTER(1,9) = 1.0/g" OPTIONS.nam_neige
sed -i "s/XUNIF_GNDLITTER(1,10) = 1.0E+20/XUNIF_GNDLITTER(1,10) = 1.0/g" OPTIONS.nam_neige
sed -i "s/XUNIF_GNDLITTER(1,11) = 1.0E+20/XUNIF_GNDLITTER(1,11) = 1.0/g" OPTIONS.nam_neige
sed -i "s/XUNIF_GNDLITTER(1,12) = 1.0E+20/XUNIF_GNDLITTER(1,12) = 1.0/g" OPTIONS.nam_neige

sed -i "s/XUNIF_Z0LITTER(1,1) = 1.0E+20/XUNIF_Z0LITTER(1,1) = 0.025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0LITTER(1,2) = 1.0E+20/XUNIF_Z0LITTER(1,2) = 0.025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0LITTER(1,3) = 1.0E+20/XUNIF_Z0LITTER(1,3) = 0.025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0LITTER(1,4) = 1.0E+20/XUNIF_Z0LITTER(1,4) = 0.025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0LITTER(1,5) = 1.0E+20/XUNIF_Z0LITTER(1,5) = 0.025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0LITTER(1,6) = 1.0E+20/XUNIF_Z0LITTER(1,6) = 0.025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0LITTER(1,7) = 1.0E+20/XUNIF_Z0LITTER(1,7) = 0.025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0LITTER(1,8) = 1.0E+20/XUNIF_Z0LITTER(1,8) = 0.025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0LITTER(1,9) = 1.0E+20/XUNIF_Z0LITTER(1,9) = 0.025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0LITTER(1,10) = 1.0E+20/XUNIF_Z0LITTER(1,10) = 0.025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0LITTER(1,11) = 1.0E+20/XUNIF_Z0LITTER(1,11) = 0.025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_Z0LITTER(1,12) = 1.0E+20/XUNIF_Z0LITTER(1,12) = 0.025/g" OPTIONS.nam_neige

sed -i "s/XUNIF_H_TREE(1) = 1.0E+20/XUNIF_H_TREE(1) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(2) = 1.0E+20/XUNIF_H_TREE(2) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(3) = 1.0E+20/XUNIF_H_TREE(3) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(4) = 1.0E+20/XUNIF_H_TREE(4) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(5) = 1.0E+20/XUNIF_H_TREE(5) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(6) = 1.0E+20/XUNIF_H_TREE(6) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(7) = 1.0E+20/XUNIF_H_TREE(7) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(8) = 1.0E+20/XUNIF_H_TREE(8) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(9) = 1.0E+20/XUNIF_H_TREE(9) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(10) = 1.0E+20/XUNIF_H_TREE(10) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(11) = 1.0E+20/XUNIF_H_TREE(11) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(12) = 1.0E+20/XUNIF_H_TREE(12) = 20./g" OPTIONS.nam_neige

sed -i "s/XUNIF_Z0_O_Z0H(1) = 1.0E+20/XUNIF_Z0_O_Z0H(1) = 10./g" OPTIONS.nam_neige
sed -i "s/XUNIF_RSMIN(1) = 1.0E+20/XUNIF_RSMIN(1) = 40./g" OPTIONS.nam_neige
sed -i "s/XUNIF_GAMMA(1) = 1.0E+20/XUNIF_GAMMA(1) = 0.04/g" OPTIONS.nam_neige
sed -i "s/XUNIF_WRMAX_CF(1) = 1.0E+20/XUNIF_WRMAX_CF(1) = 0.2/g" OPTIONS.nam_neige
sed -i "s/XUNIF_RGL(1) = 1.0E+20/XUNIF_RGL(1) = 30./g" OPTIONS.nam_neige
sed -i "s/XUNIF_CV(1) = 1.0E+20/XUNIF_CV(1) = 0.00002/g" OPTIONS.nam_neige

sed -i "s/XUNIF_GMES(1) = 1.0E+20/XUNIF_GMES(1) = 0.001/g" OPTIONS.nam_neige
sed -i "s/XUNIF_RE25(1) = 1.0E+20/XUNIF_RE25(1) = 1.8D-7/g" OPTIONS.nam_neige
sed -i "s/XUNIF_BSLAI(1) = 1.0E+20/XUNIF_BSLAI(1) = 0.08/g" OPTIONS.nam_neige
sed -i "s/XUNIF_LAIMIN(1) = 1.0E+20/XUNIF_LAIMIN(1) = 0.3/g" OPTIONS.nam_neige
sed -i "s/XUNIF_SEFOLD(1) = 1.0E+20/XUNIF_SEFOLD(1) = 12960000./g" OPTIONS.nam_neige
sed -i "s/XUNIF_GC(1) = 1.0E+20/XUNIF_GC(1) = 0.00025/g" OPTIONS.nam_neige
sed -i "s/XUNIF_DMAX(1) = 1.0E+20/XUNIF_DMAX(1) = 0.05/g" OPTIONS.nam_neige
sed -i "s/XUNIF_F2I(1) = 1.0E+20/XUNIF_F2I(1) = 0.3/g" OPTIONS.nam_neige
sed -i "s/XUNIF_H_TREE(1) = 1.0E+20/XUNIF_H_TREE(1) = 20./g" OPTIONS.nam_neige
sed -i "s/XUNIF_CE_NITRO(1) = 1.0E+20/XUNIF_CE_NITRO(1) = 5.56/g" OPTIONS.nam_neige
sed -i "s/XUNIF_CF_NITRO(1) = 1.0E+20/XUNIF_CF_NITRO(1) = 6.73/g" OPTIONS.nam_neige
sed -i "s/XUNIF_CNA_NITRO(1) = 1.0E+20/XUNIF_CNA_NITRO(1) = 1.3/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ROOT_DEPTH(1) = 1.0E+20/XUNIF_ROOT_DEPTH(1) = 1./g" OPTIONS.nam_neige
sed -i "s/XUNIF_GROUND_DEPTH(1) = 1.0E+20/XUNIF_GROUND_DEPTH(1) = 10./g" OPTIONS.nam_neige
sed -i "s/XUNIF_ROOT_LIN(1) = 1.0E+20/XUNIF_ROOT_LIN(1) = 0.05/g" OPTIONS.nam_neige
sed -i "s/XUNIF_ROOT_EXTINCTION(1) = 1.0E+20/XUNIF_ROOT_EXTINCTION(1) = 0.943/g" OPTIONS.nam_neige

sed -i "s/LECOCLIMAP = T/LECOCLIMAP = F/g" OPTIONS.nam_neige
sed -i "s/XUNIF_SEA = 1.0E+20/XUNIF_SEA = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_WATER = 1.0E+20/XUNIF_WATER = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_TOWN = 1.0E+20/XUNIF_TOWN = 0./g" OPTIONS.nam_neige
sed -i "s/XUNIF_NATURE = 1.0E+20/XUNIF_NATURE = 1./g" OPTIONS.nam_neige

sed -i "s/NGROUND_LAYER = 8/NGROUND_LAYER = 20/g" OPTIONS.nam_neige
sed -i "s/XSOILGRID = 0.01,0.1,1.1,2.4,3.9,5.2,6.1,7.9/XSOILGRID = 0.01,0.03,0.06,0.1,0.2,0.3,0.45,0.6,0.8,1.,1.25,1.5,2.,2.5,3.,4.,5.,6.5,8.,10./g" OPTIONS.nam_neige

sed -i "s/NYEAR = 1986/NYEAR = 2017/g" OPTIONS.nam_neige
sed -i "s/NMONTH = 1/NMONTH = 2/g" OPTIONS.nam_neige
sed -i "s/NDAY = 29/NDAY = 1/g" OPTIONS.nam_neige
sed -i "s/XTIME = 0./XTIME = 21600./g" OPTIONS.nam_neige
sed -i "s/NTIME = 36/NTIME = 12/g" OPTIONS.nam_neige

sed -i "s/CSCOND = \"NP89\"/CSCOND = \"PL98\"/g" OPTIONS.nam_neige
sed -i "s/XTSTEP_OUTPUT = 10800./XTSTEP_OUTPUT = 3600./g" OPTIONS.nam_neige

sed -i "s/CGRID = \"CONF PROJ\"/CGRID = \"LONLATVAL\"/g" OPTIONS.nam_neige
sed -i "s/NPOINTS = 127/NPOINTS = 21/g" OPTIONS.nam_neige
sed -i "s/XX = -88.23, -82.21, -87.09,-121.30,-114.37, -81.16, -97.25, -77.77, -85.36, -109.96,/XX = 6.64493, 6.64493,/g" OPTIONS.nam_neige
sed -i "s/-64.31, -102.27, -100.05, -70.28, -98.80, -88.55, -94.91, -101.44, -71.13, -100.14, -154.90,/6.64493, 6.64493,/g" OPTIONS.nam_neige
sed -i "s/-103.33, -73.81, -99.40, -95.28, -117.91, -97.66, -80.86, -66.32, -102.27, -112.08, -98.07,/6.64493, 6.64493,/g" OPTIONS.nam_neige
sed -i "s/-104.83, -74.40, -92.97, -106.64, -164.22, -97.73, -156.40, -82.73, -115.49, -132.76, -95.08,/6.64493, 6.64493,/g" OPTIONS.nam_neige
sed -i "s/-98.59, -103.05, -72.98, -72.02,-110.66,-73.27,-104.15,-160.73,-94.21,-94.70,-86.35,-113.17,/6.64493, 6.64493, 6.64493,/g" OPTIONS.nam_neige
sed -i "s/-111.30,-118.44,-60.83,-115.83,-79.92,-108.92,-153.60,-97.75,-108.55,-74.71,-100.04,-108.46,-99.40,-113.84,-97.57,/6.64493, 6.64493,/g" OPTIONS.nam_neige
sed -i "s/-79.42,-89.11,-94.90,-97.58,-110.38,-76.28,-101.31,-83.85,-101.73,-99.48,-77.02,-107.45,-90.02,-107.28,-76.75,/6.64493, 6.64493,/g" OPTIONS.nam_neige
sed -i "s/-90.81,-100.44,-104.68,-80.64,-115.49,-100.16,-133.75,-106.07,-105.27,-117.64,-97.29,-95.79,-108.29,-100.04,/6.64493, 6.64493,/g" OPTIONS.nam_neige
sed -i "s/-103.45,-105.73,-107.82,-70.94,-88.42,-66.14,-83.67,-118.83,-99.77,-73.82,-93.03,-73.42,-100.21,-102.94,/6.64493, 6.64493,/g" OPTIONS.nam_neige
sed -i "s/-93.65,-121.13,-70.99,-98.98,-155.67,-120.04,-97.66,-108.22,-119.55,-71.58,-118.96,-118.71,-122.77,-121.19,/6.64493, 6.64493,/g" OPTIONS.nam_neige

sed -i "s/XY = 47.72, 44.78, 43.86, 65.91, 62.09, 42.25, 52.12, 43.85, 11.57, 59.10, 54.19, 57.19, 52.37, 66.42, 50.99,/XY = 46.176849, 46.176849,/g" OPTIONS.nam_neige
sed -i "s/49.80, 49.38, 63.13, 64.99, 53.33, 59.56, 58.30, 50.82, 60.25, 64.13, 63.33, 24.64, 26.95, 54.38,/46.176849, 46.176849, 46.176849,/g" OPTIONS.nam_neige
sed -i "s/60.34, 58.59, 62.69, 55.11, 56.15, 48.61, 57.47, 60.79, 67.79, 57.85, 42.50, 55.43, 69.10, 48.04,/46.176849, 46.176849,/g" OPTIONS.nam_neige
sed -i "s/64.55, 20.21, 55.05, 48.66, 65.59, 44.45, 60.69, 66.51, 54.62, 53.85, 12.32, 66.28, 63.96, 64.95,/46.176849, 46.176849,/g" OPTIONS.nam_neige
sed -i "s/45.95, 33.30, 46.24, 61.82, 70.59, 54.07, 55.84, 57.34, 53.83, 64.15, 65.95, 65.31, 21.66, 44.47, 15.57,/46.176849, 46.176849, 46.176849,/g" OPTIONS.nam_neige
sed -i "s/62.20, 54.71, 64.54, 56.33, 60.96, 15.35, 62.28, 62.99, 50.97, 63.94, 53.77, 54.76, 53.22, 51.04,/46.176849, 46.176849,/g" OPTIONS.nam_neige
sed -i "s/55.83, 60.00, 28.24, 60.22, 54.05, 59.57, 60.02, 69.41, 60.54, 61.37, 62.96, 55.96, 62.27, 54.78,/46.176849, 46.176849,/g" OPTIONS.nam_neige
sed -i "s/59.05, 63.17, 81.80, 44.02, 52.69, 12.54, 43.34, 51.27, 59.40, 53.00, 56.49, 56.40, 60.64, 46.24,/46.176849, 46.176849,/g" OPTIONS.nam_neige
sed -i "s/60.58, 50.76, 63.76, 58.64, 39.09, 63.93, 56.37, 40.03, 18.49, 38.01, 38.70,39.02,40.26,/46.176849, 46.176849, 46.176849,/g" OPTIONS.nam_neige

sed -i "s/XDX = 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,/XDX = 0.5, 0.5, 0.5, 0.5, 0.5, 0.5/g" OPTIONS.nam_neige
sed -i "s/XDY = 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,/XDY = 0.5, 0.5, 0.5, 0.5, 0.5, 0.5/g" OPTIONS.nam_neige
sed -i "s/0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,/0.5, 0.5, 0.5,/g" OPTIONS.nam_neige

sed -i "s/XTG_SURF = 290.16/XTG_SURF = 273.15/g" OPTIONS.nam_neige
sed -i "s/XTG_ROOT = 285.52/XTG_ROOT = 273.15/g" OPTIONS.nam_neige

sed -i "s/LNOSOF = T/LNOSOF = F/g" OPTIONS.nam_neige

sed -i "s/XWSNOW = 220.//g" OPTIONS.nam_neige
sed -i "s/XZSNOW = 1.0E+20//g" OPTIONS.nam_neige
sed -i "s/XTSNOW = 273.15//g" OPTIONS.nam_neige
sed -i "s/XLWCSNOW = 0.//g" OPTIONS.nam_neige
sed -i "s/XRSNOW = 260.//g" OPTIONS.nam_neige
sed -i "s/XASNOW = 0.95//g" OPTIONS.nam_neige
sed -i "s/XSG1SNOW = 1.2//g" OPTIONS.nam_neige
sed -i "s/XSG2SNOW = 2.2//g" OPTIONS.nam_neige
sed -i "s/XHISTSNOW = 3.//g" OPTIONS.nam_neige
sed -i "s/XAGESNOW = 10.//g" OPTIONS.nam_neige
sed -i "s/XSWEMAX = 500.//g" OPTIONS.nam_neige

echo "############# $2_SNCRO8"
sed -i "s/CFORCING\_FILETYPE = \"ASCII\"/CFORCING\_FILETYPE = \"NETCDF\"/g" OPTIONS.nam_neige
sed -i "s/YZS = \"gtopo30\"/YZS = \"FORCING.nc\"/g" OPTIONS.nam_neige
sed -i "s/YZSFILETYPE = \"DIRECT\"/YZSFILETYPE = \"NETCDF\"/g" OPTIONS.nam_neige
sed -i "s/YSLOPE = \"\"/YSLOPE = \"FORCING.nc\"/g" OPTIONS.nam_neige
sed -i "s/YSLOPEFILETYPE = \"\"/YSLOPEFILETYPE = \"NETCDF\"/g" OPTIONS.nam_neige
sed -i "s/CSNOW = \"D95\"/CSNOW = \"CRO\"/g" OPTIONS.nam_neige
sed -i "s/CSNOW = \"3-L\"/CSNOW = \"CRO\"/g" OPTIONS.nam_neige
sed -i "s/NSNOW\_LAYER = 12/NSNOW\_LAYER = 8/g" OPTIONS.nam_neige
sed -i "s/LGLACIER = F/LGLACIER = T/g" OPTIONS.nam_neige
sed -i "s/LPROSNOW = T/LPROSNOW = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWDIMNC = T/LSNOWDIMNC = F/g" OPTIONS.nam_neige
sed -i "s/LRESETCUMUL = T/LRESETCUMUL = F/g" OPTIONS.nam_neige
sed -i "s/LVOLUMETRIC_SNOWLIQ = T/LVOLUMETRIC_SNOWLIQ = F/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f TESTS/FORCAGES/NETCDF/NEIGE/FORCING_CROCUS.nc FORCING.nc
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO8"
cp -f TESTS/PREP/FILES/NEIGE/PREP_08.nc  PREP_semi_restart.nc
./script_exec_semi_restart.sh "$2_SNCRO8" $3 $4 $5

ln -sf TESTS/ISBA/drdt_bst_fit_60.nc

echo "############# $2_SNCRO50_01"
sed -i "s/NSNOW\_LAYER = 8/NSNOW\_LAYER = 50/g" OPTIONS.nam_neige
sed -i "s/LGLACIER = T/LGLACIER = F/g" OPTIONS.nam_neige
sed -i "s/LPROSNOW = F/LPROSNOW = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWDIMNC = F/LSNOWDIMNC = T/g" OPTIONS.nam_neige
sed -i "s/LRESETCUMUL = F/LRESETCUMUL = T/g" OPTIONS.nam_neige
sed -i "s/LVOLUMETRIC_SNOWLIQ = F/LVOLUMETRIC_SNOWLIQ = T/g" OPTIONS.nam_neige
sed -i "s/CFILE_SNOW = \"\"/CFILE_SNOW = \"neige_20170201.nc\"/g" OPTIONS.nam_neige
sed -i "s/CTYPE_SNOW = \"\"/CTYPE_SNOW = \"NETCDF\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWMETAMO = \"B92\"/CSNOWMETAMO = \"C13\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWFALL = \"V12\"/CSNOWFALL = \"S02\"/g" OPTIONS.nam_neige
sed -i "s/LSNOWSYTRON = F/LSNOWSYTRON = T/g" OPTIONS.nam_neige
sed -i "s/CSNOWRES = \"DEF\"/CSNOWRES = \"RIL\"/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = F/LSLOPE = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWCOMPACT_BOOL = F/LSNOWCOMPACT_BOOL = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWTILLER = F/LSNOWTILLER = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_BOOL = F/LSNOWMAK_BOOL = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_PROP = F/LSNOWMAK_PROP = T/g" OPTIONS.nam_neige
sed -i "s/LSELF_PROD = F/LSELF_PROD = T/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_01"
cp -f TESTS/PREP/FILES/NEIGE/PREP_50.nc  PREP_semi_restart.nc
./script_exec_semi_restart.sh "$2_SNCRO50_01" $3 $4 $5

echo "############# $2_SNCRO50_02"
sed -i "s/CSNOWFALL = \"S02\"/CSNOWFALL = \"A76\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOMP = \"B92\"/CSNOWCOMP = \"T11\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOND = \"Y81\"/CSNOWCOND = \"I02\"/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = T/LSLOPE = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_BOOL = T/LSNOWMAK_BOOL = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_PROP = T/LSNOWMAK_PROP = F/g" OPTIONS.nam_neige
sed -i "s/LSELF_PROD = T/LSELF_PROD = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWDRIFT_SUBLIM = F/LSNOWDRIFT_SUBLIM = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWSYTRON = T/LSNOWSYTRON = F/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_02"
./script_exec_semi_restart.sh "$2_SNCRO50_02" $3 $4 $5

echo "############# $2_SNCRO50_03"
sed -i "s/CSNOWRAD = \"B92\"/CSNOWRAD = \"T17\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOMP = \"T11\"/CSNOWCOMP = \"S14\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOND = \"I02\"/CSNOWCOND = \"C11\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWHOLD = \"B92\"/CSNOWHOLD = \"SPK\"/g" OPTIONS.nam_neige
sed -i "s/LSNOWDRIFT_SUBLIM = T/LSNOWDRIFT_SUBLIM = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWSYTRON = F/LSNOWSYTRON = T/g" OPTIONS.nam_neige
sed -i "s/XCVHEATF = 0.2/XCVHEATF = 0.6/g" OPTIONS.nam_neige
sed -i "s/NIMPUR = 0/NIMPUR = 2/g" OPTIONS.nam_neige
sed -i "s/NIMPUROF = 0/NIMPUROF = 2/g" OPTIONS.nam_neige
sed -i "s/CSPECSNOW = F/CSPECSNOW = T/g" OPTIONS.nam_neige
sed -i "s/LFORCIMP = F/LFORCIMP = T/g" OPTIONS.nam_neige
sed -i "s/LPROBANDS = F/LPROBANDS = T/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = F/LSLOPE = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWCOMPACT_BOOL = T/LSNOWCOMPACT_BOOL = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWTILLER = T/LSNOWTILLER = F/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_03"
#./script_exec_semi_restart.sh "$2_SNCRO50_03" $3 $4 $5

echo "############# $2_SNCRO50_04"
sed -i "s/CSNOWRAD = \"T17\"/CSNOWRAD = \"B92\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWFALL = \"A76\"/CSNOWFALL = \"S02\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOMP = \"S14\"/CSNOWCOMP = \"T11\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOND = \"C11\"/CSNOWCOND = \"Y81\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWHOLD = \"SPK\"/CSNOWHOLD = \"B92\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWRES = \"RIL\"/CSNOWRES = \"M98\"/g" OPTIONS.nam_neige
sed -i "s/NIMPUR = 2/NIMPUR = 0/g" OPTIONS.nam_neige
sed -i "s/NIMPUROF = 2/NIMPUROF = 0/g" OPTIONS.nam_neige
sed -i "s/CSPECSNOW = T/CSPECSNOW = F/g" OPTIONS.nam_neige
sed -i "s/LPROBANDS = T/LPROBANDS = F/g" OPTIONS.nam_neige
#sed -i "s/LWRITE_TOPO = F/LWRITE_TOPO = T/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = T/LSLOPE = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWCOMPACT_BOOL = F/LSNOWCOMPACT_BOOL = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWTILLER = F/LSNOWTILLER = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_BOOL = F/LSNOWMAK_BOOL = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_PROP = F/LSNOWMAK_PROP = T/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_04"
./script_exec_semi_restart.sh "$2_SNCRO50_04" $3 $4 $5

echo "############# $2_SNCRO50_05"
sed -i "s/CSNOWCOMP = \"T11\"/CSNOWCOMP = \"S14\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOND = \"Y81\"/CSNOWCOND = \"I02\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWHOLD = \"B92\"/CSNOWHOLD = \"SPK\"/g" OPTIONS.nam_neige
sed -i "s/LSNOWSYTRON = T/LSNOWSYTRON = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWDRIFT_SUBLIM = F/LSNOWDRIFT_SUBLIM = T/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = F/LSLOPE = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWCOMPACT_BOOL = T/LSNOWCOMPACT_BOOL = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWTILLER = T/LSNOWTILLER = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_BOOL = T/LSNOWMAK_BOOL = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_PROP = T/LSNOWMAK_PROP = F/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_05"
./script_exec_semi_restart.sh "$2_SNCRO50_05" $3 $4 $5

echo "############# $2_SNCRO50_06"
sed -i "s/CSNOWRAD = \"B92\"/CSNOWRAD = \"T17\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWMETAMO = \"C13\"/CSNOWMETAMO = \"F06\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWFALL = \"S02\"/CSNOWFALL = \"V12\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOMP = \"S14\"/CSNOWCOMP = \"T11\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOND = \"I02\"/CSNOWCOND = \"C11\"/g" OPTIONS.nam_neige
sed -i "s/XCVHEATF = 0.6/XCVHEATF = 0.2/g" OPTIONS.nam_neige
sed -i "s/NIMPUR = 0/NIMPUR = 1/g" OPTIONS.nam_neige
sed -i "s/NIMPUROF = 0/NIMPUROF = 1/g" OPTIONS.nam_neige
sed -i "s/CSPECSNOW = F/CSPECSNOW = T/g" OPTIONS.nam_neige
sed -i "s/LPROBANDS = F/LPROBANDS = T/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = T/LSLOPE = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWCOMPACT_BOOL = F/LSNOWCOMPACT_BOOL = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWTILLER = F/LSNOWTILLER = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_BOOL = F/LSNOWMAK_BOOL = T/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_06"
#./script_exec_semi_restart.sh "$2_SNCRO50_06" $3 $4 $5

# OPTIONS FOR MEB
cp -f TESTS/PREP/FILES/NEIGE/OPTIONS.nam_neige_meb OPTIONS.nam_neige
cp -f TESTS/PREP/FILES/NEIGE/PREP_50_MEB.nc  PREP_semi_restart.nc
cp -f OPTIONS.nam_neige OPTIONS.nam

echo "############# $2_SNCRO50_07"
sed -i "s/CSNOWRAD = \"T17\"/CSNOWRAD = \"B92\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWMETAMO = \"F06\"/CSNOWMETAMO = \"T07\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWFALL = \"V12\"/CSNOWFALL = \"S02\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOMP = \"T11\"/CSNOWCOMP = \"S14\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOND = \"C11\"/CSNOWCOND = \"Y81\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWHOLD = \"SPK\"/CSNOWHOLD = \"B02\"/g" OPTIONS.nam_neige
sed -i "s/LSNOWDRIFT_SUBLIM = T/LSNOWDRIFT_SUBLIM = F/g" OPTIONS.nam_neige
sed -i "s/CSNOWRES = \"M98\"/CSNOWRES = \"RIL\"/g" OPTIONS.nam_neige
sed -i "s/NIMPUR = 1/NIMPUR = 0/g" OPTIONS.nam_neige
sed -i "s/NIMPUROF = 1/NIMPUROF = 0/g" OPTIONS.nam_neige
sed -i "s/CSPECSNOW = T/CSPECSNOW = F/g" OPTIONS.nam_neige
sed -i "s/LPROBANDS = T/LPROBANDS = F/g" OPTIONS.nam_neige
#sed -i "s/LWRITE_TOPO = T/LWRITE_TOPO = F/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = F/LSLOPE = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWCOMPACT_BOOL = T/LSNOWCOMPACT_BOOL = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWTILLER = T/LSNOWTILLER = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_BOOL = T/LSNOWMAK_BOOL = F/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_07"
./script_exec_semi_restart.sh "$2_SNCRO50_07" $3 $4 $5

echo "############# $2_SNCRO50_08"
sed -i "s/CSNOWMETAMO = \"T07\"/CSNOWMETAMO = \"C13\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWFALL = \"S02\"/CSNOWFALL = \"V12\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOMP = \"S14\"/CSNOWCOMP = \"T11\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOND = \"Y81\"/CSNOWCOND = \"I02\"/g" OPTIONS.nam_neige
sed -i "s/XCVHEATF = 0.2/XCVHEATF = 0.6/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = T/LSLOPE = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWCOMPACT_BOOL = F/LSNOWCOMPACT_BOOL = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_BOOL = F/LSNOWMAK_BOOL = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_PROP = F/LSNOWMAK_PROP = T/g" OPTIONS.nam_neige
sed -i "s/LSELF_PROD = F/LSELF_PROD = T/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_08"
./script_exec_semi_restart.sh "$2_SNCRO50_08" $3 $4 $5

echo "############# $2_SNCRO50_09"
sed -i "s/CSNOWRAD = \"B92\"/CSNOWRAD = \"T17\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOMP = \"T11\"/CSNOWCOMP = \"B92\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOND = \"I02\"/CSNOWCOND = \"C11\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWHOLD = \"B02\"/CSNOWHOLD = \"SPK\"/g" OPTIONS.nam_neige
sed -i "s/LSNOWSYTRON = F/LSNOWSYTRON = T/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = F/LSLOPE = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_PROP = T/LSNOWMAK_PROP = F/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_09"
#./script_exec_semi_restart.sh "$2_SNCRO50_09" $3 $4 $5

echo "############# $2_SNCRO50_10"
sed -i "s/CSNOWRAD = \"T17\"/CSNOWRAD = \"B92\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWFALL = \"V12\"/CSNOWFALL = \"S02\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOMP = \"B92\"/CSNOWCOMP = \"T11\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOND = \"C11\"/CSNOWCOND = \"Y81\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWHOLD = \"SPK\"/CSNOWHOLD = \"B02\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWRES = \"RIL\"/CSNOWRES = \"DEF\"/g" OPTIONS.nam_neige
#sed -i "s/LWRITE_TOPO = F/LWRITE_TOPO = T/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = T/LSLOPE = F/g" OPTIONS.nam_neige
sed -i "s/LSELF_PROD = T/LSELF_PROD = F/g" OPTIONS.nam_neige
sed -i "s/LMEB_LITTER = T/LMEB_LITTER = F/g" OPTIONS.nam_neige
sed -i "s/LMEB_GNDRES = F/LMEB_GNDRES = T/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_10"
./script_exec_semi_restart.sh "$2_SNCRO50_10" $3 $4 $5

echo "############# $2_SNCRO50_11"
sed -i "s/CSNOWCOMP = \"T11\"/CSNOWCOMP = \"B92\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOND = \"Y81\"/CSNOWCOND = \"I02\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWHOLD = \"B02\"/CSNOWHOLD = \"SPK\"/g" OPTIONS.nam_neige
sed -i "s/LSNOWSYTRON = T/LSNOWSYTRON = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWDRIFT_SUBLIM = F/LSNOWDRIFT_SUBLIM = T/g" OPTIONS.nam_neige
sed -i "s/XCVHEATF = 0.6/XCVHEATF = 0.2/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = F/LSLOPE = T/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_PROP = F/LSNOWMAK_PROP = T/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_11"
./script_exec_semi_restart.sh "$2_SNCRO50_11" $3 $4 $5

echo "############# $2_SNCRO50_12"
sed -i "s/CSNOWRAD = \"B92\"/CSNOWRAD = \"T17\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWFALL = \"S02\"/CSNOWFALL = \"A76\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOMP = \"B92\"/CSNOWCOMP = \"T11\"/g" OPTIONS.nam_neige
sed -i "s/CSNOWCOND = \"I02\"/CSNOWCOND = \"C11\"/g" OPTIONS.nam_neige
sed -i "s/NIMPUR = 0/NIMPUR = 2/g" OPTIONS.nam_neige
sed -i "s/NIMPUROF = 0/NIMPUROF = 2/g" OPTIONS.nam_neige
sed -i "s/CSPECSNOW = F/CSPECSNOW = T/g" OPTIONS.nam_neige
sed -i "s/LPROBANDS = F/LPROBANDS = T/g" OPTIONS.nam_neige
sed -i "s/LSLOPE = T/LSLOPE = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_BOOL = T/LSNOWMAK_BOOL = F/g" OPTIONS.nam_neige
sed -i "s/LSNOWMAK_PROP = T/LSNOWMAK_PROP = F/g" OPTIONS.nam_neige
cp -f OPTIONS.nam_neige OPTIONS.nam
cp -f OPTIONS.nam_neige OPTIONS.nam_"$2_SNCRO50_12"
#./script_exec_semi_restart.sh "$2_SNCRO50_12" $3 $4 $5

cp -f OPTIONS.nam_save OPTIONS.nam

rm -rf FORCING.nc
rm -rf PREP_semi_restart.nc
rm -rf drdt_bst_fit_60.nc
