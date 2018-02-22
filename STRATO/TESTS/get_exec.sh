### to get all cases names
#set -x

#n1
for file in script_pgd*.sh script_prep*.sh script_sea*.sh script_1d*.sh script_water*.sh script_isba*.sh script_teb*.sh script_csts*.sh script_offline*.sh script_cas_parts*.sh script_flake*.sh
do
awk '{if (/script_exec/) {print $0}}' $file > _$file
sed -e '/^cp/d' _$file > __$file 
sed -e 's/\"\ $fname\ $2\ $3//g' __$file > _$file
sed -e 's/\.\/script_exec_omp_pgd.sh\ \"/0-/g' _$file > __$file
sed -e 's/\.\/script_exec_parall.sh\ \"/0-/g' __$file > _$file
sed -e 's/\.\/script_exec_restart_parall.sh\ \"/0-/g' _$file > __$file
sed -e 's/\.\/script_exec_restart.sh\ \"/0-/g' __$file > _$file
sed -e 's/\.\/script_exec.sh\ \"/0-/g' _$file > __$file
sed -e 's/\.\/script_exec_soda_oi.sh\ \"/0-/g' __$file > _$file
sed -e 's/\.\/script_exec_soda_ekf.sh\ \"/0-/g' _$file > __$file
sed -e 's/\"\ $3\ $4\ $5//g' __$file > _$file
sed -e 's/$2//g' _$file > __$file
mv -f __$file _$file
done

## get_exec_01.sh sert au premier niveau d'appel de sous-scripts
## get_exec_02.sh sert aux autres niveaux d'appels de sous-scripts
## doivent être appelés pour chaque chemin, en dissociant les fichiers pères
## quand un fichier père appelle plusieurs scripts fils
## afin de ne pas mélanger les lignes d'appels
######################## SEA ################################
#set -x
./get_exec_01.sh script_sea.sh script_sea_canopy.sh

./get_exec_02.sh script_sea_canopy.sh script_1d_ocean.sh

######################## WATER ##########################

./get_exec_01.sh script_water.sh script_water_canopy.sh

###################### FLAKE ##############################

./get_exec_01.sh script_flake.sh script_flake_canopy.sh

###################### ISBA ###############################
#isba/photo: a
cp -f script_isba.sh script_isba_a.sh
./get_exec_01.sh script_isba_a.sh script_isba_photo.sh
#isba/photo_dif: b
cp -f script_isba.sh script_isba_b.sh
./get_exec_01.sh script_isba_b.sh script_isba_photo_dif.sh
#isba/phys: c
cp -f script_isba.sh script_isba_c.sh
./get_exec_01.sh script_isba_c.sh script_isba_phys.sh

#photo/neige: a 
cp -f script_isba_photo.sh script_isba_photo_a.sh
./get_exec_02.sh script_isba_photo_a.sh script_isba_neige.sh
#photo/canopy_parall: b
cp -f script_isba_photo.sh script_isba_photo_b.sh
./get_exec_02.sh script_isba_photo_b.sh script_isba_canopy_parall.sh

#photo_dif/neige_dif: a
cp -f script_isba_photo_dif.sh script_isba_photo_dif_a.sh
./get_exec_02.sh script_isba_photo_dif_a.sh script_isba_neige_dif.sh
#photo_dif/canopy_parall: b
cp -f script_isba_photo_dif.sh script_isba_photo_dif_b.sh
./get_exec_02.sh script_isba_photo_dif_b.sh script_isba_canopy_parall.sh

./get_exec_02.sh script_isba_neige.sh script_isba_canopy_parall.sh

./get_exec_02.sh script_isba_neige_dif.sh script_isba_canopy_parall.sh

#phys/canopy: a
cp -f script_isba_phys.sh script_isba_phys_a.sh
./get_exec_02.sh script_isba_phys_a.sh script_isba_canopy.sh
#phys/canopy_parall: b
cp -f script_isba_phys.sh script_isba_phys_b.sh
./get_exec_02.sh script_isba_phys_b.sh script_isba_canopy_parall.sh

######################### TEB ##############################

#teb/teb_prep_simple: a
cp -f script_teb.sh script_teb_a.sh
./get_exec_01.sh script_teb_a.sh script_teb_prep_simple.sh
#teb/teb_pgd: b
cp -f script_teb.sh script_teb_b.sh
./get_exec_01.sh script_teb_b.sh script_teb_pgd.sh

./get_exec_02.sh script_teb_pgd.sh script_teb_prep.sh

./get_exec_02.sh script_teb_prep.sh script_teb_run.sh

######################### CSTS ##############################

./get_exec_01.sh script_csts.sh script_csts_canopy.sh

###########################################################

for file in script_pgd*.sh script_prep*.sh script_sea*.sh script_1d*.sh script_water*.sh script_isba*.sh script_teb*.sh script_csts*.sh script_offline*.sh script_cas_parts*.sh
do
sed -e '/^# /d' _$file > __$file
sed -e '/^#/d' __$file > _$file
done

rm -f __*.sh


###########################################################

#dissociation, pour tracer les différents chemins
cp -f _script_isba_a.sh _script_isba_a1.sh
cp -f _script_isba_a.sh _script_isba_a2.sh
cp -f _script_isba_b.sh _script_isba_b1.sh
cp -f _script_isba_b.sh _script_isba_b2.sh
cp -f _script_isba_c.sh _script_isba_c1.sh
cp -f _script_isba_c.sh _script_isba_c2.sh

rm -f _script_pgd.sh

#chemins à 2 étages
./get_exec_4.sh script_water.sh script_water_canopy.sh
./get_exec_4.sh script_flake.sh script_flake_canopy.sh
./get_exec_4.sh script_csts.sh script_csts_canopy.sh
./get_exec_4.sh script_teb_a.sh script_teb_prep_simple.sh
./get_exec_4.sh script_isba_a.sh script_isba_photo.sh
./get_exec_4.sh script_isba_b.sh script_isba_photo_dif.sh
./get_exec_4.sh script_isba_c.sh script_isba_phys.sh

rm -f _script_water_canopy.sh _script_flake_canopy.sh _script_csts_canopy.sh _script_teb_prep_simple.sh

#chemins à 3 étages
./get_exec_4.sh script_sea.sh script_sea_canopy.sh script_1d_ocean.sh

rm -f _script_sea_canopy.sh _script_1d_ocean.sh

./get_exec_4.sh script_isba_a1.sh script_isba_photo_b.sh script_isba_canopy_parall.sh
./get_exec_4.sh script_isba_b1.sh script_isba_photo_dif_b.sh script_isba_canopy_parall.sh
./get_exec_4.sh script_isba_c1.sh script_isba_phys_a.sh script_isba_canopy.sh
./get_exec_4.sh script_isba_c2.sh script_isba_phys_b.sh script_isba_canopy_parall.sh

#chemins à 4 étages
./get_exec_4.sh script_isba_a2.sh script_isba_photo_a.sh script_isba_neige.sh script_isba_canopy_parall.sh
./get_exec_4.sh script_isba_b2.sh script_isba_photo_dif_a.sh script_isba_neige_dif.sh script_isba_canopy_parall.sh

rm -f _script_isba_photo.sh _script_isba_photo_dif.sh _script_isba_phys.sh _script_isba_neige.sh _script_isba_neige_dif.sh \
	_script_isba_canopy.sh _script_isba_canopy_parall.sh 

./get_exec_4.sh script_teb_b.sh script_teb_pgd.sh script_teb_prep.sh script_teb_run.sh

rm -f _script_teb_pgd.sh _script_teb_prep.sh _script_teb_run.sh

#regroupements des scripts pères dissociés
mv -f _script_isba.sh _script_isba_0.sh
mv -f _script_teb.sh _script_teb_0.sh

cat _script_isba_0.sh _script_isba_a.sh _script_isba_a1.sh _script_isba_a2.sh \
	_script_isba_b.sh _script_isba_b1.sh _script_isba_b2.sh \
	_script_isba_c.sh _script_isba_c1.sh _script_isba_c2.sh > _script_isba.sh

cat _script_teb_0.sh _script_teb_a.sh _script_teb_b.sh > _script_teb.sh

#suppression des fichiers intermédiaires
rm -f _*_0.sh
rm -f *_a.sh *_a1.sh *_a2.sh
rm -f *_b.sh *_b1.sh *_b2.sh
rm -f *_c.sh *_c1.sh *_c2.sh
