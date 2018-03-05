# 2014-04 Update for ISBA-TOP B. Vincendon
dir_run=$1
exec_new=$2
exec_old=$3

fname="INFO_CAS_PARTS.txt"

rm -f $fname
touch $fname

cd $1

cp -f TESTS/script_prep_exte.sh .
cp -f TESTS/script_offline_exte.sh .
cp -f TESTS/script_exec.sh .
cp -f TESTS/script_exec_parall.sh .
cp -f TESTS/script_to_old.sh .
cp -f TESTS/script_exec_soda_*.sh .

rm -f Forc_*.txt
rm -f Params_config.txt

cp -f TESTS/PREP/make*.sh .
cp -f TESTS/PGD/make*.sh .

ln -s TESTS/FORCAGES/HIVER/Forc*.txt .

cp -f TESTS/FORCAGES/HIVER/Params_config.txt .

####################NAMELISTS INITIALES##########################################

#namelists avec PREP uniformes
cp -f TESTS/OPTIONS.nam_000 OPTIONS.nam_PHYSIO_UNIF
./make_prep_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_isba_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_isba_snow_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_teb_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_gd_gr_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_teb_snow_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
./make_prep_flake_unif.sh "OPTIONS.nam_PHYSIO_UNIF"
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1.0E+20/XSST_UNIF = 283.15/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/XTS_WATER_UNIF = 1.0E+20/XTS_WATER_UNIF = 282./g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF

cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YCLAY = \"\"/YCLAY = \"clay_fao\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YSAND = \"\"/YSAND = \"sand_fao\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YSOC_TOP = \"\"/YSOC_TOP = \"soc_top\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YSOC_SUB = \"\"/YSOC_SUB = \"soc_sub\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YCTI = \"\"/YCTI = \"topo_index\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YPERM = \"\"/YPERM = \"perm_glo_10km\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF

#recup d'un cas ECOCLIMAP 
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam
sed -e "s/LECOCLIMAP\ \=\ F/LECOCLIMAP\ \=\ T/g" OPTIONS.nam > OPTIONS.nam_ECOCLIMAP

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_base

mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/LSET\_FORC\_ZS\ \=\ F/LSET\_FORC\_ZS\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_base

##############################################################
#TEB: on active tout
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/NTEB_PATCH\ =\ 1/NTEB_PATCH\ =\ 4/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CROAD\_DIR =\ \"UNIF\"/CROAD\_DIR\ =\ \"ORIE\"/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CWALL\_OPT\ =\ \"UNIF\"/CWALL\_OPT\ =\ \"TWO\"/g" OPTIONS.nam_save > OPTIONS.nam_base
mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CBEM\ =\ \"DEF\"/CBEM\ =\ \"BEM\"/g" OPTIONS.nam_save > OPTIONS.nam_base
#NAM_DATA_ISBA est incompatible avec l'utilisation de GARDEN-ECOCLIMAP
#mv -f OPTIONS.nam_base OPTIONS.nam_save
#sed -e "s/LGARDEN\ =\ F/LGARDEN\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam_base
#mv -f OPTIONS.nam_base OPTIONS.nam_save
#sed -e "s/LGREENROOF\ =\ F/LGREENROOF\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam_base
#mv -f OPTIONS.nam_base OPTIONS.nam_save
#sed -e "s/XUNIF\_GREENROOF\ =\ 1.0E+20/XUNIF\_GREENROOF\ =\ 0.7/g" OPTIONS.nam_save > OPTIONS.nam_base

##############################################################

cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"LONLAT REG\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLONMIN = -10./XLONMIN = -180./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLONMAX = 80./XLONMAX = 180./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATMIN = 10./XLATMIN = 45./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATMAX = 55./XLATMAX = 46./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NLON = 90/NLON = 360/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NLAT = 45/NLAT = 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRM_COVER = 1.E-6/XRM_COVER = 0.01/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRM_COAST = 1./XRM_COAST = 0.6/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRM_LAKE = 0./XRM_LAKE = 0.50/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRM_SEA = 0./XRM_SEA = 0.20/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LTOWN\_TO\_ROCK\ \=\ F/LTOWN\_TO\_ROCK\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LWATER\_TO\_NATURE\ \=\ F/LWATER\_TO\_NATURE\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/YCOVER = \"ECOCLIMAP_II_EUROP\"/YCOVER = \"ECOCLIMAP_I_GLOBAL\"/g" OPTIONS.nam_save > OPTIONS.nam_landuse
cp -f OPTIONS.nam_landuse OPTIONS.nam_save
sed -e "s/LLAND\_USE\ \=\ F/LLAND\_USE\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_landuse
cp -f OPTIONS.nam_landuse OPTIONS.nam_save
sed -e "s/NHALO\ \=\ 2/NHALO\ \=\ 180/g" OPTIONS.nam_save > OPTIONS.nam_landuse

#NETCDF
echo "########### NAM_LANDUSE_NETCDF"
ln -s TESTS/CAS_PART/LANDUSE/LAND_USE.nc .
cp -f OPTIONS.nam_landuse OPTIONS.nam_save
sed -e "s/CFNAM_VEGTYPE = \"\"/CFNAM_VEGTYPE = \"LAND_USE.nc\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFTYP_VEGTYPE = \"\"/CFTYP_VEGTYPE = \"NETCDF\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LANDUSE_NETCDF
#non compatible avec la parallélisation à cause des interpolations
#./script_exec_all_parall.sh "LANDUSE_NETCDF" $fname $2 $3
rm -f LAND_USE.nc

#formats PREP


cp -f OPTIONS.nam_landuse OPTIONS.nam_save 
sed -e "s/XUNIF_VEGTYPE(1) = 1.0E+20/XUNIF_VEGTYPE(1) = 0.06/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(2) = 1.0E+20/XUNIF_VEGTYPE(2) = 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(3) = 1.0E+20/XUNIF_VEGTYPE(3) = 0.05/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(4) = 1.0E+20/XUNIF_VEGTYPE(4) = 0.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(5) = 1.0E+20/XUNIF_VEGTYPE(5) = 0.07/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(6) = 1.0E+20/XUNIF_VEGTYPE(6) = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(7) = 1.0E+20/XUNIF_VEGTYPE(7) = 0.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(8) = 1.0E+20/XUNIF_VEGTYPE(8) = 0.07/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(9) = 1.0E+20/XUNIF_VEGTYPE(9) = 0.06/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(10) = 1.0E+20/XUNIF_VEGTYPE(10) = 0.24/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(11) = 1.0E+20/XUNIF_VEGTYPE(11) = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(12) = 1.0E+20/XUNIF_VEGTYPE(12) = 0.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(13) = 1.0E+20/XUNIF_VEGTYPE(13) = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(14) = 1.0E+20/XUNIF_VEGTYPE(14) = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(15) = 1.0E+20/XUNIF_VEGTYPE(15) = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(16) = 1.0E+20/XUNIF_VEGTYPE(16) = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(17) = 1.0E+20/XUNIF_VEGTYPE(17) = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(18) = 1.0E+20/XUNIF_VEGTYPE(18) = 0.00/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_VEGTYPE(19) = 1.0E+20/XUNIF_VEGTYPE(19) = 0.00/g" OPTIONS.nam_save > OPTIONS.nam_avant
#
cp -f OPTIONS.nam_avant OPTIONS.nam
#
#./script_prep_exte.sh "NAM_LANDUSE_PREP_NC" $fname $exec_new $exec_old "nc" 
#
cp -f OPTIONS.nam_avant OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ \=\ \"NC\"/CSURF\_FILETYPE\ \=\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
#./script_prep_exte.sh "NAM_LANDUSE_PREP_ASCII" $fname $exec_new $exec_old "txt" 
#c	
cp -f OPTIONS.nam_avant OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ \=\ \"NC\"/CSURF\_FILETYPE\ \=\ \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
#./script_prep_exte.sh "NAM_LANDUSE_PREP_LFI" $fname $exec_new $exec_old "lfi" 
#
#NC
echo "########### NAM_LANDUSE_PREP_NC"
cp -f OPTIONS.nam_landuse OPTIONS.nam_save
sed -e "s/CFNAM_VEGTYPE = \"\"/CFNAM_VEGTYPE = \"PREP_BASE.nc\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFTYP_VEGTYPE = \"\"/CFTYP_VEGTYPE = \"NC\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LANDUSE_PREP_NC
#non compatible avec la parallélisation à cause des interpolations
#./script_exec_parall.sh "NAM_LANDUSE_PREP_NC" $fname $2 $3
rm -f *_BASE*.nc *_BASE*.nc

#ASCII
echo "########### NAM_LANDUSE_PREP_ASCII"
cp -f OPTIONS.nam_landuse OPTIONS.nam_save
sed -e "s/CFNAM_VEGTYPE = \"\"/CFNAM_VEGTYPE = \"PREP_BASE.txt\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFTYP_VEGTYPE = \"\"/CFTYP_VEGTYPE = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LANDUSE_PREP_ASCII
#./script_exec.sh "NAM_LANDUSE_PREP_ASCII" $fname $2 $3
rm -f *_BASE*.txt *_BASE*.txt

#LFI
echo "########### NAM_LANDUSE_PREP_LFI"
cp -f OPTIONS.nam_landuse OPTIONS.nam_save
sed -e "s/CFNAM_VEGTYPE = \"\"/CFNAM_VEGTYPE = \"PREP_BASE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFTYP_VEGTYPE = \"\"/CFTYP_VEGTYPE = \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LANDUSE_PREP_LFI
#./script_exec.sh "NAM_LANDUSE_PREP_LFI" $fname $2 $3
rm -f *_BASE*.lfi *_BASE*.lfi

#######################CHIMIE####################################

#namelists avec PHYSIO non uniformes
cp -f TESTS/OPTIONS.nam_000 OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YCLAY = \"\"/YCLAY = \"clay_fao\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YSAND = \"\"/YSAND = \"sand_fao\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF

#recup d'un cas ECOCLIMAP 
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam
sed -e "s/LECOCLIMAP\ \=\ F/LECOCLIMAP\ \=\ T/g" OPTIONS.nam > OPTIONS.nam_ECOCLIMAP

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_base

cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/LSET\_FORC\_ZS\ \=\ F/LSET\_FORC\_ZS\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_base

#CHIMIE
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CFILE = \"\"/CFILE = \"ecmwf.OD.20040810.18-V2\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILETYPE = \"\"/CFILETYPE = \"GRIB\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLON0 = 3./XLON0 = -80./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRPK = 0./XRPK = 1.0/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATCEN = 45./XLATCEN = 42./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLONCEN = 3./XLONCEN = -86/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NIMAX = 12/NIMAX = 25/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NJMAX = 8/NJMAX = 25/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDX = 60000./XDX = 15000./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDY = 60000./XDY = 15000./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NEMIS_PGD_NBR = 0/NEMIS_PGD_NBR = 61/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CCHEM\_SURF\_FILE\ \=/CCHEM\_SURF\_FILE\ \=\ \"ReLACS\_poet\.nam\"/g" OPTIONS.nam_save > OPTIONS.nam_chimie
ln -s TESTS/CAS_PART/CHIMIE/files/ReLACS_poet.nam .
ln -s TESTS/CAS_PART/CHIMIE/files/EMISSIONS/* .
ln -s TESTS/CAS_PART/CHIMIE/files/ecmwf.OD.20040810.18-V2 .
cp -f TESTS/CAS_PART/CHIMIE/Params_config.txt .

cp -f OPTIONS.nam_chimie OPTIONS.nam_save
sed -e "s/NHALO\ \=\ 2/NHALO\ \=\ 18/g" OPTIONS.nam_save > OPTIONS.nam_chimie
cp -f OPTIONS.nam_chimie OPTIONS.nam_save
sed -e "s/NHALO_PREP\ \=\ 2/NHALO_PREP\ \=\ 5/g" OPTIONS.nam_save > OPTIONS.nam_chimie
cp -f OPTIONS.nam_chimie OPTIONS.nam_save
sed -e "s/NSCAL\ \=\ 0/NSCAL\ \=\ 59/g" OPTIONS.nam_save > OPTIONS.nam_chimie
cp -f OPTIONS.nam_chimie OPTIONS.nam_save
sed -e "s/LCH\_SURF\_EMIS\ \=\ F/LCH\_SURF\_EMIS\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_chimie
cp -f OPTIONS.nam_chimie OPTIONS.nam_save
sed -e "s/LCH\_BIO\_FLUX\ \=\ F/LCH\_BIO\_FLUX\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_chimie
cp -f OPTIONS.nam_chimie OPTIONS.nam_save
sed -e "s/LCH\_NO\_FLUX\ \=\ F/LCH\_NO\_FLUX\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_chimie

#CHIMIE NC NETCDF
echo "############# CHIMIE_NC_NETCDF"
cp -f OPTIONS.nam_chimie OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_CHIMIE_NC_NETCDF
#./script_exec_all_parall.sh "CHIMIE_NC_NETCDF" $fname $2 $3

#CHIMIE ASCII NETCDF
echo "############# CHIMIE_ASCII_NETCDF"
cp -f OPTIONS.nam_chimie OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ \=\ \"NC\"/CSURF\_FILETYPE\ \=\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_CHIMIE_ASCII_NETCDF
#./script_exec.sh "CHIMIE_ASCII_NETCDF" $fname $2 $3

#CHIMIE LFI NETCDF
echo "############# CHIMIE_LFI_NETCDF"
cp -f OPTIONS.nam_chimie OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ \=\ \"NC\"/CSURF\_FILETYPE\ \=\ \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_CHIMIE_LFI_NETCDF
#./script_exec.sh "CHIMIE_LFI_NETCDF" $fname $2 $3

#CHIMIE FA NETCDF
echo "############# CHIMIE_FA_NETCDF"
rm -f *.fa
cp -f OPTIONS.nam_chimie OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ \=\ \"NC\"/CSURF\_FILETYPE\ \=\ \"FA\"/g" OPTIONS.nam_save > OPTIONS.nam
rm -f PREP.fa
cp -f OPTIONS.nam OPTIONS.nam_CHIMIE_FA_NETCDF
#./script_exec.sh "CHIMIE_FA_NETCDF" $fname $2 $3
#
for file in TESTS/CAS_PART/CHIMIE/files/*
do
	file2=$(basename $file)
	rm -f $file2
done
for file in TESTS/CAS_PART/CHIMIE/files/EMISSIONS/*
do
	file2=$(basename $file)
	rm -f $file2
done

#SNAP
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CFILE = \"\"/CFILE = \"arome.AN.20120215.00\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILETYPE = \"\"/CFILETYPE = \"GRIB\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLAT0 = 45./XLAT0 = 49./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLON0 = 3./XLON0 = 2.5/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRPK = 0./XRPK = 1.0/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATCEN = 45./XLATCEN = 49./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLONCEN = 3./XLONCEN = 2.5/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NIMAX = 12/NIMAX = 242/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NJMAX = 8/NJMAX = 250/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDX = 60000./XDX = 2000./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDY = 60000./XDY = 2000./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NEMIS_NBR = 0/NEMIS_NBR = 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CCH_EMIS = \"AGGR\"/CCH_EMIS = \"SNAP\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CCHEM\_SURF\_FILE\ \=/CCHEM\_SURF\_FILE\ \=\ \"ReLACS\_snap\.nam\"/g" OPTIONS.nam_save > OPTIONS.nam_snap

ln -s TESTS/CAS_PART/SNAP/files/* .
cp -f TESTS/CAS_PART/SNAP/Params_config.txt .

cp -f OPTIONS.nam_snap OPTIONS.nam_save
sed -e "s/NSCAL\ \=\ 0/NSCAL\ \=\ 2/g" OPTIONS.nam_save > OPTIONS.nam_snap
cp -f OPTIONS.nam_snap OPTIONS.nam_save
sed -e "s/NHALO\ \=\ 2/NHALO\ \=\ 15/g" OPTIONS.nam_save > OPTIONS.nam_snap
cp -f OPTIONS.nam_snap OPTIONS.nam_save
sed -e "s/NHALO_PREP\ \=\ 2/NHALO_PREP\ \=\ 80/g" OPTIONS.nam_save > OPTIONS.nam_snap
cp -f OPTIONS.nam_snap OPTIONS.nam_save
sed -e "s/LCH\_SURF\_EMIS\ \=\ F/LCH\_SURF\_EMIS\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_snap

#SNAP NC NETCDF
echo "############# SNAP_NC_NETCDF"
cp -f OPTIONS.nam_snap OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SNAP_NC_NETCDF
#./script_exec.sh "SNAP_NC_NETCDF" $fname $2 $3

#SNAP ASCII NETCDF
echo "############# SNAP_ASCII_NETCDF"
cp -f OPTIONS.nam_snap OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ \=\ \"NC\"/CSURF\_FILETYPE\ \=\ \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SNAP_ASCII_NETCDF
#./script_exec.sh "SNAP_ASCII_NETCDF" $fname $2 $3

#SNAP LFI NETCDF
echo "############# SNAP_LFI_NETCDF"
cp -f OPTIONS.nam_snap OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ \=\ \"NC\"/CSURF\_FILETYPE\ \=\ \"LFI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SNAP_LFI_NETCDF
#./script_exec_all_parall.sh "SNAP_LFI_NETCDF" $fname $2 $3

#SNAP FA NETCDF
echo "############# SNAP_FA_NETCDF"
rm -f *.fa
cp -f OPTIONS.nam_snap OPTIONS.nam_save
sed -e "s/CSURF\_FILETYPE\ \=\ \"NC\"/CSURF\_FILETYPE\ \=\ \"FA\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SNAP_FA_NETCDF
#./script_exec.sh "SNAP_FA_NETCDF" $fname $2 $3

for file in TESTS/CAS_PART/SNAP/files/*
do
	file2=$(basename $file)
	rm -f $file2
done

#######################TOPMODEL####################################


#namelists avec PHYSIO non uniformes
cp -f TESTS/OPTIONS.nam_000 OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YCLAY = \"\"/YCLAY = \"clay_fao\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YSAND = \"\"/YSAND = \"sand_fao\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF

#recup d'un cas ECOCLIMAP 
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam
sed -e "s/LECOCLIMAP\ \=\ F/LECOCLIMAP\ \=\ T/g" OPTIONS.nam > OPTIONS.nam_ECOCLIMAP

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_base

cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/CTOWN = \"TEB\"/CTOWN = \"NONE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSEA = \"SEAFLX\"/CSEA = \"NONE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CWATER = \"WATFLX\"/CWATER = \"NONE\"/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CISBA\ \=\ \"2-L\"/CISBA\ \=\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND\_LAYER\ \=\ 2/NGROUND\_LAYER\ \=\ 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CSNOW\ \=\ \"D95\"/CSNOW\ \=\ \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NSNOW\_LAYER\ \=\ 1/NSNOW\_LAYER\ \=\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CKSAT\ \=\ \"DEF\"/CKSAT\ \=\ \"EXP\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CRUNOFF\ \=\ \"WSAT\"/CRUNOFF\ \=\ \"TOPD\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFORCING\_FILETYPE\ \=\ \"ASCII\"/CFORCING\_FILETYPE\ \=\ \"NETCDF\"/g" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSET\_FORC\_ZS\ \=\ F/LSET\_FORC\_ZS\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LCOUPL\_TOPD\ \=\ F/LCOUPL\_TOPD\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTSTEP\ \=\ 300./XTSTEP\ \=\ 900. /g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTSTEP\_SURF\ \=\ 300./XTSTEP\_SURF\ \=\ 900./g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_base_topd


cp -f OPTIONS.nam_base_topd OPTIONS.nam_save
sed -e "s/NYEAR = 1000000000/NYEAR = 2007/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH = 1000000000/NMONTH = 11/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY = 1000000000/NDAY = 19/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTIME = 1.0E+20/XTIME = 0./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_HUG_SURF = \"\"/CFILE_HUG_SURF = \"hug1.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_HUG_ROOT = \"\"/CFILE_HUG_ROOT = \"hug2.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_HUG_DEEP = \"\"/CFILE_HUG_DEEP = \"hug3.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG_SURF = \"\"/CFILE_TG_SURF = \"tg1.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG_ROOT = \"\"/CFILE_TG_ROOT = \"tg2.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_TG_DEEP = \"\"/CFILE_TG_DEEP = \"tg3.dat\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_HUG = \"\"/CTYPE_HUG = \"ASCLLV\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CTYPE_TG = \"\"/CTYPE_TG = \"ASCLLV\"/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLAT0 = 45./XLAT0 = 45.89891944/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLON0 = 3./XLON0 = 2.337229166/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRPK = 0./XRPK = 0.8/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATCEN = 45./XLATCEN = 44.09557/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLONCEN = 3./XLONCEN = 4.14765/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NIMAX = 12/NIMAX = 90/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NJMAX = 8/NJMAX = 60/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDX = 60000./XDX = 1000./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDY = 60000./XDY = 1000./g" OPTIONS.nam_save > OPTIONS.nam_topd

cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/NHALO\ \=\ 2/NHALO\ \=\ 15/g" OPTIONS.nam_save > OPTIONS.nam_topd

cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/XTSTEP\_OUTPUT\ \=\ 10800./XTSTEP\_OUTPUT\ \=\ 3600./g" OPTIONS.nam_save > OPTIONS.nam_topd

rm -f FORCING.nc
ln -s TESTS/CAS_PART/TOPD/files/FORCING.nc FORCING.nc
ln -s TESTS/CAS_PART/TOPD/files/boucoiran*.* .
ln -s TESTS/CAS_PART/TOPD/files/carte_f_dc.txt .
ln -s TESTS/CAS_PART/TOPD/*.dat .

echo "############# TOPD_TOPD_EXP"
cp -f OPTIONS.nam_topd OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_TOPD_TOPD_EXP
#./script_exec.sh "TOPD_TOPD_EXP" $fname $2 $3

echo "############ TOPD_TOPD_EXP_RESTART"

#./script_offline_exte.sh "TOPD_TOPD_EXP_RESTART" $fname $exec_new $exec_old "nc" 

#cas restart topmodel
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/LSTOCK\_TOPD\ \=\ F/LSTOCK\_TOPD\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NNB\_STP\_STOCK\ \=\ 0/NNB\_STP\_STOCK\ \=\ 15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY\ \=\ 19/NDAY\ \=\ 20/g" OPTIONS.nam_save > OPTIONS.nam
cp -f boucoiran_xwtop_sav.map$3 boucoiran_xwtop_init.map$3
rm -f stocks_init.txt*
cp -f stocks_sav.txt$3 stocks_init.txt$3
cp -f surfcont_sav.map$3 surfcont_init.map$3
cp -f boucoiran_xwtop_sav.map$2 boucoiran_xwtop_init.map$2
cp -f stocks_sav.txt$2 stocks_init.txt$2
cp -f surfcont_sav.map$2 surfcont_init.map$2
cp -f OPTIONS.nam OPTIONS.nam_TOPD_TOPD_EXP_RESTART
rm -f FORCING.nc
ln -s TESTS/CAS_PART/TOPD/files/FORCING.nc.RESTART FORCING.nc

#./script_exec_restart.sh "TOPD_TOPD_EXP_RESTART" $fname $2 $3
#rm -f *_OLD*.nc *_NEW*.nc


rm -f FORCING.nc
ln -s TESTS/CAS_PART/TOPD/files/FORCING.nc FORCING.nc

echo "############# TOPD_TOPD_EXP_12P"
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/NPATCH\ \=\ 1/NPATCH\ \=\ 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_TOPD_TOPD_EXP_12P
#./script_exec.sh "TOPD_TOPD_EXP_12P" $fname $2 $3

echo "############# TOPD_DIF_TOPD_DEF"
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/CKSAT\ \=\ \"EXP\"/CKSAT\ \=\ \"DEF\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CISBA\ \=\ \"3-L\"/CISBA\ \=\ \"DIF\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND\_LAYER\ \=\ 3/NGROUND\_LAYER\ \=\ 6/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_TOPD_DIF_TOPD_DEF
#./script_exec.sh "TOPD_DIF_TOPD_DEF" $fname $2 $3


echo "########### TOPD_WSAT_EXP"
#CRUNOFF/=TOPD
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/LSTOCK\_TOPD\ \=\ T/LSTOCK\_TOPD\ \=\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NNB\_STP\_STOCK\ \=\ 15/NNB\_STP\_STOCK\ \=\ 0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY = 20/NDAY = 19/g" OPTIONS.nam_save > OPTIONS.nam


cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/CRUNOFF\ \=\ \"TOPD\"/CRUNOFF\ \=\ \"WSAT\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_TOPD_WSAT_EXP
#./script_exec.sh "TOPD_WSAT_EXP" $fname $2 $3

for file in TESTS/CAS_PART/TOPD/files/*
do
	file2=$(basename $file)
	rm -f $file2
done
for file in TESTS/CAS_PART/TOPD/*
do
	file2=$(basename $file)
	rm -f $file2
done


echo "############# TOPD_TOPD_EXP_IGN"
cp -f OPTIONS.nam_base_topd OPTIONS.nam_save
sed -e "s/NYEAR = 1000000000/NYEAR = 2014/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH = 1000000000/NMONTH = 2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY = 1000000000/NDAY = 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTIME = 1.0E+20/XTIME = 28800./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUG_SURF = 1.0E+20/XHUG_SURF = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUG_ROOT = 1.0E+20/XHUG_ROOT = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUG_DEEP = 1.0E+20/XHUG_DEEP = 0.10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_SURF = 1.0E+20/XTG_SURF = 288.71/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_ROOT = 1.0E+20/XTG_ROOT = 288.71/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_DEEP = 1.0E+20/XTG_DEEP = 288.71/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"IGN\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NPOINTS = 500/NPOINTS = 142/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/620000.,588000.,596000.,604000.,612000.,620000.,628000.,636000.,644000.,652000.,660000.,564000.,572000.,580000.,/\
584000.,592000.,600000.,608000.,616000.,576000.,584000.,592000.,600000.,608000.,616000.,624000.,632000.,640000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/588000.,596000.,604000.,612000.,620000.,628000.,636000.,644000.,652000.,660000.,556000.,564000.,572000.,580000.,/\
544000.,552000.,560000.,576000.,584000.,592000.,600000.,608000.,616000.,624000.,632000.,640000.,528000.,536000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/588000.,596000.,604000.,612000.,620000.,628000.,636000.,644000.,548000.,556000.,564000.,572000.,580000.,588000.,/\
544000.,552000.,560000.,568000.,576000.,584000.,592000.,600000.,608000.,616000.,624000.,632000.,640000.,528000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/596000.,604000.,612000.,620000.,628000.,636000.,644000.,548000.,556000.,564000.,572000.,580000.,588000.,596000.,/\
536000.,544000.,552000.,560000.,568000.,576000.,584000.,592000.,600000.,608000.,616000.,624000.,632000.,640000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/604000.,612000.,620000.,628000.,636000.,644000.,652000.,660000.,548000.,556000.,564000.,572000.,580000.,588000.,/\
648000.,656000.,664000.,672000.,680000.,528000.,536000.,544000.,552000.,560000.,568000.,576000.,584000.,592000.,/g" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/596000.,604000.,612000.,620000.,628000.,636000.,644000.,652000.,660000.,668000.,548000.,556000.,564000.,572000.,/\
600000.,608000.,616000.,624000.,632000.,640000.,648000.,656000.,664000.,672000.,680000.,688000.,696000.,536000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/580000.,588000.,596000.,604000.,612000.,620000.,628000.,636000.,644000.,652000.,660000.,668000.,548000.,556000.,/\
544000.,552000.,560000.,568000.,576000.,584000.,592000.,600000.,608000.,616000.,624000.,632000.,640000.,648000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/564000.,572000.,580000.,588000.,596000.,604000.,612000.,620000.,628000.,636000.,644000.,652000.,660000.,668000.,/\
656000.,664000.,672000.,680000.,688000.,696000.,544000.,552000.,560000.,568000.,576000.,584000.,592000.,600000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/676000.,684000.,692000.,700000.,708000.,716000.,724000.,732000.,548000.,556000.,564000.,572000.,580000.,588000.,/\
608000.,616000.,624000.,632000.,640000.,648000.,656000.,664000.,672000.,680000.,688000.,696000.,552000.,560000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/596000.,604000.,612000.,620000.,628000.,636000.,644000.,652000.,660000.,668000.,676000.,684000.,692000.,700000.,/\
568000.,576000.,584000.,592000.,600000.,608000.,616000.,624000.,632000.,640000.,648000.,656000.,568000.,576000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/708000.,716000.,724000.,732000.,740000.,548000.,556000.,564000.,572000.,580000.,588000.,596000.,604000.,612000.,/\
584000.,592000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2681000.,2673000.,2673000.,2673000.,2673000.,2673000.,2673000.,2673000.,2673000.,2673000.,2673000.,2665000.,2665000.,/\
1576000.,1576000.,1576000.,1576000.,1576000.,1584000.,1584000.,1584000.,1584000.,1584000.,1584000.,1584000.,1584000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2665000.,2665000.,2665000.,2665000.,2665000.,2665000.,2665000.,2665000.,2665000.,2665000.,2665000.,2657000.,2657000.,/\
1584000.,1592000.,1592000.,1592000.,1592000.,1592000.,1592000.,1592000.,1592000.,1592000.,1592000.,1592000.,1592000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2657000.,2657000.,2657000.,2657000.,2657000.,2657000.,2657000.,2657000.,2657000.,2657000.,2649000.,2649000.,2649000.,/\
1600000.,1600000.,1600000.,1600000.,1600000.,1600000.,1600000.,1600000.,1600000.,1600000.,1600000.,1600000.,1600000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2649000.,2649000.,2649000.,2649000.,2649000.,2649000.,2649000.,2649000.,2649000.,2649000.,2641000.,2641000.,2641000.,/\
1600000.,1600000.,1608000.,1608000.,1608000.,1608000.,1608000.,1608000.,1608000.,1608000.,1608000.,1608000.,1608000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2641000.,2641000.,2641000.,2641000.,2641000.,2641000.,2641000.,2641000.,2641000.,2641000.,2641000.,2641000.,2633000.,/\
1608000.,1608000.,1608000.,1608000.,1608000.,1608000.,1608000.,1608000.,1608000.,1616000.,1616000.,1616000.,1616000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2633000.,2633000.,2633000.,2633000.,2633000.,2633000.,2633000.,2633000.,2633000.,2633000.,2633000.,2633000.,2633000.,/\
1616000.,1616000.,1616000.,1616000.,1616000.,1616000.,1616000.,1616000.,1616000.,1616000.,1616000.,1616000.,1616000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2633000.,2633000.,2625000.,2625000.,2625000.,2625000.,2625000.,2625000.,2625000.,2625000.,2625000.,2625000.,2625000.,/\
1616000.,1616000.,1616000.,1616000.,1616000.,1624000.,1624000.,1624000.,1624000.,1624000.,1624000.,1624000.,1624000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2625000.,2625000.,2625000.,2625000.,2625000.,2617000.,2617000.,2617000.,2617000.,2617000.,2617000.,2617000.,2617000.,/\
1624000.,1624000.,1624000.,1624000.,1624000.,1624000.,1624000.,1624000.,1624000.,1624000.,1624000.,1624000.,1624000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2617000.,2617000.,2617000.,2617000.,2617000.,2617000.,2617000.,2617000.,2617000.,2617000.,2617000.,2617000.,2617000.,/\
1632000.,1632000.,1632000.,1632000.,1632000.,1632000.,1632000.,1632000.,1632000.,1632000.,1632000.,1632000.,1632000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2617000.,2617000.,2617000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,/\
1632000.,1632000.,1632000.,1632000.,1632000.,1632000.,1632000.,1640000.,1640000.,1640000.,1640000.,1640000.,1640000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,2609000.,/\
1640000.,1640000.,1640000.,1640000.,1640000.,1640000.,1640000.,1640000.,1648000.,1648000.,1648000.,/" \
OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/2609000.,2609000.,2601000.,2601000.,2601000.,2601000.,2601000.,2601000.,2601000.,2601000.,2601000.,2601000.,2601000.,/\
1648000.,/" \
OPTIONS.nam_save > OPTIONS.nam_topd

cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/CCAT(1)\ \=\ \"boucoiran\"/CCAT(1)\ \=\ \"Gruevo\"/g" OPTIONS.nam_save > OPTIONS.nam_topd
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/CCAT(2)\ \=\ \"\"/CCAT(2)\ \=\ \"Vehtino\"/g" OPTIONS.nam_save > OPTIONS.nam_topd
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/CCAT(3)\ \=\ \"\"/CCAT(3)\ \=\ \"GKula\"/g" OPTIONS.nam_save > OPTIONS.nam_topd
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/XF_PARAM_BV(1)\ \=\ 2.0/XF_PARAM_BV(1)\ \=\ 5.5/g" OPTIONS.nam_save > OPTIONS.nam_topd
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/XC_DEPTH_RATIO_BV(1)\ \=\ 1.00/XC_DEPTH_RATIO_BV(1)\ \=\ 0.5/g" OPTIONS.nam_save > OPTIONS.nam_topd
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/XSPEEDR(1)\ \=\ 1.0/XSPEEDR(1)\ \=\ 3.500/g" OPTIONS.nam_save > OPTIONS.nam_topd
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/XSPEEDG(1)\ \=\ 0.3/XSPEEDG(1)\ \=\ 1.5/g" OPTIONS.nam_save > OPTIONS.nam_topd
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/XQINIT(1)\ \=\ 0.0/XQINIT(1)\ \=\ 0.7/g" OPTIONS.nam_save > OPTIONS.nam_topd
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/LSTOCK_TOPD\ \=\ F/LSTOCK_TOPD\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_topd
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/NNB_STP_STOCK\ \=\ 0/NNB_STP_STOCK\ \=\ 10/g" OPTIONS.nam_save > OPTIONS.nam_topd
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/NNB_STP_RESTART\ \=\ 15/NNB_STP_RESTART\ \=\ 10/g" OPTIONS.nam_save > OPTIONS.nam_topd
cp -f OPTIONS.nam_topd OPTIONS.nam_save
sed -e "s/NNB_TOPD\ \=\ 4/NNB_TOPD\ \=\ 12/g" OPTIONS.nam_save > OPTIONS.nam_topd

rm -f FORCING.nc
ln -s TESTS/CAS_PART/TOPD2/files/FORCING.nc FORCING.nc
ln -s TESTS/CAS_PART/TOPD2/files/GKula*.* .
ln -s TESTS/CAS_PART/TOPD2/files/Gruevo*.* .
ln -s TESTS/CAS_PART/TOPD2/files/Vehtino*.* .
ln -s TESTS/CAS_PART/TOPD2/files/*.map .
rm -f stocks_init.txt stocks_init.txt_new stocks_init.txt_old
ln -s TESTS/CAS_PART/TOPD2/files/stocks_init.txt stocks_init.txt_new
ln -s stocks_init.txt_new stocks_init.txt_old
ln -s stocks_init.txt_new stocks_init.txt

cp -f OPTIONS.nam_topd OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_TOPD_TOPD_EXP_IGN
./script_exec.sh "TOPD_TOPD_EXP_IGN" $fname $2 $3


for file in TESTS/CAS_PART/TOPD2/files/*
do
	file2=$(basename $file)
	rm -f $file2
done

rm -f bilan_*
rm -f boucoiran*
rm -f GKula*
rm -f Gruevo* 
rm -f Vehtino*
rm -f *.map 
rm -f stocks_*

#######################TSZ0####################################

cp -f TESTS/OPTIONS.nam_000 OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/XUNIF_CLAY = 0.33/XUNIF_CLAY = 0.5/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/XUNIF_SAND = 0.33/XUNIF_SAND = 0.5/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF

#recup d'un cas ECOCLIMAP 
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam
sed -e "s/LECOCLIMAP\ \=\ F/LECOCLIMAP\ \=\ T/g" OPTIONS.nam > OPTIONS.nam_ECOCLIMAP

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_save
sed -e "s/XUNIF_COVER(4) = 0./XUNIF_COVER(4) = 1./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF_ZS = 1.0E+20/XUNIF_ZS = 0./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NYEAR = 1000000000/NYEAR = 2007/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH = 1000000000/NMONTH = 9/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY = 1000000000/NDAY = 18/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTIME = 1.0E+20/XTIME = 0./g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUG_SURF = 1.0E+20/XHUG_SURF = -10./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUG_ROOT = 1.0E+20/XHUG_ROOT = -10./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUG_DEEP = 1.0E+20/XHUG_DEEP = -10./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_SURF = 1.0E+20/XHUGI_SURF = 0.001/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_ROOT = 1.0E+20/XHUGI_ROOT = 0.001/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHUGI_DEEP = 1.0E+20/XHUGI_DEEP = 0.001/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_SURF = 1.0E+20/XTG_SURF = 265.9952087402344/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_ROOT = 1.0E+20/XTG_ROOT = 265.9952087402344/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTG_DEEP = 1.0E+20/XTG_DEEP = 265.9952087402344/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"CARTESIAN\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLAT0 = 0./XLAT0 = 73./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NIMAX = 20/NIMAX = 4/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NJMAX = 20/NJMAX = 4/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDX = 250000./XDX = 1500./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDY = 250000./XDY = 1500./g" OPTIONS.nam_save > OPTIONS.nam_tsz0


cp -f OPTIONS.nam_tsz0 OPTIONS.nam_save
sed -e "s/CNATURE\ \=\ \"ISBA\"/CNATURE\ \=\ \"TSZ0\"/g" OPTIONS.nam_save > OPTIONS.nam_tsz0
cp -f OPTIONS.nam_tsz0 OPTIONS.nam_save
sed -e "s/CSEA\ \=\ \"SEAFLX\"/CSEA\ \=\ \"NONE\"/g" OPTIONS.nam_save > OPTIONS.nam_tsz0
cp -f OPTIONS.nam_tsz0 OPTIONS.nam_save
sed -e "s/CWATER\ \=\ \"WATFLX\"/CWATER\ \=\ \"NONE\"/g" OPTIONS.nam_save > OPTIONS.nam_tsz0
cp -f OPTIONS.nam_tsz0 OPTIONS.nam_save
sed -e "s/CTOWN\ \=\ \"TEB\"/CTOWN\ \=\ \"NONE\"/g" OPTIONS.nam_save > OPTIONS.nam_tsz0

echo "########## TSZ0"
cp -f OPTIONS.nam_tsz0 OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_TSZ0
rm -f Forc_*.txt
rm -f Params_config.txt
ln -s TESTS/CAS_PART/TSZ0/Forc_*.txt .
cp -f TESTS/CAS_PART/TSZ0/Params_config.txt .
#./script_exec_all_parall.sh "TSZ0" $fname $2 $3


#######################IDEAL####################################


cp -f OPTIONS.nam_TSZ0 OPTIONS.nam_save
sed -e "s/NIMAX = 1/NIMAX = 5/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NJMAX = 1/NJMAX = 5/g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDX = 1500./XDX = 250000./g" OPTIONS.nam_save > OPTIONS.nam 
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDY = 1500./XDY = 250000./g" OPTIONS.nam_save > OPTIONS.nam_ideal

cp -f OPTIONS.nam_ideal OPTIONS.nam_save
sed -e "s/CNATURE\ \=\ \"ISBA\"/CNATURE\ \=\ \"FLUX\"/g" OPTIONS.nam_save > OPTIONS.nam_ideal
cp -f OPTIONS.nam_ideal OPTIONS.nam_save
sed -e "s/CSEA\ \=\ \"SEAFLX\"/CSEA\ \=\ \"FLUX\"/g" OPTIONS.nam_save > OPTIONS.nam_ideal
cp -f OPTIONS.nam_ideal OPTIONS.nam_save
sed -e "s/CWATER\ \=\ \"WATFLX\"/CWATER\ \=\ \"FLUX\"/g" OPTIONS.nam_save > OPTIONS.nam_ideal
cp -f OPTIONS.nam_ideal OPTIONS.nam_save
sed -e "s/CTOWN\ \=\ \"TEB\"/CTOWN\ \=\ \"FLUX\"/g" OPTIONS.nam_save > OPTIONS.nam_ideal

echo "########## IDEAL"
cp -f OPTIONS.nam_ideal OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_IDEAL
rm -f Forc_*.txt
rm -f Params_config.txt
ln -s TESTS/CAS_PART/IDEAL/Forc_*.txt .
cp -f TESTS/CAS_PART/IDEAL/Params_config.txt .
#./script_exec_all_parall.sh "IDEAL" $fname $2 $3



#######################SODA####################################

#namelists avec PHYSIO non uniformes
cp -f TESTS/OPTIONS.nam_000 OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YCLAY = \"\"/YCLAY = \"clay_fao\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/YSAND = \"\"/YSAND = \"sand_fao\"/g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF

#recup d'un cas ECOCLIMAP 
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam
sed -e "s/LECOCLIMAP\ \=\ F/LECOCLIMAP\ \=\ T/g" OPTIONS.nam > OPTIONS.nam_ECOCLIMAP

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_base
#PGD
cp -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/XLAT0 = 45./XLAT0 = 59.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLON0 = 3./XLON0 = -10.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRPK = 1.0/XRPK = 0.858959896930664/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NIMAX = 12/NIMAX = 39/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NJMAX = 8/NJMAX = 39/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDX = 60000./XDX = 11000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDY = 60000./XDY = 11000./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLATCEN = 45./XLATCEN = 59.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XLONCEN = 3./XLONCEN = 14.0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CISBA = \"2-L\"/CISBA = \"3-L\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NGROUND_LAYER = 2/NGROUND_LAYER = 3/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LRESTART_2M = F/LRESTART_2M = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTSTEP_OUTPUT = 10800./XTSTEP_OUTPUT = 3600./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NHALO_PREP = 2/NHALO_PREP = 80/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_assim1

#links
ln -s TESTS/CAS_PART/SODA/PREP/* .
rm -f Forc_*.txt
ln -s TESTS/CAS_PART/SODA/FORC/Forc_* .
rm -f Params_config.txt
cp -f TESTS/CAS_PART/SODA/FORC/Params_config.txt .
ln -s TESTS/CAS_PART/SODA/ASSIM/* .

ln -s TESTS/CAS_PART/SODA/OI/* .

#OI1
echo "############ SODA_OI1"
cp -f OPTIONS.nam_assim1 OPTIONS.nam_save
sed -e "s/LOBSWG = T/LOBSWG = F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE = \"\"/CFILE = \"ECMWF_GRIB_FILE\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILETYPE = \"\"/CFILETYPE = \"GRIB\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SODA_OI1
cp -f OPTIONS.nam_SODA_OI1 OPTIONS.nam
#./script_exec_soda_oi.sh "SODA_OI1" $fname $2 $3

#OI2
cp -f OPTIONS.nam_SODA_OI1 OPTIONS.nam_save
sed -e "s/LOBSWG = T/LOBSWG = F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LOBS2M = F/LOBS2M = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LIMVEG = T/LIMVEG = F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LEXTRAP_NATURE = F/LEXTRAP_NATURE = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LREAD_SST_FROM_FILE = F/LREAD_SST_FROM_FILE = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LAESNM = F/LAESNM = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LECSST = F/LECSST = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LWATERTG2 = F/LWATERTG2 = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NPRINTLEV = 0/NPRINTLEV = 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_assim2

echo "############ SODA_OI2"
cp -f OPTIONS.nam_assim2 OPTIONS.nam_SODA_OI2
cp -f OPTIONS.nam_SODA_OI2 OPTIONS.nam
#./script_exec_soda_oi.sh "SODA_OI2" $fname $2 $3

#EKF
cp -f OPTIONS.nam_assim1 OPTIONS.nam_save
sed -e "s/CASSIM_ISBA = \"OI   \"/CASSIM_ISBA = \"EKF  \"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSET_FORC_ZS = F/LSET_FORC_ZS = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NYEAR = 1000000000/NYEAR = 2007/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NMONTH = 1000000000/NMONTH = 7/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NDAY = 1000000000/NDAY = 10/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTIME = 1.0E+20/XTIME = 21600./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XSST_UNIF = 1.0E+20/XSST_UNIF = 283.15/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTS_WATER_UNIF = 1.0E+20/XTS_WATER_UNIF = 282./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_assim3
./make_prep_unif.sh "OPTIONS.nam_assim3"
./make_prep_isba_unif.sh "OPTIONS.nam_assim3"
./make_prep_isba_snow_unif.sh "OPTIONS.nam_assim3"
./make_prep_teb_unif.sh "OPTIONS.nam_assim3"
./make_prep_gd_gr_unif.sh "OPTIONS.nam_assim3"
./make_prep_teb_snow_unif.sh "OPTIONS.nam_assim3"
./make_prep_gd_gr_snow_unif.sh "OPTIONS.nam_assim3"
./make_prep_flake_unif.sh "OPTIONS.nam_assim3"

echo "############ SODA_EKF1"
cp -f OPTIONS.nam_assim3 OPTIONS.nam_save
sed -e "s/CGRID = \"CONF PROJ\"/CGRID = \"IGN\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NPOINTS = 500/NPOINTS = 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NPATCH = 1/NPATCH = 12/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CPHOTO = \"NON\"/CPHOTO = \"NIT\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFORCING_FILETYPE = \"ASCII\"/CFORCING_FILETYPE = \"NETCDF\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LEXTRAP_SEA = T/LEXTRAP_SEA = F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LEXTRAP_WATER = T/LEXTRAP_WATER = F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_FORMAT_OBS = \"FA   \"/CFILE_FORMAT_OBS = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LOBSNAT = F/LOBSNAT = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SODA_EKF1
rm -f BGROUND*
./script_exec_soda_ekf.sh "SODA_EKF1" $fname $2 $3

echo "############ SODA_EKF2"
cp -f OPTIONS.nam_assim3 OPTIONS.nam_save
sed -e "s/LBEV = T/LBEV = F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LBFIXED = F/LBFIXED = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LOBSFILE = T/LOBSFILE = F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CBIO = \"BIOMA1\"/CBIO = \"LAI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NNCO = 0, 0, 1, 1/NNCO = 1, 1, 0, 0/g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_save
sed -e "s/LEXTRAP_SEA = T/LEXTRAP_SEA = F/g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_save
sed -e "s/LEXTRAP_WATER = T/LEXTRAP_WATER = F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SODA_EKF2
rm -f BGROUND*
#./script_exec_soda_ekf.sh "SODA_EKF2" $fname $2 $3

#ENKF
echo "############ SODA_ENKF1"
cp -f OPTIONS.nam_SODA_EKF1 OPTIONS.nam_save
sed -e "s/CASSIM_ISBA = \"EKF  \"/CASSIM_ISBA = \"ENKF \"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NVAR = 2/NVAR = 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NNCV = 1, 0, 1, 0, 0/NNCV = 1, 1, 0, 1, 0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NOBSTYPE = 2/NOBSTYPE = 1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NNCO = 0, 0, 1, 1/NNCO = 0, 0, 1, 0/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/NENS = 1/NENS = 3/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XADDINFL_M = 0., 0., 0., 0., 0./XADDINFL_M = 0.3, 0.5, 0., 0., 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XINFL_M = 0., 0., 0., 0., 0./XINFL_M = 0.4, 0.7, 0., 1., 0./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XADDTIMECORR_M = 0., 0., 0., 0., 0./XADDTIMECORR_M = 3., 8., 0., 0., 0./g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_save
sed -e "s/CFILE_FORMAT_OBS = \"FA   \"/CFILE_FORMAT_OBS = \"ASCII\"/g" OPTIONS.nam_save > OPTIONS.nam
cp OPTIONS.nam OPTIONS.nam_save
sed -e "s/LOBSNAT = F/LOBSNAT = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SODA_ENKF1
./script_exec_soda_enkf.sh "SODA_ENKF1" $fname $2 $3 3


echo "############ SODA_ENKF2"
cp -f OPTIONS.nam_SODA_ENKF1 OPTIONS.nam_save
sed -e "s/LPB_CORRELATIONS = F/LPB_CORRELATIONS = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LBIAS_CORRECTION = F/LBIAS_CORRECTION = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LENKF = T/LENKF = F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LDENKF = F/LDENKF = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SODA_ENKF2
./script_exec_soda_enkf.sh "SODA_ENKF2" $fname $2 $3 4


echo "############ SODA_ENKF3"
cp -f OPTIONS.nam_SODA_ENKF1 OPTIONS.nam_save
sed -e "s/LPERTURBATION_RUN = F/LPERTURBATION_RUN = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SODA_ENKF3
./script_exec_soda_enkf.sh "SODA_ENKF3" $fname $2 $3 3


for file in `ls TESTS/CAS_PART/SODA/PREP`
do
	rm -f $file
done
for file in `ls TESTS/CAS_PART/SODA/ASSIM`
do
	rm -f $file
done
for file in `ls TESTS/CAS_PART/SODA/OI`
do
	rm -f $file
done

