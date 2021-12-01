dir_run=$1
exec_new=$2
exec_old=$3

fname="INFO_CSTS.txt"

rm -f $fname
touch $fname

cd $1

cp -f TESTS/script_exec.sh .
cp -f TESTS/script_to_old.sh .
cp -f TESTS/script_exec_parall.sh .

rm -f Forc_*.txt
rm -f Params_config.txt

cp -f TESTS/PREP/make*.sh .
cp -f TESTS/PGD/make*.sh .

ln -s TESTS/FORCAGES/HIVER/Params_config.txt .
ln -s TESTS/FORCAGES/HIVER/Forc*.txt .

####################NAMELISTS INITIALES##########################################

#namelists avec PREP uniformes
cp TESTS/OPTIONS.nam_000 OPTIONS.nam_PHYSIO_UNIF
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

#ZS uniforme
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam_save
sed -e "s/XUNIF_ZS = 1.0E+20/XUNIF_ZS = 250./g" OPTIONS.nam_save > OPTIONS.nam_PHYSIO_UNIF

#recup d'un cas ECOCLIMAP 
cp -f OPTIONS.nam_PHYSIO_UNIF OPTIONS.nam
sed -e "s/LECOCLIMAP\ \=\ F/LECOCLIMAP\ \=\ T/g" OPTIONS.nam > OPTIONS.nam_ECOCLIMAP

cp -f OPTIONS.nam_ECOCLIMAP OPTIONS.nam_base

mv -f OPTIONS.nam_base OPTIONS.nam_save
sed -e "s/LSET\_FORC\_ZS\ \=\ F/LSET\_FORC\_ZS\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_csts

##############################################################

#isba : DIF , CSNOW = CROCUS
mv -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/CSNOW\ =\ \"D95\"/CSNOW\ =\ \"EBA\"/g" OPTIONS.nam_save > OPTIONS.nam_csts

#TEB: on active tout
mv -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/NTEB_PATCH\ =\ 1/NTEB_PATCH\ =\ 4/g" OPTIONS.nam_save > OPTIONS.nam_csts
mv -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/CROAD\_DIR =\ \"UNIF\"/CROAD\_DIR\ =\ \"ORIE\"/g" OPTIONS.nam_save > OPTIONS.nam_csts
mv -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/CWALL\_OPT\ =\ \"UNIF\"/CWALL\_OPT\ =\ \"TWO\"/g" OPTIONS.nam_save > OPTIONS.nam_csts
mv -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/CBEM\ =\ \"DEF\"/CBEM\ =\ \"BEM\"/g" OPTIONS.nam_save > OPTIONS.nam_csts
mv -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LGARDEN\ =\ F/LGARDEN\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam_csts
mv -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LGREENROOF\ =\ F/LGREENROOF\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam_csts
mv -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/XUNIF\_GREENROOF\ =\ 1.0E+20/XUNIF\_GREENROOF\ =\ 0.7/g" OPTIONS.nam_save > OPTIONS.nam_csts
./make_greenroof_unif1.sh "OPTIONS.nam_csts"


##############################################################


#on commence par tester NAM_SURF_ATM

#aldthres
echo "########### NAM_SURF_ATM_ALDTHRES"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LALDTHRES\ \=\ F/LALDTHRES\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_ALDTHRES" $fname $2 $3


#Aldz0h
echo "########### NAM_SURF_ATM_ALDZ0H"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LALDZ0H\ \=\ F/LALDZ0H\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_ALDZ0H" $fname $2 $3


#variables qui n'ont pas le même rôle avec ou sans drag_coef_arp
echo "########### NAM_SURF_ATM_VCHRNK_VZ0CM"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/XVCHRNK\ \=\ 0\.015/XVCHRNK\ \=\ 0\.022/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XVZ0CM\ \=\ 0\./XVZ0CM\ \=\ 0\.01/g" OPTIONS.nam_save > OPTIONS.nam
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_VCHRNK_VZ0CM" $fname $2 $3


#drag_coef_arp 
echo "########### NAM_SURF_ATM_DRAG_COEF_ARP"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LDRAG_COEF_ARP\ \=\ F/LDRAG_COEF_ARP\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_drag_coef_arp
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_DRAG_COEF_ARP" $fname $2 $3

echo "########### NAM_SURF_ATM_DRAG_COEF_ARP_VCHRNK_VZ0CM"
cp -f OPTIONS.nam_drag_coef_arp OPTIONS.nam_save
sed -e "s/XVCHRNK\ \=\ 0\.015/XVCHRNK\ \=\ 0\.022/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XVZ0CM\ \=\ 0\./XVZ0CM\ \=\ 0\.01/g" OPTIONS.nam_save > OPTIONS.nam
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_DRAG_COEF_ARP_VCHRNK_VZ0CM" $fname $2 $3

echo "########### NAM_SURF_ATM_DRAG_COEF_ARP_VALEURS"
cp -f OPTIONS.nam_drag_coef_arp OPTIONS.nam_save
sed -e "s/XEDB\ \=\ 5\./XEDB\ \=\ 4\.5/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XEDC\ \=\ 5\./XEDC\ \=\ 4\./g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XEDD\ \=\ 5\./XEDD\ \=\ 5\.5/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XEDK\ \=\ 1\./XEDK\ \=\ 1\.2/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUSURIC\ \=\ 1\./XUSURIC\ \=\ 0\.7/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUSURID\ \=\ 0\.035/XUSURID\ \=\ 0\.020/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUSURICL\ \=\ 4\./XUSURICL\ \=\ 4\.3/g" OPTIONS.nam_save > OPTIONS.nam
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_DRAG_COEF_ARP_VALEURS" $fname $2 $3


#+vziustar0_arp
echo "########### NAM_SURF_ATM_DRAG_COEF_ARP_VALEURS"
cp -f OPTIONS.nam_drag_coef_arp OPTIONS.nam_save
sed -e "s/LVZIUSTAR0_ARP\ \=\ F/LVZIUSTAR0_ARP\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_vziustar0_arp
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_DRAG_COEF_ARP_VZIUSTAR0_ARP" $fname $2 $3

echo "########### NAM_SURF_ATM_DRAG_COEF_ARP_VZIUSTAR0_ARP_0.75"
cp -f OPTIONS.nam_vziustar0_arp OPTIONS.nam_save
sed -e "s/XVZIUSTAR0\ \=\ 0\./XVZIUSTAR0\ \=\ 0\.75/g" OPTIONS.nam_save > OPTIONS.nam
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_DRAG_COEF_ARP_VZIUSTAR0_ARP_0.75" $fname $2 $3

echo "########### NAM_SURF_ATM_DRAG_COEF_ARP_VZIUSTAR0_ARP_0.92_RZHZ0M_1.1"
cp -f OPTIONS.nam_vziustar0_arp OPTIONS.nam_save
sed -e "s/XRZHZ0M\ \=\ 1\./XRZHZ0M\ \=\ 0\.81/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XVZIUSTAR0\ \=\ 0\./XVZIUSTAR0\ \=\ 0\.92/g" OPTIONS.nam_save > OPTIONS.nam
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_DRAG_COEF_ARP_VZIUSTAR0_ARP_0.92_RZHZ0M_1.1" $fname $2 $3


#RRGUST_ARP
echo "########### NAM_SURF_ATM_RRGUST_ARP"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LRRGUST_ARP\ \=\ F/LRRGUST_ARP\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_rrgust
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_RRGUST_ARP" $fname $2 $3


#on bouge les constantes RRGUST_ARP 
echo "########### NAM_SURF_ATM_RRGUST_ARP_VALEURS"
cp -f OPTIONS.nam_rrgust OPTIONS.nam_save
sed -e "s/XRRSCALE\ \=\ 1\.15E\-4/XRRSCALE\ \=\ 1\.82E\-4/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRRGAMMA\ \=\ 0\.8/XRRGAMMA\ \=\ 0\.7/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUTILGUST\ \=\ 0\.125/XUTILGUST\ \=\ 0\.096/g" OPTIONS.nam_save > OPTIONS.nam
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_RRGUST_ARP_VALEURS" $fname $2 $3


#CPL_ARP
echo "########### NAM_SURF_ATM_CPL_ARP"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/CCPSURF\ \=\ \"DRY\"/CCPSURF\ \=\ \"HUM\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LCPL_ARP\ \=\ F/LCPL_ARP\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
#attention, CPL_ARP n'est pas compatible avec GARDEN
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LGARDEN\ \=\ T/LGARDEN\ \=\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LGREENROOF\ \=\ T/LGREENROOF\ \=\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XUNIF\_GREENROOF\ =\ 0.7/XUNIF\_GREENROOF\ =\ 1.0E+20/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_cpl
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_CPL_ARP" $fname $2 $3
echo "########### NAM_SURF_ATM_CPL_ARP_QVNPLUS"
cp -f OPTIONS.nam_cpl OPTIONS.nam_save
sed -e "s/LQVNPLUS\ \=\ F/LQVNPLUS\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_CPL_ARP_QVNPLUS" $fname $2 $3

#LARP_PN : activate ARPEGE PN values for CV and TAU_ICE
echo "########### NAM_SURF_ATM_LARP_PN"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LARP_PN\ \=\ F/LARP_PN\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_LARP
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_LARP_PN" $fname $2 $3

#on teste les autres booléens

echo "########### NAM_REPROD_OPER"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LREPROD_OPER = F/LREPROD_OPER = T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XEVERG_RSMIN = 175./XEVERG_RSMIN = 250./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XEVERG_VEG = 1./XEVERG_VEG = 0.99/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CDGAVG = \"INV\"/CDGAVG = \"ARI\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CDGDIF = \"ROOT\"/CDGDIF = \"SOIL\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CIMPLICIT_WIND = \"NEW\"/CIMPLICIT_WIND = \"OLD\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CQSAT = \"NEW\"/CQSAT = \"OLD\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CCHARNOCK = \"NEW\"/CCHARNOCK = \"OLD\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NAM_REPROD_OPER
./script_csts_canopy.sh "OPTIONS.nam" "NAM_REPROD_OPER" $fname $2 $3

echo "########### NAM_SURF_ATM_VERTSHIFT" 
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LVERTSHIFT\ \=\ T/LVERTSHIFT\ \=\ F/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NAM_SURF_ATM_VERTSHIFT
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_VERTSHIFT" $fname $2 $3

echo "########### NAM_SURF_ATM_NOSOF"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LNOSOF\ \=\ F/LNOSOF\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NAM_SURF_ATM_NOSOF
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_NOSOF" $fname $2 $3

echo "########### NAM_SURF_ATM_LCPL_GCM"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LCPL\_GCM\ \=\ F/LCPL\_GCM\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NAM_SURF_ATM_LCPL_GCM
./script_exec.sh "NAM_SURF_ATM_LCPL_GCM" $fname $2 $3

#on teste les autres constantes
echo "########## NAM_SURF_ATM_CSTS"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/XCISMIN\ \=\ 6\.7E\-5/XCISMIN\ \=\ 6\.2E\-5/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XVMODMIN\ \=\ 0\./XVMODMIN\ \=\ 0\.02/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XDELTA\_MAX\ \=\ 1\./XDELTA\_MAX\ \=\ 0\.94/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRIMAX\ \=\ 0\.2/XRIMAX\ \=\ 0\.32/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NAM_SURF_ATM_CSTS
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_ATM_CSTS" $fname $2 $3



#NAM_SURF_CSTS
echo "########## NAM_SURF_CSTS"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/CSEA\_ALB\ \=\ \"UNIF\"/CSEA\_ALB\ \=\ \"TA96\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LGLACIER\ \=\ F/LGLACIER\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam

cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XEMISSN\ \=\ 0\.99/XEMISSN\ \=\ 1\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XANSMIN\ \=\ 0\.50/XANSMIN\ \=\ 0\.41/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XANSMAX\ \=\ 0\.85/XANSMAX\ \=\ 0\.89/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XAGLAMIN\ \=\ 0\.8/XAGLAMIN\ \=\ 0\.77/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XAGLAMAX\ \=\ 0\.85/XAGLAMAX\ \=\ 0\.87/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XALBWAT\ \=\ 0\.065/XALBWAT\ \=\ 0\.135/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XALBCOEF\_TA96\ \=\ 0\.037/XALBCOEF\_TA96\ \=\ 0\.052/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XALBSCA\_WAT\ \=\ 0\.06/XALBSCA\_WAT\ \=\ 0\.1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XEMISWAT\ \=\ 0\.96/XEMISWAT\ \=\ 0\.98/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XALBWATICE\ \=\ 0\.40/XALBWATICE\ \=\ 0\.85/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XALBSEAICE\ \=\ 0\.85/XALBSEAICE\ \=\ 0\.71/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XALBWATSNOW\ \=\ 0\.60/XALBWATSNOW\ \=\ 0\.85/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XEMISWATICE\ \=\ 0\.97/XEMISWATICE\ \=\ 1\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XHGLA\ \=\ 33\.3/XHGLA\ \=\ 29\.1/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XWSNV\ \=\ 5\.0/XWSNV\ \=\ 4\.2/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XCFFV\ \=\ 3\.0/XCFFV\ \=\ 2\.4/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XZ0SN\ \=\ 0\.001/XZ0SN\ \=\ 0\.003/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XZ0HSN\ \=\ 0\.0001/XZ0HSN\ \=\ 0.0005/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XTAU\_SMELT\ \=\ 300\./XTAU\_SMELT\ \=\ 220\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XZ0FLOOD\ \=\ 0\.0002/XZ0FLOOD\ \=\ 0\.0005/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NAM_SURF_CSTS
./script_csts_canopy.sh "OPTIONS.nam" "NAM_SURF_CSTS" $fname $2 $3



#NAM_CROCUSn

echo "########## NAM_CROCUSN_CSTS"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LGLACIER\ \=\ F/LGLACIER\ \=\ T/g" OPTIONS.nam_save > OPTIONS.nam_crocus

#les constantes
cp -f OPTIONS.nam_crocus OPTIONS.nam_save
sed -e "s/XZ0ICEZ0SNOW\ \=\ 10\./XZ0ICEZ0SNOW\ \=\ 15\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XRHOTHRESOLD\_ICE\ \=\ 850\./XRHOTHRESOLD\_ICE\ \=\ 790\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XALBICE1\ \=\ 0\.38/XALBICE1\ \=\ 0\.52/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XALBICE2\ \=\ 0\.23/XALBICE2\ \=\ 0\.35/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XALBICE3\ \=\ 0\.08/XALBICE3\ \=\ 0\.13/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XVAGING\_NOGLACIER\ \=\ 60\./XVAGING\_NOGLACIER\ \=\ 75\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XVAGING\_GLACIER\ \=\ 900\./XVAGING\_GLACIER\ \=\ 840\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_NAM_CROCUSN_CSTS
./script_exec.sh "NAM_CROCUSN_CSTS" $fname $2 $3


#NAM_SSOn
#Z01D
echo "######### SSO_Z01D"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/CROUGH\ \=\ \"NONE\"/CROUGH\ \=\ \"Z01D\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SSO_Z01D
./script_exec_parall.sh "SSO_Z01D" $fname $2 $3
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XFRACZ0\ \=\ 2\./XFRACZ0\ \=\ 5\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SSO_Z01D_FRACZ0
#Z04D
echo "######### SSO_Z04D"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/CROUGH\ \=\ \"NONE\"/CROUGH\ \=\ \"Z04D\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SSO_Z04D
./script_exec_parall.sh "SSO_Z04D" $fname $2 $3
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XFRACZ0\ \=\ 2\./XFRACZ0\ \=\ 5\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SSO_Z04D_FRACZ0
#BE04
echo "######### SSO_BE04"
cp -f OPTIONS.nam_csts OPTIONS.nam_save
sed -e "s/LISBA\_CANOPY\ =\ F/LISBA\_CANOPY\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LTEB\_CANOPY\ =\ F/LTEB\_CANOPY\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LWAT\_SBL\ =\ F/LWAT\_SBL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/LSEA\_SBL\ =\ F/LSEA\_SBL\ =\ T/g" OPTIONS.nam_save > OPTIONS.nam
mv -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/CROUGH\ \=\ \"NONE\"/CROUGH\ \=\ \"BE04\"/g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SSO_BE04
./script_exec.sh "SSO_BE04" $fname $2 $3
echo "######### SSO_BE04_COEFBE"
cp -f OPTIONS.nam OPTIONS.nam_save
sed -e "s/XCOEFBE\ \=\ 2\./XCOEFBE\ \=\ 5\./g" OPTIONS.nam_save > OPTIONS.nam
cp -f OPTIONS.nam OPTIONS.nam_SSO_BE04_COEFBE
./script_exec_parall.sh "SSO_BE04_COEFBE" $fname $2 $3

for file in TESTS/PREP/make*.sh
do
	file2=$(basename $file)
	rm $file2
done

for file in TESTS/PGD/make*.sh
do
	file2=$(basename $file)
	rm $file2
done

rm -f Params_config.txt .
rm -f Forc*.txt .

