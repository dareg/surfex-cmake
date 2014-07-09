#!/bin/ksh
########################################################################################
#                                                                                      #
#   Extended Kalman Filter for Land Data Assimilation                                  #
#   ISBA land surface scheme (SURFEX version)                                          #
#   Analysis of Ts, T2, wg and w2 (chosen in namelist OPTIONS.nam)                     #
#   Observations : screen-level parameters T2m et HU2m  + Superficial soil moisture    #
#                                                                                      #
#   Jean-Francois Mahfouf (20 November 2008)                                           #
#                                                                                      #
#   C. Rüdiger: adaptation for running in ASCII environment                            #
#   A. Barbu : Analysis of w2 and LAI                                                  #
#              Observations wg (swi2) and LAI                                          #
########################################################################################
#
# Initialisation forcing, observations, and execution directories
#
#set -e
#set -x
echo 'Initialisation ...'
#
expid=test_varassim
exp_options=V8
exp_prep=V8
exp_pgd=V8
#
# change to local path
SURFEX_SVN=$HOME/SVN/surfex/branches/V8_FOR_UNIFIED_ASSIMILATION
WORKDIR_svn=$SURFEX_SVN/MY_RUN/KTEST/TEST-VAR
#
repforcing=$WORKDIR_svn/forcing/
represults=$WORKDIR_svn/results/${expid}/
reprun=$WORKDIR_svn/rundir/work
repobs=$WORKDIR_svn/data/                
repname=$WORKDIR_svn/namelist/
repinput=$WORKDIR_svn/input

mkdir -p $represults
#----------------------------------
#  Cleaning of useless files
#---------------------------------
rm -f $reprun/*

ln -s $SURFEX_SVN/MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.bin $reprun/ecoclimapI_covers_param.bin
ln -s $SURFEX_SVN/MY_RUN/ECOCLIMAP/ecoclimapII_eu_covers_param.bin $reprun/ecoclimapII_eu_covers_param.bin
ln -s $SURFEX_SVN/MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.dat $reprun/ecoclimapI_covers_param.dat
ln -s $SURFEX_SVN/MY_RUN/ECOCLIMAP/ecoclimapII_eu_covers_param.dat $reprun/ecoclimapII_eu_covers_param.dat

repbin=/home/faroux/SVN/surfex/branches/V8_FOR_UNIFIED_ASSIMILATION/exe

ln -s    $repbin/VARASSIM-LXgfortran-SFX-V7-3-0-MYSRC-MPIAUTO-DEBUG  $reprun/VARASSIM
ln -s    $repbin/OFFLINE-LXgfortran-SFX-V7-3-0-MYSRC-MPIAUTO-DEBUG  $reprun/OFFLINE

#
cp -f $repname'OPTIONS_'${exp_options}.nam $reprun/OPTIONS_init.nam
cp -f $repname'OPTIONS_'${exp_options}.nam $represults/OPTIONS_init.nam
#
cd $reprun

# set the number of control variables (NVAR) in namelist
NVAR=2
#
sed "s/NVAR         = 1/NVAR         = $NVAR/"  OPTIONS_init.nam > OPTIONS_init.nam_
mv -f OPTIONS_init.nam_ OPTIONS_init.nam

cp -f $repinput/'PGD_'${exp_pgd}.txt $reprun/PGD.txt
cp -f $repinput/'PREP_'${exp_prep}.txt $reprun/PREP_REF.txt
#---------------------------------------------------
# Define the experiment identifier and beginning date
#----------------------------------------------------
AAAAMMJJRR_start=2007071006
AAAAMMJJRR_end=2007071506

AAAAMMJJRR=$AAAAMMJJRR_start

typeset -Z2 RR
typeset -Z2 RRobs

echo 'Running EKF from '${AAAAMMJJRR_start}' to '${AAAAMMJJRR_end}
#----------------------------------
# Loop over assimilation windows
#-------------------------------------
while [ $AAAAMMJJRR  -le $AAAAMMJJRR_end ] ; do
  echo 'Starting period: ' $AAAAMMJJRR

  cp -f OPTIONS_init.nam OPTIONS.nam

  aa=`echo $AAAAMMJJRR | cut -c3-4`
  mm=`echo $AAAAMMJJRR | cut -c5-6`
  jj=`echo $AAAAMMJJRR | cut -c7-8`
  RR=`echo $AAAAMMJJRR | cut -c9-10`

  case $RR in
#   00) NT=00;;
    06) NT=06;;
#   12) NT=12;;
#   18) NT=18;;
    *) echo "Error in date $1 !!"; exit 1;;
  esac

  AAAAMMJJRRobs=`$repinput/smsdate $AAAAMMJJRR 24 `
  aaobs=`echo $AAAAMMJJRRobs | cut -c3-4`
  mmobs=`echo $AAAAMMJJRRobs | cut -c5-6`
  jjobs=`echo $AAAAMMJJRRobs | cut -c7-8`
  RRobs=`echo $AAAAMMJJRRobs | cut -c9-10`

#-------------------------------------------
# Get atmospheric forcing data 
#--------------------------------------------
  cp -f $repforcing'20'$aa$mm${jj}${NT}'_FORCING_7_3.nc' $reprun/'FORCING.nc'

# loop on control variables`
 echo 'Generating perturbations ...'

  vv=1 
  while [ $vv -le $NVAR ] ; do
 
#------------------------------------------------
# VARASSIM - define perturbed initial conditions
#--------------------------------------------------
    cp -f PREP_REF.txt PREP.txt
	sed "s/LPRT        =  F/LPRT        =  T/" OPTIONS.nam > OPTIONS.nam_
	mv -f OPTIONS.nam_ OPTIONS.nam
#------------------------------------------------
# Reset variable name to be read/perturbed/written
#-------------------------------------------------
	sed "s/NIVAR         = 1/NIVAR         = $vv/" OPTIONS.nam > OPTIONS.nam_
	mv -f OPTIONS.nam_ OPTIONS.nam
	
	echo 'define perturbation' $vv
	rm -f stderr

	./VARASSIM 2>stderr
cp -f SURFOUT.txt SURFOUT$vv.txt
    mv -f PERTURB                PERTURB_${vv}_$aaobs$mmobs${jjobs}'H'$RRobs.DAT

    if [ $vv -eq 1 ] ; then     
         sed -i "s/WG2/oldWG2/" PREP.txt 
    else 
         sed -i "s/LAI/oldLAI/" PREP.txt 
         sed -i "s/BIOMA1/oldBIOMA1/" PREP.txt
    fi

    cat PREP.txt >> SURFOUT.txt                   # changed for ASCII treatment of PREP
    mv -f SURFOUT.txt PREP.txt                    # changed for ASCII treatment of PREP
    cp -f PREP.txt PREP$vv.txt
	mv -f stderr stderr1$vv
	rm -f fort.10
   
    sed "s/NIVAR         = $vv/NIVAR         = 1/" OPTIONS.nam > OPTIONS.nam_
	mv -f OPTIONS.nam_ OPTIONS.nam
	sed "s/LPRT        =  T/LPRT        =  F/" OPTIONS.nam > OPTIONS.nam_
	mv -f OPTIONS.nam_ OPTIONS.nam

#-------------------------------------------------
#  SURFEX run with perturbed initial conditions
#------------------------------------------------
	echo 'SURFEX run with perturbation' $vv

	./OFFLINE 2>stderr 
         
	mv -f stderr stderr2$vv
	rm -f *.tex LISTING.txt *TXT
        mv SURFOUT.txt PREP.txt

#------------------------------------------------------
# VARASSIM - storage of perturbed simulated observations
#------------------------------------------------------- 
	sed "s/LSIM        =  F/LSIM        =  T/" OPTIONS.nam > OPTIONS.nam_
	mv -f OPTIONS.nam_ OPTIONS.nam
	echo 'store perturbed run' $vv

	./VARASSIM 2>stderr_assim

	mv -f stderr_assim stderr3$vv
	sed "s/LSIM        =  T/LSIM        =  F/" OPTIONS.nam > OPTIONS.nam_
	mv -f OPTIONS.nam_ OPTIONS.nam

	
    mv -f OBSIMU               OBSIMU_PERT_${vv}_$aaobs$mmobs${jjobs}'H'$RRobs.DAT
    mv -f MDSIMU               MDSIMU_PERT_${vv}_$aaobs$mmobs${jjobs}'H'$RRobs.DAT

	rm -f fort.10

	(( vv = vv + 1));

  done  # var loop

#--------------------------------------------
# SURFEX run with reference initial conditions
#---------------------------------------------
  echo 'SURFEX run with reference vars...'
  cp -f PREP_REF.txt PREP.txt

  ./OFFLINE 2>stderr

  mv -f stderr stderr4

  for jvar in *.BIN
  do mv $jvar ${represults}${jvar}_${aa}${mm}${jj} 
     gzip  ${represults}${jvar}_${aa}${mm}${jj} ${represults}${jvar}_${aa}${mm}${jj}'.gz'
  done

  rm -f *.tex LISTING.txt *TXT
  cp -f SURFOUT.txt PREP.txt                   

if [ $AAAAMMJJRR -eq $AAAAMMJJRR_end ]; then
  cp -f SURFOUT.txt  $represults'SURFOUT_'$aaobs$mmobs${jjobs}'_'${expid}'.txt'
  gzip  ${represults}'SURFOUT_'$aaobs$mmobs${jjobs}'_'${expid}'.txt' ${represults}'SURFOUT_'$aaobs$mmobs${jjobs}'_'${expid}'.txt.gz'
fi
 
#-----------------------------------------------------------
# VARASSIM - storage of reference simulated observations
#-----------------------------------------------------------

  sed "s/LSIM        =  F/LSIM        =  T/" OPTIONS.nam > OPTIONS.nam_
  mv -f OPTIONS.nam_ OPTIONS.nam
  echo 'store references vars'

  ./VARASSIM > LISTING_varassim.txt 2>stderr 

  mv -f stderr stderr5
  sed "s/LSIM        =  T/LSIM        =  F/" OPTIONS.nam > OPTIONS.nam_
  mv -f OPTIONS.nam_ OPTIONS.nam

  mv -f OBSIMU               OBSIMU_REFR_$aaobs$mmobs${jjobs}'H'$RRobs.DAT
  mv -f MDSIMU               MDSIMU_REFR_$aaobs$mmobs${jjobs}'H'$RRobs.DAT

#---------------------------------------------------------------------------
# VARASSIM - evolve B - could be reset to BO periodically (test to be defined)
#----------------------------------------------------------
if [ $AAAAMMJJRR -eq $AAAAMMJJRR_start ] ; then
    mv -f  BGROUNDin0 BGROUNDin
  fi

#
  sed "s/LBEV        =  F/LBEV        =  T/" OPTIONS.nam > OPTIONS.nam_
  mv -f OPTIONS.nam_ OPTIONS.nam

  ./VARASSIM 2>stderr

  mv -f stderr stderr7
  sed "s/LBEV        =  T/LBEV        =  F/" OPTIONS.nam > OPTIONS.nam_
  mv -f OPTIONS.nam_ OPTIONS.nam
  mv -f BGROUNDout BGROUNDin


#----------------------------------------------------------
# VARASSIM - soil analysis
#-------------------------------
  echo 'Soil analysis ...'
  if [ -e $repobs'CANARI_NATURE_'$aaobs$mmobs${jjobs}'H'$RRobs.DAT ] ; then 
  echo 'File CANARI exists'
  cp -f $repobs'CANARI_NATURE_'$aaobs$mmobs${jjobs}'H'$RRobs.DAT $reprun/'CANARI_NATURE_'$aaobs$mmobs${jjobs}'H'$RRobs.DAT
  else
  echo "File CANARI does not exist"
  cp -f $repobs'CANARI_NATURE_NULL.DAT' $reprun/'CANARI_NATURE_'$aaobs$mmobs${jjobs}'H'$RRobs.DAT
  fi 

   ./VARASSIM 2>stderr
 
  mv -f stderr stderr6
  cp -f $reprun/'ANAL_INCR' $represults'ANAL_INCR_'${aa}${mm}${jj}.'txt'
  cp -f $reprun/'INNOV' $represults'INNOV_'${aa}${mm}${jj}.'txt'
  cp -f $reprun/'HO_WG2_v1' $represults'HO_WG2_WG1_'$aaobs$mmobs${jjobs}.txt
  cp -f $reprun/'HO_WG2_v2' $represults'HO_WG2_LAI_'$aaobs$mmobs${jjobs}.txt
  cp -f $reprun/'HO_LAI_v1' $represults'HO_LAI_LAI_'$aaobs$mmobs${jjobs}.txt
  cp -f $reprun/'HO_LAI_v2' $represults'HO_LAI_LAI_'$aaobs$mmobs${jjobs}.txt


  sed -i "s/WG2/oldWG2/" PREP.txt 
  sed -i "s/LAI/oldLAI/" PREP.txt 
  sed -i "s/BIOMA1/oldBIOMA1/" PREP.txt 

  cat PREP.txt >> SURFOUT.txt                   # changed for ASCII treatment of PREP
  mv -f SURFOUT.txt PREP.txt                    # changed for ASCII treatment of PREP
 
#----------------------------------------------
# Save PREP file for restarting (if needed)
#----------------------------------------------
if [ $jj -eq  1 ] ; then 
    cp -f  PREP.txt $represults'PREP_'$aaobs$mmobs${jj}'_'${expid}'.txt'
    gzip  ${represults}'PREP_'$aaobs$mmobs${jj}'_'${expid}'.txt' ${represults}'PREP_'$aaobs$mmobs${jj}'_'${expid}'.txt.gz'
fi

  cp -f PREP.txt PREP_REF.txt                   # this is the analysis ready for next cycle
  mv -f BGROUNDout BGROUNDin



#-------------------------
# Advance the time
#--------------------------
  echo 'end of period '$AAAAMMJJRR
  AAAAMMJJRR=$AAAAMMJJRRobs
#
rm -f ${reprun}/*.DAT
rm -f ${reprun}/*.nc

  done

mv -f $reprun/'OPTIONS.nam' $represults'OPTIONS.nam'
rm -f VARASSIM OFFLINE

echo '############################################'
echo 'RUN EKF EXITED AT TIME: '$AAAAMMJJRR
echo '############################################'

 




