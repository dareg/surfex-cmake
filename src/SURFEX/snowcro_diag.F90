!     #########

SUBROUTINE SNOWCRO_DIAG(HSNOWMETAMO, &
                        PSNOWDZ, PSNOWSWE, PSNOWRHO, PSNOWGRAN1, PSNOWGRAN2, PSNOWAGE, &
                        PSNOWHIST, PSNOWTEMP, PSNOWLIQ, PDIRCOSZW, PSNOWDEND, PSNOWSPHER, &
                        PSNOWSIZE, PSNOWSSA, PSNOWTYPEMEPRA, PSNOWRAM, PSNOWSHEAR, &
                        PACC_RAT,&
                        PSNOWDEPTH_1DAYS, PSNOWDEPTH_3DAYS, PSNOWDEPTH_5DAYS, &
                        PSNOWDEPTH_7DAYS, PSNOWSWE_1DAYS, PSNOWSWE_3DAYS, PSNOWSWE_5DAYS,&
                        PSNOWSWE_7DAYS, PSNOWRAM_SONDE, PSNOW_WETTHICKNESS, PSNOW_REFROZENTHICKNESS)

! Diagnostics of Crocus snowpack model
! Authors: P. Hagenmuller, Meteo-France, July 2016
! Modified Summer 2017 (P. Hagenmuller)
!
! Note that the Mepra diagnosis is the exact copy of the original Mepra (version in snowtools)
! and that this version explicitely contains incoherences (see comments in code and list below).
! In consequence, Mepra results should be considered for what it is worth.
!
!########################Mepra overall organization################################################!
!   0) Initialization of working variables
!   1) Loop on points and layers, to compute layer properties:
!       a) grain morphology: size (PSNOWSIZE), dendricity (PSNOWDEND), sphericity (PSNOWSPHER) and
!          snow type (PSNOWTYPEMEPRA)
!       b) mechanical properties: ram strength and shear strength

USE MODD_SURF_PAR,      ONLY : XUNDEF
USE MODD_CSTS,          ONLY : XRHOLI, XRHOLW

USE MODD_SNOW_PAR,ONLY :&
XX,XD1,XD2,XD3,&
!
JPTAB_DEND,JPTAB_NODEND,&
JP_PP_PP,JP_PP_DF,JP_DF_DF,JP_DF_RG,JP_DF_FC,JP_RG_RG,JP_RG_MF,&
JP_RG_FC,JP_FC_FC,JP_FC_DH,JP_DH_DH,JP_MF_MF,JP_MF_DH,JP_MF_FC,&
!
JPACC_HIG,JPACC_MOD,JPACC_LOW,JPACC_NUL,JPACC_NAN,&
XACC_RAT_HIG,XACC_RAT_MOD,XACC_SLA_STR,&
!
JPNAT_VLO,JPNAT_LOW,JPNAT_MOA,JPNAT_MOD,JPNAT_HIG,JPNAT_VHI,JPNAT_NAN,&
JPNAT_TAB,JPNAT_ACT,&
XNAT_RAT_HIG,XNAT_RAT_MOD,&
XNAT_HEI_HIG,XNAT_HEI_MOD,XNAT_HEI_LOW, XNAT_HEI_MIN,&
!
JPPRO_SUP_NAN,JPPRO_SUP_NEW,JPPRO_SUP_WET,JPPRO_SUP_FRO,&
XPRO_SUP_CRU,XPRO_SUP_DEP,&
!
JPPRO_INF_NAN,JPPRO_INF_HAR,JPPRO_INF_SOF,&
XPRO_INF_RAM,XPRO_INF_COE,&
!
JPAVA_NEW_DRY,JPAVA_NEW_WET,JPAVA_NEW_MIX,JPAVA_SLA_SUR,&
JPAVA_MEL_SUR,JPAVA_MEL_GRO,JPAVA_NAN
IMPLICIT NONE

! Options
CHARACTER(3), INTENT(IN)         :: HSNOWMETAMO ! metamorphism option

! Prognostic of snowcro (layer variables)
REAL, DIMENSION(:,:), INTENT(IN) :: PSNOWDZ          ! slope vertical layer height (projection in diag_misc_isban)(m)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWSWE         ! mass (kg/m2)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO         ! density (kg/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWGRAN1       ! grain morphology parameter 1 (-)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWGRAN2       ! grain morphology parameter 2 (-)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWAGE         ! age since snowfall (day)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWHIST        ! historical parameter (-) in {0,1,2,3,4,5} ######Why REAL ?
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWTEMP        ! temperature (K)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWLIQ         ! liquid water content (m) WARNING not in kg/m3???
!
! Characteristics of slope
REAL, DIMENSION(:),   INTENT(IN) :: PDIRCOSZW           ! cosine of slope angle (-)
!
! Diagnostic variables of Mepra and snowpro
! Layer variables
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWDEND          ! dendricity (-)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSPHER         ! sphericity (-)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSIZE          ! grain size (m)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSSA           ! specific surface area (m2/kg)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWTYPEMEPRA     ! snow type (-) INTEGER*1 is enough
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWRAM           ! ram penetration strength (kgf = 9.81 N)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSHEAR         ! shear strength (kgf/dm2 = 0.981 kPa)
REAL, DIMENSION(:,:), INTENT(OUT) :: PACC_RAT           ! accidental ratio shear strength/stress
!Layer-integrated variables
REAL, DIMENSION(:),   INTENT(OUT) :: PSNOWDEPTH_1DAYS   ! height of snow with age <= 1 day  (m)
REAL, DIMENSION(:),   INTENT(OUT) :: PSNOWDEPTH_3DAYS   ! height of snow with age <= 3 days (m)
REAL, DIMENSION(:),   INTENT(OUT) :: PSNOWDEPTH_5DAYS   ! height of snow with age <= 5 days (m)
REAL, DIMENSION(:),   INTENT(OUT) :: PSNOWDEPTH_7DAYS   ! height of snow with age <= 7 days (m)
REAL, DIMENSION(:),   INTENT(OUT) :: PSNOWSWE_1DAYS     ! swe with age <= 1 day  (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT) :: PSNOWSWE_3DAYS     ! swe with age <= 3 days (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT) :: PSNOWSWE_5DAYS     ! swe with age <= 5 days (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT) :: PSNOWSWE_7DAYS     ! swe with age <= 7 days (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT) :: PSNOWRAM_SONDE     ! 1st top penetration depth of ramsonde (m)
REAL, DIMENSION(:),   INTENT(OUT) :: PSNOW_WETTHICKNESS !
REAL, DIMENSION(:),   INTENT(OUT) :: PSNOW_REFROZENTHICKNESS  !
!
!
!
! Declarations of local variables used for calculations
! scalar
INTEGER :: ICLASS_DEND, ICLASS_SPHER, ICLASS_SIZE, ICLASS_HIST, ICLASS
REAL    :: ZDIAM                      ! for SSA calculations
REAL    :: ZRAM_FIN,ZRAM_DEN,ZRAM_ANG   ! for ram strength calculations
REAL    :: ZSHE_SPH,ZSHE_DEN,ZSHE_MTS,& ! for shear strength calculations
ZSHE_FE,ZSHE_FRE,ZSHE_FRE_SEC
REAL::ZEC,Z_TGM

LOGICAL,DIMENSION(SIZE(PSNOWSWE,1)) :: GRAM, GWET, GREFROZEN

LOGICAL :: LTHERM

LOGICAL :: GTHERMSTATE !thermal_state <= 2
REAL    :: ZSCW        !liquid water content (kg/m3)
REAL::EPSI!!!!!!!!!!!!VERRUE A REMPLACER PAR EPSI DE SURFEX QUAND TERMINE
INTEGER :: JJ       ! Array index of point
INTEGER :: JST      ! Array index of layer

! PRINT*,ASSOCIATED(PSNOWDEPTH_1DAYS),ASSOCIATED(PSNOWDEPTH_3DAYS)
! PRINT*,ASSOCIATED(PSNOWDEPTH_5DAYS),ASSOCIATED(PSNOWDEPTH_7DAYS)
! PRINT*,"in snowcrodiag"
! PRINT*,ALLOCATED(PSNOWDEPTH_1DAYS),ALLOCATED(PSNOWDEPTH_3DAYS)
! PRINT*,ALLOCATED(PSNOWDEPTH_5DAYS),ALLOCATED(PSNOWDEPTH_7DAYS)
! PRINT*,SIZE(PSNOWDEPTH_1DAYS),SIZE(PSNOWDEPTH_3DAYS)
! PRINT*,SIZE(PSNOWDEPTH_5DAYS),SIZE(PSNOWDEPTH_7DAYS)

! Initializations
! Scalar
GRAM      = .TRUE.
GWET      = .TRUE.
GREFROZEN = .TRUE.
EPSI      = 1e-16 !To change to correct XUNDEF

! One dimensional variables
PSNOWDEPTH_1DAYS= 0.
PSNOWDEPTH_3DAYS= 0.
PSNOWDEPTH_5DAYS= 0.
PSNOWDEPTH_7DAYS= 0.
PSNOWSWE_1DAYS  = 0.
PSNOWSWE_3DAYS  = 0.
PSNOWSWE_5DAYS  = 0.
PSNOWSWE_7DAYS  = 0.
PSNOWRAM_SONDE  = 0.
PSNOW_WETTHICKNESS      = 0.
PSNOW_REFROZENTHICKNESS = 0.

! Two dimensional variables
PSNOWDEND       = XUNDEF
PSNOWSPHER      = XUNDEF
PSNOWSIZE       = XUNDEF
PSNOWSSA        = XUNDEF
PSNOWTYPEMEPRA  = XUNDEF
PSNOWRAM        = XUNDEF
PSNOWSHEAR      = XUNDEF
PACC_RAT        = XUNDEF

PRINT *, 'TOTO'

DO JST=1,SIZE(PSNOWSWE,2)
  DO JJ=1,SIZE(PSNOWSWE,1)
    
    IF (PSNOWSWE(JJ,JST)>0) THEN
      !Do something only in case of non-empty layer

      !WARNING:
      ! PSNOWLIQ is not well exported before
      ! RIGHT PLACE TO DO THAT?
      PSNOWLIQ(JJ,JST) = PSNOWLIQ(JJ,JST) * PDIRCOSZW(JJ)

      !########################Grain morphology######################################################!
      ! Computes dendricity, sphericity, grain size, optical dimater and snow type from variables
      ! grain1 and grain2.

      IF (PSNOWGRAN1(JJ,JST) < 0) THEN
        ! Dendritic case

        ! Dendricity, sphericity and grain size
        PSNOWSIZE (JJ,JST) =   XUNDEF  !Grain size not defined for dendritic snow
        PSNOWDEND (JJ,JST) = - PSNOWGRAN1(JJ,JST) / XX
        PSNOWSPHER(JJ,JST) =   PSNOWGRAN2(JJ,JST) / XX

        ! Optical diameter for SSA diagnostic
        ZDIAM = PSNOWDEND(JJ,JST) * XD1 + (1 - PSNOWDEND(JJ,JST)) * &
        (PSNOWSPHER(JJ,JST) * XD2 + (1 - PSNOWSPHER(JJ,JST)) * XD3)
        ZDIAM = ZDIAM / 10000.

        ! 10 classes of dendricity 0:[0,0.1[, ..., 9:[0.9,1.0[ (value 1.0 does not exist)
        ICLASS_DEND = INT(10 * PSNOWDEND(JJ,JST))

        ! 10 classes of sphericity 0:[0,0.05[, 1:[0.05,0.15[, ..., 9:[0.85,1.0]
        !###########Strange way of defining sphericity classes -> Check with very old versions
        ICLASS_SPHER = MIN(INT(10 * PSNOWSPHER(JJ,JST) + 0.05),9)

        ! Overall 10x10 classes from 1 to 100 (included)
        ICLASS = 1 + ICLASS_DEND + ICLASS_SPHER * 10

        ! Snow type obtained in table JPTAB_DEND
        PSNOWTYPEMEPRA(JJ,JST) = JPTAB_DEND(ICLASS)


      ELSE
        ! Non dendritic case

        ! Dendricity,sphericty and grain size
        PSNOWSIZE (JJ,JST) = PSNOWGRAN2(JJ,JST)
        PSNOWDEND (JJ,JST) = 0
        PSNOWSPHER(JJ,JST) = PSNOWGRAN1(JJ,JST) / XX

        ! Optical diameter for SSA diagnostic
        ZDIAM = PSNOWSIZE(JJ,JST) * PSNOWSPHER(JJ,JST) + &
        MAX( 0.0004, 0.5*PSNOWSIZE(JJ,JST) ) * ( 1.-PSNOWSPHER(JJ,JST) )

        ! 10 classes of sphericity 0:[0,0.05[, 1:[0.05,0.15[, ..., 9:[0.85,1.0]
        !###########Strange way of defining sphericity classes -> Check with very old versions
        ICLASS_SPHER = MIN(INT(10 * PSNOWSPHER(JJ,JST) + 0.05),9)

        ! 6 classes of historical variable {0,1,...,5}
        ICLASS_HIST = NINT(PSNOWHIST(JJ,JST))

        ! 3 classes of grain size in mm 0:[0,0.55[, 1:[0.55,1.05[, 2:[1.05, +inf[
        !#########Strange +0.05
        IF     (PSNOWSIZE(JJ,JST) < 0.00055) THEN
          ICLASS_SIZE = 0
        ELSEIF (PSNOWSIZE(JJ,JST) < 0.00105) THEN
          ICLASS_SIZE = 1
        ELSE
          ICLASS_SIZE = 2
        ENDIF

        ! Overall 10x3x6 classes from 1 to 180 (included)
        ICLASS = 1 + ICLASS_SPHER + ICLASS_SIZE * 10 + ICLASS_HIST * 30

        !Snow type obtained in table JPTAB_NODEND
        PSNOWTYPEMEPRA(JJ,JST) = JPTAB_NODEND(ICLASS)

      ENDIF

      ! Additional constrain to define JP_RG_MF (MF/RG)
      IF ((PSNOWLIQ(JJ,JST) > EPSI) .AND. (PSNOWHIST(JJ,JST) < 2)) THEN !CHANGED FROM > 0 TO > XUEPSI
        PSNOWTYPEMEPRA(JJ,JST) = JP_RG_MF
      ENDIF


      ! For all snow types (dendritic and non-dendritic)


      !########################Thermal state######################################################!
      ! Table of thermal state definition
      !# T(deg C) tel(kg/m3)# tel < 5      # 5 <= tel < 50     # 50 <= tel#
      !####################################################################
      !#        T < -2      # 0            # 0                 # 0        #
      !# -2  <= T < -0.2    # 1            # 1                 # 1        #
      !#-0.2 <= T           # 2            # 3                 # 4        #

      ! In practice, only the distinction between thermal state >2 or <= 2 is considered here, i.e.
      ZSCW = XRHOLW * PSNOWLIQ(JJ,JST) / (PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ)) !ZSCW = liquid water content (kg/m3)
      GTHERMSTATE = (PSNOWTEMP(JJ,JST) < 272.96) .OR. (ZSCW < 5)           !True if thermstate <= 2


      !########################Ram strength######################################################!
      ! Computes penetration resistance force (PSNWORAM, kgf) of ramsonde as a function of snow type
      ! (PSNOWTYPEMEPRA), dendricity (PSNWODEND), sphericity (PSNOWSPHER), grain diameter (PSNWOSIZE,m)
      ! and density (PSNOWRHO,kg/m3)
      !
      ! WARNING: ram_strength expressed in kgf (i.e. weight of 1kg = 9.81 N) corresponding to a
      ! standard ramsonde. Roughly, the higher the ram strength, the higher the penetration resistance.
      ! But this is not an intrinsic material property, it is dependent on the measurement instrument.
      ! All threshold on ram strength in Mepra are thus also expressed in kgf.

      ! Variable ZRAM_FIN used several times
      ZRAM_FIN = 0.17 * PSNOWRHO(JJ,JST) - 31

      IF (PSNOWDEND(JJ,JST) > 0) THEN
        ! Dendritic case
        ZRAM_DEN = MAX(1., 0.018 * PSNOWRHO(JJ,JST) - 1.363)
        ZRAM_ANG = MAX(2., ZRAM_FIN * PSNOWSPHER(JJ,JST) + (1 - PSNOWSPHER(JJ,JST)) * (0.5 * ZRAM_FIN + 0.6))
        PSNOWRAM(JJ,JST) = PSNOWDEND(JJ,JST) * ZRAM_DEN + (1 - PSNOWDEND(JJ,JST)) * ZRAM_ANG

      ELSE
        ! Non-dendritic cases
        SELECT CASE (INT(PSNOWTYPEMEPRA(JJ,JST)))

        ! RG type
        CASE (JP_RG_RG)
          IF (PSNOWRHO(JJ,JST) < 200) THEN
            PSNOWRAM(JJ,JST) = 3
          ELSE
            PSNOWRAM(JJ,JST) = ZRAM_FIN
          ENDIF

        ! RG/FC type
        CASE (JP_RG_FC)
          IF (PSNOWRHO(JJ,JST) < 200) THEN
            PSNOWRAM(JJ,JST) = 2
          ELSE
            PSNOWRAM(JJ,JST) = PSNOWSPHER(JJ,JST)  * ZRAM_FIN +&
            (1- PSNOWSPHER(JJ,JST)) * (ZRAM_FIN * (0.8-1000*PSNOWSIZE(JJ,JST)) + 2000*PSNOWSIZE(JJ,JST))
          ENDIF

        ! FC and FC/DH types
        CASE (JP_FC_FC,JP_FC_DH)
          IF(PSNOWSIZE(JJ,JST) > 0.0008) THEN
            PSNOWRAM(JJ,JST) = 2
          ELSEIF (PSNOWRHO(JJ,JST) < 200) THEN
            PSNOWRAM(JJ,JST) = 3        * (0.8 - 1000 * PSNOWSIZE(JJ,JST)) + 2000 * PSNOWSIZE(JJ,JST)
          ELSE
            PSNOWRAM(JJ,JST) = ZRAM_FIN * (0.8 - 1000 * PSNOWSIZE(JJ,JST)) + 2000 * PSNOWSIZE(JJ,JST)
          ENDIF

        ! MF/RG, MF, MF/DH, MF/FC types
        CASE (JP_RG_MF,JP_MF_MF,JP_MF_DH,JP_MF_FC)
          IF (GTHERMSTATE) THEN
            PSNOWRAM(JJ,JST) = MAX(10., 0.103 * PSNOWRHO(JJ,JST) - 19.666)
          ELSEIF (PSNOWRHO(JJ,JST) < 250) THEN
            PSNOWRAM(JJ,JST) = 1
          ELSEIF (PSNOWRHO(JJ,JST) < 350) THEN
            PSNOWRAM(JJ,JST) = 2
          ELSE
            PSNOWRAM(JJ,JST) = 0.16 * PSNOWRHO(JJ,JST) - 54
          ENDIF

        ! All other cases
        CASE DEFAULT
          PSNOWRAM(JJ,JST) = 2

        END SELECT
      ENDIF

      !########################Shear strength#######################################################
      ! Return the shear strength (kgf/dm2= 0.981 kPa) as a function
      ! of snow type, dendricity, sphericity, grain diameter (m), historical variable, density (rho, kg/m3),
      ! liquid water content (kg/m3) and thermal state.
      !
      ! WARNING: One of the most time-consuming calculation, here. The evaluation of the paramaterization
      ! by J.P. Navarre is elusive. The parameterization is threshold sensitive. Moreover, the complexity
      ! of the parameterization appears, in some cases (e.g. C_den, C_spher), to be useless.

      IF ((.NOT.GTHERMSTATE).AND.((PSNOWTYPEMEPRA(JJ,JST)==JP_RG_FC).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_RG_RG).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_RG).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_RG_MF).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_DH).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_FC).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_MF))) THEN

        ! Case of wet snow of type RG, DF+RG, MF+RG, MF+DH, MF+FC, MF and RG+FC
        IF     (PSNOWRHO(JJ,JST) < 200) THEN
          PSNOWSHEAR(JJ,JST) = 0.1
        ELSEIF (PSNOWRHO(JJ,JST) < 320) THEN
          PSNOWSHEAR(JJ,JST) = 0.02  * PSNOWRHO(JJ,JST) - 3.9
        ELSE
          PSNOWSHEAR(JJ,JST) = 0.068 * PSNOWRHO(JJ,JST) - 18.64
        ENDIF

      ELSE
        ! Other case

        ! General formula for not wet snow
        ! shear_strength = ZSHE_SPH * ZSHE_DEN * ZSHE_MTS * ZSHE_FE * ZSHE_FRE) * (rho^2/10^4 - 0.6) + 0.12 with
        ! ZSHE_SPH coefficient of sphericity/faceting effect,
        ! ZSHE_DEN coefficient of dendrites effect,
        ! ZSHE_MTS coefficient of size effect
        ! ZSHE_FE  coefficient of water content
        ! ZSHE_FRE coefficient of melt/freeze cycles


        ! ZEC corresponds to the ratio between the liquid water content (kg/m3) and
        ! the theoretical maximum liquid water content retaind by the snow layer,
        ! defined as 5% of solid ice volumetric content (m3/m3)
        ZEC = ZSCW / (50. * (1. + (ZSCW - PSNOWRHO(JJ,JST))/XRHOLI))

        ! ZSHE_SPH sphericity
        ! Qualitatively ZSHE_SPH increases almost linearly from (spher,ZSHE_SPH)=(0,0.45),
        ! i.e. faceted snow, to (spher,ZSHE_SPH)=(1.0,1.15), i.e. roundish snow.
        ! The special case for histo = 3 or 5 and spher = 0.8 is to constrained the
        ! sphericity below 0.8 for snow that has been angular and humid (relatively
        ! resistant, even angular).

        IF ((PSNOWSPHER(JJ,JST) > 0.8).AND.&
           ((PSNOWHIST(JJ,JST)==3).OR.(PSNOWHIST(JJ,JST)==5))) THEN
          ZSHE_SPH = 1.05
        ELSEIF (PSNOWSPHER(JJ,JST) <= 0.25) THEN
          ZSHE_SPH = 0.45  + 0.7  * PSNOWSPHER(JJ,JST)
        ELSEIF (PSNOWSPHER(JJ,JST) <= 0.50) THEN
          ZSHE_SPH = 0.625 + 1.0 * (PSNOWSPHER(JJ,JST) - 0.25)
        ELSEIF (PSNOWSPHER(JJ,JST) <= 0.75) THEN
          ZSHE_SPH = 0.875 + 0.6 * (PSNOWSPHER(JJ,JST) - 0.50)
        ELSE
          ZSHE_SPH = 1.025 + 0.5 * (PSNOWSPHER(JJ,JST) - 0.75)
        ENDIF

        ! ZSHE_DEN dendricity
        ! Qualitatively, ZSHE_DEN decreases almost linearly from (dendr,ZSHE_DEN) = (0,1),
        ! i.e. old snow, to (dendr,ZSHE_DEN) = (1.0,0.45), i.e. recent snow.

        IF (PSNOWDEND(JJ,JST) <= 0.25) THEN
          ZSHE_DEN = 1.0 - 0.4 *  PSNOWDEND(JJ,JST)
        ELSEIF (PSNOWDEND(JJ,JST) <= 0.50) THEN
          ZSHE_DEN = 0.9 - 0.4 * (PSNOWDEND(JJ,JST) - 0.25)
        ELSEIF (PSNOWDEND(JJ,JST) <= 0.75) THEN
          ZSHE_DEN = 0.8 - 0.8 * (PSNOWDEND(JJ,JST) - 0.50)
        ELSE
          ZSHE_DEN = 0.6 - 0.6 * (PSNOWDEND(JJ,JST) - 0.75)
        ENDIF

        ! ZSHE_MTS
        ! As expected, the grain size has no effect on dendritic snow.
        ! Qualitatively, on non dendritic snow with a grain size larger than a kind of simplified optical radius (Z_TGM),
        ! ZSHE_MTS accounts for the fact that the larger the grains, the smaller the strength, especially in case
        ! of faceted snow types.

        Z_TGM = 0.0004 - 0.0001 * PSNOWSPHER(JJ,JST)

        IF ((PSNOWDEND(JJ,JST) > 0).OR.(PSNOWSIZE(JJ,JST) <= Z_TGM)) THEN
          ZSHE_MTS = 1
        ELSE
          ZSHE_MTS = 1 - (0.8 - 0.2 * PSNOWSPHER(JJ,JST)) * 530 * (PSNOWSIZE(JJ,JST) - Z_TGM)
        ENDIF

        ! ZSHE_FE Wetting of snow
        ! For low liquid water content (lower than thresold depending on density),
        ! C_fe increases linearly with the liquid water content (max = 1.1). For larger values,
        ! C_fe decreases linearly with the liquid water content until C_fe reaches about 0.62,
        ! then it continues to decrease but more slowly.

        IF (ZEC <= 0.1) THEN
          ZSHE_FE = 1 + ZEC
        ELSEIF (ZEC <= 0.3) THEN
          ZSHE_FE = 1.335 - 2.35 * ZEC
        ELSEIF (ZEC <  0.9) THEN
          ZSHE_FE = 0.750 - 0.4  * ZEC
        ELSE
          ZSHE_FE = MAX(0.15,MIN(0.35,(PSNOWRHO(JJ,JST)-ZSCW)/1000))
        ENDIF

        ! ZSHE_FRE takes into account the strengthening effect of melt/freeze cycles.
        ! In case of never humid snow (histo = 0 or 1), C_fre = 1.0
        ! In case of "has been humid snow" (histo = 2 or 3) (but no cycle),
        ! C_fre = 1.0 if the snow is still humid, and C_fre about 2.0 if completely dry (ec=0).
        ! In case of melt/freeze cycle snow, C_fre decreases from about 1.95 to 1 with
        ! liquid water content, with 1/scw trend.

        IF ((PSNOWHIST(JJ,JST) <= 1).OR.(ZEC > 0.5)) THEN
          ZSHE_FRE = 1.
        ELSEIF (ZEC<EPSI) THEN !Changed from==0 to <1e-16
          ! C_fre_sec = 1.5 * ((1.15 + 0.2 * (1.-spher)) / 1.15) * (1 + 0.2/C_mts)
          ZSHE_FRE = (1.7608695652 - 0.2608695652 * PSNOWSPHER(JJ,JST)) * (1. + 0.2/ZSHE_MTS)
        ELSEIF (PSNOWHIST(JJ,JST) <= 3) THEN
          ZSHE_FRE = 1
        ELSEIF (ZEC <= 0.1) THEN
          ZSHE_FRE = -2. * ZEC + 1.5
        ELSE
          ZSHE_FRE = -0.75 * ZEC + 1.375
        ENDIF

        PSNOWSHEAR(JJ,JST) = MAX(0.05, &
        ZSHE_SPH * ZSHE_DEN * ZSHE_MTS * ZSHE_FE * ZSHE_FRE *&
        (PSNOWRHO(JJ,JST)*PSNOWRHO(JJ,JST) / 10000. -0.6) + 0.12)

      ENDIF
!
!
!
!
!
!
!
!
      ! All cases
      ! Compute depth and SWE of recent snow
      IF(PSNOWAGE(JJ,JST)<=1)THEN
        PSNOWDEPTH_1DAYS(JJ) = PSNOWDEPTH_1DAYS(JJ) + PSNOWDZ (JJ,JST)
        PSNOWSWE_1DAYS  (JJ) = PSNOWSWE_1DAYS  (JJ) + PSNOWSWE(JJ,JST)
      ENDIF

      IF(PSNOWAGE(JJ,JST)<=3)THEN
        PSNOWDEPTH_3DAYS(JJ) = PSNOWDEPTH_3DAYS(JJ) + PSNOWDZ (JJ,JST)    
        PSNOWSWE_3DAYS  (JJ) = PSNOWSWE_3DAYS  (JJ) + PSNOWSWE(JJ,JST)
      ENDIF
    
      IF(PSNOWAGE(JJ,JST)<=5)THEN
        PSNOWDEPTH_5DAYS(JJ) = PSNOWDEPTH_5DAYS(JJ) + PSNOWDZ (JJ,JST)    
        PSNOWSWE_5DAYS  (JJ) = PSNOWSWE_5DAYS  (JJ) + PSNOWSWE(JJ,JST)
      ENDIF
  
      IF(PSNOWAGE(JJ,JST)<=7)THEN
        PSNOWDEPTH_7DAYS(JJ) = PSNOWDEPTH_7DAYS(JJ) + PSNOWDZ (JJ,JST)    
        PSNOWSWE_7DAYS  (JJ) = PSNOWSWE_7DAYS  (JJ) + PSNOWSWE(JJ,JST)
      END IF
    
      ! Ram sonde penetration
      IF ((GRAM(JJ)).AND.(PSNOWRAM(JJ,JST)<=2.)) THEN
        PSNOWRAM_SONDE(JJ)=PSNOWRAM_SONDE(JJ)+PSNOWDZ(JJ,JST)
      ELSE
        GRAM(JJ)=.FALSE.
      ENDIF

      ! Depth of wet snow
      IF ((GWET(JJ)).AND.(PSNOWLIQ(JJ,JST)>0.)) THEN
        PSNOW_WETTHICKNESS(JJ)=PSNOW_WETTHICKNESS(JJ)+PSNOWDZ(JJ,JST)
      ELSE
        GWET(JJ)=.FALSE.
      ENDIF
      ! Depth of refrozen snow
      IF ((GREFROZEN(JJ)).AND.(PSNOWHIST(JJ,JST)>=2).AND.(PSNOWTEMP(JJ,JST)<273.15)) THEN
        PSNOW_REFROZENTHICKNESS(JJ)=PSNOW_REFROZENTHICKNESS(JJ)+PSNOWDZ(JJ,JST)
      ELSE
        GREFROZEN(JJ)=.FALSE.
      ENDIF    
    
      ! Specific surface area
      IF ( HSNOWMETAMO=='B92' ) THEN
        PSNOWSSA(JJ,JST) = 6. / (XRHOLI*ZDIAM)
      ELSE
        PSNOWSSA(JJ,JST) = 6. / (XRHOLI*PSNOWGRAN1(JJ,JST))
      END IF

    END IF
  END DO
END DO
  
END SUBROUTINE SNOWCRO_DIAG
