!     #########

SUBROUTINE SNOWCRO_DIAG(HSNOWMETAMO, &
                        PSNOWDZ, PSNOWSWE, PSNOWRHO, PSNOWGRAN1, PSNOWGRAN2, PSNOWAGE, &
                        PSNOWHIST, PSNOWTEMP, PSNOWLIQ, PDIRCOSZW, PSNOWDEND, PSNOWSPHER, &
                        PSNOWSIZE, PSNOWSSA, PSNOWTYPEMEPRA, PSNOWRAM, PSNOWSHEAR, &
                        PACC_RAT, PNAT_RAT, &
                        PSNOWDEPTH_1DAYS, PSNOWDEPTH_3DAYS, PSNOWDEPTH_5DAYS, PSNOWDEPTH_7DAYS,&
                        PSNOWSWE_1DAYS, PSNOWSWE_3DAYS, PSNOWSWE_5DAYS,PSNOWSWE_7DAYS,&
                        PSNOWRAM_SONDE, PSNOW_WETTHICKNESS, PSNOW_REFTHICKNESS,&
                        PDEP_HIG, PDEP_MOD, PDEP_SUP, PDEP_TOT, PDEP_HUM,&
                        PACC_LEV, PNAT_LEV, PPRO_SUP_TYP, PPRO_INF_TYP, PAVA_TYP)
!
! Diagnostics of Crocus snowpack model
! Authors: P. Hagenmuller, Meteo-France, July 2016
! Modified Summer 2017 (P. Hagenmuller)
!
!
! Note that the Mepra diagnosis is the exact copy of the original Mepra (version in snowtools)
! and that this version explicitely contains incoherences (see comments in code and list below).
! In consequence, Mepra results should be considered for what they are worth.
!
!
!########################Mepra overall organization################################################!
!   0) Initialization of working variables
!   1) Loop on points and layers, to compute layer properties:
!       a) grain morphology: size (PSNOWSIZE), dendricity (PSNOWDEND), sphericity (PSNOWSPHER) and
!          snow type (PSNOWTYPEMEPRA)
!       b) mechanical properties: ram strength and shear strength
!
!
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
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDZ             ! slope-parallel layer height (projection in diag_misc_isban)(m)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWSWE            ! mass (kg/m2)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO            ! density (kg/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWGRAN1          ! grain morphology parameter 1 (-)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWGRAN2          ! grain morphology parameter 2 (-)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWAGE            ! age since snowfall (day)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWHIST           ! historical parameter (-) in {0-5} ######Why REAL ?
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWTEMP           ! temperature (K)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWLIQ            ! vertical liquid water content (m)
!
! Characteristics of slope
REAL, DIMENSION(:),   INTENT(IN)    :: PDIRCOSZW           ! cosine of slope angle (-)
!
! Diagnostic variables of Mepra and snowpro
! Layer variables
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWDEND           ! dendricity (-)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSPHER          ! sphericity (-)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSIZE           ! grain size (m)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSSA            ! specific surface area (m2/kg)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWTYPEMEPRA      ! snow type (-) INTEGER*1 is enough
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWRAM            ! ram penetration strength (kgf = 9.81 N)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSHEAR          ! shear strength (kgf/dm2 = 0.981 kPa)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PACC_RAT            ! accidental ratio shear strength/stress
REAL, DIMENSION(:,:), INTENT(OUT)   :: PNAT_RAT            ! natural ratio shear strength/stress
!Layer-integrated variables
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWDEPTH_1DAYS    ! height of snow with age <= 1 day  (m)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWDEPTH_3DAYS    ! height of snow with age <= 3 days (m)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWDEPTH_5DAYS    ! height of snow with age <= 5 days (m)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWDEPTH_7DAYS    ! height of snow with age <= 7 days (m)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWSWE_1DAYS      ! swe with age <= 1 day  (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWSWE_3DAYS      ! swe with age <= 3 days (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWSWE_5DAYS      ! swe with age <= 5 days (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWSWE_7DAYS      ! swe with age <= 7 days (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWRAM_SONDE      ! top penetration depth of ramsonde (m)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOW_WETTHICKNESS  !
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOW_REFTHICKNESS  !
REAL, DIMENSION(:),   INTENT(OUT)   :: PDEP_HIG            ! depth of high instability (m)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDEP_MOD            ! depth of moderate instability (m)
REAL, DIMENSION(:),   INTENT(INOUT) :: PDEP_SUP            ! snow depth in superior profile(m)
REAL, DIMENSION(:),   INTENT(INOUT) :: PDEP_TOT            ! total snow depth (m)
REAL, DIMENSION(:),   INTENT(INOUT) :: PDEP_HUM            ! height of the uppest continuous block of humid snow in the sup
REAL, DIMENSION(:),   INTENT(OUT)   :: PACC_LEV            ! accidental risk index (0-4)
REAL, DIMENSION(:),   INTENT(INOUT) :: PNAT_LEV            ! natural risk index (0-6)
REAL, DIMENSION(:),   INTENT(INOUT) :: PPRO_SUP_TYP        ! type of superior profile (0, 4, 5, 6)
REAL, DIMENSION(:),   INTENT(OUT)   :: PPRO_INF_TYP        ! type of inferior profile (0, 1, 6)
REAL, DIMENSION(:),   INTENT(INOUT) :: PAVA_TYP            ! type of avalanche (0-6)
!
!
!
!
! Declarations of local variables used for calculations
! scalar
INTEGER :: ICLASS_DEND, ICLASS_SPHER,&                     ! for snow type calculations
           ICLASS_SIZE, ICLASS_HIST, ICLASS,&              ! for snow type calculations
           JJ, JST,&                                       ! array indices (point, layer)
           IPRO_CLASS,&                                    ! for natural risk index calculations
           IVAL_HIG, IVAL_MOD, IVAL_LOW,&                  ! for natural risk index calculations
           IACC_FROMNAT
REAL    :: ZDIAM,&                                         ! for SSA calculations
           ZRAM_FIN, ZRAM_DEN, ZRAM_ANG,&                  ! for ram strength calculations
           ZSHE_SPH, ZSHE_DEN, ZSHE_MTS,&                  ! for shear strength calculations
           ZSHE_FE , ZSHE_FRE, ZSHE_FRE_SEC,&              ! for shear strength calculations
           ZEC,Z_TGM,&                                     ! for XXXX
           ZSCW,&                                          ! liquid water content (kg/m3)
           ZSNOWDZ,&
           ZSNOWSWE,&
           ZSNOW_DEPTH_TOP,&
           ZSNOW_HEIGHT_TOP,&
           ZNAT_LEV_TMP
LOGICAL :: GMF, GPP,&                                      ! combination of snow types
           GTHERMSTATE                                     ! thermal_state <= 2




REAL::EPSI!!!!!!!!!!!!VERRUE A REMPLACER PAR EPSI DE SURFEX QUAND TERMINE


!1d (location dimension)
REAL     , DIMENSION(SIZE(PSNOWSWE,1)) :: ZWEIGH_STRESS,&
                                          ZSKIER_STRESS,&
                                          ZBETA,&
                                          ZSNOW_THICK,&
                                          ZCRUST_THICKNESS,&    !Total thickness of frozen MF-like
                                          ZSNOW_DEPTH,&         !Depth of current layer
                                          ZPRO_CRUST,&          !Thickness of continuous crust
                                          ZDEP_INF,&
                                          ZHUMTHICK,&
                                          ZACC_HIG_DEP,&
                                          ZACC_MOD_DEP,&
                                          ZNAT_HIG_DEP,&
                                          ZNAT_MOD_DEP,&
                                          ZDEP_SUP_PRE,&
                                          ZDEP_TOT_PRE,&
                                          ZDEP_HUM_PRE,&
                                          ZAVA_TYP_PRE,&
                                          ZNAT_LEV_PRE,&
                                          ZPRO_SUP_TYP_PRE,&
                                          ZTHICK_BOT            ! Thickness of the bottom layer
!
INTEGER  , DIMENSION(SIZE(PSNOWSWE,1)) :: IPRO_SUP_LIM
LOGICAL,   DIMENSION(SIZE(PSNOWSWE,1)) :: GRAM,&
                                          GWET,&
                                          GREFROZEN,&
                                          GBELOW_SLAB,&
                                          GPREV_ISNOT_MF,&
                                          GACC_MOD_TEMP,&
                                          GCOULD_BE_NEW,&
                                          GIKNOTMF,&
                                          GFOUND,&
                                          GHUM,&                !Does all layers in sup pro and heigth > thr have thermal state > 2?
                                          GDRY,&                !Does all layers in sup pro and height > thr have thermal state <=2?
                                          GMEL_GRO,&            !Is there one layer in sup pro between depth > thr and height >thr with thermal state <= 2
                                          GTHERMSTATE_BOT       !thermal_state <= 2 of bottom layer
!
!
!
! Initializations
EPSI      = 1e-16 !To change to correct XUNDEF

!Saving previous time step variables
DO JJ=1,SIZE(PSNOWSWE,1)
  ZDEP_SUP_PRE    (JJ) = PDEP_SUP    (JJ)
  ZDEP_TOT_PRE    (JJ) = PDEP_TOT    (JJ)
  ZDEP_HUM_PRE    (JJ) = PDEP_HUM    (JJ)
  ZAVA_TYP_PRE    (JJ) = PAVA_TYP    (JJ)
  ZNAT_LEV_PRE    (JJ) = PNAT_LEV    (JJ)
  ZPRO_SUP_TYP_PRE(JJ) = PPRO_SUP_TYP(JJ)
ENDDO

! Tests
GRAM            = .TRUE.
GWET            = .TRUE.
GREFROZEN       = .TRUE.
GACC_MOD_TEMP   = .FALSE.
GBELOW_SLAB     = .FALSE.
GPREV_ISNOT_MF  = .TRUE. !To be changed in GPREV_IS_MF ?
GIKNOTMF        = .FALSE.
GHUM            = .TRUE. !OK
GDRY            = .TRUE.
GMEL_GRO        = .TRUE.
GCOULD_BE_NEW   = .FALSE.
GFOUND          = .FALSE.

GTHERMSTATE_BOT = .FALSE. ! not necessary
ZTHICK_BOT      = 0.  !not necessary
!Local cumulative


! One dimensional variables (cumulative!!!)
PSNOWDEPTH_1DAYS        = 0.
PSNOWDEPTH_3DAYS        = 0.
PSNOWDEPTH_5DAYS        = 0.
PSNOWDEPTH_7DAYS        = 0.
PSNOWSWE_1DAYS          = 0.
PSNOWSWE_3DAYS          = 0.
PSNOWSWE_5DAYS          = 0.
PSNOWSWE_7DAYS          = 0.
PSNOWRAM_SONDE          = 0.
PSNOW_WETTHICKNESS      = 0.
PSNOW_REFTHICKNESS      = 0.
ZWEIGH_STRESS           = 0.
ZSKIER_STRESS           = 0.
ZBETA                   = 0.
ZSNOW_THICK             = 0.
ZCRUST_THICKNESS        = 0.
ZSNOW_DEPTH             = 0.
ZHUMTHICK               = 0.

ZPRO_CRUST      = 0.
PDEP_HUM        = 0.
PDEP_TOT = 0.
PDEP_SUP = 0.

! One dimensional (non cumulative)
ZNAT_HIG_DEP            = 0
ZNAT_MOD_DEP            = 0
ZACC_HIG_DEP            = XUNDEF
ZACC_MOD_DEP            = XUNDEF
IPRO_SUP_LIM = 0


! Two dimensional variables (intitalization not absolutely necessary)
PSNOWDEND       = XUNDEF
PSNOWSPHER      = XUNDEF
PSNOWSIZE       = XUNDEF
PSNOWSSA        = XUNDEF
PSNOWTYPEMEPRA  = XUNDEF
PSNOWRAM        = XUNDEF
PSNOWSHEAR      = XUNDEF
PACC_RAT        = XUNDEF
PNAT_RAT        = XUNDEF
PDEP_HIG        = 0
PDEP_MOD        = 0
PACC_LEV        = JPACC_LOW
PPRO_INF_TYP    = JPPRO_INF_NAN
PPRO_SUP_TYP    = JPPRO_SUP_NAN
PAVA_TYP        = JPAVA_NAN!in case no real layers
PNAT_LEV        = JPNAT_NAN!JPNAT_VLO

! merge ZSNOW_THICK and ZSNOw_DEPTH

PRINT *, 'TOTO'

DO JST=1,SIZE(PSNOWSWE,2)
  DO JJ=1,SIZE(PSNOWSWE,1)
!
    IF (PSNOWSWE(JJ,JST)>0) THEN
!     Do something only in case of non-empty layer
!
!     WARNING:
!     PSNOWLIQ is not well exported before
!     RIGHT PLACE TO DO THAT?
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
        ZDIAM = PSNOWDEND(JJ,JST) * XD1 +(1 - PSNOWDEND(JJ,JST)) * &
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
      ZSCW = XRHOLW * PSNOWLIQ(JJ,JST) / (PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ)) !ZSCW = LWC (kg/m3)
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


      !########################Strength stress ratio######################################################!
      ! Computes the ratio PACC_RAT between the shear strength (PSNOWSHEAR kgf/dm2 = 0.981 kPa) of the
      ! current layer and the shear stress due to a skier and to the overlying snow weight at the TOP of the layer
      ! (= bottom of just above layer).
      ! and the ratio PNAT_RAT between the shear strength (PSNOWSHEAR kgf/dm2 = 0.981 kPa) and the shear
      ! stress due only to the weight of the layers at the BOTTOM of the current layer.
      !
      ! Stress skier takes into account that the additional stress induced by the skier
      ! is more distributed (so lower) in the snowpack far below the skier than just below the skier.
      ! The calculation is based on the simplification of the semi-infinite uniform elastic layer
      ! theory developed by Boussinesq. Expert rules are used to account for
      ! bridging effects of different snow types (calculation of beta).
      !
      ! Low value of beta indicates high bridging effect (i.e. decrease of max shear stress)
      ! slab_thr = 1.5          #Threshold on shear strength (kgf/dm2)
      ! beta_refrozen = 0.5     #(MF/RG, MF, MF/DH or MF/FC) and thermal_state <= 2
      ! beta_humid = 1.1        #(MF/RG, MF, MF/DH or MF/FC) and thermal_state >  2
      ! beta_evolved = 1.0      #Other snow types and shear_strength >  slab_thr
      ! beta_recent_dry = 1.2   #Other snow types and shear_strength <= slab_thr
      !
      ! WARNING: two slightly different versions in "python snowtools" or "MEPRA fortran" of PACC_RAT.
      ! The version in MEPRA fortran is used to calculate the accidental risk and the version
      ! in snowtools just provides the value of 'MEPRA_ACCIDENTAL_RATIO' but is not used for
      ! other calculations. Differences are observed on the calculation of contps
      ! (e.g. -0.015 *snow_depth + 1.45). And the indexing of beta is unclear.
      ! Corrected in latest surfex version (24.11.2015) ?
      !
      ! WARNING: not defined for slope angle = 0 and first top layer
      ! WARNING !!!!!!!!!!!!!! the division by PDIRCOSZW(JJ) does not make sense but agrees with snowtools (check snowtools)


      ! Update of beta (bridging factor)
      IF (JST==1) THEN
        ZBETA(JJ) = 0
      ELSE
        IF ((PSNOWTYPEMEPRA(JJ,JST-1)==JP_MF_MF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST-1)==JP_RG_MF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST-1)==JP_MF_FC).OR.&
            (PSNOWTYPEMEPRA(JJ,JST-1)==JP_MF_DH)) THEN
          IF (PSNOWTEMP(JJ,JST-1) < 272.96) THEN
            ZBETA(JJ) = ZBETA(JJ) + 0.5 * PSNOWDZ(JJ,JST-1) / PDIRCOSZW(JJ)
          ELSE
            ZBETA(JJ) = ZBETA(JJ) + 1.1 * PSNOWDZ(JJ,JST-1) / PDIRCOSZW(JJ)
          ENDIF
        ELSE
          IF (PSNOWSHEAR(JJ,JST-1) > 1.5) THEN
            ZBETA(JJ) = ZBETA(JJ) + 1.0 * PSNOWDZ(JJ,JST-1) / PDIRCOSZW(JJ)
          ELSE
            ZBETA(JJ) = ZBETA(JJ) + 1.2 * PSNOWDZ(JJ,JST-1) / PDIRCOSZW(JJ)
          ENDIF
        ENDIF
      ENDIF

      IF((PDIRCOSZW(JJ)<1).AND.(JST>1)) THEN
        ! ACC_RAT is defined only for slope > 25, i.e 40 and not for the top layer
        ! total_stress(jst) = weight_stress(jst-1) + beta(jst) * skier_stress(jst-1)
        PACC_RAT(JJ,JST) = ZWEIGH_STRESS(JJ) + ZBETA(JJ) / ZSNOW_THICK(JJ) * ZSKIER_STRESS(JJ)
        PACC_RAT(JJ,JST) = PSNOWSHEAR(JJ,JST) / PACC_RAT(JJ,JST)
      ELSE
        PACC_RAT(JJ,JST) = XUNDEF !-1 in previous versions
      ENDIF

      ! Update of snow thickness
      ! ZSNOW_THICK = thickness of top of current JST layer -> thickness of bottom of current JST layer
      ZSNOW_THICK(JJ) = ZSNOW_THICK(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)

      ! Update of skier_stress
      IF (    ZSNOW_THICK(JJ) < 0.10) THEN
        ZSKIER_STRESS(JJ) = -15.0 * ZSNOW_THICK(JJ) + 4.00
      ELSEIF (ZSNOW_THICK(JJ) < 0.15) THEN
        ZSKIER_STRESS(JJ) = -10.0 * ZSNOW_THICK(JJ) + 3.50
      ELSEIF (ZSNOW_THICK(JJ) < 0.20) THEN
        ZSKIER_STRESS(JJ) = - 8.0 * ZSNOW_THICK(JJ) + 3.20
      ELSEIF (ZSNOW_THICK(JJ) < 0.35) THEN
        ZSKIER_STRESS(JJ) = - 4.0 * ZSNOW_THICK(JJ) + 2.40
      ELSEIF (ZSNOW_THICK(JJ) < 0.50) THEN
        ZSKIER_STRESS(JJ) = - 2.0 * ZSNOW_THICK(JJ) + 1.70
      ELSEIF (ZSNOW_THICK(JJ) < 0.80) THEN
        ZSKIER_STRESS(JJ) = - 1.5 * ZSNOW_THICK(JJ) + 1.45
      ELSE
        ZSKIER_STRESS(JJ) = 0
      ENDIF

      ZSKIER_STRESS(JJ) = 1.4 * ZSKIER_STRESS(JJ) !Could be incorporated above

      ! Update of shear stress due to overlying layers
      ZWEIGH_STRESS(JJ) = ZWEIGH_STRESS(JJ) + &
      PSNOWRHO(JJ,JST) * PSNOWDZ(JJ,JST) * SQRT(1-PDIRCOSZW(JJ)*PDIRCOSZW(JJ)) / 100.

      IF((PDIRCOSZW(JJ)<1).AND.(JST>1)) THEN
        ! ACC_RAT is defined only for slope > 10 degrees and not for the top layer
        PNAT_RAT(JJ,JST) = PSNOWSHEAR(JJ,JST) / ZWEIGH_STRESS(JJ)
      ELSE
        PNAT_RAT(JJ,JST) = XUNDEF !-1 in previous versions
      ENDIF

      !##########################################Accidental risk###################################
      !
      ! Update of snow depth (check with ZSNOW_THICK ??)
      ZSNOW_DEPTH(JJ) = ZSNOW_DEPTH(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)

      ! Depth of top of current layer
      ZSNOW_DEPTH_TOP = ZSNOW_DEPTH(JJ) - PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)

      ! Update of total thickness of crusts above current layer (included)
      ! In practice, it is not necessary to calculate it for the whole snowpack
      ! when a high instability is already found
      IF((PSNOWTEMP(JJ,JST) < 272.96).AND.&
         ((PSNOWTYPEMEPRA(JJ,JST)==JP_RG_MF).OR.&
          (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_MF).OR.&
          (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_DH).OR.&
          (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_FC))) THEN

          ZCRUST_THICKNESS(JJ) = ZCRUST_THICKNESS(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
      ENDIF

      !
      !If no slab structure was found yet (GBELOW_SLAB = to be under a slab)
      IF(.NOT.GBELOW_SLAB(JJ)) THEN

        !Update of GPREV_ISNOT_MF(JJ) = True if layer just above is MF-like snow
        IF((PSNOWTYPEMEPRA(JJ,JST)==JP_RG_MF).OR.&
           (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_MF).OR.&
           (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_DH).OR.&
           (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_FC)) THEN

          GPREV_ISNOT_MF(JJ) = .FALSE.

        !Update of GBELOW_SLAB(JJ) = True if current layer is below a slab
        ELSEIF((PSNOWSHEAR(JJ,JST).GT.XACC_SLA_STR).AND.&
               (GPREV_ISNOT_MF(JJ))                .AND.&
               ((PSNOWTYPEMEPRA(JJ,JST)==JP_DF_RG).OR.&
                (PSNOWTYPEMEPRA(JJ,JST)==JP_RG_RG).OR.&
                (PSNOWTYPEMEPRA(JJ,JST)==JP_RG_FC).OR.&
                (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_DF))) THEN

          GBELOW_SLAB(JJ) = .TRUE.

        ELSE
          !Update of GPREV_ISNOT_MF(JJ) = True if layer just above is MF-like snow
          GPREV_ISNOT_MF(JJ) = .TRUE.
        ENDIF

      !If a slab structure was found and but no high instability -> searching for a weak layer
      ELSEIF ((ZSNOW_DEPTH_TOP>= 0.01).AND.&
              (ZSNOW_DEPTH_TOP<  1.00).AND.&
              (ZACC_HIG_DEP(JJ) == XUNDEF)) THEN

        !Permanent weak layer composed of facets or depth hoar
        IF((PSNOWTYPEMEPRA(JJ,JST)==JP_FC_FC).OR.&
           (PSNOWTYPEMEPRA(JJ,JST)==JP_FC_DH).OR.&
           (PSNOWTYPEMEPRA(JJ,JST)==JP_DH_DH)) THEN

          !If skier_ratio < threshold1
          IF((PACC_RAT(JJ,JST).LT.XACC_RAT_HIG)) THEN
            ZACC_HIG_DEP(JJ) = ZSNOW_DEPTH(JJ) - 0.5 * PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)

            !Condition on the crust thickness above weak layer
            IF(ZCRUST_THICKNESS(JJ) <= 0.01) THEN
              PACC_LEV(JJ) = JPACC_HIG
            ELSE
              PACC_LEV(JJ) = JPACC_MOD
            ENDIF

          !Else if no moderate was found yet and skier_ratio < XACC_RAT_MOD
          ELSEIF((ZACC_MOD_DEP(JJ) == XUNDEF).AND.(PACC_RAT(JJ,JST).LT.XACC_RAT_MOD)) THEN
            ZACC_MOD_DEP(JJ) = ZSNOW_DEPTH(JJ) - 0.5 * PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)

            !Condition on the crust thickness above weak layer
            IF(ZCRUST_THICKNESS(JJ) <= 0.01) THEN
              PACC_LEV(JJ) = JPACC_MOD
            ELSE
              PACC_LEV(JJ) = JPACC_LOW
            ENDIF

          ENDIF

        !Else if not already found looking for a
        !temporary weak layer composed of precipitation particles or decomposed snow
        !Note that there are no condition on ratio_acc
        ELSEIF((.NOT.GACC_MOD_TEMP(JJ)).AND.&
               ((PSNOWTYPEMEPRA(JJ,JST)==JP_PP_PP).OR.&
                (PSNOWTYPEMEPRA(JJ,JST)==JP_PP_DF).OR.&
                (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_DF))) THEN

          ZACC_MOD_DEP(JJ) = ZSNOW_DEPTH(JJ) - 0.5 * PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
          GACC_MOD_TEMP(JJ) = .TRUE.

          !Condition on the crust thickness above weak layer
          IF(ZCRUST_THICKNESS(JJ) <= 0.01) THEN
            PACC_LEV(JJ) = JPACC_MOD
          ELSE
            PACC_LEV(JJ) = JPACC_LOW
          ENDIF

        ENDIF
      ENDIF


      !The accidental risk is later (see below) combined with the natural (spontaneous) risk.
!
!
!
      !#######################Superior profile detection and classification########################
      !
      ! Note that there are more superior profile classes in original Mepra but they are never used
      ! in practice. In consequence, they are not considered here. Three types of sup. profile are
      ! considered: NEW (new snow), FRO (top frozen), WET (top wet).
      !
      ! Profile type NEW
      ! The superior profile, in case of the presence recent snow, is composed of
      ! the recent snow layers (even not at surface) and all above layers and possibly directly
      ! below humid layers (thermal_state>2), where thin ice crusts (connected thickness < 1 cm and
      ! composed of dry MF-like snow layers) are "ignored".

      ! Profile types WET or FRO
      ! Roughly, these profiles are composed of connected layers of MF-like snow if the top of
      ! this block of MF layers is close to the snow surface. There are additional constraints to
      ! define these profiles but they do not make clear sense to me.

      !!!!TODO!!!!!
      !PPRO_SUP_TYP when no snow = JPPRO_SUP_NAN = 6
      !HSUP HINF in case of flat terrain ???

      ! GMF = the snowtype is MF-like (MF/RG, MF, MF/DH, MF/FC)
      GMF = (PSNOWTYPEMEPRA(JJ,JST)==JP_RG_MF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_MF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_DH).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_FC)

      ! GPP = the snowtype is new-like, i.e. being of type PP, PP/DF, DF, DF/RG, DF/FC
      ! not strictly equivalent of being dendritic
      GPP = (PSNOWTYPEMEPRA(JJ,JST)==JP_PP_PP).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_PP_DF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_DF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_RG).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_FC)

      IF(GPP) THEN
        ! When recent snow is found, then
        ! the superior profile type is, in all cases, NEW,
        ! the superior profile reaches at least this layer,
        ! what is direclty below potentially belongs to the superior profile and
        ! the thickness of the connex crusts is re-initiliazed to zero.
        PPRO_SUP_TYP (JJ) = JPPRO_SUP_NEW
        IPRO_SUP_LIM (JJ) = JST
        PDEP_SUP     (JJ) = ZSNOW_DEPTH(JJ)
        GCOULD_BE_NEW(JJ) = .TRUE.
        ZPRO_CRUST   (JJ) = 0

      ELSEIF(GCOULD_BE_NEW(JJ)) THEN
        ! if the current layer potentially belongs to the sup. profile

        IF((.NOT.GTHERMSTATE).AND.(ZPRO_CRUST(JJ) < XPRO_SUP_CRU)) THEN
          ! if the current layer is humid and the above crust thickness < 1 cm
          ! then the current layer and potentially what is below, belongs to the sup. profile and
          ! the crust thickness is re-initialized to zero.
          IPRO_SUP_LIM(JJ) = JST
          PDEP_SUP    (JJ) = ZSNOW_DEPTH(JJ)
          ZPRO_CRUST  (JJ) = 0

        ELSEIF(GMF.AND.GTHERMSTATE) THEN
          ! elseif the current layer is a crust
          ! then we increase the crust thickness accordingly
          ZPRO_CRUST(JJ) = ZPRO_CRUST(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)

        ELSE
          ! in all other cases, the layer and what is below do not belong to the sup. profile.
          GCOULD_BE_NEW(JJ) = .FALSE.

        ENDIF

      ELSEIF(PPRO_SUP_TYP(JJ).NE.JPPRO_SUP_NEW) THEN
      ! if a sup. profile of type NEW was not found yet.
      ! COULD BE IMPROVED FOR CLARITY ...

        IF(.NOT.GFOUND(JJ)) THEN
          ! if not found yet, we are looking for the first layer of type MF or being at a
          ! depth > 3 cm. This layer thermal state determines whether the profile is FRO or WET
          ! (not consistent...). In case of this layer being not MF, we are not sure yet that the
          ! sup. profile will be FRO/WET or NAN.

          IF(GMF.OR.(ZSNOW_DEPTH(JJ)>XPRO_SUP_DEP)) THEN
            GFOUND      (JJ) = .TRUE.
            IPRO_SUP_LIM(JJ) = JST
            PDEP_SUP    (JJ) = ZSNOW_DEPTH(JJ)

            IF(GTHERMSTATE) THEN
              PPRO_SUP_TYP(JJ) = JPPRO_SUP_FRO
            ELSE
              PPRO_SUP_TYP(JJ) = JPPRO_SUP_WET
            ENDIF

            GIKNOTMF(JJ) = .NOT.GMF

            IF(GIKNOTMF(JJ)) THEN
              PPRO_SUP_TYP(JJ) = JPPRO_SUP_NAN
            ENDIF

          ENDIF
        ENDIF

        IF((GFOUND(JJ)).AND.(.NOT.GIKNOTMF(JJ))) THEN
          IF(GMF) THEN
!            !if the current layer is of type MF
!            !then we increase the sup. profile and
!            !we do not care anymore on the type of the "found" layer
            IPRO_SUP_LIM(JJ) = JST
            PDEP_SUP    (JJ) = ZSNOW_DEPTH(JJ)


          ELSE
            GIKNOTMF(JJ)=.TRUE.

          ENDIF
        ENDIF
      ENDIF

!
!
      !!!!!!!!!!!!!!!!!!!!!Cumulative SWE and DEPTH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !Compute depth and SWE of snow with age < X days (X in {1,3,5,7})
      !
      !WARNING: there is some rounding, age could be expressed as int in seconds or minutes
      !but there would be a loss of genericity of the possible crocus timesetps
      !to get equality with snow_pro.py, replace age<duration by age<=duration in python

      ! PSNOWSWE is slope perpendicular, it is here porjected on the vertical axis,
      ! as it is later done in diag_misc_isban.
      ZSNOWSWE = PSNOWSWE(JJ,JST) / PDIRCOSZW(JJ)
      ZSNOWDZ  = PSNOWDZ (JJ,JST) / PDIRCOSZW(JJ)

      IF(    PSNOWAGE(JJ,JST) <= 7) THEN
        PSNOWDEPTH_7DAYS(JJ) = PSNOWDEPTH_7DAYS(JJ) + ZSNOWDZ
        PSNOWSWE_7DAYS  (JJ) = PSNOWSWE_7DAYS  (JJ) + ZSNOWSWE
      ELSEIF(PSNOWAGE(JJ,JST) <= 5) THEN
        PSNOWDEPTH_5DAYS(JJ) = PSNOWDEPTH_5DAYS(JJ) + ZSNOWDZ
        PSNOWSWE_5DAYS  (JJ) = PSNOWSWE_5DAYS  (JJ) + ZSNOWSWE
      ELSEIF(PSNOWAGE(JJ,JST) <= 3) THEN
        PSNOWDEPTH_3DAYS(JJ) = PSNOWDEPTH_3DAYS(JJ) + ZSNOWDZ
        PSNOWSWE_3DAYS  (JJ) = PSNOWSWE_3DAYS  (JJ) + ZSNOWSWE
      ELSEIF(PSNOWAGE(JJ,JST) <= 1) THEN
        PSNOWDEPTH_1DAYS(JJ) = PSNOWDEPTH_1DAYS(JJ) + ZSNOWDZ
        PSNOWSWE_1DAYS  (JJ) = PSNOWSWE_1DAYS  (JJ) + ZSNOWSWE
      ENDIF

      ! Ramsonde top penetration
      IF ((GRAM(JJ)).AND.(PSNOWRAM(JJ,JST)<=2.)) THEN
        PSNOWRAM_SONDE(JJ)   = PSNOWRAM_SONDE(  JJ) + ZSNOWDZ
      ELSE
        GRAM(JJ)=.FALSE.
      ENDIF

      ! Depth of top wet snow
      IF ((GWET(JJ)).AND.(PSNOWLIQ(JJ,JST)>0)) THEN
        PSNOW_WETTHICKNESS(JJ) = PSNOW_WETTHICKNESS(JJ) + PSNOWDZ(JJ,JST)
      ELSE
        GWET(JJ)=.FALSE.
      ENDIF

      ! Depth of top refrozen snow
      IF (GREFROZEN(JJ).AND.(PSNOWHIST(JJ,JST)>=2).AND.(PSNOWTEMP(JJ,JST)<273.15)) THEN
        PSNOW_REFTHICKNESS(JJ) = PSNOW_REFTHICKNESS(JJ) + PSNOWDZ(JJ,JST)
      ELSE
        GREFROZEN(JJ)=.FALSE.
      ENDIF


      ! Specific surface area
      ! in snowpro the density of ice (XRHOLI) is 900 kg/m3 instead of 917 kg/m3 (SURFEX)
      ! which needs to be changed in snow_pro.py
      ! Note: only tested with B92
      IF ( HSNOWMETAMO=='B92' ) THEN
        PSNOWSSA(JJ,JST) = 6. / (XRHOLI*ZDIAM)
      ELSE
        PSNOWSSA(JJ,JST) = 6. / (XRHOLI*PSNOWGRAN1(JJ,JST))
      END IF

    ENDIF
  END DO
END DO

!RE-initialization for new calculations
DO JJ=1,SIZE(PSNOWSWE,1)
  PDEP_TOT   (JJ) = ZSNOW_DEPTH(JJ)
  ZSNOW_DEPTH(JJ) = 0
  ZDEP_INF   (JJ) = PDEP_TOT(JJ) - PDEP_SUP(JJ)
  GWET       (JJ) = .TRUE.
END DO
!
!
DO JST=1,SIZE(PSNOWSWE,2)
  DO JJ=1,SIZE(PSNOWSWE,1)
!
    IF (PSNOWSWE(JJ,JST) > 0) THEN
      ! Do something only in case of non-empty layer
      !
      ! Update of snow depth (depth of bottom of current layer)
      !!!!!!!!!!!!!!!!!ALREADY CALCULATED
      ZSNOW_HEIGHT_TOP = PDEP_TOT(JJ) - ZSNOW_DEPTH(JJ)
      ZSNOW_DEPTH(JJ)  = ZSNOW_DEPTH(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)

!
      IF(JST<=IPRO_SUP_LIM(JJ)) THEN
        ! Inside superior profile
!
        ZSCW = XRHOLW * PSNOWLIQ(JJ,JST) / (PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ)) !liquid water content
        GTHERMSTATE = (PSNOWTEMP(JJ,JST) < 272.96).OR.(ZSCW < 5)!True if thermstate <= 2
!
        ! Calculations of the height of the uppest continuous block of humid snow in the sup.
        ! profile (slightly different from SNOWWETTHICKNESS)
        IF(GWET(JJ).AND.(.NOT.GTHERMSTATE)) THEN
          PDEP_HUM(JJ) = PDEP_HUM(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
        ELSEIF(PDEP_HUM(JJ).GT.0) THEN
          GWET(JJ) = .FALSE.
        ENDIF
!
        ! Determination of natural risk based on the stress/strength ratio
        IF(ZSNOW_HEIGHT_TOP.GE.XNAT_HEI_MIN) THEN
        !Only searching above a certain height

          IF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NEW) THEN
          ! Case of NEW sup. profile

            IF(PNAT_RAT(JJ,JST).LE.XNAT_RAT_HIG) THEN
              ZNAT_HIG_DEP(JJ) = ZSNOW_DEPTH(JJ)
            ELSEIF(PNAT_RAT(JJ,JST).LE.(XNAT_RAT_MOD + 0.05)) THEN !verrue d'origine inconnue
              ZNAT_MOD_DEP(JJ) = ZSNOW_DEPTH(JJ)
            ENDIF

          ELSE
          ! Case of WET or FRO profiles (code not accesible for nan profile)

            IF(PNAT_RAT(JJ,JST).LE.XNAT_RAT_HIG) THEN
              ZNAT_HIG_DEP(JJ) = ZSNOW_DEPTH(JJ)
            ELSEIF(PNAT_RAT(JJ,JST).LE.XNAT_RAT_MOD) THEN
              ZNAT_MOD_DEP(JJ) = ZSNOW_DEPTH(JJ)
            ENDIF


          ENDIF


          ! Calculations used for avalanche type determination

          ! GHUM and GDRY used only for profile sup. NEW
          IF(GTHERMSTATE) THEN
            GHUM(JJ) = .FALSE.
          ELSE
            GDRY(JJ) = .FALSE.
          ENDIF

          ! GMEL_GRO used only for profiles sup. WET and NAN
          IF((ZSNOW_DEPTH(JJ).GT.XNAT_HEI_MIN).AND.GTHERMSTATE) THEN
            GMEL_GRO(JJ) = .FALSE.
          ENDIF

          !ZHUMTHICK used for profiles sup. FRO and NAN
          !Non-sense ...
          ZHUMTHICK(JJ) = ZHUMTHICK(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)



        ENDIF
!
      ELSE
        !inside inferior profile to determine its type
        IF(PPRO_INF_TYP(JJ).NE.JPPRO_INF_HAR) THEN
          IF(ZSNOW_HEIGHT_TOP.GT.(XPRO_INF_COE * ZDEP_INF(JJ))) THEN
          !NOOOOOOOOOOOOOOT SURE of top
            IF(PSNOWRAM(JJ,JST).LT.XPRO_INF_RAM) THEN
              PPRO_INF_TYP(JJ) = JPPRO_INF_SOF
            ELSE
              PPRO_INF_TYP(JJ) = JPPRO_INF_HAR
            ENDIF
          ENDIF
        ENDIF

      ENDIF

    ! Storing the current GTHERMSTATE to get the one of the botom layer
    GTHERMSTATE_BOT(JJ) = GTHERMSTATE

    ! Storing the current layer thickness to get the one of the botom layer
    ZTHICK_BOT(JJ) = PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)

    ENDIF

  ENDDO
ENDDO

!Loop only on points
DO JJ=1,SIZE(PSNOWSWE,1)

  !!!!!!!!Natural risks index determination!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SELECT CASE(INT(PPRO_SUP_TYP(JJ)))
    CASE(JPPRO_SUP_NEW)
      IPRO_CLASS = 0
    CASE(JPPRO_SUP_WET)
      IPRO_CLASS = 15
    CASE(JPPRO_SUP_FRO)
      IPRO_CLASS = 30
    CASE DEFAULT
      IPRO_CLASS = 0 !CHECKKKKKKKKKKKKKK
  END SELECT

  IF    (ZNAT_HIG_DEP(JJ) > XNAT_HEI_HIG) THEN
    IVAL_HIG = JPNAT_TAB(15 + IPRO_CLASS)
  ELSEIF(ZNAT_HIG_DEP(JJ) > XNAT_HEI_MOD) THEN
    IVAL_HIG = JPNAT_TAB(12 + IPRO_CLASS)
  ELSEIF(ZNAT_HIG_DEP(JJ) > XNAT_HEI_LOW) THEN
    IVAL_HIG = JPNAT_TAB(9  + IPRO_CLASS)
  ELSEIF(ZNAT_HIG_DEP(JJ) > 0           ) THEN
    IVAL_HIG = JPNAT_TAB(6  + IPRO_CLASS)
  ELSE
    IVAL_HIG = JPNAT_TAB(3  + IPRO_CLASS)
  ENDIF

  IF    (ZNAT_MOD_DEP(JJ) > XNAT_HEI_HIG) THEN
    IVAL_MOD = JPNAT_TAB(14 + IPRO_CLASS)
  ELSEIF(ZNAT_MOD_DEP(JJ) > XNAT_HEI_MOD) THEN
    IVAL_MOD = JPNAT_TAB(11 + IPRO_CLASS)
  ELSEIF(ZNAT_MOD_DEP(JJ) > XNAT_HEI_LOW) THEN
    IVAL_MOD = JPNAT_TAB(8  + IPRO_CLASS)
  ELSEIF(ZNAT_MOD_DEP(JJ) > 0)            THEN
    IVAL_MOD = JPNAT_TAB(5  + IPRO_CLASS)
  ELSE
    IVAL_MOD = JPNAT_TAB(2  + IPRO_CLASS)
  ENDIF

  IF    (PDEP_SUP(JJ) > XNAT_HEI_HIG) THEN
    IVAL_LOW = JPNAT_TAB(13 + IPRO_CLASS)
  ELSEIF(PDEP_SUP(JJ) > XNAT_HEI_MOD) THEN
    IVAL_LOW = JPNAT_TAB(10 + IPRO_CLASS)
  ELSEIF(PDEP_SUP(JJ) > XNAT_HEI_LOW) THEN
    IVAL_LOW = JPNAT_TAB(7  + IPRO_CLASS)
  ELSEIF(PDEP_SUP(JJ) > 0)            THEN
    IVAL_LOW = JPNAT_TAB(4  + IPRO_CLASS)
  ELSE
    IVAL_LOW = JPNAT_TAB(1  + IPRO_CLASS)
  ENDIF

  IF(IVAL_HIG.NE.JPNAT_MOA) THEN
    PNAT_LEV(JJ) = MAX(INT(IVAL_HIG),MAX(INT(IVAL_MOD),INT(IVAL_LOW)))
  ELSE
    PNAT_LEV(JJ) = JPNAT_MOA
  ENDIF

  !!!!!!!!!!!!!!!!!Avalanche type!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  !Does not make sense, whatever...
  ZHUMTHICK(JJ) = ZHUMTHICK(JJ) + ZTHICK_BOT(JJ)

  IF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NEW) THEN
    !Profile sup. NEW

    IF(GHUM(JJ)) THEN
      PAVA_TYP(JJ) = JPAVA_NEW_WET
    ELSEIF(GDRY(JJ)) THEN
      PAVA_TYP(JJ) = JPAVA_NEW_DRY
    ELSE
      PAVA_TYP(JJ) = JPAVA_NEW_MIX
    ENDIF

  ELSEIF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NAN) THEN
    !Profile sup. NAN
    PAVA_TYP(JJ) = JPAVA_NAN

  ELSEIF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_WET) THEN
    !Profile sup. WET

    IF(PPRO_INF_TYP(JJ).EQ.JPPRO_INF_NAN) THEN
      !Profile inf NAN
      IF(GMEL_GRO(JJ)) THEN
        PAVA_TYP(JJ) = JPAVA_MEL_GRO
      ELSE
        PAVA_TYP(JJ) = JPAVA_MEL_SUR
      ENDIF

    ELSEIF(PPRO_INF_TYP(JJ).EQ.JPPRO_INF_HAR) THEN
      !Profile inf HAR
      PAVA_TYP(JJ) = JPAVA_MEL_SUR
    ELSE
      !Profile inf SOF
      PAVA_TYP(JJ) = JPAVA_MEL_GRO
    ENDIF

  ELSEIF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_FRO) THEN
    !Profile sup. WET
    IF(PPRO_INF_TYP(JJ).EQ.JPPRO_INF_NAN) THEN
      IF((.NOT.GTHERMSTATE_BOT(JJ)).AND.(ZHUMTHICK(JJ).GT.PDEP_TOT(JJ))) THEN
        PAVA_TYP(JJ) = JPAVA_MEL_SUR
      ELSE
        PAVA_TYP(JJ) = JPAVA_MEL_GRO
      ENDIF
    ELSEIF(PPRO_INF_TYP(JJ).EQ.JPPRO_INF_HAR) THEN
      !Profile inf HAR
      PAVA_TYP(JJ) = JPAVA_MEL_SUR
    ELSE
      !Profile inf SOF
      PAVA_TYP(JJ) = JPAVA_MEL_GRO
    ENDIF

    IF(PNAT_LEV(JJ).EQ.JPNAT_VLO) PAVA_TYP(JJ) = JPAVA_NAN

  ENDIF


  !!!!!Risk nat actualize
  ZNAT_LEV_TMP = PNAT_LEV(JJ)
  ! only in case, there is not a NAN risk in the previous step
  IF(ZNAT_LEV_PRE(JJ).NE.JPNAT_NAN) THEN

    !For type sup. NEW
    IF((  PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NEW   )   .AND.&
       (  PAVA_TYP    (JJ).EQ.ZAVA_TYP_PRE(JJ))   .AND.&
       (  PDEP_TOT    (JJ).LT.ZDEP_TOT_PRE(JJ))   .AND.&
       (  PDEP_SUP    (JJ).LE.ZDEP_SUP_PRE(JJ))   .AND.&
       (((PAVA_TYP(JJ).EQ.JPAVA_NEW_MIX).AND.(PDEP_HUM(JJ).LE.ZDEP_HUM_PRE(JJ))).OR.&
        (PAVA_TYP(JJ).NE.JPAVA_NEW_MIX))) THEN

      ZNAT_LEV_TMP = JPNAT_ACT(INT(PNAT_LEV(JJ) + ZNAT_LEV_PRE(JJ)*7 + 1))
      !
      !Weird addtional instruction
      IF(( PAVA_TYP(JJ).NE.JPAVA_NEW_MIX).AND.&
         ((PNAT_LEV(JJ).EQ.JPNAT_HIG).OR.(PNAT_LEV(JJ).EQ.JPNAT_VHI)).AND.&
         ((ZNAT_LEV_PRE(JJ).EQ.JPNAT_MOD).OR.(ZNAT_LEV_PRE(JJ).EQ.JPNAT_HIG).OR.&
          (ZNAT_LEV_PRE(JJ).EQ.JPNAT_VHI))) THEN
         ZNAT_LEV_TMP = JPNAT_MOD
      ENDIF

    ENDIF

    IF((PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NEW   )   .AND.&
       (PAVA_TYP    (JJ).EQ.ZAVA_TYP_PRE(JJ))   .AND.&
       (PDEP_TOT    (JJ).LT.ZDEP_TOT_PRE(JJ))   .AND.&
       (PNAT_LEV    (JJ).EQ.JPNAT_MOA)) THEN

        IF(ZNAT_LEV_PRE(JJ).EQ.JPNAT_HIG) ZNAT_LEV_TMP = JPNAT_MOD
        IF(ZNAT_LEV_PRE(JJ).EQ.JPNAT_VHI) ZNAT_LEV_TMP = JPNAT_HIG

    ENDIF

    !PACC_LEV(JJ) = PNAT_LEV(JJ) only for testing
    PNAT_LEV(JJ) = ZNAT_LEV_TMP

!
    !For type sup. WET or FRO
    IF(((PPRO_SUP_TYP    (JJ).EQ.JPPRO_SUP_WET).OR.(PPRO_SUP_TYP    (JJ).EQ.JPPRO_SUP_FRO)).AND.&
       ((ZPRO_SUP_TYP_PRE(JJ).EQ.JPPRO_SUP_WET).OR.(ZPRO_SUP_TYP_PRE(JJ).EQ.JPPRO_SUP_FRO)).AND.&
       (ZDEP_INF(JJ).GE.(ZDEP_TOT_PRE(JJ)-ZDEP_SUP_PRE(JJ)-0.05))) THEN

      IF(PNAT_LEV(JJ).EQ.JPNAT_MOD) THEN
        PNAT_LEV(JJ) = JPNAT_LOW
      ENDIF

      IF((PPRO_SUP_TYP    (JJ).EQ.JPPRO_SUP_WET).AND.&
         ((PNAT_LEV(JJ).EQ.JPNAT_HIG).OR.(PNAT_LEV(JJ).EQ.JPNAT_VHI))) THEN

        IF((ZNAT_LEV_PRE(JJ).EQ.JPNAT_MOD).OR.&
           (ZNAT_LEV_PRE(JJ).EQ.JPNAT_HIG).OR.&
           (ZNAT_LEV_PRE(JJ).EQ.JPNAT_VHI)) THEN
          PNAT_LEV(JJ) = JPNAT_MOD
        ELSEIF(ZNAT_LEV_PRE(JJ).EQ.JPNAT_LOW) THEN
          PNAT_LEV(JJ) = JPNAT_LOW
        ENDIF
      ENDIF
    ENDIF

    !For no profil sup. (JPPRO_SUP_NAN)
    IF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NAN) THEN
      PNAT_LEV(JJ) = JPNAT_VLO
      !PAVA_TYP(JJ) = JPAVA_MEL_SUR COMPREND PAS D OU CA VIENT?
    ENDIF

  ENDIF

  !Case with no snow (most of the time we do not even enter snowcrodiag)
  IF(PDEP_TOT(JJ).LE.0.0000001) THEN
    PNAT_LEV(JJ) = JPNAT_NAN
  ENDIF

! Combination of natural and accidental risk levels
!
  IF    ((PNAT_LEV(JJ).EQ.JPNAT_VHI).OR.(PNAT_LEV(JJ).EQ.JPNAT_HIG)) THEN
    IACC_FROMNAT = JPACC_HIG
  ELSEIF((PNAT_LEV(JJ).EQ.JPNAT_MOD).OR.(PNAT_LEV(JJ).EQ.JPNAT_MOA)) THEN
    IACC_FROMNAT = JPACC_MOD
  ELSE
    IACC_FROMNAT = JPACC_LOW !!!!!!!!!!!!!!!!!!JP_ACC_LOW ???
  ENDIF

  IF(PDEP_TOT(JJ).LE.XNAT_HEI_LOW) THEN
    PACC_LEV(JJ) = JPACC_LOW
  ENDIF


  !Accidental risk
  IF(ZACC_HIG_DEP(JJ) == XUNDEF) THEN
    IF(ZACC_MOD_DEP(JJ) == XUNDEF) THEN
      ZACC_HIG_DEP(JJ) = -1
      ZACC_MOD_DEP(JJ) = -1
    ELSE
      ZACC_HIG_DEP(JJ) = -1
      ZACC_MOD_DEP(JJ) = ZACC_MOD_DEP(JJ)
    ENDIF
  ELSE
    ZACC_HIG_DEP(JJ) = ZACC_HIG_DEP(JJ)
    ZACC_MOD_DEP(JJ) = -1
  ENDIF

  IF((IACC_FROMNAT.GE.PACC_LEV(JJ)).OR.(PDIRCOSZW(JJ).GT.0.8)) THEN
    PACC_LEV(JJ) = IACC_FROMNAT
    PDEP_HIG(JJ) = ZNAT_HIG_DEP(JJ)
    PDEP_MOD(JJ) = ZNAT_MOD_DEP(JJ)
  ELSE
    PDEP_HIG(JJ) = ZACC_HIG_DEP(JJ)
    PDEP_MOD(JJ) = ZACC_MOD_DEP(JJ)
  ENDIF

  IF((PDIRCOSZW(JJ).GT.0.8).AND.(PACC_LEV(JJ).EQ.JPACC_LOW)) THEN
    PACC_LEV(JJ) = JPACC_NUL
  ENDIF

  !PDEP_HIG(JJ) = ZNAT_HIG_DEP(JJ) !Only for testing
  !PDEP_MOD(JJ) = ZNAT_MOD_DEP(JJ) !Only for testing

  !!!!!!!!!!!!!!!!!!!!
  !Verrue additionelle pour le cas sans pente
  !Does not make sense, since these variables could be also defined for flat terrain
  IF(PDIRCOSZW(JJ).GT.0.99) THEN
  !!!!!!!!!!!!!!??replace 0 by XUNDEF
    PDEP_SUP(JJ) = 0
    PDEP_HUM(JJ) = 0
    PDEP_HIG(JJ) = XUNDEF
    PDEP_MOD(JJ) = XUNDEF
    PAVA_TYP(JJ) = 0
    PNAT_LEV(JJ) = JPNAT_VLO
    PACC_LEV(JJ) = JPACC_NUL

  ENDIF

  IF(PDEP_HIG(JJ).EQ.0) THEN
    PDEP_HIG(JJ) = XUNDEF
  ENDIF

  IF(PDEP_MOD(JJ).EQ.0) THEN
    PDEP_MOD(JJ) = XUNDEF
  ENDIF


ENDDO


END SUBROUTINE SNOWCRO_DIAG
