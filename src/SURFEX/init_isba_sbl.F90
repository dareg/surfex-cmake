!     #########
    SUBROUTINE INIT_ISBA_SBL(IO, IP, IMX, IMT, IR, ICP, PTSTEP, PPA, PPS, PTA, PQA, PRHOA, PU, PV, &
                             PDIR_SW, PSCA_SW, PSW_BANDS, PRAIN, PSNOW, PZREF, PUREF, PSSO_SLOPE )  
!     #################################################################################
!
!!****  *INIT_WATER_SBL* - inits water SBL profiles
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     S. Riette
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    03/2010
!!------------------------------------------------------------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_CANOPY_n, ONLY : CANOPY_t
!
USE MODD_TYPE_SNOW
!
USE MODD_CSTS,             ONLY : XCPD, XRD, XP00, XG, XLVTT
USE MODD_SURF_ATM,         ONLY : LNOSOF
USE MODD_CANOPY_TURB,      ONLY : XALPSBL
!
USE MODI_CLS_TQ
USE MODI_ISBA_SNOW_FRAC
USE MODI_WET_LEAVES_FRAC
USE MODI_VEG
USE MODI_DRAG
USE MODI_CLS_WIND
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: IMT
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(CANOPY_t), INTENT(INOUT) :: ICP
!
REAL,               INTENT(IN)   :: PTSTEP   ! timestep of the integration
REAL, DIMENSION(:), INTENT(IN)  :: PPA       ! pressure at forcing level             (Pa)
REAL, DIMENSION(:), INTENT(IN)  :: PPS       ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(:), INTENT(IN)  :: PTA       ! air temperature forcing               (K)
REAL, DIMENSION(:), INTENT(IN)  :: PQA       ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(:), INTENT(IN)  :: PRHOA     ! air density                           (kg/m3)
REAL, DIMENSION(:), INTENT(IN)  :: PU        ! zonal wind                            (m/s)
REAL, DIMENSION(:), INTENT(IN)  :: PV        ! meridian wind                         (m/s)
REAL, DIMENSION(:,:),INTENT(IN) :: PDIR_SW   ! direct  solar radiation (on horizontal surf.)
!                                            !                                       (W/m2)
REAL, DIMENSION(:,:),INTENT(IN) :: PSCA_SW   ! diffuse solar radiation (on horizontal surf.)
!                                            !                                        (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(:), INTENT(IN)  :: PSNOW     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)  :: PRAIN     ! liquid precipitation                  (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)  :: PZREF     ! height of T,q forcing                 (m)
REAL, DIMENSION(:), INTENT(IN)  :: PUREF     ! height of wind forcing                (m)
REAL, DIMENSION(:), INTENT(IN)  :: PSSO_SLOPE! slope of S.S.O.                         (-)
!                                            ! canopy       
!
!*      0.2    declarations of local variables
!
!* forcing variables
!
REAL, DIMENSION(SIZE(PTA))   :: ZWIND    ! lowest atmospheric level wind speed           (m/s)
REAL, DIMENSION(SIZE(PTA))   :: ZEXNA    ! Exner function at lowest SBL scheme level     (-)
REAL, DIMENSION(SIZE(PTA))   :: ZQA      ! specific humidity                             (kg/m3)
!
! SBL turbulence scheme
!
REAL, DIMENSION(SIZE(PTA))   ::ZRI
REAL, DIMENSION(SIZE(PTA))   ::ZCD
REAL, DIMENSION(SIZE(PTA))   ::ZCDN
REAL, DIMENSION(SIZE(PTA))   ::ZCH
REAL, DIMENSION(SIZE(PTA))   ::ZTNM
REAL, DIMENSION(SIZE(PTA))   ::ZQNM
REAL, DIMENSION(SIZE(PTA))   ::ZHUNM
REAL, DIMENSION(SIZE(PTA))   ::ZP_SLOPE_COS
REAL, DIMENSION(SIZE(PTA))   ::ZZ0
REAL, DIMENSION(SIZE(PTA))   ::ZZ0H
REAL, DIMENSION(SIZE(PTA))   ::ZEXNS
REAL, DIMENSION(SIZE(PTA))   ::ZTS
REAL, DIMENSION(SIZE(PTA))   ::ZHU
REAL, DIMENSION(SIZE(PTA))   ::ZQS
REAL, DIMENSION(SIZE(PTA))   ::ZZ0EFF
REAL, DIMENSION(SIZE(PTA))   ::ZWG
REAL, DIMENSION(SIZE(PTA))   ::ZWGI
REAL, DIMENSION(SIZE(PTA))   ::ZVEG
REAL, DIMENSION(SIZE(PTA))   ::ZRESA
REAL, DIMENSION(SIZE(PTA))   ::ZHUG
REAL, DIMENSION(SIZE(PTA))   ::ZHUGI
REAL, DIMENSION(SIZE(PTA))   ::ZHV
REAL, DIMENSION(SIZE(PTA))   ::ZCPS
REAL, DIMENSION(SIZE(PTA))   ::ZWRMAX_CF
REAL, DIMENSION(SIZE(PTA))   ::ZWR
REAL, DIMENSION(SIZE(PTA))   ::ZZ0_WITH_SNOW
REAL, DIMENSION(SIZE(PTA))   ::ZPSNG
REAL, DIMENSION(SIZE(PTA))   ::ZPSNV
REAL, DIMENSION(SIZE(PTA))   ::ZPSNV_A
REAL, DIMENSION(SIZE(PTA))   ::ZPSN
REAL, DIMENSION(SIZE(PTA))   ::ZSNOWALB
REAL, DIMENSION(SIZE(PTA),SIZE(IR%TSNOW%WSNOW,2)) ::ZSNOWSWE
REAL, DIMENSION(SIZE(PTA),SIZE(IR%TSNOW%WSNOW,2)) ::ZSNOWRHO
REAL, DIMENSION(SIZE(PTA))   ::ZFFG
REAL, DIMENSION(SIZE(PTA))   ::ZFFGNOS
REAL, DIMENSION(SIZE(PTA))   ::ZFFV
REAL, DIMENSION(SIZE(PTA))   ::ZFFVNOS
REAL, DIMENSION(SIZE(PTA))   ::ZFF
REAL, DIMENSION(SIZE(PTA))   ::ZRS
REAL, DIMENSION(SIZE(PTA))   ::ZP_GLOBAL_SW
REAL, DIMENSION(SIZE(PTA))   ::ZF2
REAL, DIMENSION(SIZE(PTA))   ::ZF5
REAL, DIMENSION(SIZE(PTA))   ::ZLAI
REAL, DIMENSION(SIZE(PTA))   ::ZGAMMA
REAL, DIMENSION(SIZE(PTA))   ::ZRGL
REAL, DIMENSION(SIZE(PTA))   ::ZRSMIN
REAL, DIMENSION(SIZE(PTA))   ::ZDELTA
REAL, DIMENSION(SIZE(PTA))   ::ZWRMAX
REAL, DIMENSION(SIZE(PTA))   ::ZCLS_WIND_ZON
REAL, DIMENSION(SIZE(PTA))   ::ZCLS_WIND_MER
REAL, DIMENSION(SIZE(PTA),SIZE(IR%TSNOW%WSNOW,2))   ::ZSUM_LAYER
REAL, DIMENSION(SIZE(PTA))   ::ZSUM
REAL, DIMENSION(SIZE(PTA))   :: ZLEG_DELTA  ! soil evaporation delta fn
REAL, DIMENSION(SIZE(PTA))   :: ZLEGI_DELTA ! soil sublimation delta fn
REAL, DIMENSION(SIZE(PTA))   :: ZLVTT
!
INTEGER                     :: JSWB
INTEGER                     :: JLAYER
INTEGER                     :: JPATCH
!
REAL, DIMENSION(SIZE(PTA),SIZE(IP%XPATCH,2)) ::ZWSNOW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_ISBA_SBL',0,ZHOOK_HANDLE)
!    
!Means over patches
ZTS     = SUM(IR%XTG (:,1,:)*IP%XPATCH(:,:),DIM=2)
ZWG     = SUM(IR%XWG (:,1,:)*IP%XPATCH(:,:),DIM=2)
ZWGI    = SUM(IR%XWGI(:,1,:)*IP%XPATCH(:,:),DIM=2)
ZZ0     = SUM(IP%XPATCH(:,:)*IMT%XZ0(:,:)  ,DIM=2)
!
!We choose to set ZZ0EFF and ZZ0_WITH_SNOW equal to ZZ0
ZZ0EFF        = ZZ0
ZZ0_WITH_SNOW = ZZ0
ZZ0H(:) = SUM(IP%XPATCH(:,:) * IMT%XZ0(:,:)/IMX%XZ0_O_Z0H(:,:),DIM=2)
ZVEG(:) = SUM(IP%XPATCH(:,:) * IMT%XVEG(:,:)                  ,DIM=2)
!
ZP_SLOPE_COS(:) = 1./SQRT(1.+PSSO_SLOPE(:)**2)
IF (LNOSOF) ZP_SLOPE_COS(:) = 1.0
!
ZRESA(:) = SUM(IP%XPATCH(:,:)*IR%XRESA(:,:),DIM=2)
WHERE(ZVEG(:)>0)
  ZLAI     (:)= SUM(IP%XPATCH(:,:)*IMT%XVEG(:,:)*IMT%XLAI(:,:)     ,&
                DIM=2,MASK=IMT%XVEG(:,:)>0) / ZVEG(:)
  ZWRMAX_CF(:)= SUM(IP%XPATCH(:,:)*IMT%XVEG(:,:)*IMT%XWRMAX_CF(:,:),&
                DIM=2,MASK=IMT%XVEG(:,:)>0) / ZVEG(:)
  ZWR      (:)= SUM(IP%XPATCH(:,:)*IMT%XVEG(:,:)*IR%XWR(:,:)       ,&
                DIM=2,MASK=IMT%XVEG(:,:)>0) / ZVEG(:)
ELSEWHERE
  ZLAI     (:) = IMT%XLAI     (:,1)
  ZWRMAX_CF(:) = IMT%XWRMAX_CF(:,1)
  ZWR      (:) = IR%XWR       (:,1)
ENDWHERE
!
ZSUM_LAYER(:,:) = 0.
ZSUM        (:) = 0.
!
DO JLAYER=1,SIZE(IR%TSNOW%WSNOW,2)
  ZSNOWSWE  (:,JLAYER) = SUM(IP%XPATCH(:,:)*IR%TSNOW%WSNOW(:,JLAYER,:),DIM=2)
  ZSUM_LAYER(:,JLAYER) = SUM(IP%XPATCH(:,:),DIM=2,MASK=IR%TSNOW%WSNOW(:,JLAYER,:)>0)
  WHERE(ZSUM_LAYER(:,JLAYER)>0)      
    ZSNOWRHO(:,JLAYER)= SUM( IP%XPATCH(:,:)*IR%TSNOW%RHO(:,JLAYER,:), DIM=2, &
                             MASK=IR%TSNOW%WSNOW(:,JLAYER,:)>0) / ZSUM_LAYER(:,JLAYER)
  ELSEWHERE
    ZSNOWRHO(:,JLAYER)=IR%TSNOW%RHO(:,JLAYER,1)
  ENDWHERE
END DO
!
ZSUM(:)=SUM(ZSUM_LAYER(:,:),DIM=2)
!
ZWSNOW(:,:) = 0.
DO JPATCH=1,SIZE(IR%TSNOW%WSNOW,3)
  DO JLAYER=1,SIZE(IR%TSNOW%WSNOW,2)
    ZWSNOW(:,JPATCH) = ZWSNOW(:,JPATCH) + IR%TSNOW%WSNOW(:,JLAYER,JPATCH)
  ENDDO
ENDDO    
!
WHERE(ZSUM(:)>0)         
  ZSNOWALB(:) = SUM(IP%XPATCH(:,:)*IR%TSNOW%ALB(:,:),DIM=2,MASK=ZWSNOW(:,:)>0) / ZSUM(:)      
ELSEWHERE
  ZSNOWALB(:) = IR%TSNOW%ALB(:,1)
ENDWHERE
!
ZRGL  (:) = SUM(IP%XPATCH(:,:) * IMT%XRGL  (:,:),DIM=2)
ZRSMIN(:) = SUM(IP%XPATCH(:,:) * IMT%XRSMIN(:,:),DIM=2)
ZGAMMA(:) = SUM(IP%XPATCH(:,:) * IMT%XGAMMA(:,:),DIM=2)
!
ZEXNA(:) = (PPA(:)/XP00)**(XRD/XCPD)
ZEXNS(:) = (PPS(:)/XP00)**(XRD/XCPD)
ZQA  (:) = PQA(:) / PRHOA(:)
ZWIND(:) = SQRT(PU**2+PV**2)
!
!We compute the snow fractions
 CALL ISBA_SNOW_FRAC(IR%TSNOW%SCHEME, ZSNOWSWE, ZSNOWRHO, ZSNOWALB,   &
                     ZVEG, ZLAI, ZZ0, ZPSN, ZPSNV_A, ZPSNG, ZPSNV   )  
!
!We compute total shortwave incoming radiation needed by veg
ZP_GLOBAL_SW(:) = 0.
DO JSWB=1,SIZE(PSW_BANDS)
  ZP_GLOBAL_SW(:)   = ZP_GLOBAL_SW(:) + (PDIR_SW(:,JSWB) + PSCA_SW(:,JSWB))
END DO
!
!We choose the case HPHOTO=='NON' and a humid soil (ZF2=1) to compute ZRS
ZF2(:)=1.0
 CALL VEG(ZP_GLOBAL_SW, PTA, ZQA, PPS, ZRGL, ZLAI, ZRSMIN, ZGAMMA, ZF2, ZRS)
!Calculation of ZDELTA
 CALL WET_LEAVES_FRAC(ZWR, ZVEG, ZWRMAX_CF, ZZ0_WITH_SNOW, ZLAI, ZWRMAX, ZDELTA)
!
!We choose the case LFLOOD=false
ZFFG   (:) = 0.0
ZFFGNOS(:) = 0.0
ZFFV   (:) = 0.0
ZFFVNOS(:) = 0.0
ZFF    (:) = 0.0
!
ZF5    (:) = 1.0
ZLVTT  (:) = XLVTT
!We compute ZCD, ZCH and ZRI
 CALL DRAG(IO%CISBA, IR%TSNOW%SCHEME, IO%CCPSURF,  PTSTEP, ZTS, ZWG, ZWGI, &
           ZEXNS, ZEXNA, PTA, ZWIND, ZQA, PRAIN, PSNOW, PPS, ZRS, ZVEG,    &
           ZZ0, ZZ0EFF, ZZ0H, IP%XWFC(:,1), IP%XWSAT(:,1), ZPSNG, ZPSNV,   &
           PZREF, PUREF, ZP_SLOPE_COS, ZDELTA, ZF5, ZRESA, ZCH, ZCD, ZCDN, &
           ZRI, ZHUG, ZHUGI, ZHV, ZHU, ZCPS, ZQS, ZFFG, ZFFV, ZFF, ZFFGNOS,&
           ZFFVNOS, ZLEG_DELTA, ZLEGI_DELTA, ZWR, PRHOA, ZLVTT            )  
!
!Initialisation of T, Q, Wind and TKE on all canopy levels
DO JLAYER=1,ICP%NLVL
  !
  CALL CLS_TQ(PTA, ZQA, PPA, PPS, PZREF, ZCD, ZCH, ZRI, ZTS, ZHU, ZZ0H, &
              ICP%XZ(:,JLAYER), ZTNM, ZQNM, ZHUNM           ) 
  ! 
  ICP%XT(:,JLAYER)=ZTNM
  ICP%XQ(:,JLAYER)=ZQNM
  !
  CALL CLS_WIND(PU, PV, PUREF, ZCD, ZCDN, ZRI, ICP%XZ(:,JLAYER), &
                ZCLS_WIND_ZON, ZCLS_WIND_MER                 )
  !
  ICP%XU(:,JLAYER) = SQRT( ZCLS_WIND_ZON(:)**2 + ZCLS_WIND_MER(:)**2 )
  ICP%XTKE (:,JLAYER) = XALPSBL * ZCD(:) * ( PU(:)**2 + PV(:)**2 )
  ICP%XP   (:,JLAYER) = PPA(:) + XG * PRHOA(:) * (ICP%XZ(:,ICP%NLVL) - ICP%XZ(:,JLAYER))
  !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('INIT_ISBA_SBL',1,ZHOOK_HANDLE) 
!
END SUBROUTINE INIT_ISBA_SBL
