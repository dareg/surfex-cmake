!     ###############################################################################
SUBROUTINE COUPLING_TSZ0_n (DTCO, UG, U, USS, ICP, AG, CHI, DTI, DGI, GB, ISS, IG,  &
                            I, PKCI, PKD, PK, DTZ, DST, SLT, HPROGRAM, HCOUPLING,   &
                            PTSTEP, KYEAR, KMONTH, KDAY, PTIME, KI, KSV, KSW, PTSUN,&
                            PZENITH, PZENITH2, PAZIM, PZREF, PUREF, PZS, PU, PV,    &
                            PQA, PTA, PRHOA, PSV, PCO2, HSV, PRAIN, PSNOW, PLW,     &
                            PDIR_SW, PSCA_SW, PSW_BANDS, PPS, PPA, PSFTQ, PSFTH,    &
                            PSFTS, PSFCO2, PSFU, PSFV, PTRAD, PDIR_ALB, PSCA_ALB,   &
                            PEMIS, PTSURF, PZ0, PZ0H, PQSURF, PPEW_A_COEF,          &
                            PPEW_B_COEF, PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,     &
                            PPEQ_B_COEF, HTEST      )  
!     ###############################################################################
!
!!****  *COUPLING_TSZ0_n * - Call of fluxes from vegetation scheme ISBA but 
!!        without temporal evolution of the soil/vegetation.
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
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified    09/2012 : J. Escobar , SIZE(PTA) not allowed without-interface , replace by KI
!!      B. Decharme 04/2013 new coupling variables
!!      P. LeMoigne 12/2014 bug in "implicit" coefficients 
!!------------------------------------------------------------------
!
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_SURFEX_n, ONLY : ISBA_DIAG_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_PACK_CH_ISBA, ONLY : PACK_CH_ISBA_t
USE MODD_PACK_DIAG_ISBA, ONLY : PACK_DIAG_ISBA_t
USE MODD_PACK_ISBA, ONLY : PACK_ISBA_t
USE MODD_CANOPY_n, ONLY : CANOPY_t
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_SLT_n, ONLY : SLT_t
!
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS,   ONLY : XP00, XRD, XCPD
!
USE MODI_TSZ0
USE MODI_COUPLING_ISBA_OROGRAPHY_n
! 
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(AGRI_t), INTENT(INOUT) :: AG
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(ISBA_DIAG_t), INTENT(INOUT) :: DGI
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(SSO_t), INTENT(INOUT) :: ISS 
TYPE(GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(PACK_CH_ISBA_t), INTENT(INOUT) :: PKCI
TYPE(PACK_DIAG_ISBA_t), INTENT(INOUT) :: PKD
TYPE(PACK_ISBA_t), INTENT(INOUT) :: PK
TYPE(CANOPY_t), INTENT(INOUT) :: ICP
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SSO_t), INTENT(INOUT) :: USS
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(SLT_t), INTENT(INOUT) :: SLT
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=1),    INTENT(IN)  :: HCOUPLING ! type of coupling
                                              ! 'E' : explicit
                                              ! 'I' : implicit
INTEGER,             INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,             INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,             INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                INTENT(IN)  :: PTIME     ! current time since midnight (UTC, s)
INTEGER,             INTENT(IN)  :: KI        ! number of points
INTEGER,             INTENT(IN)  :: KSV       ! number of scalars
INTEGER,             INTENT(IN)  :: KSW       ! number of short-wave spectral bands
REAL, DIMENSION(KI), INTENT(IN)  :: PTSUN     ! solar time                    (s from midnight)
REAL,                INTENT(IN)  :: PTSTEP    ! atmospheric time-step                 (s)
REAL, DIMENSION(KI), INTENT(IN)  :: PZREF     ! height of T,q forcing                 (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PUREF     ! height of wind forcing                (m)
!
REAL, DIMENSION(KI), INTENT(IN)  :: PTA       ! air temperature forcing               (K)
REAL, DIMENSION(KI), INTENT(IN)  :: PQA       ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(KI), INTENT(IN)  :: PRHOA     ! air density                           (kg/m3)
REAL, DIMENSION(KI,KSV),INTENT(IN) :: PSV     ! scalar variables
!                                             ! chemistry:   first char. in HSV: '#'  (molecule/m3)
!                                             !
 CHARACTER(LEN=6), DIMENSION(KSV),INTENT(IN):: HSV  ! name of all scalar variables
REAL, DIMENSION(KI), INTENT(IN)  :: PU        ! zonal wind                            (m/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PV        ! meridian wind                         (m/s)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PDIR_SW ! direct  solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PSCA_SW ! diffuse solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KSW),INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH   ! zenithal angle at t      (radian from the vertical)
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH2  ! zenithal angle at t+1    (radian from the vertical)
REAL, DIMENSION(KI), INTENT(IN)  :: PAZIM     ! azimuthal angle      (radian from North, clockwise)
REAL, DIMENSION(KI), INTENT(IN)  :: PLW       ! longwave radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI), INTENT(IN)  :: PPS       ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PPA       ! pressure at forcing level             (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PZS       ! atmospheric model orography           (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PCO2      ! CO2 concentration in the air          (kg/m3)
REAL, DIMENSION(KI), INTENT(IN)  :: PSNOW     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PRAIN     ! liquid precipitation                  (kg/m2/s)
!
!
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFU      ! zonal momentum flux                   (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFV      ! meridian momentum flux                (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFCO2    ! flux of CO2                           (m/s*kg_CO2/kg_air)
REAL, DIMENSION(KI,KSV),INTENT(OUT):: PSFTS   ! flux of scalar var.                   (kg/m2/s)
!
REAL, DIMENSION(KI), INTENT(OUT) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PDIR_ALB! direct albedo for each spectral band  (-)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PSCA_ALB! diffuse albedo for each spectral band (-)
REAL, DIMENSION(KI), INTENT(OUT) :: PEMIS     ! emissivity                            (-)
!
REAL, DIMENSION(KI), INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
REAL, DIMENSION(KI), INTENT(OUT) :: PZ0       ! roughness length for momentum         (m)
REAL, DIMENSION(KI), INTENT(OUT) :: PZ0H      ! roughness length for heat             (m)
REAL, DIMENSION(KI), INTENT(OUT) :: PQSURF    ! specific humidity at surface          (kg/kg)
!
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_A_COEF! implicit coefficients
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_B_COEF! needed if HCOUPLING='I'
REAL, DIMENSION(KI), INTENT(IN) :: PPET_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPET_B_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_B_COEF
 CHARACTER(LEN=2),    INTENT(IN) :: HTEST ! must be equal to 'OK'

!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(KI,I%O%NGROUND_LAYER,I%O%NPATCH) :: ZTG   ! soil temperature
REAL, DIMENSION(KI,I%O%NGROUND_LAYER,I%O%NPATCH) :: ZWG   ! soil water content
REAL, DIMENSION(KI,I%O%NGROUND_LAYER,I%O%NPATCH) :: ZWGI  ! soil ice content
REAL, DIMENSION(KI,I%O%NPATCH) :: ZWR   ! interception reservoir
REAL, DIMENSION(KI,I%O%NPATCH) :: ZRESA ! aerodynamical resistance
REAL, DIMENSION(KI,I%R%TSNOW%NLAYER,I%O%NPATCH) :: ZWSNOW! snow reservoir
REAL, DIMENSION(KI,I%R%TSNOW%NLAYER,I%O%NPATCH) :: ZRHOSN! snow density
REAL, DIMENSION(KI,I%R%TSNOW%NLAYER,I%O%NPATCH) :: ZHEASN! snow heat content
REAL, DIMENSION(KI,I%O%NPATCH) :: ZALBSN! snow albedo
REAL, DIMENSION(KI,I%O%NPATCH) :: ZEMISN! snow emissivity
!
REAL, DIMENSION(KI)     :: ZPEW_A_COEF ! implicit coefficients
REAL, DIMENSION(KI)     :: ZPEW_B_COEF ! needed if HCOUPLING='I'
REAL, DIMENSION(KI)     :: ZPET_A_COEF
REAL, DIMENSION(KI)     :: ZPEQ_A_COEF
REAL, DIMENSION(KI)     :: ZPET_B_COEF
REAL, DIMENSION(KI)     :: ZPEQ_B_COEF
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('COUPLING_TSZ0_N',0,ZHOOK_HANDLE)
!
!*      1.     Specified evolution of ISBA prognostic variables
!              ------------------------------------------------
!
 CALL TSZ0(DTZ, &
           PTIME, PTSTEP, I%IP%XWFC, I%R%XTG, I%R%XWG)
!
!
!*      2.     Saves the prognostic variables
!              ------------------------------
!
ZTG  (:,:,:) = I%R%XTG        (:,:,:)
ZWG  (:,:,:) = I%R%XWG        (:,:,:)
ZWGI (:,:,:) = I%R%XWGI       (:,:,:)
ZWR  (:,:)   = I%R%XWR        (:,:)
ZRESA(:,:)   = I%R%XRESA      (:,:)
ZWSNOW(:,:,:)= I%R%TSNOW%WSNOW(:,:,:)
ZRHOSN(:,:,:)= I%R%TSNOW%RHO  (:,:,:)
ZALBSN(:,:)  = I%R%TSNOW%ALB  (:,:)
IF (I%R%TSNOW%SCHEME=='3-L' .OR. I%R%TSNOW%SCHEME=='CRO') THEN
  ZHEASN(:,:,:)= I%R%TSNOW%HEAT (:,:,:)
  ZEMISN(:,:)  = I%R%TSNOW%EMIS (:,:)
END IF
!
!
!*      3.     Call to surface scheme
!              ----------------------
!
 CALL COUPLING_ISBA_OROGRAPHY_n(DTCO, UG, U, USS, ICP, AG, CHI, DTI, DGI, GB, ISS, IG, I, &
                                PKCI, PKD, PK, DST, SLT, HPROGRAM, 'E', 0.001, KYEAR,     &
                                KMONTH, KDAY, PTIME,  KI, KSV, KSW, PTSUN, PZENITH,       &
                                PZENITH2, PAZIM, PZREF, PUREF, PZS, PU, PV, PQA, PTA,     &
                                PRHOA, PSV, PCO2, HSV, PRAIN, PSNOW, PLW, PDIR_SW,        &
                                PSCA_SW, PSW_BANDS, PPS, PPA, PSFTQ, PSFTH, PSFTS, PSFCO2,&
                                PSFU, PSFV, PTRAD, PDIR_ALB, PSCA_ALB, PEMIS, PTSURF, PZ0,&
                                PZ0H, PQSURF, PPEW_A_COEF, PPEW_B_COEF, PPET_A_COEF,      &
                                PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF, 'OK'  )  
!
!
!*      4.     Removes temporal evolution of ISBA variables
!              --------------------------------------------
!
!
I%R%XTG  (:,:,:) = ZTG
I%R%XWG  (:,:,:) = ZWG
I%R%XWGI (:,:,:) = ZWGI
I%R%XWR  (:,:)   = ZWR
I%R%XRESA(:,:)   = ZRESA
I%R%TSNOW%WSNOW(:,:,:) = ZWSNOW
I%R%TSNOW%RHO  (:,:,:) = ZRHOSN
I%R%TSNOW%ALB  (:,:)   = ZALBSN
IF (I%R%TSNOW%SCHEME=='3-L' .OR. I%R%TSNOW%SCHEME=='CRO') THEN
  I%R%TSNOW%HEAT (:,:,:) = ZHEASN
  I%R%TSNOW%EMIS (:,:)   = ZEMISN
END IF
!
IF (LHOOK) CALL DR_HOOK('COUPLING_TSZ0_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE COUPLING_TSZ0_n
