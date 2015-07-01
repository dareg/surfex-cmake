!     ###############################################################################
SUBROUTINE COUPLING_WATFLUX_OROG_n (CHF, CHS, CHW, DTS, DGF, DGMF, DGO, DGS, DGSI, DST, F, O, OR, SG, S, SLT, &
                                     DGW, W, WSB, &
                                    HPROGRAM, HCOUPLING, PTIMEC,                              &
                 PTSTEP, KYEAR, KMONTH, KDAY, PTIME, KI, KSV, KSW, PTSUN, PZENITH, PZENITH2, &
                 PAZIM, PZREF, PUREF, PZS, PU, PV, PQA, PTA, PRHOA, PSV, PCO2, HSV,          &
                 PRAIN, PSNOW, PLW, PDIR_SW, PSCA_SW, PSW_BANDS, PPS, PPA,                   &
                 PSFTQ, PSFTH, PSFTS, PSFCO2, PSFU, PSFV,                                    &
                 PTRAD, PDIR_ALB, PSCA_ALB, PEMIS, PTSURF, PZ0, PZ0H, PQSURF,                &
                 PPEW_A_COEF, PPEW_B_COEF,                                                   &
                 PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                         &
                 HTEST                                                                       )  
!     ###############################################################################
!
!!****  *COUPLING_WATFLUX_OROG_n * - Modifies the input forcing if not
!!           initially at lake level
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
!!      B. Decharme   2008   reset the subgrid topographic effect on the forcing
!!      J. Escobar    09/2012  SIZE(PTA) not allowed without-interface , replace by KI
!!      B. Decharme  04/2013 new coupling variables
!!                           improve forcing vertical shift
!!-------------------------------------------------------------
!
!
!
!
USE MODD_CH_FLAKE_n, ONLY : CH_FLAKE_t
USE MODD_CH_SEAFLUX_n, ONLY : CH_SEAFLUX_t
USE MODD_CH_WATFLUX_n, ONLY : CH_WATFLUX_t
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_t
USE MODD_DIAG_FLAKE_n, ONLY : DIAG_FLAKE_t
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DIAG_MISC_FLAKE_t
USE MODD_DIAG_OCEAN_n, ONLY : DIAG_OCEAN_t
USE MODD_DIAG_SEAFLUX_n, ONLY : DIAG_SEAFLUX_t
USE MODD_DIAG_SEAICE_n, ONLY : DIAG_SEAICE_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_OCEAN_n, ONLY : OCEAN_t
USE MODD_OCEAN_REL_n, ONLY : OCEAN_REL_t
USE MODD_SEAFLUX_GRID_n, ONLY : SEAFLUX_GRID_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SLT_n, ONLY : SLT_t
!
USE MODD_DIAG_WATFLUX_n, ONLY : DIAG_WATFLUX_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
USE MODD_WATFLUX_SBL_n, ONLY : WATFLUX_SBL_t
!
USE MODD_SURF_PAR,         ONLY : XUNDEF
USE MODD_CSTS,             ONLY : XCPD, XRD, XP00
!
USE MODD_SURF_ATM, ONLY : LVERTSHIFT
!
USE MODI_FORCING_VERT_SHIFT
!
USE MODI_COUPLING_SEAWAT_SBL_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
TYPE(CH_FLAKE_t), INTENT(INOUT) :: CHF
TYPE(CH_SEAFLUX_t), INTENT(INOUT) :: CHS
TYPE(CH_WATFLUX_t), INTENT(INOUT) :: CHW
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: DTS
TYPE(DIAG_FLAKE_t), INTENT(INOUT) :: DGF
TYPE(DIAG_MISC_FLAKE_t), INTENT(INOUT) :: DGMF
TYPE(DIAG_OCEAN_t), INTENT(INOUT) :: DGO
TYPE(DIAG_SEAFLUX_t), INTENT(INOUT) :: DGS
TYPE(DIAG_SEAICE_t), INTENT(INOUT) :: DGSI
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(OCEAN_t), INTENT(INOUT) :: O
TYPE(OCEAN_REL_t), INTENT(INOUT) :: OR
TYPE(SEAFLUX_GRID_t), INTENT(INOUT) :: SG
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SLT_t), INTENT(INOUT) :: SLT
!
TYPE(DIAG_WATFLUX_t), INTENT(INOUT) :: DGW
TYPE(WATFLUX_t), INTENT(INOUT) :: W
TYPE(WATFLUX_SBL_t), INTENT(INOUT) :: WSB
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=1),    INTENT(IN)  :: HCOUPLING ! type of coupling
                                              ! 'E' : explicit
                                              ! 'I' : implicit
REAL,                INTENT(IN)  :: PTIMEC    ! cumulated time since beginning of simulation
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
REAL, DIMENSION(KI)  ::  ZPEQ_B_COEF   ! 1st explicit coefficient
REAL, DIMENSION(KI)  ::  ZPET_B_COEF   ! 2nd explicit coefficient
!
REAL, DIMENSION(KI)  :: ZTA    ! Temperature at forcing height above surface orography
REAL, DIMENSION(KI)  :: ZPA    ! Pressure    at forcing height above surface orography
REAL, DIMENSION(KI)  :: ZPS    ! Pressure    at surface orography
REAL, DIMENSION(KI)  :: ZQA    ! Humidity    at forcing height above surface orography
REAL, DIMENSION(KI)  :: ZRHOA  ! Density     at forcing height above surface orography
REAL, DIMENSION(KI)  :: ZLW    ! LW rad      at forcing height above surface orography
REAL, DIMENSION(KI)  :: ZRAIN  ! Rainfall    at forcing height above surface orography
REAL, DIMENSION(KI)  :: ZSNOW  ! Snowfall    at forcing height above surface orography
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
! Preliminaries:
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('COUPLING_WATFLUX_OROG_N',0,ZHOOK_HANDLE)
!
ZPEQ_B_COEF(:) = PPEQ_B_COEF(:)
ZPET_B_COEF(:) = PPET_B_COEF(:)
!
IF(LVERTSHIFT)THEN
!
  ZTA  (:) = XUNDEF
  ZQA  (:) = XUNDEF
  ZPS  (:) = XUNDEF
  ZPA  (:) = XUNDEF
  ZRHOA(:) = XUNDEF
  ZLW  (:) = XUNDEF
  ZRAIN(:) = XUNDEF
  ZSNOW(:) = XUNDEF
!       
   CALL FORCING_VERT_SHIFT(PZS,W%XZS,PTA,PQA,PPA,PRHOA,PLW,PRAIN,PSNOW,&
                           ZTA,ZQA,ZPA,ZRHOA,ZLW,ZRAIN,ZSNOW         )
!
   ZPS(:) = ZPA(:) + (PPS(:) - PPA(:))
!
  IF (HCOUPLING=='I') THEN
    ZPEQ_B_COEF = PPEQ_B_COEF + ZQA - PQA
    ZPET_B_COEF = PPET_B_COEF + ZTA/(ZPA/XP00)**(XRD/XCPD) - PTA/(PPA/XP00)**(XRD/XCPD)
  ENDIF
!
ELSE
!
  ZTA  (:) = PTA  (:)
  ZQA  (:) = PQA  (:)
  ZPS  (:) = PPS  (:)
  ZPA  (:) = PPA  (:)
  ZRHOA(:) = PRHOA(:)
  ZLW  (:) = PLW  (:)
  ZRAIN(:) = PRAIN(:)
  ZSNOW(:) = PSNOW(:)
!
ENDIF
!
 CALL COUPLING_SEAWAT_SBL_n(CHF, CHS, CHW, DTS, DGF, DGMF, DGO, DGS, DGSI, DGW, DST, F, O, OR, SG, S, SLT, W, &
                            HPROGRAM, HCOUPLING, PTIMEC, 'W',                    &
               PTSTEP, KYEAR, KMONTH, KDAY, PTIME,                               &
               KI, KSV, KSW,                                                     &
               PTSUN, PZENITH, PZENITH2, PAZIM,                                  &
               PZREF, PUREF, W%XZS, PU, PV, ZQA, ZTA, ZRHOA, PSV, PCO2, HSV,       &
               ZRAIN, ZSNOW, ZLW, PDIR_SW, PSCA_SW, PSW_BANDS, ZPS, ZPA,         &
               PSFTQ, PSFTH, PSFTS, PSFCO2, PSFU, PSFV, W%LSBL, W%XTS, W%XZ0,          &
               WSB%XZ, WSB%XU, WSB%NLVL, WSB%XTKE, WSB%XT, WSB%XQ, WSB%XLMO, WSB%XZF, WSB%XDZ, WSB%XDZF, WSB%XP,             &
               DGW%N2M, DGW%XT2M, DGW%XQ2M, DGW%XHU2M, DGW%XZON10M, DGW%XMER10M, DGW%XWIND10M, DGW%XWIND10M_MAX, &
               DGW%XT2M_MIN, DGW%XT2M_MAX, DGW%XHU2M_MIN, DGW%XHU2M_MAX,                         &
               PTRAD, PDIR_ALB, PSCA_ALB, PEMIS, PTSURF, PZ0, PZ0H, PQSURF,      &
               PPEW_A_COEF, PPEW_B_COEF,                                         &
               PPET_A_COEF, PPEQ_A_COEF, ZPET_B_COEF, ZPEQ_B_COEF,               &
               'OK'                                                              )
!
IF (LHOOK) CALL DR_HOOK('COUPLING_WATFLUX_OROG_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE COUPLING_WATFLUX_OROG_n
