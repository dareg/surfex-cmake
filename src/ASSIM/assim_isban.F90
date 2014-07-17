!     ###############################################################################
SUBROUTINE ASSIM_ISBA_n(HPROGRAM,KI,                                   &
                        PCON_RAIN, PSTRAT_RAIN, PCON_SNOW, PSTRAT_SNOW,&
                        PCLOUDS,   PLSM,        PEVAPTR,   PEVAP,      &
                        PSWEC,     PTSC,                               &
                        PTS,       PT2M,        PHU2M,     PSWE,       &
                        HTEST, OD_MASKEXT, PLON_IN, PLAT_IN )

!     ###############################################################################
!
!!****  *ASSIM_ISBA_n * - Chooses the surface assimilation schemes for ISBA
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     T. Aspelien
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2012
!!      Trygve Aspelien, Separating IO  06/2013
!!--------------------------------------------------------------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_ASSIM,          ONLY : CASSIM_ISBA,LAESNM,LEXTRAP_NATURE,NPRINTLEV
!
USE MODD_SURF_ATM_n,     ONLY : NR_NATURE
USE MODD_ISBA_GRID_n,    ONLY : XLAT, XLON
USE MODD_ISBA_n,         ONLY : TSNOW,XWG,XWGI,XTG,XZS
!
USE YOMHOOK,             ONLY : LHOOK,   DR_HOOK
USE PARKIND1,            ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_OI_HOR_EXTRAPOL_SURF
USE MODI_ASSIM_ISBA_UPDATE_SNOW
USE MODI_ASSIM_NATURE_ISBA_EKF
USE MODI_ASSIM_NATURE_ISBA_OI
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),    INTENT(IN) :: HPROGRAM  ! program calling surf. schemes
INTEGER,             INTENT(IN) :: KI
REAL, DIMENSION(KI), INTENT(IN) :: PCON_RAIN
REAL, DIMENSION(KI), INTENT(IN) :: PSTRAT_RAIN
REAL, DIMENSION(KI), INTENT(IN) :: PCON_SNOW
REAL, DIMENSION(KI), INTENT(IN) :: PSTRAT_SNOW
REAL, DIMENSION(KI), INTENT(IN) :: PCLOUDS
REAL, DIMENSION(KI), INTENT(IN) :: PLSM
REAL, DIMENSION(KI), INTENT(IN) :: PEVAPTR
REAL, DIMENSION(KI), INTENT(IN) :: PEVAP
REAL, DIMENSION(KI), INTENT(IN) :: PSWEC
REAL, DIMENSION(KI), INTENT(IN) :: PTSC
REAL, DIMENSION(KI), INTENT(IN) :: PTS
REAL, DIMENSION(KI), INTENT(IN) :: PT2M
REAL, DIMENSION(KI), INTENT(IN) :: PHU2M
REAL, DIMENSION(KI), INTENT(IN) :: PSWE
CHARACTER(LEN=2),    INTENT(IN) :: HTEST ! must be equal to 'OK'
LOGICAL,  DIMENSION (KI), INTENT(IN) ::  OD_MASKEXT
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLON_IN
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLAT_IN
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
LOGICAL, DIMENSION(:), ALLOCATABLE :: GINTERP_NATURE
LOGICAL, DIMENSION(:), ALLOCATABLE :: GINTERP_SN
REAL,    DIMENSION(:), ALLOCATABLE :: ZTS_EP,ZTS_EP0
REAL,    DIMENSION(:), ALLOCATABLE :: ZTP_EP,ZTP_EP0
REAL,    DIMENSION(:), ALLOCATABLE :: ZT3_EP,ZT3_EP0
REAL,    DIMENSION(:), ALLOCATABLE :: ZWS_EP,ZWS_EP0
REAL,    DIMENSION(:), ALLOCATABLE :: ZWP_EP,ZWP_EP0
REAL,    DIMENSION(:), ALLOCATABLE :: ZTL_EP,ZTL_EP0
REAL,    DIMENSION(:), ALLOCATABLE :: ZSWE_EP,ZSWE_EP0
REAL,    DIMENSION(:), ALLOCATABLE :: ZSNR_EP,ZSNR_EP0
REAL,    DIMENSION(:), ALLOCATABLE :: ZSNA_EP,ZSNA_EP0
REAL,    DIMENSION(KI) :: ZSWE
REAL,    DIMENSION(KI) :: ZSWE_ORIG
INTEGER :: I,IL,IP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_N',0,ZHOOK_HANDLE)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_ISBA_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
! Set snow layers and patches
IP = 1
IL = 1
!
ZSWE = PSWE
!
! Soil assimilation
IF ( CASSIM_ISBA == 'EKF  ' ) THEN
  !
  ! Snow analysis/update
  IF (LAESNM) THEN
    WRITE(*,*) 'UPDATE SNOW FROM ANALYSED CANARI VALUES'
    CALL ASSIM_ISBA_UPDATE_SNOW(HPROGRAM,KI,ZSWE,ZSWE_ORIG,.TRUE.,.TRUE.,HTEST)
  ELSE
    WRITE(*,*) 'SNOW IS NOT UPDATED FROM ANALYSED CANARI VALUES'
  ENDIF
  !
  ! Run EKF for soil
  CALL ASSIM_NATURE_ISBA_EKF(HPROGRAM, KI, PT2M, PHU2M, HTEST)
  !
ELSEIF ( CASSIM_ISBA == 'OI   ' ) THEN
  !
  ! Snow analysis/update. Store the original field in the surfex file
  IF (LAESNM) THEN
    WRITE(*,*) 'UPDATE SNOW FROM ANALYSED CANARI VALUES'
    CALL ASSIM_ISBA_UPDATE_SNOW(HPROGRAM,KI,ZSWE,ZSWE_ORIG,.TRUE.,.FALSE.,HTEST)
  ELSE
    WRITE(*,*) 'SNOW IS NOT UPDATED FROM ANALYSED CANARI VALUES'
  ENDIF
  !
  ! Run OI for soil
  CALL ASSIM_NATURE_ISBA_OI(HPROGRAM, KI,                                  &
                            PCON_RAIN, PSTRAT_RAIN, PCON_SNOW, PSTRAT_SNOW,&
                            PCLOUDS,   PLSM,        PEVAPTR,   PEVAP,      &
                            PSWEC,     PTSC,                               &
                            PTS,       PT2M,        PHU2M,     ZSWE,       &
                            HTEST, OD_MASKEXT, PLON_IN, PLAT_IN )
  !
  ! Snow analysis/update (changed in oi_cacsts). Get the full increment
  IF (LAESNM) THEN
    WRITE(*,*) 'UPDATE SNOW FROM ANALYSED OI_CACSTS VALUES'
    CALL ASSIM_ISBA_UPDATE_SNOW(HPROGRAM,KI,ZSWE,ZSWE_ORIG,.FALSE.,.TRUE.,HTEST)
  ELSE
    WRITE(*,*) 'SNOW IS NOT UPDATED FROM ANALYSED OI_CACSTS VALUES'
  ENDIF
  !
ELSE
  CALL ABOR1_SFX(CASSIM_ISBA//' is not a defined scheme for ASSIM_ISBA_N')
ENDIF

! Extrapolation if requested
IF ( LEXTRAP_NATURE ) THEN
  !
  ALLOCATE(ZWS_EP(KI),ZWP_EP(KI),ZTS_EP(KI),ZTP_EP(KI),ZT3_EP(KI),&
           ZTL_EP(KI),ZSWE_EP(KI),ZSNR_EP(KI),ZSNA_EP(KI))
  !  
  ZWS_EP  = XWG(:,1,IP)
  ZWP_EP  = XWG(:,2,IP)
  ZTS_EP  = XTG(:,1,IP)
  ZTP_EP  = XTG(:,2,IP)
  ZT3_EP  = XTG(:,3,IP)
  ZTL_EP  = XWGI(:,2,IP)
  ZSWE_EP = TSNOW%WSNOW(:,IL,IP)
  ZSNR_EP = TSNOW%RHO  (:,IL,IP)
  ZSNA_EP = TSNOW%ALB  (:,   IP)
  !
  ALLOCATE(GINTERP_NATURE(KI),GINTERP_SN(KI))
  !
  ! Search for the nearest grid point values for land surface fields
  ! at locations where the CANARI land fraction is less than 50%
  ! and therefore useless values MIGTH be given
  GINTERP_NATURE = .FALSE.
  GINTERP_SN     = .FALSE.
  !
  ! Snow albedo and density are also extrapolated in points 
  ! which get initial snow in the snow analysis
  WHERE ( ZSWE_EP(:) < 1.0E-10 .AND. PSWE(:)>= 1.0E-10 )
    GINTERP_SN(:) = .TRUE.
    ZSNA_EP(:)    = XUNDEF
    ZSNR_EP(:)    = XUNDEF
  END WHERE
  ZSWE_EP(:) = PSWE(:)
  !
  WHERE ( PLSM(:) < 0.5 )
    GINTERP_NATURE(:) = .TRUE.
    GINTERP_SN(:) = .TRUE.
    ZTS_EP(:)     = XUNDEF
    ZTP_EP(:)     = XUNDEF
    ZT3_EP(:)     = XUNDEF
    ZWS_EP(:)     = XUNDEF
    ZWP_EP(:)     = XUNDEF
    ZTL_EP(:)     = XUNDEF
    ZSWE_EP(:)    = XUNDEF
    ZSNA_EP(:)    = XUNDEF
    ZSNR_EP(:)    = XUNDEF
  END WHERE
  !
  ALLOCATE(ZWS_EP0(KI),ZWP_EP0(KI),ZTS_EP0(KI),ZTP_EP0(KI),ZT3_EP0(KI),&
           ZTL_EP0(KI),ZSWE_EP0(KI),ZSNR_EP0(KI),ZSNA_EP0(KI))
  !
  ZWS_EP0(:) = ZWS_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,XLAT,XLON,ZWS_EP0,XLAT,XLON,ZWS_EP,GINTERP_NATURE)
  ZWP_EP0(:) = ZWP_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,XLAT,XLON,ZWP_EP0,XLAT,XLON,ZWP_EP,GINTERP_NATURE)
  ZTS_EP0(:) = ZTS_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,XLAT,XLON,ZTS_EP0,XLAT,XLON,ZTS_EP,GINTERP_NATURE,XZS)
  ZTP_EP0(:) = ZTP_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,XLAT,XLON,ZTP_EP0,XLAT,XLON,ZTP_EP,GINTERP_NATURE,XZS)
  ZT3_EP0(:) = ZT3_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,XLAT,XLON,ZT3_EP0,XLAT,XLON,ZT3_EP,GINTERP_NATURE,XZS)
  ZTL_EP0(:) = ZTL_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,XLAT,XLON,ZTL_EP0,XLAT,XLON,ZTL_EP,GINTERP_NATURE)
  ZSWE_EP0(:) = ZSWE_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,XLAT,XLON,ZSWE_EP0,XLAT,XLON,ZSWE_EP,GINTERP_SN)
  ZSNR_EP0(:) = ZSNR_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,XLAT,XLON,ZSNR_EP0,XLAT,XLON,ZSNR_EP,GINTERP_SN)
  ZSNA_EP0(:) = ZSNA_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,XLAT,XLON,ZSNA_EP0,XLAT,XLON,ZSNA_EP,GINTERP_SN)  
  !
  DEALLOCATE(ZWS_EP0,ZWP_EP0,ZTS_EP0,ZTP_EP0,ZT3_EP0,ZTL_EP0,ZSWE_EP0,ZSNR_EP0,ZSNA_EP0)
  !
  ! PRINT values produced by OI_HO_EXTRAPOL_SURF for TS
  IF ( NPRINTLEV > 2 ) THEN
    DO I=1,KI
     IF (GINTERP_NATURE(I)) THEN
       PRINT *,'Surface temperature set to ',ZTS_EP(I),'from nearest neighbour at I=',NR_NATURE(I)
     ENDIF
    ENDDO
  ENDIF
  !
  DEALLOCATE(GINTERP_NATURE,GINTERP_SN)
  !
  ! Set extrpolated fields to global
  XWG (:,1,IP) = ZWS_EP(:)
  XWG (:,2,IP) = ZWP_EP(:)
  XTG (:,1,IP) = ZTS_EP(:)
  XTG (:,2,IP) = ZTP_EP(:)
  XTG (:,3,IP) = ZT3_EP(:)
  XWGI(:,2,IP) = ZTL_EP(:)
  TSNOW%WSNOW(:,IL,IP) = ZSWE_EP(:)
  TSNOW%RHO  (:,IL,IP) = ZSNR_EP(:)
  TSNOW%ALB  (:,   IP) = ZSNA_EP(:)
  !
  DEALLOCATE(ZWS_EP,ZWP_EP,ZTS_EP,ZTP_EP,ZT3_EP,ZTL_EP,ZSWE_EP,ZSNR_EP,ZSNA_EP)
  !
ENDIF

! Snow analysis/update security
IF (LAESNM) THEN

  ! removes very small values due to computation precision
  WHERE( TSNOW%WSNOW(:,IL,IP) < 1.0E-10 ) TSNOW%WSNOW(:,IL,IP) = 0.0

  ! No SNOW
  WHERE ( TSNOW%WSNOW(:,IL,IP) == 0.0 )
    TSNOW%RHO(:,IL,IP) = XUNDEF
    TSNOW%ALB(:,IP)    = XUNDEF
  END WHERE
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_ISBA_n
