!     #################################################################################
SUBROUTINE ASSIM_SURF_ATM_n(HPROGRAM, KI,                                               &
                            PCON_RAIN, PSTRAT_RAIN, PCON_SNOW, PSTRAT_SNOW,             &
                            PCLOUDS,   PLSM,        PEVAPTR,   PEVAP,                   &
                            PSWEC,     PTSC,                                            &
                            PTS,       PT2M,        PHU2M,     PSWE,                    &
                            PSST,      PSIC,                                            &
                            HTEST , ODINLINE, OD_MASKEXT, PLON, PLAT)
!     #################################################################################
!
!
!!****  *ASSIM_SURF_ATM_n * - Driver to call the schemes for the 
!!       four surface types (SEA, WATER, NATURE, TOWN)
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
!!-------------------------------------------------------------
!
USE MODD_SURF_CONF,      ONLY : CPROGNAME
USE MODD_SURF_ATM_n,     ONLY : NSIZE_SEA, NSIZE_WATER, NSIZE_TOWN, NSIZE_NATURE, &
                                NR_SEA,    NR_WATER,    NR_TOWN,    NR_NATURE
USE MODD_ISBA_n,      ONLY : XTG, TSNOW, XWG, XWGI
USE MODD_SEAFLUX_n,   ONLY : XSST
USE MODD_WATFLUX_n,   ONLY : XTS
USE MODD_TEB_n,       ONLY : XT_ROAD
!
USE MODD_ASSIM,          ONLY : XAT2M_ISBA, XAHU2M_ISBA, XAZON10M_ISBA, XAMER10M_ISBA, XAT2M_TEB, LAROME
!
#ifdef ARO 
USE MODD_IO_SURF_ARO,ONLY : LWRITE, LCOUNTW, LFMWRIT, XGPGW, YSURFEX_CACHE_OUT,      &
                            SURFEX_FIELD_BUF_PREALLOC, SURFEX_FIELD_BUF_SET_RECORD,  &
                            NCOUNTW, NCOUNTW_TOT 
#endif
!
USE MODI_ABOR1_SFX
USE MODI_ASSIM_SEA_n
USE MODI_ASSIM_INLAND_WATER_n
USE MODI_ASSIM_NATURE_n
USE MODI_ASSIM_TOWN_n
!
USE YOMHOOK,             ONLY : LHOOK,   DR_HOOK
USE PARKIND1,            ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),    INTENT(IN) :: HPROGRAM     ! program calling surf. schemes
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
REAL, DIMENSION(KI), INTENT(IN) :: PSST
REAL, DIMENSION(KI), INTENT(IN) :: PSIC
CHARACTER(LEN=2),   INTENT(IN)  :: HTEST        ! must be equal to 'OK'
LOGICAL, INTENT (IN) :: ODINLINE
LOGICAL,  DIMENSION (KI), INTENT(IN) ::  OD_MASKEXT
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLON
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLAT
!
!*      0.2    declarations of local variables
!
INTEGER :: JTILE                        ! loop on type of surface
LOGICAL :: GNATURE, GTOWN, GWATER, GSEA ! .T. if the corresponding surface is represented
LOGICAL :: GLKEEPEXTZONE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ASSIM_SURF_ATM_N',0,ZHOOK_HANDLE)
!
CPROGNAME = HPROGRAM
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_SURF_ATMN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
GLKEEPEXTZONE = .FALSE.
!
!-------------------------------------------------------------------------------------
! Preliminaries: Tile related operations
!-------------------------------------------------------------------------------------

! FLAGS for the various surfaces:
!
GSEA      = NSIZE_SEA    >0
GWATER    = NSIZE_WATER  >0
GTOWN     = NSIZE_TOWN   >0
GNATURE   = NSIZE_NATURE >0
!
! Tile counter:
!
JTILE     = 0 
!
!--------------------------------------------------------------------------------------
! Call interfaces for sea, water, nature and town here...
!--------------------------------------------------------------------------------------
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! SEA Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE = JTILE + 1
!
IF(GSEA)THEN
!
  CALL ASSIM_TREAT_SURF(JTILE,NSIZE_SEA,NR_SEA)
!
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! INLAND WATER Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE = JTILE + 1
!
IF(GWATER)THEN
!
  CALL ASSIM_TREAT_SURF(JTILE,NSIZE_WATER,NR_WATER)
!
ENDIF 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! NATURAL SURFACE Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE = JTILE + 1
!
IF(GNATURE)THEN
!
  CALL ASSIM_TREAT_SURF(JTILE,NSIZE_NATURE,NR_NATURE)

  IF ( ALLOCATED(XAT2M_ISBA))    DEALLOCATE(XAT2M_ISBA)
  IF ( ALLOCATED(XAHU2M_ISBA))   DEALLOCATE(XAHU2M_ISBA)
  IF ( ALLOCATED(XAZON10M_ISBA)) DEALLOCATE(XAZON10M_ISBA)
  IF ( ALLOCATED(XAMER10M_ISBA)) DEALLOCATE(XAMER10M_ISBA)
!
ENDIF 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! URBAN Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE = JTILE + 1
!
IF(GTOWN)THEN
!
  CALL ASSIM_TREAT_SURF(JTILE,NSIZE_TOWN,NR_TOWN)

  IF ( ALLOCATED(XAT2M_TEB))    DEALLOCATE(XAT2M_TEB)
!
ENDIF
!
IF (ODINLINE) THEN
#ifdef ARO
! Count 2D fields in MSE
  NCOUNTW_TOT = 0
  LWRITE      = .FALSE.
  LCOUNTW     = .TRUE.
  NCOUNTW     = 0
! NINDX1, NINDX2, NKPROMA already set 
  IF (.NOT. LFMWRIT) CALL SURFEX_FIELD_BUF_SET_RECORD(YSURFEX_CACHE_OUT, .FALSE.)
  !
  CALL WRITE
  !  
  IF (LFMWRIT) THEN
    ALLOCATE (XGPGW(IGPCOMP, NCOUNTW_TOT))
    XGPGW = XUNDEF
  ELSE
    CALL SURFEX_FIELD_BUF_PREALLOC  (YSURFEX_CACHE_OUT)
    CALL SURFEX_FIELD_BUF_SET_RECORD(YSURFEX_CACHE_OUT, .TRUE.)
  ENDIF
  !
  LWRITE      = .TRUE.
  LCOUNTW     = .FALSE.
  NCOUNTW     = 0

  IF (LFMWRIT) DEALLOCATE (XGPGW)
#endif
ENDIF

IF (LHOOK) CALL DR_HOOK('ASSIM_SURF_ATM_N',1,ZHOOK_HANDLE)
!
!=======================================================================================
CONTAINS

SUBROUTINE WRITE 

USE MODI_WRITE_SURF

CHARACTER(LEN=10)  :: YVAR
CHARACTER(LEN=100) :: YPREFIX ! Prefix of the prognostic variable  (in LFI file)
INTEGER :: IRESP

IF (LHOOK) CALL DR_HOOK ('ASSIM_SURF_ATM_n:WRITE', 0, ZHOOK_HANDLE)

CALL DD ('WG1', XWG(:,1,1))

YVAR='WG1'
YPREFIX='X_Y_WG1 (m3/m3)                                   '
CALL WRITE_SURF(HPROGRAM,YVAR,XWG(:,1,1),IRESP,HCOMMENT=YPREFIX)

CALL DD ('WG2', XWG(:,2,1))

YVAR='WG2'
YPREFIX='X_Y_WG2 (m3/m3)                                   '
CALL WRITE_SURF(HPROGRAM,YVAR,XWG(:,2,1),IRESP,HCOMMENT=YPREFIX)

CALL DD ('WGI2', XWGI(:,2,1))

YVAR='WGI2'
YPREFIX='X_Y_WGI2 (m3/m3)                                  '
CALL WRITE_SURF(HPROGRAM,YVAR,XWGI(:,2,1),IRESP,HCOMMENT=YPREFIX)

CALL DD ('TG1', XTG(:,1,1))

YVAR='TG1'
YPREFIX='X_Y_TG1 (K)                                       '
CALL WRITE_SURF(HPROGRAM,YVAR,XTG(:,1,1),IRESP,HCOMMENT=YPREFIX)

CALL DD ('TG2', XTG(:,2,1))

YVAR='TG2'
YPREFIX='X_Y_TG2 (K)                                       '
CALL WRITE_SURF(HPROGRAM,YVAR,XTG(:,2,1),IRESP,HCOMMENT=YPREFIX)

CALL DD ('SST', XSST(:))

YVAR='SST'
YPREFIX='X_Y_SST (K)                                       '
CALL WRITE_SURF(HPROGRAM,YVAR,XSST(:),IRESP,HCOMMENT=YPREFIX)

CALL DD ('TS_WATER', XTS(:))

YVAR='TS_WATER'
YPREFIX='X_Y_TS_WATER (K)                                  '
CALL WRITE_SURF(HPROGRAM,YVAR,XTS(:),IRESP,HCOMMENT=YPREFIX)

IF (NSIZE_TOWN > 0 .AND. LAROME) THEN
  CALL DD ('T_ROAD3', XT_ROAD(:,3))

  YVAR='TROAD3'
  YPREFIX='X_Y_T_ROAD3 (K)                                   '
  CALL WRITE_SURF(HPROGRAM,YVAR,XT_ROAD(:,3),IRESP,HCOMMENT=YPREFIX)
ENDIF

CALL DD ('WSNOW_VEG1', TSNOW%WSNOW(:,1,1))

YVAR='WSN_VEG1'
YPREFIX='X_Y_WSNOW_VEG1 (kg/m2)                            '
CALL WRITE_SURF(HPROGRAM,YVAR,TSNOW%WSNOW(:,1,1),IRESP,HCOMMENT=YPREFIX)

IF (LHOOK) CALL DR_HOOK ('ASSIM_SURF_ATM_n:WRITE', 1, ZHOOK_HANDLE)

END SUBROUTINE WRITE

SUBROUTINE DD (CDN, PX)
CHARACTER(LEN=*), INTENT (IN) :: CDN
REAL, INTENT (IN) :: PX (:)

REAL :: ZX (SIZE (PX))
INTEGER :: I, N
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('ASSIM_SURF_ATM_n:DD', 0, ZHOOK_HANDLE)

IF (ODINLINE) THEN
#ifdef ARO
  IF (.NOT.LWRITE.AND.LHOOK) CALL DR_HOOK ('ASSIM_SURF_ATM_n:DD', 1, ZHOOK_HANDLE)
  IF (.NOT.LWRITE) RETURN
#endif
  N = COUNT (.NOT. OD_MASKEXT)
  ZX (1:N) = PACK (PX, .NOT. OD_MASKEXT)
ELSE
  ZX = PX
  N = SIZE (PX)
ENDIF

WRITE (0, *) TRIM(CDN)//" = " 
WRITE (0, *) N, MINVAL(ZX(1:N)), MAXVAL(ZX(1:N))
!WRITE (0, '(10(E14.6,", "))') ZX (1:N)

IF (LHOOK) CALL DR_HOOK ('ASSIM_SURF_ATM_n:DD', 1, ZHOOK_HANDLE)

END SUBROUTINE DD
!
!=======================================================================================
SUBROUTINE ASSIM_TREAT_SURF(KTILE,KSIZE,KMASK)
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                   :: KTILE
INTEGER, INTENT(IN)                   :: KSIZE
INTEGER, INTENT(IN), DIMENSION(KSIZE) :: KMASK
!
REAL,DIMENSION(KSIZE)                 :: ZP_PCON_RAIN
REAL,DIMENSION(KSIZE)                 :: ZP_PSTRAT_RAIN
REAL,DIMENSION(KSIZE)                 :: ZP_PCON_SNOW
REAL,DIMENSION(KSIZE)                 :: ZP_PSTRAT_SNOW
REAL,DIMENSION(KSIZE)                 :: ZP_PCLOUDS
REAL,DIMENSION(KSIZE)                 :: ZP_PLSM
REAL,DIMENSION(KSIZE)                 :: ZP_PEVAPTR
REAL,DIMENSION(KSIZE)                 :: ZP_PEVAP
REAL,DIMENSION(KSIZE)                 :: ZP_PSWEC
REAL,DIMENSION(KSIZE)                 :: ZP_PTSC
REAL,DIMENSION(KSIZE)                 :: ZP_PTS
REAL,DIMENSION(KSIZE)                 :: ZP_PT2M
REAL,DIMENSION(KSIZE)                 :: ZP_PHU2M
REAL,DIMENSION(KSIZE)                 :: ZP_PSWE
REAL,DIMENSION(KSIZE)                 :: ZP_PSST
REAL,DIMENSION(KSIZE)                 :: ZP_PSIC
REAL,DIMENSION(KSIZE)                 :: ZP_LON
REAL,DIMENSION(KSIZE)                 :: ZP_LAT
LOGICAL,DIMENSION(KSIZE)              :: GD_MASKEXT
INTEGER                               :: JJ,JI
!
DO JJ=1,KSIZE
  JI=KMASK(JJ)
  ZP_PLSM(JJ)        = PLSM(JI)  
  ZP_PCON_RAIN(JJ)   = PCON_RAIN(JI)
  ZP_PSTRAT_RAIN(JJ) = PSTRAT_RAIN(JI)
  ZP_PCON_SNOW(JJ)   = PCON_SNOW(JI)
  ZP_PSTRAT_SNOW(JJ) = PSTRAT_SNOW(JI)
  ZP_PCLOUDS(JJ)     = PCLOUDS(JI)
  ZP_PEVAPTR(JJ)     = PEVAPTR(JI)
  ZP_PEVAP(JJ)       = PEVAP(JI)
  ZP_PSWE(JJ)        = PSWE(JI)  
  ZP_PSWEC(JJ)       = PSWEC(JI)
  ZP_PTSC(JJ)        = PTSC(JI)
  ZP_PTS(JJ)         = PTS(JI) 
  ZP_PT2M(JJ)        = PT2M(JI)
  ZP_PHU2M(JJ)       = PHU2M(JI)
  ZP_PSST(JJ)        = PSST(JI)
  ZP_PSIC(JJ)        = PSIC(JI)
  ZP_LON(JJ)         = PLON(JI)
  ZP_LAT(JJ)         = PLAT(JI)
  GD_MASKEXT(JJ)     = OD_MASKEXT(JI)
ENDDO

IF (KTILE==1) THEN
 
  WRITE(*,*) '*********************************************'
  WRITE(*,*) '*      ASSIMILATIONS FOR SEA POINTS         *'
  WRITE(*,*) '*********************************************'
 
  CALL  ASSIM_SEA_n(HPROGRAM,KSIZE,ZP_PTS,ZP_PSST,ZP_PSIC,ZP_PLSM,HTEST,&
              ODINLINE,GLKEEPEXTZONE,GD_MASKEXT)

ELSEIF (KTILE==2) THEN
  
  WRITE(*,*) '*********************************************'  
  WRITE(*,*) '*      ASSIMILATIONS FOR WATER POINTS       *'
  WRITE(*,*) '*********************************************'
  CALL ASSIM_INLAND_WATER_n(HPROGRAM,KSIZE,ZP_PTS,ZP_PLSM,HTEST,&
                            ODINLINE,GLKEEPEXTZONE,GD_MASKEXT)

ELSEIF (KTILE==3) THEN
  
  WRITE(*,*) '*********************************************'  
  WRITE(*,*) '*      ASSIMILATIONS FOR NATURE POINTS      *'
  WRITE(*,*) '*********************************************'
  CALL ASSIM_NATURE_n(HPROGRAM,KSIZE,                                             &
                      ZP_PCON_RAIN, ZP_PSTRAT_RAIN, ZP_PCON_SNOW, ZP_PSTRAT_SNOW, &
                      ZP_PCLOUDS,   ZP_PLSM,        ZP_PEVAPTR,   ZP_PEVAP,       & 
                      ZP_PSWEC,     ZP_PTSC,                                      &
                      ZP_PTS,       ZP_PT2M,        ZP_PHU2M,     ZP_PSWE,        & 
                      HTEST, ODINLINE, GD_MASKEXT, ZP_LON, ZP_LAT )
  
ELSEIF (KTILE==4) THEN
  
  WRITE(*,*) '*********************************************'  
  WRITE(*,*) '*      ASSIMILATIONS FOR URBAN POINTS       *'
  WRITE(*,*) '*********************************************'
  CALL ASSIM_TOWN_n(HPROGRAM,KSIZE,ZP_PT2M,HTEST)
  
ENDIF

END SUBROUTINE ASSIM_TREAT_SURF
!=======================================================================================
END SUBROUTINE ASSIM_SURF_ATM_n
!=======================================================================================

