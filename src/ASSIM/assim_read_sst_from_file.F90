!     ###############################################################################
SUBROUTINE ASSIM_READ_SST_FROM_FILE(KI,PITM,PSST,PSIC,HTEST)

!     ###############################################################################
!
!!****  *ASSIM_READ_SST_FROM_FILE * - Reads SST from file
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
!!--------------------------------------------------------------------
!
USE MODD_SURF_ATM_n, ONLY : NR_SEA, NSIZE_SEA
USE MODD_SEAFLUX_n, ONLY : XSST
!
USE MODD_ASSIM,         ONLY : LECSST, LREAD_SST_FROM_FILE
USE MODD_SURF_PAR,      ONLY : XUNDEF
!
#ifdef FA
USE MODD_IO_SURF_FA,    ONLY : CFILEIN_FA, CDNOMC
#endif
!
USE MODI_ABOR1_SFX
USE MODI_INIT_IO_SURF_n
USE MODI_READ_SURF
USE MODI_END_IO_SURF_n
USE MODI_IO_BUFF_CLEAN_n
USE MODI_UNPACK_SAME_RANK
!
USE YOMHOOK,            ONLY : LHOOK,DR_HOOK
USE PARKIND1,           ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER,            INTENT(IN)  :: KI
REAL,DIMENSION(KI), INTENT(IN)  :: PITM
REAL,DIMENSION(KI), INTENT(OUT) :: PSST
REAL,DIMENSION(KI), INTENT(OUT) :: PSIC  ! Not used at the moment
CHARACTER(LEN=2),   INTENT(IN)  :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
CHARACTER(LEN=6)     :: YPROGRAM2 = 'FA    '
REAL, DIMENSION(SIZE(PSST)) :: ZSST
REAL                 :: ZFMAX, ZFMIN, ZFMEAN
INTEGER              :: IRESP
REAL(KIND=JPRB)      :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ASSIM_SEA_N',0,ZHOOK_HANDLE)

IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_SEA_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
PSIC(:) = 0.
!
IF ( LREAD_SST_FROM_FILE ) THEN
  !
  !  Read SST from boundaries when SST analysis NOT is performed in CANARI
  !
  !  Define FA file name for SST analysis interpolated from boundary file 
  !
#ifdef FA
  CFILEIN_FA = 'SST_SIC'        ! input SST and SIC analysis  
  CDNOMC     = 'CADRE SST'      ! new frame name 
  WRITE(*,*) 'READING SST FROM ',TRIM(CFILEIN_FA)
#endif
!
!  Open FA file
!
  CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')
!
!  Read SST_SIC 
!
  IF ( LECSST ) THEN
  ! SST field interpolated from ECMWF SST ANALYSIS to model domain
    CALL READ_SURF(YPROGRAM2,'SURFSEA.TEMPERA',PSST,IRESP)
  ELSE
  ! Surface temperature from boundary in SST_SIC
    CALL READ_SURF(YPROGRAM2,'SURFTEMPERATURE',PSST,IRESP)
  ENDIF
!
!  Close SST_SIC file
!
  CALL END_IO_SURF_n(YPROGRAM2)
  CALL IO_BUFF_CLEAN_n
  WRITE(*,*) 'READ SST_SIC OK'

ELSE
  !
  IF ( NSIZE_SEA>0 ) THEN
    CALL UNPACK_SAME_RANK(NR_SEA,XSST,ZSST)
    PSST(:) = ZSST(:)
  ELSE
    PSST(:) = XUNDEF
  ENDIF
  !
ENDIF
!
ZFMIN = MINVAL(PSST)
ZFMAX = MAXVAL(PSST)
ZFMEAN = SUM(PSST)/FLOAT(KI)
!
IF ( LECSST ) THEN
  WRITE(*,*) '  ECMWF_SST_SIC'
  WRITE(*,'("  SURFSEA.TEMPERA - min, mean, max: ",3E13.4)') ZFMIN, ZFMEAN, ZFMAX
  ! Replace -9999. with UNDEF
  WHERE ( PSST(:)< 0. )
    PSST(:) = XUNDEF
  ENDWHERE
ELSE
  WRITE(*,*) '  Boundary file'
  WRITE(*,'("  SURFTEMPERATURE - min, mean, max: ",3E13.4)') ZFMIN, ZFMEAN, ZFMAX
  ! To avoid surface temperatures influenced by land, NATURE points are replaced with UNDEF
  WHERE ( PITM(:)>0.5 )
    PSST(:) = XUNDEF
  ENDWHERE
ENDIF
!
ZFMIN = MINVAL(PSST)
ZFMAX = MAXVAL(PSST)
ZFMEAN = SUM(PSST)/FLOAT(KI)
WRITE(*,*) '  Replaced land by UNDEF '
WRITE(*,'("  SST            - min, mean, max: ",3E13.4)') ZFMIN, ZFMEAN, ZFMAX

IF (LHOOK) CALL DR_HOOK('ASSIM_READ_SST_FROM_FILE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_READ_SST_FROM_FILE
