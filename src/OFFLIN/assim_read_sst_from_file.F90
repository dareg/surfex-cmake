!     ###############################################################################
SUBROUTINE ASSIM_READ_SST_FROM_FILE(YPROGRAM,KI,PSST,PITM,HTEST)

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
USE MODD_ASSIM,         ONLY : LECSST
USE MODD_SURF_PAR,      ONLY : XUNDEF
!
#ifdef FA
USE MODD_IO_SURF_FA,    ONLY : CFILEIN_FA,CDNOMC
#endif
!
USE YOMHOOK,            ONLY : LHOOK,DR_HOOK
USE PARKIND1,           ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_INIT_IO_SURF_n
USE MODI_READ_SURF
USE MODI_END_IO_SURF_n
USE MODI_IO_BUFF_CLEAN_n
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN)  :: YPROGRAM  ! program calling surf. schemes
INTEGER,            INTENT(IN)  :: KI
REAL,DIMENSION(KI), INTENT(INOUT) :: PSST
REAL,DIMENSION(KI), INTENT(IN)  :: PITM
CHARACTER(LEN=2),   INTENT(IN)  :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
CHARACTER(LEN=6)     :: YPROGRAM2 = 'FA    '
INTEGER              :: IRESP
REAL                 :: ZFMAX,ZFMIN,ZFMEAN
REAL, DIMENSION (KI) :: ZSST
REAL, DIMENSION (KI) :: PTS
REAL(KIND=JPRB)      :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ASSIM_SEA_N',0,ZHOOK_HANDLE)

IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_SEA_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF

!  Read SST from boundaries when SST analysis NOT is performed in CANARI
!
!  Define FA file name for SST analysis interpolated from boundary file 
!
#ifdef FA
CFILEIN_FA = 'SST_SIC'        ! input SST and SIC analysis  
CDNOMC     = 'CADRE SST'      ! new frame name 
#endif
!
WRITE(*,*) 'READING SST FROM ',TRIM(CFILEIN_FA)
!
!  Open FA file
!
CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')
!
!  Read SST_SIC 
!
IF ( LECSST ) THEN
  ! SST field interpolated from ECMWF SST ANALYSIS to model domain
  CALL READ_SURF(YPROGRAM2,'SURFSEA.TEMPERA',ZSST,IRESP)
ELSE
  ! Surface temperature from boundary in SST_SIC
  CALL READ_SURF(YPROGRAM2,'SURFTEMPERATURE',ZSST,IRESP)
ENDIF
!
!  Close SST_SIC file
!
CALL END_IO_SURF_n(YPROGRAM2)
CALL IO_BUFF_CLEAN_n
WRITE(*,*) 'READ SST_SIC OK'

ZFMIN = MINVAL(ZSST)
ZFMAX = MAXVAL(ZSST)
ZFMEAN = SUM(ZSST)/FLOAT(KI)

IF ( LECSST ) THEN
  WRITE(*,*) '  ECMWF_SST_SIC'
  WRITE(*,'("  SURFSEA.TEMPERA - min, mean, max: ",3E13.4)') ZFMIN, ZFMEAN, ZFMAX
  ! Replace -9999. with UNDEF
  WHERE ( ZSST(:)< 0. )
    ZSST(:) = XUNDEF
  END WHERE
ELSE
  WRITE(*,*) '  Boundary file'
  WRITE(*,'("  SURFTEMPERATURE - min, mean, max: ",3E13.4)') ZFMIN, ZFMEAN, ZFMAX
  ! To avoid surface temperatures influenced by land, NATURE points are replaced with UNDEF
  WHERE ( PTS(:)/=XUNDEF .OR. PITM(:)>0.5 )
    ZSST(:) = XUNDEF
  END WHERE
ENDIF

ZFMIN = MINVAL(ZSST)
ZFMAX = MAXVAL(ZSST)
ZFMEAN = SUM(ZSST)/FLOAT(KI)
WRITE(*,*) '  Replaced land by UNDEF '
WRITE(*,'("  ZSST            - min, mean, max: ",3E13.4)') ZFMIN, ZFMEAN, ZFMAX

IF (LHOOK) CALL DR_HOOK('ASSIM_READ_SST_FROM_FILE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_READ_SST_FROM_FILE
