!     #########
      SUBROUTINE WRITESURF_SEAFLUX_n(HPROGRAM)
!     ########################################
!
!!****  *WRITE_SEAFLUX_n* - writes SEAFLUX fields
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003 
!!      Modified    01/2014 : S. Senesi : handle seaice scheme
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SEAFLUX_n, ONLY : XSST, XZ0, TTIME, &
                           LINTERPOL_SST,    &
                           XSST_MTH,         &
                           LHANDLE_SIC
!
USE MODD_DIAG_SEAICE_n, ONLY: LDIAG_SEAICE
!
USE MODI_WRITE_SURF
USE MODI_WRITESURF_OCEAN_n
USE MODI_WRITESURF_SEAICE_N
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: JMTH, INMTH
 CHARACTER(LEN=2 ) :: YMTH
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_SEAFLUX_N',0,ZHOOK_HANDLE)
!
CALL WRITESURF_OCEAN_n(HPROGRAM)
!
!*       2.     Sea-ice prognostic fields:
!               --------------------------
!
!* flag to tell if Sea Ice model is used
!
YCOMMENT='flag to handle sea ice cover'
CALL WRITE_SURF(HPROGRAM,'HANDLE_SIC',LHANDLE_SIC,IRESP,YCOMMENT)
!
IF (LHANDLE_SIC) CALL WRITESURF_SEAICE_n(HPROGRAM)
!
!
!*       3.     Prognostic fields:
!               -----------------
!
!* water temperature
!
IF(LINTERPOL_SST)THEN
!
  INMTH=SIZE(XSST_MTH,2)
!
  DO JMTH=1,INMTH
     WRITE(YMTH,'(I2)') (JMTH-1)
     YRECFM='SST_MTH'//ADJUSTL(YMTH(:LEN_TRIM(YMTH)))
     YCOMMENT='SST month t'//ADJUSTL(YMTH(:LEN_TRIM(YMTH)))
     CALL WRITE_SURF(HPROGRAM,YRECFM,XSST_MTH(:,JMTH),IRESP,HCOMMENT=YCOMMENT)
  ENDDO
!
ENDIF
!
YRECFM='SST'
YCOMMENT='SST'
 CALL WRITE_SURF(HPROGRAM,YRECFM,XSST(:),IRESP,HCOMMENT=YCOMMENT)  
!
!-------------------------------------------------------------------------------
!
!*       4.     Semi-prognostic fields:
!               ----------------------
!
!* roughness length
!
YRECFM='Z0SEA'
YCOMMENT='Z0SEA (m)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,XZ0(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!*       5.  Time
!            ----
!
YRECFM='DTCUR'
YCOMMENT='s'
 CALL WRITE_SURF(HPROGRAM,YRECFM,TTIME,IRESP,HCOMMENT=YCOMMENT)
IF (LHOOK) CALL DR_HOOK('WRITESURF_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_SEAFLUX_n
