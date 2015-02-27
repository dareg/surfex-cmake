!     #########
      SUBROUTINE WRITESURF_CPL_GCM_n(HPROGRAM)
!     #######################################
!
!!****  *WRITESURF_CPL_GCM_n* - routine to write physical fields into
!!                              the restart file for ARPEGE/ALADIN run
!!
!!    PURPOSE
!!    -------
!!       The purpose of this routine is to store the 
!!       physical fields into the restart file . Indeed, 
!!       when ARPEGE/ALADIN is used, theses fields 
!!       are not initialized at the begin of a run.
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_ATM,      ONLY : LCPL_GCM
USE MODD_SURF_ATM_n,    ONLY : NSIZE_FULL, XRAIN, XSNOW, XZ0, XZ0H, XQSURF
!
USE MODI_WRITE_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! calling program
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
INTEGER           :: IRESP          ! Error code after redding
CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WRITESURF_PRECIP_N',0,ZHOOK_HANDLE)
!
YRECFM='LCPL_GCM'
YCOMMENT='flag to store physical fields in restart file'
CALL WRITE_SURF(HPROGRAM,YRECFM,LCPL_GCM,IRESP,HCOMMENT=YCOMMENT)
!
IF(LCPL_GCM)THEN
!
   YRECFM='RAIN_GCM'
   YCOMMENT='RAINFALL FOR RESTART (kg/m2/s)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XRAIN(:),IRESP,HCOMMENT=YCOMMENT)
!
   YRECFM='SNOW_GCM'
   YCOMMENT='SNOWFALL FOR RESTART (kg/m2/s)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XSNOW(:),IRESP,HCOMMENT=YCOMMENT)
!
   YRECFM='Z0_GCM'
   YCOMMENT='Z0 FOR RESTART (m)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XZ0(:),IRESP,HCOMMENT=YCOMMENT)
!
   YRECFM='Z0H_GCM'
   YCOMMENT='Z0H FOR RESTART (m)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XZ0H(:),IRESP,HCOMMENT=YCOMMENT)
!
   YRECFM='QS_GCM'
   YCOMMENT='QS FOR RESTART (kg/kg)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XQSURF(:),IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PRECIP_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_CPL_GCM_n
