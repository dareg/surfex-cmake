!     #########
      SUBROUTINE INIT_CPL_GCM_n(HPROGRAM,HINIT)
!     ########################################
!
!!****  *INIT_CPL_GCM_n* - routine to read  physical fields into  
!!                         the restart file for ARPEGE/ALADIN run
!!
!!    PURPOSE
!!    -------
!!       The purpose of this routine is to initialise some 
!!       physical fields. Indeed, when ARPEGE/ALADIN is used, 
!!       these field is not initialize at the begin of a run.
!!
!!
!!**  METHOD
!!    ------
!!      The data are read in the initial surface file :
!!        - 2D data fields
!!          
!!      It does not read the grid definition. This should have been
!!      read already.
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
!!	B. Decharme   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_ATM_n,      ONLY : NSIZE_FULL, XRAIN, XSNOW, XZ0, XZ0H, XQSURF
!
USE MODI_READ_SURF
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
 CHARACTER(LEN=3),  INTENT(IN)  :: HINIT    ! choice of fields to initialize
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP      ! Error code after redding
! 
CHARACTER(LEN=12) :: YRECFM     ! Name of the article to be read
!
INTEGER           :: IVERSION   ! surface version
INTEGER           :: IBUG       ! surface bugfix
!
LOGICAL           :: LREAD      ! work key
LOGICAL           :: LCPL_GCM   ! work key
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_CPL_GCM_N',0,ZHOOK_HANDLE)
!
LCPL_GCM = .FALSE.
!
IF (HINIT=='PGD') THEN
!     
   ALLOCATE(XRAIN (0))
   ALLOCATE(XSNOW (0))           
   ALLOCATE(XZ0   (0))           
   ALLOCATE(XZ0H  (0))           
   ALLOCATE(XQSURF(0))           
!     
ELSE
!
   ALLOCATE(XRAIN (NSIZE_FULL))
   ALLOCATE(XSNOW (NSIZE_FULL))
   ALLOCATE(XZ0   (NSIZE_FULL))
   ALLOCATE(XZ0H  (NSIZE_FULL))
   ALLOCATE(XQSURF(NSIZE_FULL))
!
   XRAIN (:) = 0.0
   XSNOW (:) = 0.0
   XZ0   (:) = 0.001
   XZ0H  (:) = 0.001
   XQSURF(:) = 0.0
!
ENDIF
!
YRECFM='VERSION'
CALL READ_SURF(HPROGRAM,YRECFM,IVERSION,IRESP)
YRECFM='BUG'
CALL READ_SURF(HPROGRAM,YRECFM,IBUG,IRESP)
!
LREAD=(HINIT/='PGD'.AND.HINIT/='PRE'.AND.(IVERSION>7.OR.(IVERSION==7.AND.IBUG>=3)))
!
IF (LREAD) THEN
   YRECFM='LCPL_GCM'
   CALL READ_SURF(HPROGRAM,YRECFM,LCPL_GCM,IRESP)
ENDIF  
!
IF (LREAD.AND.LCPL_GCM) THEN
!
   YRECFM='RAIN_GCM'
   CALL READ_SURF(HPROGRAM,YRECFM,XRAIN(:),IRESP)
!
   YRECFM='SNOW_GCM'
   CALL READ_SURF(HPROGRAM,YRECFM,XSNOW(:),IRESP)
!
   YRECFM='Z0_GCM'
   CALL READ_SURF(HPROGRAM,YRECFM,XZ0(:),IRESP)
!
   YRECFM='Z0H_GCM'
   CALL READ_SURF(HPROGRAM,YRECFM,XZ0H(:),IRESP)
!
   YRECFM='QS_GCM'
   CALL READ_SURF(HPROGRAM,YRECFM,XQSURF(:),IRESP)
!
ENDIF        
!
IF (LHOOK) CALL DR_HOOK('INIT_CPL_GCM_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_CPL_GCM_n
