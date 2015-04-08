!     #########
      SUBROUTINE WRITE_DIAG_SEB_OCEAN_n(HPROGRAM)
!     #################################
!
!!****  *WRITE_DIAG_SEB_OCEAN_n* - write the oceanic diagnostic fields
!!
!!    PURPOSE
!!    -------
!!
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
!!      C. Lebeaupin Brossier   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2007
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_OCEAN_n, ONLY : O => OCEAN
USE MODD_DIAG_OCEAN_n, ONLY : DGO => DIAG_OCEAN
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
!
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
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_OCEAN_N',0,ZHOOK_HANDLE)
 CALL INIT_IO_SURF_n(HPROGRAM,'SEA   ','SEAFLX','WRITE')
!
!
!*       2.     Mean values in OML :
!               --------------------
!
  YRECFM='TOML'
  YCOMMENT='X_Y_'//YRECFM
!
  CALL WRITE_SURF(HPROGRAM,YRECFM,DGO%XTOCMOY(:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='SOML'
  YCOMMENT='X_Y_'//YRECFM
!
  CALL WRITE_SURF(HPROGRAM,YRECFM,DGO%XSOCMOY(:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='UOML'
  YCOMMENT='X_Y_'//YRECFM
!
  CALL WRITE_SURF(HPROGRAM,YRECFM,DGO%XUOCMOY(:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='VOML'
  YCOMMENT='X_Y_'//YRECFM
!
  CALL WRITE_SURF(HPROGRAM,YRECFM,DGO%XVOCMOY(:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='DOML'
  YCOMMENT='X_Y_'//YRECFM
!
  CALL WRITE_SURF(HPROGRAM,YRECFM,DGO%XDOCMOY(:),IRESP,HCOMMENT=YCOMMENT)
!------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_OCEAN_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE WRITE_DIAG_SEB_OCEAN_n
