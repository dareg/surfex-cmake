!     #########
      SUBROUTINE WRITE_DIAG_SEB_SEAICE_n(HPROGRAM)
!     #################################
!
!!****  *WRITE_DIAG_SEB_SEAICE_n* - write the seaice diagnostic fields
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
!!	S.Senesi                *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2014
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIAG_SEAICE_n
USE MODD_SEAFLUX_n,  ONLY : XTICE, XICE_ALB, CSEAICE_SCHEME
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
REAL(KIND=JPRB)   :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_SEAICE_N',0,ZHOOK_HANDLE)
 CALL INIT_IO_SURF_n(HPROGRAM,'SEA   ','SEAFLX','WRITE')
!
!
!
!
 YCOMMENT='Sea-ice temperature (K)'
 CALL WRITE_SURF(HPROGRAM,'TSICE',XTICE(:),IRESP,YCOMMENT)
!
 YCOMMENT='Sea-ice albedo (-)'
 CALL WRITE_SURF(HPROGRAM,'IALB',XICE_ALB(:),IRESP,YCOMMENT)
!
 IF (TRIM(CSEAICE_SCHEME) == 'GELATO') THEN 
    YCOMMENT='Sea-ice thickness (m)'
    CALL WRITE_SURF(HPROGRAM,'SIT',XSIT(:),IRESP,YCOMMENT)
    !
    YCOMMENT='Sea-ice snow depth (m)'
    CALL WRITE_SURF(HPROGRAM,'SND',XSND(:),IRESP,YCOMMENT)
    !
    YCOMMENT='Sea mixed layer temp for Glt (K)'
    CALL WRITE_SURF(HPROGRAM,'SIMLT',XMLT(:),IRESP,YCOMMENT)
    !
 ENDIF
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)

IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_SEAICE_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE WRITE_DIAG_SEB_SEAICE_n
