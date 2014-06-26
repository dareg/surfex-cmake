!     #########
      SUBROUTINE DEFAULT_ASSIM(OASSIM,HASSIM,HASSIM_ISBA,NPRINTLEV,     &
                               OAROME,OECSST,OAESST,OAESNM,             &
                               OALADSURF,OREAD_SST_FROM_FILE,           &
                               OEXTRAP_SEA,OEXTRAP_WATER,OEXTRAP_NATURE,&
                               OWATERTG2,OBFIXED,PERROBS_M,KNCO,OOBSFILE)
!     ########################################################################
!
!!****  *DEFAULT_ISBA* - routine to set default values for the configuration for ISBA assimilation scheme
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
!!	L. Jarlan  *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2005
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_ASSIM, ONLY : NOBSMAX
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
LOGICAL,           INTENT(OUT) :: OASSIM        ! assimilation or not
CHARACTER(LEN=5),  INTENT(OUT) :: HASSIM        ! type of corrections PLUS/2DVAR
CHARACTER(LEN=5),  INTENT(OUT) :: HASSIM_ISBA
INTEGER,           INTENT(OUT) :: NPRINTLEV
LOGICAL,           INTENT(OUT) :: OAROME
LOGICAL,           INTENT(OUT) :: OECSST
LOGICAL,           INTENT(OUT) :: OAESST
LOGICAL,           INTENT(OUT) :: OAESNM
LOGICAL,           INTENT(OUT) :: OALADSURF
LOGICAL,           INTENT(OUT) :: OREAD_SST_FROM_FILE
LOGICAL,           INTENT(OUT) :: OEXTRAP_SEA
LOGICAL,           INTENT(OUT) :: OEXTRAP_WATER
LOGICAL,           INTENT(OUT) :: OEXTRAP_NATURE
LOGICAL,           INTENT(OUT) :: OWATERTG2
LOGICAL,           INTENT(OUT) :: OBFIXED
REAL, DIMENSION(NOBSMAX), INTENT(OUT) :: PERROBS_M
INTEGER, DIMENSION(NOBSMAX), INTENT(OUT) :: KNCO
LOGICAL, INTENT(OUT) :: OOBSFILE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DEFAULT_ASSIM',0,ZHOOK_HANDLE)
!
OASSIM    = .FALSE.
HASSIM    = "PLUS "
HASSIM_ISBA = "OI" 
NPRINTLEV = 0
OAROME    = .TRUE.
OECSST    = .FALSE.
OAESST    = .FALSE.
OAESNM    = .FALSE.
OALADSURF = .TRUE.
OREAD_SST_FROM_FILE=.FALSE.
OEXTRAP_SEA    = .TRUE.
OEXTRAP_WATER  = .TRUE.
OEXTRAP_NATURE = .FALSE.
OWATERTG2      = .FALSE.
OBFIXED = .FALSE.
PERROBS_M = (/1.0,0.1,0.4,0.2/)
KNCO = (/1,1,0,0/)
OOBSFILE = .FALSE.
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_ASSIM',1,ZHOOK_HANDLE)
!
END SUBROUTINE DEFAULT_ASSIM
