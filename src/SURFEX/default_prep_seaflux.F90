!     #########
      SUBROUTINE DEFAULT_PREP_SEAFLUX
!     ###########################
!
!!****  *DEFAULT_PREP_SEAFLUX* - routine to set default values for the configuration for SEAFLUX field preparation
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
!!	S. Malardel   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2003 
!!      01/2008     C. Lebeaupin Brossier ! initialization of oceanic var. 
!!                                        ! from MERCATOR analyses types
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PREP_SEAFLUX,   ONLY : CFILE_SEAFLX, CTYPE, CFILEPGD_SEAFLX, CTYPEPGD, XSST_UNIF
!
USE MODN_PREP_SEAFLUX,   ONLY : LSEA_SBL, LOCEAN_MERCATOR, LOCEAN_CURRENT

USE MODD_SURF_PAR,   ONLY : XUNDEF
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
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DEFAULT_PREP_SEAFLUX',0,ZHOOK_HANDLE)
CFILE_SEAFLX = '                          '
CTYPE        = 'GRIB  '
!
CFILEPGD_SEAFLX = '                          '
CTYPEPGD        = '      '
!
XSST_UNIF = XUNDEF
!
LSEA_SBL = .FALSE.
LOCEAN_MERCATOR = .FALSE.
LOCEAN_CURRENT = .FALSE.
IF (LHOOK) CALL DR_HOOK('DEFAULT_PREP_SEAFLUX',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_PREP_SEAFLUX
