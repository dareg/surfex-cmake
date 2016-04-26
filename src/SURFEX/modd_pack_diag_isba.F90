!######################
MODULE MODD_PACK_DIAG_ISBA
!######################
!
!!****  *MODD_PACK_DIAG_ISBA - declaration of packed diagnostics for ISBA scheme
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!!      Modified       10/2004 by P. Le Moigne: add Halstead coefficient
!!      Modified       11/2009 by S. Senesi: add precipitation intercepted by the vegetation (XP_RRVEG)
!!      Modified       04-09 by A.L. Gibelin  : Add carbon diagnostics
!!      Modified       10-14 by P. Samuelsson: MEB
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!------------------------------------------------------------------------------
!
TYPE PACK_DIAG_ISBA_t

  INTEGER :: NSIZE_SIMPLE
  INTEGER :: NSIZE_GROUND
  INTEGER :: NSIZE_SNOW
  INTEGER :: NSIZE_KSW
  INTEGER :: NSIZE_ABC
  INTEGER :: NSIZE_0
  INTEGER :: NSIZE_00
  REAL, POINTER, DIMENSION(:,:) :: XBLOCK_SIMPLE
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_GROUND
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_SNOW
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_KSW
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_ABC
  REAL, POINTER, DIMENSION(:,:) :: XBLOCK_0
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_00
!
TYPE(DIAG_t) :: D
TYPE(DIAG_t) :: DP
!
TYPE(DIAG_MISC_ISBA_t) :: DMI
!
TYPE(DIAG_EVAP_ISBA_t) :: DE
TYPE(DIAG_EVAP_ISBA_t) :: DEP
!
TYPE(GR_BIOG_t) :: GB
!
END TYPE PACK_DIAG_ISBA_t
!
!-------------------------------------------------------------------------------
!


CONTAINS

!
!




SUBROUTINE PACK_DIAG_ISBA_INIT(YPACK_DIAG_ISBA)
TYPE(PACK_DIAG_ISBA_t), INTENT(INOUT) :: YPACK_DIAG_ISBA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YPACK_DIAG_ISBA%XBLOCK_SIMPLE)
  NULLIFY(YPACK_DIAG_ISBA%XBLOCK_GROUND)
  NULLIFY(YPACK_DIAG_ISBA%XBLOCK_SNOW)
  NULLIFY(YPACK_DIAG_ISBA%XBLOCK_KSW)
  NULLIFY(YPACK_DIAG_ISBA%XBLOCK_ABC)
  NULLIFY(YPACK_DIAG_ISBA%XBLOCK_0)
  NULLIFY(YPACK_DIAG_ISBA%XBLOCK_00)
IF (LHOOK) CALL DR_HOOK("MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE PACK_DIAG_ISBA_INIT


END MODULE MODD_PACK_DIAG_ISBA
