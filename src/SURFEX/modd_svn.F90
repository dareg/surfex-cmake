!        ###############
         MODULE MODD_SV_n
!        ###############
!
!!****  *MODD_NSV* - declaration of scalar variables numbers
!!
!!    PURPOSE
!!    -------
!!       Arrays to store the per-model NSV_* values number (suffix _A denote an array)
!!
!!    AUTHOR
!!    ------
!!      P. Tulet   Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original  01/2004
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE SV_t
!
!###############################################################################
!
! variables updated for the current model
!
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSV ! name of the scalar variables
  INTEGER    :: NSV_CHSBEG, NSV_CHSEND    !  index of first and last gas chemistry related scalar variable
  INTEGER    :: NBEQ                      ! number of chemical gas species in the surface scheme
  INTEGER    :: NSV_DSTBEG, NSV_DSTEND    ! index of first and last dust related scalar variable
  INTEGER    :: NDSTEQ                    ! number of dust related species in scalar variables list
  INTEGER    :: NSV_SLTBEG, NSV_SLTEND    ! index of first and last sea salt related scalar variable
  INTEGER    :: NSLTEQ                    ! number of sea salt related species in scalar variables list
  INTEGER    :: NSV_AERBEG, NSV_AEREND    ! index of first and last aerosol related scalar variabl
  INTEGER    :: NAEREQ                    ! number of aerosols variables

!
!
END TYPE SV_t

TYPE(SV_t), ALLOCATABLE, TARGET, SAVE :: SV_MODEL(:)

TYPE(SV_t), POINTER :: SV => NULL()
!$OMP THREADPRIVATE(SV)

CONTAINS

SUBROUTINE SV_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_SV_N:SV_GOTO_MODEL',0,ZHOOK_HANDLE)

SV => SV_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_SV_N:SV_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE SV_GOTO_MODEL

SUBROUTINE SV_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SV_N:SV_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(SV_MODEL(KMODEL))
SV => SV_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(SV_MODEL(J)%CSV)
ENDDO
SV_MODEL(:)%NBEQ=0
SV_MODEL(:)%NSV_CHSBEG=0
SV_MODEL(:)%NSV_CHSEND=0
SV_MODEL(:)%NSV_DSTBEG=0
SV_MODEL(:)%NSV_DSTEND=0
SV_MODEL(:)%NDSTEQ=0
SV_MODEL(:)%NSV_SLTBEG=0
SV_MODEL(:)%NSV_SLTEND=0
SV_MODEL(:)%NSLTEQ=0
SV_MODEL(:)%NSV_AERBEG=0
SV_MODEL(:)%NSV_AEREND=0
SV_MODEL(:)%NAEREQ=0
IF (LHOOK) CALL DR_HOOK("MODD_SV_N:SV_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE SV_ALLOC

SUBROUTINE SV_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SV_N:SV_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(SV_MODEL)) DEALLOCATE(SV_MODEL)
IF (ASSOCIATED(SV)) NULLIFY(SV)
IF (LHOOK) CALL DR_HOOK("MODD_SV_N:SV_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE SV_DEALLO

END MODULE MODD_SV_n
