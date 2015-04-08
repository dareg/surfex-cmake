!     #####################
      MODULE MODD_CH_SEAFLUX_n
!     ######################
!
!!
!!    PURPOSE
!!    -------
!     
!   
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!
!!    AUTHOR
!!    ------
!!  P. Tulet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!  16/07/03 (P. Tulet)  restructured for externalization
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE CH_SEAFLUX_t
!
  CHARACTER(LEN=6)              :: CCH_DRY_DEP              !  deposition scheme
  REAL, DIMENSION(:,:), POINTER :: XDEP =>NULL()                    ! final dry deposition
                                                            ! velocity  for sea
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSV =>NULL()           ! name of the scalar var.
  INTEGER    :: NSV_CHSBEG, NSV_CHSEND                      ! chemical begin and ending
                                                            ! index of the HSV/CSV array
  INTEGER    :: NSV_DSTBEG, NSV_DSTEND                      ! dust begin and ending
  INTEGER    :: NSV_SLTBEG, NSV_SLTEND                      ! sea salt begin and ending
  INTEGER    :: NSV_AERBEG, NSV_AEREND                      ! aerosol begin and ending
  INTEGER    :: NBEQ                                        ! number of chemical species
  INTEGER    :: NDSTEQ                                      ! number of dust species
  INTEGER    :: NSLTEQ                                      ! number of sea salt species
  INTEGER    :: NAEREQ                                      ! number of aerosol species
                                                            ! in the surface scheme
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES=>NULL()      ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES=>NULL()
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES=>NULL()
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES=>NULL()
!
END TYPE CH_SEAFLUX_t

TYPE(CH_SEAFLUX_t), ALLOCATABLE, TARGET, SAVE :: CH_SEAFLUX_MODEL(:)

TYPE(CH_SEAFLUX_t), POINTER :: CH_SEAFLUX => NULL()
!$OMP THREADPRIVATE(CH_SEAFLUX)

CONTAINS

SUBROUTINE CH_SEAFLUX_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_CH_SEAFLUX_N:CH_SEAFLUX_GOTO_MODEL',0,ZHOOK_HANDLE)

CH_SEAFLUX => CH_SEAFLUX_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_CH_SEAFLUX_N:CH_SEAFLUX_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE CH_SEAFLUX_GOTO_MODEL

SUBROUTINE CH_SEAFLUX_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_SEAFLUX_N:CH_SEAFLUX_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(CH_SEAFLUX_MODEL(KMODEL))
CH_SEAFLUX_MODEL(:)%CCH_DRY_DEP=' '
CH_SEAFLUX_MODEL(:)%NBEQ=0
CH_SEAFLUX_MODEL(:)%NDSTEQ=0
CH_SEAFLUX_MODEL(:)%NSLTEQ=0
CH_SEAFLUX_MODEL(:)%NAEREQ=0
CH_SEAFLUX_MODEL(:)%NSV_CHSBEG=0
CH_SEAFLUX_MODEL(:)%NSV_CHSEND=0
CH_SEAFLUX_MODEL(:)%NSV_DSTBEG=0
CH_SEAFLUX_MODEL(:)%NSV_DSTEND=0
CH_SEAFLUX_MODEL(:)%NSV_SLTBEG=0
CH_SEAFLUX_MODEL(:)%NSV_SLTEND=0
CH_SEAFLUX_MODEL(:)%NSV_AERBEG=0
CH_SEAFLUX_MODEL(:)%NSV_AEREND=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_SEAFLUX_N:CH_SEAFLUX_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE CH_SEAFLUX_ALLOC

SUBROUTINE CH_SEAFLUX_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_SEAFLUX_N:CH_SEAFLUX_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(CH_SEAFLUX_MODEL)) DEALLOCATE(CH_SEAFLUX_MODEL)
IF (ASSOCIATED(CH_SEAFLUX)) NULLIFY(CH_SEAFLUX)
IF (LHOOK) CALL DR_HOOK("MODD_CH_SEAFLUX_N:CH_SEAFLUX_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE CH_SEAFLUX_DEALLO

END MODULE MODD_CH_SEAFLUX_n
