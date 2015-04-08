!     #####################
      MODULE MODD_CH_FLAKE_n
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
!!      Modified    04/2013, P. Le Moigne: FLake chemistry
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

TYPE CH_FLAKE_t
!
  CHARACTER(LEN=6)              :: CCH_DRY_DEP              ! deposition scheme
  REAL, DIMENSION(:,:), POINTER :: XDEP =>NULL()                    ! final dry deposition
                                                            ! velocity  for lakes
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSV =>NULL()           ! name of the scalar var.
  INTEGER    :: NSV_CHSBEG, NSV_CHSEND                      ! chemical begin and ending
                                                            ! index of the HSV/CSV array
  INTEGER    :: NBEQ                                        ! number of chemical species
                                                            ! in the surface scheme

  INTEGER    :: NSV_DSTBEG, NSV_DSTEND                      ! dust begin and ending
  INTEGER    :: NSV_SLTBEG, NSV_SLTEND                      ! sea salt begin and ending
  INTEGER    :: NSV_AERBEG, NSV_AEREND                      ! aerosol begin and ending
  INTEGER    :: NDSTEQ                                      ! number of dust species
  INTEGER    :: NSLTEQ                                      ! number of sea salt species
  INTEGER    :: NAEREQ                                      ! number of aerosol species

  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES=>NULL()      ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES=>NULL()
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES=>NULL()
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES=>NULL()


!
END TYPE CH_FLAKE_t

TYPE(CH_FLAKE_t), ALLOCATABLE, TARGET, SAVE :: CH_FLAKE_MODEL(:)

TYPE(CH_FLAKE_t), POINTER :: CH_FLAKE => NULL()
!$OMP THREADPRIVATE(CH_FLAKE)

CONTAINS

SUBROUTINE CH_FLAKE_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!

! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_CH_FLAKE_N:CH_FLAKE_GOTO_MODEL',0,ZHOOK_HANDLE)

CH_FLAKE => CH_FLAKE_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_CH_FLAKE_N:CH_FLAKE_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE CH_FLAKE_GOTO_MODEL

SUBROUTINE CH_FLAKE_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_FLAKE_N:CH_FLAKE_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(CH_FLAKE_MODEL(KMODEL))
CH_FLAKE_MODEL(:)%CCH_DRY_DEP=' '
CH_FLAKE_MODEL(:)%NSV_CHSBEG=0
CH_FLAKE_MODEL(:)%NSV_CHSEND=0
CH_FLAKE_MODEL(:)%NSV_DSTBEG=0
CH_FLAKE_MODEL(:)%NSV_DSTEND=0
CH_FLAKE_MODEL(:)%NBEQ=0
CH_FLAKE_MODEL(:)%NSV_SLTBEG=0
CH_FLAKE_MODEL(:)%NSV_SLTEND=0
CH_FLAKE_MODEL(:)%NSV_AERBEG=0
CH_FLAKE_MODEL(:)%NSV_AEREND=0
CH_FLAKE_MODEL(:)%NDSTEQ=0
CH_FLAKE_MODEL(:)%NSLTEQ=0
CH_FLAKE_MODEL(:)%NAEREQ=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_FLAKE_N:CH_FLAKE_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE CH_FLAKE_ALLOC

SUBROUTINE CH_FLAKE_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_FLAKE_N:CH_FLAKE_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(CH_FLAKE_MODEL)) DEALLOCATE(CH_FLAKE_MODEL)
IF (ASSOCIATED(CH_FLAKE)) NULLIFY(CH_FLAKE)
IF (LHOOK) CALL DR_HOOK("MODD_CH_FLAKE_N:CH_FLAKE_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE CH_FLAKE_DEALLO

END MODULE MODD_CH_FLAKE_n
