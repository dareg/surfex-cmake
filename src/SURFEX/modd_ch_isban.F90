!     #####################
      MODULE MODD_CH_ISBA_n
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

TYPE CH_ISBA_t
!
  CHARACTER(LEN=28)  :: CCHEM_SURF_FILE  ! name of general (chemical) purpose ASCII input file
  CHARACTER(LEN=6)                :: CCH_DRY_DEP            !  deposition scheme
  REAL, DIMENSION(:,:,:), POINTER :: XDEP =>NULL()                  ! final dry deposition
                                                            ! velocity  for nature
  REAL, DIMENSION(:,:),   POINTER :: XSOILRC_SO2 =>NULL()           ! for SO2
  REAL, DIMENSION(:,:),   POINTER :: XSOILRC_O3 =>NULL()            ! for O3
  LOGICAL                         :: LCH_BIO_FLUX           ! flag for the calculation of
                                                            ! biogenic fluxes
  LOGICAL                         :: LCH_NO_FLUX            ! flag for the calculation of
                                                            ! biogenic NO fluxes
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSV =>NULL()           ! name of the scalar var.
  INTEGER    :: NSV_CHSBEG, NSV_CHSEND                      ! chemical begin and ending
                                                            ! index of the HSV/CSV array
  INTEGER    :: NBEQ                                        ! number of chemical species
                                                            ! in the surface scheme
  INTEGER    :: NSV_DSTBEG, NSV_DSTEND                      ! index of first and last dust
!                                                           ! related scalar variable
  INTEGER    :: NSV_SLTBEG, NSV_SLTEND                      ! index of first and last sea salt
!                                                           ! related scalar variable
  INTEGER    :: NSV_AERBEG, NSV_AEREND                      ! index of first and last dust
!                                                           ! related scalar variable
  INTEGER    :: NDSTEQ                                      ! number of dust related species
!                                                           ! in scalar variables list
  INTEGER    :: NAEREQ                                      ! number of aerosol species
  INTEGER    :: NSLTEQ                                      ! number of sea salt related species
!                                                           ! in scalar variables list
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES =>NULL()     ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES =>NULL()     ! NAME OF CHEMICAL SPECIES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES =>NULL()     ! NAME OF CHEMICAL SPECIES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES =>NULL()     ! NAME OF CHEMICAL SPECIES

!
END TYPE CH_ISBA_t

TYPE(CH_ISBA_t), ALLOCATABLE, TARGET, SAVE :: CH_ISBA_MODEL(:)

TYPE(CH_ISBA_t), POINTER :: CH_ISBA => NULL()
!$OMP THREADPRIVATE(CH_ISBA)

CONTAINS

SUBROUTINE CH_ISBA_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!

! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_CH_ISBA_N:CH_ISBA_GOTO_MODEL',0,ZHOOK_HANDLE)

CH_ISBA => CH_ISBA_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_CH_ISBA_N:CH_ISBA_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE CH_ISBA_GOTO_MODEL

SUBROUTINE CH_ISBA_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_ISBA_N:CH_ISBA_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(CH_ISBA_MODEL(KMODEL))
CH_ISBA_MODEL(:)%CCHEM_SURF_FILE=' '
CH_ISBA_MODEL(:)%CCH_DRY_DEP=' '
CH_ISBA_MODEL(:)%LCH_BIO_FLUX=.FALSE.
CH_ISBA_MODEL(:)%LCH_NO_FLUX=.FALSE.
CH_ISBA_MODEL(:)%NBEQ=0
CH_ISBA_MODEL(:)%NDSTEQ=0
CH_ISBA_MODEL(:)%NAEREQ=0
CH_ISBA_MODEL(:)%NSLTEQ=0
CH_ISBA_MODEL(:)%NSV_CHSBEG=0
CH_ISBA_MODEL(:)%NSV_CHSEND=0
CH_ISBA_MODEL(:)%NSV_DSTBEG=0
CH_ISBA_MODEL(:)%NSV_DSTEND=0
CH_ISBA_MODEL(:)%NSV_SLTBEG=0
CH_ISBA_MODEL(:)%NSV_SLTEND=0
CH_ISBA_MODEL(:)%NSV_AERBEG=0
CH_ISBA_MODEL(:)%NSV_AEREND=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_ISBA_N:CH_ISBA_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE CH_ISBA_ALLOC

SUBROUTINE CH_ISBA_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_ISBA_N:CH_ISBA_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(CH_ISBA_MODEL)) DEALLOCATE(CH_ISBA_MODEL)
IF (ASSOCIATED(CH_ISBA)) NULLIFY(CH_ISBA)
IF (LHOOK) CALL DR_HOOK("MODD_CH_ISBA_N:CH_ISBA_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE CH_ISBA_DEALLO

END MODULE MODD_CH_ISBA_n
