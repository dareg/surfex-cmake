!     #####################
      MODULE MODD_CH_SURF_n
!     #####################
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
!!
!!    AUTHOR
!!    ------
!!  P. Tulet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!  16/07/03 (P. Tulet)  restructured for externalization
!!   10/2011 (S. Queguiner) Add CCH_EMIS
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

TYPE CH_SURF_t
!
  CHARACTER(LEN=4)              :: CCH_EMIS            ! Option for chemical emissions
                                                       ! 'NONE' : no emission
                                                       ! 'AGGR' : one aggregated value
                                                       !    for each specie and hour
                                                       ! 'SNAP' : from SNAP data using
                                                       !    potential emission & temporal profiles
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES ! NAME OF CHEMICAL
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES ! NAME OF AEROSOL SPECIES
                                                       ! SPECIES (FOR DIAG ONLY)
  CHARACTER(LEN=28)             :: CCHEM_SURF_FILE     ! name of general 
                                                       ! (chemical) purpose
                                                       ! ASCII input file
  REAL, DIMENSION(:), POINTER   :: XCONVERSION         ! emission unit 
                                                       ! conversion factor
  LOGICAL  :: LCH_SURF_EMIS                            ! T : chemical emissions
                                                       ! are used
  LOGICAL  :: LCH_EMIS                                 ! T : chemical emissions
                                                       ! are present in the file
!
END TYPE CH_SURF_t

TYPE(CH_SURF_t), ALLOCATABLE, TARGET, SAVE :: CH_SURF_MODEL(:)

TYPE(CH_SURF_t), POINTER :: CH_SURF => NULL()
!$OMP THREADPRIVATE(CH_SURF)

CONTAINS

SUBROUTINE CH_SURF_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_CH_SURF_N:CH_SURF_GOTO_MODEL',0,ZHOOK_HANDLE)

CH_SURF => CH_SURF_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_CH_SURF_N:CH_SURF_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE CH_SURF_GOTO_MODEL

SUBROUTINE CH_SURF_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_SURF_N:CH_SURF_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(CH_SURF_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(CH_SURF_MODEL(J)%CCH_NAMES)
  NULLIFY(CH_SURF_MODEL(J)%CAER_NAMES)
  NULLIFY(CH_SURF_MODEL(J)%XCONVERSION)
ENDDO
CH_SURF_MODEL(:)%CCH_EMIS=' '
CH_SURF_MODEL(:)%CCHEM_SURF_FILE=' '
CH_SURF_MODEL(:)%LCH_SURF_EMIS=.FALSE.
CH_SURF_MODEL(:)%LCH_EMIS=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_CH_SURF_N:CH_SURF_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE CH_SURF_ALLOC

SUBROUTINE CH_SURF_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_SURF_N:CH_SURF_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(CH_SURF_MODEL)) DEALLOCATE(CH_SURF_MODEL)
IF (ASSOCIATED(CH_SURF)) NULLIFY(CH_SURF)
IF (LHOOK) CALL DR_HOOK("MODD_CH_SURF_N:CH_SURF_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE CH_SURF_DEALLO

END MODULE MODD_CH_SURF_n
