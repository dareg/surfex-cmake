!     ###########################
      MODULE MODD_CH_EMIS_FIELD_n
!     ###########################
!
!!****  *MODD_CH_EMIS_FIELD_n* - declaration of chemical emission data arrays
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     chemical emission data arrays.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!      D. Gazen   *L.A.*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/03/2001                      
!!      01/12/03    (D.Gazen) change emissions handling for surf. externalization
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_EFUTIL
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE CH_EMIS_FIELD_t
!
  REAL               :: XTIME_SIMUL  = 0.
  INTEGER            :: NTIME_MAX
  INTEGER            :: NEMIS_NBR
!                          ! number of chemical pgd fields chosen by user
  CHARACTER(LEN=3) , DIMENSION(:), POINTER :: CEMIS_AREA
!                          ! areas where chemical pgd fields are defined
!                          ! 'ALL' : everywhere
!                          ! 'SEA' : where sea exists
!                          ! 'LAN' : where land exists
!                          ! 'WAT' : where inland water exists
!                          ! 'NAT' : where natural or agricultural areas exist
!                          ! 'TWN' : where town areas exist
!                          ! 'STR' : where streets are present
!                          ! 'BLD' : where buildings are present
!                          !
  CHARACTER(LEN=40), DIMENSION(:), POINTER :: CEMIS_COMMENT ! comment
  CHARACTER(LEN=40), DIMENSION(:), POINTER :: CEMIS_NAME
!                          ! name of the chemical pgd fields (emitted species)
!
  INTEGER,           DIMENSION(:), POINTER :: NEMIS_TIME   ! emission time
!
  REAL,              DIMENSION(:,:), POINTER:: XEMIS_FIELDS ! emission pgd fields values
!
  INTEGER                                          :: NEMISPEC_NBR ! Number of chemical species
!
  TYPE(EMISSVAR_T),  DIMENSION(:), POINTER :: TSEMISS      ! Offline emission struct array
!
  TYPE(PRONOSVAR_T),               POINTER     :: TSPRONOSLIST ! Head pointer on pronostic
!                                                              variables list
!-------------------------------------------------------------------------------
!
END TYPE CH_EMIS_FIELD_t

TYPE(CH_EMIS_FIELD_t), ALLOCATABLE, TARGET, SAVE :: CH_EMIS_FIELD_MODEL(:)

TYPE(CH_EMIS_FIELD_t), POINTER :: CH_EMIS_FIELD => NULL()
!$OMP THREADPRIVATE(CH_EMIS_FIELD)

CONTAINS

SUBROUTINE CH_EMIS_FIELD_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_CH_EMIS_FIELD_N:CH_EMIS_FIELD_GOTO_MODEL',0,ZHOOK_HANDLE)

CH_EMIS_FIELD => CH_EMIS_FIELD_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_CH_EMIS_FIELD_N:CH_EMIS_FIELD_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE CH_EMIS_FIELD_GOTO_MODEL

SUBROUTINE CH_EMIS_FIELD_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_EMIS_FIELD_N:CH_EMIS_FIELD_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(CH_EMIS_FIELD_MODEL(KMODEL))
CH_EMIS_FIELD => CH_EMIS_FIELD_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(CH_EMIS_FIELD_MODEL(J)%CEMIS_AREA)
  NULLIFY(CH_EMIS_FIELD_MODEL(J)%CEMIS_COMMENT)
  NULLIFY(CH_EMIS_FIELD_MODEL(J)%CEMIS_NAME)
  NULLIFY(CH_EMIS_FIELD_MODEL(J)%NEMIS_TIME)
  NULLIFY(CH_EMIS_FIELD_MODEL(J)%XEMIS_FIELDS)
  NULLIFY(CH_EMIS_FIELD_MODEL(J)%TSEMISS)
ENDDO
CH_EMIS_FIELD_MODEL(:)%XTIME_SIMUL=0.
CH_EMIS_FIELD_MODEL(:)%NEMIS_NBR=0
CH_EMIS_FIELD_MODEL(:)%NTIME_MAX=-1
CH_EMIS_FIELD_MODEL(:)%NEMISPEC_NBR=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_EMIS_FIELD_N:CH_EMIS_FIELD_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE CH_EMIS_FIELD_ALLOC

SUBROUTINE CH_EMIS_FIELD_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_EMIS_FIELD_N:CH_EMIS_FIELD_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(CH_EMIS_FIELD_MODEL)) DEALLOCATE(CH_EMIS_FIELD_MODEL)
IF (ASSOCIATED(CH_EMIS_FIELD)) NULLIFY(CH_EMIS_FIELD)
IF (LHOOK) CALL DR_HOOK("MODD_CH_EMIS_FIELD_N:CH_EMIS_FIELD_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE CH_EMIS_FIELD_DEALLO

END MODULE MODD_CH_EMIS_FIELD_n

