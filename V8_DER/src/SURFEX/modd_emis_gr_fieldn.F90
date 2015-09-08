!     ###########################
      MODULE MODD_EMIS_GR_FIELD_n
!     ###########################
!
!!****  *MODD_EMIS_GR_FIELD_n* - declaration of chemical emission data arrays
!!                               for model n
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     chemical emission data arrays for model n.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_EMIS_GR_FIELD)
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

TYPE EMIS_GR_FIELD_t
!
!
  INTEGER                                     :: NEMIS_GR_NBR
!                          ! number of chemical pgd fields chosen by user
  CHARACTER(LEN=3) , DIMENSION(:), POINTER:: CEMIS_GR_AREA
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
  CHARACTER(LEN=40), DIMENSION(:), POINTER:: CEMIS_GR_NAME
!                          ! name of the chemical pgd fields (emitted species)
!
  INTEGER,       DIMENSION(:),     POINTER:: NEMIS_GR_TIME   ! emission time
!
  REAL,          DIMENSION(:,:,:), POINTER:: XEMIS_GR_FIELDS ! emission pgd fields values
!
  INTEGER                                          :: NEMISPEC_NBR ! Number of chemical species
!
  TYPE(EMISSVAR_T),  DIMENSION(:), POINTER :: TSEMISS      ! Offline emission struct array
!
  TYPE(PRONOSVAR_T),               POINTER     :: TSPRONOSLIST ! Head pointer on pronostic
!                                                              variables list
!-------------------------------------------------------------------------------
!
END TYPE EMIS_GR_FIELD_t

TYPE(EMIS_GR_FIELD_t), ALLOCATABLE, TARGET, SAVE :: EMIS_GR_FIELD_MODEL(:)

TYPE(EMIS_GR_FIELD_t), POINTER :: EMIS_GR_FIELD => NULL()
!$OMP THREADPRIVATE(EMIS_GR_FIELD)

CONTAINS

SUBROUTINE EMIS_GR_FIELD_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_EMIS_GR_FIELD_N:EMIS_GR_FIELD_GOTO_MODEL',0,ZHOOK_HANDLE)

EMIS_GR_FIELD => EMIS_GR_FIELD_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_EMIS_GR_FIELD_N:EMIS_GR_FIELD_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE EMIS_GR_FIELD_GOTO_MODEL

SUBROUTINE EMIS_GR_FIELD_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_EMIS_GR_FIELD_N:EMIS_GR_FIELD_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(EMIS_GR_FIELD_MODEL(KMODEL))
EMIS_GR_FIELD => EMIS_GR_FIELD_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(EMIS_GR_FIELD_MODEL(J)%CEMIS_GR_AREA)
  NULLIFY(EMIS_GR_FIELD_MODEL(J)%CEMIS_GR_NAME)
  NULLIFY(EMIS_GR_FIELD_MODEL(J)%NEMIS_GR_TIME)
  NULLIFY(EMIS_GR_FIELD_MODEL(J)%XEMIS_GR_FIELDS)
ENDDO
EMIS_GR_FIELD_MODEL(:)%NEMIS_GR_NBR=0
EMIS_GR_FIELD_MODEL(:)%NEMISPEC_NBR=0
IF (LHOOK) CALL DR_HOOK("MODD_EMIS_GR_FIELD_N:EMIS_GR_FIELD_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE EMIS_GR_FIELD_ALLOC

SUBROUTINE EMIS_GR_FIELD_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_EMIS_GR_FIELD_N:EMIS_GR_FIELD_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(EMIS_GR_FIELD_MODEL)) DEALLOCATE(EMIS_GR_FIELD_MODEL)
IF (ASSOCIATED(EMIS_GR_FIELD)) NULLIFY(EMIS_GR_FIELD)
IF (LHOOK) CALL DR_HOOK("MODD_EMIS_GR_FIELD_N:EMIS_GR_FIELD_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE EMIS_GR_FIELD_DEALLO

END MODULE MODD_EMIS_GR_FIELD_n
