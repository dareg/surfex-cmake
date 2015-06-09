!     ###########################
      MODULE MODD_CH_SNAP_n
!     ###########################
!
!!****  *MODD_CH_SNAP_n* - declaration of chemical emission data arrays
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
!!      M.Leriche 04/2014  change length of CHARACTER for emission 6->12
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_EFUTIL
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE CH_EMIS_SNAP_t
!
  INTEGER            :: NEMIS_NBR
!                          ! number of chemical pgd fields chosen by user
  CHARACTER(LEN=3)                         :: CCONVERSION
!                          ! Unit conversion code
  CHARACTER(LEN=5)                         :: CSNAP_TIME_REF
!                          ! Reference time for Snap temporal profiles
!                          !  'UTC  ' : UTC   time
!                          !  'SOLAR' : SOLAR time
!                          !  'LEGAL' : LEGAL time
!                          !

  CHARACTER(LEN=12), DIMENSION(:), POINTER :: CEMIS_NAME
!                          ! name of the chemical fields (emitted species)
  CHARACTER(LEN=40), DIMENSION(:), POINTER :: CEMIS_COMMENT
!                          ! comment on the chemical fields (emitted species)
!
  REAL,     DIMENSION(:,:,:), POINTER:: XEMIS_FIELDS_SNAP ! Emission factor for
!                                                         ! each chemical specie and
!                                                         ! each snap
  REAL,     DIMENSION(:,:),   POINTER:: XEMIS_FIELDS      ! Emission for each specie
!                                                         ! (at a given time taking into 
!                                                         ! account all snaps)
  REAL,     DIMENSION(:),     POINTER:: XDELTA_LEGAL_TIME ! Difference (in hours)) between
!                                                         ! Legal time and UTC time
  INTEGER            :: NEMIS_SNAP                        ! number of snaps
  INTEGER            :: NSNAP_M                           ! number of months
  INTEGER            :: NSNAP_D                           ! number of days
  INTEGER            :: NSNAP_H                           ! number of hours
  REAL,              DIMENSION(:,:,:), POINTER:: XSNAP_MONTHLY
  REAL,              DIMENSION(:,:,:), POINTER:: XSNAP_DAILY
  REAL,              DIMENSION(:,:,:), POINTER:: XSNAP_HOURLY
  REAL,              DIMENSION(:),     POINTER:: XCONVERSION ! conversion factor
!
  TYPE(PRONOSVAR_T),               POINTER     :: TSPRONOSLIST ! Head pointer on pronostic
!                                                              variables list
!-------------------------------------------------------------------------------
!
END TYPE CH_EMIS_SNAP_t

TYPE(CH_EMIS_SNAP_t), ALLOCATABLE, TARGET, SAVE :: CH_EMIS_SNAP_MODEL(:)

TYPE(CH_EMIS_SNAP_t), POINTER :: CH_EMIS_SNAP => NULL()
!$OMP THREADPRIVATE(CH_EMIS_SNAP)

CONTAINS

SUBROUTINE CH_EMIS_SNAP_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_CH_SNAP_n:CH_EMIS_SNAP_GOTO_MODEL',0,ZHOOK_HANDLE)

CH_EMIS_SNAP => CH_EMIS_SNAP_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_CH_SNAP_n:CH_EMIS_SNAP_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE CH_EMIS_SNAP_GOTO_MODEL

SUBROUTINE CH_EMIS_SNAP_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_SNAP_n:CH_EMIS_FIELD_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(CH_EMIS_SNAP_MODEL(KMODEL))
CH_EMIS_SNAP => CH_EMIS_SNAP_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(CH_EMIS_SNAP_MODEL(J)%CEMIS_COMMENT)
  NULLIFY(CH_EMIS_SNAP_MODEL(J)%CEMIS_NAME)
  NULLIFY(CH_EMIS_SNAP_MODEL(J)%XDELTA_LEGAL_TIME)
  NULLIFY(CH_EMIS_SNAP_MODEL(J)%XEMIS_FIELDS)
  NULLIFY(CH_EMIS_SNAP_MODEL(J)%XEMIS_FIELDS_SNAP)
  NULLIFY(CH_EMIS_SNAP_MODEL(J)%XSNAP_DAILY)
  NULLIFY(CH_EMIS_SNAP_MODEL(J)%XSNAP_HOURLY)
  NULLIFY(CH_EMIS_SNAP_MODEL(J)%XSNAP_MONTHLY)
  NULLIFY(CH_EMIS_SNAP_MODEL(J)%XCONVERSION)
ENDDO
CH_EMIS_SNAP_MODEL(:)%CCONVERSION=' '
CH_EMIS_SNAP_MODEL(:)%CSNAP_TIME_REF=' '
CH_EMIS_SNAP_MODEL(:)%NEMIS_NBR=0
CH_EMIS_SNAP_MODEL(:)%NEMIS_SNAP=0
CH_EMIS_SNAP_MODEL(:)%NSNAP_M=0
CH_EMIS_SNAP_MODEL(:)%NSNAP_D=0
CH_EMIS_SNAP_MODEL(:)%NSNAP_H=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_SNAP_n:CH_EMIS_FIELD_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE CH_EMIS_SNAP_ALLOC

SUBROUTINE CH_EMIS_SNAP_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_SNAP_n:CH_EMIS_FIELD_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(CH_EMIS_SNAP_MODEL)) DEALLOCATE(CH_EMIS_SNAP_MODEL)
IF (ASSOCIATED(CH_EMIS_SNAP)) NULLIFY(CH_EMIS_SNAP)
IF (LHOOK) CALL DR_HOOK("MODD_CH_SNAP_n:CH_EMIS_FIELD_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE CH_EMIS_SNAP_DEALLO

END MODULE MODD_CH_SNAP_n

