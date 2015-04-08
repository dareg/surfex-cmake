!     ####################
      MODULE MODD_SURF_ATM_n
!     ######################
!
!!****  *MODD_SURF_ATM - declaration of surface parameters
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
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
!!      V. Masson and A. Boone   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE SURF_ATM_t
!
!-----------------------------------------------------------------------------------------------------
!
! Type of each surface scheme
!
  CHARACTER(LEN=6) :: CTOWN      ! name of the urban surface scheme
  CHARACTER(LEN=6) :: CNATURE    ! name of the soil&vegetation surface scheme
  CHARACTER(LEN=6) :: CWATER     ! name of the scheme for inland water
  CHARACTER(LEN=6) :: CSEA       ! name for the ocean scheme
!
!-----------------------------------------------------------------------------------------------------
!
! Surface/Tile Fractions:
!
  REAL, POINTER, DIMENSION(:)   :: XTOWN     ! urban surface fraction of the grid box   (-)
  REAL, POINTER, DIMENSION(:)   :: XNATURE   ! natural surface fraction of the grid box (-)
  REAL, POINTER, DIMENSION(:)   :: XWATER    ! inland water fraction of the grid box    (-)
  REAL, POINTER, DIMENSION(:)   :: XSEA      ! sea/ocean fraction of the grid box       (-)
!
!-------------------------------------------------------------------------------
!
! type of initialization of vegetation: from cover types (ecoclimap) or parameters prescribed
!
  LOGICAL                        :: LECOCLIMAP ! T: parameters computed from ecoclimap
!                                              ! F: they are read in the file
!
!-------------------------------------------------------------------------------
!
! change water (not lake) to nature and/or town to rock : arrange cover properly
!
  LOGICAL                        :: LWATER_TO_NATURE ! T: Change Wetland treated as inland water into nature 
  LOGICAL                        :: LTOWN_TO_ROCK    ! T: Change Town into Rock
!
!-------------------------------------------------------------------------------
!
! include urban green areas for urbanized covers
!
  LOGICAL                        :: LGARDEN    ! T: define urban green areas
!                                              ! F: no urban green areas
!
!-----------------------------------------------------------------------------------------------------
!
! Masks and number of grid elements for each tile surface
!
! Sea/Ocean:
!
  INTEGER                               :: NSIZE_SEA    ! number of grid points by proc containing a
!                                                     ! sea surface                              (-)
  INTEGER                               :: NDIM_SEA     ! total number of grid points containing a
!                                                     ! sea surface                             (-)
  INTEGER, POINTER, DIMENSION(:)    :: NR_SEA       ! sea/ocean surface mask                  (-)
!
! Inland Water:
!
  INTEGER                               :: NSIZE_WATER  ! number of grid points containing an 
!                                                     ! inland water surface                    (-)
  INTEGER                               :: NDIM_WATER   ! total number of grid points by proc containing an
!                                                     ! inland surface
  INTEGER, POINTER, DIMENSION(:)    :: NR_WATER
!
! Town:
!
  INTEGER                               :: NSIZE_TOWN   ! number of grid points by proc containing an 
!                                                     ! urban surface                           (-)
  INTEGER                               :: NDIM_TOWN    ! total number of grid points containing an
!                                                     ! urban surface
  INTEGER, POINTER, DIMENSION(:)    :: NR_TOWN      ! urban surface mask                      (-)
!
! Natural surface:
!
  INTEGER                               :: NSIZE_NATURE ! number of grid points by proc containing a 
!                                                     ! natural surface                         (-)
  INTEGER                               :: NDIM_NATURE  ! total number of grid points containing a
!                                                     ! natural surface                         (-)
  INTEGER, POINTER, DIMENSION(:)    :: NR_NATURE    ! natural surface mask                    (-)
!
! All surfaces:
!
  INTEGER                               :: NSIZE_FULL   ! total number of grid points by proc     (-)
  INTEGER                               :: NDIM_FULL    ! total number of grid points             (-)
!
!-----------------------------------------------------------------------------------------------------
!
! Surface fields (only 1 horizontal dimension)
!
  REAL, POINTER, DIMENSION(:,:) :: XCOVER    ! fraction of each ecosystem for each grid box (-)
  LOGICAL, POINTER, DIMENSION(:):: LCOVER    ! GCOVER(i)=T --> ith cover field is not 0.
  REAL, POINTER, DIMENSION(:)   :: XZS       ! orography                                    (m)
!
!-------------------------------------------------------------------------------
!
  TYPE (DATE_TIME)                      :: TTIME            ! current date and time
!
  REAL                                  :: XOUT_TSTEP       ! output writing time step
!
!-----------------------------------------------------------------------------------------------------
!
! physical fields need into the restart file for ARPEGE/ALADIN run
!
  REAL, POINTER, DIMENSION(:)   :: XRAIN    ! Rainfall rate at surface               (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XSNOW    ! snowfall rate at surface               (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XZ0      ! surface roughness length for momentum  (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0H     ! surface roughness length for heat      (m)
  REAL, POINTER, DIMENSION(:)   :: XQSURF   ! specific humidity at surface           (kg/kg)
!
!-----------------------------------------------------------------------------------------------------
!
!
END TYPE SURF_ATM_t
!
TYPE(SURF_ATM_t), ALLOCATABLE, TARGET, SAVE :: SURF_ATM_MODEL(:)

TYPE(SURF_ATM_t), POINTER :: SURF_ATM => NULL()
!$OMP THREADPRIVATE(SURF_ATM)

CONTAINS

SUBROUTINE SURF_ATM_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_SURF_ATM_N:SURF_ATM_GOTO_MODEL',0,ZHOOK_HANDLE)

SURF_ATM => SURF_ATM_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_SURF_ATM_N:SURF_ATM_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE SURF_ATM_GOTO_MODEL
!
SUBROUTINE SURF_ATM_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_N:SURF_ATM_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(SURF_ATM_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(SURF_ATM_MODEL(J)%XTOWN)
  NULLIFY(SURF_ATM_MODEL(J)%XNATURE)
  NULLIFY(SURF_ATM_MODEL(J)%XWATER)
  NULLIFY(SURF_ATM_MODEL(J)%XSEA)
  NULLIFY(SURF_ATM_MODEL(J)%NR_SEA)
  NULLIFY(SURF_ATM_MODEL(J)%NR_WATER)
  NULLIFY(SURF_ATM_MODEL(J)%NR_TOWN)
  NULLIFY(SURF_ATM_MODEL(J)%NR_NATURE)
  NULLIFY(SURF_ATM_MODEL(J)%XCOVER)
  NULLIFY(SURF_ATM_MODEL(J)%LCOVER)
  NULLIFY(SURF_ATM_MODEL(J)%XZS)
  NULLIFY(SURF_ATM_MODEL(J)%XRAIN)
  NULLIFY(SURF_ATM_MODEL(J)%XSNOW)
  NULLIFY(SURF_ATM_MODEL(J)%XZ0)
  NULLIFY(SURF_ATM_MODEL(J)%XZ0H)
  NULLIFY(SURF_ATM_MODEL(J)%XQSURF)
ENDDO
SURF_ATM_MODEL(:)%CTOWN=' '
SURF_ATM_MODEL(:)%CNATURE=' '
SURF_ATM_MODEL(:)%CWATER=' '
SURF_ATM_MODEL(:)%CSEA=' '
SURF_ATM_MODEL(:)%LECOCLIMAP=.FALSE.
SURF_ATM_MODEL(:)%LWATER_TO_NATURE=.FALSE.
SURF_ATM_MODEL(:)%LTOWN_TO_ROCK=.FALSE.
SURF_ATM_MODEL(:)%LGARDEN=.FALSE.
SURF_ATM_MODEL(:)%NSIZE_SEA=0
SURF_ATM_MODEL(:)%NDIM_SEA=0
SURF_ATM_MODEL(:)%NSIZE_WATER=0
SURF_ATM_MODEL(:)%NDIM_WATER=0
SURF_ATM_MODEL(:)%NSIZE_TOWN=0
SURF_ATM_MODEL(:)%NDIM_TOWN=0
SURF_ATM_MODEL(:)%NSIZE_NATURE=0
SURF_ATM_MODEL(:)%NDIM_NATURE=0
SURF_ATM_MODEL(:)%NSIZE_FULL=0
SURF_ATM_MODEL(:)%NDIM_FULL=0
SURF_ATM_MODEL(:)%XOUT_TSTEP=0.
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_N:SURF_ATM_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE SURF_ATM_ALLOC

SUBROUTINE SURF_ATM_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_N:SURF_ATM_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(SURF_ATM_MODEL)) DEALLOCATE(SURF_ATM_MODEL)
IF (ASSOCIATED(SURF_ATM)) NULLIFY(SURF_ATM)
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_N:SURF_ATM_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE SURF_ATM_DEALLO

END MODULE MODD_SURF_ATM_n

