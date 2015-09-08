!     ####################
      MODULE MODD_DATA_COVER_n
!     ######################
!
!!****  *MODD_DATA_COVER_n* - declaration of correspondances between surface
!!                            classes and variables, for parameters that
!!                            can change as function of physical options
!!                            (GARDENs or not).
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
!!      V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       09/2011
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DATA_COVER_t
!
!-----------------------------------------------------------------------------------------------------
!
REAL, DIMENSION(:,:,:), POINTER :: XDATA_WEIGHT
!
REAL, DIMENSION(:),   POINTER :: XDATA_TOWN   ! artificial surfaces fraction
REAL, DIMENSION(:),   POINTER :: XDATA_NATURE ! natural and cul. fraction
REAL, DIMENSION(:),   POINTER :: XDATA_SEA    ! sea fraction
REAL, DIMENSION(:),   POINTER :: XDATA_WATER  ! inland water fraction
REAL, DIMENSION(:,:), POINTER :: XDATA_VEGTYPE! vegetation types fractions
REAL, DIMENSION(:),   POINTER :: XDATA_GARDEN ! garden fraction
REAL, DIMENSION(:),   POINTER :: XDATA_BLD    ! building fraction in
                                              ! artificial areas
REAL, DIMENSION(:),   POINTER :: XDATA_WALL_O_HOR  ! ratio of vert. surf.
!                                                  ! over hor. surf.
!
LOGICAL                           :: LGARDEN      ! T: define urban green areas
!                                                 ! F: no urban green areas
!
INTEGER :: NYEAR        ! current year for ecoclimap2
!
!-----------------------------------------------------------------------------------------------------
!
END TYPE DATA_COVER_t

TYPE(DATA_COVER_t), ALLOCATABLE, TARGET, SAVE :: DATA_COVER_MODEL(:)

TYPE(DATA_COVER_t), POINTER :: DATA_COVER => NULL()
!$OMP THREADPRIVATE(DATA_COVER)

CONTAINS

SUBROUTINE DATA_COVER_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DATA_COVER_N:DATA_COVER_GOTO_MODEL',0,ZHOOK_HANDLE)

DATA_COVER => DATA_COVER_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DATA_COVER_N:DATA_COVER_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE DATA_COVER_GOTO_MODEL

SUBROUTINE DATA_COVER_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_COVER_N:DATA_COVER_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DATA_COVER_MODEL(KMODEL))
DATA_COVER => DATA_COVER_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(DATA_COVER_MODEL(J)%XDATA_WEIGHT)
  NULLIFY(DATA_COVER_MODEL(J)%XDATA_TOWN)
  NULLIFY(DATA_COVER_MODEL(J)%XDATA_NATURE)
  NULLIFY(DATA_COVER_MODEL(J)%XDATA_SEA)
  NULLIFY(DATA_COVER_MODEL(J)%XDATA_WATER)
  NULLIFY(DATA_COVER_MODEL(J)%XDATA_VEGTYPE)
  NULLIFY(DATA_COVER_MODEL(J)%XDATA_GARDEN)
  NULLIFY(DATA_COVER_MODEL(J)%XDATA_BLD)
  NULLIFY(DATA_COVER_MODEL(J)%XDATA_WALL_O_HOR)
ENDDO
DATA_COVER_MODEL(:)%LGARDEN=.FALSE.
DATA_COVER_MODEL(:)%NYEAR=9999
IF (LHOOK) CALL DR_HOOK("MODD_DATA_COVER_N:DATA_COVER_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_COVER_ALLOC

SUBROUTINE DATA_COVER_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_COVER_N:DATA_COVER_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DATA_COVER_MODEL)) DEALLOCATE(DATA_COVER_MODEL)
IF (ASSOCIATED(DATA_COVER)) NULLIFY(DATA_COVER)
IF (LHOOK) CALL DR_HOOK("MODD_DATA_COVER_N:DATA_COVER_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_COVER_DEALLO

END MODULE MODD_DATA_COVER_n

