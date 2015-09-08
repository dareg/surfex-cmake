!     #################
      MODULE MODD_OCEAN_n
!     #################
!
!!****  *MODD_OCEAN_n - declaration of ocean varaiables 
!!                          for 1D oceanic model
!!
!!    PURPOSE
!!    -------
!     Declaration of ocean varaiables
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
!!      C. Lebeaupin   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       04/2006
!!      Modified       07/2012, P. Le Moigne : CMO1D phasing
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE OCEAN_t
!
!
!   Switche for interactive coupling with oceanic model
LOGICAL:: LMERCATOR   !set to .true. to initialize oceanic var. from Mercator
LOGICAL:: LCURRENT    !set to .true. to make initialize ocean state with current      
LOGICAL:: LPROGSST    !set to .true. to make SST evolve with tendance
INTEGER:: NTIME_COUPLING! coupling time frequency 
INTEGER:: NOCTCOUNT   !oceanic model counter
! General surface: 
!
REAL, POINTER, DIMENSION(:,:) :: XSEAT  ! oceanic temperature profiles 
REAL, POINTER, DIMENSION(:,:) :: XSEAS  ! oceanic salinity profiles 
REAL, POINTER, DIMENSION(:,:) :: XSEAU  ! oceanic zonal current profiles 
REAL, POINTER, DIMENSION(:,:) :: XSEAV  ! oceanic meridian current profiles 
REAL, POINTER, DIMENSION(:,:) :: XSEAE  ! oceanic kinetic turbulent energy profiles (^(1/2))
REAL, POINTER, DIMENSION(:,:) :: XSEABATH !bathymetry indice
                                          !=1 for free sea water
                                          !=0 for sea-bed
REAL, POINTER, DIMENSION(:) ::   XSEAHMO! oceanic mixing lengths
!
REAL, POINTER, DIMENSION(:,:) :: XLE,XLK! oceanic mixing lengths
REAL, POINTER, DIMENSION(:,:) :: XKMEL,XKMELM  ! oceanic mixing coefficients
!
REAL, POINTER, DIMENSION(:) ::   XSEATEND! SST tendance
!
REAL, POINTER, DIMENSION(:,:) ::   XDTFSOL ! Temp tendancy due to solar flux
REAL, POINTER, DIMENSION(:)   ::   XDTFNSOL! -------------------- non solar flux
!
END TYPE OCEAN_t
!
TYPE(OCEAN_t), ALLOCATABLE, TARGET, SAVE :: OCEAN_MODEL(:)

TYPE(OCEAN_t), POINTER :: OCEAN => NULL()
!$OMP THREADPRIVATE(OCEAN)

CONTAINS

SUBROUTINE OCEAN_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_OCEAN_N:OCEAN_GOTO_MODEL',0,ZHOOK_HANDLE)

OCEAN => OCEAN_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_OCEAN_N:OCEAN_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE OCEAN_GOTO_MODEL

SUBROUTINE OCEAN_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_OCEAN_N:OCEAN_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(OCEAN_MODEL(KMODEL))
OCEAN => OCEAN_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(OCEAN_MODEL(J)%XSEAT)
  NULLIFY(OCEAN_MODEL(J)%XSEAS)
  NULLIFY(OCEAN_MODEL(J)%XSEAU)
  NULLIFY(OCEAN_MODEL(J)%XSEAV)
  NULLIFY(OCEAN_MODEL(J)%XSEAE)
  NULLIFY(OCEAN_MODEL(J)%XSEABATH)
  NULLIFY(OCEAN_MODEL(J)%XSEAHMO)
  NULLIFY(OCEAN_MODEL(J)%XLE)
  NULLIFY(OCEAN_MODEL(J)%XLK)
  NULLIFY(OCEAN_MODEL(J)%XKMEL)
  NULLIFY(OCEAN_MODEL(J)%XKMELM)
  NULLIFY(OCEAN_MODEL(J)%XSEATEND)
  NULLIFY(OCEAN_MODEL(J)%XDTFNSOL)
  NULLIFY(OCEAN_MODEL(J)%XDTFSOL)
ENDDO
OCEAN_MODEL(:)%LMERCATOR=.FALSE.
OCEAN_MODEL(:)%LCURRENT=.FALSE.
OCEAN_MODEL(:)%LPROGSST=.FALSE.
OCEAN_MODEL(:)%NTIME_COUPLING=0
OCEAN_MODEL(:)%NOCTCOUNT=0
IF (LHOOK) CALL DR_HOOK("MODD_OCEAN_N:OCEAN_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE OCEAN_ALLOC

SUBROUTINE OCEAN_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_OCEAN_N:OCEAN_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(OCEAN_MODEL)) DEALLOCATE(OCEAN_MODEL)
IF (ASSOCIATED(OCEAN)) NULLIFY(OCEAN)
IF (LHOOK) CALL DR_HOOK("MODD_OCEAN_N:OCEAN_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE OCEAN_DEALLO

END MODULE MODD_OCEAN_n
