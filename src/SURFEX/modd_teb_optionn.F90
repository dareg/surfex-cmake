!     ################
      MODULE MODD_TEB_OPTION_n
!     ################
!
!!****  *MODD_TEB_n - declaration of surface parameters for urban surface
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!!      A. Lemonsu      07/2012         Key for urban hydrology
!!      V. Masson       06/2013         splits module
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


TYPE TEB_OPTIONS_t
! TEB scheme option
!
  LOGICAL                        :: LCANOPY      ! T: SBL scheme within the canopy
                                                 ! F: no atmospheric layers below forcing level      
  LOGICAL                        :: LGARDEN      ! T: Urban green areas (call ISBA from TEB)
                                                 ! F: No urban green areas
  CHARACTER(LEN=4)               :: CROAD_DIR    ! TEB option for road directions
                                                 ! 'UNIF' : no specific direction
                                                 ! 'ORIE' : many road ORIEntations
                                                 ! ( one per TEB patch)
  CHARACTER(LEN=4)               :: CWALL_OPT    ! TEB option for walls
                                                 ! 'UNIF' : uniform walls
                                                 ! 'TWO ' : two separated walls
  CHARACTER(LEN=3)               :: CBLD_ATYPE   ! Type of averaging for walls
                                                 ! 'ARI'  : Characteristics are
                                                 !          linearly averaged
                                                 ! 'MAJ ' : Majoritary building in
                                                 !          grid mesh is chosen
  CHARACTER(LEN=6)               :: CZ0H         ! TEB option for z0h roof & road
                                                 ! 'MASC95' : Mascart et al 1995
                                                 ! 'BRUT82' : Brustaert     1982
                                                 ! 'KAND07' : Kanda         2007
  CHARACTER(LEN=5)               :: CCH_BEM      ! BEM option for roof/wall outside convective coefficient
                                                 ! 'DOE-2' : DOE-2 model from
                                                 ! EnergyPlus Engineering reference, p65
  CHARACTER(LEN=3)               :: CBEM         ! TEB option for the building energy model
                                                 ! 'DEF':  DEFault version force-restore model from Masson et al. 2002
                                                 ! 'BEM':  Building Energy Model Bueno et al. 2011

  CHARACTER(LEN=3)               :: CTREE        ! TEB option for the high vegetation
                                                 ! 'DEF':  DEFault version without radiative, dynamic effects or turbulent fluxes
                                                 ! 'RAD':  only RADiative effects 
                                                 ! 'DYN':  radiative and DYNamic effects 
                                                 ! 'FLX':  radiative, dynamic effects, and turbulent fluxes 
  LOGICAL                        :: LGREENROOF   ! T: green roofs (call ISBA from TEB)
  LOGICAL                        :: LHYDRO       ! T: urban subsoil and hydrology processes
  LOGICAL                        :: LSOLAR_PANEL ! T: solar panels on roofs
! 
! type of initialization of vegetation: from cover types (ecoclimap) or parameters prescribed
!
  LOGICAL                        :: LECOCLIMAP   ! T: parameters computed from ecoclimap
!                                                ! F: they are read in the file
!
! General surface: 
!
  REAL, POINTER, DIMENSION(:)   :: XZS           ! orography                        (m)
  REAL, POINTER, DIMENSION(:,:) :: XCOVER        ! fraction of each ecosystem       (-)
  LOGICAL, POINTER, DIMENSION(:):: LCOVER        ! GCOVER(i)=T --> ith cover field is not 0.
  INTEGER                       :: NTEB_PATCH    ! number of TEB patches
  REAL, POINTER, DIMENSION(:,:) :: XTEB_PATCH    ! fraction of each TEB patch
!
! Number of layers
!
  INTEGER                       :: NROOF_LAYER   ! number of layers in roofs
  INTEGER                       :: NROAD_LAYER   ! number of layers in roads
  INTEGER                       :: NWALL_LAYER   ! number of layers in walls
!
! Date:
!
  TYPE (DATE_TIME)              :: TTIME         ! current date and time
!
! Time-step:
!
  REAL                          :: XTSTEP        ! time step for TEB
!
  REAL                          :: XOUT_TSTEP    ! TEB output writing time step
!
END TYPE TEB_OPTIONS_t

TYPE(TEB_OPTIONS_t), ALLOCATABLE, TARGET, SAVE :: TEB_OPTIONS_MODEL(:)

TYPE(TEB_OPTIONS_t), POINTER :: TEB_OPTIONS => NULL()
!$OMP THREADPRIVATE(TEB_OPTIONS)

CONTAINS
!----------------------------------------------------------------------------

SUBROUTINE TEB_OPTIONS_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_N:TEB_OPTIONS_GOTO_MODEL',0,ZHOOK_HANDLE)

TEB_OPTIONS => TEB_OPTIONS_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_N:TEB_OPTIONS_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE TEB_OPTIONS_GOTO_MODEL

SUBROUTINE TEB_OPTIONS_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J,JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_OPTIONS_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_OPTIONS_MODEL(KMODEL))
TEB_OPTIONS => TEB_OPTIONS_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(TEB_OPTIONS_MODEL(J)%XZS)
  NULLIFY(TEB_OPTIONS_MODEL(J)%XCOVER)
  NULLIFY(TEB_OPTIONS_MODEL(J)%LCOVER)
  NULLIFY(TEB_OPTIONS_MODEL(J)%XTEB_PATCH)
ENDDO
TEB_OPTIONS_MODEL(:)%LCANOPY=.FALSE.
TEB_OPTIONS_MODEL(:)%LGARDEN=.FALSE.
TEB_OPTIONS_MODEL(:)%CROAD_DIR=' '
TEB_OPTIONS_MODEL(:)%CWALL_OPT=' '
TEB_OPTIONS_MODEL(:)%CBLD_ATYPE=' '
TEB_OPTIONS_MODEL(:)%CZ0H=' '
TEB_OPTIONS_MODEL(:)%CCH_BEM=' '
TEB_OPTIONS_MODEL(:)%CBEM=' '
TEB_OPTIONS_MODEL(:)%CTREE=' '
TEB_OPTIONS_MODEL(:)%LGREENROOF=.FALSE.
TEB_OPTIONS_MODEL(:)%LHYDRO=.FALSE.
TEB_OPTIONS_MODEL(:)%LSOLAR_PANEL=.FALSE.
TEB_OPTIONS_MODEL(:)%LECOCLIMAP=.FALSE.
TEB_OPTIONS_MODEL(:)%NTEB_PATCH=0
TEB_OPTIONS_MODEL(:)%NROOF_LAYER=0
TEB_OPTIONS_MODEL(:)%NROAD_LAYER=0
TEB_OPTIONS_MODEL(:)%NWALL_LAYER=0
TEB_OPTIONS_MODEL(:)%XTSTEP=0.
TEB_OPTIONS_MODEL(:)%XOUT_TSTEP=0.
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_OPTIONS_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_OPTIONS_ALLOC

SUBROUTINE TEB_OPTIONS_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_OPTIONS_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_OPTIONS_MODEL)) DEALLOCATE(TEB_OPTIONS_MODEL)
IF (ASSOCIATED(TEB_OPTIONS)) NULLIFY(TEB_OPTIONS)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_OPTIONS_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_OPTIONS_DEALLO

!----------------------------------------------------------------------------
END MODULE MODD_TEB_OPTION_n
