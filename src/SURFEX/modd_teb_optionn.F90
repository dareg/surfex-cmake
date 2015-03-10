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

LOGICAL, POINTER :: LCANOPY=>NULL()
!$OMP THREADPRIVATE(LCANOPY)
LOGICAL, POINTER :: LGARDEN=>NULL()
!$OMP THREADPRIVATE(LGARDEN)
 CHARACTER(LEN=4), POINTER :: CROAD_DIR=>NULL()
!$OMP THREADPRIVATE(CROAD_DIR)
 CHARACTER(LEN=4), POINTER :: CWALL_OPT=>NULL()
!$OMP THREADPRIVATE(CWALL_OPT)
 CHARACTER(LEN=3), POINTER :: CBLD_ATYPE=>NULL()
!$OMP THREADPRIVATE(CBLD_ATYPE)
 CHARACTER(LEN=6), POINTER :: CZ0H=>NULL()
!$OMP THREADPRIVATE(CZ0H)
 CHARACTER(LEN=5), POINTER :: CCH_BEM=>NULL()
!$OMP THREADPRIVATE(CCH_BEM)
 CHARACTER(LEN=3), POINTER :: CBEM=>NULL()
!$OMP THREADPRIVATE(CBEM)
 CHARACTER(LEN=3), POINTER :: CTREE=>NULL()
!$OMP THREADPRIVATE(CTREE)
LOGICAL, POINTER :: LGREENROOF=>NULL()
!$OMP THREADPRIVATE(LGREENROOF)
LOGICAL, POINTER :: LHYDRO=>NULL()
!$OMP THREADPRIVATE(LHYDRO)
LOGICAL, POINTER :: LSOLAR_PANEL=>NULL()
!$OMP THREADPRIVATE(LSOLAR_PANEL)
LOGICAL, POINTER :: LECOCLIMAP=>NULL()
!$OMP THREADPRIVATE(LECOCLIMAP)
REAL, POINTER, DIMENSION(:)   :: XZS=>NULL()
!$OMP THREADPRIVATE(XZS)
REAL, POINTER, DIMENSION(:,:) :: XCOVER=>NULL()
!$OMP THREADPRIVATE(XCOVER)
LOGICAL, POINTER, DIMENSION(:):: LCOVER=>NULL()
!$OMP THREADPRIVATE(LCOVER)
INTEGER, POINTER :: NTEB_PATCH=>NULL()
!$OMP THREADPRIVATE(NTEB_PATCH)
REAL, POINTER, DIMENSION(:,:) :: XTEB_PATCH=>NULL()
!$OMP THREADPRIVATE(XTEB_PATCH)
INTEGER, POINTER :: NROOF_LAYER=>NULL()
!$OMP THREADPRIVATE(NROOF_LAYER)
INTEGER, POINTER :: NROAD_LAYER=>NULL()
!$OMP THREADPRIVATE(NROAD_LAYER)
INTEGER, POINTER :: NWALL_LAYER=>NULL()
!$OMP THREADPRIVATE(NWALL_LAYER)
TYPE (DATE_TIME), POINTER :: TTIME=>NULL()
!$OMP THREADPRIVATE(TTIME)
REAL, POINTER :: XTSTEP=>NULL()
!$OMP THREADPRIVATE(XTSTEP)
REAL, POINTER :: XOUT_TSTEP=>NULL()
!$OMP THREADPRIVATE(XOUT_TSTEP)

!----------------------------------------------------------------------------
CONTAINS
!----------------------------------------------------------------------------

SUBROUTINE TEB_OPTIONS_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Save current state for allocated arrays
IF (LKFROM) THEN
TEB_OPTIONS_MODEL(KFROM)%XZS=>XZS
TEB_OPTIONS_MODEL(KFROM)%XCOVER=>XCOVER
TEB_OPTIONS_MODEL(KFROM)%LCOVER=>LCOVER
TEB_OPTIONS_MODEL(KFROM)%XTEB_PATCH=>XTEB_PATCH
ENDIF
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_N:TEB_OPTIONS_GOTO_MODEL',0,ZHOOK_HANDLE)
LCANOPY=>TEB_OPTIONS_MODEL(KTO)%LCANOPY
LGARDEN=>TEB_OPTIONS_MODEL(KTO)%LGARDEN
CZ0H=>TEB_OPTIONS_MODEL(KTO)%CZ0H
CCH_BEM=>TEB_OPTIONS_MODEL(KTO)%CCH_BEM
CBEM=>TEB_OPTIONS_MODEL(KTO)%CBEM
CTREE=>TEB_OPTIONS_MODEL(KTO)%CTREE
CROAD_DIR=>TEB_OPTIONS_MODEL(KTO)%CROAD_DIR
CWALL_OPT=>TEB_OPTIONS_MODEL(KTO)%CWALL_OPT
CBLD_ATYPE=>TEB_OPTIONS_MODEL(KTO)%CBLD_ATYPE
LGREENROOF=>TEB_OPTIONS_MODEL(KTO)%LGREENROOF
LHYDRO=>TEB_OPTIONS_MODEL(KTO)%LHYDRO
LSOLAR_PANEL=>TEB_OPTIONS_MODEL(KTO)%LSOLAR_PANEL
LECOCLIMAP=>TEB_OPTIONS_MODEL(KTO)%LECOCLIMAP
XZS=>TEB_OPTIONS_MODEL(KTO)%XZS
XCOVER=>TEB_OPTIONS_MODEL(KTO)%XCOVER
LCOVER=>TEB_OPTIONS_MODEL(KTO)%LCOVER
NTEB_PATCH=>TEB_OPTIONS_MODEL(KTO)%NTEB_PATCH
XTEB_PATCH=>TEB_OPTIONS_MODEL(KTO)%XTEB_PATCH
NROOF_LAYER=>TEB_OPTIONS_MODEL(KTO)%NROOF_LAYER
NROAD_LAYER=>TEB_OPTIONS_MODEL(KTO)%NROAD_LAYER
NWALL_LAYER=>TEB_OPTIONS_MODEL(KTO)%NWALL_LAYER
TTIME=>TEB_OPTIONS_MODEL(KTO)%TTIME
XTSTEP=>TEB_OPTIONS_MODEL(KTO)%XTSTEP
XOUT_TSTEP=>TEB_OPTIONS_MODEL(KTO)%XOUT_TSTEP
IF (LHOOK) CALL DR_HOOK('MODD_TEB_N:TEB_OPTIONS_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE TEB_OPTIONS_GOTO_MODEL

SUBROUTINE TEB_OPTIONS_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J,JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_OPTIONS_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_OPTIONS_MODEL(KMODEL))
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
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_OPTIONS_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_OPTIONS_DEALLO

!----------------------------------------------------------------------------
END MODULE MODD_TEB_OPTION_n
