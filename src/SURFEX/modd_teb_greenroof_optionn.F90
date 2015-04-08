!##################
MODULE MODD_TEB_GREENROOF_OPTION_n
!##################
!
!!****  *MODD_TEB_GREENROOF - declaration of ISBA scheme packed surface parameters for urban green roofs
!!
!!    PURPOSE
!!    -------
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
!!      A. Lemonsu *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       09/2009
!!      C. de Munck     06/2011 
!!      V. Masson       06/2013 splits module in 4
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_SNOW
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE




TYPE TEB_GREENROOF_OPTIONS_t
!-------------------------------------------------------------------------------
!
! type of initialization : from cover types (ecoclimap) or parameters prescribed
!
  LOGICAL                         :: LPAR_GREENROOF ! T: parameters computed from ecoclimap
!                                                   ! F: they are read in the file
!
! ISBA Scheme Options specific to urban green roofs:
!
  CHARACTER(LEN=3)                :: CISBA_GR       ! type of ISBA version ('2-L' = default, '3-L', 'DIF')
  CHARACTER(LEN=4)                :: CSCOND_GR      ! Thermal conductivity ('DEF '= NP89 implicit method , 
                                                    ! 'PL98' = Peters-Lidard et al. 1998 used for explicit computation of CG)
!
  LOGICAL                          :: LTR_ML_GR
!-------------------------------------------------------------------------------
!
! type of initialization of vegetation: from cover types (ecoclimap) or parameters prescribed
!
  INTEGER                         :: NLAYER_GR       ! number of ground layers
  INTEGER                         :: NTIME_GR        ! number of time data : for VEG, LAI, EMIS, Z0
!
  INTEGER                              :: NLAYER_HORT_GR
  INTEGER                              :: NLAYER_DUN_GR
!
  REAL, POINTER, DIMENSION(:)          :: XSOILGRID_GR        ! Soil layer grid as reference for DIF
!-------------------------------------------------------------------------------
!
! - SGH scheme
!                                                     
  CHARACTER(LEN=4)                :: CRUNOFF_GR      ! surface runoff formulation for green roofs
!                                                    ! 'WSAT'
!                                                    ! 'DT92'
!                                                    ! 'SGH ' Topmodel
!
!SGH scheme and vertical hydrology
!
  CHARACTER(LEN=3)                :: CKSAT_GR        ! ksat
!                                                    ! 'DEF' = default value 
!                                                    ! 'SGH' = profil exponentiel
  CHARACTER(LEN=3)                :: CHORT_GR        ! Horton runoff
!                                                    ! 'DEF' = no Horton runoff
!                                                    ! 'SGH' = Horton runoff
  LOGICAL                         :: LSOC_GR         ! soil organic carbon effect
!                                                    ! False = default value 
!                                                    ! True = SOC profil
!
!-------------------------------------------------------------------------------
!                                 
! Type of green roof (characterization of green roof structure based on GR vegetation)
!
  CHARACTER(LEN=5)                :: CTYP_GR         ! type of green roof
!
!-------------------------------------------------------------------------------
!
END TYPE TEB_GREENROOF_OPTIONS_t

TYPE(TEB_GREENROOF_OPTIONS_t), ALLOCATABLE, TARGET, SAVE :: TEB_GREENROOF_OPTIONS_MODEL(:)

TYPE(TEB_GREENROOF_OPTIONS_t), POINTER :: TEB_GREENROOF_OPTIONS => NULL()
!$OMP THREADPRIVATE(TEB_GREENROOF_OPTIONS)

CONTAINS

SUBROUTINE TEB_GREENROOF_OPTIONS_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_GREENROOF_N:TEB_GREENROOF_OPTIONS_GOTO_MODEL',0,ZHOOK_HANDLE)

TEB_GREENROOF_OPTIONS => TEB_GREENROOF_OPTIONS_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_GREENROOF_N:TEB_GREENROOF_OPTIONS_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE TEB_GREENROOF_OPTIONS_GOTO_MODEL

SUBROUTINE TEB_GREENROOF_OPTIONS_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_GREENROOF_OPTIONS_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(TEB_GREENROOF_OPTIONS_MODEL(J)%XSOILGRID_GR)
ENDDO
TEB_GREENROOF_OPTIONS_MODEL(:)%LPAR_GREENROOF=.TRUE.
TEB_GREENROOF_OPTIONS_MODEL(:)%CISBA_GR=' '
TEB_GREENROOF_OPTIONS_MODEL(:)%LTR_ML_GR=.FALSE.
TEB_GREENROOF_OPTIONS_MODEL(:)%LSOC_GR=.FALSE.
TEB_GREENROOF_OPTIONS_MODEL(:)%CRUNOFF_GR=' '
TEB_GREENROOF_OPTIONS_MODEL(:)%CSCOND_GR=' '
TEB_GREENROOF_OPTIONS_MODEL(:)%CKSAT_GR=' '
TEB_GREENROOF_OPTIONS_MODEL(:)%CHORT_GR=' '
TEB_GREENROOF_OPTIONS_MODEL(:)%CTYP_GR=' '
TEB_GREENROOF_OPTIONS_MODEL(:)%NLAYER_GR=0
TEB_GREENROOF_OPTIONS_MODEL(:)%NLAYER_HORT_GR=0
TEB_GREENROOF_OPTIONS_MODEL(:)%NLAYER_DUN_GR=0
TEB_GREENROOF_OPTIONS_MODEL(:)%NTIME_GR=0
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_OPTIONS_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GREENROOF_OPTIONS_ALLOC

SUBROUTINE TEB_GREENROOF_OPTIONS_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_OPTIONS_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_GREENROOF_OPTIONS_MODEL)) DEALLOCATE(TEB_GREENROOF_OPTIONS_MODEL)
IF (ASSOCIATED(TEB_GREENROOF_OPTIONS)) NULLIFY(TEB_GREENROOF_OPTIONS)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_OPTIONS_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GREENROOF_OPTIONS_DEALLO

END MODULE MODD_TEB_GREENROOF_OPTION_n
