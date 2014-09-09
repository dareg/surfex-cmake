!     #########
      SUBROUTINE GET_SFX_LAND(OCPL_GW,OCPL_FLOOD,OCPL_CALVING,  &
                              PRUNOFF,PDRAIN,PCALVING,PRECHARGE, &
                              PPFLOOD,PEFLOOD,PIFLOOD            )  
!     ###############################################################################
!
!!****  *GET_SFX_LAND* - routine to get some land surface variables from surfex
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	B. Decharme      *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODD_SURF_ATM_n, ONLY : NR_NATURE, NSIZE_NATURE
!
USE MODD_ISBA_n,     ONLY : LGLACIER,                    &
                            XCPL_DRAIN, XCPL_RUNOFF,     & 
                            XCPL_ICEFLUX, XCPL_RECHARGE, &
                            XCPL_PFLOOD, XCPL_EFLOOD,    &
                            XCPL_IFLOOD
!
USE MODI_UNPACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
LOGICAL,            INTENT(IN)  :: OCPL_GW     ! groundwater/surface key
LOGICAL,            INTENT(IN)  :: OCPL_FLOOD   ! flood key
LOGICAL,            INTENT(IN)  :: OCPL_CALVING ! calving key
!
REAL, DIMENSION(:), INTENT(OUT) :: PRUNOFF    ! Cumulated Surface runoff             (kg/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PDRAIN     ! Cumulated Deep drainage              (kg/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PCALVING   ! Cumulated Calving flux               (kg/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PRECHARGE  ! Cumulated Recharge to groundwater    (kg/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PPFLOOD    ! Cumulated flood precip interception  (kg/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PEFLOOD    ! Cumulated flood evaporation          (kg/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PIFLOOD    ! Cumulated flood infiltration         (kg/m2)
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(SIZE(PCALVING)) :: ZCALVING
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_SFX_LAND',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.0   Initialization
!              --------------
!
PRUNOFF  (:) = XUNDEF
PDRAIN   (:) = XUNDEF
PCALVING (:) = XUNDEF
PRECHARGE(:) = XUNDEF
PPFLOOD  (:) = XUNDEF
PEFLOOD  (:) = XUNDEF
PIFLOOD  (:) = XUNDEF
!
!*       2.0   Get variable over nature
!              ------------------------
!
IF(NSIZE_NATURE>0)THEN
!
! * surface runoff
!
  CALL UNPACK_SAME_RANK(NR_NATURE,XCPL_RUNOFF(:),PRUNOFF(:),XUNDEF)
  XCPL_RUNOFF (:) = 0.0
!
! * deep drainage
!
  CALL UNPACK_SAME_RANK(NR_NATURE,XCPL_DRAIN(:),PDRAIN(:),XUNDEF)
  XCPL_DRAIN(:) = 0.0
!
! * Calving flux
!
  IF(OCPL_CALVING)THEN
    CALL UNPACK_SAME_RANK(NR_NATURE,XCPL_ICEFLUX(:),PCALVING(:),XUNDEF)
    XCPL_ICEFLUX(:) = 0.0
  ELSEIF(LGLACIER)THEN
    XCPL_DRAIN  (:) = XCPL_DRAIN(:) + XCPL_ICEFLUX(:)
    XCPL_ICEFLUX(:) = 0.0
  ENDIF
!
! * groundwater recharge 
!
  IF(OCPL_GW)THEN
    CALL UNPACK_SAME_RANK(NR_NATURE,XCPL_RECHARGE(:),PRECHARGE(:),XUNDEF)
    XCPL_RECHARGE(:)=0.0
  ENDIF
!
! * floodplain source terms
!
  IF(OCPL_FLOOD)THEN
    CALL UNPACK_SAME_RANK(NR_NATURE,XCPL_PFLOOD(:),PPFLOOD(:),XUNDEF)
    CALL UNPACK_SAME_RANK(NR_NATURE,XCPL_EFLOOD(:),PEFLOOD(:),XUNDEF)
    CALL UNPACK_SAME_RANK(NR_NATURE,XCPL_IFLOOD(:),PIFLOOD(:),XUNDEF)
    XCPL_PFLOOD(:) = 0.0
    XCPL_EFLOOD(:) = 0.0
    XCPL_IFLOOD(:) = 0.0
  ENDIF
!
ENDIF
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_SFX_LAND',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_SFX_LAND
