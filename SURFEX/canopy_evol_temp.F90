!     #########
      SUBROUTINE CANOPY_EVOL_TEMP(KI,KLVL,PTSTEP,KIMPL,PTHA,PK,PDKDDVDZ,           &
                                   PSFLUX_T,PFORC_T,PDFORC_TDT,PDZ,PDZF,PTH,PWTH,   &
                                   PALFA,PBETA                                      ) 
!     #########################################
!
!!****  *CANOPY_EVOL_TEMP* - evolution of wind in canopy 
!!                        
!!
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2006 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_TRIDIAG_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,                  INTENT(IN)    :: KI        ! number of horizontal points
INTEGER,                  INTENT(IN)    :: KLVL      ! number of levels in canopy
REAL,                     INTENT(IN)    :: PTSTEP    ! time-step                             (s)
INTEGER,                  INTENT(IN)    :: KIMPL     ! implicitation code: 
!                                                    ! 1 : computes only alfa and beta coupling
!                                                    !     coefficients for all variables
!                                                    ! 2 : computes temporal evolution of the
!                                                    !     variables
REAL, DIMENSION(KI),      INTENT(IN)    :: PTHA      ! potential temp.  at forcing level     (K)
REAL, DIMENSION(KI,KLVL), INTENT(IN)    :: PK        ! mixing exchange coefficient           (m2/s)
REAL, DIMENSION(KI,KLVL), INTENT(IN)    :: PDKDDVDZ  ! derivative of mixing coefficient as a
!                                                    ! function of vertical gradient of wind
!                                                    ! (at mid levels)                       (m2)
REAL, DIMENSION(KI),      INTENT(IN)    :: PSFLUX_T  ! surface flux w'Th'                    (mK/s)
REAL, DIMENSION(KI,KLVL), INTENT(IN)    :: PFORC_T   ! tendency of wind due to canopy drag   (K/s)
REAL, DIMENSION(KI,KLVL), INTENT(IN)    :: PDFORC_TDT! formal derivative of the tendency of
!                                                    ! wind due to canopy drag               (1/s)
REAL, DIMENSION(KI,KLVL), INTENT(IN)    :: PDZ       ! deltaZ between canopy half levels,
!                                                    ! located at full levels                (m)
REAL, DIMENSION(KI,KLVL), INTENT(IN)    :: PDZF      ! deltaZ between canopy (full) levels,
!                                                    ! located at half levels                (m)
REAL, DIMENSION(KI,KLVL), INTENT(INOUT) :: PTH       ! pot. temp. at canopy levels           (K)
REAL, DIMENSION(KI,KLVL), INTENT(OUT)   :: PWTH      ! turbulent flux (at half levels)       (mK/s)
REAL, DIMENSION(KI),      INTENT(OUT)   :: PALFA     !  V+(1) = alfa F(1) + beta
REAL, DIMENSION(KI),      INTENT(OUT)   :: PBETA     !  V+(1) = alfa F(1) + beta
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER                     :: JLAYER   ! loop counter on layers
!
REAL, DIMENSION(KI,KLVL)   :: ZDTHDZ   ! dTh/dz at mid levels
REAL, DIMENSION(KI,KLVL)   :: ZF       ! turbulent flux at mid levels
REAL, DIMENSION(KI,KLVL)   :: ZDFDDVDZ ! derivative of turbulent flux as a
!                                       ! function of vertical gradient of wind variable
!                                       ! (at mid levels)
REAL, DIMENSION(KI,KLVL)   :: ZEXT     ! external forcing at full levels
REAL, DIMENSION(KI,KLVL)   :: ZDEXTDV  ! derivative of external forcing as a
!                                       ! function of vertical variable
!                                       ! (at full levels)
REAL, DIMENSION(KI,KLVL)   :: ZTH      ! work variable : pot. temp at futur instant 
!                                       ! (or past at the end of the routine) 
REAL, DIMENSION(KI)         :: ZDTHADT  ! dTHa/dt   at forcing level
REAL, DIMENSION(KI)         :: ZDWTHDZ  ! dw'Th'/dz at forcing level
LOGICAL                     :: LIMPL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!*    1. initializations
!        ---------------
!
!* external forces
!
IF (LHOOK) CALL DR_HOOK('CANOPY_EVOL_TEMP',0,ZHOOK_HANDLE)
ZEXT = 0.
ZDEXTDV = 0.
!
!* turbulent flux
!
ZF       = 0.
ZDFDDVDZ = 0.
!-------------------------------------------------------------------------------
!
!*    2. pot. temp. vertical derivative (at half levels below full levels)
!        ------------------------------
!
ZDTHDZ(:,:) = -999.
DO JLAYER=2,KLVL
  ZDTHDZ(:,JLAYER) = (PTH(:,JLAYER) - PTH(:,JLAYER-1)) / PDZF(:,JLAYER)
END DO
!
!-------------------------------------------------------------------------------
!
!*    3. turbulent flux (at half levels below full levels)
!        --------------
!
ZF(:,:) = -999.
!
!* surface flux
ZF(:,1) = PSFLUX_T(:)
!
!* flux at other levels
DO JLAYER=2,KLVL
  ZF(:,JLAYER) = -PK (:,JLAYER) * ZDTHDZ(:,JLAYER)
END DO
!
!-------------------------------------------------------------------------------
!
!*    4. formal derivative of turbulent flux for variable X=(dU/dz)
!        ----------------------------------------------------------
!
!* no implicitation of surface flux
!
ZDFDDVDZ(:,1) = 0.
!
!* other levels (at half levels below full levels)
!
ZDFDDVDZ(:,2:KLVL) = - PK(:,2:KLVL) - PDKDDVDZ(:,2:KLVL) * ZDTHDZ(:,2:KLVL)
!
!
!-------------------------------------------------------------------------------
!
!*    5. Forcing due to drag (at full levels)
!        -------------------
!
ZEXT = ZEXT + PFORC_T
ZDEXTDV = ZDEXTDV + PDFORC_TDT
!
!-------------------------------------------------------------------------------
!
!*    6. External forcing due to large-scale forces (at full levels)
!        ------------------------------------------
!
!
!* forces due to large-scale forcing
!
! These are computed from wind evolution equation at forcing level :
!
! dUa/dt = Large_Scale_Forcing - d(u'w')/dz|z=forcing_level
!
! because vertical derivative of turbulent flux is not available at forcing
! level, one must make the assumption that the turbulent flux is uniform
! between the forcing level and the level just below.
! This means that one assume that the forcing layer is in the inertial sublayer
! (where turbulent fluxes are constant).
!
ZDTHADT(:) = ( PTHA(:) - PTH(:,KLVL) ) /PTSTEP
!
ZDWTHDZ(:) = 0.
!
DO JLAYER=1,KLVL
  ZEXT(:,JLAYER) = ZEXT(:,JLAYER) + ZDWTHDZ(:) + ZDTHADT(:)
END DO
!
!-------------------------------------------------------------------------------
!
!*    7. adds Forces & divergence of turbulent flux to dU/dt (at full levels)
!        ---------------------------------------------------
!
LIMPL=(KIMPL==1)
CALL TRIDIAG_SURF(PTH,ZF,ZDFDDVDZ,ZEXT,ZDEXTDV,PTSTEP, &
                   PDZF,PDZ,ZTH,LIMPL,PALFA,PBETA    ) 
!
!
!-------------------------------------------------------------------------------
!
!*    8. updated turbulent flux (at half levels below full levels)
!        ----------------------
!
PWTH(:,:) = -999.
PWTH(:,1) = PSFLUX_T(:) 
DO JLAYER=2,KLVL
  PWTH(:,JLAYER) = PWTH(:,JLAYER-1)                                                         &
    +   ( PFORC_T(:,JLAYER-1) + PDFORC_TDT(:,JLAYER-1) * (ZTH(:,JLAYER-1)-PTH(:,JLAYER-1)) ) &
      * PDZ(:,JLAYER-1) 
END DO
!
!-------------------------------------------------------------------------------
!
!*    8. New value of potential temperature (at full levels)
!        ----------------------------------
!
!
PTH(:,:) = ZTH(:,:)
IF (LHOOK) CALL DR_HOOK('CANOPY_EVOL_TEMP',1,ZHOOK_HANDLE)
!
!----------------------------------------------------------------
!
END SUBROUTINE CANOPY_EVOL_TEMP
