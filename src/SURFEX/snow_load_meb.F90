!   ############################################################################
SUBROUTINE SNOW_LOAD_MEB(IR, DGEIP, PTSTEP, PSR, PWRVNMAX, PKVN, PCHEATV, PMELTVN, &
                         PVELC, PSUBVCOR)
!   ############################################################################
!
!!****  *SNOW_LOAD_MEB*
!!
!!    PURPOSE
!!    -------
!
!     Calculate temporal evolution of canopy-intercepted intercepted snow
!     
!!**  METHOD
!!    ------
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      P. Samuelsson           * SMHI *
!!      A. Boone                * CNRM-GAME, Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2011
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
!
USE MODD_CSTS,     ONLY : XTT, XLMTT, XLVTT, XLSTT
!
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    Declaration of Arguments
!
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIP
!
REAL,               INTENT(IN)    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)    :: PSR, PCHEATV, PVELC, PMELTVN, PWRVNMAX, PKVN
!
REAL, DIMENSION(:), INTENT(OUT)   :: PSUBVCOR
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSR))        :: ZSRINT, ZUNLOAD, ZWRVN, ZSUB
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.3    declarations of local parameters
!
! Snow unloading parameters (Roesch el al., Clim. Dyn., 2001)
!
REAL, PARAMETER                   :: ZUNLOAD_T     = 1.5E+5   ! K s
REAL, PARAMETER                   :: ZUNLOAD_TT    = 270.15   ! K
REAL, PARAMETER                   :: ZUNLOAD_V     = 1.87E+5  ! m 
!
!-------------------------------------------------
! 0) Initialization
!
IF (LHOOK) CALL DR_HOOK('SNOW_LOAD_MEB',0,ZHOOK_HANDLE)
!
ZSRINT(:)      = 0.0
ZWRVN(:)       = 0.0
ZSUB(:)        = 0.0
ZUNLOAD(:)     = 0.0
!
!
! 1) First consider the case when maximum interception is zero...
! this only occurs when vegetation canopy is *totally* buried. The follwing line
! results in non-zero snow loading (total removal of intercepted snow) 
! only during the timestep when vegetation has just been buried:
!
!
WHERE(PWRVNMAX(:) == 0.0)
!
   DGEIP%XSR_GN(:) = IR%XWR(:,1)/PTSTEP    ! kg m-2 s-1
   IR%XWR(:,1)       = 0.0

! for a totally buried canopy, the following are zero:

   DGEIP%XMELTCV(:)     = 0.0
   DGEIP%XFRZCV(:)      = 0.0
   PSUBVCOR(:)    = 0.0
!
!
ELSEWHERE
!
!
! 2) Case for snow beneath or only partially covering the vegetation canopy:
!
!
! The following are computed as steps to ensure mass conservation.
!
! Interception: gain

   ZSRINT(:)      = MAX(0.0,PWRVNMAX(:)-IR%XWR(:,1))*(1.0-EXP(-PKVN(:)*PSR(:)*PTSTEP)) ! kg m-2
   ZSRINT(:)      = MIN(PSR(:)*PTSTEP, ZSRINT(:))  ! kg m-2 
   ZWRVN(:)       = IR%XWR(:,1) + ZSRINT(:)           ! kg m-2 

   DGEIP%XSR_GN(:) = MAX(0.0, PSR(:) - ZSRINT(:)/PTSTEP) ! kg m-2 s-1

! Sublimation: gain or loss
! NOTE for the rare case that sublimation exceeds snow mass (possible as traces of snow disappear)
! compute a mass correction to be removed from soil (to conserve mass): PSUBVCOR

   ZSUB(:)        = DGEIP%XLESC(:)*(PTSTEP/XLSTT)              ! kg m-2
   PSUBVCOR(:)    = MAX(0.0, ZSUB(:) - ZWRVN(:))/PTSTEP  ! kg m-2 s-1
   ZWRVN(:)       = MAX(0.0, ZWRVN(:) - ZSUB(:))         ! kg m-2

! Phase change: loss (melt of snow mass)

   DGEIP%XMELTCV(:)     = PTSTEP*MAX(0.0, PMELTVN(:))         ! kg m-2  
   DGEIP%XMELTCV(:)     = MIN(DGEIP%XMELTCV(:), ZWRVN(:))
   ZWRVN(:)       = ZWRVN(:) - DGEIP%XMELTCV(:)
   IR%XWR(:,1)        = IR%XWR(:,1)  + DGEIP%XMELTCV(:)               ! NOTE...liq reservoir can exceed maximum holding
                                                        !        capacity here, but this is accounted for
                                                        !        in main prognostic PWRV routine.

! Phase change: gain (freeze of intercepted water) 
! Note, to get a better estimate of water available for freezing, remove Er in 
! estimation of water for freezing:
! Also, update liquid water stored on the canopy here:

   DGEIP%XFRZCV(:)      = PTSTEP*MAX(0.0, -PMELTVN(:))        ! kg m-2  
   DGEIP%XFRZCV(:)      = MIN(DGEIP%XFRZCV(:), MAX(0.0,IR%XWR(:,1)-DGEIP%XLERCV(:)*(PTSTEP/XLVTT)))
   ZWRVN(:)       = ZWRVN(:) + DGEIP%XFRZCV(:)
   IR%XWR(:,1)        = IR%XWR(:,1)  - DGEIP%XFRZCV(:)

! Unloading (falling off branches, etc...): loss
! Note, the temperature effect is assumed to vanish for cold temperatures.

   ZUNLOAD(:)     = MIN(ZWRVN(:), IR%XWR(:,1)*( PVELC(:)*(PTSTEP/ZUNLOAD_V)          &
                     + MAX(0.0, IR%XTV(:,1)-ZUNLOAD_TT)*(PTSTEP/ZUNLOAD_T) ))            ! kg m-2 
   ZWRVN(:)       = ZWRVN(:) - ZUNLOAD(:)                                           ! kg m-2 
   DGEIP%XSR_GN(:) = DGEIP%XSR_GN(:) + ZUNLOAD(:)/PTSTEP

! Diagnostic updates:
! final phase change (units)

   DGEIP%XMELTCV(:)     = DGEIP%XMELTCV(:)/PTSTEP ! kg m-2 s-1
   DGEIP%XFRZCV(:)      = DGEIP%XFRZCV(:) /PTSTEP ! kg m-2 s-1

! Prognostic Updates:

   IR%XWR(:,1)       = ZWRVN(:)

   IR%XTV(:,1)         = IR%XTV(:,1) + (DGEIP%XFRZCV(:) - DGEIP%XMELTCV(:))*(XLMTT*PTSTEP)/PCHEATV(:) ! K

END WHERE
!
IF (LHOOK) CALL DR_HOOK('SNOW_LOAD_MEB',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOW_LOAD_MEB
