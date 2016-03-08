!     #########################################
      SUBROUTINE CANOPY_GRID(KI,CP)
!     #########################################
!
!!****  *CANOPY_GRID* - computation of vertical grid coordinatesa at 
!!                      half levels and grid depths at half and full
!!                      levels
!!                        
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!
!
!  --------------------------------- PZ(k+1)                     PDZ(k+1)
!                                                                           ^
!                                                                           |
!                                                                           |
!  - - - - - - - - - - - - - - - - - PZf(k+1)                               | PDZf(k+1)
!                                                              ^            |
!                                                              |            |
!  --------------------------------- PZ(k), XU, XT, XQ, XTKE   | PDZ(k)     V
!                                                              |            ^
!  - - - - - - - - - - - - - - - - - PZf(k)                    V            | PDZf(k)
!  --------------------------------- PZ(k-1)                     PDZ(k-1)   V
!  - - - - - - - - - - - - - - - - - PZf(k-1)
!

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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2006 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CANOPY_n, ONLY : CANOPY_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,                  INTENT(IN)    :: KI     ! number of horizontal points
!
TYPE(CANOPY_t), INTENT(INOUT) :: CP
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JLAYER                 ! loop counter on layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CANOPY_GRID',0,ZHOOK_HANDLE)
!
!*    1. Geometric computations
!        ----------------------
!
!
!*    1.1 layer depths (variable located at half levels below full levels)
!         ------------
!
CP%XDZF(:,:) = -999.
CP%XDZF(:,1) = 2.*CP%XZ(:,1)
DO JLAYER=2,CP%NLVL
  CP%XDZF(:,JLAYER) = CP%XZ(:,JLAYER) - CP%XZ(:,JLAYER-1)
END DO
!
!*    1.2 Layer heights (variable located at half levels below full levels)
!         -------------
!
CP%XZF(:,:) = -999.
CP%XZF(:,1) = 0.
DO JLAYER=2,CP%NLVL
  CP%XZF(:,JLAYER) = 2.*CP%XZ(:,JLAYER-1) - CP%XZF(:,JLAYER-1)
END DO
!
!
!*    1.3 layer depths (variable located at full levels)
!         ------------
!
CP%XDZ(:,:) = -999.
DO JLAYER=1,CP%NLVL-1
  CP%XDZ(:,JLAYER) = CP%XZF(:,JLAYER+1) - CP%XZF(:,JLAYER)
END DO
!
CP%XDZ(:,CP%NLVL) = 2.*(CP%XZ(:,CP%NLVL)-CP%XZF(:,CP%NLVL))
!
IF (LHOOK) CALL DR_HOOK('CANOPY_GRID',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE CANOPY_GRID
