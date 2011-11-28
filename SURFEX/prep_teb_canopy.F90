!     #########
SUBROUTINE PREP_TEB_CANOPY()
!     #################################################################################
!
!!****  *PREP_TEB_CANOPY* - prepares TEB canopy fields
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2006
!!      S. Riette   06/2009 XT, XU, XQ, XTKE are set to XUNDEF
!!                          No more argument needed
!!------------------------------------------------------------------
!
!
USE MODD_TEB_GRID_n,     ONLY : NDIM
USE MODD_TEB_CANOPY_n,   ONLY : NLVL, XZ, XU, XT, XQ, XTKE, XLMO, XP, XQ0
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!*      0.2    declarations of local variables
!
INTEGER :: JLAYER
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZF ! altitudes at half levels
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!
!*      1.    number of levels (MUST be at least equal to 2)
!             ----------------
!
NLVL = 6
!NLVL = 7
IF (LHOOK) CALL DR_HOOK('PREP_TEB_CANOPY',0,ZHOOK_HANDLE)
!NLVL = 9
!
!*      2.    height of half levels (where turbulent fluxes will be)
!             ---------------------
!
!* Warning :   ZZF(:,1)   MUST BE ZERO
ALLOCATE(ZZF(NDIM,NLVL))
ZZF(:,1) = 0.
ZZF(:,2) = 1.
ZZF(:,3) = 3.
ZZF(:,4) = 5.
ZZF(:,5) = 8.
ZZF(:,6) = 12.

!ZZF(:,1) = 0.
!ZZF(:,2) = 3.
!ZZF(:,3) = 4.
!ZZF(:,4) = 10.
!ZZF(:,5) = 12.
!ZZF(:,6) = 14.
!ZZF(:,7) = 16.
!ZZF(:,8) = 20.
!ZZF(:,9) = 25.


!ZZF(:,1) = 0.
!ZZF(:,2) = 1.
!ZZF(:,3) = 2.
!ZZF(:,4) = 3.
!ZZF(:,5) = 5.
!ZZF(:,6) = 8.
!ZZF(:,7) = 12.

ALLOCATE(XZ(NDIM,NLVL))
DO JLAYER=1,NLVL-1
  XZ(:,JLAYER) = 0.5 * (ZZF(:,JLAYER)+ZZF(:,JLAYER+1))
END DO
XZ(:,NLVL) = 1.5 * ZZF(:,NLVL) - 0.5 * ZZF(:,NLVL-1)
!
DEALLOCATE(ZZF)
!
!
!*      3.    wind in canopy (m/s)
!             --------------
!
ALLOCATE(XU(NDIM,NLVL))
XU(:,:) = XUNDEF
!
!*      4.    temperature in canopy (K)
!             ---------------------
!
ALLOCATE(XT(NDIM,NLVL))
XT(:,:) = XUNDEF
!
!*      5.    humidity in canopy (kg/m3)
!             ------------------
!
ALLOCATE(XQ(NDIM,NLVL))
XQ(:,:) = XUNDEF
!
!*      6.    Tke in canopy (m2/s2)
!             -------------
!
ALLOCATE(XTKE(NDIM,NLVL))
XTKE(:,:) = XUNDEF
!
!*      7.    Monin-Obhukov length (m)
!             -------------
!
ALLOCATE(XLMO(NDIM,NLVL))
XLMO(:,:) = XUNDEF
!
!*      8.    pressure (Pa)
!             --------
!
ALLOCATE(XP(NDIM,NLVL))
XP(:,:) = XUNDEF
!
!*      9.    Surface flux (mK/s)
!             ------------
!
ALLOCATE(XQ0(NDIM))
XQ0(:) = XUNDEF
IF (LHOOK) CALL DR_HOOK('PREP_TEB_CANOPY',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_TEB_CANOPY
