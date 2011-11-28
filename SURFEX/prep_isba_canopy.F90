!     #########
SUBROUTINE PREP_ISBA_CANOPY()
!     #################################################################################
!
!!****  *PREP_ISBA_CANOPY* - prepares ISBA canopy fields
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
USE MODD_ISBA_GRID_n,     ONLY : NDIM
USE MODD_ISBA_CANOPY_n,   ONLY : NLVL, XZ, XU, XT, XQ, XTKE, XLMO, XP
USE MODD_SURF_PAR,        ONLY : XUNDEF
!
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
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZF    ! altitudes at half levels
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!
!*      1.    number of levels (MUST be at least equal to 2)
!             ----------------
!
IF (LHOOK) CALL DR_HOOK('PREP_ISBA_CANOPY',0,ZHOOK_HANDLE)
NLVL = 6
!NLVL = 13
!
!*      2.    height of half levels (where turbulent fluxes will be)
!             ---------------------
!
!* Warning :   ZZF(:,1)   MUST BE ZERO
ALLOCATE(ZZF(NDIM,NLVL))
ZZF(:,1) = 0.
ZZF(:,2) = 1
ZZF(:,3) = 3.
ZZF(:,4) = 5.
ZZF(:,5) = 8.
ZZF(:,6) = 12.
!ZZF(:,1) = 0.
!ZZF(:,2) = 0.2
!ZZF(:,3) = 1.
!ZZF(:,4) = 3.
!ZZF(:,5) = 6.
!ZZF(:,6) = 12.

!ZZF(:,1) = 0.
!ZZF(:,2) = 1.
!ZZF(:,3) = 2.
!ZZF(:,4) = 3.
!ZZF(:,5) = 4.
!ZZF(:,6) = 5.
!ZZF(:,7) = 6.
!ZZF(:,8) = 7.
!ZZF(:,9) = 8.
!ZZF(:,10) = 9.
!ZZF(:,11) = 10.
!ZZF(:,12) = 11.
!ZZF(:,13) = 12.


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
ALLOCATE(XLMO(NDIM))
XLMO(:) = XUNDEF
!
!*      8.    pressure (Pa)
!             --------
!
ALLOCATE(XP(NDIM,NLVL))
XP(:,:) = XUNDEF
IF (LHOOK) CALL DR_HOOK('PREP_ISBA_CANOPY',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_ISBA_CANOPY
