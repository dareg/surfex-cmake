!     #########
SUBROUTINE PREP_FLAKE_SBL()
!     #################################################################################
!
!!****  *PREP_FLAKE_SBL* - prepares FLAKE SBL fields
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
!!      E. Martin   01/2012 XUNDEF fields are no more written in PREP file
!!------------------------------------------------------------------
!
!
USE MODD_FLAKE_GRID_n, ONLY : FG => FLAKE_GRID
USE MODD_FLAKE_SBL_n, ONLY : FSB => FLAKE_SBL
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
IF (LHOOK) CALL DR_HOOK('PREP_FLAKE_SBL',0,ZHOOK_HANDLE)
FSB%NLVL = 6
!
!*      2.    height of half levels (where turbulent fluxes will be)
!             ---------------------
!
!* Warning :   ZZF(:,1)   MUST BE ZERO
ALLOCATE(ZZF(FG%NDIM,FSB%NLVL))
ZZF(:,1) = 0.
ZZF(:,2) = 1
ZZF(:,3) = 3.
ZZF(:,4) = 5.
ZZF(:,5) = 8.
ZZF(:,6) = 12.

ALLOCATE(FSB%XZ(FG%NDIM,FSB%NLVL))
DO JLAYER=1,FSB%NLVL-1
  FSB%XZ(:,JLAYER) = 0.5 * (ZZF(:,JLAYER)+ZZF(:,JLAYER+1))
END DO
FSB%XZ(:,FSB%NLVL) = 1.5 * ZZF(:,FSB%NLVL) - 0.5 * ZZF(:,FSB%NLVL-1)
!
DEALLOCATE(ZZF)
!
IF (LHOOK) CALL DR_HOOK('PREP_FLAKE_SBL',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_FLAKE_SBL
