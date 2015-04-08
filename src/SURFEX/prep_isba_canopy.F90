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
!!      E. Martin   01/2012 XUNDEF fields are no more written in PREP file
!!------------------------------------------------------------------
!
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
USE MODD_ISBA_CANOPY_n, ONLY : ICP => ISBA_CANOPY
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
ICP%NLVL = 6
!
!*      2.    height of half levels (where turbulent fluxes will be)
!             ---------------------
!
!* Warning :   ZZF(:,1)   MUST BE ZERO
ALLOCATE(ZZF(IG%NDIM,ICP%NLVL))
ZZF(:,1) = 0.
ZZF(:,2) = 1
ZZF(:,3) = 3.
ZZF(:,4) = 5.
ZZF(:,5) = 8.
ZZF(:,6) = 12.

ALLOCATE(ICP%XZ(IG%NDIM,ICP%NLVL))
DO JLAYER=1,ICP%NLVL-1
  ICP%XZ(:,JLAYER) = 0.5 * (ZZF(:,JLAYER)+ZZF(:,JLAYER+1))
END DO
ICP%XZ(:,ICP%NLVL) = 1.5 * ZZF(:,ICP%NLVL) - 0.5 * ZZF(:,ICP%NLVL-1)
!
DEALLOCATE(ZZF)
!
IF (LHOOK) CALL DR_HOOK('PREP_ISBA_CANOPY',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_ISBA_CANOPY
