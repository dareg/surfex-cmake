!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###############################################
      SUBROUTINE GET_GRID_COORD_LONLATVAL(KGRID_PAR,KL,PGRID_PAR,PX,PY)
!     ###############################################
!
!!****  *GET_GRID_COORD_LONLATVAL* - computes X and Y of ALL points
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
!!      E. Martin   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2007 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODE_GRIDTYPE_LONLATVAL
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,                    INTENT(IN)  :: KGRID_PAR  ! size of PGRID_PAR
INTEGER,                    INTENT(IN)  :: KL         ! number of points
REAL, DIMENSION(KGRID_PAR), INTENT(IN)  :: PGRID_PAR  ! parameters defining this grid
REAL, DIMENSION(KL),        INTENT(OUT) :: PX         ! X (m)
REAL, DIMENSION(KL),        INTENT(OUT) :: PY         ! Y (m)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!---------------------------------------------------------------------------
!
!*       1.    Projection and 2D grid parameters
!              ---------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_GRID_COORD_LONLATVAL',0,ZHOOK_HANDLE)
 CALL GET_GRIDTYPE_LONLATVAL(PGRID_PAR,PX=PX,PY=PY)
IF (LHOOK) CALL DR_HOOK('GET_GRID_COORD_LONLATVAL',1,ZHOOK_HANDLE)
!
!---------------------------------------------------------------------------
!
END SUBROUTINE GET_GRID_COORD_LONLATVAL
