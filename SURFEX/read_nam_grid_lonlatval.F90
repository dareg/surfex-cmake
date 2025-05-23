!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################################################################
      SUBROUTINE READ_NAM_GRID_LONLATVAL(PGRID_FULL_PAR,KDIM_FULL,HPROGRAM,KGRID_PAR,KL,PGRID_PAR,HDIR)
!     ################################################################
!
!!****  *READ_NAM_GRID_LONLATVAL* - routine to read in namelist the horizontal grid
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
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NSIZE_TASK
!
USE MODE_POS_SURF
!
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_GET_LUOUT
!
USE MODE_GRIDTYPE_LONLATVAL
USE MODI_READ_AND_SEND_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
REAL, DIMENSION(:), POINTER :: PGRID_FULL_PAR
INTEGER, INTENT(IN) :: KDIM_FULL
!
 CHARACTER(LEN=6),           INTENT(IN)    :: HPROGRAM   ! calling program
INTEGER,                    INTENT(INOUT) :: KGRID_PAR  ! size of PGRID_PAR
INTEGER,                    INTENT(OUT)   :: KL         ! number of points
REAL, DIMENSION(KGRID_PAR), INTENT(OUT)   :: PGRID_PAR  ! parameters defining this grid
 CHARACTER(LEN=1), INTENT(IN) :: HDIR
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: ILUOUT ! output listing logical unit
INTEGER :: ILUNAM ! namelist file  logical unit
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZX, ZX0       ! X conformal coordinate of grid mesh
REAL, DIMENSION(:),   ALLOCATABLE :: ZY, ZY0       ! Y conformal coordinate of grid mesh
REAL, DIMENSION(:),   ALLOCATABLE :: ZDX, ZDX0      ! X grid mesh size
REAL, DIMENSION(:),   ALLOCATABLE :: ZDY, ZDY0      ! Y grid mesh size
!
!*       0.3   Declarations of namelist
!              ------------------------
!
INTEGER :: NPOINTS  ! number of points
REAL, DIMENSION(100000) :: XX  ! X coordinate of grid mesh center (in meters)
REAL, DIMENSION(100000) :: XY  ! Y coordinate of grid mesh center (in meters)
REAL, DIMENSION(100000) :: XDX ! X mesh size (in meters)
REAL, DIMENSION(100000) :: XDY ! Y mesh size (in meters)
!
REAL, DIMENSION(:), POINTER :: ZGRID_PAR
!
LOGICAL :: GFOUND
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
NAMELIST/NAM_LONLATVAL/NPOINTS,XX,XY,XDX,XDY
!
!------------------------------------------------------------------------------
!
!*       1.    opening of namelist
! 
IF (LHOOK) CALL DR_HOOK('READ_NAM_GRID_LONLATVAL',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IF (HDIR/='H') THEN
  !
  CALL OPEN_NAMELIST(HPROGRAM,ILUNAM)
  !
  !---------------------------------------------------------------------------
  !
  !*       2.    Reading of projection parameters
  !              --------------------------------
  !
  CALL POSNAM(ILUNAM,'NAM_LONLATVAL',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_LONLATVAL)
  !
  !---------------------------------------------------------------------------
  CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
  !---------------------------------------------------------------------------
  !
  !*       3.    Number of points
  !              ----------------
  !
  KL = NPOINTS
  !
  !---------------------------------------------------------------------------
  !
  !*       3.    Array of X and Y coordinates
  !              ----------------------------
  !
  !
  ALLOCATE(ZX(KL))
  ALLOCATE(ZY(KL))
  ZX(:) = XX(:KL)
  ZY(:) = XY(:KL)
  !
  !---------------------------------------------------------------------------
  !
  !*       4.    Array of X and Y increments
  !              ---------------------------
  !
  ALLOCATE(ZDX(KL))
  ALLOCATE(ZDY(KL))
  ZDX(:) = XDX(:KL)
  ZDY(:) = XDY(:KL)
  !
ELSE
  !
  ALLOCATE(ZX0(KDIM_FULL),ZY0(KDIM_FULL),ZDX0(KDIM_FULL),ZDY0(KDIM_FULL))
  !
  CALL GET_GRIDTYPE_LONLATVAL(PGRID_FULL_PAR,&
                              PX=ZX0,PY=ZY0,PDX=ZDX0,PDY=ZDY0)
  !
  KL = NSIZE_TASK(NRANK)
  ALLOCATE(ZX(KL),ZY(KL),ZDX(KL),ZDY(KL))
  !
  CALL READ_AND_SEND_MPI(ZX0,ZX)
  CALL READ_AND_SEND_MPI(ZY0,ZY)
  CALL READ_AND_SEND_MPI(ZDX0,ZDX)
  CALL READ_AND_SEND_MPI(ZDY0,ZDY)
  !
  DEALLOCATE(ZX0,ZY0,ZDX0,ZDY0)  
  !
ENDIF
!---------------------------------------------------------------------------
!
!*       8.    All this information stored into pointer PGRID_PAR
!              --------------------------------------------------
!
 CALL PUT_GRIDTYPE_LONLATVAL(ZGRID_PAR,ZX,ZY,ZDX,ZDY)
!
!---------------------------------------------------------------------------
DEALLOCATE(ZX)
DEALLOCATE(ZY)
DEALLOCATE(ZDX)
DEALLOCATE(ZDY)
!---------------------------------------------------------------------------
!
!* 1st call : initializes dimension
!
IF (KGRID_PAR==0) THEN
  KGRID_PAR = SIZE(ZGRID_PAR)
!
ELSE
!
!* 2nd call : initializes grid array
!
  PGRID_PAR(:) = 0.
  PGRID_PAR(:) = ZGRID_PAR
END IF
!
DEALLOCATE(ZGRID_PAR)
IF (LHOOK) CALL DR_HOOK('READ_NAM_GRID_LONLATVAL',1,ZHOOK_HANDLE)
!
!---------------------------------------------------------------------------
!
END SUBROUTINE READ_NAM_GRID_LONLATVAL
