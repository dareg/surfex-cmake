!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################################################################
      SUBROUTINE READ_NAM_GRID_LONLAT_ROT(PGRID_FULL_PAR,KDIM_FULL,HPROGRAM,KGRID_PAR,KL,PGRID_PAR,HDIR)
!     ################################################################
!
!!****  *READ_NAM_GRID_LONLAT_ROT* - routine to read in namelist the horizontal grid
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
!!      P. Samuelsson   SMHI
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/2012 
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
USE MODE_GRIDTYPE_LONLAT_ROT
!
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
REAL, DIMENSION(:), ALLOCATABLE :: ZLAT, ZLAT0 ! latitude  of all points
REAL, DIMENSION(:), ALLOCATABLE :: ZLON, ZLON0 ! longitude of all points
!
REAL,    DIMENSION(:), POINTER     :: ZGRID_PAR
!
LOGICAL :: GFOUND
!
!
!*       0.3   Declarations of namelist
!              ------------------------
!
REAL    :: XWEST   ! West longitude in rotated grid (degrees)
REAL    :: XSOUTH  ! South latitude in rotated grid  (degrees)
REAL    :: XDLON   ! Longitudal grid spacing  (degrees)
REAL    :: XDLAT   ! Latitudal grid spacing  (degrees)
REAL    :: XPOLON  ! Longitude of rotated pole (degrees)
REAL    :: XPOLAT  ! Latitude of rotated pole  (degrees)
INTEGER :: NLON    ! number of points in longitude
INTEGER :: NLAT    ! number of points in latitude
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
NAMELIST/NAM_LONLAT_ROT/XWEST,XSOUTH,XDLON,XDLAT,XPOLON,XPOLAT,NLON,NLAT
!
!------------------------------------------------------------------------------
!
!*       1.    opening of namelist
! 
IF (LHOOK) CALL DR_HOOK('READ_NAM_GRID_LONLAT_ROT',0,ZHOOK_HANDLE)
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
  CALL POSNAM(ILUNAM,'NAM_LONLAT_ROT',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_LONLAT_ROT)
  !
  !---------------------------------------------------------------------------
  !
  !*       3.    Number of points
  !              ----------------
  !
  KL = NLON * NLAT
  !
  !---------------------------------------------------------------------------
  CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
  !---------------------------------------------------------------------------
  !
  !*       4.    All this information stored into pointer PGRID_PAR
  !              --------------------------------------------------
  !
  ALLOCATE(ZLAT(KL))
  ALLOCATE(ZLON(KL))
  !
  CALL LATLON_LONLAT_ROT(XWEST,XSOUTH,XDLON,XDLAT,XPOLON,XPOLAT, &
                           NLON,NLAT,ZLON,ZLAT                   )  
  !
ELSE
  !
  ALLOCATE(ZLON0(KDIM_FULL),ZLAT0(KDIM_FULL))
  !
  CALL GET_GRIDTYPE_LONLAT_ROT(PGRID_FULL_PAR,PWEST=XWEST,PSOUTH=XSOUTH,&
                       PDLON=XDLON,PDLAT=XDLAT,PPOLON=XPOLON,PPOLAT=XPOLAT,&
                       KLON=NLON,KLAT=NLAT,KL=KL,PLON=ZLON0,PLAT=ZLAT0)
  !
  ALLOCATE(ZLON(KL),ZLAT(KL))
  !
  CALL READ_AND_SEND_MPI(ZLON0,ZLON)
  CALL READ_AND_SEND_MPI(ZLAT0,ZLAT)
  !
  DEALLOCATE(ZLON0,ZLAT0)
  !
ENDIF
!
 CALL PUT_GRIDTYPE_LONLAT_ROT(ZGRID_PAR,                                 &
                               XWEST,XSOUTH,XDLON,XDLAT,XPOLON,XPOLAT,  &
                               NLON,NLAT,KL,ZLON,ZLAT                   )  
!
DEALLOCATE(ZLAT)
DEALLOCATE(ZLON)
!
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
IF (LHOOK) CALL DR_HOOK('READ_NAM_GRID_LONLAT_ROT',1,ZHOOK_HANDLE)
!
!---------------------------------------------------------------------------
!
END SUBROUTINE READ_NAM_GRID_LONLAT_ROT
