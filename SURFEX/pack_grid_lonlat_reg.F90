!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##############################################################
      SUBROUTINE PACK_GRID_LONLAT_REG(KMASK_SIZE,KMASK,KGRID_PAR1,PGRID_PAR1,KGRID_PAR2,OPACK,PGRID_PAR2)
!     ##############################################################
!
!!**** *PACK_GRID_LONLAT_REG* packs the grid definition vector
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson         Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    03/2004
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODI_PACK_SAME_RANK
USE MODE_GRIDTYPE_LONLAT_REG
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
INTEGER,                        INTENT(IN)    :: KMASK_SIZE ! size of mask
INTEGER, DIMENSION(KMASK_SIZE), INTENT(IN)    :: KMASK      ! mask used
INTEGER,                        INTENT(IN)    :: KGRID_PAR1 ! size of input grid vector
REAL,    DIMENSION(KGRID_PAR1), INTENT(IN)    :: PGRID_PAR1 ! parameters of input grid
INTEGER,                        INTENT(INOUT) :: KGRID_PAR2 ! size of output grid vector
LOGICAL,                        INTENT(IN)    :: OPACK      ! flag to pack the grid vector
REAL,    DIMENSION(KGRID_PAR2), INTENT(OUT)   :: PGRID_PAR2 ! parameters of output grid
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
REAL    :: ZLONMIN ! minimum longitude (degrees)
REAL    :: ZLONMAX ! maximum longitude (degrees)
REAL    :: ZLATMIN ! minimum latitude  (degrees)
REAL    :: ZLATMAX ! maximum latitude  (degrees)
INTEGER :: ILON    ! number of points in longitude
INTEGER :: ILAT    ! number of points in latitude
INTEGER :: IL      ! number of points used
REAL, DIMENSION(:), ALLOCATABLE    :: ZLAT1     ! latitude of all grid points
REAL, DIMENSION(:), ALLOCATABLE    :: ZLON1     ! longitude of all grid points
REAL, DIMENSION(:), ALLOCATABLE    :: ZLAT2     ! latitude of subset of grid points
REAL, DIMENSION(:), ALLOCATABLE    :: ZLON2     ! longitude of subset of grid points

!
REAL, DIMENSION(:), POINTER       :: ZGRID_PAR2 ! parameters of output grid
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
!*    1.     Computes grid parameters
!            ------------------------
!
IF (LHOOK) CALL DR_HOOK('PACK_GRID_LONLAT_REG',0,ZHOOK_HANDLE)
 CALL GET_GRIDTYPE_LONLAT_REG(PGRID_PAR1,ZLONMIN,ZLONMAX,     &
                               ZLATMIN,ZLATMAX,ILON,ILAT, IL   )  
ALLOCATE(ZLAT1(IL))
ALLOCATE(ZLON1(IL))
!
 CALL GET_GRIDTYPE_LONLAT_REG(PGRID_PAR1,PLON=ZLON1,PLAT=ZLAT1)
!----------------------------------------------------------------------------
!
!*    2.     Packs latitude and longitude arrays
!            -----------------------------------
!
!
ALLOCATE(ZLAT2(KMASK_SIZE))
ALLOCATE(ZLON2(KMASK_SIZE))
!
 CALL PACK_SAME_RANK(KMASK,ZLAT1,ZLAT2)
 CALL PACK_SAME_RANK(KMASK,ZLON1,ZLON2)
!
DEALLOCATE(ZLAT1)
DEALLOCATE(ZLON1)

!----------------------------------------------------------------------------
!
!*    3.     Stores data in new grid vector
!            ------------------------------
!
 CALL PUT_GRIDTYPE_LONLAT_REG(ZGRID_PAR2,ZLONMIN,ZLONMAX,                     &
                               ZLATMIN,ZLATMAX,ILON,ILAT,KMASK_SIZE,ZLON2,ZLAT2)  

DEALLOCATE(ZLAT2)
DEALLOCATE(ZLON2)
!----------------------------------------------------------------------------
!
IF (OPACK) THEN
  PGRID_PAR2(:) = ZGRID_PAR2(:)
ELSE
  KGRID_PAR2    = SIZE(ZGRID_PAR2(:))
END IF
!
DEALLOCATE(ZGRID_PAR2)
IF (LHOOK) CALL DR_HOOK('PACK_GRID_LONLAT_REG',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PACK_GRID_LONLAT_REG
