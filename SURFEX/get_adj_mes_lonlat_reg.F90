!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE GET_ADJ_MES_LONLAT_REG(KGRID_PAR,KL,PGRID_PAR,KLEFT,KRIGHT,KTOP,KBOTTOM)
!     ##############################################################
!
!!**** *GET_ADJ_MES_LONLAT_REG* get the near grid mesh indices
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
INTEGER,                         INTENT(IN)    :: KGRID_PAR ! size of PGRID_PAR
INTEGER,                         INTENT(IN)    :: KL        ! number of points
REAL,    DIMENSION(KGRID_PAR),   INTENT(IN)    :: PGRID_PAR ! grid parameters
INTEGER, DIMENSION(KL),          INTENT(OUT)   :: KLEFT     ! left   mesh index
INTEGER, DIMENSION(KL),          INTENT(OUT)   :: KRIGHT    ! right  mesh index
INTEGER, DIMENSION(KL),          INTENT(OUT)   :: KTOP      ! top    mesh index
INTEGER, DIMENSION(KL),          INTENT(OUT)   :: KBOTTOM   ! bottom mesh index
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
INTEGER :: JLAT, JLON
INTEGER :: JL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_ADJ_MES_LONLAT_REG',0,ZHOOK_HANDLE)
 CALL GET_GRIDTYPE_LONLAT_REG(PGRID_PAR,ZLONMIN,ZLONMAX, &
                               ZLATMIN,ZLATMAX,ILON,ILAT  )  
!
KLEFT  (:) = 0
KRIGHT (:) = 0
KTOP   (:) = 0
KBOTTOM(:) = 0
!
IF (ILON*ILAT==KL) THEN
  DO JLAT=1,ILAT
    DO JLON=1,ILON
      JL = JLON + ILON * (JLAT-1)
      IF (JLON>1   ) KLEFT  (JL) = JL-1
      IF (JLON<ILON) KRIGHT (JL) = JL+1
      IF (JLAT>1   ) KBOTTOM(JL) = JL-ILON
      IF (JLAT<ILAT) KTOP   (JL) = JL+ILON
      IF (JLON==1    .AND. ZLONMIN + 360. == ZLONMAX) KLEFT  (JL) = JL+ILON-1
      IF (JLON==ILON .AND. ZLONMIN + 360. == ZLONMAX) KRIGHT (JL) = JL-ILON+1
    END DO
  END DO
END IF
IF (LHOOK) CALL DR_HOOK('GET_ADJ_MES_LONLAT_REG',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_ADJ_MES_LONLAT_REG
