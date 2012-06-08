!     #########
      SUBROUTINE GET_ADJ_MES_LONLATVAL(KGRID_PAR,KL,PGRID_PAR,KLEFT,KRIGHT,KTOP,KBOTTOM)
!     ##############################################################
!
!!**** *GET_ADJ_MES_LONLATVAL* get the near grid mesh indices
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
USE MODE_GRIDTYPE_LONLATVAL
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
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
REAL,DIMENSION(KL)    :: ZX
REAL,DIMENSION(KL)    :: ZY
REAL,DIMENSION(KL)    :: ZDX
REAL,DIMENSION(KL)    :: ZDY
INTEGER :: JLAT, JLON
INTEGER :: IL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_ADJ_MES_LONLATVAL',0,ZHOOK_HANDLE)
CALL GET_GRIDTYPE_LONLATVAL(PGRID_PAR,IL,ZX,ZY,ZDX,ZDY)
!
KLEFT  (:) = 0
KRIGHT (:) = 0
KTOP   (:) = 0
KBOTTOM(:) = 0
!
DO JLAT=1,KL
  DO JLON=1,KL
    IF (ZX(JLON)==ZX(JLAT)) THEN
      IF (ABS(ZY(JLON)-ZY(JLAT))==(ZDY(JLON)+ZDY(JLAT))/2.) THEN !si les points sont bien adjacents
        IF (ZY(JLON)<ZY(JLAT)) THEN
          KBOTTOM(JLAT)=JLON
        ELSEIF (ZY(JLON)>ZY(JLAT)) THEN
          KTOP(JLAT)=JLON
        ENDIF
      ENDIF
    ELSEIF (ZY(JLON)==ZY(JLAT)) THEN
      IF (ABS(ZX(JLON)-ZX(JLAT))==(ZDX(JLON)+ZDX(JLAT))/2.) THEN !si les points sont bien adjacents
        IF (ZX(JLON)<ZX(JLAT)) THEN
          KLEFT(JLAT)=JLON
        ELSEIF (ZX(JLON)>ZX(JLAT)) THEN
          KRIGHT(JLAT)=JLON
        ENDIF
      ENDIF
    ENDIF
  ENDDO
ENDDO  
IF (LHOOK) CALL DR_HOOK('GET_ADJ_MES_LONLATVAL',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_ADJ_MES_LONLATVAL
