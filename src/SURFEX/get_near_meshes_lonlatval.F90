!     #########
      SUBROUTINE GET_NEAR_MESHES_LONLATVAL(KGRID_PAR,KL,PGRID_PAR,KNEAR_NBR,OLIST,KLIST,KNEAR)
!     ##############################################################
!
!!**** *GET_NEAR_MESHES_LONLATVAL* get the near grid mesh indices
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
INTEGER,                         INTENT(IN)    :: KNEAR_NBR ! number of nearest points wanted
REAL,    DIMENSION(KGRID_PAR),   INTENT(IN)    :: PGRID_PAR ! grid parameters
LOGICAL, DIMENSION(KL),          INTENT(IN)    :: OLIST     ! position in complete array of points
!                                                           ! for which one wants near mesh indices
INTEGER,                         INTENT(IN)    :: KLIST     ! number of points for which one 
                                                            ! wants near mesh indices
INTEGER, DIMENSION(KLIST,KNEAR_NBR),INTENT(OUT)   :: KNEAR     ! near mesh indices
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
REAL,DIMENSION(KL)    :: ZX
REAL,DIMENSION(KL)    :: ZY
REAL,DIMENSION(KL)    :: ZDX
REAL,DIMENSION(KL)    :: ZDY
INTEGER :: JLAT, JLON, JJ, JLIST
REAL, DIMENSION(KL,KL) :: XDIS
INTEGER :: IDIST
INTEGER :: ICOUNT
INTEGER :: IL
INTEGER, DIMENSION(1) :: ID0
REAL :: D0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_LONLATVAL',0,ZHOOK_HANDLE)
CALL GET_GRIDTYPE_LONLATVAL(PGRID_PAR,IL,ZX,ZY,ZDX,ZDY)
!
KNEAR  (:,:) = 0
!
IDIST = INT(SQRT(FLOAT(KNEAR_NBR)))
!
! calcul de la distance de tous les points 2 à 2
!
XDIS = 1.E20
DO JLON=1,KL
  IF (.NOT. OLIST(JLON)) CYCLE
  DO JLAT=1,KL
    XDIS(JLON,JLAT)=SQRT((ZX(JLON)-ZX(JLAT))**2+(ZY(JLON)-ZY(JLAT))**2)
  ENDDO
ENDDO
!
! on prend les knear_nbr premiers, pour chaque
JLIST = 0
DO JLON=1,KL
  IF (.NOT. OLIST(JLON)) CYCLE
  ICOUNT = 0
  JLIST = JLIST + 1
  DO JLAT=1,KNEAR_NBR
    ID0=MINLOC(XDIS(JLON,:))
    D0=MINVAL(XDIS(JLON,:))
    DO JJ=1,KL !on récupère celui qui a les indices minimaux
      IF (XDIS(JLON,JJ)==D0 .AND. &
            (ZX(JJ).LT.ZX(ID0(1)) .OR. &
             (ZX(JJ)==ZX(ID0(1)) .AND. ZY(JJ).LT.ZY(ID0(1))))) THEN  
        ID0(1)=JJ
      ENDIF
    ENDDO    
    !on le garde dans le tableau s'il répond au critère
    IF (abs(ZX(JLON)-ZX(ID0(1))).LE.(IDIST/2)*(ZDX(JLON)+ZDX(ID0(1)))/2 .AND. &
          abs(ZY(JLON)-ZY(ID0(1))).LE.(IDIST/2)*(ZDY(JLON)+ZDY(ID0(1)))/2) THEN  
      ICOUNT = ICOUNT + 1
      KNEAR(JLIST,ICOUNT)=ID0(1)
    ENDIF
    XDIS(JLON,ID0(1))=MAXVAL(XDIS(JLON,:))+1
    XDIS(JLON,ID0(1))=MAXVAL(XDIS(JLON,:))+1
  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_LONLATVAL',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_NEAR_MESHES_LONLATVAL
