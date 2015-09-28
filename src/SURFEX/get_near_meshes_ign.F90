!     #########
      SUBROUTINE GET_NEAR_MESHES_IGN(KGRID_PAR,KL,PGRID_PAR,KNEAR_NBR,KNEAR)
!     ##############################################################
!
!!**** *GET_NEAR_MESHES_IGN* get the near grid mesh indices
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
USE MODE_GRIDTYPE_IGN
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
INTEGER, DIMENSION(:,:),POINTER :: KNEAR    ! near mesh indices
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
REAL, DIMENSION(KL,KNEAR_NBR) :: ZNEAR
REAL,DIMENSION(KL)  :: ZX
REAL,DIMENSION(KL)  :: ZY
REAL,DIMENSION(KL)  :: ZDX
REAL,DIMENSION(KL)  :: ZDY
REAL :: ZDIS
INTEGER :: JP1, JP2, JN1, JN2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_IGN_1',0,ZHOOK_HANDLE)
!
 CALL GET_GRIDTYPE_IGN(PGRID_PAR,PX=ZX,PY=ZY,PDX=ZDX,PDY=ZDY)
!
KNEAR(:,:) = 0
ZNEAR(:,:) = 0.
!
! calcul de la distance de tous les points 2 à 2
!
!
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_IGN_1',1,ZHOOK_HANDLE)
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_IGN_2',0,ZHOOK_HANDLE)
!
!$OMP PARALLEL DO PRIVATE(JP1,JP2,ZDIS)
DO JP1=1,KL-1
  !
  DO JP2=JP1+1,KL
    !
    ZDIS = SQRT((ZX(JP1)-ZX(JP2))**2+(ZY(JP1)-ZY(JP2))**2)
    !
    IF (JP1==1) THEN 
      !
      ZNEAR(JP2,1) = ZDIS
      KNEAR(JP2,1) = 1
      !
      IF (JP2==2) THEN
        !
        ZNEAR(1,1) = ZDIS
        KNEAR(1,1) = 2
        !
      ELSE
        !
        CALL GET_NEAR_POINTS(JP1,JP2,ZDIS)
        !
      ENDIF
      !
    ELSE
      !
      CALL GET_NEAR_POINTS(JP1,JP2,ZDIS)
      CALL GET_NEAR_POINTS(JP2,JP1,ZDIS)
      !
    ENDIF
  ENDDO
  !
ENDDO
!$OMP END PARALLEL DO
!
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_IGN_2',1,ZHOOK_HANDLE)
!
CONTAINS
!
SUBROUTINE GET_NEAR_POINTS(KP1,KP2,PDIS)
!
INTEGER, INTENT(IN) :: KP1
INTEGER, INTENT(IN) :: KP2
REAL, INTENT(IN) :: PDIS
!
INTEGER :: JN1,JN2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
DO JN1 = 1,KNEAR_NBR
  !
  IF (PDIS<ZNEAR(KP1,JN1) .OR. KNEAR(KP1,JN1)==0) THEN
    !
    IF (JN1<KNEAR_NBR) THEN
      !
      DO JN2=KNEAR_NBR,JN1+1,-1
        !
        ZNEAR(KP1,JN2) = ZNEAR(KP1,JN2-1)
        KNEAR(KP1,JN2) = KNEAR(KP1,JN2-1)
        !
      ENDDO
      !
    ENDIF
    !
    ZNEAR(KP1,JN1) = PDIS
    KNEAR(KP1,JN1) = KP2
    !
    EXIT
    !
  ENDIF
  !
ENDDO
!
END SUBROUTINE GET_NEAR_POINTS
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_NEAR_MESHES_IGN
