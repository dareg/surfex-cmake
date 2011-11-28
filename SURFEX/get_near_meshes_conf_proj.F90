!     #########
      SUBROUTINE GET_NEAR_MESHES_CONF_PROJ(KGRID_PAR,KL,PGRID_PAR,KNEAR_NBR,OLIST,KLIST,KNEAR)
!     ##############################################################
!
!!**** *GET_NEAR_MESHES_CONF_PROJ* get the near grid mesh indices
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
USE MODE_GRIDTYPE_CONF_PROJ
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
INTEGER, DIMENSION(KLIST,KNEAR_NBR),INTENT(OUT)   :: KNEAR  ! near mesh indices
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
INTEGER                            :: IIMAX, IJMAX
INTEGER                            :: JI, JJ
INTEGER                            :: JX, JY
INTEGER                            :: JL
INTEGER                            :: JLIST
INTEGER                            :: IDIST
INTEGER                            :: ICOUNT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_CONF_PROJ',0,ZHOOK_HANDLE)
CALL GET_GRIDTYPE_CONF_PROJ(PGRID_PAR,KIMAX=IIMAX,KJMAX=IJMAX)
!
KNEAR  (:,:) = 0
!
IDIST = INT(SQRT(FLOAT(KNEAR_NBR)))
!
JLIST = 0
!
IF (IIMAX*IJMAX==KL) THEN
  DO JJ=1,IJMAX
    DO JI=1,IIMAX
      ICOUNT = 0
      JL = JI + IIMAX * (JJ-1)
      IF (.NOT. OLIST(JL)) CYCLE
      JLIST = JLIST + 1
      KNEAR(JLIST,:) = 0      
      DO JX=-(IDIST-1)/2,IDIST/2
        DO JY=-(IDIST-1)/2,IDIST/2
          IF (JI+JX>0 .AND. JI+JX<IIMAX+1 .AND. JJ+JY>0 .AND. JJ+JY<IJMAX+1) THEN
            ICOUNT = ICOUNT + 1
            KNEAR(JLIST,ICOUNT) = (JI+JX) + IIMAX * (JJ+JY-1)
          END IF
        END DO
      END DO
    END DO
  END DO
END IF
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_CONF_PROJ',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_NEAR_MESHES_CONF_PROJ
