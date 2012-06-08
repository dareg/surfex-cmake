!     #########
      SUBROUTINE GET_NEAR_MESHES_IGN(KGRID_PAR,KL,PGRID_PAR,KNEAR_NBR,OLIST,KLIST,KNEAR)
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
INTEGER, DIMENSION(KLIST,KNEAR_NBR),INTENT(OUT) :: KNEAR    ! near mesh indices
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_IGN',0,ZHOOK_HANDLE)
KNEAR  (:,:) = 0
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_IGN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_NEAR_MESHES_IGN
