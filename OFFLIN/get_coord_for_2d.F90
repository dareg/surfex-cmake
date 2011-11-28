!     #########
      SUBROUTINE GET_COORD_FOR_2D(KLUOUT,IL,ZXX,ZYY)
!     #######################################################
!!****  *GET_FOR_2D* - 
!!
!!    PURPOSE
!!    -------
!!
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
!!	S. Faroux   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2010 
!-------------------------------------------------------------------------------
!
USE MODD_SURF_ATM_GRID_n, ONLY : CGRID
!
USE MODI_GET_GRID_COORD
!
INTEGER, INTENT(IN)  :: KLUOUT   ! output listing
INTEGER, INTENT(IN) :: IL
REAL, DIMENSION(:), INTENT(OUT) :: ZXX
REAL, DIMENSION(:), INTENT(OUT) :: ZYY
!
REAL, DIMENSION(IL)   :: ZX           
REAL, DIMENSION(IL)   :: ZY        
INTEGER :: J
!
!
IF (CGRID.EQ.'CONF PROJ '.OR.CGRID.EQ.'CARTESIAN '&
        .OR.CGRID.EQ.'LONLAT REG') THEN
        CALL GET_GRID_COORD(KLUOUT,PX=ZX,PY=ZY)
        DO J=1,SIZE(ZXX)
          ZXX(J)=ZX(J)
        ENDDO
        DO J=1,SIZE(ZYY)
          ZYY(J)=ZY((J-1)*SIZE(ZXX)+1)
        ENDDO
ENDIF
!
!
END SUBROUTINE GET_COORD_FOR_2D
