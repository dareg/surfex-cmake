!     #########
      SUBROUTINE GET_DIM_FOR_2D(KDIM1,KDIM2,HTYPE)
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
USE MODD_SURF_ATM_GRID_n, ONLY : CGRID, XGRID_PAR, NGRID_PAR
!
USE MODI_GET_GRID_DIM
!
INTEGER, INTENT(OUT) :: KDIM1     ! 1st dimension
INTEGER, INTENT(OUT) :: KDIM2     ! 2nd dimension
CHARACTER(LEN=3), INTENT(OUT) :: HTYPE
!
LOGICAL   :: ORECT     ! T if rectangular grid     
INTEGER :: J
!
KDIM1=0
KDIM2=0
HTYPE='XY '
IF (CGRID.EQ.'CONF PROJ '.OR.CGRID.EQ.'CARTESIAN '&
        .OR.CGRID.EQ.'LONLAT REG') THEN
        IF (CGRID.EQ.'LONLAT REG') HTYPE='LON'
        CALL GET_GRID_DIM(CGRID,NGRID_PAR,XGRID_PAR,ORECT,KDIM1,KDIM2)
ENDIF
!
END SUBROUTINE GET_DIM_FOR_2D
