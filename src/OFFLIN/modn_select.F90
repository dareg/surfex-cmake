!!
!!    #####################
      MODULE MODN_SELECT
!!    #####################
!!
!!*** *MODN_SELECT*
!!
!!    PURPOSE
!!    -------
!       namelist for output writing selection
!!
!!**  AUTHOR
!!    ------
!!    P. Le Moigne      *CNRM*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 3/2009
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_SELECT, ONLY : LSELECT,CNAME_SELECT
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
SAVE
NAMELIST /NAM_SELECT/LSELECT,CNAME_SELECT
!
END MODULE MODN_SELECT
