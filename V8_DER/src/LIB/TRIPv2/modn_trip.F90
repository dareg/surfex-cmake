!#######################
MODULE  MODN_TRIP
!#######################
!
!!****  *MODN_TRIP* define the variables and namelist for TRIP
!                       driver programs
!!
!!    PURPOSE
!!    -------
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
CHARACTER(LEN=3) :: CGROUNDW = 'DEF'  !Use groundwater scheme
                                      !'DEF' = No groundwater scheme
                                      !'CST' = Constant transfert time
                                      !'DIF' = Groundwater diffusive scheme 
CHARACTER(LEN=3) :: CVIT  = 'DEF'     !Type of stream flow velocity
                                      !'DEF' = Constant velocit = 0.5m/s
                                      !'VAR' = variable velocity
LOGICAL          :: LFLOOD = .FALSE.  !if true, use TRIP-FLOOD
!
REAL             :: XTAUG_UNIF = 30.0 ! Constant transfert time value
REAL             :: XCVEL      = 0.5  ! Constant velocity value
REAL             :: XRATMED    = 1.4  ! Meandering ratio
REAL             :: XTSTEP     = 3600.
!
!-------------------------------------------------------------------------------
!
!*       1.    NAMELISTS
!              ---------
!
NAMELIST/NAM_TRIP/CGROUNDW,CVIT,LFLOOD,XTAUG_UNIF,XCVEL,XRATMED,XTSTEP
!
!-------------------------------------------------------------------------------
END MODULE MODN_TRIP
