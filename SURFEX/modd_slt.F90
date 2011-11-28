MODULE MODD_SLT

!Purpose: 
!Declare variables and constants necessary to do the sea salt calculations
!Here are only the variables which depend on the grid!

!Author: Alf Grini / Pierre Tulet

IMPLICIT NONE

REAL,    PARAMETER  :: XDENSITY_SLT=2100.0 ![kg/m3] density of sea salt
REAL,    PARAMETER  :: XMOLARWEIGHT_SLT=58.d-3 ![kg/mol] molar weight of sea salt
!Set emission related parameters
INTEGER, PARAMETER  :: NEMISMODES_MAX=3
INTEGER, DIMENSION(NEMISMODES_MAX), PARAMETER :: JORDER=(/2,1,3/)        !Sea salt modes in order of importance

  !REAL,DIMENSION(NEMISMODES_MAX)   :: XEMISRADIUS_INI=(/0.2, 2.0, 12.0/) ! number madian radius initialization for sea salt mode (um) values from Vignatti for a distribution in number
  REAL,DIMENSION(NEMISMODES_MAX)   :: XEMISRADIUS_INI=0.5*(/0.28, 2.25, 15.32/) ! number madian radius initialization for sea salt mode (um)
  !REAL,DIMENSION(NEMISMODES_MAX)   :: XEMISSIG_INI=(/1.9, 2.0, 3.0/)     ! dispersion initialization for sea salt mode; values from Vignatti
  REAL,DIMENSION(NEMISMODES_MAX)   :: XEMISSIG_INI=(/1.58, 2.0, 2.0/)     ! dispersion initialization for sea salt mode

END MODULE MODD_SLT
