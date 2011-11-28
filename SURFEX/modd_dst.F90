!     ################
MODULE MODD_DST

!Purpose: 
!Declare variables and constants necessary to do the dust calculations
!Here are only the variables which depend on the grid!

!Author: Alf Grini <alf.grini@cnrm.meteo.fr>

IMPLICIT NONE

REAL,    PARAMETER  :: XDENSITY_DST=2500.0 ![kg/m3] density of dust
REAL,    PARAMETER  :: XMOLARWEIGHT_DST=100.d-3 ![kg/mol] molar weight of dust
!Set emission related parameters
INTEGER, PARAMETER  :: NEMISMODES_MAX=3
INTEGER, DIMENSION(NEMISMODES_MAX), PARAMETER :: JORDER=(/3,2,1/)        !Dust modes in order of importance

  REAL,DIMENSION(NEMISMODES_MAX)   :: XEMISRADIUS_INI=0.5*(/1.5, 6.7, 14.2/) ! mass madian radius initialization for dust mode (um)
  REAL,DIMENSION(NEMISMODES_MAX)   :: XEMISSIG_INI=(/1.70, 1.60, 1.50/)    ! dispersion initialization for dust mode
  REAL,DIMENSION(NEMISMODES_MAX)   :: XMSS_FRC_SRC_INI=(/0.09, 0.45, 0.46/) ! Mass fraction from each mode


END MODULE MODD_DST
