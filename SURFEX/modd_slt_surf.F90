MODULE MODD_SLT_SURF

  USE MODD_CSTS, only :    &
         XPI                 &!Definition of pi
         ,XMD                &![kg/mol] molar weight of air
         ,XAVOGADRO           ![molec/mol] avogadros number  

  IMPLICIT NONE

  real, parameter  :: xdensity_salt = 2.1e3         ![kg/m3] density of sea salt
  real, parameter  :: xmolarweight_salt = 58.e-3   ![kg/mol] molar weight sea salt
  real, parameter  :: xm3toum3          = 1.d18     ![um3/m3] conversion factor
  real, parameter  :: xum3tom3          = 1.d-18    ![m3/um3] conversion factor
  real, parameter  :: xsixth            = 1./6.     ![-] one sixth
  INTEGER      :: JPMODE_SLT= 3  ! number of sea salt modes (max 3; default = 1)
  REAL, PARAMETER                        :: XDENSITY_SLT=2100.0 ![kg/m3] density of sea salt
  REAL, PARAMETER                        :: XMOLARWEIGHT_SLT=58.d-3 ![kg/mol] molar weight of sea salt
  !
  INTEGER                                :: NSIZE_LARGEST_SLT !Size of largest sea salt emitter vector in any patch
  INTEGER                                :: NVEGNO_SLT        !Number of vegetation classes considered sea salt emitters
  CHARACTER(LEN=6), DIMENSION(:), ALLOCATABLE         :: CSV_SLT           !Name of scalar variables 
  INTEGER                                :: NSLT_MDEBEG        !Index of mass flux in first sea salt mode in scalar list
  INTEGER                                :: NSLTMDE           !Number of sea salt modes emitted
  CHARACTER(LEN=1)                       :: CSLTYN            !Do we calculate sea salt ('Y' or 'N') depending CSV list
  CHARACTER(LEN=6), DIMENSION(:), ALLOCATABLE         :: CSLTNAMES         !Names of sea salt related scalars
    !Set emission related parameters
  INTEGER, PARAMETER  :: NEMISMODES_MAX=3
  CHARACTER(LEN=4)  :: CRGUNITS   = 'MASS'  ! type of log-normal geometric mean radius
  INTEGER, DIMENSION(NEMISMODES_MAX), PARAMETER :: JORDER_SLT=(/3,2,1/) !Dust modes in order of importance
  CHARACTER(LEN=5)    :: CEMISPARAM_SLT='Sig01'        ! Reference to paper where emission parameterization is proposed
  REAL,DIMENSION(NEMISMODES_MAX)   :: XEMISRADIUS_INI_SLT=(/0.2, 2.0, 12.0/)         ! number madian radius initialization for sea salt mode (um)
  REAL,DIMENSION(NEMISMODES_MAX)   :: XEMISSIG_INI_SLT=(/1.9, 2.0, 3.0/)             ! dispersion initialization for sea salt mode
   LOGICAL      :: LVARSIG_SLT    = .FALSE.   ! switch to active pronostic dispersion for all modes
   LOGICAL      :: LRGFIX_SLT     = .TRUE.   ! switch to active pronostic mean radius for all modes


END MODULE MODD_SLT_SURF
