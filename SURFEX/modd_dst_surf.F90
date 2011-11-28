!     ################
MODULE MODD_DST_SURF

  USE MODD_CSTS, only :    &
         XPI                 &!Definition of pi
         ,XMD                &![kg/mol] molar weight of air
         ,XAVOGADRO           ![molec/mol] avogadros number  

  IMPLICIT NONE

  real, parameter  :: xdensity_dust = 2.5e3         ![kg/m3] density of dust
  real, parameter  :: xmolarweight_dust = 100.e-3   ![kg/mol] molar weight dust
  real, parameter  :: xm3toum3          = 1.d18     ![um3/m3] conversion factor
  real, parameter  :: xum3tom3          = 1.d-18    ![m3/um3] conversion factor
  real, parameter  :: xsixth            = 1./6.     ![-] one sixth
  INTEGER      :: JPMODE_DST= 3  ! number of dust modes (max 3; default = 1)
  REAL, PARAMETER                        :: XDENSITY_DST=2500.0 ![kg/m3] density of dust
  REAL, PARAMETER                        :: XMOLARWEIGHT_DST=100.d-3 ![kg/mol] molar weight of dust
  !
  INTEGER                                :: NSIZE_LARGEST_DST !Size of largest dust emitter vector in any patch
  INTEGER                                :: NVEGNO_DST        !Number of vegetation classes considered dust emitters
  CHARACTER(LEN=6), DIMENSION(:), ALLOCATABLE         :: CSV_DST           !Name of scalar variables 
  INTEGER                                :: NDST_MDEBEG        !Index of mass flux in first dust mode in scalar list
  INTEGER                                :: NDSTMDE           !Number of dust modes emitted
  CHARACTER(LEN=1)                       :: CDSTYN            !Do we calculate dust ('Y' or 'N') depending CSV list
  CHARACTER(LEN=6), DIMENSION(:), ALLOCATABLE         :: CDSTNAMES         !Names of dust related scalars
    !Set emission related parameters
  INTEGER, PARAMETER  :: NEMISMODES_MAX=3
  INTEGER, DIMENSION(NEMISMODES_MAX), PARAMETER :: JORDER=(/3,2,1/) !Dust modes in order of importance
  CHARACTER(LEN=5)    :: CEMISPARAM='AMMA'        ! Reference to paper where emission parameterization is proposed
  CHARACTER(LEN=6)    :: CVERMOD='NONE'            ! CVERMOD='CMDVER' Reference pour activer la version modifieé
  CHARACTER(LEN=4)  :: CRGUNITD   = 'NUMB'             ! type of log-normal geometric mean radius
   LOGICAL      :: LVARSIG    = .FALSE.   ! switch to active pronostic dispersion for all modes
   LOGICAL      :: LRGFIX_DST = .TRUE.   ! switch to active pronostic mean radius for all modes
  ! XEMISRADIUS_INI, XEMISSIG_INI, XMSS_FRC_SRC_INI are defined in INIT_DST_n
  REAL,DIMENSION(NEMISMODES_MAX)   :: XEMISRADIUS_INI  ! mass madian radius initialization for dust mode (um)
  REAL,DIMENSION(NEMISMODES_MAX)   :: XEMISSIG_INI     ! dispersion initialization for dust mode
  REAL,DIMENSION(NEMISMODES_MAX)   :: XMSS_FRC_SRC_INI ! Mass fraction from each mode
  REAL :: XFLX_MSS_FDG_FCT = 8.E-4


END MODULE MODD_DST_SURF
