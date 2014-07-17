!     ##################
      MODULE MODD_ASSIM
!     ##################
!
!!****  *MODD_ASSIM - declaration of keys for assimilation schemes
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	L. Jarlan   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       23/02/05
! 
!       Inclusion of OI constants 21/05/09 (J.-F. Mahfouf)  
!!       Add all assim keys         04/2012  T.Aspelien
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE

!-------------------------------------------------------------------------------
!
! Assimilation Scheme Options:
!
 LOGICAL                               :: LASSIM               ! Assimilation or not
                                                               !'.TRUE.'
                                                               !'.FALSE.'
 LOGICAL                               :: LAROME               ! If reading AROME fields
 LOGICAL                               :: LECSST               ! Use ECMWF SST
 LOGICAL                               :: LAESST               ! SST analysis performed
 LOGICAL                               :: LAESNM               ! Update snow analysis
 LOGICAL                               :: LALADSURF            
 LOGICAL                               :: LREAD_SST_FROM_FILE  ! Read SST from file
 LOGICAL                               :: LWATERTG2            ! Use deep soil temperature as lake temp.
 LOGICAL                               :: LEXTRAP_SEA          ! Extrapolation of sea points
 LOGICAL                               :: LEXTRAP_WATER        ! Extrapolation of inland water  points
 LOGICAL                               :: LEXTRAP_NATURE       ! Extrapolation of nature points
 LOGICAL                               :: LPRT                 ! Running VARASSIM in a perturbation mode
 LOGICAL                               :: LSIM                 ! Running VARASSIM in a perturbation mode 
 LOGICAL                               :: LBEV                 ! Running VARASSIM to evolve B
 LOGICAL                               :: LBFIXED    
 LOGICAL                               :: LOBSFILE 

 INTEGER, PARAMETER                    :: NOBSMAX = 4
 INTEGER, PARAMETER                    :: NVARMAX = 5
 INTEGER,DIMENSION(NOBSMAX)            :: NNCO                 ! Select the type of observations to be assimilated 
 INTEGER,DIMENSION(NVARMAX)            :: NNCV                 ! Select the control variables to be used 
 INTEGER                               :: NOBSTYPE
 INTEGER                               :: NOBS
 INTEGER                               :: NIPERT 
 INTEGER                               :: NIVAR                ! counter for ctnrl vars
 INTEGER                               :: NVAR                 ! number of cntrl vars
 INTEGER                               :: NBOUTPUT  
 INTEGER                               :: NPRINTLEV            ! Verbosity 

 CHARACTER(LEN=12)                     :: CBIO                 ! Name of Biomass variable
 CHARACTER(LEN=100)                    :: CPREFIX_BIO          ! The prefix of the Biomass variable 
 CHARACTER(LEN=5)                      :: CASSIM_ISBA          ! OI/EKF
 CHARACTER(LEN=5)                      :: CASSIM               ! type of correction
 CHARACTER(LEN=3),DIMENSION(NVARMAX)   :: CVAR_M               ! X is ctrl
                                                               ! 'PLUS ' (default)
                                                               ! 'AVERA'            
                                                               ! '2DVAR'
 CHARACTER(LEN=100),DIMENSION(NVARMAX) :: CPREFIX_M            ! The prefix of the control variables (in PREP.txt file) (max dim)                          
 CHARACTER(LEN=10),DIMENSION(:), ALLOCATABLE  :: COBS          ! Identifier for simulated observations
 CHARACTER(LEN=3),DIMENSION(:), ALLOCATABLE   :: CVAR          ! Identifier for control variable

REAL,DIMENSION(12)                     :: XALPH
 REAL,DIMENSION(NVARMAX)               :: XTPRT_M              ! The perturbation amplitude (max dim)
 REAL,DIMENSION(NVARMAX)               :: XSIGMA_M             ! covariance of background errors if B is fixed (max dim)
!                                                              ! covariance of model errors if B evolving (max dim)
 REAL,DIMENSION(NOBSMAX)               :: XERROBS_M            ! Observational standard deviation
 REAL,DIMENSION(:,:,:,:),ALLOCATABLE   :: XF_PATCH             ! vector of model observations (for each pacth)
 REAL,DIMENSION(:,:,:,:),ALLOCATABLE   :: XF                   ! Vector of forecast control variables 
 REAL,DIMENSION(:,:,:),ALLOCATABLE     :: XI 
 REAL,DIMENSION(:,:), ALLOCATABLE      :: XYO                  ! vector of observations
 REAL,DIMENSION(:,:),ALLOCATABLE       :: XLAI_PASS
 REAL,DIMENSION(:,:),ALLOCATABLE       :: XBIO_PASS
 REAL,DIMENSION(:,:),ALLOCATABLE       :: XAT2M_ISBA
 REAL,DIMENSION(:,:),ALLOCATABLE       :: XAHU2M_ISBA
 REAL,DIMENSION(:,:),ALLOCATABLE       :: XAZON10M_ISBA
 REAL,DIMENSION(:,:),ALLOCATABLE       :: XAMER10M_ISBA
 REAL,DIMENSION(:),ALLOCATABLE         :: XAT2M_TEB
 REAL,DIMENSION(:),ALLOCATABLE         :: XTPRT           ! The perturbation amplitude
 REAL,DIMENSION(:),ALLOCATABLE         :: XSIGMA          ! covariance of background errors if B is fixed
                                                          ! covariance of model errors if B evolving  
 REAL,DIMENSION(:),ALLOCATABLE         :: XERROBS
 REAL                                  :: XSCALE_Q        ! scaling factor of Q matrix w.r.t. the initial B
 REAL                                  :: XSCALE_QLAI

!
! Constants and options of the soil OI analysis
!
 LOGICAL ::  LHUMID,  LIMVEG, LISSEW,   L_SM_WP, LFGEL,      LCLIM,   LLDHMT,  &
             LOBSWG,  LOBS2M 
 INTEGER ::  NMINDJ,   NNEBUL, NNEIGT,   NNEIGW,  NR_SM_WP,  NECHGU,  NTVGLA,  &
             NSEAICE, NLISSEW, NIDJ,    NITRAD 
 REAL    ::  XANEBUL,  XRCLIMN, XRCLIMTP,  XRCLIMTS, XRCLIMV, XRCLIMWP, XRCLIMWS, &
             XSCOEFH,  XSCOEFT, XSEVAP,    XSIGH2MO, XSIGT2MO,    XSNEIGT,  XSNEIGW,  &
             XSPRECIP, XSWFC,   XV10MX,    XRD1,     XRTINER,     XWCRIN,   XWPMX,    &
             XWSMX,    XTMERGL, XRZHZ0G,   XRCLIMCA, XRCLISST,    XRWPIA,   XRWPIB,   &
             XRSNSA,   XRSNSB,  XSALBM,    XSALBB,   XSEMIB,      XSZZ0B,   XSMU0,    &
             XSICE,    XSEMIM,  XRA_SM_WP, XRSCALDW, XSPRECIP2,                     &
             XREPSM,   XRCDTR,  XSIGHP1,   XSIGT2MR, XSIGH2MR,    XRSABR,            &
             XRARGR,   XGWFC,   XEWFC,     XGWWILT,  XEWWILT,     XG1WSAT,  XG2WSAT,  &
             XREPS1,   XREPS2,  XREPS3,    XADWR,    XSODELX(0:9),                  &
             XSIGWGO,  XSIGWGB, XSIGW2B,   XRTHR_QC, XSIGWGO_MAX, XRSCAL_JAC
!
END MODULE MODD_ASSIM
