!     #########
      SUBROUTINE READ_NAM_PGD_ISBA(HPROGRAM, KPATCH, KGROUND_LAYER,                         &
                                   HISBA, HPEDOTF, HPHOTO, OTR_ML,                          &
                                   HCLAY, HCLAYFILETYPE, PUNIF_CLAY, OIMP_CLAY,             &
                                   HSAND, HSANDFILETYPE, PUNIF_SAND, OIMP_SAND,             &
                                   HSOM_TOP, HSOM_SUB, HSOMFILETYPE, PUNIF_SOM, OIMP_SOM,   &
                                   HCTI, HCTIFILETYPE, OIMP_CTI,                            &
                                   HRUNOFFB, HRUNOFFBFILETYPE, PUNIF_RUNOFFB,               &
                                   HWDRAIN,  HWDRAINFILETYPE , PUNIF_WDRAIN, PSOILGRID,     &
                                   HPH, HPHFILETYPE, PUNIF_PH, HFERT, HFERTFILETYPE,        &
                                   PUNIF_FERT      )  
!     ##############################################################
!
!!**** *READ_NAM_PGD_ISBA* reads namelist for ISBA
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
!!    EXTERNAL
!!    --------
!OTR_ML,                        !
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    01/2005
!!       2008 B. Decharme : uniform value of subgrid drainage coefficient
!!    12/2008 E. Martin   : files of data for subgrid drainage 
!!                          and subgridrunoff
!!    06/2009 B. Decharme : files of data for topographic index
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_PAR, ONLY : XUNDEF, NUNDEF
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
!
USE MODE_POS_SURF
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
CHARACTER(LEN=6),    INTENT(IN)    :: HPROGRAM      ! Type of program
INTEGER,             INTENT(OUT)   :: KPATCH        ! number of patches
INTEGER,             INTENT(OUT)   :: KGROUND_LAYER ! number of soil layers
CHARACTER(LEN=3),    INTENT(OUT)   :: HISBA         ! ISBA option
CHARACTER(LEN=4),    INTENT(OUT)   :: HPEDOTF       ! Pedo-transfert function for DIF
CHARACTER(LEN=3),    INTENT(OUT)   :: HPHOTO        ! photosynthesis option
LOGICAL,             INTENT(OUT)   :: OTR_ML        ! new radiative transfert
CHARACTER(LEN=28),   INTENT(OUT)   :: HSAND         ! file name for sand fraction
CHARACTER(LEN=28),   INTENT(OUT)   :: HCLAY         ! file name for clay fraction
CHARACTER(LEN=28),   INTENT(OUT)   :: HCTI          ! file name for topographic index
CHARACTER(LEN=28),   INTENT(OUT)   :: HRUNOFFB      ! file name for runoffb parameter
CHARACTER(LEN=28),   INTENT(OUT)   :: HWDRAIN       ! file name for wdrain parameter
CHARACTER(LEN=6),    INTENT(OUT)   :: HSANDFILETYPE ! sand data file type
CHARACTER(LEN=6),    INTENT(OUT)   :: HCLAYFILETYPE ! clay data file type
CHARACTER(LEN=6),    INTENT(OUT)   :: HCTIFILETYPE  ! topographic index data file type
CHARACTER(LEN=6),    INTENT(OUT)   :: HRUNOFFBFILETYPE ! subgrid runoff data file type
CHARACTER(LEN=6),    INTENT(OUT)   :: HWDRAINFILETYPE  ! subgrid drainage data file type
REAL,                INTENT(OUT)   :: PUNIF_SAND    ! uniform value of sand fraction
REAL,                INTENT(OUT)   :: PUNIF_CLAY    ! uniform value of clay fraction
REAL,                INTENT(OUT)   :: PUNIF_RUNOFFB ! uniform value of subgrid runoff coefficient
REAL,                INTENT(OUT)   :: PUNIF_WDRAIN  ! uniform value of subgrid drainage coefficient
LOGICAL,             INTENT(OUT)   :: OIMP_SAND     ! Imposed values for Sand
LOGICAL,             INTENT(OUT)   :: OIMP_CLAY     ! Imposed values for Clay
LOGICAL,             INTENT(OUT)   :: OIMP_CTI      ! Imposed values for topographic index statistics
CHARACTER(LEN=28),   INTENT(OUT)   :: HSOM_TOP      ! file name for organic matter
CHARACTER(LEN=28),   INTENT(OUT)   :: HSOM_SUB      ! file name for organic matter
CHARACTER(LEN=6),    INTENT(OUT)   :: HSOMFILETYPE  ! organic matter data file type
REAL,                INTENT(OUT)   :: PUNIF_SOM     ! uniform value of organic matter (%)
LOGICAL,             INTENT(OUT)   :: OIMP_SOM      ! Imposed maps of organic matter
REAL, DIMENSION(:),  INTENT(OUT)   :: PSOILGRID     ! Soil layer thickness for DIF
CHARACTER(LEN=28),   INTENT(OUT)   :: HPH           ! file name for pH
CHARACTER(LEN=28),   INTENT(OUT)   :: HFERT         ! file name for fertilisation rate
CHARACTER(LEN=6),    INTENT(OUT)   :: HPHFILETYPE   ! pH data file type
CHARACTER(LEN=6),    INTENT(OUT)   :: HFERTFILETYPE ! fertilisation data file type
REAL,                INTENT(OUT)   :: PUNIF_PH      ! uniform value of pH
REAL,                INTENT(OUT)   :: PUNIF_FERT    ! uniform value of fertilisation rate
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER                           :: ILUOUT    ! output listing logical unit
INTEGER                           :: ILUNAM    ! namelist file logical unit
LOGICAL                           :: GFOUND    ! flag when namelist is present
!
!*    0.3    Declaration of namelists
!            ------------------------
!
INTEGER                  :: NPATCH           ! number of patches
INTEGER                  :: NGROUND_LAYER    ! number of soil layers
CHARACTER(LEN=3)         :: CISBA            ! ISBA option
CHARACTER(LEN=4)         :: CPEDO_FUNCTION   ! Pedo-transfert function for DIF
CHARACTER(LEN=3)         :: CPHOTO           ! photosynthesis option
LOGICAL                  :: LTR_ML           ! new radiative transfert
CHARACTER(LEN=28)        :: YSAND            ! file name for sand fraction
CHARACTER(LEN=28)        :: YCLAY            ! file name for clay fraction
CHARACTER(LEN=28)        :: YCTI             ! file name for topographic index
CHARACTER(LEN=28)        :: YRUNOFFB         ! file name for runoffb parameter
CHARACTER(LEN=28)        :: YWDRAIN          ! file name for wdrain parameter
CHARACTER(LEN=28)        :: YPH              ! file name for pH
CHARACTER(LEN=28)        :: YFERT            ! file name for fertilisation rate
CHARACTER(LEN=6)         :: YSANDFILETYPE    ! sand data file type
CHARACTER(LEN=6)         :: YCLAYFILETYPE    ! clay data file type
CHARACTER(LEN=6)         :: YCTIFILETYPE     ! topographic index data file type
CHARACTER(LEN=6)         :: YRUNOFFBFILETYPE ! subgrid runoff data file type
CHARACTER(LEN=6)         :: YWDRAINFILETYPE  ! subgrid drainage data file type
CHARACTER(LEN=6)         :: YPHFILETYPE      ! pH data file type
CHARACTER(LEN=6)         :: YFERTFILETYPE    ! fertilisation data file type
LOGICAL                  :: LIMP_SAND        ! Imposed maps of Sand from another PGD file
LOGICAL                  :: LIMP_CLAY        ! Imposed maps of Clay from another PGD file
LOGICAL                  :: LIMP_CTI         ! Imposed values for topographic index statistics from another PGD file
REAL                     :: XUNIF_SAND    ! uniform value of sand fraction
REAL                     :: XUNIF_CLAY    ! uniform value of clay fraction
REAL                     :: XUNIF_RUNOFFB ! uniform value of subgrid runoff coefficient
REAL                     :: XUNIF_WDRAIN  ! uniform value of subgrid drainage coefficient
REAL                     :: XUNIF_PH      ! uniform value of pH
REAL                     :: XUNIF_FERT    ! uniform value of fertilisation rate
!
REAL, DIMENSION(150)     :: XSOILGRID     ! Soil layer thickness for DIF
!
CHARACTER(LEN=28)        :: YSOM_TOP      ! file name for organic matter
CHARACTER(LEN=28)        :: YSOM_SUB      ! file name for organic matter
CHARACTER(LEN=6)         :: YSOMFILETYPE  ! organic matter data file type
REAL                     :: XUNIF_SOM     ! uniform value of organic matter (%)
LOGICAL                  :: LIMP_SOM      ! Imposed maps of organic matter
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
NAMELIST/NAM_ISBA/ NPATCH, NGROUND_LAYER, CISBA, CPEDO_FUNCTION, CPHOTO,   &
                   LTR_ML, YCLAY, YCLAYFILETYPE, XUNIF_CLAY, LIMP_CLAY,    &
                   YSAND, YSANDFILETYPE, XUNIF_SAND, LIMP_SAND,            &
                   YSOM_TOP, YSOM_SUB, YSOMFILETYPE, XUNIF_SOM, LIMP_SOM,  &
                   YCTI, YCTIFILETYPE, LIMP_CTI,                           &
                   YRUNOFFB, YRUNOFFBFILETYPE, XUNIF_RUNOFFB,              &
                   YWDRAIN,  YWDRAINFILETYPE,  XUNIF_WDRAIN, XSOILGRID,    &
                   YPH, YPHFILETYPE, XUNIF_PH, YFERT, YFERTFILETYPE,       &
                   XUNIF_FERT   
!
!-------------------------------------------------------------------------------
!
!*    1.      Initializations of defaults
!             ---------------------------
!
!#####################
IF (LHOOK) CALL DR_HOOK('READ_NAM_PGD_ISBA',0,ZHOOK_HANDLE)
NPATCH         = 1
NGROUND_LAYER  = NUNDEF
CISBA          = '3-L'
CPEDO_FUNCTION = 'CH78'
CPHOTO         = 'NON'
LTR_ML         = .FALSE.
XSOILGRID(:)   = XUNDEF
!#####################
!
XUNIF_CLAY       = 0.33
XUNIF_SAND       = 0.33
XUNIF_SOM        = XUNDEF
XUNIF_RUNOFFB    = 0.5
XUNIF_WDRAIN     = 0.
XUNIF_PH        = XUNDEF
XUNIF_FERT      = XUNDEF
!
YCLAY            = '                          '
YSAND            = '                          '
YSOM_TOP         = '                          '
YSOM_SUB         = '                          '
YCTI             = '                          '
YRUNOFFB         = '                          '
YWDRAIN          = '                          '
YPH              = '                          '
YFERT            = '                          '
!
YCLAYFILETYPE    = '      '
YSANDFILETYPE    = '      '
YSOMFILETYPE     = '      '
YCTIFILETYPE     = '      '
YRUNOFFBFILETYPE = '      '
YWDRAINFILETYPE  = '      ' 
YPHFILETYPE      = '      '
YPHFILETYPE      = '      '
!
LIMP_CLAY        = .FALSE.
LIMP_SAND        = .FALSE.
LIMP_SOM         = .FALSE.
LIMP_CTI         = .FALSE.
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*    2.      Reading of namelist
!             -------------------
!
CALL OPEN_NAMELIST(HPROGRAM,ILUNAM)
!
CALL POSNAM(ILUNAM,'NAM_ISBA',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_ISBA)
!
CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
!
!-------------------------------------------------------------------------------
!
KPATCH           = NPATCH           ! number of patches
KGROUND_LAYER    = NGROUND_LAYER    ! number of soil layers
PSOILGRID        = XSOILGRID        ! soil layer tickness for DIF
HISBA            = CISBA            ! ISBA option
HPEDOTF          = CPEDO_FUNCTION   ! Pedo-transfert function for DIF
HPHOTO           = CPHOTO           ! photosynthesis option
OTR_ML           = LTR_ML           ! new radiative transfert
HSAND            = YSAND            ! file name for sand fraction
HCLAY            = YCLAY            ! file name for clay fraction
HSOM_TOP         = YSOM_TOP         ! file name for organic matter
HSOM_SUB         = YSOM_SUB         ! file name for organic matter
HCTI             = YCTI             ! file name for topographic index
HRUNOFFB         = YRUNOFFB         ! file name for subgrid runoff
HWDRAIN          = YWDRAIN          ! file name for subgrid drainage
HSANDFILETYPE    = YSANDFILETYPE    ! sand data file type
HCLAYFILETYPE    = YCLAYFILETYPE    ! clay data file type
HSOMFILETYPE     = YSOMFILETYPE     ! organic matter data file type
HCTIFILETYPE     = YCTIFILETYPE     ! topographic index data file type
HRUNOFFBFILETYPE = YRUNOFFBFILETYPE ! subgrid runoff data file type
HWDRAINFILETYPE  = YWDRAINFILETYPE  ! subgrid drainage data file type
PUNIF_SAND       = XUNIF_SAND       ! uniform value of sand fraction
PUNIF_CLAY       = XUNIF_CLAY       ! uniform value of clay fraction
PUNIF_SOM        = XUNIF_SOM        ! uniform value of organic matter
PUNIF_RUNOFFB    = XUNIF_RUNOFFB    ! uniform value of subgrid runoff coefficient
PUNIF_WDRAIN     = XUNIF_WDRAIN     ! uniform value of subgrid drainage coefficient
OIMP_SAND        = LIMP_SAND        ! Imposed values for SAND
OIMP_CLAY        = LIMP_CLAY        ! Imposed values for CLAY
OIMP_SOM         = LIMP_SOM         ! Imposed values for organic matter
OIMP_CTI         = LIMP_CTI         ! Imposed values for topographic index statistics
!
HPH           = YPH           ! file name for pH value
HFERT         = YFERT         ! file name for fertilisation data
HPHFILETYPE   = YPHFILETYPE   ! pH data file type
HFERTFILETYPE = YFERTFILETYPE ! Fertilisation data file type
PUNIF_PH      = XUNIF_PH      ! uniform value of pH
PUNIF_FERT    = XUNIF_FERT    ! uniform value of fertilisation rate
!
IF (LHOOK) CALL DR_HOOK('READ_NAM_PGD_ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_NAM_PGD_ISBA
