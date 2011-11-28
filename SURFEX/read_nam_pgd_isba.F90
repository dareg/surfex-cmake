!     #########
      SUBROUTINE READ_NAM_PGD_ISBA(HPROGRAM,KPATCH, KGROUND_LAYER, HISBA, HPEDOTF,HPHOTO,   &
                                   HCLAY, HCLAYFILETYPE, PUNIF_CLAY, OIMP_CLAY,             &
                                   HSAND, HSANDFILETYPE, PUNIF_SAND, OIMP_SAND,             &
                                   HORGMAT, HORGMATFILETYPE, PUNIF_ORGMAT, OIMP_ORGMAT,     &
                                   HDENSITY, HDENSITYFILETYPE, PUNIF_DENSITY, OIMP_DENSITY, &         
                                   HCTI, HCTIFILETYPE, OIMP_CTI,                            &
                                   HRUNOFFB, HRUNOFFBFILETYPE, PUNIF_RUNOFFB,               &
                                   HWDRAIN,  HWDRAINFILETYPE , PUNIF_WDRAIN                 )  
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
!!
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
CHARACTER(LEN=28),   INTENT(OUT)   :: HORGMAT          ! file name for organic matter
CHARACTER(LEN=28),   INTENT(OUT)   :: HDENSITY         ! file name for soil density
CHARACTER(LEN=6),    INTENT(OUT)   :: HORGMATFILETYPE  ! organic matter data file type
CHARACTER(LEN=6),    INTENT(OUT)   :: HDENSITYFILETYPE ! soil density data file type
REAL,                INTENT(OUT)   :: PUNIF_ORGMAT     ! uniform value of organic matter (%)
REAL,                INTENT(OUT)   :: PUNIF_DENSITY    ! uniform value of soil density
LOGICAL,             INTENT(OUT)   :: OIMP_ORGMAT      ! Imposed maps of organic matter
LOGICAL,             INTENT(OUT)   :: OIMP_DENSITY     ! Imposed maps of soil density
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
CHARACTER(LEN=28)        :: YSAND            ! file name for sand fraction
CHARACTER(LEN=28)        :: YCLAY            ! file name for clay fraction
CHARACTER(LEN=28)        :: YCTI             ! file name for topographic index
CHARACTER(LEN=28)        :: YRUNOFFB         ! file name for runoffb parameter
CHARACTER(LEN=28)        :: YWDRAIN          ! file name for wdrain parameter
CHARACTER(LEN=6)         :: YSANDFILETYPE    ! sand data file type
CHARACTER(LEN=6)         :: YCLAYFILETYPE    ! clay data file type
CHARACTER(LEN=6)         :: YCTIFILETYPE     ! topographic index data file type
CHARACTER(LEN=6)         :: YRUNOFFBFILETYPE ! subgrid runoff data file type
CHARACTER(LEN=6)         :: YWDRAINFILETYPE  ! subgrid drainage data file type
LOGICAL                  :: LIMP_SAND        ! Imposed maps of Sand from another PGD file
LOGICAL                  :: LIMP_CLAY        ! Imposed maps of Clay from another PGD file
LOGICAL                  :: LIMP_CTI         ! Imposed values for topographic index statistics from another PGD file
REAL                     :: XUNIF_SAND    ! uniform value of sand fraction
REAL                     :: XUNIF_CLAY    ! uniform value of clay fraction
REAL                     :: XUNIF_RUNOFFB ! uniform value of subgrid runoff coefficient
REAL                     :: XUNIF_WDRAIN  ! uniform value of subgrid drainage coefficient
!
CHARACTER(LEN=28)        :: YORGMAT          ! file name for organic matter
CHARACTER(LEN=28)        :: YDENSITY         ! file name for soil density
CHARACTER(LEN=6)         :: YORGMATFILETYPE  ! organic matter data file type
CHARACTER(LEN=6)         :: YDENSITYFILETYPE ! soil density data file type
REAL                     :: XUNIF_ORGMAT     ! uniform value of organic matter (%)
REAL                     :: XUNIF_DENSITY    ! uniform value of soil density
LOGICAL                  :: LIMP_ORGMAT      ! Imposed maps of organic matter
LOGICAL                  :: LIMP_DENSITY     ! Imposed maps of soil density
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
NAMELIST/NAM_ISBA/ NPATCH, NGROUND_LAYER, CISBA, CPEDO_FUNCTION, CPHOTO,   &
                   YCLAY, YCLAYFILETYPE, XUNIF_CLAY, LIMP_CLAY,            &
                   YSAND, YSANDFILETYPE, XUNIF_SAND, LIMP_SAND,            &
                   YORGMAT, YORGMATFILETYPE, XUNIF_ORGMAT, LIMP_ORGMAT,    &
                   YDENSITY, YDENSITYFILETYPE, XUNIF_DENSITY, LIMP_DENSITY,&                     
                   YCTI, YCTIFILETYPE, LIMP_CTI,                           &
                   YRUNOFFB, YRUNOFFBFILETYPE, XUNIF_RUNOFFB,              &
                   YWDRAIN,  YWDRAINFILETYPE,  XUNIF_WDRAIN   
!
!-------------------------------------------------------------------------------
!
!*    1.      Initializations of defaults
!             ---------------------------
!
!#####################
IF (LHOOK) CALL DR_HOOK('READ_NAM_PGD_ISBA',0,ZHOOK_HANDLE)
NPATCH         = 1
NGROUND_LAYER  = 3
CISBA          = '3-L'
CPEDO_FUNCTION = 'CH78'
CPHOTO         = 'NON'
!#####################
!
XUNIF_CLAY       = 0.33
XUNIF_SAND       = 0.33
XUNIF_ORGMAT     = XUNDEF
XUNIF_DENSITY    = XUNDEF
XUNIF_RUNOFFB    = 0.5
XUNIF_WDRAIN     = 0.
!
YCLAY            = '                          '
YSAND            = '                          '
YORGMAT          = '                          '
YDENSITY         = '                          '
YCTI             = '                          '
YRUNOFFB         = '                          '
YWDRAIN          = '                          '
!
YCLAYFILETYPE    = '      '
YSANDFILETYPE    = '      '
YORGMATFILETYPE  = '      '
YDENSITYFILETYPE = '      '
YCTIFILETYPE     = '      '
YRUNOFFBFILETYPE = '      '
YWDRAINFILETYPE  = '      ' 
!
LIMP_CLAY        = .FALSE.
LIMP_SAND        = .FALSE.
LIMP_ORGMAT      = .FALSE.
LIMP_DENSITY     = .FALSE.
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
HISBA            = CISBA            ! ISBA option
HPEDOTF          = CPEDO_FUNCTION   ! Pedo-transfert function for DIF
HPHOTO           = CPHOTO           ! photosynthesis option
HSAND            = YSAND            ! file name for sand fraction
HCLAY            = YCLAY            ! file name for clay fraction
HORGMAT          = YORGMAT          ! file name for organic matter
HDENSITY         = YDENSITY         ! file name for soil density
HCTI             = YCTI             ! file name for topographic index
HRUNOFFB         = YRUNOFFB         ! file name for subgrid runoff
HWDRAIN          = YWDRAIN          ! file name for subgrid drainage
HSANDFILETYPE    = YSANDFILETYPE    ! sand data file type
HCLAYFILETYPE    = YCLAYFILETYPE    ! clay data file type
HORGMATFILETYPE  = YORGMATFILETYPE  ! organic matter data file type
HDENSITYFILETYPE = YDENSITYFILETYPE ! soil density data file type
HCTIFILETYPE     = YCTIFILETYPE     ! topographic index data file type
HRUNOFFBFILETYPE = YRUNOFFBFILETYPE ! subgrid runoff data file type
HWDRAINFILETYPE  = YWDRAINFILETYPE  ! subgrid drainage data file type
PUNIF_SAND       = XUNIF_SAND       ! uniform value of sand fraction
PUNIF_CLAY       = XUNIF_CLAY       ! uniform value of clay fraction
PUNIF_ORGMAT     = XUNIF_ORGMAT     ! uniform value of organic matter
PUNIF_DENSITY    = XUNIF_DENSITY    ! uniform value of soil density
PUNIF_RUNOFFB    = XUNIF_RUNOFFB    ! uniform value of subgrid runoff coefficient
PUNIF_WDRAIN     = XUNIF_WDRAIN     ! uniform value of subgrid drainage coefficient
OIMP_SAND        = LIMP_SAND        ! Imposed values for SAND
OIMP_CLAY        = LIMP_CLAY        ! Imposed values for CLAY
OIMP_ORGMAT      = LIMP_ORGMAT      ! Imposed values for organic matter
OIMP_DENSITY     = LIMP_DENSITY     ! Imposed values for soil density
OIMP_CTI         = LIMP_CTI         ! Imposed values for topographic index statistics
IF (LHOOK) CALL DR_HOOK('READ_NAM_PGD_ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_NAM_PGD_ISBA
