!     #########
      SUBROUTINE PGD_TEB_GARDEN(HPROGRAM,OECOCLIMAP)
!     ##############################################################
!
!!**** *PGD_TEB_GARDEN* monitor for averaging and interpolations of TEB physiographic fields
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
!!    Original    03/2010
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_PGD_GRID,       ONLY : NL
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_TEB_n,          ONLY : XCOVER, LCOVER, XZS
USE MODD_TEB_GARDEN_n,   ONLY : NPATCH, NGROUND_LAYER, NNBIOMASS,        &
                                  CISBA, CPHOTO, CPEDOTF, LTR_ML,        &
                                  XCLAY, XSAND, XRUNOFFB, XWDRAIN,       &
                                  XZ0EFFJPDIR,                           &
                                  XAOSIP, XAOSIM, XAOSJP, XAOSJM,        &
                                  XHO2IP, XHO2IM, XHO2JP, XHO2JM,        &
                                  XSSO_SLOPE, XSOILGRID  
USE MODD_TEB_GRID_n,     ONLY : CGRID, XGRID_PAR, XLAT, XLON, XMESH_SIZE, NDIM
USE MODD_DATA_TEB_GARDEN_n, ONLY : NTIME
!
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_ISBA_PAR,       ONLY : NOPTIMLAYER, XOPTIMGRID
!
USE MODI_GET_LUOUT
USE MODI_READ_NAM_PGD_ISBA
USE MODI_PGD_FIELD
USE MODI_TEST_NAM_VAR_SURF
!
USE MODI_PGD_TEB_GARDEN_PAR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=6), INTENT(IN)  :: HPROGRAM   ! Type of program
LOGICAL,          INTENT(IN)  :: OECOCLIMAP ! T if parameters are computed with ecoclimap
!                                           ! F if all parameters must be specified
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER                           :: ILUOUT    ! output listing logical unit
INTEGER                           :: JLAYER    ! loop counter
!
!*    0.3    Declaration of namelists
!            ------------------------
!
INTEGER                  :: IPATCH           ! number of patches
INTEGER                  :: IGROUND_LAYER    ! number of soil layers
CHARACTER(LEN=3)         :: YISBA            ! ISBA option
CHARACTER(LEN=4)         :: YPEDOTF          ! Pedo-transfert function for DIF
CHARACTER(LEN=3)         :: YPHOTO           ! photosynthesis option
LOGICAL                  :: GTR_ML           ! new radiative transfert
REAL                     :: ZRM_PATCH        ! threshold to remove little fractions of patches
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
REAL                     :: XUNIF_SAND       ! uniform value of sand fraction
REAL                     :: XUNIF_CLAY       ! uniform value of clay fraction
REAL                     :: XUNIF_RUNOFFB    ! uniform value of subgrid runoff coefficient
REAL                     :: XUNIF_WDRAIN     ! uniform subgrid drainage parameter
LOGICAL                  :: LIMP_SAND        ! Imposed maps of Sand
LOGICAL                  :: LIMP_CLAY        ! Imposed maps of Clay
LOGICAL                  :: LIMP_CTI         ! Imposed maps of topographic index statistics
REAL, DIMENSION(150)     :: ZSOILGRID        ! Soil layer thickness for DIF
!
! Not used in TEB garden
!
CHARACTER(LEN=28)        :: YSOM_TOP      ! file name for organic carbon
CHARACTER(LEN=28)        :: YSOM_SUB      ! file name for organic carbon
CHARACTER(LEN=6)         :: YSOMFILETYPE  ! organic carbon data file type
REAL                     :: XUNIF_SOM     ! uniform value of organic matter (kg/m2)
LOGICAL                  :: LIMP_SOM      ! Imposed maps of organic carbon
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PGD_TEB_GARDEN',0,ZHOOK_HANDLE)
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*    1.      Reading of namelist 
!             -------------------
!
NPATCH        = 0
NGROUND_LAYER = 0
CISBA         = '   '
CPEDOTF       = '   '
CPHOTO        = '   '
!
CALL READ_NAM_PGD_ISBA(HPROGRAM, IPATCH, IGROUND_LAYER,                         &
                       YISBA, YPEDOTF, YPHOTO,  GTR_ML, ZRM_PATCH,              &
                       YCLAY, YCLAYFILETYPE, XUNIF_CLAY, LIMP_CLAY,             &
                       YSAND, YSANDFILETYPE, XUNIF_SAND, LIMP_SAND,             &
                       YSOM_TOP, YSOM_SUB, YSOMFILETYPE, XUNIF_SOM, LIMP_SOM,   &
                       YCTI, YCTIFILETYPE, LIMP_CTI,                            &                       
                       YRUNOFFB, YRUNOFFBFILETYPE, XUNIF_RUNOFFB,               &
                       YWDRAIN,  YWDRAINFILETYPE , XUNIF_WDRAIN, ZSOILGRID      )  
!
NPATCH        = IPATCH
NGROUND_LAYER = IGROUND_LAYER
CISBA         = YISBA
CPEDOTF       = YPEDOTF
CPHOTO        = YPHOTO
LTR_ML        = GTR_ML
!
!-------------------------------------------------------------------------------
!
!*    2.      Coherence of options
!             --------------------
!
  CALL TEST_NAM_VAR_SURF(ILUOUT,'CISBA',CISBA,'2-L','3-L','DIF')
  CALL TEST_NAM_VAR_SURF(ILUOUT,'CPEDOTF',CPEDOTF,'CH78','CO84')  
  CALL TEST_NAM_VAR_SURF(ILUOUT,'CPHOTO',CPHOTO,'NON','AGS','LAI','AST','LST','NIT','NCB')
!
  SELECT CASE (CISBA)
!  
    CASE ('2-L')
!            
      NGROUND_LAYER = 2
      CPEDOTF       ='CH78'       
      WRITE(ILUOUT,*) '*****************************************'
      WRITE(ILUOUT,*) '* With option CISBA = ',CISBA,'         *'
      WRITE(ILUOUT,*) '* the number of soil layers is set to 2 *'
      WRITE(ILUOUT,*) '* theta(psi) function = Brook and Corey *'
      WRITE(ILUOUT,*) '* Pedo transfert function = CH78        *'          
      WRITE(ILUOUT,*) '*****************************************'
!      
    CASE ('3-L')
!            
      NGROUND_LAYER = 3
      CPEDOTF       ='CH78'         
      WRITE(ILUOUT,*) '*****************************************'
      WRITE(ILUOUT,*) '* With option CISBA = ',CISBA,'         *'
      WRITE(ILUOUT,*) '* the number of soil layers is set to 3 *'
      WRITE(ILUOUT,*) '* theta(psi) function = Brook and Corey *'
      WRITE(ILUOUT,*) '* Pedo transfert function = CH78        *'        
      WRITE(ILUOUT,*) '*****************************************'
!
    CASE ('DIF')
!
    IF(NGROUND_LAYER==NUNDEF)THEN
      IF(OECOCLIMAP)THEN
        NGROUND_LAYER=NOPTIMLAYER
      ELSE
        WRITE(ILUOUT,*) '****************************************'
        WRITE(ILUOUT,*) '* Number of ground layer not specified *'
        WRITE(ILUOUT,*) '****************************************'
        CALL ABOR1_SFX('PGD_TEB_GARDEN: NGROUND_LAYER MUST BE DONE IN NAM_ISBA')
      ENDIF
    ENDIF
! 
    ALLOCATE(XSOILGRID(NGROUND_LAYER))
    XSOILGRID(:)=XUNDEF
    XSOILGRID(:)=ZSOILGRID(1:NGROUND_LAYER) 
    IF(ALL(ZSOILGRID(:)==XUNDEF))THEN
      IF(OECOCLIMAP) XSOILGRID(1:NGROUND_LAYER)=XOPTIMGRID(1:NGROUND_LAYER)
    ELSEIF(COUNT(XSOILGRID/=XUNDEF)/=NGROUND_LAYER)THEN
      WRITE(ILUOUT,*) '********************************************************'
      WRITE(ILUOUT,*) '* Soil grid reference values /= number of ground layer *'
      WRITE(ILUOUT,*) '********************************************************'
      CALL ABOR1_SFX('PGD_TEB_GARDEN: XSOILGRID must be coherent with NGROUND_LAYER in NAM_ISBA')            
    ENDIF
!
    WRITE(ILUOUT,*) '*****************************************'
    WRITE(ILUOUT,*) '* Option CISBA            = ',CISBA
    WRITE(ILUOUT,*) '* Pedo transfert function = ',CPEDOTF    
    WRITE(ILUOUT,*) '* Number of soil layers   = ',NGROUND_LAYER
    IF(OECOCLIMAP)THEN
      WRITE(ILUOUT,*) '* Soil layers grid (m)    = ',XSOILGRID(1:NGROUND_LAYER)
    ENDIF
    WRITE(ILUOUT,*) '*****************************************'
!      
  END SELECT
!
  SELECT CASE (CPHOTO)
    CASE ('AGS','LAI','AST','LST')
      NNBIOMASS = 1
    CASE ('NIT')
      NNBIOMASS = 3
    CASE ('NCB')
      NNBIOMASS = 6
  END SELECT
  WRITE(ILUOUT,*) '*****************************************'
  WRITE(ILUOUT,*) '* With option CPHOTO = ',CPHOTO,'               *'
  WRITE(ILUOUT,*) '* the number of biomass pools is set to ', NNBIOMASS
  WRITE(ILUOUT,*) '*****************************************'
!
  IF (NPATCH<1 .OR. NPATCH>NVEGTYPE) THEN
    WRITE(ILUOUT,*) '*****************************************'
    WRITE(ILUOUT,*) '* Number of patch must be between 1 and ', NVEGTYPE
    WRITE(ILUOUT,*) '* You have chosen NPATCH = ', NPATCH
    WRITE(ILUOUT,*) '*****************************************'
    CALL ABOR1_SFX('PGD_TEB_GARDEN: NPATCH MUST BE BETWEEN 1 AND NVEGTYPE')
  END IF
!
  IF ( CPHOTO/='NON' .AND. NPATCH/=12 ) THEN
    WRITE(ILUOUT,*) '*****************************************'
    WRITE(ILUOUT,*) '* With option CPHOTO = ', CPHOTO
    WRITE(ILUOUT,*) '* Number of patch must be equal to 12 '
    WRITE(ILUOUT,*) '* But you have chosen NPATCH = ', NPATCH
    WRITE(ILUOUT,*) '*****************************************'
    CALL ABOR1_SFX('PGD_TEB_GARDEN: CPHOTO='//CPHOTO//' REQUIRES NPATCH=12')
  END IF
!
!-------------------------------------------------------------------------------
!
!*    3.      Sand fraction
!             -------------
!
ALLOCATE(XSAND(NDIM,NGROUND_LAYER))
!
IF(LIMP_SAND)THEN
!
  CALL ABOR1_SFX('PGD_TEB_GARDEN: LIMP_SAND IS NOT CONSISTENT WITH TEB_GARDEN')
!
ELSE
!
CALL PGD_FIELD(HPROGRAM,'sand fraction','TWN',YSAND,YSANDFILETYPE,XUNIF_SAND,XSAND(:,1))
ENDIF
!
DO JLAYER=1,NGROUND_LAYER
  XSAND(:,JLAYER) = XSAND(:,1)
END DO
!-------------------------------------------------------------------------------
!
!*    4.      Clay fraction
!             -------------
!
ALLOCATE(XCLAY(NDIM,NGROUND_LAYER))
!
IF(LIMP_CLAY)THEN
!
  CALL ABOR1_SFX('PGD_TEB_GARDEN: LIMP_SAND IS NOT CONSISTENT WITH TEB_GARDEN')
!
ELSE
CALL PGD_FIELD(HPROGRAM,'clay fraction','TWN',YCLAY,YCLAYFILETYPE,XUNIF_CLAY,XCLAY(:,1))
ENDIF
!
DO JLAYER=1,NGROUND_LAYER
  XCLAY(:,JLAYER) = XCLAY(:,1)
END DO
!-------------------------------------------------------------------------------
!
!*    5.      Subgrid runoff 
!             --------------
!
ALLOCATE(XRUNOFFB(NDIM))
CALL PGD_FIELD                                                                              &
       (HPROGRAM,'subgrid runoff','TWN',YRUNOFFB,YRUNOFFBFILETYPE,XUNIF_RUNOFFB,XRUNOFFB(:))  
!
!-------------------------------------------------------------------------------
!
!*    6.      Drainage coefficient
!             --------------------
!
ALLOCATE(XWDRAIN(NDIM))
CALL PGD_FIELD                                                                              &
       (HPROGRAM,'subgrid drainage','TWN',YWDRAIN,YWDRAINFILETYPE,XUNIF_WDRAIN,XWDRAIN(:))  
!
!-------------------------------------------------------------------------------
!
!*    7.      Fileds prescribed to 0
!             ----------------------
!
ALLOCATE(XZ0EFFJPDIR(NDIM))
ALLOCATE(XAOSIP     (NDIM))
ALLOCATE(XAOSIM     (NDIM))
ALLOCATE(XAOSJP     (NDIM))
ALLOCATE(XAOSJM     (NDIM))
ALLOCATE(XHO2IP     (NDIM))
ALLOCATE(XHO2IM     (NDIM))
ALLOCATE(XHO2JP     (NDIM))
ALLOCATE(XHO2JM     (NDIM))
ALLOCATE(XSSO_SLOPE (NDIM))
!
XZ0EFFJPDIR(:) = 0.
XAOSIP     (:) = 0.
XAOSIM     (:) = 0.
XAOSJP     (:) = 0.
XAOSJM     (:) = 0.
XHO2IP     (:) = 0.
XHO2IM     (:) = 0.
XHO2JP     (:) = 0.
XHO2JM     (:) = 0.
XSSO_SLOPE (:) = 0.
!
!-------------------------------------------------------------------------------
!
!*    8.      ISBA specific fields
!             --------------------
!
NTIME = 12
CALL PGD_TEB_GARDEN_PAR(HPROGRAM)
!
IF (LHOOK) CALL DR_HOOK('PGD_TEB_GARDEN',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE PGD_TEB_GARDEN
