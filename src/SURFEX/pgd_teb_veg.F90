!     #########
      SUBROUTINE PGD_TEB_VEG (DTCO, DTGD, DTGR, UG, U, USS, TGDO, TGDP, TGRO, TGRP, TG, &
                              TIR, TOP, TVG, &
                              HPROGRAM)
!     ##############################################################
!
!!**** *PGD_TEB_VEG* monitor for averaging and interpolations of physiographic fields
!!                   for natural covers of TEB
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
!!    J.Escobar   11/2013   Add USE MODI_PGD_TEB_GREENROOF
!!    V. Masson   04/2014   Adds Irrigation
!!    P. Samuelsson 02/2014 Introduced dummy variable in call to READ_NAM_PGD_ISBA for MEB
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
!
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_TEB_GARDEN_n, ONLY : DATA_TEB_GARDEN_t
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DATA_TEB_GREENROOF_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TEB_GREENROOF_PGD_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
!
USE MODD_PGD_GRID,          ONLY : NL
USE MODD_DATA_COVER_PAR,    ONLY : NVEGTYPE
!
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_ISBA_PAR,       ONLY : NOPTIMLAYER, XOPTIMGRID
!
USE MODI_GET_LUOUT
USE MODI_READ_NAM_PGD_ISBA
USE MODI_PGD_FIELD
USE MODI_TEST_NAM_VAR_SURF
!
USE MODI_PGD_TEB_GREENROOF
USE MODI_PGD_TEB_GARDEN_PAR
USE MODI_PGD_TEB_IRRIG
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
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_TEB_GARDEN_t), INTENT(INOUT) :: DTGD
TYPE(DATA_TEB_GREENROOF_t), INTENT(INOUT) :: DTGR
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_GREENROOF_PGD_t), INTENT(INOUT) :: TGRP
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
!
 CHARACTER(LEN=6), INTENT(IN)  :: HPROGRAM   ! Type of program
!                                           ! F if all parameters must be specified
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER                  :: ILUOUT    ! output listing logical unit
INTEGER                  :: JLAYER    ! loop counter
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
 CHARACTER(LEN=28)        :: YSOC_TOP      ! file name for organic carbon
 CHARACTER(LEN=28)        :: YSOC_SUB      ! file name for organic carbon
 CHARACTER(LEN=28)        :: YPERM         ! file name for permafrost distribution
 CHARACTER(LEN=28)        :: YGW               ! file name for groundwater map
 CHARACTER(LEN=6)         :: YSOCFILETYPE  ! organic carbon data file type
 CHARACTER(LEN=6)         :: YPERMFILETYPE ! permafrost distribution data file type
 CHARACTER(LEN=6)         :: YGWFILETYPE      ! groundwater distribution data file type
REAL                     :: XUNIF_SOC_TOP ! uniform value of organic carbon top soil (kg/m2)
REAL                     :: XUNIF_SOC_SUB ! uniform value of organic carbon sub soil (kg/m2)
REAL                     :: XUNIF_PERM    ! uniform permafrost distribution
REAL                     :: XUNIF_GW         ! uniform groundwater distribution
LOGICAL                  :: LIMP_SOC      ! Imposed maps of organic carbon
LOGICAL                  :: LIMP_PERM     ! Imposed maps of permafrost distribution
LOGICAL                  :: GMEB          ! Multi-energy balance (MEB)
LOGICAL                  :: LIMP_GW          ! Imposed maps of groundwater distribution
 CHARACTER(LEN=28)        :: YPH           ! file name for pH
 CHARACTER(LEN=28)        :: YFERT         ! file name for fertilisation rate
 CHARACTER(LEN=6)         :: YPHFILETYPE   ! pH data file type
 CHARACTER(LEN=6)         :: YFERTFILETYPE ! fertilisation data file type
REAL                     :: XUNIF_PH      ! uniform value of pH
REAL                     :: XUNIF_FERT    ! uniform value of fertilisation rate
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PGD_TEB_VEG',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*    1.      Reading of namelist NAM_ISBA for general options of vegetation
!             --------------------------------------------------------------
!
TGDO%NGROUND_LAYER = 0
TVG%CISBA         = '   '
TVG%CPEDOTF       = '   '
TVG%CPHOTO        = '   '
!
CALL READ_NAM_PGD_ISBA(HPROGRAM, IPATCH, IGROUND_LAYER,                         &
                       YISBA, YPEDOTF, YPHOTO,  GTR_ML, ZRM_PATCH,              &
                       YCLAY, YCLAYFILETYPE, XUNIF_CLAY, LIMP_CLAY,             &
                       YSAND, YSANDFILETYPE, XUNIF_SAND, LIMP_SAND,             &
                       YSOC_TOP, YSOC_SUB, YSOCFILETYPE, XUNIF_SOC_TOP,         &
                       XUNIF_SOC_SUB, LIMP_SOC, YCTI, YCTIFILETYPE, LIMP_CTI,   &
                       YPERM, YPERMFILETYPE, XUNIF_PERM, LIMP_PERM, GMEB,       &                     
                       YGW, YGWFILETYPE, XUNIF_GW, LIMP_GW,                     &                     
                       YRUNOFFB, YRUNOFFBFILETYPE, XUNIF_RUNOFFB,               &
                       YWDRAIN,  YWDRAINFILETYPE , XUNIF_WDRAIN, ZSOILGRID,     &
                       YPH, YPHFILETYPE, XUNIF_PH, YFERT, YFERTFILETYPE,        &
                       XUNIF_FERT                                               )  
!
TGDO%NGROUND_LAYER = IGROUND_LAYER
TVG%CISBA         = YISBA
TVG%CPEDOTF       = YPEDOTF
TVG%CPHOTO        = YPHOTO
TVG%LTR_ML        = GTR_ML
!
!-------------------------------------------------------------------------------
!
!*    2.      Coherence of options
!             --------------------
!
  CALL TEST_NAM_VAR_SURF(ILUOUT,'CISBA',TVG%CISBA,'2-L','3-L','DIF')
  CALL TEST_NAM_VAR_SURF(ILUOUT,'CPEDOTF',TVG%CPEDOTF,'CH78','CO84')
  CALL TEST_NAM_VAR_SURF(ILUOUT,'CPHOTO',TVG%CPHOTO,'NON','AGS','LAI','AST','LST','NIT','NCB')
  !
  IF (TVG%CPHOTO=='NCB') THEN
    TVG%CPHOTO = 'NIT'
    WRITE(ILUOUT,*) '****************************************************************'
    WRITE(ILUOUT,*) '* FOR GARDENS, AGS OPTION HAS BEEN CHANGED FROM "NCB" TO "NIT" *'
    WRITE(ILUOUT,*) '****************************************************************'
  END IF
!
  SELECT CASE (TVG%CISBA)
    CASE ('2-L')
      TGDO%NGROUND_LAYER = 2
      TVG%CPEDOTF       ='CH78'       
      WRITE(ILUOUT,*) '*****************************************'
      WRITE(ILUOUT,*) '* With option CISBA = ',TVG%CISBA,'         *'
      WRITE(ILUOUT,*) '* the number of soil layers is set to 2 *'
      WRITE(ILUOUT,*) '* theta(psi) function = Brook and Corey *'
      WRITE(ILUOUT,*) '* Pedo transfert function = CH78        *'          
      WRITE(ILUOUT,*) '*****************************************'
    CASE ('3-L')
      TGDO%NGROUND_LAYER = 3
      TVG%CPEDOTF       ='CH78'         
      WRITE(ILUOUT,*) '*****************************************'
      WRITE(ILUOUT,*) '* With option CISBA = ',TVG%CISBA,'         *'
      WRITE(ILUOUT,*) '* the number of soil layers is set to 3 *'
      WRITE(ILUOUT,*) '* theta(psi) function = Brook and Corey *'
      WRITE(ILUOUT,*) '* Pedo transfert function = CH78        *'        
      WRITE(ILUOUT,*) '*****************************************'
    CASE ('DIF')
      IF(TGDO%NGROUND_LAYER==NUNDEF)THEN
        IF(TOP%LECOCLIMAP)THEN
          TGDO%NGROUND_LAYER=NOPTIMLAYER
        ELSE
          WRITE(ILUOUT,*) '****************************************'
          WRITE(ILUOUT,*) '* Number of ground layer not specified *'
          WRITE(ILUOUT,*) '****************************************'
          CALL ABOR1_SFX('PGD_TEB_GARDEN: NGROUND_LAYER MUST BE DONE IN NAM_ISBA')
        ENDIF
      ENDIF
! 
      ALLOCATE(TGDO%XSOILGRID(TGDO%NGROUND_LAYER))
      TGDO%XSOILGRID(:)=XUNDEF
      TGDO%XSOILGRID(:)=ZSOILGRID(1:TGDO%NGROUND_LAYER) 
      IF(ALL(ZSOILGRID(:)==XUNDEF))THEN
        IF(TOP%LECOCLIMAP) TGDO%XSOILGRID(1:TGDO%NGROUND_LAYER)=XOPTIMGRID(1:TGDO%NGROUND_LAYER)
      ELSEIF(COUNT(TGDO%XSOILGRID/=XUNDEF)/=TGDO%NGROUND_LAYER)THEN
        WRITE(ILUOUT,*) '********************************************************'
        WRITE(ILUOUT,*) '* Soil grid reference values /= number of ground layer *'
        WRITE(ILUOUT,*) '********************************************************'
        CALL ABOR1_SFX('PGD_TEB_GARDEN: XSOILGRID must be coherent with NGROUND_LAYER in NAM_ISBA')            
      ENDIF
!
      WRITE(ILUOUT,*) '*****************************************'
      WRITE(ILUOUT,*) '* Option CISBA            = ',TVG%CISBA
      WRITE(ILUOUT,*) '* Pedo transfert function = ',TVG%CPEDOTF    
      WRITE(ILUOUT,*) '* Number of soil layers   = ',TGDO%NGROUND_LAYER
      IF(TOP%LECOCLIMAP)THEN
        WRITE(ILUOUT,*) '* Soil layers grid (m)    = ',TGDO%XSOILGRID(1:TGDO%NGROUND_LAYER)
      ENDIF
      WRITE(ILUOUT,*) '*****************************************' 

  END SELECT
!
  SELECT CASE (TVG%CPHOTO)
    CASE ('AGS','LAI','AST','LST')
      TVG%NNBIOMASS = 1
    CASE ('NIT')
      TVG%NNBIOMASS = 3
  END SELECT
  WRITE(ILUOUT,*) '*****************************************'
  WRITE(ILUOUT,*) '* With option CPHOTO = ',TVG%CPHOTO,'               *'
  WRITE(ILUOUT,*) '* the number of biomass pools is set to ', TVG%NNBIOMASS
  WRITE(ILUOUT,*) '*****************************************'
!
!-------------------------------------------------------------------------------
!
!*    3.      Sand fraction
!             -------------
!
ALLOCATE(TGDP%XSAND(TG%NDIM,TGDO%NGROUND_LAYER))
!
IF(LIMP_SAND)THEN
!
  CALL ABOR1_SFX('PGD_TEB_VEG: LIMP_SAND IS NOT CONSISTENT WITH TEB_GARDEN')
!
ELSE
!
 CALL PGD_FIELD(DTCO, UG, U, USS, &
                HPROGRAM,'sand fraction','TWN',YSAND,YSANDFILETYPE,XUNIF_SAND,TGDP%XSAND(:,1))
ENDIF
!
DO JLAYER=1,TGDO%NGROUND_LAYER
  TGDP%XSAND(:,JLAYER) = TGDP%XSAND(:,1)
END DO
!-------------------------------------------------------------------------------
!
!*    4.      Clay fraction
!             -------------
!
ALLOCATE(TGDP%XCLAY(TG%NDIM,TGDO%NGROUND_LAYER))
!
IF(LIMP_CLAY)THEN
!
  CALL ABOR1_SFX('PGD_TEB_VEG: LIMP_SAND IS NOT CONSISTENT WITH TEB_GARDEN')
!
ELSE
 CALL PGD_FIELD(DTCO, UG, U, USS, &
                HPROGRAM,'clay fraction','TWN',YCLAY,YCLAYFILETYPE,XUNIF_CLAY,TGDP%XCLAY(:,1))
ENDIF
!
DO JLAYER=1,TGDO%NGROUND_LAYER
  TGDP%XCLAY(:,JLAYER) = TGDP%XCLAY(:,1)
END DO
!-------------------------------------------------------------------------------
!
!*    5.      Subgrid runoff 
!             --------------
!
ALLOCATE(TGDP%XRUNOFFB(TG%NDIM))
 CALL PGD_FIELD(DTCO, UG, U, USS, &
                HPROGRAM,'subgrid runoff','TWN',YRUNOFFB,YRUNOFFBFILETYPE,XUNIF_RUNOFFB,TGDP%XRUNOFFB(:))
!
!-------------------------------------------------------------------------------
!
!*    6.      Drainage coefficient
!             --------------------
!
ALLOCATE(TGDP%XWDRAIN(TG%NDIM))
 CALL PGD_FIELD(DTCO, UG, U, USS, &
                HPROGRAM,'subgrid drainage','TWN',YWDRAIN,YWDRAINFILETYPE,XUNIF_WDRAIN,TGDP%XWDRAIN(:))
!
!-------------------------------------------------------------------------------
!
!*    7.      Interpolation of GARDEN physiographic fields
!             --------------------------------------------
!
DTGD%NTIME = 12
 CALL PGD_TEB_GARDEN_PAR(DTCO, DTGD, UG, U, USS, TGDO, TGDP, TG, &
                         HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*    8.      Case of greenroofs
!             ------------------
!
IF (TOP%LGREENROOF) CALL PGD_TEB_GREENROOF(DTCO, DTGR, UG, U, USS, TGRO, TGRP, TG, &
                                           HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*    9.      Irrigation of gardens and greenroofs
!             ------------------------------------
!
CALL PGD_TEB_IRRIG(DTCO, UG, U, USS, TG, TIR, &
                   HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*    9.      Case of urban hydrology
!             -----------------------
!
IF (TOP%LHYDRO) print*," CALL PGD_TEB_URBHYDRO(HPROGRAM,LECOCLIMAP)"
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PGD_TEB_GARDEN',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE PGD_TEB_VEG
