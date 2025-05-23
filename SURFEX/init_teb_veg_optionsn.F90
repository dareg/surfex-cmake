!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!#############################################################
SUBROUTINE INIT_TEB_VEG_OPTIONS_n (CHT, OSURF_DIAG_ALBEDO, OGREENROOF, GDO, GRO, HPROGRAM)
!#############################################################
!
!!****  *INIT_TEB_TEB_VEG_n* - routine to initialize ISBA
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
!!      A. Lemonsu  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2009
!!      B. Decharme 07/2011 : read pgd+prep
!!      B. Decharme 04/2013 : delete CTOPREG option (never used)
!!                            water table / surface coupling
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
USE MODD_CH_TEB_n, ONLY : CH_TEB_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
!
USE MODD_READ_NAMELIST,   ONLY : LNAM_READ
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
USE MODD_TYPE_SNOW, ONLY: SURF_SNOW
!
!


USE MODD_DATA_COVER_PAR,  ONLY: NVEGTYPE
USE MODD_SURF_PAR,        ONLY: XUNDEF, NUNDEF
USE MODD_ISBA_ALBEDO,     ONLY: CALBEDO
!
USE MODD_ISBA_PAR,        ONLY : XOPTIMGRID
!
USE MODN_TEB_n,           ONLY : XTSTEP
!
USE MODI_READ_NAM_PGD_ISBA
USE MODI_DEFAULT_ISBA
USE MODI_DEFAULT_CH_DEP
USE MODI_DEFAULT_CH_BIO_FLUX
USE MODI_DEFAULT_CROCUS
USE MODI_READ_DEFAULT_TEB_VEG_n
USE MODI_READ_TEB_VEG_CONF_n
USE MODI_GET_LUOUT
USE MODI_READ_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
!
!
TYPE(CH_TEB_t), INTENT(INOUT) :: CHT
LOGICAL, INTENT(OUT) :: OSURF_DIAG_ALBEDO
LOGICAL, INTENT(IN) :: OGREENROOF 
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GDO
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GRO
!
 CHARACTER(LEN=6),                   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IVERSION, IBUGFIX  ! surface version
INTEGER           :: ILUOUT   ! unit of output listing file
INTEGER           :: IRESP    ! Error code after redding
CHARACTER(LEN=12) :: YRECFM   ! Name of the article to be read
CHARACTER(LEN=4 ) :: YLVL
!
INTEGER :: JLAYER ! loop counter on layers
!
REAL                              :: ZOUT_TSTEP
CHARACTER(LEN=3)                  :: YRAIN 
LOGICAL                           :: GCANOPY_DRAG
LOGICAL                           :: GGLACIER
LOGICAL                           :: GFLOOD
LOGICAL                           :: GWTD
LOGICAL                           :: GVEGUPD
LOGICAL                           :: GSPINUPCARBS
LOGICAL                           :: GSPINUPCARBW
REAL                              :: ZSPINMAXS
REAL                              :: ZSPINMAXW
REAL                              :: ZCO2_START
REAL                              :: ZCO2_END
INTEGER                           :: INBYEARSPINS
INTEGER                           :: INBYEARSPINW
REAL:: ZCVHEATF
INTEGER                  :: IPATCH           ! number of patches
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_VEG_OPTIONS_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!               Other little things
!
OSURF_DIAG_ALBEDO = .FALSE.
!
IF (LNAM_READ) THEN
 !
 !*       1.     Defaults
 !               --------
 !
 !        1.1. Hard defaults
 !      
 !       Definition of default options for ISBA (in MODD_TEB_VEG_n)
 !       REM - TSTEP, OUT_TSTEP, CANOPY_DRAG are defined as local variables
 !             because they are already in init_teb.f90 (these options are 
 !             forced to the same values for TEB and urban green areas)
 !
 CALL DEFAULT_ISBA(GDO%XTSTEP, GDO%XOUT_TSTEP, GDO%CRUNOFF, &
                   GDO%CSCOND, GDO%CC1DRY, GDO%CSOILFRZ,   &
                   GDO%CDIFSFCOND, GDO%CSNOWRES, GDO%CCPSURF, GDO%XCGMAX, &
                   GDO%XCDRAG, GDO%CKSAT, GDO%LSOC, GDO%CRAIN, GDO%CHORT, &
                   GDO%LGLACIER, GDO%LCANOPY_DRAG, GDO%LVEGUPD, GDO%LSPINUPCARBS, &
                   GDO%LSPINUPCARBW, GDO%XSPINMAXS, GDO%XSPINMAXW, GDO%XCO2_START,&
                   GDO%XCO2_END, GDO%NNBYEARSPINS, GDO%NNBYEARSPINW, &
                   GDO%LPERTSURF, GDO%XPERT_LOW, GDO%XPERT_HIGH, &
                   GDO%LNITRO_DILU, GDO%XCSMAX, ZCVHEATF )
 !
 CALL DEFAULT_CH_BIO_FLUX(CHT%LCH_BIO_FLUX)
 !
 CALL DEFAULT_CROCUS(GDO%CSNOWDRIFT,GDO%LSNOWDRIFT_SUBLIM,GDO%LSNOW_ABS_ZENITH,&
                     GDO%CSNOWMETAMO,GDO%CSNOWRAD, GDO%LATMORAD, GDO%LSNOWSYTRON, &
                     GDO%CSNOWFALL, GDO%CSNOWCOND , GDO%CSNOWHOLD, GDO%CSNOWCOMP, &
                     GDO%CSNOWZREF, GDO%LSNOWCOMPACT_BOOL, GDO%LSNOWMAK_BOOL,     &
                     GDO%LPRODSNOWMAK, GDO%LSNOWMAK_PROP, GDO%LSNOWTILLER, GDO%LSELF_PROD) 
 !
ENDIF
!        1.2. Defaults from file header
!    
 CALL READ_DEFAULT_TEB_VEG_n(CHT, GDO, HPROGRAM)
!
 CALL READ_TEB_VEG_CONF_n(CHT, GDO, HPROGRAM)
!
!-------------------------------------------------------------------------------
GDO%NPATCH = 1
GDO%CRESPSL = 'DEF'
GDO%LMEB_GNDRES = .FALSE.
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
!*       2.     Definition of version
!               ---------------------
!
YRECFM='VERSION'
 CALL READ_SURF(HPROGRAM,YRECFM,IVERSION,IRESP)
!
YRECFM='BUG'
 CALL READ_SURF(HPROGRAM,YRECFM,IBUGFIX,IRESP)
!
!*       2.     Initialisation of ISBA options
!               ------------------------------
!
!
YRECFM='TWN_ISBA'
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM='GD_ISBA'
 CALL READ_SURF(HPROGRAM,YRECFM,GDO%CISBA,IRESP)
!
IF (IVERSION>=7) THEN
  !
  YRECFM='TWN_PEDOTF'
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM='GD_PEDOTF'
  CALL READ_SURF(HPROGRAM,YRECFM,GDO%CPEDOTF,IRESP)
  !
ELSE
  GDO%CPEDOTF = 'CH78'
ENDIF
!
YRECFM='TWN_PHOTO'
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM='GD_PHOTO'
 CALL READ_SURF(HPROGRAM,YRECFM,GDO%CPHOTO,IRESP)
!
YRECFM='TWN_LAYER'
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM='GD_LAYER'
 CALL READ_SURF(HPROGRAM,YRECFM,GDO%NGROUND_LAYER,IRESP)
!
ALLOCATE(GDO%LMEB_PATCH(1))
GDO%LMEB_PATCH(:) = .FALSE.
!
!* new radiative transfert
!
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=2) THEN
  !
  YRECFM='TWN_TR_ML'
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM='GD_TR_ML'
  CALL READ_SURF(HPROGRAM,YRECFM,GDO%LTR_ML,IRESP)
  !
ELSE 
  GDO%LTR_ML = .FALSE.
ENDIF
!
IF (IVERSION>8 .OR. IVERSION==8 .AND. IBUGFIX>=1) THEN
  !
  YRECFM='GD_ALBEDO'
  CALL READ_SURF(HPROGRAM,YRECFM,GDO%CALBEDO,IRESP)
  !
ELSE
  !  
  GDO%CALBEDO = CALBEDO
  !
ENDIF
!
!* Reference grid for DIF
!
IF(GDO%CISBA=='DIF') THEN
  ALLOCATE(GDO%XSOILGRID(GDO%NGROUND_LAYER))
  GDO%XSOILGRID=XUNDEF
  IF (IVERSION>=8) THEN
     DO JLAYER=1,GDO%NGROUND_LAYER
        WRITE(YLVL,'(I4)') JLAYER
        YRECFM='GD_SGRID'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
        CALL READ_SURF(HPROGRAM,YRECFM,GDO%XSOILGRID(JLAYER),IRESP)
     ENDDO
  ELSEIF (IVERSION==7 .AND. IBUGFIX>=2) THEN
    YRECFM='TWN_SOILGRID'
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM='GD_SOILGRID'
    CALL READ_SURF(HPROGRAM,YRECFM,GDO%XSOILGRID,IRESP,HDIR='-')
  ELSE
    GDO%XSOILGRID(1:GDO%NGROUND_LAYER)=XOPTIMGRID(1:GDO%NGROUND_LAYER)
  ENDIF
ELSE
  ALLOCATE(GDO%XSOILGRID(0))
ENDIF
!
!* number of biomass pools
!
GDO%NNBIOMASS=1
IF (GDO%CPHOTO=='NIT') GDO%NNBIOMASS=3
!
!-------------------------------------------------------------------------------
!
IF (OGREENROOF) THEN
  !
  GRO%NPATCH = 1
  ALLOCATE(GRO%LMEB_PATCH(1))
  GRO%LMEB_PATCH(:) = .FALSE.
  GRO%CRESPSL   = GDO%CRESPSL  
  !
  YRECFM='GR_ISBA'
  CALL READ_SURF(HPROGRAM,YRECFM,GRO%CISBA,IRESP)
  !
  GRO%CPEDOTF   = GDO%CPEDOTF  
  GRO%CPHOTO    = GDO%CPHOTO
  !
  YRECFM='GR_LAYER'
  CALL READ_SURF( HPROGRAM,YRECFM,GRO%NGROUND_LAYER,IRESP)
  !
  GRO%LTR_ML = .FALSE.
  !
  GRO%NNBIOMASS = GDO%NNBIOMASS
  !
  YRECFM='GR_SCOND'
  CALL READ_SURF(HPROGRAM,YRECFM,GRO%CSCOND,IRESP)
  !
  GRO%XTSTEP = GDO%XTSTEP
  GRO%XOUT_TSTEP = GDO%XOUT_TSTEP
  GRO%CRUNOFF = GDO%CRUNOFF
  GRO%CALBEDO = GDO%CALBEDO
  GRO%CC1DRY = GDO%CC1DRY
  GRO%CSOILFRZ = GDO%CSOILFRZ
  GRO%CDIFSFCOND = GDO%CDIFSFCOND 
  GRO%CSNOWRES = GDO%CSNOWRES
  GRO%CCPSURF = GDO%CCPSURF
  GRO%XCGMAX = GDO%XCGMAX 
  GRO%XCSMAX = GDO%XCSMAX 
  GRO%XCDRAG = GDO%XCDRAG 
  GRO%CKSAT = GDO%CKSAT
  GRO%LSOC = GDO%LSOC
  GRO%CRAIN = GDO%CRAIN
  GRO%CHORT = GDO%CHORT
  GRO%LGLACIER = GDO%LGLACIER
  GRO%LCANOPY_DRAG = GDO%LCANOPY_DRAG
  GRO%LVEGUPD = GDO%LVEGUPD
  GRO%LSPINUPCARBS = GDO%LSPINUPCARBS
  GRO%LSPINUPCARBW = GDO%LSPINUPCARBW
  GRO%XSPINMAXS = GDO%XSPINMAXS
  GRO%XSPINMAXW = GDO%XSPINMAXW
  GRO%XCO2_START = GDO%XCO2_START
  GRO%XCO2_END = GDO%XCO2_END
  GRO%NNBYEARSPINS = GDO%NNBYEARSPINS
  GRO%NNBYEARSPINW = GDO%NNBYEARSPINW
  GRO%LNITRO_DILU = GDO%LNITRO_DILU
  !
!  GRO%CSNOWDRIFT = GDO%CSNOWDRIFT
!  GRO%LSNOWDRIFT_SUBLIM = GDO%LSNOWDRIFT_SUBLIM
!  GRO%LSNOW_ABS_ZENITH = GDO%LSNOW_ABS_ZENITH
!  GRO%CSNOWMETAMO = GDO%CSNOWMETAMO
!  GRO%CSNOWRAD = GDO%CSNOWRAD
  !
  GRO%LMEB_GNDRES = GRO%LMEB_GNDRES
  !
  CALL DEFAULT_CROCUS(GRO%CSNOWDRIFT,GRO%LSNOWDRIFT_SUBLIM,GRO%LSNOW_ABS_ZENITH,&
                     GRO%CSNOWMETAMO,GRO%CSNOWRAD, GRO%LATMORAD, GRO%LSNOWSYTRON, &
										 GRO%CSNOWFALL, GRO%CSNOWCOND , GRO%CSNOWHOLD, GRO%CSNOWCOMP, &
										 GRO%CSNOWZREF, GRO%LSNOWCOMPACT_BOOL, GRO%LSNOWMAK_BOOL,     &
										 GRO%LPRODSNOWMAK, GRO%LSNOWMAK_PROP, GRO%LSNOWTILLER, GRO%LSELF_PROD) 
ENDIF
!
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_VEG_OPTIONS_N',1,ZHOOK_HANDLE)
END SUBROUTINE INIT_TEB_VEG_OPTIONS_n
