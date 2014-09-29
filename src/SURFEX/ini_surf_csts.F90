SUBROUTINE INI_SURF_CSTS 
!     ##################
!
!!****  *INI_SURF_CSTS * - routine to initialize all surface parameter as
!!                         emissivity and anbedo
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!      The physical constants are set to their default numerical values 
!!      or specified in namelist NAM_SURF_CSTS
!!     
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
!!    AUTHOR
!!    ------
!!  	B. Decharme       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2009
!!      M Lafaysse 05/2014 : snow parameters
!!      B. Decharme    05/13 : Add NAM_SURF_REPROD_OPER for versions reproductibility
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,  ONLY : XUNDEF
USE MODD_SURF_CONF, ONLY : CPROGNAME
!
USE MODD_WATER_PAR
USE MODD_FLOOD_PAR
USE MODD_SNOW_PAR,  ONLY : XEMISSN, XANSMIN, XANSMAX, &
                           XAGLAMIN, XAGLAMAX, XHGLA, &
                           XWSNV, XZ0SN, XZ0HSN,      &
                           XTAU_SMELT,                &
                           XALBICE1, XALBICE2, XALBICE3,&
                           XRHOTHRESHOLD_ICE, XZ0ICEZ0SNOW, &
                           XVAGING_NOGLACIER, XVAGING_GLACIER, &
                           XPERCENTAGEPORE
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODE_POS_SURF
!
USE MODD_REPROD_OPER,  ONLY : XEVERG_RSMIN, XEVERG_VEG, &
                                   CDGAVG, CIMPLICIT_WIND,   &
                                   CQSAT, CCHARNOCK, CDGDIF
USE MODI_TEST_NAM_VAR_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER               :: ILUOUT    ! unit of output listing file
INTEGER               :: ILUNAM    ! namelist file  logical unit
LOGICAL               :: GFOUND    ! true if namelist is found
!
LOGICAL               :: LREPROD_OPER
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
NAMELIST/NAM_SURF_CSTS/ XEMISSN, XANSMIN, XANSMAX, XAGLAMIN, XAGLAMAX, &
                        XALBWAT, XALBCOEF_TA96, XALBSCA_WAT, XEMISWAT, &
                        XALBWATICE, XEMISWATICE, XHGLA, XWSNV, XCFFV,  &
                        XZ0SN, XZ0HSN, XTAU_SMELT, XALBSEAICE,         &
                        XZ0FLOOD, XALBWATSNOW
!
NAMELIST/NAM_REPROD_OPER/ LREPROD_OPER, XEVERG_RSMIN, XEVERG_VEG, &
                          CDGAVG, CDGDIF, CIMPLICIT_WIND, CQSAT,  &
                          CCHARNOCK
!
NAMELIST/NAM_SURF_SNOW_CSTS/ XZ0ICEZ0SNOW, XRHOTHRESHOLD_ICE,          &
                             XALBICE1, XALBICE2, XALBICE3,             &
                             XVAGING_NOGLACIER, XVAGING_GLACIER,       &
                             XPERCENTAGEPORE
!
!-------------------------------------------------------------------------------
!*	 0. INIT
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INI_SURF_CSTS',0,ZHOOK_HANDLE)
!
XALBWAT     = XUNDEF
XALBSCA_WAT = XUNDEF
XALBSEAICE  = XUNDEF
XALBWATICE  = XUNDEF
XALBWATSNOW = XUNDEF
XEMISWAT    = XUNDEF
XEMISWATICE = XUNDEF
XEMISSN     = XUNDEF
!
!-------------------------------------------------------------------------------
!*	 1. Default values
!-------------------------------------------------------------------------------
!
! Minimum and maximum values of the albedo of snow:
!
XANSMIN = 0.50 ! (-)
XANSMAX = 0.85 ! (-)
!
! Minimum and maximum values of the albedo of permanet snow/ice:
!
XAGLAMIN = 0.8 ! (-)
XAGLAMAX = 0.85 ! (-)
! 
! Height of aged snow in glacier case (allows Pn=1)
!
XHGLA    = 33.3 !(m)
! 
! Coefficient for calculation of snow fraction over vegetation
!
XWSNV = 5.0 !(-)
!
! Water direct albedo coefficient (option "TA96")
!
XALBCOEF_TA96 =  0.037
!
! Coefficient for calculation of floodplain fraction over vegetation
!
XCFFV = 3.0
!
! Roughness length of pure snow surface (m)
!
XZ0SN = 0.001
!
! Roughness length for heat of pure snow surface (m)
!
XZ0HSN = 0.0001
!
! Snow Melt timescale with D95 (s): needed to prevent time step 
! dependence of melt when snow fraction < unity.
!
XTAU_SMELT = 300.
!
! NAM_SURF_SNOW_CSTS
!
! Roughness length ratio between ice and snow
XZ0ICEZ0SNOW = 10.
!
! 3 bands spectral albedo for glacier ice (CROCUS)
! Default values from Lejeune et al 2009 (Zongo, Bolivia)
XALBICE1 = 0.38
XALBICE2 = 0.23
XALBICE3 = 0.08
!
! Gerbaux et al 2005 (Saint Sorlin)
! PALBICE1=0.23
! PALBICE2=0.16
! PALBICE3=0.05
!
! Density threshold for ice detection kg.m-3
XRHOTHRESHOLD_ICE = 850.
!
! Parameters for ageing effect on albedo
XVAGING_NOGLACIER = 60.
XVAGING_GLACIER   = 900.

! percentage of the total pore volume to compute the max liquid water holding capacity   !Pahaut 1976
XPERCENTAGEPORE = 0.05
!
! Roughness length for flood (m)
!
XZ0FLOOD = 0.0002
!-------------------------------------------------------------------------------
!
! * Reproductibility for SURFEX OPER
!
LREPROD_OPER = .FALSE. ! default
!
! * Vegetation parameters for tropical forest
!
!XEVERG_RSMIN : old = 250. (Manzi 1993) but observations range 
!               from 140 to 180. According to Delire et al. (1997) and 
!               new tests over 6 local sites, 175. is recommended
!               Should be the default after check with AROME/ALADIN
!
XEVERG_RSMIN = 175.  !Rsmin
!
!XEVERG_VEG : old = 0.99 (Manzi 1993) but according to Delire et al. (1997) and 
!             new tests over 6 local sites, 1.0 is recommended because 0.99
!             induces unrealistic bare soil evaporation for Tropical forest
!             Should be the default after check with AROME/ALADIN
!
XEVERG_VEG   = 1.0  !Veg fraction
!
! * Soil depth average
!
CDGAVG = 'INV'
!
! * Soil depth with ISBA-DF
!
CDGDIF = 'ROOT'
!
! * wind implicitation option
!
CIMPLICIT_WIND = 'NEW'
!
! * qsat computation
!
CQSAT = 'NEW'
!
! * Charnock parameter
!
CCHARNOCK = 'NEW'
!
!-------------------------------------------------------------------------------
!*	 2. User values
!-------------------------------------------------------------------------------
!
 CALL GET_LUOUT(CPROGNAME,ILUOUT)
!
 CALL OPEN_NAMELIST(CPROGNAME,ILUNAM)
!
 CALL POSNAM(ILUNAM,'NAM_SURF_CSTS',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_SURF_CSTS)
!
 CALL POSNAM(ILUNAM,'NAM_SURF_SNOW_CSTS',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_SURF_SNOW_CSTS)
!
!-------------------------------------------------------------------------------
!*	 3. For Reproductibility
!-------------------------------------------------------------------------------
!
 CALL POSNAM(ILUNAM,'NAM_REPROD_OPER',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_REPROD_OPER)
!
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CDGAVG',CDGAVG,'ARI','INV')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CDGDIF',CDGDIF,'SOIL','ROOT')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CIMPLICIT_WIND',CIMPLICIT_WIND,'NEW','OLD')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CQSAT',CIMPLICIT_WIND,'NEW','OLD')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CCHARNOCK',CIMPLICIT_WIND,'NEW','OLD')
!
 CALL TEST_NAM_VAR_SURF(ILUOUT,'XEVERG_RSMIN',XEVERG_RSMIN,175.0,250.0)
 CALL TEST_NAM_VAR_SURF(ILUOUT,'XEVERG_VEG',XEVERG_VEG,1.0,0.99) 
!
IF(LREPROD_OPER)THEN
  XEVERG_RSMIN   = 250.
  XEVERG_VEG     = 0.99
  CDGAVG         = 'ARI'
  CDGDIF         = 'SOIL'
  CQSAT          = 'OLD'
  CCHARNOCK      = 'OLD'
ENDIF
!
! Water diffuse albedo
IF(XALBSCA_WAT==XUNDEF)THEN
  IF(LREPROD_OPER)THEN
    XALBSCA_WAT =  0.06
  ELSE
    XALBSCA_WAT =  0.065
  ENDIF
ENDIF
!
! Water global albedo (option "UNIF")
!
IF(XALBWAT==XUNDEF)THEN
  IF(LREPROD_OPER)THEN
    XALBWAT =  0.135
  ELSE
    XALBWAT =  0.065
  ENDIF
ENDIF
!
! Sea ice albedo
!
IF(XALBSEAICE==XUNDEF)THEN
  IF(LREPROD_OPER)THEN
    XALBSEAICE =  0.85
  ELSE
    XALBSEAICE =  0.71
  ENDIF
ENDIF
!
! water ice and snow albedo
!
IF(XALBWATICE==XUNDEF)THEN
  IF(LREPROD_OPER)THEN
    XALBWATICE =  0.85
  ELSE
    XALBWATICE =  0.40
  ENDIF
ENDIF
!
IF(XALBWATSNOW==XUNDEF)THEN
  IF(LREPROD_OPER)THEN
    XALBWATSNOW =  0.85
  ELSE
    XALBWATSNOW =  0.60
  ENDIF
ENDIF
!                   
! Water emissivity
!
IF(XEMISWAT==XUNDEF)THEN
  IF(LREPROD_OPER)THEN
    XEMISWAT =  0.98
  ELSE
    XEMISWAT =  0.96
  ENDIF
ENDIF
!
! Sea ice emissivity
!
IF(XEMISWATICE==XUNDEF)THEN
  IF(LREPROD_OPER)THEN
    XEMISWATICE =  1.0
  ELSE
    XEMISWATICE =  0.97
  ENDIF
ENDIF
!
!
! Snow emissivity:
!
IF(XEMISSN==XUNDEF)THEN
  IF(LREPROD_OPER)THEN
    XEMISSN =  1.0
  ELSE
    XEMISSN =  0.99
  ENDIF
ENDIF
!
!-------------------------------------------------------------------------------
!
 CALL CLOSE_NAMELIST(CPROGNAME,ILUNAM)
!
IF (LHOOK) CALL DR_HOOK('INI_SURF_CSTS',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_SURF_CSTS 
