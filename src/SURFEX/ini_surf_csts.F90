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
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_CONF, ONLY : CPROGNAME
!
USE MODD_WATER_PAR
USE MODD_FLOOD_PAR
USE MODD_SNOW_PAR,  ONLY : XEMISSN, XANSMIN, XANSMAX, &
                           XAGLAMIN, XAGLAMAX, XHGLA, &
                           XWSNV, XZ0SN, XZ0HSN, &
                           XTAU_SMELT,&
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
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER               :: ILUOUT    ! unit of output listing file
INTEGER               :: ILUNAM    ! namelist file  logical unit
LOGICAL               :: GFOUND    ! true if namelist is found
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
NAMELIST/NAM_SURF_CSTS/ XEMISSN, XANSMIN, XANSMAX, XAGLAMIN, XAGLAMAX, &
                        XALBWAT, XALBCOEF_TA96, XALBSCA_WAT, XEMISWAT, &
                        XALBWATICE, XEMISWATICE, XHGLA, XWSNV, XCFFV,  &
                        XZ0SN, XZ0HSN, XTAU_SMELT  
NAMELIST/NAM_SURF_SNOW_CSTS/ XZ0ICEZ0SNOW, XRHOTHRESHOLD_ICE,          &
                             XALBICE1, XALBICE2, XALBICE3,             &
                             XVAGING_NOGLACIER, XVAGING_GLACIER,       &
                             XPERCENTAGEPORE
!
!-------------------------------------------------------------------------------
!*	 1. Default values
!-------------------------------------------------------------------------------
!
! Snow emissivity:
!
IF (LHOOK) CALL DR_HOOK('INI_SURF_CSTS',0,ZHOOK_HANDLE)
XEMISSN = 1.0  ! (-)
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
! Water global albedo (option "UNIF")
!
XALBWAT =  0.135
!
! Water direct albedo coefficient (option "TA96")
!
XALBCOEF_TA96 =  0.037
!
! Water diffuse albedo
!
XALBSCA_WAT =  0.06
!                        
! Water emissivity
!
XEMISWAT =  0.98                      
!
! Sea ice albedo
!
XALBWATICE = 0.85
!
! Sea ice emissivity
!
XEMISWATICE = 1.0
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
 CALL CLOSE_NAMELIST(CPROGNAME,ILUNAM)

IF (LHOOK) CALL DR_HOOK('INI_SURF_CSTS',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_SURF_CSTS 
