!     #########
SUBROUTINE COTWORES(PTSTEP, PANF, PABC, PPOI, PAN, PANDAY,                     &
            PANFM,PCSP, PTG,                                                     &
            PSW_RAD, PRA, PQA, PQSAT, PLE, PPSNV, PDELTA, PF2,                   &
            PLAI, PRS, PZENITH, PGC, PRHOA,                                      &
            PFZERO, PEPSO, PGAMM, PQDGAMM, PGMES, PQDGMES,                       &
            PT1GMES, PT2GMES, PAMAX, PQDAMAX, PT1AMAX, PT2AMAX, PIACAN,          &
            PGPP, PRESP_LEAF)  
!   #########################################################################
!
!!****  *COTWORES*  
!!
!!    PURPOSE
!!    -------
!!
!!    Calculates net assimilation of CO2 and leaf conductance.
!!              
!!**  METHOD
!!    ------
!!    Calvet et al. 1998 Forr. Agri. Met. [from model of Jacobs(1994)]
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    USE MODD_CST
!!    USE MODD_CO2V_PAR
!!    USE MODI_COTWO
!!    USE MODI_CCETR
!!    USE MODE_THERMOS
!!
!!    REFERENCE
!!    ---------
!!
!!    Calvet et al. 1998 Forr. Agri. Met. 
!!      
!!    AUTHOR
!!    ------
!!
!!	A. Boone           * Meteo-France *
!!      (following Belair)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/10/97 
!!      V. Masson and V. Rivailland 12/2003 modificatino of ISBA routines order
!!      L. Jarlan   27/10/04 : add of T2 to calculate soil respiration and use
!!                              of CRESPSL key to manage the calculation of soil
!!                              respiration
!!                             PAN et PPST in kgCO2 m-2 s-1 to be fully
!!                              compatible with vegetation growth module (lailoss.f90)
!!      P Le Moigne 09/2005 AGS modifs of L. Jarlan
!!      S. Lafont      03/09 : change units of EPSO GAMM ANDAY
!!      A.L. Gibelin   06/09 : suppress evolution of [CO2] in canopy
!!      A.L. Gibelin   06/09 : move calculations of some CO2 fluxes
!!      A.L. Gibelin   06/09 : add RESP_LEAF
!!      A.L. Gibelin   07/09 : ensure coherence between cotwores and cotworestress
!!      A.L. Gibelin   07/09 : Suppress PPST and PPSTF as outputs, and diagnose GPP
!!
!-------------------------------------------------------------------------------
!
!
USE MODD_CSTS,     ONLY : XMD, XTT, XRHOLW, XLVTT
USE MODD_CO2V_PAR, ONLY : XPARCF, XRDCF, XCONDCTMIN, XDMAX_AGS, XMCO2
!
USE MODI_CCETR
USE MODI_COTWO
!
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
REAL,              INTENT(IN)  :: PTSTEP     ! time step
REAL, DIMENSION(:),INTENT(IN)  :: PABC       ! abscissa needed for integration
!                                            ! of net assimilation and stomatal
!                                            ! conductance over canopy depth
REAL, DIMENSION(:),INTENT(IN)  :: PPOI       ! Gaussian weights (as above)
!
REAL,DIMENSION(:),INTENT(IN):: PCSP, PTG, PF2, PSW_RAD, PRA 
!                                    PCSP  = atmospheric concentration of CO2
!                                    PTG   = updated leaf temperature
!                                    PF2   = normalized soil water stress factor
!                                    PSW_RAD = incident solar radiation
!                                    PRA   = aerodynamic resistance
!
REAL,DIMENSION(:),INTENT(IN):: PQA,  PQSAT, PLE, PPSNV, PDELTA, PLAI, PRHOA
!                                    PQA   = atmospheric mixing ratio
!                                    PQSAT = surface saturation mixing ratio
!                                    PLE   = evapotranspiration (kgH2O kgAir-1 m s-1)
!                                    PPSNV = snow cover fraction
!                                    PDELTA= fraction of the foliage covered
!                                        by intercepted water
!                                    PLAI  = leaf area index
!                                    PRHOA = air density
!
REAL,DIMENSION(:),INTENT(IN):: PZENITH
!	                             PZENITH = solar zenith angle needed 
!                                    for computation of diffusuion of solar
!                                    radiation: for CO2 model.
!
REAL,DIMENSION(:),INTENT(IN):: PFZERO, PEPSO, PGAMM, PQDGAMM, PGMES, PGC,           &
                                 PQDGMES, PT1GMES, PT2GMES, PAMAX, PQDAMAX,           &
                                 PT1AMAX, PT2AMAX  
!                                    PFZERO    = ideal value of F, no photorespiration or 
!                                                saturation deficit
!                                    PEPSO     = maximum initial quantum use efficiency 
!                                                (kgCO2 J-1 PAR)
!                                    PGAMM     = CO2 conpensation concentration (ppmv)
!                                    PQDGAMM   = Q10 function for CO2 conpensation 
!                                                concentration
!                                    PGMES     = mesophyll conductance (m s-1)
!                                    PGC       = cuticular conductance (m s-1)
!                                    PQDGMES   = Q10 function for mesophyll conductance 
!                                    PT1GMES   = reference temperature for computing 
!                                                compensation concentration function for 
!                                                mesophyll conductance: minimum temperature 
!                                    PT2GMES   = reference temperature for computing 
!                                                compensation concentration function for 
!                                                mesophyll conductance: maximum temperature
!                                    PAMAX     = leaf photosynthetic capacity (kg m-2 s-1)
!                                    PQDAMAX   = Q10 function for leaf photosynthetic capacity
!                                    PT1AMAX   = reference temperature for computing 
!                                                compensation concentration function for leaf 
!                                                photosynthetic capacity: minimum temperature
!                                    PT2AMAX   = reference temperature for computing 
!                                                compensation concentration function for leaf 
!                                                photosynthetic capacity: maximum temperature
!
REAL,DIMENSION(:),INTENT(INOUT):: PAN, PANDAY, PRS, PANFM, PGPP
!                                    PAN   = Net assimilation of CO2
!                                    PANDAY= cumulated daily net assimilation of CO2 (kgCO2/m2/day)
!                                    PRS   = stomatal resistance
!                                    PANFM = maximum leaf assimilation
!                                    PGPP  = Gross Primary Production
!
REAL,DIMENSION(:),INTENT(OUT):: PANF
!	                             PANF  = total assimilation over canopy
!
REAL, DIMENSION(:,:), INTENT(OUT) :: PIACAN ! PAR in the canopy at different gauss level
!
REAL,DIMENSION(:),INTENT(OUT):: PRESP_LEAF
!	                             PRESP_LEAF = dark respiration over canopy
!
!*      0.2    declarations of local variables
!
REAL, PARAMETER                :: ZDENOM_MIN  = 1.E-6 ! minimum denominator:
!                                                     ! numerical factor to prevent division by 0
REAL, PARAMETER                :: ZRS_MAX     = 5000. ! maximum canopy resistance (s m-1)
REAL, PARAMETER                :: ZRS_MIN     = 1.E-4 ! minimum canopy resistance (s m-1)
!
INTEGER                     :: JINT ! index for loops
!
REAL, DIMENSION(SIZE(PLAI)) :: ZGAMMT, ZDSP, ZANMAX, ZGMEST, ZDMAX
!                                  ZGAMMT  = compensation point  
!	                           ZDSP    = saturation deficit of atmosphere 
!                                            verses the leaf surface (with correction)
!
!
REAL, DIMENSION(SIZE(PLAI)) :: ZIA, ZTSPC, ZXMUS, ZTAN, ZTGS,    &
                                 ZXIA, ZAN0, ZGS0, ZXTGS, ZRDK  
!	                                    ZIA   = absorbed PAR
!                                           ZTSPC = temperature conversion (K to C) 
!                                           ZXMUS = cosine of solar zenith angle
!                                           ZTAN  = sum for integrated net assimilation 
!                                           ZTGS  = sum for integrated leaf conductance
!                                           ZXIA  = incident radiation after diffusion
!                                           ZAN0  = net asimilation at each interval
!                                                   in the canopy
!                                                   in the canopy
!                                           ZGS0  = leaf conductance at each interval
!                                                   in the canopy        
!                                           ZXTGS = total canopy conductance
!                                           ZRDK = dark respiration
!
REAL, DIMENSION(SIZE(PLAI)) :: ZCONVE1, ZEPSO
!                                           ZEPSO conversion of PEPSO in kgCO2/kgair m/s
!
REAL ZABC
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                           ZABC    = abscissa needed for integration
!                                                     of net assimilation and stomatal 
!                                                     conductance over canopy depth 
!                                                     (working scalar)
!-------------------------------------------------------------------------------
!
! STOMATAL RESISTANCE: ENTRY VARIABLES TO CO2 ROUTINE:
!   CS        = CO2 concentration (kgCO2 kgair-1) cs
!   DSP       = specific humidity deficit (kgH2O kgair-1) ds
!   TSM       = surface temperature (C) ts
!   RG        = global radiation (W m-2) rg
!
! initialisation: convert from ppm to mg/m-3

IF (LHOOK) CALL DR_HOOK('COTWORES',0,ZHOOK_HANDLE)
ZCONVE1(:) = XMCO2*PRHOA(:)/XMD
!
! initialisation: convert from K to C
!
ZTSPC(:)   = PTG(:) - XTT               
!
! absorbed PAR
!
ZIA(:)     = PSW_RAD(:)*XPARCF               
!

!    See (Varlet-Granchet C., M. Chartier, G. Gosse,  and R. Bonhomme, 1981: 
!    Rayonnement utilise pour la photosynthese des vegetaux en
!    conditions naturelles: caracterisation et variations. 
!    Oecol. Plant. 2(16), 189-202.)
!
!
! compensation point (ppm): temperature response
!
ZGAMMT(:) = PGAMM(:)*PQDGAMM(:)**(0.1*(ZTSPC(:)-25.0))
!
! specific humidity deficit (kg kg-1)
!
ZDSP(:)   = MAX( 0.0, PQSAT(:) - PQA(:) - PLE(:)*PRA(:)/(PRHOA*XLVTT) )

! cosine of solar zenith angle 
!
ZXMUS(:) = MAX(COS(PZENITH(:)),0.01)
!
!
! Compute temperature response functions:
!
! kg/m2/s
ZANMAX(:) = ( PAMAX(:)*PQDAMAX(:)**(0.1*(ZTSPC(:)-25.0)) )    &
                 /( (1.0+EXP(0.3*(PT1AMAX(:)-ZTSPC(:))))*       &
                    (1.0+EXP(0.3*(ZTSPC(:)-PT2AMAX(:)))) )  
! m/s
ZGMEST(:) = ( PGMES(:)*PQDGMES(:)**(0.1*(ZTSPC(:)-25.0)) )    &
               /( (1.0+EXP(0.3*(PT1GMES(:)-ZTSPC(:))))*               &
                  (1.0+EXP(0.3*(ZTSPC(:)-PT2GMES(:)))) )  
!
! Add soil moisture stress effect to leaf conductance:
!
ZGMEST(:) = ZGMEST(:)*PF2(:)
!
!
! Integration over the canopy: SIZE(PABC) increments
! are used to approximate the integral.
!
ZTAN(:) = 0.0
ZTGS(:) = 0.0
!
! Unit conversion
! ZANMAX and ZEPSO from kgCO2/m2/s to kgCO2/kgair m/s by dividing by RHOA (kgair/m3)
! ZGAMMT from ppm to kgCO2/kgair
ZANMAX(:)=ZANMAX(:)/PRHOA(:)
ZEPSO(:)=PEPSO(:)/PRHOA(:)
ZGAMMT(:)=ZGAMMT(:)*XMCO2/XMD*1e-6
!
DO JINT = 1, SIZE(PABC)
!
!  Diffusion of incident radiation:
!
   ZABC = PABC(JINT)
!
   CALL CCETR(ZXIA,ZIA,ZXMUS,ZABC,PLAI)
!
!  PAR at different Gauss  level in micmolphot/m2/s
!
   PIACAN(:,JINT)= ZXIA(:)
!
!  Compute conductance and assimilation of CO2: 
!
   ZDMAX = XDMAX_AGS
!
   CALL COTWO(ZAN0, ZGS0, PGC, PCSP,                         &
                ZDSP, ZDMAX, ZXIA, PF2, ZGAMMT,                &
                PFZERO, ZGMEST, ZEPSO, ZANMAX, ZRDK            )  
!
! kgCO2/kgair m/s
  ZTAN(:) = ZTAN(:) + ZAN0(:)*PPOI(JINT) 
  ZTGS(:) = ZTGS(:) + ZGS0(:)*PPOI(JINT) 

END DO

!
!
! Total assimilation
!
PANF(:)= ZTAN(:)
!
! Net assimilation over canopy
!
PAN(:) = (1.0-PDELTA(:))*(1.0-PPSNV(:))*PANF(:)*PLAI(:)
!
! Dark respiration over canopy (does not depend on radiation, 
! no need to integrate over vertical dimension)
!
PRESP_LEAF(:) = (1.0-PDELTA(:))*(1.0-PPSNV(:))*ZRDK(:)*PLAI(:)
!
! Gross primary production over canopy
!
PGPP(:) = PAN(:) + PRESP_LEAF(:)
!
! Daily Net assimilation over canopy (kgCO2/m2/s)
!
PANDAY(:) = PANDAY(:) + PAN(:) * PTSTEP * PRHOA(:)
!
! Adjust maximum leaf assimilation:
!
PANFM(:) = MAX( PANF(:), PANFM(:) )
!
! Total conductance over canopy 
!
ZXTGS(:) = ZTGS(:)*PLAI(:)
!
! Canopy resistance from Ags:
!
PRS(:) = MIN( 1.0/(ZXTGS(:)+ZDENOM_MIN), ZRS_MAX)
PRS(:) = MAX( PRS(:), ZRS_MIN)
IF (LHOOK) CALL DR_HOOK('COTWORES',1,ZHOOK_HANDLE)
!
END SUBROUTINE COTWORES
