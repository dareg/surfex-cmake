!     #########
SUBROUTINE COTWORESTRESS(PTSTEP, PVEGTYPE, OSTRESSDEF, PAH, PBH, PF2I,         &
            PANF, PABC, PPOI, PAN, PANDAY, PANFM, PCSP, PTG,                     &
            PSW_RAD, PRA, PQA, PQSAT, PLE, PPSNV, PDELTA, PF2,                   &
            PLAI, PRS, PZENITH, PGC, PRHOA, PDMAX,                               &
            PFZERO, PEPSO, PGAMM, PQDGAMM, PGMES, PQDGMES,                       &
            PT1GMES, PT2GMES, PAMAX, PQDAMAX, PT1AMAX, PT2AMAX, PIACAN,           &
            PGPP, PRESP_LEAF)  
!   #########################################################################
!
!!****  *COTWORESTRESS*  
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
!!        S. Lafont    03/11 : Correct a bug fopr grassland below wilting point
!!        A. Boone     11/11 : Add rsmax to MODD_ISBA_PA
!!
!-------------------------------------------------------------------------------
!
USE MODD_CSTS,           ONLY : XMD, XTT, XRHOLW, XLVTT
USE MODD_ISBA_PAR,       ONLY : XRS_MAX
USE MODD_CO2V_PAR,       ONLY : XPARCF, XDMAXX, XDMAXN, XAW, XBW, XASW, &
                                  XCONDCTMIN, XMCO2  
USE MODD_DATA_COVER_PAR, ONLY : NVT_TREE, NVT_EVER, NVT_CONI
USE MODD_SURF_PAR,       ONLY : XUNDEF
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
REAL,                INTENT(IN) :: PTSTEP    ! time step
REAL, DIMENSION(:,:),INTENT(IN) :: PVEGTYPE
!                                PVEGTYPE  = type de vegetation (1 a 9)
!
REAL, DIMENSION(:),INTENT(IN)  :: PABC       ! abscissa needed for integration
!                                            ! of net assimilation and stomatal
!                                            ! conductance over canopy depth
REAL, DIMENSION(:),INTENT(IN)  :: PPOI       ! Gaussian weights (as above)
!
LOGICAL,DIMENSION(:),INTENT(IN):: OSTRESSDEF
REAL,DIMENSION(:),INTENT(IN):: PAH, PBH, PF2I
!                                    PAH, PBH  = coefficients for universal herbaceous
!                                                stress relation 
!                                    OSTRESSDEF   = water stress vegetation comportement 
!                                                (true:defensif false:offensif)
!                                    PF2I      = critical normalized soil water stress
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
REAL,DIMENSION(:),INTENT(IN):: PFZERO, PEPSO, PGAMM, PQDGAMM, PGMES, PGC, PDMAX,    &
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
!                                    PDMAX     = maximum saturation deficit of 
!                                                  atmosphere tolerate by vegetation
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
!
!*      0.2    declarations of local variables
!
REAL, PARAMETER                :: ZDENOM_MIN  = 1.E-6 ! minimum denominator:
!                                                     ! numerical factor to prevent division by 0
REAL, PARAMETER                :: ZRS_MIN     = 1.E-4 ! minimum canopy resistance (s m-1)
!
INTEGER                     :: JINT ! index for loops
!
REAL, DIMENSION(SIZE(PLAI,1)) :: ZDMAXT, ZDMAXSTAR, ZFZEROSTAR, &
                                   ZFZERON, ZGMESTN, ZFZEROT  
!                                 ZDMAXT     = maximum saturation deficit of atmosphere
!                                              tolerate by vegetation
!                                 ZDMAXSTAR  = maximum saturation deficit of atmosphere
!                                              tolerate by vegetation without soil water stress
!                                 ZFZEROSTAR = initial optimal ratio Ci/Cs for woody vegetation
!                                 ZFZERON  = minimum value for "fzero" in defensive woody strategy
!                                 ZGMESTN  = gmest value at zf2=zf2i in offensive woody strategy
!                                 ZFZEROT  = ideal value of F, no photorespiration or 
!                                            saturation deficit
!
REAL, DIMENSION(SIZE(PLAI,1)) :: ZGAMMT, ZDSP, ZANMAX, ZGMEST
!                                  ZGAMMT  = compensation point 
!	                           ZDSP    = saturation deficit of atmosphere 
!                                            verses the leaf surface (with correction)
!
REAL, DIMENSION(SIZE(PLAI,1)) :: ZIA, ZTSPC, ZXMUS, ZTAN, ZTGS,    &
                                   ZXIA, ZAN0, ZGS0, ZXTGS, ZRDK  
!	                                    ZIA   = absorbed PAR
!                                           ZTSPC = temperature conversion (K to C) 
!                                           ZXMUS = cosine of solar zenith angle
!                                           ZTAN  = sum for integrated net assimilation 
!                                           ZTGS  = sum for integrated leaf conductance
!                                           ZXIA  = incident radiation after diffusion
!                                           ZAN0  = net asimilation at each interval
!                                                   in the canopy
!                                           ZGS0  = leaf conductance at each interval
!                                                   in the canopy        
!                                           ZXTGS = total canopy conductance
!                                           ZRDK = dark respiration
!
REAL, DIMENSION(SIZE(PLAI)) :: ZCONVE1, ZEPSO, ZLAI
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

IF (LHOOK) CALL DR_HOOK('COTWORESTRESS',0,ZHOOK_HANDLE)
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
ZLAI(:)=PLAI(:)
WHERE (PLAI(:)==XUNDEF) ZLAI(:)=0.0

!    See (Varlet-Granchet C., M. Chartier, G. Gosse,  and R. Bonhomme, 1981: 
!    Rayonnement utilise pour la photosynthese des vegetaux en
!    conditions naturelles: caracterisation et variations. 
!    Oecol. Plant. 2(16), 189-202.)
!
!-------------------------------------
! Add soil moisture stress effect to leaf conductance:
! OFFENSIVE and DEFENSIVE water stress response
!
ZGMEST(:)  = PGMES(:)
ZFZEROT(:) = PFZERO(:)
ZDMAXT(:)  = PDMAX(:)
!-------------------
!
! -HERBACEOUS-
!
! DEFENSIVE soil water stress response: 
! PAH and PBH are original coefficients of Calvet 2000
!
ZDMAXSTAR(:) = EXP((LOG(ZGMEST(:)*1000.)-PAH(:))/PBH(:))/1000.
!
WHERE( OSTRESSDEF(:) .AND. (PF2(:) >= PF2I(:))&
   .AND.(PVEGTYPE(:,NVT_TREE) + PVEGTYPE(:,NVT_EVER) + PVEGTYPE(:,NVT_CONI)<0.5))  
  ZDMAXT(:) = ZDMAXSTAR(:) - (ZDMAXSTAR(:)-XDMAXN)* &
                (1.-PF2(:))/(1.-PF2I(:))   
  ZGMEST(:) = EXP(PAH(:)+PBH(:)*LOG(ZDMAXT(:)*1000.))/1000.
END WHERE
WHERE( OSTRESSDEF(:) .AND. (PF2(:) < PF2I(:))&
   .AND.(PVEGTYPE(:,NVT_TREE) + PVEGTYPE(:,NVT_EVER) + PVEGTYPE(:,NVT_CONI)<0.5))  
  ZDMAXT(:)  = XDMAXN
  ZGMEST(:) = (EXP(PAH(:)+PBH(:)*LOG(XDMAXN*1000.))/1000.)*PF2(:)/PF2I(:)
END WHERE
! 
! OFFENSIVE soil water stress response:
!
WHERE ( (.NOT. OSTRESSDEF(:)) .AND. (PF2(:) >= PF2I(:))&
   .AND.(PVEGTYPE(:,NVT_TREE) + PVEGTYPE(:,NVT_EVER) + PVEGTYPE(:,NVT_CONI)<0.5))  
  ZDMAXT(:) = XDMAXX+(ZDMAXSTAR(:)-XDMAXX)*&
                (PF2(:)-PF2I(:))/(1.-PF2I(:))   
  ZGMEST(:) = EXP(PAH(:)+PBH(:)*LOG(ZDMAXT(:)*1000.))/1000.
END WHERE
WHERE ( (.NOT. OSTRESSDEF(:)) .AND. (PF2(:) < PF2I(:))&
   .AND.(PVEGTYPE(:,NVT_TREE) + PVEGTYPE(:,NVT_EVER) + PVEGTYPE(:,NVT_CONI)<0.5))  
  ZDMAXT(:)  = XDMAXX*PF2(:)/PF2I(:)
  ZGMEST(:)  = EXP(PAH(:)+PBH(:)*LOG(XDMAXX*1000.))/1000.
END WHERE
!
!to limit photosynthesis under wilting point
WHERE ( (.NOT. OSTRESSDEF(:)) .AND. (ZDMAXT(:)<=XDMAXN) &
 .AND.(PVEGTYPE(:,NVT_TREE) + PVEGTYPE(:,NVT_EVER) + PVEGTYPE(:,NVT_CONI)<0.5))
  ZDMAXT(:) = XDMAXN
  ZGMEST(:) = (EXP(PAH(:)+PBH(:)*LOG(XDMAXN*1000.))/1000.)*PF2(:)/PF2I(:)
END WHERE
!
!-------------------
!
! -WOODY-
!
! DEFENSIVE soil water stress response:
!
ZFZEROSTAR(:) = ( XAW - LOG(ZGMEST(:)*1000.) )/XBW
ZFZERON(:)    = ( XASW - LOG(ZGMEST(:)*1000.) )/XBW
ZGMESTN(:)    = EXP(XASW - XBW*ZFZEROSTAR(:))/1000.
!
WHERE ( OSTRESSDEF(:) .AND. (PF2(:) >= PF2I(:)) &
   .AND. (PVEGTYPE(:,NVT_TREE) + PVEGTYPE(:,NVT_EVER) + PVEGTYPE(:,NVT_CONI)>0.5))  
  ZFZEROT(:) = ZFZEROSTAR(:) - (ZFZEROSTAR(:)-ZFZERON(:))* &
                 (1.-PF2(:))/(1.-PF2I(:))  
END WHERE
WHERE ( OSTRESSDEF(:) .AND. (PF2(:) < PF2I(:)) &
   .AND. (PVEGTYPE(:,NVT_TREE) + PVEGTYPE(:,NVT_EVER) + PVEGTYPE(:,NVT_CONI)>0.5))  
  ZGMEST(:)  = MAX( 1.0E-10,ZGMEST(:)*PF2(:)/PF2I(:) )
  ZFZEROT(:) = MIN( .95,( XASW-LOG(ZGMEST(:)*1000.) )/XBW )
END WHERE
!
! OFFENSIVE soil water stress response:
!
WHERE ( (.NOT. OSTRESSDEF(:)) .AND. (PF2(:) >= PF2I(:))&
   .AND. (PVEGTYPE(:,NVT_TREE) + PVEGTYPE(:,NVT_EVER) + PVEGTYPE(:,NVT_CONI)>0.5))  
 ZGMEST(:) = ZGMEST(:) - (ZGMEST(:)-ZGMESTN(:))* &
               (1.-PF2(:))/(1.-PF2I(:))  
  ZFZEROT(:) = ZFZEROSTAR(:)
END WHERE 
WHERE ( (.NOT. OSTRESSDEF(:)) .AND. (PF2(:) < PF2I(:))&
   .AND. (PVEGTYPE(:,NVT_TREE) + PVEGTYPE(:,NVT_EVER) + PVEGTYPE(:,NVT_CONI)>0.5))  
  ZGMEST(:) = ZGMESTN(:)*PF2(:)/PF2I(:)
  ZFZEROT(:) = MIN( .95,( XASW-LOG(ZGMEST(:)*1000.) )/XBW )
END WHERE
! 
!-------------------------
!
! compensation point (ppm): temperature response
!
ZGAMMT(:) = PGAMM(:)*PQDGAMM(:)**(0.1*(ZTSPC(:)-25.0))
!
! specific humidity deficit (kg kg-1)
!
ZDSP(:)   = MAX( 0.0, PQSAT(:) - PQA(:) - PLE(:)*PRA(:)/(PRHOA*XLVTT) )
!
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
ZGMEST(:) = ( ZGMEST(:)*PQDGMES(:)**(0.1*(ZTSPC(:)-25.0)) )    &
               /( (1.0+EXP(0.3*(PT1GMES(:)-ZTSPC(:))))*          &
                  (1.0+EXP(0.3*(ZTSPC(:)-PT2GMES(:)))) )  
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
   CALL CCETR(ZXIA,ZIA,ZXMUS,ZABC,ZLAI)
!
!  PAR at different Gauss  level in micmolphot/m2/s
!
   PIACAN(:,JINT)= ZXIA(:)
!
!  Compute conductance and assimilation of CO2: 
!
   CALL COTWO(ZAN0, ZGS0, PGC, PCSP,                         &
                ZDSP, ZDMAXT, ZXIA, PF2, ZGAMMT,               &
                ZFZEROT, ZGMEST, ZEPSO, ZANMAX, ZRDK           )  
 !
! kgCO2/kgair m/s
   ZTAN(:) = ZTAN(:) + ZAN0(:)*PPOI(JINT) 
   ZTGS(:) = ZTGS(:) + ZGS0(:)*PPOI(JINT)
!
END DO

!
!
! Total assimilation
!
PANF(:)= ZTAN(:)
!
! Net assimilation over canopy
!
PAN(:) = (1.0-PDELTA(:))*(1.0-PPSNV(:))*PANF(:)*ZLAI(:)
!
! Dark respiration over canopy (does not depend on radiation, 
! no need to integrate over vertical dimension)
!
PRESP_LEAF(:) = (1.0-PDELTA(:))*(1.0-PPSNV(:))*ZRDK(:)*ZLAI(:)
!
! Gross primary production over canopy
!
PGPP(:) = PAN(:) + PRESP_LEAF(:)
!
! Cumulated daily net assimilation over canopy (kgCO2/m2/day)
!
PANDAY(:) = PANDAY(:) + PAN(:) * PTSTEP * PRHOA(:)
!
! Adjust maximum leaf assimilation:
!
PANFM(:) = MAX( PANF(:), PANFM(:) )
!
! Total conductance over canopy 
!
ZXTGS(:) = ZTGS(:)*ZLAI(:)
!
! Canopy resistance from Ags:
!
PRS(:) = MIN( 1.0/(ZXTGS(:)+ZDENOM_MIN), XRS_MAX)
PRS(:) = MAX( PRS(:), ZRS_MIN)
IF (LHOOK) CALL DR_HOOK('COTWORESTRESS',1,ZHOOK_HANDLE)
!
END SUBROUTINE COTWORESTRESS
