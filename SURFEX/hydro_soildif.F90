!     #########
      SUBROUTINE HYDRO_SOILDIF(HDIF,PTSTEP,                                         &
                               PBCOEF,PWSAT,PCONDSAT,PCONDSAT_EXP,PMPOTSAT,         &
                               PWFC,PDG,PDZG,PPG,PLETR,PLEG,PLEGI,PF2WGHT,PWG,PWGI, &
                               PTG,PPS,PF2,PWDRAIN,PDRAIN,PRUNOFF,HKSAT,PWWILT,     &
                               HHORT,PFSAT,PEXP_DIF,PALPHA,PN,PM,PL,PWRES           )
!     ##########################################################################
!
!
!!****  *HYDRO_SOILDIF*  
!
!!    PURPOSE
!!    -------
!     This subroutine solves the 1-D (z) mass conservation equation
!     (mixed-form Richard's equation: tendency in terms of volumetric
!     water content, gradient in terms of matric potential)
!     for liquid water (using Darcy's Law for the vertical flux)
!     together with the Clapp and Hornberger (1978) simplification to
!     the Brooks and Corey (1966) empirical model for relating matric
!     potential and hydraulic conductivity to soil water content.
!     Any set of parameters can be used (eg. Clapp and Hornberger 1978;
!     Cosby et al. 1984; etc.) Modifications to the equations is also
!     made for the Van Genucten model/relationships. The equations
!     also incorporate vapor transfer for dry soils. The soil porosity
!     is modified in the presence of soil ice. Soil ice content
!     is also updated herein due to changes resulting from sublimation.
!     The layer averaged set of equations are time differenced
!     using an implicit time scheme. 
!     The equations/model *assume* a heterogeneous soil texture profile,
!     however, if the soil properties are homogeneous, the equations
!     collapse into the standard homogeneous approach (i.e. give the
!     same results as). The eqs are solved rapidly by taking advantage of the
!     fact that the matrix is tridiagonal. 
!     Note that the exponential profile of hydraulic conductivity with soil
!     depth is applied to interfacial conductivity (or interblock)
!
!     
!!**  METHOD
!!    ------
!
!     Direct calculation
!
!!    EXTERNAL
!!    --------
!
!     None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    USE MODD_CST
!!    USE MODI_TRIDIAG_GROUND
!!      
!!    REFERENCE
!!    ---------
!!
!!    Boone (2000)
!!    Boone et al. (2000)
!!      
!!    AUTHOR
!!    ------
!!	A. Boone          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    16/02/00 Boone
!!      Modif       04/2010  B.Decharme, geometric mean for interfacial conductivity
!!                                       Brook and Corey or Van Genuchten 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,     ONLY : XLSTT, XRHOLW, XLVTT
USE MODD_ISBA_PAR, ONLY : XWGMIN
!
USE MODE_HYDRO_DIF
USE MODI_TRIDIAG_GROUND
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
CHARACTER(LEN=*),  INTENT(IN)       :: HDIF   ! NOTE: Only used when HISBA = DIF
!                                             ! 'BC' = Brook and Corey
!                                             ! 'VG' = Van Genuchten
!
REAL, INTENT(IN)                    :: PTSTEP ! Model time step (s)
!
REAL, DIMENSION(:), INTENT(IN)      :: PPS, PF2, PPG, PLETR, PLEG, PLEGI, PWDRAIN, PFSAT
!                                      PPS    = surface pressure (Pa)
!                                      PF2    = total water stress factor (-)
!                                      PPG    = throughfall rate: 
!                                               rate at which water reaches the surface
!                                               of the soil (from snowmelt, rain, canopy
!                                               drip, etc...) (kg m-2 s-1)
!                                      PLETR  = transpiration rate (W m-2)
!                                      PLEG   = bare-soil evaporation rate (W m-2)
!                                      PLEGI  = surface layer sublimation rate (W m-2)
!                                      PWDRAIN= minimum Wg for drainage (m3 m-3)
!                                      PFSAT  = Saturated fraction
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PTG, PDG, PDZG
!                                      PTG = layer-average soil temperatures (K)
!                                      PDG = soil layer depth       (m)
!                                      PDZG= soil layer thicknesses (m)
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PWSAT, PCONDSAT, PCONDSAT_EXP, &
                                       PF2WGHT, PWFC, PEXP_DIF
!                                      PWSAT        = porosity profile (m3 m-3)
!                                      PCONDSAT     = hydraulic conductivity at saturation (m s-1)
!                                      PCONDSAT_EXP = exponential hydraulic conductivity at saturation (m s-1)
!                                      PF2WGHT      = root-zone transpiration weights (-)
!                                      PWFC         = field capacity water content (m3 m-3)
!                                      PEXP_DIF     = exponential coef for interfacial hydraulic conductivity
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PBCOEF,PMPOTSAT
!                                      PMPOTSAT = matric potential at saturation (m) (BC parameters)
!                                      PBCOEF   = slope of the retention curve (-) (BC parameters)
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PALPHA, PN, PM, PL, PWRES
!                                      Van Genuchten parameter
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWG, PWGI
!                                      PWG  = volumetric liquid water content (m3 m-3) 
!                                      PWGI = volumetric ice content (m3 m-3)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PDRAIN, PRUNOFF
!                                      PDRAIN   = drainage (flux out of model base) (kg m-2 s-1)
!                                      PRUNOFF  = runoff (due to saturation (lateral) (kg m-2 s-1)
!
CHARACTER(LEN=*),     INTENT(IN)    :: HHORT    ! Hortonian runoff
!
CHARACTER(LEN=*),     INTENT(IN)    :: HKSAT    ! soil hydraulic profil option
!                                               ! 'DEF'  = ISBA homogenous soil
!                                               ! 'SGH'  = ksat exponential decay
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PWWILT
!                                      PWWILT = wilting point volumetric water
!                                             content (m3 m-3)
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ        ! loop control
!
INTEGER                             :: INLVLD    ! Number of grid layers
!
REAL, DIMENSION(SIZE(PDZG,1))       :: ZINFILTMAX, ZINFILT, ZEXCESS, ZLETR
!                                      ZINFILTMAX = maximum allowable infiltration rate
!                                                   (from Darcy's Law) (m s-1)
!                                      ZEXCESS    = working variable: excess soil water
!                                                   which is used as a constraint 
!                                      ZLETR      = adjusted transpiration (normalized) (W m-2)
!                                      
!
REAL, DIMENSION(SIZE(PDZG,1),SIZE(PDZG,2)) :: ZWFLUX, ZDFLUXDT1, ZDFLUXDT2, ZWFLUXN
!                                      ZWFLUX    = vertical soil water flux (+ up) (m s-1)
!                                      ZDFLUXDT  = total vertical flux derrivative
!                                      ZDFLUXDT1 = vertical flux derrivative: dF_j/dw_j
!                                      ZDFLUXDT2 = vertical flux derrivative: dF_j/dw_j+1
!                                      ZWFLUXN   = vertical soil water flux at end of time 
!
REAL, DIMENSION(SIZE(PDZG,1))     :: ZWFLUX0, ZDFLUXDT20, ZWFLUXN0
!                                                  As above but at z=0 (surface)
!
REAL, DIMENSION(SIZE(PDZG,1))     :: ZDEPTH, ZDINFILT, ZCAPACITY
!                                    ZDEPTH      = total soil (depth) (m)
!                                    ZDINFILT    = maximum potential infiltration depth during a 
!                                                  time step (m)
!                                    ZCAPACITY   = simple volumetric water holding capacity estimate for
!                                                  wetting front penetration (-) 
!
REAL, DIMENSION(SIZE(PDZG,1),SIZE(PDZG,2)) :: ZPSI, ZK, ZKI, ZNU, ZWSAT, ZWFC, ZHEAD, &
                                              ZVAPCOND, ZFRZ
!                                      ZNU       = interfacial total conductivity (m s-1)
!                                                  at level z_j
!                                      ZKI       = interfacial hydraulic conductivity (m s-1)
!                                                  at level z_j
!                                      ZK        = hydraulic conductivity (m s-1)
!                                      ZHEAD     = matric potential gradient (-)
!                                      ZWSAT     = ice modified soil porosity (m3 m-3)
!                                      ZWFC      = ice modified soil field capacity (m3 m-3)
!                                      ZFRZ      = diffusion coefficient for freezing (-)
!                                      ZVAPCOND  = vapor conductivity (m s-1) 
!
REAL, DIMENSION(SIZE(PDZG,1),SIZE(PDZG,2)) :: ZDKDT1, ZDKDT2, ZDHEADDT1, ZDHEADDT2
!                                      ZDKDT1    = hydraulic conductivity derrivative w/r/t
!                                                  upper layer water content
!                                      ZDKDT2    = "" lower layer water content
!                                      ZDHEADDT1 = matric potential gradient derrivative w/r/t
!                                                  upper layer water content
!                                      ZDHEADDT2 = "" lower layer water content
!                                                  
!
REAL, DIMENSION(SIZE(PDZG,1),SIZE(PDZG,2)) :: ZAMTRX, ZBMTRX, ZCMTRX, ZFRC, ZSOL, &
                                              ZDZDIF, ZSGDRAIN, ZWORK
!                                      ZAMTRX    = leftmost diagonal element of tri-diagonal
!                                                  coefficient matrix 
!                                      ZBMTRX    = center diagonal element (vector)
!                                      ZCMTRX    = rightmost diagonal element (vector)
!                                      ZFRC      = forcing function (vector)
!                                      ZSOL      = solution vector
!                                      ZDZDIF    = distance between the midpoints of 
!                                                  consecutive layers (m)
!                                      ZSGDRAIN  = sub-grid drainge (m s-1): calibration parameter
!                                      ZWORK     = working variable for thickness-weighted averages
!
REAL, DIMENSION(SIZE(PDZG,1),SIZE(PDZG,2))  :: ZWLIM, ZWWILT, ZWRES
!
REAL, PARAMETER                     :: ZWGHT = 0.5  ! time scheme weight for calculating flux.
!                                                     varies from 0 (explicit time scheme)
!                                                     to 1 (backward difference implicit)
!                                                     Default is 1/2 (Crank-Nicholson)
!
REAL, PARAMETER                     :: ZEICE = 6.0  ! Ice vertical diffusion impedence factor 
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
! 0. Initialization:
!    ---------------
!
IF (LHOOK) CALL DR_HOOK('HYDRO_SOILDIF',0,ZHOOK_HANDLE)
INLVLD = SIZE(PDZG(:,:),2)
!
PDRAIN (:) = 0.0
PRUNOFF(:) = 0.0
!
ZFRC(:,:) = 0.0
!
! 1. First adjust surface soil ice content for sublimation
!    -----------------------------------------------------
!
PWGI(:,1)    = PWGI(:,1) - PLEGI(:)*PTSTEP/(XLSTT*XRHOLW*PDZG(:,1))
!
! The remaining code in this block are merely constraints to ensure a highly
! accurate water budget: most of the time this code will not have any
! effect on the soil water profile.
!
! If sublimation causes all of the remaining ice to be extracted, remove
! some of the liquid (a correction): NOTE that latent heating already accounted
! for in sublimation term, so no need to alter soil temperature.
!
ZEXCESS(:)  = MAX(0.0,  - PWGI(:,1))
PWG(:,1)    = PWG(:,1)  - ZEXCESS(:)
PWGI(:,1)   = PWGI(:,1) + ZEXCESS(:)
!
! If sublimation is negative (condensation), make sure ice does not
! exceed maximum possible. If it does, then put excess ice into layer below:
! This correction should rarely if ever cause any ice accumulation in the
! sub-surface layer: this is especially true of deeper layers but it is
! accounted for none-the-less.
!
DO JJ=1,INLVLD-1
   ZEXCESS(:)   = MAX(0.0, PWGI(:,JJ) - (PWSAT(:,JJ)-XWGMIN) )
   PWGI(:,JJ)   = PWGI(:,JJ)   - ZEXCESS(:)
   PWGI(:,JJ+1) = PWGI(:,JJ+1) + ZEXCESS(:)*(PDZG(:,JJ)/PDZG(:,JJ+1))
ENDDO
ZEXCESS(:)      = MAX(0.0, PWGI(:,INLVLD) - (PWSAT(:,INLVLD)-XWGMIN) )
PWGI(:,INLVLD)  = PWGI(:,INLVLD)   - ZEXCESS(:)
PDRAIN(:)       = ZEXCESS(:)*PDZG(:,INLVLD)*XRHOLW/PTSTEP
!
!
! 2. Modification/addition of frozen soil parameters
!    -----------------------------------------------
!
! Modify soil porosity as ice assumed to become part
! of solid soil matrix (with respect to liquid flow):
!
ZWSAT(:,:) = MAX(XWGMIN, PWSAT(:,:)-PWGI(:,:))
ZWFC(:,:)  = PWFC(:,:)*ZWSAT(:,:)/PWSAT(:,:)
ZWWILT(:,:)= PWWILT(:,:)*ZWSAT(:,:)/PWSAT(:,:)
!
IF(HDIF=='VG')THEN
  ZWRES(:,:)=PWRES(:,:)*ZWSAT(:,:)/PWSAT(:,:)
ELSE
  ZWRES(:,:)=0.0
ENDIF
!
IF(HKSAT=='SGH')THEN
  ZWLIM(:,:)=ZWWILT(:,:)
ELSE
  IF(HDIF=='BC')THEN
    ZWLIM(:,:)=XWGMIN
  ELSE
    ZWLIM(:,:)=MAX(XWGMIN,ZWRES(:,:))
  ENDIF
ENDIF
!
! Factor from (Johnsson and Lundin 1991), except here it is normalized so that it
! goes to zero in the limit as all available pore space is filled up with ice.
! For now, a simple constant is used for all soils. Further modifications
! will be made as research warrents.
!
ZFRZ(:,:)   = 10.**(-ZEICE*(PWGI(:,:)/(PWGI(:,:)+PWG(:,:))))
!
! 3.) Linear (in time) sub-grid drainage term (input in m s-1 herein)
!     (very simple way to take into account sub-grid
!      water table effects, pores, etc...):
!      Can be calibrated. Also, make this effect vanish/cease for 
!      extremely dry soil layers. This term/option is OFF if WDRAIN = 0.
!
DO JJ=1,INLVLD
   ZSGDRAIN(:,JJ) = PWDRAIN(:)
ENDDO
!
ZSGDRAIN(:,:)=WDRAIN_FUNC(HDIF,PWG,ZWFC,ZWSAT,ZWLIM,ZSGDRAIN,ZFRZ,PCONDSAT_EXP,PBCOEF,PM,PL,ZWRES)
!
! 4. Linearized water flux: values at "t"
!    ------------------------------------
!     calculate flux at the beginning of the time step:
!     first some preliminary computations:
!
! Distance between consecuative layer mid-points (m):
!
DO JJ=1,INLVLD-1
   ZDZDIF(:,JJ)  = 0.5*(PDZG(:,JJ)+PDZG(:,JJ+1))
ENDDO
ZDZDIF(:,INLVLD) = 0.5*PDZG(:,INLVLD)
!
! Matric potential (m) :
!
ZPSI(:,:) = PSI_FUNC(HDIF,PWG,ZWSAT,PMPOTSAT,PBCOEF,PALPHA,PN,PM,ZWRES)
!
! Hydraulic conductivity (m s-1):
!
ZK(:,:) = K_FUNC(HDIF,ZPSI,ZFRZ,PCONDSAT,PMPOTSAT,PBCOEF,PALPHA,PN,PM,PL)
!
! Vapor diffusion conductivity (m s-1)
!
DO JJ=1,INLVLD
   ZVAPCOND(:,JJ) = VAPCONDCF(PTG(:,JJ),PWG(:,JJ),PWGI(:,JJ),ZPSI(:,JJ),PPS(:),PWSAT(:,JJ),PWFC(:,JJ))
ENDDO
!
! Total interfacial conductivity (m s-1)
!
DO JJ=1,INLVLD-1
   ZKI(:,JJ) = PEXP_DIF(:,JJ) * SQRT(ZK(:,JJ)* ZK(:,JJ+1))
   ZNU(:,JJ) = ZKI(:,JJ) + SQRT(ZVAPCOND(:,JJ)*ZVAPCOND(:,JJ+1))
ENDDO
ZKI(:,INLVLD) = PEXP_DIF(:,INLVLD) * ZK(:,INLVLD)
ZNU(:,INLVLD) = 0.0
!
! Potential gradient (dimensionless):
!
DO JJ=1,INLVLD-1
   ZHEAD(:,JJ)  = (ZPSI(:,JJ)-ZPSI(:,JJ+1))/ZDZDIF(:,JJ)
ENDDO
ZHEAD(:,INLVLD) = 0.0
!
! Total Sub-surface soil water fluxes (m s-1): (+ up, - down) using Darcy's
! Law with an added linear drainage term:
!
ZWFLUX(:,:)     = -ZNU(:,:) *ZHEAD(:,:) - ZKI(:,:) - ZSGDRAIN(:,:)
!
! 5. Infiltration at "t"
!    -------------------
!
! Surface flux term (m s-1): Infiltration (and surface runoff)
! Surface fluxes are limited a Green-Ampt approximation from Abramopoulos et al
! (1988) and Entekhabi and Eagleson (1989).
! Note: A modification for frozen soils: rainfall or snowmelt over very 
! frozen upper layers can result in too much infiltration/excess.
! As the maximum rate is explicit, there are no Jacobian elements to be added.
! The maximum possible rate is the minimum of this and the saturated hydraulic
! conductivity.
!
!Doit être encore retravailler si VIC ou TOPMODEL (1-FSAT)
IF(HHORT/='SGH')THEN
   !Green-Ampt approximation for maximum infiltration over ~20cm
   ZINFILTMAX(:) = INFMAX_FUNC(HDIF,ZPSI,ZFRZ,PCONDSAT_EXP,PMPOTSAT,PDZG,PDG)
   ZINFILT   (:) = MIN(PPG(:)/XRHOLW,ZINFILTMAX(:))
ELSE
   ! Infiltration already calculated in hydro_sgh when Horton option is used
   ZINFILT(:) = PPG(:)/XRHOLW
ENDIF
!
!Fast(temporal)-response runoff (surface excess) (kg m-2 s-1):
PRUNOFF(:)   = MAX(0.0,PPG(:)/XRHOLW-ZINFILT(:))*XRHOLW
!
! A simple way to test capacity of layers to store water from top-down:
! This infiltration fix is mainly for large time steps
! and high potential infiltration rates...which is a fairly
! rare combination but nonetheless accounted for here.
! Reason: for a large time step and high rain/melt rates, one
! could have a wetting front progressing deep into the soil
! during a single time step.
! This permits a wetting to enter the layer via the source/sink
! fn below the surface layer.
!
ZDINFILT (:) =ZINFILT(:)*PTSTEP
ZCAPACITY(:) =(ZWSAT(:,1)-PWG(:,1))*PDZG(:,1)
ZDEPTH   (:) =PDG(:,1)
!
ZWORK(:,:)=0.0
ZWORK(:,1)=PDG(:,1)
DO JJ=2,INLVLD
   WHERE(ZDINFILT (:)>=ZCAPACITY(:))
         ZWORK (:,JJ)=PDZG(:,JJ)
         ZDEPTH(:   )=PDG (:,JJ)
   ENDWHERE
   ZCAPACITY(:) = ZCAPACITY(:)+(ZWSAT(:,JJ)-PWG(:,JJ))*PDZG(:,JJ)
ENDDO
!
ZWFLUX0(:) = -ZINFILT(:)*PDG(:,1)/ZDEPTH(:)
DO JJ=2,INLVLD
   ZFRC(:,JJ)  = ZINFILT(:) * ZWORK(:,JJ)/ZDEPTH(:)
ENDDO
!
! 6. Linearized water flux: values at "t+dt"
!    ---------------------------------------
! Flux Derrivative terms, see A. Boone thesis (Annexe E):
!
! Matric potential and Hydraulic conductivity:
! (as vapor conductivity has no bearing on numerical stability
! for the time steps considered (i.e. it changes slowly relative to
! the typical time steps used herein), and because its derrivative
! is complex, its derivative is neglected).
!
DO JJ=1,INLVLD-1
!
   ZDHEADDT1(:,JJ) = DPSI_FUNC(HDIF,ZPSI(:,JJ  ),PWG(:,JJ  ),ZWSAT(:,JJ  ),PBCOEF(:,JJ  ),PN(:,JJ  ),PM(:,JJ  ),ZWRES(:,JJ  )) &
                   / ZDZDIF(:,JJ)
   ZDHEADDT2(:,JJ) = DPSI_FUNC(HDIF,ZPSI(:,JJ+1),PWG(:,JJ+1),ZWSAT(:,JJ+1),PBCOEF(:,JJ+1),PN(:,JJ+1),PM(:,JJ+1),ZWRES(:,JJ+1)) &
                   / ZDZDIF(:,JJ)
!   
   ZDKDT1(:,JJ) = PEXP_DIF(:,JJ) * DK_FUNC(HDIF,PWG(:,JJ  ),ZWSAT(:,JJ  ),ZK(:,JJ),ZK(:,JJ+1), &
                                           PBCOEF(:,JJ  ),PM(:,JJ  ),PL(:,JJ  ),ZWRES(:,JJ  )  )
   ZDKDT2(:,JJ) = PEXP_DIF(:,JJ) * DK_FUNC(HDIF,PWG(:,JJ+1),ZWSAT(:,JJ+1),ZK(:,JJ),ZK(:,JJ+1), &
                                           PBCOEF(:,JJ+1),PM(:,JJ+1),PL(:,JJ+1),ZWRES(:,JJ+1)  )
!   
ENDDO
ZDHEADDT1(:,INLVLD) = 0.0
ZDHEADDT2(:,INLVLD) = 0.0
ZDKDT1   (:,INLVLD) = PEXP_DIF(:,INLVLD) * 2.0*DK_FUNC(HDIF,PWG(:,INLVLD),ZWSAT(:,INLVLD),ZK(:,INLVLD),ZK(:,INLVLD), &
                                                       PBCOEF(:,INLVLD),PM(:,INLVLD),PL(:,INLVLD),ZWRES(:,INLVLD)    )
ZDKDT2   (:,INLVLD) = 0.0
!
! Total Flux derrivative terms:
!
! Interior and at model base:
!
ZDFLUXDT1(:,:)     = -ZDKDT1(:,:)*ZHEAD(:,:) - ZNU(:,:)*ZDHEADDT1(:,:) - ZDKDT1(:,:)
ZDFLUXDT2(:,:)     = -ZDKDT2(:,:)*ZHEAD(:,:) + ZNU(:,:)*ZDHEADDT2(:,:) - ZDKDT2(:,:)
!
! At surface:
!
ZDFLUXDT20(:)      = 0.0
!
!
! 7. Jacobian Matrix coefficients
!    ----------------------------
!
ZAMTRX(:,1)      = 0.0
ZBMTRX(:,1)      = (PDZG(:,1)/PTSTEP) - ZWGHT*(ZDFLUXDT1(:,1)-ZDFLUXDT20(:))
ZCMTRX(:,1)      = -ZWGHT*ZDFLUXDT2(:,1)
DO JJ=2,INLVLD-1
   ZAMTRX(:,JJ)  = ZWGHT*ZDFLUXDT1(:,JJ-1)
   ZBMTRX(:,JJ)  = (PDZG(:,JJ)/PTSTEP) - ZWGHT*(ZDFLUXDT1(:,JJ)-ZDFLUXDT2(:,JJ-1)) 
   ZCMTRX(:,JJ)  = -ZWGHT*ZDFLUXDT2(:,JJ)
ENDDO
ZAMTRX(:,INLVLD) = ZWGHT*ZDFLUXDT1(:,INLVLD-1)
ZBMTRX(:,INLVLD) = (PDZG(:,INLVLD)/PTSTEP) - ZWGHT*(ZDFLUXDT1(:,INLVLD)-ZDFLUXDT2(:,INLVLD-1)) 
ZCMTRX(:,INLVLD) = 0.0 
!     
! Forcing function:
! Soil water sink terms: convert from (W m-2) and (kg m-2 s-1) to (m s-1)
!
ZLETR(:)      = PLETR(:)/(PF2(:)*XLVTT*XRHOLW)
!
! surface layer:
!
ZFRC(:,1)     = ZWFLUX(:,1) - ZWFLUX0(:) - PLEG(:)/(XLVTT*XRHOLW) - PF2WGHT(:,1)*ZLETR(:)
!
! sub-surface: vertical distribution of soil water extraction:
!
DO JJ=2,INLVLD
   ZFRC(:,JJ) = ZFRC(:,JJ) + ZWFLUX(:,JJ) - ZWFLUX(:,JJ-1) - PF2WGHT(:,JJ)*ZLETR(:)
ENDDO
!
! 
! Solve Matrix Equation: tridiagonal system: solve for soil
! water (volumetric water content) tendencies:
!
CALL TRIDIAG_GROUND(ZAMTRX,ZBMTRX,ZCMTRX,ZFRC,ZSOL)
!
! Update liquid water content (m3 m-3):
!
PWG(:,:) = PWG(:,:) + ZSOL(:,:)
!
!
! 8. Final calculations and diagnostics:
!    -----------------------------------
!
! supersaturated drainage (kg m-2 s-1):
! Should be improved for ice/permafrost regime
DO JJ=1,INLVLD-1
   ZEXCESS(:)  = MAX(0.0, PWG(:,JJ) - ZWSAT(:,JJ))
   PWG(:,JJ+1) = PWG(:,JJ+1) + ZEXCESS(:)*(PDZG(:,JJ)/PDZG(:,JJ+1))
ENDDO
ZEXCESS(:)     = MAX(0.0, PWG(:,INLVLD) - ZWSAT(:,INLVLD))
PDRAIN(:)      = PDRAIN(:) + ZEXCESS(:)*PDZG(:,INLVLD)*XRHOLW/PTSTEP
!
PWG(:,:) = MIN(PWG(:,:),ZWSAT(:,:))
!
! Outputs: 
!   
! First, final fluxes (at end of time step) (m s-1):
!
ZWFLUXN0(:)       = ZWFLUX0(:) + ZDFLUXDT20(:)*ZSOL(:,1)
DO JJ=1,INLVLD-1
   ZWFLUXN(:,JJ)  = ZWFLUX(:,JJ) + (ZDFLUXDT1(:,JJ)*ZSOL(:,JJ)+ZDFLUXDT2(:,JJ)*ZSOL(:,JJ+1))
ENDDO
ZWFLUXN(:,INLVLD) = ZWFLUX(:,INLVLD) + ZDFLUXDT1(:,INLVLD)*ZSOL(:,INLVLD)
!
! Drainage or baseflow out of bottom of model (slow time response) (kg m-2 s-1):
! Final fluxes (if needed) over the time step (kg m-2 s-1)
! would be calculated as (for water budget checks) as F = [ wgt*F(t+dt) + (1.-wgt)*F(t) ]*XRHOLW
!
PDRAIN(:)    = PDRAIN(:) - (ZWGHT*ZWFLUXN(:,INLVLD) + (1.-ZWGHT)*ZWFLUX(:,INLVLD))*XRHOLW
!
! Possible correction/Constraint application: 
! if the soil water happens to fall below the minimum, then
! extract needed water from the layer below: this should
! generally be non-existant: but added to ensure conservation
! even for the most extreme events.
!
DO JJ=1,INLVLD-1
   ZEXCESS(:)  = MAX(0., XWGMIN  - PWG(:,JJ))
   PWG(:,JJ)   = PWG(:,JJ)   + ZEXCESS(:) 
   PWG(:,JJ+1) = PWG(:,JJ+1) - ZEXCESS(:)*PDZG(:,JJ)/PDZG(:,JJ+1)
ENDDO
!
! NOTE, the following can permit a very small mass loss/conservation error.
! But negative moisture can arise for *completely* dry/dessicated soils 
! owing to the above check because vertical fluxes
! can be *very* small but nonzero. Here correct owing to
! small numerical drainage.
!
WHERE(PWG(:,INLVLD) < XWGMIN)
   PWG(:,INLVLD)   = XWGMIN
   PDRAIN(:)       = 0.0
END WHERE
!
IF (LHOOK) CALL DR_HOOK('HYDRO_SOILDIF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!
!
!####################################################################
!####################################################################
!####################################################################
FUNCTION VAPCONDCF(PTG,PWG,PWGI,PPSIA,PPS,PWSAT,PWFC) RESULT(PVAPCOND)
!
! Uses method of Braud et al. (1993) for
! vapor conductivity (m s-1)
!
USE MODD_CSTS,       ONLY : XMV, XMD, XTT, XP00, XG, XRV, XRHOLW
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODE_THERMOS
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)      :: PTG,PWG,PWGI,PPSIA, &
                                       PPS,PWSAT,PWFC
!
REAL, DIMENSION(SIZE(PPS))          :: PVAPCOND
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PPS))          :: ZDVA, ZFVA, ZCHI, ZHUM, ZPV,   &
                                       ZESAT, ZQSAT, ZESATI,          &
                                       ZQSATI, ZWG
!
! Parameters:
!
REAL, PARAMETER                     :: ZTORTY = 0.66         ! (-)
REAL, PARAMETER                     :: ZNV    = 1.88         ! (-)
REAL, PARAMETER                     :: ZCV    = 2.17e-5      ! (m2/s)
REAL, PARAMETER                     :: ZWK    = 0.05         ! (m3 m-3)
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('HYDRO_SOILDIF:VAPCONDCF',0,ZHOOK_HANDLE)
!
ZDVA(:)     = XUNDEF
ZFVA(:)     = XUNDEF
ZCHI(:)     = XUNDEF
ZHUM(:)     = XUNDEF
ZPV(:)      = XUNDEF
ZQSAT(:)    = XUNDEF
ZQSATI(:)   = XUNDEF
ZESAT(:)    = XUNDEF
ZESATI(:)   = XUNDEF
!
PVAPCOND(:) = 0.0
!
! Only perform this computation if the soil is sufficiently
! dry (as otherwise the hydraulic conductivity dominates
! the diffusion coefficient). Arbitrarily base threshold on field
! capacity water content:
!
ZWG(:) = PWG(:) + PWGI(:)
WHERE(ZWG(:)< PWFC(:) .AND. ZWG(:) > XWGMIN)
!
! Vapor pressure over liquid and solid water surfaces (Pa), respectively:
!
   ZQSAT(:)    = QSAT(PTG(:), PPS(:))
   ZESAT(:)    = ZQSAT(:)* PPS(:)/((XMV/XMD)+ZQSAT(:) *(1.-(XMV/XMD)))
!
   ZQSATI(:)   = QSATI(MIN(XTT,PTG(:)),PPS(:))
   ZESATI(:)   = ZQSATI(:)*PPS(:)/((XMV/XMD)+ZQSATI(:)*(1.-(XMV/XMD)))
!
! molecular diffusivity of water vapor (m2 s-1):
!
   ZDVA(:)     = ZCV*(XP00/PPS(:))*((PTG(:)/XTT)**ZNV)
!
! function of pore space: 
!
   ZFVA(:)     = (PWSAT(:) - ZWG(:))*(1.+(ZWG(:)/(PWSAT(:)-ZWK)))
   ZFVA(:)     = MIN(ZFVA(:),PWSAT(:))
!
! relative humidity of air in soil pores:
!
   ZHUM(:)     = MAX(0.001, EXP(PPSIA(:)*XG/(XRV*PTG(:))) )
!
! fraction of frozen water:
!
   ZCHI(:)     = PWGI(:)/ZWG(:)
!
! vapor pressure within pore space (Pa):
!
   ZPV(:)      = ZHUM(:)*(ZCHI(:)*ZESAT(:) + (1.-ZCHI(:))*ZESATI(:))
!
! vapor conductivity (kg m-2 s-1)
!
   PVAPCOND(:) = ZTORTY*PPS(:)*ZDVA(:)*ZFVA(:)*XG*ZPV(:)/                  &
                 ((PPS(:)-ZPV(:))*(XRV*XRV*PTG(:)*PTG(:)))
!
! vapor conductivity (m s-1)
!
   PVAPCOND(:) = PVAPCOND(:)/XRHOLW
!
END WHERE
!
IF (LHOOK) CALL DR_HOOK('HYDRO_SOILDIF:VAPCONDCF',1,ZHOOK_HANDLE)
!
END FUNCTION VAPCONDCF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE HYDRO_SOILDIF 






