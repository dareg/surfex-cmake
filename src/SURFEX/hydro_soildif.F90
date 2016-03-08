!     #########
      SUBROUTINE HYDRO_SOILDIF(IO, IP, INI, IMX, IR, PTSTEP, PPG, PLETR, PLEG, PEVAPCOR, &
                               PF2WGHT, PPS, PQSAT, PQSATI, PDRAIN, PHORTON, KMAX_LAYER, PQSB)
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
!!      Modif       04/2010  B.Decharme: geometric mean for interfacial conductivity
!!                                       Brook and Corey 
!!      Modif       08/2011  B.Decharme: Optimization using global loops
!!                  10/12    B.Decharme: EVAPCOR snow correction in DIF
!!                  04/13    B.Decharme: Subsurface runoff if SGH (DIF option only)
!!                  07/2013  B.Decharme: Surface / Water table depth coupling
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODD_SURF_PAR, ONLY : XUNDEF, NUNDEF
USE MODD_CSTS,     ONLY : XRHOLW
USE MODD_ISBA_PAR, ONLY : XWGMIN
!
USE MODE_HYDRO_DIF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMX
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
REAL, INTENT(IN)                    :: PTSTEP ! Model time step (s)
!
REAL, DIMENSION(:), INTENT(IN)      :: PPS, PPG, PLETR, PLEG, PEVAPCOR
!                                      PPS    = surface pressure (Pa)
!                                      PPG    = throughfall rate: 
!                                               rate at which water reaches the surface
!                                               of the soil (from snowmelt, rain, canopy
!                                               drip, etc...) (m/s)
!                                      PLETR  = transpiration rate (m/s)
!                                      PLEG   = bare-soil evaporation rate (m/s)
!                                      PEVAPCOR = correction for any excess evaporation 
!                                                from snow as it completely ablates (m/s)
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PQSAT,PQSATI
!                                      specific humidity at saturation
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PF2WGHT
!                                      PF2WGHT      = root-zone transpiration weights (-)
!
INTEGER,               INTENT(IN)   :: KMAX_LAYER  
!                                      KMAX_LAYER = Max number of soil moisture layers (DIF option)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PDRAIN, PHORTON
!                                      PDRAIN   = drainage (flux out of model base) (kg m-2 s-1)
!
REAL, DIMENSION(:),   INTENT(OUT)   :: PQSB     !Lateral subsurface flow [kg/m²/s]
!
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ, JL    ! loop control
!
INTEGER                             :: INJ, INL, IDEPTH ! Number of point and grid layers
!
REAL, DIMENSION(SIZE(IP%XDZG(:,:,1),1))       :: ZINFILTMAX, ZINFILTC, ZEXCESS, ZDGN, ZWGTOT, ZPSIWTD, ZWTD
!                                      ZINFILTMAX = maximum allowable infiltration rate
!                                                   (from Darcy's Law) (m s-1)
!                                      ZEXCESS    = working variable: excess soil water
!                                                   which is used as a constraint 
!                                      ZDGN       = Depth of the last node (m)
!                                                   and the water table (m s-1)
!                                      ZWGTOT    = total soil moisture for ZINFNEG computation
!                                      ZPSIWTD   = matric potential at saturation for water table depth coupling
!                                      ZWTD      = water table depth positive below soil surface (m)
!
REAL, DIMENSION(SIZE(IP%XDZG(:,:,1),1),SIZE(IP%XDZG(:,:,1),2)) :: ZWFLUX, ZDFLUXDT1, ZDFLUXDT2, ZWFLUXN
!                                      ZWFLUX    = vertical soil water flux (+ up) (m s-1)
!                                      ZDFLUXDT  = total vertical flux derrivative
!                                      ZDFLUXDT1 = vertical flux derrivative: dF_j/dw_j
!                                      ZDFLUXDT2 = vertical flux derrivative: dF_j/dw_j+1
!                                      ZWFLUXN   = vertical soil water flux at end of time 
!
REAL, DIMENSION(SIZE(IP%XDZG(:,:,1),1),SIZE(IP%XDZG(:,:,1),2)) :: ZPSI, ZK, ZNU, ZWSAT, ZHEAD, &
                                              ZVAPCOND, ZFRZ, ZKI,         &
                                              ZCAPACITY, ZINFNEG
!                                      ZNU       = interfacial total conductivity (m s-1)
!                                                  at level z_j
!                                      ZK        = hydraulic conductivity (m s-1)
!                                      ZHEAD     = matric potential gradient (-)
!                                      ZWSAT     = ice modified soil porosity (m3 m-3)
!                                      ZFRZ      = diffusion coefficient for freezing (-)
!                                      ZVAPCOND  = vapor conductivity (m s-1) 
!                                      ZKI       = interfacial hydraulic conductivity (m s-1) at level z_j
!                                      ZCAPACITY = simple volumetric water holding capacity estimate for
!                                                  wetting front penetration (-) 
!                                      ZINFNEG   = Negative infiltration (m s-1)
!
REAL, DIMENSION(SIZE(IP%XDZG(:,:,1),1),SIZE(IP%XDZG(:,:,1),2)) :: ZAMTRX, ZBMTRX, ZCMTRX, ZFRC, ZSOL, &
                                              ZTOPQS
!                                      ZAMTRX    = leftmost diagonal element of tri-diagonal
!                                                  coefficient matrix 
!                                      ZBMTRX    = center diagonal element (vector)
!                                      ZCMTRX    = rightmost diagonal element (vector)
!                                      ZFRC      = forcing function (vector)
!                                      ZSOL      = solution vector
!
REAL, DIMENSION(SIZE(IP%XDZG(:,:,1),1),SIZE(IP%XDZG(:,:,1),2))  :: ZINFLAYER
!
REAL, PARAMETER                     :: ZWGHT = 0.5  ! time scheme weight for calculating flux.
!                                                     varies from 0 (explicit time scheme)
!                                                     to 1 (backward difference implicit)
!                                                     Default is 1/2 (Crank-Nicholson)
!
REAL, PARAMETER                     :: ZEICE = 6.0  ! Ice vertical diffusion impedence factor 
!
REAL                                :: ZLOG10, ZS, ZLOG, ZWDRAIN
!
REAL                                :: ZDKDT1, ZDKDT2, ZDHEADDT1, ZDHEADDT2
!                                      ZDKDT1    = hydraulic conductivity derrivative w/r/t upper layer water content
!                                      ZDKDT2    = "" lower layer water content
!                                      ZDHEADDT1 = matric potential gradient derrivative w/r/t upper layer water content
!                                      ZDHEADDT2 = "" lower layer water content
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialization:
!    ---------------
!
IF (LHOOK) CALL DR_HOOK('HYDRO_SOILDIF',0,ZHOOK_HANDLE)
!
INJ = SIZE(IP%XDZG(:,:,1),1)
INL = KMAX_LAYER
!
ZLOG10 = LOG(10.0)
!
PDRAIN    (:) = 0.0
PQSB      (:) = 0.0
PHORTON   (:) = 0.0
ZINFILTC  (:) = 0.0
ZEXCESS   (:) = 0.0
ZINFILTMAX(:) = 0.0
!
ZDGN      (:) = XUNDEF
ZPSIWTD   (:) = XUNDEF
ZWTD      (:) = XUNDEF
!
ZINFNEG  (:,:) = 0.0
ZINFLAYER(:,:) = 0.0
ZFRC     (:,:) = 0.0
ZFRZ     (:,:) = 0.0
!
ZWSAT    (:,:) = XUNDEF
ZCAPACITY(:,:) = XUNDEF
ZPSI     (:,:) = XUNDEF
ZK       (:,:) = XUNDEF
ZVAPCOND (:,:) = XUNDEF
ZKI      (:,:) = XUNDEF
ZSOL     (:,:) = XUNDEF
ZNU      (:,:) = XUNDEF
ZHEAD    (:,:) = XUNDEF
ZWFLUX   (:,:) = XUNDEF
ZWFLUXN  (:,:) = XUNDEF
ZDFLUXDT1(:,:) = XUNDEF
ZDFLUXDT2(:,:) = XUNDEF
ZAMTRX   (:,:) = XUNDEF
ZBMTRX   (:,:) = XUNDEF
ZCMTRX   (:,:) = XUNDEF
!
! Modification/addition of frozen soil parameters
! -----------------------------------------------
!
DO JL=1,INL
  DO JJ=1,INJ
!
    IDEPTH=IMX%NWG_LAYER(JJ,1)
    IF(JL<=IDEPTH)THEN
!
!     Modify soil porosity as ice assumed to become part
!     of solid soil matrix (with respect to liquid flow):
      ZWSAT (JJ,JL) = MAX(XWGMIN, IP%XWSAT(JJ,JL)-IR%XWGI(JJ,JL,1))   
!
!     Factor from (Johnsson and Lundin 1991), except here it is normalized so that it
!     goes to zero in the limit as all available pore space is filled up with ice.
!     For now, a simple constant is used for all soils. Further modifications
!     will be made as research warrents.
!     Old : 10.**(-ZEICE*(IR%XWGI(JJ,JL,1)/(IR%XWGI(JJ,JL,1)+IR%XWG(JJ,JL,1))))
      ZFRZ(JJ,JL) = EXP(ZLOG10*(-ZEICE*(IR%XWGI(JJ,JL,1)/(IR%XWGI(JJ,JL,1)+IR%XWG(JJ,JL,1)))))
!
    ENDIF
!
  ENDDO
ENDDO
!
! Lateral sub-surface flow (m s-1) if Topmodel
! --------------------------------------------
!
IF(IO%CRUNOFF=='SGH')THEN
  ZTOPQS(:,:)=ZFRZ(:,:)*INI%XTOPQS(:,:,1)
ELSE
  ZTOPQS(:,:)=0.0
ENDIF
!
! 1. Infiltration at "t"
!    -------------------
!
! Surface flux term (m s-1): Infiltration (and surface runoff)
! Surface fluxes are limited a Green-Ampt approximation from Abramopoulos et al
! (1988) and Entekhabi and Eagleson (1989).
! Note : when Horton option is used, infiltration already calculated in hydro_sgh
IF(IO%CHORT/='SGH')THEN
!  Green-Ampt approximation for maximum infiltration (derived form)
  ZINFILTMAX(:) = INFMAX_FUNC(IR%XWG(:,:,1), ZWSAT, ZFRZ, IP%XCONDSAT(:,:,1), IP%XMPOTSAT, &
                               IP%XBCOEF, IP%XDZG(:,:,1), IMX%XDG(:,:,1), IO%NLAYER_HORT)
!  Fast(temporal)-response runoff (surface excess) (m s-1):
  PHORTON   (:) = (1.-INI%XFSAT(:)) * MAX(0.0,PPG(:)-ZINFILTMAX(:))
ENDIF
!
!
! 2. Initialise soil moisture profile according to infiltration terms at "t"
!    ----------------------------------------------------------------------
!
!Surface cumulative infiltration  (m)
ZINFILTC(:) = MAX(0.0,PPG(:)-PHORTON(:))*PTSTEP
!
DO JL=1,INL
  DO JJ=1,INJ
    IDEPTH=IMX%NWG_LAYER(JJ,1)
    IF(JL<=IDEPTH)THEN
!     Simple volumetric water holding capacity estimate for wetting front penetration
      ZCAPACITY(JJ,JL) = MAX(0.0,ZWSAT(JJ,JL)-IR%XWG(JJ,JL,1))*IP%XDZG(JJ,JL,1)
!     Infiltration terms (m) :
      ZINFLAYER(JJ,JL) = MIN(ZINFILTC(JJ),ZCAPACITY(JJ,JL))
!     Soil moisture (m3/m3) :
      IR%XWG(JJ,JL,1) = IR%XWG(JJ,JL,1)+ZINFLAYER(JJ,JL)/IP%XDZG(JJ,JL,1)
!     Put remainding infiltration into the next layer (m)
      ZINFILTC(JJ) = ZINFILTC(JJ) - ZINFLAYER(JJ,JL)
    ENDIF
  ENDDO
ENDDO
!
! 3. Compute Fast(temporal)-response runoff and Possible negative infiltration
!    -------------------------------------------------------------------------
!
!Possible negative infiltration  (m s-1)
ZWGTOT(:)=0.0 
DO JL=1,INL
  DO JJ=1,INJ
    IDEPTH=IMX%NWG_LAYER(JJ,1)
    IF(JL<IDEPTH)THEN
      ZINFNEG(JJ,JL) = (MIN(0.0,PPG(JJ))-PEVAPCOR(JJ))*IP%XDZG(JJ,JL,1)*IR%XWG(JJ,JL,1)
      ZWGTOT (JJ   ) = ZWGTOT(JJ)+IP%XDZG(JJ,JL,1)*IR%XWG(JJ,JL,1)
    ENDIF
  ENDDO
ENDDO
DO JL=1,INL
  DO JJ=1,INJ
    ZINFNEG(JJ,JL) = ZINFNEG(JJ,JL)/ZWGTOT(JJ)
  ENDDO
ENDDO 
!
!Fast(temporal)-response runoff (surface excess) (kg m2 s-1):
!special case : if infiltration > total soil capacity
PHORTON(:)=(PHORTON(:)+ZINFILTC(:)/PTSTEP)*XRHOLW
!
!
! 4. Initialise matric potential and hydraulic conductivity at "t"
!    -------------------------------------------------------------
!
DO JL=1,INL
  DO JJ=1,INJ    
    IDEPTH=IMX%NWG_LAYER(JJ,1)
    IF(JL<=IDEPTH)THEN
!     Matric potential (m) :
!     psi=mpotsat*(w/wsat)**(-bcoef)
      ZS          = MIN(1.0,IR%XWG(JJ,JL,1)/ZWSAT(JJ,JL))
      ZLOG        = IP%XBCOEF(JJ,JL)*LOG(ZS)
      ZPSI(JJ,JL) = IP%XMPOTSAT(JJ,JL)*EXP(-ZLOG)
!     Hydraulic conductivity from matric potential (m s-1):
!     k=frz*condsat*(psi/mpotsat)**(-2-3/bcoef)
      ZLOG      = -ZLOG*(2.0+3.0/IP%XBCOEF(JJ,JL))
      ZK(JJ,JL) = ZFRZ(JJ,JL)*IP%XCONDSAT(JJ,JL,1)*EXP(-ZLOG)
    ENDIF
  ENDDO
ENDDO    
!
! Prepare water table depth coupling
! ----------------------------------
!
DO JJ=1,INJ
  IDEPTH=IMX%NWG_LAYER(JJ,1)   
! Depth of the last node
  ZDGN   (JJ) = 0.5*(IMX%XDG(JJ,IDEPTH,1)+IMX%XDG(JJ,IDEPTH-1,1))
  ZPSIWTD(JJ) = ZPSI(JJ,IDEPTH)
  IF(IP%XWTD(JJ)/=XUNDEF)THEN  
!   Water table depth
    ZWTD(JJ)    = MAX(IMX%XDG(JJ,IDEPTH,1),-IP%XWTD(JJ))
!   Modify matric potential at saturation for water table coupling
    ZS          = MIN(1.0,ZWSAT(JJ,IDEPTH)/IP%XWSAT(JJ,IDEPTH))
    ZLOG        = IP%XBCOEF(JJ,IDEPTH)*LOG(ZS)
    ZPSIWTD(JJ) = IP%XMPOTSAT(JJ,IDEPTH)*EXP(-ZLOG)
  ENDIF
ENDDO
!
! 5. Vapor diffusion conductivity (m s-1)
!    ------------------------------------
!
ZVAPCOND(:,:) = VAPCONDCF(IR%XTG(:,:,1), PPS, IR%XWG(:,:,1), IR%XWGI(:,:,1), ZPSI,&
                          IP%XWSAT, IP%XWFC, PQSAT, PQSATI, IMX%NWG_LAYER(:,1), INL)
ZVAPCOND(:,:) = ZFRZ(:,:)*ZVAPCOND(:,:)
!
! 6. Linearized water flux: values at "t"
!    ------------------------------------
!    calculate flux at the beginning of the time step:
!
DO JL=1,INL
  DO JJ=1,INJ
!
    IDEPTH=IMX%NWG_LAYER(JJ,1)
    IF(JL<IDEPTH)THEN
!
!     Total interfacial conductivity (m s-1) And Potential gradient (dimensionless):
      ZKI  (JJ,JL) = SQRT(ZK(JJ,JL)*ZK(JJ,JL+1))
      ZNU  (JJ,JL) = ZKI(JJ,JL) + SQRT(ZVAPCOND(JJ,JL)*ZVAPCOND(JJ,JL+1))
      ZHEAD(JJ,JL) = (ZPSI(JJ,JL)-ZPSI(JJ,JL+1))/IP%XDZDIF(JJ,JL,1)
!
!     Total Sub-surface soil water fluxes (m s-1): (+ up, - down) using Darcy's
!     Law with an added linear drainage term:
      ZWFLUX(JJ,JL) = -ZNU(JJ,JL) * ZHEAD(JJ,JL) - ZKI(JJ,JL)
!
    ELSEIF(JL==IDEPTH)THEN !Last layers   
!        
!     Total interfacial conductivity (m s-1) And Potential gradient (dimensionless):
      ZKI  (JJ,IDEPTH) = ZK(JJ,IDEPTH)
      ZNU  (JJ,IDEPTH) = ZK(JJ,IDEPTH) * IP%XFWTD(JJ)
      ZHEAD(JJ,IDEPTH) = (ZPSI(JJ,IDEPTH)-ZPSIWTD(JJ))/(ZWTD(JJ)-ZDGN(JJ))
!
!     Total Sub-surface soil water fluxes (m s-1): (+ up, - down) using Darcy's
!     Law with an added linear drainage term:
      ZWFLUX(JJ,IDEPTH) = -ZNU(JJ,IDEPTH) * ZHEAD(JJ,IDEPTH) - ZKI(JJ,IDEPTH)
!
    ENDIF
!
  ENDDO
ENDDO
!
!
! 7. Linearized water flux: values at "t+dt"
!    ---------------------------------------
! Flux Derrivative terms, see A. Boone thesis (Annexe E).
!
DO JL=1,INL
  DO JJ=1,INJ
    IDEPTH=IMX%NWG_LAYER(JJ,1)        
    IF(JL<IDEPTH)THEN                
      ZDHEADDT1 = -IP%XBCOEF(JJ,JL  )*ZPSI(JJ,JL  )/(IR%XWG(JJ,JL  ,1)*IP%XDZDIF(JJ,JL,1))
      ZDHEADDT2 = -IP%XBCOEF(JJ,JL+1)*ZPSI(JJ,JL+1)/(IR%XWG(JJ,JL+1,1)*IP%XDZDIF(JJ,JL,1))
      ZDKDT1    = (2.*IP%XBCOEF(JJ,JL  )+3.)*ZKI(JJ,JL)/(2.0*IR%XWG(JJ,JL  ,1))
      ZDKDT2    = (2.*IP%XBCOEF(JJ,JL+1)+3.)*ZKI(JJ,JL)/(2.0*IR%XWG(JJ,JL+1,1))
!     Total Flux derrivative terms:
      ZDFLUXDT1(JJ,JL) = -ZDKDT1*ZHEAD(JJ,JL) - ZNU(JJ,JL)*ZDHEADDT1 - ZDKDT1
      ZDFLUXDT2(JJ,JL) = -ZDKDT2*ZHEAD(JJ,JL) + ZNU(JJ,JL)*ZDHEADDT2 - ZDKDT2  
    ELSEIF(JL==IDEPTH)THEN !Last layers
      ZDHEADDT1 = -IP%XBCOEF(JJ,IDEPTH)*ZPSI   (JJ,IDEPTH)/(IR%XWG  (JJ,IDEPTH,1)*(ZWTD(JJ)-ZDGN(JJ))) &
                  +IP%XBCOEF(JJ,IDEPTH)*ZPSIWTD(JJ       )/(ZWSAT(JJ,IDEPTH)*(ZWTD(JJ)-ZDGN(JJ)))
      ZDHEADDT2 = 0.0
      ZDKDT1    = (2.*IP%XBCOEF(JJ,IDEPTH)+3.)*ZK(JJ,IDEPTH)/IR%XWG(JJ,IDEPTH,1)
      ZDKDT2    = 0.0                
!     Total Flux derrivative terms:
      ZDFLUXDT1(JJ,IDEPTH) = -ZDKDT1*ZHEAD(JJ,IDEPTH)*IP%XFWTD(JJ) - ZNU(JJ,IDEPTH)*ZDHEADDT1 - ZDKDT1
      ZDFLUXDT2(JJ,IDEPTH) = -ZDKDT2*ZHEAD(JJ,IDEPTH)*IP%XFWTD(JJ) + ZNU(JJ,IDEPTH)*ZDHEADDT2 - ZDKDT2  
    ENDIF
  ENDDO
ENDDO
!
! 8. Jacobian Matrix coefficients and Forcing function
!    -------------------------------------------------
!     
!surface layer:
ZFRC  (:,1) = ZWFLUX(:,1) - PLEG(:) - PF2WGHT(:,1)*PLETR(:) + ZINFNEG(:,1) - ZTOPQS(:,1)
ZAMTRX(:,1) = 0.0
ZBMTRX(:,1) = (IP%XDZG(:,1,1)/PTSTEP) - ZWGHT*ZDFLUXDT1(:,1)
ZCMTRX(:,1) = -ZWGHT*ZDFLUXDT2(:,1)
!
!Other sub-surface layers:       
DO JL=2,INL
  DO JJ=1,INJ   
    IDEPTH=IMX%NWG_LAYER(JJ,1)
    IF(JL<=IDEPTH)THEN
      ZFRC  (JJ,JL) = ZWFLUX (JJ,JL) - ZWFLUX(JJ,JL-1) - PF2WGHT(JJ,JL)*PLETR(JJ) + ZINFNEG(JJ,JL) - ZTOPQS(JJ,JL)
      ZAMTRX(JJ,JL) = ZWGHT*ZDFLUXDT1(JJ,JL-1)
      ZBMTRX(JJ,JL) = (IP%XDZG(JJ,JL,1)/PTSTEP) - ZWGHT*(ZDFLUXDT1(JJ,JL)-ZDFLUXDT2(JJ,JL-1))       
      ZCMTRX(JJ,JL) = -ZWGHT*ZDFLUXDT2(JJ,JL)
    ENDIF
  ENDDO
ENDDO
!
!Solve Matrix Equation: tridiagonal system: solve for soil
!water (volumetric water content) tendencies:
!
CALL TRIDIAG_DIF(ZAMTRX,ZBMTRX,ZCMTRX,ZFRC,IMX%NWG_LAYER(:,1),INL,ZSOL)
!
! 9. Final calculations and diagnostics:
!    -----------------------------------
!
!
DO JL=1,INL
  DO JJ=1,INJ
!   
    IDEPTH=IMX%NWG_LAYER(JJ,1)
    IF(JL<IDEPTH)THEN
! 
!     Update liquid water content (m3 m-3):
      IR%XWG(JJ,JL,1)   = IR%XWG(JJ,JL,1) + ZSOL(JJ,JL)    
!
!     Supersaturated drainage (kg m-2 s-1):
      ZEXCESS(JJ)  = MAX(0.0, IR%XWG(JJ,JL,1) - ZWSAT(JJ,JL))
      IR%XWG(JJ,JL  ,1) = MIN(IR%XWG(JJ,JL,1),ZWSAT(JJ,JL))
      IR%XWG(JJ,JL+1,1) = IR%XWG(JJ,JL+1,1) + ZEXCESS(JJ)*(IP%XDZG(JJ,JL,1)/IP%XDZG(JJ,JL+1,1))
!
!     final fluxes (at end of time step) (m s-1):
      ZWFLUXN(JJ,JL) = ZWFLUX(JJ,JL) + ZDFLUXDT1(JJ,JL)*ZSOL(JJ,JL) + ZDFLUXDT2(JJ,JL)*ZSOL(JJ,JL+1)
!
!     Total topmodel subsurface flow
      PQSB (JJ) = PQSB(JJ) + ZTOPQS(JJ,JL)*XRHOLW
!
    ELSEIF(JL==IDEPTH)THEN
! 
!     Update liquid water content (m3 m-3):
      IR%XWG(JJ,IDEPTH,1)   = IR%XWG(JJ,IDEPTH,1) + ZSOL(JJ,IDEPTH) 
!        
!     Supersaturated drainage (kg m-2 s-1):
      ZEXCESS(JJ)    = MAX(0.0, IR%XWG(JJ,IDEPTH,1) - ZWSAT(JJ,IDEPTH))
      IR%XWG(JJ,IDEPTH,1) = MIN(IR%XWG(JJ,IDEPTH,1),ZWSAT(JJ,IDEPTH))   
      PDRAIN (JJ)    = PDRAIN(JJ) + ZEXCESS(JJ)*IP%XDZG(JJ,IDEPTH,1)*XRHOLW/PTSTEP
!   
!     final fluxes (at end of time step) (m s-1):
      ZWFLUXN(JJ,IDEPTH) = ZWFLUX(JJ,IDEPTH) + ZDFLUXDT1(JJ,IDEPTH)*ZSOL(JJ,IDEPTH)
!   
!     Drainage or baseflow out of bottom of model (slow time response) (kg m-2 s-1):
!     Final fluxes (if needed) over the time step (kg m-2 s-1)
!     would be calculated as (for water budget checks) as F = [ wgt*F(t+dt) + (1.-wgt)*F(t) ]*XRHOLW
      PDRAIN (JJ) = PDRAIN(JJ)-(ZWGHT*ZWFLUXN(JJ,IDEPTH)+(1.-ZWGHT)*ZWFLUX(JJ,IDEPTH))*XRHOLW
!
!     Total topmodel subsurface flow
      PQSB (JJ) = PQSB(JJ) + ZTOPQS(JJ,IDEPTH)*XRHOLW
!
    ENDIF
!
  ENDDO
ENDDO
!
! Possible correction/Constraint application: 
!
DO JL=1,INL
  DO JJ=1,INJ   
    IDEPTH=IMX%NWG_LAYER(JJ,1)
    IF(JL<IDEPTH)THEN
!     if the soil water happens to fall below the minimum, then
!     extract needed water from the layer below: this should
!     generally be non-existant: but added to ensure conservation
!     even for the most extreme events.              
      ZEXCESS(JJ)  = MAX(0., XWGMIN  - IR%XWG(JJ,JL,1))
      IR%XWG(JJ,JL,1)   = IR%XWG(JJ,JL,1)   + ZEXCESS(JJ) 
      IR%XWG(JJ,JL+1,1) = IR%XWG(JJ,JL+1,1) - ZEXCESS(JJ)*IP%XDZG(JJ,JL,1)/IP%XDZG(JJ,JL+1,1)
    ELSEIF(JL==IDEPTH.AND.IR%XWG(JJ,IDEPTH,1)<XWGMIN)THEN
!     NOTE, negative moisture can arise for *completely* dry/dessicated soils 
!     owing to the above check because vertical fluxes
!     can be *very* small but nonzero. Here correct owing to
!     small numerical drainage.
      PDRAIN(JJ)     = PDRAIN(JJ) + MIN(0.0,IR%XWG(JJ,IDEPTH,1)-XWGMIN)*IP%XDZG(JJ,IDEPTH,1)*XRHOLW/PTSTEP
      IR%XWG(JJ,IDEPTH,1) = XWGMIN
    ENDIF
  ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('HYDRO_SOILDIF',1,ZHOOK_HANDLE)
!
END SUBROUTINE HYDRO_SOILDIF 






