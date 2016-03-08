!     #########
      SUBROUTINE HYDRO_SGH(IO, P, IP, INI, IMX, IR, DGEIP, DGMI, PTSTEP, PPG, PPG_MELT, PDUNNE )  
!
!     #####################################################################
!
!!****  *HYDRO_SGH*  
!!
!!    PURPOSE
!!    =======
!
!     1. Determine the Horton runoff that take account of a spatial subgri
!        exponential distribution of the precipitation and of the surface ksat.
!     1. Determine the surface saturated fraction (dt92 or Topmodel).
!     3. Determine the Dunne runoff (dt92 or Topmodel).
!     4. Determine the infiltration rate.
!     5. Determine the flooplains interception and infiltration rate.
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!        ===================
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
!
USE MODD_CSTS,      ONLY : XRHOLW, XDAY, XCL, XCI, XRHOLI
USE MODD_ISBA_PAR,  ONLY : XWGMIN, XSPHSOIL, XDRYWGHT
USE MODD_SURF_PAR,  ONLY : XUNDEF
USE MODD_SGH_PAR,   ONLY : XHORT_DEPTH
!
#ifdef TOPD
USE MODD_COUPLING_TOPD, ONLY : LCOUPL_TOPD, XAS_NATURE, XATOP, NMASKT_PATCH
#endif
!
USE MODI_HYDRO_DT92
!
USE MODE_HYDRO_DIF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMX
TYPE(ISBA_PGD_t), INTENT(INOUT) :: P
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIP
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
!
REAL, INTENT(IN)                 :: PTSTEP
!                                   timestep of the integration
!
REAL, DIMENSION(:), INTENT(INOUT):: PPG
REAL, DIMENSION(:), INTENT(IN)   :: PPG_MELT
!                                   PPG      = water reaching the ground
!                                   PPG_MELT = snowmelt reaching the ground
!
REAL, DIMENSION(:), INTENT(OUT)  :: PDUNNE
!                                   PDUNNE  = Dunne runoff
!
!*      0.2    declarations of local variables
!
REAL, PARAMETER                            :: ZEICE = 6.0  ! Ice vertical diffusion impedence factor 
!
REAL, DIMENSION(SIZE(PPG))                 :: ZPG_INI, ZFROZEN, ZIMAX_ICE, ZIMAX, &
                                              ZHORT_R, ZHORT_M, ZSOILMAX, ZIF_MAX
!                                             ZFROZEN  = frozen soil fraction for runoff
!                                             ZIMAX_ICE    = maximum infiltration rate for frozen soil
!                                             ZIMAX     = maximum infiltration rate for unfrozen soil
REAL, DIMENSION(SIZE(PPG))                 :: ZWG2_AVG, ZWGI2_AVG, ZWSAT_AVG, ZWWILT_AVG
!                                             Average water and ice content
!                                             values over the soil depth D2 (for calculating surface runoff)
!
REAL, DIMENSION(SIZE(IMX%XDG(:,:,1),1),SIZE(IMX%XDG(:,:,1),2)) :: ZWSAT, ZFRZ
!
REAL, DIMENSION(SIZE(PPG))                 :: ZPG_WORK, ZRUISDT, ZNL_HORT, ZDEPTH
!
REAL, DIMENSION(SIZE(PPG))                 :: ZRUNOFF_TOPD
!
REAL                                       :: ZEFFICE, ZLOG10, ZLOG, ZS, ZD_H
!
REAL                                       :: ZTDIURN, ZSOILHEATCAP
!                                             ZTDIURN      = thermal penetration depth for restore (m)
!                                             ZSOILHEATCAP = Total soil volumetric heat capacity [J/(m3 K)]
!
INTEGER                                    :: INJ, INL, JJ, JL, IDEPTH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!allocate
!
IF (LHOOK) CALL DR_HOOK('HYDRO_SGH',0,ZHOOK_HANDLE)
!
!initialize
!
ZFROZEN  (:)  = 0.0
ZIMAX_ICE(:)  = 0.0
ZIMAX    (:)  = 0.0
!
ZWSAT  (:,:)  = 0.0
!
ZLOG10 = LOG(10.0)
!
!IO%CRUNOFF = DT92 ZFSAT calculation
ZWG2_AVG(:)   = 0.0
ZWGI2_AVG(:)  = 0.0
ZWSAT_AVG(:)  = 0.0
ZWWILT_AVG(:) = 0.0
!
!IO%CHORT=SGH
ZHORT_R(:) = 0.0
ZHORT_M(:) = 0.0
!
!DGEIP%XIFLOOD calculation
ZSOILMAX(:) = 0.0
ZIF_MAX(:)  = 0.0
!
!IO%CRUNOFF = DT92 DUNNE calculation
ZPG_WORK(:)   = 0.0
ZRUISDT(:)    = 0.0
!
!IO%CRUNOFF = TOPD DUNNE calculation
ZRUNOFF_TOPD(:) = 0.0
!
!to limit numerical artifacts
ZPG_INI(:) = PPG(:) + PPG_MELT(:)
!
!
INJ=SIZE(IMX%XDG(:,:,1),1)
INL=MAXVAL(IMX%NWG_LAYER(:,1))
!
!-------------------------------------------------------------------------------
!
!*           1. Surface saturated fraction
!            -----------------------------
!
IF( IO%CRUNOFF=='DT92' .OR. IO%CRUNOFF == 'TOPD' )THEN
!
! Calculate the layer average water content for the sub-grid
! surface runoff computation: use IP%XRUNOFFD(:,1) as the depth over which
! runoff is calculated.
!
! First, determine a weight for each layer's contribution
! to thickness averaged water content and soil properties for runoff.
!
  IF (IO%CISBA == 'DIF') THEN
!
! Vertically averaged soil properties and moisture for surface runoff computation:
!
    DO JL=1,IO%NLAYER_DUN
      DO JJ=1,INJ
        IDEPTH=IMX%NWG_LAYER(JJ,1)
        IF(JL<=IDEPTH)THEN
          ZWG2_AVG  (JJ) = ZWG2_AVG  (JJ) + IP%XSOILWGHT(JJ,JL,1)*IR%XWG (JJ,JL,1)/MAX(1.E-6,IP%XRUNOFFD(JJ,1))
          ZWGI2_AVG (JJ) = ZWGI2_AVG (JJ) + IP%XSOILWGHT(JJ,JL,1)*IR%XWGI(JJ,JL,1)/MAX(1.E-6,IP%XRUNOFFD(JJ,1))
          ZWSAT_AVG (JJ) = ZWSAT_AVG (JJ) + IP%XSOILWGHT(JJ,JL,1)*IP%XWSAT (JJ,JL)/MAX(1.E-6,IP%XRUNOFFD(JJ,1))
          ZWWILT_AVG(JJ) = ZWWILT_AVG(JJ) + IP%XSOILWGHT(JJ,JL,1)*IP%XWWILT(JJ,JL)/MAX(1.E-6,IP%XRUNOFFD(JJ,1))
        ENDIF
      ENDDO
    ENDDO
!
  ELSE
!           
    ZWG2_AVG(:)   = IR%XWG(:,2,1)
    ZWGI2_AVG(:)  = IR%XWGI(:,2,1)
    ZWSAT_AVG(:)  = IP%XWSAT(:,1)
    ZWWILT_AVG(:) = IP%XWWILT(:,1)
!      
  ENDIF
!
  IF(IO%CHORT=='SGH')THEN
    !runoff over frozen soil explicitly calculated
    ZWGI2_AVG(:)=0.0
  ENDIF
!
   DO JJ=1,INJ
     ZS = MIN(1.0,(ZWG2_AVG(JJ)+ZWGI2_AVG(JJ)-ZWWILT_AVG(JJ))/(ZWSAT_AVG(JJ)-ZWWILT_AVG(JJ)))
     INI%XFSAT(JJ) = 1.0-(1.0-MAX(0.0,ZS))**(P%XRUNOFFB(JJ)/(P%XRUNOFFB(JJ)+1.))
   ENDDO        
!
ENDIF
!
!*           2. Horton runoff
!            ----------------
!
IF(IO%CHORT=='SGH'.OR.IO%LFLOOD)THEN  
!
  IF(IO%CISBA == 'DIF')THEN
!
!   no subgrid frozen soil fraction of the grid cells
    ZFROZEN(:) = 0.0
!    
    DO JL=1,IO%NLAYER_HORT
      DO JJ=1,INJ   
!              
!       Modify soil porosity as ice assumed to become part
!       of solid soil matrix (with respect to liquid flow):                
        ZWSAT(JJ,JL) = MAX(XWGMIN, IP%XWSAT(JJ,JL)-IR%XWGI(JJ,JL,1)) 
!        
!       Impedance Factor from (Johnsson and Lundin 1991).
        ZFRZ(JJ,JL) = EXP(ZLOG10*(-ZEICE*(IR%XWGI(JJ,JL,1)/(IR%XWGI(JJ,JL,1)+IR%XWG(JJ,JL,1)))))
!
      ENDDO
    ENDDO    
!
!   Calculate infiltration MAX using green-ampt approximation (derived form)
    ZIMAX(:) = INFMAX_FUNC(IR%XWG(:,:,1), ZWSAT, ZFRZ, IP%XCONDSAT(:,:,1), IP%XMPOTSAT, IP%XBCOEF, &
                           IP%XDZG(:,:,1), IMX%XDG(:,:,1), IO%NLAYER_HORT)
!  
  ELSE
!
    DO JJ=1,INJ
!
!    Total soil volumetric heat capacity [J/(m3 K)]:
!
      ZSOILHEATCAP = XCL*XRHOLW*IR%XWG (JJ,2,1) +                           &
                     XCI*XRHOLI*IR%XWGI(JJ,2,1) +                           &
                     XSPHSOIL*XDRYWGHT*(1.0-IP%XWSAT(JJ,1))
!                     
!     Soil thickness which corresponds to the diurnal surface temperature
!     wave penetration depth as T2 is the average temperature for this layer:
!
      ZTDIURN   = MIN(IMX%XDG(JJ,2,1), 4./(ZSOILHEATCAP*DGMI%XCG(JJ)))
!    
!     Effective frozen depth penetration 
!
      ZEFFICE = IMX%XDG(JJ,2,1)*IR%XWGI(JJ,2,1)/(IR%XWGI(JJ,2,1)+IR%XWG(JJ,2,1))
!
!     Modify soil porosity as ice assumed to become part
!     of solid soil matrix (with respect to liquid flow):
!
      ZWSAT(JJ,1) = MAX(XWGMIN, IP%XWSAT(JJ,1)-IR%XWGI(JJ,2,1)) 
!
!     calculate the subgrid frozen soil fraction of the grid cells
!
      ZFROZEN (JJ) = MIN(1.,ZEFFICE/MAX(IMX%XD_ICE(JJ,1),ZTDIURN))
!
!     Impedance Factor from (Johnsson and Lundin 1991).
!
      ZFRZ(JJ,1) = EXP(ZLOG10*(-ZEICE*MIN(1.,ZEFFICE/ZTDIURN)))
!
!     Calculate infiltration MAX on frozen soil as Johnsson and Lundin (1991).
!     The max infiltration is equal to the unsaturated conductivity function at a
!     water content corresponding to the total porosity less the ice-filled volume.
!
      ZS =MIN(1.,ZWSAT(JJ,1)/IP%XWSAT(JJ,1))
      ZIMAX_ICE(JJ)=ZFRZ(JJ,1)*IP%XKSAT_ICE(JJ,1)*(ZS**(2*IP%XBCOEF(JJ,1)+3.))
!
!     Calculate infiltration MAX on unfrozen soil using green-ampt approximation
!    
      ZS   =MIN(1.,IR%XWG(JJ,2,1)/ZWSAT(JJ,1))
      ZD_H =MIN(0.10,IMX%XDG(JJ,2,1))
      ZIMAX(JJ)=IP%XCONDSAT(JJ,1,1)*(IP%XBCOEF(JJ,1)*IP%XMPOTSAT(JJ,1)*(ZS-1.0)/ZD_H+1.0)
!
    ENDDO
!
  ENDIF
!
ENDIF
!
IF(IO%CHORT=='SGH')THEN
!
! calculate the Horton runoff generated by the rainfall rate
!
  IF(IO%CRAIN=='SGH')THEN
!
    WHERE(PPG(:)>0.)
       ZHORT_R(:) = (1.- ZFROZEN(:))* PPG(:)/((ZIMAX    (:)*XRHOLW*INI%XMUF(:)/PPG(:)) + 1.) & !unfrozen soil
                  +      ZFROZEN(:) * PPG(:)/((ZIMAX_ICE(:)*XRHOLW*INI%XMUF(:)/PPG(:)) + 1.)   !frozen soil
    END WHERE        
!
  ELSE
!
    ZHORT_R(:) = (1.- ZFROZEN(:))* MAX(0.,PPG(:)-ZIMAX    (:)*XRHOLW)          & !unfrozen soil
               +      ZFROZEN(:) * MAX(0.,PPG(:)-ZIMAX_ICE(:)*XRHOLW)            !frozen soil
!
  ENDIF
!
! calculate the Horton runoff generated by the snow melt
!        
  ZHORT_M(:) = (1.- ZFROZEN(:))* MAX(0.,PPG_MELT(:)-ZIMAX    (:)*XRHOLW)          & !unfrozen soil
             +      ZFROZEN(:) * MAX(0.,PPG_MELT(:)-ZIMAX_ICE(:)*XRHOLW)            !frozen soil
!
! calculate the  total Horton runoff 
!
  WHERE(INI%XFFLOOD(:)<=INI%XFSAT(:))
        DGEIP%XHORT(:) = (1. - INI%XFSAT(:))   * (ZHORT_R(:) + ZHORT_M(:))
  ELSEWHERE
        DGEIP%XHORT(:) = (1. - INI%XFFLOOD(:)) * (ZHORT_R(:) + ZHORT_M(:))
  ENDWHERE
!
ELSE
!
  DGEIP%XHORT(:) = 0.0
!
ENDIF
!
! calculate all water reaching the ground
!
PPG  (:) = PPG(:) + PPG_MELT(:)   
!
!
!*           3. Dunne runoff and flood interception
!            --------------------------------------
!
! Interception by the flooplains
!
IF(IO%LFLOOD)THEN
  DGEIP%XPFLOOD(:)=INI%XFFLOOD(:)*MAX(0.0,PPG(:))
ELSE
  DGEIP%XPFLOOD(:)=0.0
ENDIF
!
IF(IO%CRUNOFF=='SGH ')THEN
!        
! calculate the Dunne runoff with TOPMODEL
!
  PDUNNE(:) = MAX(PPG(:),0.0) * MAX(INI%XFSAT(:)-INI%XFFLOOD(:),0.0)
!
ELSEIF (IO%CRUNOFF=='DT92' .OR. IO%CRUNOFF=='TOPD')THEN
!
!*       Dumenil et Todini (1992)  RUNOFF SCHEME
!        ---------------------------------------         
!
! surface runoff done only on the Fsat-Fflood fraction
!
  ZPG_WORK(:) = PPG(:) - DGEIP%XHORT(:) - DGEIP%XPFLOOD(:)
!
#ifdef TOPD
  IF ( LCOUPL_TOPD.AND.IO%CRUNOFF == 'TOPD' )THEN
    !
    DO JJ=1,SIZE(NMASKT_PATCH)
      IF (NMASKT_PATCH(JJ)/=0) THEN
        IF ( XATOP(NMASKT_PATCH(JJ))/=0. .AND. XAS_NATURE(NMASKT_PATCH(JJ))/=XUNDEF ) THEN
          ZRUNOFF_TOPD(JJ) = MAX(PPG(JJ),0.0) * MAX(XAS_NATURE(NMASKT_PATCH(JJ)),0.0)
        ENDIF
      ENDIF 
    ENDDO
    !
  ENDIF
#endif
  !
  CALL HYDRO_DT92(PTSTEP, P%XRUNOFFB, ZWWILT_AVG,  IP%XRUNOFFD(:,1), ZWSAT_AVG, &
                  ZWG2_AVG, ZWGI2_AVG, ZPG_WORK, ZRUISDT           )
!
  PDUNNE(:) = ZRUISDT(:)*IP%XRUNOFFD(:,1)*XRHOLW/PTSTEP
  !
#ifdef TOPD
  IF (LCOUPL_TOPD.AND.IO%CRUNOFF == 'TOPD') THEN
    PDUNNE(:) = ZRUNOFF_TOPD(:) +  PDUNNE(:)*(1-XATOP(NMASKT_PATCH(:)))
  ENDIF
#endif
  !
  IF(IO%LFLOOD)THEN
    WHERE(INI%XFFLOOD(:)>=INI%XFSAT(:).AND.INI%XFFLOOD(:)>0.0)PDUNNE(:) = 0.0
  ENDIF   
  !
ELSE
! 
! Default case (no subgrid runoff)
!
  INI%XFSAT (:) = 0.0
  PDUNNE(:) = 0.0
!
ENDIF
!
! calculate the infiltration rate after runoff
!
PPG  (:) = PPG(:) - PDUNNE(:) - DGEIP%XHORT(:) - DGEIP%XPFLOOD(:)
!
! Supress numerical artifacts:
!
WHERE (ZPG_INI(:)<0.0)
  PPG(:)     = ZPG_INI(:)
  DGEIP%XHORT(:) = 0.0
  PDUNNE (:) = 0.0
  DGEIP%XPFLOOD(:) = 0.0
ENDWHERE
!
!*           4. infiltration rate from floodplains (à revoir pour DF !!!)
!            -------------------------------------
!
IF(IO%LFLOOD)THEN
!
! calculate the maximum flood infiltration
!
  ZIF_MAX(:) = MAX(0.,(1.- ZFROZEN(:))) * ZIMAX    (:)*XRHOLW &   !unfrozen soil
             +             ZFROZEN(:)   * ZIMAX_ICE(:)*XRHOLW     !frozen soil
!
  DGEIP%XIFLOOD(:)=MAX(0.0,(INI%XFFLOOD(:)-INI%XFSAT(:)))*MIN(INI%XPIFLOOD(:)/PTSTEP,ZIF_MAX(:))
!
  IF(IO%CISBA == 'DIF')THEN
    ZDEPTH(:)=0.0
    DO JL=1,IO%NLAYER_HORT
      DO JJ=1,INJ
        IF(ZDEPTH(JJ)<XHORT_DEPTH)THEN
          ZSOILMAX(JJ) = ZSOILMAX(JJ)+MAX(0.0,ZWSAT(JJ,JL)-IR%XWG(JJ,JL,1))*IP%XDZG(JJ,JL,1)*XRHOLW/PTSTEP
          ZDEPTH  (JJ) = IMX%XDG(JJ,JL,1)
        ENDIF
      ENDDO
    ENDDO
  ELSE
    DO JJ=1,INJ
      ZWSAT(JJ,1)  = MAX(XWGMIN, IP%XWSAT(JJ,1)-IR%XWGI(JJ,2,1)) 
      ZSOILMAX(JJ) = MAX(0.0,ZWSAT(JJ,1)-IR%XWG(JJ,2,1))*IMX%XDG(JJ,2,1)*XRHOLW/PTSTEP
    ENDDO
  ENDIF
!
  DGEIP%XIFLOOD(:)=MIN(DGEIP%XIFLOOD(:),ZSOILMAX(:))
!
ELSE
!
  DGEIP%XIFLOOD(:)=0.0
!
ENDIF
!
!calculate the infiltration rate
!
PPG  (:) = PPG(:) + DGEIP%XIFLOOD(:)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('HYDRO_SGH',1,ZHOOK_HANDLE)
!
END SUBROUTINE HYDRO_SGH
