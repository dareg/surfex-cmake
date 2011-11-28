!     #########
      SUBROUTINE HYDRO_SGH  (HISBA,HRUNOFF,HRAIN,HHORT,PTSTEP,           &
                               PD_G,PDZG,PWSAT,PWWILT,PWG,PWGI,PPG,      &
                               PPG_MELT,PMUF,PKSAT,PBCOEF,PMPOTSAT,      &
                               PKSAT_ICE,PD_ICE,PFSAT,PHORTON,PDUNNE,    &
                               PFFLOOD,PPIFLOOD,PIFLOOD,PPFLOOD,         &
                               PRUNOFFB,PRUNOFFD,PTDIURN                 )  
!
!     #####################################################################
!
!!****  *HYDRO_SGH*  
!!
!!    PURPOSE
!!    =======
!
!     1. Determine the Horton runoff that take account of a spatial subgrid 
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
!
USE MODD_CSTS,      ONLY : XRHOLW, XDAY
USE MODD_ISBA_PAR,  ONLY : XWGMIN
!
USE MODI_HYDRO_DT92
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
CHARACTER(LEN=*),INTENT(IN)      :: HISBA   ! hydrology/soil:
!                                           ! '2-L'  = single column
!                                           ! '3-L'  = root zone/baseflow layer
!                                           ! 'DIF'  = N-layer diffusion: Richard's Eq.
!                                           !         (Boone and Etchevers 2001)
!
CHARACTER(LEN=*),     INTENT(IN) :: HRUNOFF ! surface runoff formulation
!                                           ! 'WSAT'
!                                           ! 'DT92'
!                                           ! 'SGH ' Topmodel
!
!
CHARACTER(LEN=*), INTENT(IN)     :: HRAIN   ! Rainfall spatial distribution
                                            ! 'DEF' = No rainfall spatial distribution
                                            ! 'SGH' = Rainfall exponential spatial distribution
                                            ! 
!
CHARACTER(LEN=*), INTENT(IN)     :: HHORT   ! Horton runoff
                                            ! 'DEF' = no Horton runoff
                                            ! 'SGH' = Horton runoff
!
REAL, INTENT(IN)                 :: PTSTEP
!                                   timestep of the integration
!
REAL, DIMENSION(:,:), INTENT(IN) :: PWG,PWGI
!                                   PWG   = layer average liquid volumetric water content (m3 m-3)
!                                   PWGI  = layer average frozen volumetric water content (m3 m-3)
!
REAL, DIMENSION(:,:), INTENT(IN) :: PD_G,PDZG,PWSAT,PWWILT
REAL, DIMENSION(:), INTENT(IN)   :: PKSAT,PBCOEF,PFSAT,PMPOTSAT
!                                   PD_G  = layer depth (m)
!                                   PDZG= layer thickness (m)
!                                   PKSAT = surface saturated hydraulic conductivity
!                                   PWSAT = soil porosity (m3 m-3)
!                                   PMPOTSAT = surface saturated matric potential
!                                   PFSAT = satured fraction
!
REAL, DIMENSION(:), INTENT(INOUT):: PPG
REAL, DIMENSION(:), INTENT(IN)   :: PPG_MELT, PMUF
!                                   PPG      = water reaching the ground
!                                   PPG_MELT = snowmelt reaching the ground
!                                   PMUF      = wet fraction reached by rain
!
REAL, DIMENSION(:), INTENT(IN)   :: PKSAT_ICE, PD_ICE
!                                   PKSAT_ICE = hydraulic conductivity at saturation (m s-1)
!                                               on frozen soil depth (Horton calculation)
!                                   PD_ICE    = depth of the soil column for
!                                               fraction of frozen soil calculation (m)
!
REAL, DIMENSION(:), INTENT(OUT)  :: PDUNNE, PHORTON
!                                   PDUNNE  = Dunne runoff
!                                   PHORTON = Horton runoff
!
REAL, DIMENSION(:), INTENT(IN   ) :: PFFLOOD
REAL, DIMENSION(:), INTENT(IN   ) :: PPIFLOOD
REAL, DIMENSION(:), INTENT(INOUT) :: PIFLOOD, PPFLOOD
!                                    PIFLOOD = Floodplain infiltration     [kg/m²/s]
!                                    PPFLOOD = Floodplain interception     [kg/m²/s]
!
REAL, DIMENSION(:), INTENT(IN)    :: PRUNOFFB ! slope of the runoff curve
REAL, DIMENSION(:), INTENT(IN)    :: PRUNOFFD
!                                    PRUNOFFD = depth over which sub-grid runoff calculated (m)
!
REAL, DIMENSION(:), INTENT(IN)    :: PTDIURN
!                                    PTDIURN      = penetration depth for restore (m)
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PPG))                 :: ZPG_INI, ZFSAT,                    &
                                              ZWSAT, ZFROZEN, ZIMAX_ICE, ZIMAX,  &
                                              ZHORT_R, ZHORT_M, ZSOILMAX, ZIF_MAX                   
!                                             ZFROZEN  = frozen soil fraction for runoff
!                                             ZIMAX_ICE    = maximum infiltration rate for frozen soil
!                                             ZIMAX     = maximum infiltration rate for unfrozen soil
REAL, DIMENSION(SIZE(PPG))                 :: ZEFFICE, ZFRZ, ZS, ZD_H
!                                             ZEFFICE   = effective soi ice penetration depth for restore (m)
!                                             ZFRZ     = diffusion impedance factor
REAL, DIMENSION(SIZE(PPG))                 :: ZWG2_AVG, ZWGI2_AVG, ZWSAT_AVG, ZWWILT_AVG
!                                             Average water and ice content
!                                             values over the soil depth D2 (for calculating surface runoff)
REAL, DIMENSION(SIZE(PD_G,1),SIZE(PD_G,2)) :: ZSOILWGHT  ! ISBA-DIF: weights for vertical
!                                                          integration of soil water and properties
REAL, DIMENSION(SIZE(PPG))                 :: ZSOILZ     ! Soil depths from surface to base of each layer (ISBA-DF)     
REAL, DIMENSION(SIZE(PPG))                 :: ZPG_WORK, ZRUISDT
!
INTEGER                                    :: JWRK, INLVLD
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
ZFSAT  (:) = 0.0
!
ZWSAT    (:)  = 0.0
ZFROZEN  (:)  = 0.0
ZIMAX_ICE(:)  = 0.0
ZIMAX    (:)  = 0.0
!
!HRUNOFF = DT92 ZFSAT calculation
ZWG2_AVG(:)   = 0.0
ZWGI2_AVG(:)  = 0.0
ZWSAT_AVG(:)  = 0.0
ZWWILT_AVG(:) = 0.0
!HRUNOFF = DT92 & HISBA = DIF ZFSAT calculation
ZSOILWGHT(:,:) = 0.0
ZSOILZ(:)      = 0.0
!
!HISBA /= DIF ZWSAT calculation
ZEFFICE  (:)   = 0.0
ZFRZ     (:)   = 0.0
ZS       (:)   = 0.0
ZD_H     (:)   = 0.0
!
!HHORT=SGH
ZHORT_R(:) = 0.0
ZHORT_M(:) = 0.0
!
!PIFLOOD calculation
ZSOILMAX(:) = 0.0
ZIF_MAX(:)  = 0.0
!
!HRUNOFF = DT92 DUNNE calculation
ZPG_WORK(:)   = 0.0
ZRUISDT(:)    = 0.0
!
!to limit numerical artifacts
ZPG_INI(:) = PPG(:) + PPG_MELT(:)
!
!-------------------------------------------------------------------------------
!
!*           1. Surface saturated fraction
!            -----------------------------
!
IF(HRUNOFF=='SGH ')THEN
!        
   ZFSAT(:) = PFSAT(:)
!   
ELSEIF(HRUNOFF=='DT92')THEN
!
! Calculate the layer average water content for the sub-grid
! surface runoff computation: use PRUNOFFD as the depth over which
! runoff is calculated.
!
! First, determine a weight for each layer's contribution
! to thickness averaged water content and soil properties for runoff.
!
   IF (HISBA == 'DIF') THEN                                   
!
      INLVLD=SIZE(PD_G,2)
!
      ZSOILWGHT(:,:) = 0.0
      ZSOILZ(:)      = 0.0
      DO JWRK=1,INLVLD
         ZSOILZ(:)         = ZSOILZ(:) + PDZG(:,JWRK)  
         ZSOILWGHT(:,JWRK) = MIN(PDZG(:,JWRK), MAX(0.0, PRUNOFFD(:)-ZSOILZ(:)+PDZG(:,JWRK)))
      ENDDO
!
! Normalize the weights:
!
      DO JWRK=1,INLVLD
         ZSOILWGHT(:,JWRK) = ZSOILWGHT(:,JWRK)/MAX(1.E-6, SUM(ZSOILWGHT(:,:),2))
      ENDDO
!
! Vertically averaged soil properties and moisture for surface runoff computation:
!
      ZWG2_AVG(:)   = SUM(ZSOILWGHT(:,:)*PWG(:,:),   2)
      ZWGI2_AVG(:)  = SUM(ZSOILWGHT(:,:)*PWGI(:,:),  2)
      ZWSAT_AVG(:)  = SUM(ZSOILWGHT(:,:)*PWSAT(:,:), 2)
      ZWWILT_AVG(:) = SUM(ZSOILWGHT(:,:)*PWWILT(:,:),2)
!
   ELSE
      ZWG2_AVG(:)   = PWG(:, 2)
      ZWGI2_AVG(:)  = PWGI(:,2)
      ZWSAT_AVG(:)  = PWSAT(:,1)
      ZWWILT_AVG(:) = PWWILT(:,1)
   ENDIF
!
   IF(HHORT=='SGH')THEN
     !runoff over frozen soil explicitly calculated
     ZWGI2_AVG(:)=0.0
   ENDIF
!
   ZS(:)=MIN(1.0,(ZWG2_AVG(:)+ZWGI2_AVG(:)-ZWWILT_AVG(:))/(ZWSAT_AVG(:)-ZWWILT_AVG(:)))
   ZS(:)=MAX(0.0,ZS(:))
   ZFSAT(:) = 1.0-(1.0-ZS(:))**(PRUNOFFB(:)/(PRUNOFFB(:)+1.))
!
ENDIF
!
!*           2. Horton runoff (à revoir pour DF !!!)
!            ----------------
!
IF(HISBA == 'DIF')THEN
!
  ZWSAT    (:) = PWSAT(:,1)
  ZFROZEN  (:) = 0.
  ZIMAX_ICE(:) = PKSAT(:)
  ZIMAX    (:) = PKSAT(:)
!  
ELSE
!
! Effective frozen depth penetration 
!
  ZEFFICE(:)=PD_G(:,2)*PWGI(:,2)/(PWGI(:,2)+PWG(:,2))
!
! Modify soil porosity as ice assumed to become part
! of solid soil matrix (with respect to liquid flow):
!
  ZWSAT(:) = MAX(XWGMIN, PWSAT(:,1)-PWGI(:,2)) 
!
! calculate the subgrid frozen soil fraction of the grid cells
!
  ZFROZEN (:) = MIN(1.,ZEFFICE(:)/MAX(PD_ICE(:),PTDIURN(:)))
!
! Impedance Factor from (Johnsson and Lundin 1991).
!
  ZFRZ(:) = 10.**(-6.*MIN(1.,ZEFFICE(:)/PTDIURN(:)))
!
! Calculate infiltration MAX on frozen soil as Johnsson and Lundin (1991).
! The max infiltration is equal to the unsaturated conductivity function at a
! water content corresponding to the total porosity less the ice-filled volume.
!
  ZS(:) =MIN(1.,(ZWSAT(:)/PWSAT(:,2)))
  ZIMAX_ICE(:)=ZFRZ(:)*PKSAT_ICE(:)*(ZS(:)**(2*PBCOEF(:)+3.))
!
! Calculate infiltration MAX on unfrozen soil using green-ampt approximation
!    
  ZS   (:)=MIN(1.,PWG(:,2)/ZWSAT(:))
  ZD_H (:)=MIN(0.10,PD_G(:,2))
  ZIMAX(:)=PKSAT(:)*(PBCOEF(:)*PMPOTSAT(:)*(ZS(:)-1.0)/ZD_H(:)+1.0)
!
ENDIF
!
IF(HHORT=='SGH')THEN
!
! calculate the Horton runoff generated by the rainfall rate
!
  IF(HRAIN=='SGH')THEN
!
    WHERE(PPG(:)>0.)
       ZHORT_R(:) = (1.- ZFROZEN(:))* PPG(:)/((ZIMAX    (:)*XRHOLW*PMUF(:)/PPG(:)) + 1.) & !unfrozen soil
                  +      ZFROZEN(:) * PPG(:)/((ZIMAX_ICE(:)*XRHOLW*PMUF(:)/PPG(:)) + 1.)   !frozen soil
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
  WHERE(PFFLOOD(:)<=ZFSAT(:))
        PHORTON(:) = (1. - ZFSAT(:)) * (ZHORT_R(:)+ZHORT_M(:))
  ELSEWHERE
        PHORTON(:) = (1. - PFFLOOD(:)) * (ZHORT_R(:) + ZHORT_M(:))
  ENDWHERE
!
! calculate all water reaching the ground
!
  PPG  (:) = PPG(:) + PPG_MELT(:)
!
ELSE
!
  PHORTON(:) = 0.0
!
! calculate all water reaching the ground
!
  PPG  (:) = PPG(:) + PPG_MELT(:)        
!
ENDIF
!
!*           3. Dunne runoff and flood interception
!            --------------------------------------
!
IF(HRUNOFF=='SGH ')THEN
!        
! calculate the Dunne runoff with TOPMODEL
!
  PDUNNE(:) = MAX(PPG(:),0.0) * MAX(ZFSAT(:)-PFFLOOD(:),0.0)
!
ELSEIF (HRUNOFF=='DT92')THEN
!
!*       Dumenil et Todini (1992)  RUNOFF SCHEME
!        ---------------------------------------         
!
! surface runoff exit only on the Fsat-Fflood fraction
!
  WHERE(PFFLOOD(:)<ZFSAT(:))
        ZPG_WORK(:) = PPG(:) - PHORTON(:) - PFFLOOD(:)*MAX(0.0,PPG(:))
  ELSEWHERE
        ZPG_WORK(:) = 0.0
  ENDWHERE
!
   CALL HYDRO_DT92(PTSTEP,                                &
                  PRUNOFFB, ZWWILT_AVG,                   &
                  PRUNOFFD, ZWSAT_AVG,                    &
                  ZWG2_AVG, ZWGI2_AVG,                    &
                  ZPG_WORK, ZRUISDT                       )
!
   PDUNNE(:) = ZRUISDT(:)*PRUNOFFD(:)*XRHOLW/PTSTEP
!
ELSE
! 
! Default case (no subgrid runoff)
!
  PDUNNE(:) = 0.0
!
ENDIF
!
! Interception by the flooplains
!
PPFLOOD(:)=PFFLOOD(:)*MAX(0.0,PPG(:))
!
! calculate the infiltration rate after runoff
!
PPG  (:) = PPG(:) - PDUNNE(:) - PHORTON(:) - PPFLOOD(:)
!
! Supress numerical artifacts:
!
WHERE (ZPG_INI(:)<0.0)
       PPG(:)     = ZPG_INI(:)
       PHORTON(:) = 0.0
       PDUNNE (:) = 0.0
       PPFLOOD(:) = 0.0
ENDWHERE
!
!*           4. infiltration rate from floodplains (à revoir pour DF !!!)
!            -------------------------------------
!
! calculate the maximum flood infiltration
!
ZIF_MAX(:) = (1.- ZFROZEN(:))* ZIMAX    (:)*XRHOLW &   !unfrozen soil
           +      ZFROZEN(:) * ZIMAX_ICE(:)*XRHOLW     !frozen soil
!
PIFLOOD(:)=MAX(0.0,(PFFLOOD(:)-ZFSAT(:)))*MIN(PPIFLOOD(:),ZIF_MAX(:))
!
IF(HISBA == 'DIF')THEN
  ZSOILMAX(:) = 0.0
ELSE
  ZSOILMAX(:) = MAX(0.0,ZWSAT(:)-PWG(:,2))*PD_G(:,2)*XRHOLW/PTSTEP
ENDIF
!
PIFLOOD(:)=MIN(PIFLOOD(:),ZSOILMAX(:))
!
!calculate the infiltration rate from floodplains
!
PPG  (:) = PPG(:) + PIFLOOD(:)
!
!-------------------------------------------------------------------------------
!
!deallocate
!
IF (LHOOK) CALL DR_HOOK('HYDRO_SGH',1,ZHOOK_HANDLE)
!
END SUBROUTINE HYDRO_SGH
