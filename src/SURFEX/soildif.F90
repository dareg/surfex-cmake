!     #########
      SUBROUTINE SOILDIF(IO, IP, IMX, IR, DGMI, PVEG, PCV, PPIFLOOD, &
                         PFROZEN1, PFFG, PFFV, PSOILCONDZ, PSOILHCAPZ   )
!     ##########################################################################
!
!!****  *SOIL*  
!!
!!    PURPOSE
!!    -------
!
!     Calculates the coefficients related to the soil (i.e., surface heat capacities, CG, CT,
!     and thermal conductivity and heat capacity profiles)
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
!!    USE MODD_PARAMETERS
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    Noilhan and Planton (1989)
!!    Belair (1995)
!!    Boone (2000)
!!      
!!    AUTHOR
!!    ------
!!      A. Boone           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    
!!                  25/03/99      (Boone)   Added Johansen (1975)/Peters-Lidard 
!!                                          option to explicitly compute CG
!!                  08/25/02      (Boone)   DIF option code only
!!                  25/05/08     (Decharme) Add Flood properties
!!                  03/08/11     (Decharme) Optimization
!!                     04/13     (Decharme) good soil moisture extrapolation computation
!!                  23/07/13     (Decharme) Surface / Water table depth coupling
!!                  23/10/14     (Decharme) revise all thermo properties
!!                                          delete NP89 option for thermal cond
!!                                          because not physical with explicit soil.
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
!
USE MODD_CSTS,       ONLY : XCL, XCI, XRHOLW, XRHOLI, XPI, XDAY, XCONDI, XTT, XLMTT, XG
USE MODD_ISBA_PAR,   ONLY : XCONDWTR, XWGMIN, XWTD_MAXDEPTH, & 
                            XOMRHO, XOMSPH, XOMCONDDRY,      &
                            XOMCONDSLD, XCVHEATF
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMX
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
!
REAL, DIMENSION(:), INTENT(IN)    :: PVEG, PCV
!                                      Soil and vegetation parameters
!                                      PVEG = fraction of vegetation
!                                      PCV  = the heat capacity of the vegetation
!
REAL, DIMENSION(:), INTENT(IN)   :: PPIFLOOD ! Floodplain potential infiltration or water mass (kg m-2)
!
REAL, DIMENSION(:), INTENT(OUT)   :: PFROZEN1
!                                      PFROZEN1 = fraction of ice in superficial soil
!
REAL, DIMENSION(:), INTENT(IN)   :: PFFV, PFFG
!                                   PFFG = Floodplain fraction over the ground
!                                   without snow (ES)
!                                   PFFV = Floodplain fraction over vegetation
!                                   without snow (ES)
!
REAL, DIMENSION(:,:), INTENT(OUT) :: PSOILCONDZ, PSOILHCAPZ
!                                    PSOILHCAP = soil heat capacity        (J m-3 K-1)
!                                    PSOILCOND = soil thermal conductivity (W m-1 K-1)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(IR%XTG,1),SIZE(IR%XTG,2)) :: ZMATPOT, ZCONDDRYZ, ZCONDSLDZ, ZVEGMULCH
!                                           ZMATPOT    = soil matric potential (m)
!
REAL                         :: ZFROZEN2DF, ZUNFROZEN2DF, ZCONDSATDF, ZLOG_CONDI, ZLOG_CONDWTR,  &
                                ZSATDEGDF, ZKERSTENDF, ZWORK1, ZWORK2, ZWORK3, ZLOG, ZWTOT, ZWL
!    
REAL, PARAMETER              :: ZCTMAX = 1.E-4 ! Maximum thermal inertia
!
REAL, PARAMETER              :: ZTHICKM = 0.04 ! Mulch thickness (m)
!
REAL, DIMENSION(SIZE(PVEG)) :: ZFF, ZCF, ZCV !Thermal inertia of the flood or vegetation
!
REAL, DIMENSION(SIZE(PVEG)) :: ZWTD ! Water table depth if no coupling (m)  
!
INTEGER :: INI, INL, JJ, JL, IDEPTH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       0.     Initialization:
!               ---------------
!
IF (LHOOK) CALL DR_HOOK('SOILDIF',0,ZHOOK_HANDLE)
!
INI=SIZE(IR%XWG,1)
INL=SIZE(IR%XWG,2)
!
ZFF (:) = 0.0
ZCF (:) = XUNDEF
!
!-------------------------------------------------------------------------------
!
!*       1.     WATER TABLE DETH ADJUSTMENT FOR ISBA (m)
!               -----------------------------------------
!
WHERE(IP%XWTD(:)==XUNDEF)
! no water table / surface coupling over some regions        
  ZWTD     (:) = XWTD_MAXDEPTH 
ELSEWHERE
  ZWTD     (:) = IP%XFWTD(:)/MAX(-IP%XWTD(:),0.001) + (1.0-IP%XFWTD(:))/MAX(-IP%XWTD(:),XWTD_MAXDEPTH)
  ZWTD     (:) = 1.0/ZWTD(:)
ENDWHERE
!
!-------------------------------------------------------------------------------
!
!*       2.     MATRIC POTENTIAL AND MOISTURE EXTRAPOLATION
!               -------------------------------------------
!
DO JL=1,INL
   DO JJ=1,INI
!   
      IDEPTH=IMX%NWG_LAYER(JJ,1)
      IF(JL>IDEPTH)THEN
!                           
!       total matric potential
        ZWORK1  = MIN(1.0,(IR%XWG(JJ,IDEPTH,1)+IR%XWGI(JJ,IDEPTH,1))/IP%XWSAT(JJ,IDEPTH))
        ZLOG    = IP%XBCOEF(JJ,IDEPTH)*LOG(ZWORK1)
        ZMATPOT(JJ,IDEPTH) = IP%XMPOTSAT(JJ,IDEPTH)*EXP(-ZLOG)

!       extrapolation of total matric potential
        ZWORK1         = 0.5*(IMX%XDG(JJ,IDEPTH,1) + IMX%XDG(JJ,IDEPTH-1,1))
        ZWORK2         = 0.5*(IMX%XDG(JJ,JL,    1) + IMX%XDG(JJ,JL-1,1))
        ZWORK3         = MAX(0.0,(ZWTD(JJ)-ZWORK2)/(ZWORK2-ZWORK1))
        ZMATPOT(JJ,JL) = (IP%XMPOTSAT(JJ,JL)+ZWORK3*ZMATPOT(JJ,IDEPTH))/(1.0+ZWORK3)
!
!       total soil water content computation
        ZWORK1      = MAX(1.0,ZMATPOT(JJ,JL)/IP%XMPOTSAT(JJ,JL))
        ZLOG        = LOG(ZWORK1)/IP%XBCOEF(JJ,JL)
        ZWTOT       = IP%XWSAT(JJ,JL)*EXP(-ZLOG)
        ZWTOT       = MAX(XWGMIN,ZWTOT)
!
!       soil liquid water content computation
        ZMATPOT(JJ,JL) = MIN(IP%XMPOTSAT(JJ,JL),XLMTT*(IR%XTG(JJ,JL,1)-XTT)/(XG*IR%XTG(JJ,JL,1)))
!        
        ZWORK1      = MAX(1.0,ZMATPOT(JJ,JL)/IP%XMPOTSAT(JJ,JL))
        ZLOG        = LOG(ZWORK1)/IP%XBCOEF(JJ,JL)
        ZWL         = IP%XWSAT(JJ,JL)*EXP(-ZLOG)
        ZWL         = MAX(ZWL,XWGMIN)
        IR%XWG (JJ,JL,1) = MIN(ZWL,ZWTOT )
!        
!       soil ice computation        
        IR%XWGI(JJ,JL,1) = MAX(0.0,ZWTOT-IR%XWG(JJ,JL,1))
!        
      ENDIF
   ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       3.     SURFACE FROZEN FRACTION
!               -----------------------
!
!
! Surface soil water reservoir frozen fraction:
!
PFROZEN1(:) = IR%XWGI(:,1,1)/(IR%XWGI(:,1,1) + MAX(IR%XWG(:,1,1),XWGMIN))
!
!-------------------------------------------------------------------------------
!
!*       4.     SIMPLE LITTER/MULCH EFFECT
!               --------------------------
!
! This takes into account the insulating effect of dead vegetation/leaf litter/mulch on
! uppermost soil layers thermal properties. Use organic matter thermal properties.
!
!
ZCONDDRYZ (:,:) = IP%XCONDDRY (:,:)
ZCONDSLDZ (:,:) = IP%XCONDSLD (:,:)
!
IF(IO%CDIFSFCOND == 'MLCH') THEN
!  
  DO JL=1,INL
     DO JJ=1,INI  
!
        ZVEGMULCH(JJ,JL) = PVEG(JJ)*MIN(IP%XDZG(JJ,JL,1),&
                MAX(0.0,ZTHICKM-IMX%XDG(JJ,JL,1)+IP%XDZG(JJ,JL,1)))/IP%XDZG(JJ,JL,1)     
!
        IF(ZVEGMULCH(JJ,JL)>0.0)THEN
           ZCONDDRYZ (JJ,JL) = 1.0/((1.0-ZVEGMULCH(JJ,JL))/IP%XCONDDRY(JJ,JL)+ZVEGMULCH(JJ,JL)/XOMCONDDRY)
           ZCONDSLDZ (JJ,JL) = 1.0/((1.0-ZVEGMULCH(JJ,JL))/IP%XCONDSLD(JJ,JL)+ZVEGMULCH(JJ,JL)/XOMCONDSLD)
        ENDIF
!
     ENDDO
  ENDDO
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       5.     THE THERMAL CONDUCTIVITY OF BARE-GROUND
!               ---------------------------------------
!
! Calculate thermal conductivity using PL98 :
!
ZLOG_CONDI   = LOG(XCONDI)
ZLOG_CONDWTR = LOG(XCONDWTR)
!
DO JL=1,INL
   DO JJ=1,INI
!     
      ZFROZEN2DF   = IR%XWGI(JJ,JL,1)/(IR%XWGI(JJ,JL,1) + MAX(IR%XWG(JJ,JL,1),XWGMIN))
      ZUNFROZEN2DF = (1.0-ZFROZEN2DF)*IP%XWSAT(JJ,JL)
!
!Old: CONDSATDF=(CONDSLDZ**(1.0-WSAT))*(CONDI**(WSAT-UNFROZEN2DF))*(CONDWTR**UNFROZEN2DF)  
      ZWORK1      = LOG(ZCONDSLDZ(JJ,JL))*(1.0-IP%XWSAT(JJ,JL))
      ZWORK2      = ZLOG_CONDI*(IP%XWSAT(JJ,JL)-ZUNFROZEN2DF)
      ZWORK3      = ZLOG_CONDWTR*ZUNFROZEN2DF
      ZCONDSATDF  = EXP(ZWORK1+ZWORK2+ZWORK3)
!
      ZSATDEGDF   = MAX(0.1, (IR%XWGI(JJ,JL,1)+IR%XWG(JJ,JL,1))/IP%XWSAT(JJ,JL))
      ZSATDEGDF   = MIN(1.0,ZSATDEGDF)
      ZKERSTENDF  = LOG10(ZSATDEGDF) + 1.0
      ZKERSTENDF  = (1.0-ZFROZEN2DF)*ZKERSTENDF + ZFROZEN2DF *ZSATDEGDF  
!
! Thermal conductivity of soil:
!
      PSOILCONDZ(JJ,JL) = ZKERSTENDF*(ZCONDSATDF-ZCONDDRYZ(JJ,JL)) + ZCONDDRYZ(JJ,JL)  
!
   ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       6.     THE HEAT CAPACITY OF BARE-GROUND
!               --------------------------------
!
! Soil Heat capacity [J/(m3 K)]
!
DO JL=1,INL
   DO JJ=1,INI
      PSOILHCAPZ(JJ,JL) = (1.0-IP%XWSAT(JJ,JL))*IP%XHCAPSOIL(JJ,JL) +         &
                               IR%XWG  (JJ,JL,1) *XCL*XRHOLW        +       &
                               IR%XWGI (JJ,JL,1) *XCI*XRHOLI    
   ENDDO
ENDDO
!
! Surface soil thermal inertia [(m2 K)/J]
!
DGMI%XCG(:) = 1.0 / ( IMX%XDG(:,1,1) * PSOILHCAPZ(:,1) )
!
DGMI%XCG(:) = MIN(ZCTMAX,DGMI%XCG(:))
!
!-------------------------------------------------------------------------------
!
!*       7.     THE HEAT CAPACITY OF VEGETATION
!               --------------------------------
!
! Vegetation thermal inertia [(m2 K)/J]
!
ZCV(:) = 1.0 / ( XCVHEATF/PCV(:) +  XCL * IR%XWR(:,1) )
!
ZCV(:) = MIN(ZCTMAX,ZCV(:))
!
!-------------------------------------------------------------------------------
!
!*       8.     THE HEAT CAPACITY OF FLOOD
!               --------------------------------
!
IF(IO%LFLOOD)THEN
!
  ZFF(:) = PVEG(:)*PFFV(:) + (1.-PVEG(:))*PFFG(:)
!
  WHERE(ZFF(:)>0.0)
    ZCF(:) = 1.0 / ( XCL * PPIFLOOD(:) )
  ENDWHERE
!
  ZCF(:) = MIN(ZCTMAX,ZCF(:))
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      9.      GRID-AVERAGED HEAT CAPACITY
!               ---------------------------
!
! With contribution from the ground, flood and vegetation for explicit
! (ISBA-ES) snow scheme option (i.e. no snow effects included here):
!
DGMI%XCT(:) = 1. / ( (1.-PVEG(:))*(1.-PFFG(:)) / DGMI%XCG(:)     &
                 +  PVEG(:) *(1.-PFFV(:)) / ZCV(:)     &
                 +  ZFF (:)               / ZCF(:)     )  
!
!-------------------------------------------------------------------------------
!
!*      10.     RESTORE DEFAULT VALUES
!               ----------------------
!
! restore default moisture and ice values under moisture soil depth
!
DO JL=1,INL
   DO JJ=1,INI
      IDEPTH=IMX%NWG_LAYER(JJ,1)
      IF(JL>IDEPTH)THEN
        IR%XWG (JJ,JL,1) = XUNDEF
        IR%XWGI(JJ,JL,1) = XUNDEF
      ENDIF
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SOILDIF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SOILDIF
