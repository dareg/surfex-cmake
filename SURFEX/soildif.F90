!     #########
      SUBROUTINE SOILDIF( HSCOND, HDIFSFCOND, HDIF,                              &
                         PVEG, PCV, PFFG, PFFV,                                  &
                         PCG, PCGMAX, PCT, PFROZEN1,                             &
                         PWG, PWGI,                                              &
                         PHCAPSOILZ, PCONDDRYZ, PCONDSLDZ,                       &
                         PBCOEF, PWSAT, PMPOTSAT, PALPHA, PN, PM,                &
                         PWRES, PSOILCONDZ, PSOILHCAPZ                           )  
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
!!	A. Boone           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    
!!                  25/03/99      (Boone)   Added Johansen (1975)/Peters-Lidard 
!!                                          option to explicitly compute CG
!!                  08/25/02      (Boone)   DIF option code only
!!                  25/05/08     (Decharme) Add Flood properties
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,       ONLY : XCL, XCI, XRHOLW, XRHOLI, XPI, XDAY, XCONDI
USE MODD_ISBA_PAR,   ONLY : XCONDWTR, XWGMIN
USE MODD_SURF_PAR,   ONLY : XUNDEF
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
!
CHARACTER(LEN=*),     INTENT(IN)   :: HSCOND  ! thermal conductivity formulation
!                                             ! 'NP89' :  Noilhan and Planton 
!                                             !  (1989: McCumber-Pielke (1981) and
!                                             !  Clapp and Hornberger (1978))
!                                             ! 'PL98' Method of Johansen (1975) as
!                                             ! presented by Peters-Lidard (JAS: 1998)
!
CHARACTER(LEN=*),     INTENT(IN)  :: HDIFSFCOND ! NOTE: Only used when HISBA = DIF
!                                               ! MLCH' = include the insulating effect of leaf
!                                               !         litter/mulch on the surface thermal cond.
!                                               ! 'DEF' = no mulch effect
!
CHARACTER(LEN=*),     INTENT(IN)  :: HDIF       ! NOTE: Only used when HISBA = DIF
!                                               ! 'BC' = Brook and Corey
!                                               ! 'VG' = Van Genuchten
!
REAL, DIMENSION(:), INTENT(IN)    :: PVEG, PCV
!                                      Soil and vegetation parameters
!                                      PVEG = fraction of vegetation
!                                      PCV  = the heat capacity of the vegetation
!
REAL, DIMENSION(:,:), INTENT(IN)  :: PHCAPSOILZ, PCONDDRYZ, PCONDSLDZ 
!                                    PHCAPSOILZ = soil heat capacity [J/(K m3)]
!                                    PCONDDRYZ  = soil dry thermal conductivity 
!                                                 [W/(m K)] 
!                                    PCONDSLDZ  = soil solids thermal conductivity 
!                                                 [W/(m K)]
!
REAL, DIMENSION(:,:), INTENT(IN)  :: PBCOEF, PWSAT, PMPOTSAT
!                                    PBCOEF   = profile of b-parameter (-)
!                                    PWSAT    = profile of porosity (m3/m3)
!                                    PMPOTSAT = profile of matric potential at saturation (m)
!
REAL, DIMENSION(:,:), INTENT(IN)  :: PWG, PWGI 
!                                    PWG    = soil liquid water content (m3/m3)
!                                    PWGI   = soil frozen water content (m3/m3)
!
REAL, DIMENSION(:), INTENT(OUT)   :: PFROZEN1, PCG, PCT
!                                      PFROZEN1 = fraction of ice in superficial soil
!                                      PCT      = averaged surface heat capacity of the grid (m2 K J-1)
!                                      PCG      = averaged surface soil heat capacity (m2 K J-1)
!
REAL,               INTENT(IN)    :: PCGMAX
!                                      Maximum soil heat capacity
!
REAL, DIMENSION(:,:), INTENT(OUT) :: PSOILCONDZ, PSOILHCAPZ
!                                    PSOILHCAP = soil heat capacity        (J m-3 K-1)
!                                    PSOILCOND = soil thermal conductivity (W m-1 K-1)
!
!
REAL, DIMENSION(:), INTENT(IN)   :: PFFV, PFFG
!                                   PFFG = Floodplain fraction over the ground
!                                   without snow (ES)
!                                   PFFV = Floodplain fraction over vegetation
!                                   without snow (ES)
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PALPHA, PN, PM, PWRES
!                                     Van Genuchten parameter
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PWG,1),SIZE(PWG,2)) :: ZMATPOT, ZFROZEN2DF, ZUNFROZEN2DF, &
                                  ZCONDSATDF, ZSATDEGDF, ZKERSTENDF  
!                               ZMATPOT    = soil matric potential (m)
!
REAL, PARAMETER              :: ZVEGMULCH  = 0.10 ! reduction factor for the surface layer
!                                                 ! thermal condictivity due to the presence
!                                                 ! of mulch/organic material (ISBA-DF)
!                                                 ! Set to 1 to remove this effect. 
!
REAL, DIMENSION(SIZE(PVEG)) :: ZFF, ZCF !heat capacity of the flood
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       0.     Initialization:
!               ---------------
!
IF (LHOOK) CALL DR_HOOK('SOILDIF',0,ZHOOK_HANDLE)
ZMATPOT(:,:)      = XUNDEF
ZFROZEN2DF(:,:)   = XUNDEF
ZUNFROZEN2DF(:,:) = XUNDEF
ZCONDSATDF(:,:)   = XUNDEF
ZSATDEGDF(:,:)    = XUNDEF
ZKERSTENDF(:,:)   = XUNDEF
!
ZCF = XUNDEF
!
!-------------------------------------------------------------------------------
!
!*       1.     THE HEAT CAPACITY OF BARE-GROUND
!               --------------------------------
!               Explicit soil thermal diffusion option:
!
IF(HSCOND == 'NP89')THEN
!
! Surface soil water reservoir frozen fraction:
!
  PFROZEN1(:) = PWGI(:,1)/(PWGI(:,1) + MAX(PWG(:,1),XWGMIN))
!
!
! Calculate the thermal conductivity [W/(m K)] using the method of McCumber and 
! Pielke (1981) used implicitly by Noilhan and Planton (1989).
! First calculate the soil water potential using Clapp and Hornberger (1978)
! (m): NOTE that this method DOES NOT explicitly account for soil ice
! so use a ice-weighted average (to prevent excessively low values
! when a soil layer totally frozen): 
!
  ZMATPOT(:,:)    = PSI_FUNC(HDIF,PWG,PWSAT,PMPOTSAT,PBCOEF,PALPHA,PN,PM,PWRES)
  PSOILCONDZ(:,:) = 418.*EXP( -(LOG10(-ZMATPOT(:,:))+4.7) )
  PSOILCONDZ(:,:) = MAX(0.171, PSOILCONDZ(:,:))
  PSOILCONDZ(:,:) = (1.0-PWSAT(:,:)+PWG(:,:)+PWGI(:,:))                      &
                         /( (PWGI(:,:)/XCONDI) +                               &
                         ((1.0-PWSAT(:,:)+PWG(:,:))/PSOILCONDZ(:,:)) )  
!
ELSE
!
! Calculate thermal conductivity using PL98, but for explicit layers:
!
  ZFROZEN2DF(:,:)   = PWGI(:,:)/(PWGI(:,:) + MAX(PWG(:,:),XWGMIN))
!
  PFROZEN1(:)       = ZFROZEN2DF(:,1)
!
  ZUNFROZEN2DF(:,:) = (1.0-ZFROZEN2DF(:,:))*PWSAT(:,:)
!
  ZCONDSATDF(:,:)   = (PCONDSLDZ(:,:)**(1.0-PWSAT(:,:)))*                 &
                        (XCONDI**(PWSAT(:,:)-ZUNFROZEN2DF(:,:)))*           &
                        (XCONDWTR**ZUNFROZEN2DF(:,:))  
!
   ZSATDEGDF(:,:)    = MAX(0.1, PWG(:,:)/PWSAT(:,:))
!
   ZKERSTENDF(:,:)   = LOG10(ZSATDEGDF(:,:)) + 1.0
!
   ZKERSTENDF(:,:)  = (1.0-ZFROZEN2DF(:,:))*ZKERSTENDF(:,:) +               &
                             ZFROZEN2DF(:,:) *ZSATDEGDF(:,:)  
!
! Thermal conductivity of soil:
!
   PSOILCONDZ(:,:)  = ZKERSTENDF(:,:)*(ZCONDSATDF(:,:)-PCONDDRYZ(:,:))      &
                        + PCONDDRYZ(:,:)  
!
ENDIF
!
!
! This takes into account the insulating effect of dead vegetation/leaf litter/mulch on
! the uppermost soil layer thermal conductivity: it is a simple modification
! of the ideas presented by Gonzalez-Sosa et al., AFM, 1999: the thermal
! conductivity is reduced by the factor 'ZVEGMULCH'. The main impact is
! to reduce the thermal coupling between the surface layer and the 
! sub-surface soil. In the limit when
! there is no vegetation, the conductivity collapses into the bare-soil value.
! If the option is not in force ( HDIFSFCOND /= 'MLCH') then use only soil
! properties (no mulch effect). 
!
IF(HDIFSFCOND == 'MLCH') PSOILCONDZ(:,1) = (1.0 - PVEG(:) * (1.0 - ZVEGMULCH)) * PSOILCONDZ(:,1) 
!
! Soil Heat capacity [J/(m3 K)]
!
PSOILHCAPZ(:,:) = (1.0-PWSAT(:,:) )*PHCAPSOILZ(:,:)     +         &
                         PWG(:,:)  *XCL*XRHOLW            +         &
                         PWGI(:,:) *XCI*XRHOLI       
!
! Surface soil thermal inertia [(m2 K)/J]
!
PCG(:) = 2.*SQRT(XPI/(PSOILCONDZ(:,1)*PSOILHCAPZ(:,1)*XDAY))
PCG(:) = MIN( PCG(:), PCGMAX )
!
!-------------------------------------------------------------------------------
!
!*       2.     THE HEAT CAPACITY OF FLOOD
!               --------------------------------
!
ZFF(:) = PVEG(:)*PFFV(:) + (1.-PVEG(:))*PFFG(:)
!
WHERE (ZFF(:) > 0.)                                                 
       ZCF(:) = 2.0 * SQRT( XPI/(XCONDWTR*XRHOLW*XCL*XDAY) )
END WHERE                 
!
!-------------------------------------------------------------------------------
!
!*      3.      GRID-AVERAGED HEAT CAPACITY
!               ---------------------------
!
! With contribution from the ground, flood and vegetation for explicit
! (ISBA-ES) snow scheme option (i.e. no snow effects included here):
!
PCT(:) = 1. / ( (1.-PVEG(:))*(1.-PFFG(:)) / PCG(:)     &
                   +  PVEG(:) *(1.-PFFV(:)) / PCV(:)     &
                   +  ZFF (:)               / ZCF(:)     )  
IF (LHOOK) CALL DR_HOOK('SOILDIF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SOILDIF
