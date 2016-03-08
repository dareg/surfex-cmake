!     #########
      SUBROUTINE COTWOINIT_n (IO,PGMES,PCO2,PGC,PDMAX, IP  )  
!     #######################################################################
!
!!****  *COTWOINIT*  
!!
!!    PURPOSE
!!    -------
!
!     Initialize model to calculate net assimilation of 
!     CO2 and leaf conductance.
!              
!!**  METHOD
!!    ------
!     Calvet at al (1998) [from model of Jacobs(1994)]
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    USE MODD_CO2V_PAR
!!    USE MODI_COTWO  
!!
!!    REFERENCE
!!    ---------
!!
!!    Calvet et al. (1998)
!!      
!!    AUTHOR
!!    ------
!!
!!      A. Boone           * Meteo-France *
!!      (following Belair)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/10/97
!!      (V. Rivalland) 10/04/02  Add: IP%XAH and IP%XBH coefficients for
!!                               herbaceous water stress response
!!      (P. LeMoigne) 03/2004:   computation of zgmest in SI units
!!      (P. LeMoigne) 10/2004:   possibility of 2 different FZERO
!!      (L. Jarlan)   10/2004:   initialization of DMAX
!!      P Le Moigne   09/2005    AGS modifs of L. Jarlan
!!      S. Lafont     03/2009    change unit of AMAX
!!      A.L. Gibelin  04/2009    TAU_WOOD for NCB option 
!!      A.L. Gibelin  04/2009    Suppress useless GPP and RDK arguments 
!!      A.L. Gibelin  07/2009    Suppress PPST and PPSTF as outputs
!!      B. Decharme   05/2012    Optimization
!!      R. Alkama     05/2012    add 7 new vegtype (19  instead 12)
!!      C. Delire     01/2014    Define a dummy LAI from top and total lai for Dark respiration 
!!
!-------------------------------------------------------------------------------
!
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE, NVT_C3, NVT_C4, NVT_IRR, NVT_TROG,     &
                                NVT_TEBD, NVT_BONE, NVT_TRBE, NVT_TRBD, NVT_TEBE,&
                                NVT_TENE, NVT_BOBD, NVT_BOND, NVT_SHRB, NVT_GRAS
USE MODD_CSTS,           ONLY : XMD
USE MODD_CO2V_PAR,       ONLY : XTOPT, XFZERO1, XFZERO2, XFZEROTROP, XEPSO, XGAMM, XQDGAMM, &
                                  XQDGMES, XT1GMES, XT2GMES, XAMAX,               &
                                  XQDAMAX, XT1AMAX, XT2AMAX, XAH, XBH,            &
                                  XDSPOPT, XIAOPT, XAW, XBW, XMCO2, XMC, XTAU_WOOD  
! 
USE MODE_COTWO,          ONLY : GAULEG
USE MODI_COTWO  
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
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
!
REAL,DIMENSION(:),INTENT(IN)  :: PGMES, PCO2
!                                     PGMES     = mesophyll conductance (m s-1)
!                                     PCO2      = atmospheric CO2 concentration
!
REAL,DIMENSION(:),INTENT(IN)   :: PDMAX, PGC    
!                                     PDMAX     = maximum air saturation deficit tolerate
!                                                 by vegetation
!                                     PGC       = cuticular conductance (m s-1)
!
!
!*      0.2    declaration of local variables
!
INTEGER                           :: JCLASS    ! indexes for loops
INTEGER                           :: ICLASS    ! indexes for loops
INTEGER                           :: ICO2TYPE  ! type of CO2 vegetation
INTEGER                           :: IRAD      ! with or without new radiative transfer
!
REAL, DIMENSION(SIZE(IP%XANMAX,1))     :: ZGS, ZGAMMT, ZTOPT, ZANMAX, ZGMEST, ZGPP, ZRDK, ZEPSO
!                                    ZTOPT     = optimum  temperature for compensation 
!                                                point
!                                    ZANMAX    = maximum photosynthesis rate
!                                    ZGS       = leaf conductance
!                                    ZGAMMT    = temperature compensation point
!                                    ZGPP      = gross primary production
!                                    ZRDK      = dark respiration
!
!
REAL, DIMENSION(SIZE(IP%XANMAX,1))     :: ZCO2INIT3, ZCO2INIT4, ZCO2INIT5, ZCO2INIT2,ZCO2INIT1
!                                    working arrays for initializing surface 
!                                    temperature, saturation deficit, global radiation,
!                                    optimum temperature for determining maximum 
!                                    photosynthesis rate, and soil water stress (none)
REAL, DIMENSION(SIZE(PDMAX))      :: ZDMAX
REAL, DIMENSION(SIZE(PDMAX))      :: ZWORK
!                                    Local variable in order to initialise DMAX
!                                    following Calvet, 2000 (AST or LST cases)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('COTWOINIT_N',0,ZHOOK_HANDLE)
!
ZTOPT  (:) = 0.
IP%XFZERO (:,1) = 0.
IP%XEPSO  (:,1) = 0.
IP%XGAMM  (:,1) = 0.
IP%XQDGAMM(:,1) = 0.
IP%XQDGMES(:,1) = 0.
IP%XT1GMES(:,1) = 0.
IP%XT2GMES(:,1) = 0.
IP%XAMAX  (:,1) = 0.
IP%XQDAMAX(:,1) = 0.
IP%XT1AMAX(:,1) = 0.
IP%XT2AMAX(:,1) = 0.
IP%XTAU_WOOD(:,1) = 0.
!
IP%XAH    (:,1) = 0.
IP%XBH    (:,1) = 0.
!
ZEPSO (:) = 0.
ZGPP (:) = 0.
ZRDK (:) = 0.
ZGAMMT (:) = 0.
ZANMAX (:) = 0.
ZGMEST (:) = 0.
ZCO2INIT3(:) = 0.
ZCO2INIT4(:) = 0.
ZCO2INIT5(:) = 0.
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
! DETERMINE GAUSSIAN WEIGHTS NEEDED FOR CO2 MODEL 
! -----------------------------------------------
!
 CALL GAULEG(0.0,1.0,IP%XABC,IP%XPOI,SIZE(IP%XABC))
!
!
! INITIALIZE VARIOUS PARAMETERS FOR CO2 MODEL:
! --------------------------------------------
! as a function of CO2 vegetation class, C3=>1, C4=>2
!
DO JCLASS=1,NVEGTYPE
  !
  IF (JCLASS==NVT_C4 .OR. JCLASS==NVT_IRR .OR. JCLASS==NVT_TROG) THEN
    ICO2TYPE = 2   ! C4 type
  ELSE
    ICO2TYPE = 1   ! C3 type
  END IF
  IF(IO%LAGRI_TO_GRASS.AND.(JCLASS==NVT_C4 .OR. JCLASS==NVT_IRR)) ICO2TYPE = 1
  IF (IO%LTR_ML) THEN
    IRAD = 1   ! running with new radiative transfer
  ELSE
    IRAD = 2
  ENDIF
  !
  ZTOPT  (:) = ZTOPT  (:) + XTOPT  (ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IF (IO%CPHOTO == 'AGS' .OR. IO%CPHOTO == 'LAI') THEN
     IP%XFZERO (:,1) = IP%XFZERO (:,1) + XFZERO1 (ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  ELSE
     IF((JCLASS==NVT_TEBD) .OR. (JCLASS==NVT_BONE) .OR.                         &
        (JCLASS==NVT_TRBD) .OR. (JCLASS==NVT_TEBE) .OR. (JCLASS==NVT_TENE) .OR. &
        (JCLASS==NVT_BOBD) .OR. (JCLASS==NVT_BOND) .OR. (JCLASS==NVT_SHRB)) THEN
        IP%XFZERO (:,1) = IP%XFZERO (:,1) + ((XAW - LOG(PGMES(:)*1000.0))/XBW)*IP%XVEGTYPE_PATCH(:,JCLASS,1)
     ELSE IF (JCLASS==NVT_TRBE) THEN
        IP%XFZERO (:,1) = IP%XFZERO (:,1) + XFZEROTROP(IRAD) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
     ELSE
        IP%XFZERO (:,1) = IP%XFZERO (:,1) + XFZERO2 (ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
     ENDIF
  ENDIF
  !
  IP%XEPSO  (:,1) = IP%XEPSO  (:,1) + XEPSO  (ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IP%XGAMM  (:,1) = IP%XGAMM  (:,1) + XGAMM  (ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IP%XQDGAMM(:,1) = IP%XQDGAMM(:,1) + XQDGAMM(ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IP%XQDGMES(:,1) = IP%XQDGMES(:,1) + XQDGMES(ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IP%XT1GMES(:,1) = IP%XT1GMES(:,1) + XT1GMES(ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IP%XT2GMES(:,1) = IP%XT2GMES(:,1) + XT2GMES(ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IP%XQDAMAX(:,1) = IP%XQDAMAX(:,1) + XQDAMAX(ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IP%XT1AMAX(:,1) = IP%XT1AMAX(:,1) + XT1AMAX(ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IP%XT2AMAX(:,1) = IP%XT2AMAX(:,1) + XT2AMAX(ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IP%XAH    (:,1) = IP%XAH    (:,1) + XAH    (ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IP%XBH    (:,1) = IP%XBH    (:,1) + XBH    (ICO2TYPE) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  !
  IF(IO%LAGRI_TO_GRASS.AND.(JCLASS==NVT_C3 .OR. JCLASS==NVT_C4 .OR. JCLASS==NVT_IRR))THEN
    ICLASS=NVT_GRAS
  ELSE
    ICLASS=JCLASS
  ENDIF    
  !
  IP%XTAU_WOOD(:,1) = IP%XTAU_WOOD(:,1) + XTAU_WOOD(ICLASS) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  IP%XAMAX    (:,1) = IP%XAMAX    (:,1) + XAMAX    (ICLASS) * IP%XVEGTYPE_PATCH(:,JCLASS,1)
  !
END DO
!
IP%XQDGAMM(:,1)=LOG(IP%XQDGAMM(:,1))
IP%XQDGMES(:,1)=LOG(IP%XQDGMES(:,1))
IP%XQDAMAX(:,1)=LOG(IP%XQDAMAX(:,1))
!
!
! INITIALIZE VARIOUS VARIABLES FOR CO2 MODEL:
! -------------------------------------------
!
!
! compute temperature responses:
!
!before optimization (with non log IP%XQDGAMM) : 
!ZGAMMT(:) = IP%XGAMM(:)*(IP%XQDGAMM(:)**(0.1*(ZTOPT(:)-25.0)))
ZWORK (:) = (0.1*(ZTOPT(:)-25.0)) * IP%XQDGAMM(:,1)
ZGAMMT(:) = IP%XGAMM(:,1)*EXP(ZWORK(:))
!
!before optimization (with non log IP%XQDAMAX) :
!ZANMAX(:) = ( IP%XAMAX(:)*IP%XQDAMAX(:)**(0.1*(ZTOPT(:)-25.0)) ) / ...
ZWORK (:) = (0.1*(ZTOPT(:)-25.0)) * IP%XQDAMAX(:,1)
ZANMAX(:) = ( IP%XAMAX(:,1)*EXP(ZWORK(:)) )                   &
               /( (1.0+EXP(0.3*(IP%XT1AMAX(:,1)-ZTOPT(:))))*  &
                  (1.0+EXP(0.3*(ZTOPT(:)-IP%XT2AMAX(:,1)))) )  
!
!before optimization (with non log IP%XQDGMES) :
!ZGMEST(:) = ( PGMES(:)*IP%XQDGMES(:)**(0.1*(ZTOPT(:)-25.0)) )    &
ZWORK (:) = (0.1*(ZTOPT(:)-25.0)) * IP%XQDGMES(:,1)
ZGMEST(:) = ( PGMES(:)*EXP(ZWORK(:)) )                   &
               /( (1.0+EXP(0.3*(IP%XT1GMES(:,1)-ZTOPT(:))))*  &
                  (1.0+EXP(0.3*(ZTOPT(:)-IP%XT2GMES(:,1)))) )  
!
!
! initialize other variables: (using optimum values for some variables)
!
ZCO2INIT3(:) = XDSPOPT
ZCO2INIT4(:) = XIAOPT
ZCO2INIT5(:) = 1.0
!
! Define a dummy LAI from top (zco2init2=0.1) and total lai (zco2init=1) for Dark respiration extinction parameterization 
!
ZCO2INIT2(:) = 0.1
ZCO2INIT1(:) = 1.0
!
! Add soil moisture stress effect to leaf conductance:
!
ZGMEST(:) = ZGMEST(:)*ZCO2INIT5(:)
!
! Initialise DMAX following Calvet (2000) in the case of 'AST' or 'LST' photosynthesis option
!
IF((IO%CPHOTO=='AST').OR.(IO%CPHOTO=='LST').OR.(IO%CPHOTO=='NIT').OR.(IO%CPHOTO=='NCB')) THEN
   ZDMAX(:) = EXP((LOG(ZGMEST(:)*1000.)-IP%XAH(:,1))/IP%XBH(:,1))/1000.
ELSE
   ZDMAX(:) = PDMAX(:)
ENDIF
!
! Compute maximum/initial/optimum net assimilation of CO2:
!
! Unit conversion with a constant value of 1.2 for PRHOA as it is not known here
! ZANMAX and ZEPSO from kgCO2/m2/s to kgCO2/kgair m/s by dividing by RHOA (kgair/m3)
! ZGAMMT from ppm to kgCO2/kgair
ZANMAX(:)=ZANMAX(:)/1.2
ZEPSO(:)=IP%XEPSO(:,1)/1.2
ZGAMMT(:)=ZGAMMT(:)*XMCO2/XMD*1e-6
!
CALL COTWO(PCO2, ZCO2INIT5, ZCO2INIT4, ZCO2INIT3, ZGAMMT, &
           IP%XFZERO(:,1), ZEPSO, ZANMAX, ZGMEST, PGC, ZDMAX,     &
           IP%XANMAX(:,1), ZGS, ZRDK, ZCO2INIT2, ZCO2INIT1        )                     
! change by sebastien IP%XEPSO change into ZEPSO for units consistency
!
!
!
IF (LHOOK) CALL DR_HOOK('COTWOINIT_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE COTWOINIT_n
