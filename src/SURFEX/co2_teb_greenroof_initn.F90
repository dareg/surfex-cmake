!     #########
      SUBROUTINE CO2_TEB_GREENROOF_INIT_n (TGRO, TGRIP, TGRM, &
                                           PCO2)
!     #####################
!
!!****  *CO2_TEB_GREENROOF_INIT_n* - routine to initialize ISBA-AGS variables
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
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
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2003 
!!      J.C. Calvet 01/2004 Externalization
!!      P Le Moigne 11/2004 cotwoinit changed into cotwoinit_n
!!      P Le Moigne 09/2005 AGS modifs of L. Jarlan
!!      S Lafont    09/2008 Add initialisation of POI and ABC (needed for TORI)
!!      A.L. Gibelin 04/2009 : TAU_WOOD for NCB option 
!!      A.L. Gibelin 04/2009 : Add carbon spinup
!!      A.L. Gibelin 07/2009 : Suppress RDK and transform GPP as a diagnostic
!!      A.L. Gibelin 07/2009 : Suppress PPST and PPSTF as outputs
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_TEB_VEG_PARAM_n, ONLY : TEB_VEG_PARAM_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
!
USE MODI_COTWOINIT_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: TGRIP
TYPE(TEB_VEG_PARAM_t), INTENT(INOUT) :: TGRM
!
REAL, DIMENSION(:), INTENT(IN) :: PCO2 ! air CO2 concentration (kg/kg)
!
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(SIZE(TGRM%X%XVEGTYPE,1)) :: ZTAU_WOOD
INTEGER :: ILU   ! size of arrays
INTEGER :: JP    ! loop on tiles
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CO2_TEB_GREENROOF_INIT_N',0,ZHOOK_HANDLE)
ILU = SIZE(TGRM%X%XVEGTYPE,1)
!
ALLOCATE(TGRIP%XANMAX        (ILU,1))
ALLOCATE(TGRIP%XFZERO        (ILU,1))
ALLOCATE(TGRIP%XEPSO         (ILU,1))
ALLOCATE(TGRIP%XGAMM         (ILU,1))
ALLOCATE(TGRIP%XQDGAMM       (ILU,1))
ALLOCATE(TGRIP%XQDGMES       (ILU,1))
ALLOCATE(TGRIP%XT1GMES       (ILU,1))
ALLOCATE(TGRIP%XT2GMES       (ILU,1))
ALLOCATE(TGRIP%XAMAX         (ILU,1))
ALLOCATE(TGRIP%XQDAMAX       (ILU,1))
ALLOCATE(TGRIP%XT1AMAX       (ILU,1))
ALLOCATE(TGRIP%XT2AMAX       (ILU,1))
ALLOCATE(TGRIP%XAH           (ILU,1))
ALLOCATE(TGRIP%XBH           (ILU,1))
!
     CALL COTWOINIT_n(TGRO%LAGRI_TO_GRASS, TGRO%LTR_ML, &
                      TGRO%CPHOTO, TGRM%X%XVEGTYPE,TGRM%T%CUR%XGMES(:,1),&
                      PCO2,TGRM%T%CUR%XGC(:,1),TGRM%X%XDMAX(:,1),TGRIP%XABC(:),&
                      TGRIP%XPOI(:),TGRIP%XANMAX(:,1), TGRIP%XFZERO(:,1),      &
                      TGRIP%XEPSO(:,1),TGRIP%XGAMM(:,1),TGRIP%XQDGAMM(:,1), &
                      TGRIP%XQDGMES(:,1),TGRIP%XT1GMES(:,1), TGRIP%XT2GMES(:,1), &
                      TGRIP%XAMAX(:,1),TGRIP%XQDAMAX(:,1),TGRIP%XT1AMAX(:,1),    &
                      TGRIP%XT2AMAX(:,1),TGRIP%XAH(:,1),TGRIP%XBH(:,1),ZTAU_WOOD        )  
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CO2_TEB_GREENROOF_INIT_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE CO2_TEB_GREENROOF_INIT_n
