MODULE MODE_EGGANGLES
!
! Version 2006.0614 by JD GRIL
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DOC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  All these functions make a package tool for angle.
!  In functions where appears DOM and UNIT indicate the domain of validity :
!         DOM     UNIT           Longitudes           Latitudes
!         "-+"    "D"          [-180.0,180.0[       [-90.0,90.0]
!         "0+"    "D"             [0.0,360.0[       [-90.0,90.0]
!         "-+"    "R"             [-pi,pi[        [-pi/2.0,pi/2.0] 
!         "0+"    "R"              [0,pi[         [-pi/2.0,pi/2.0]
!  (defaults values are DOM = "-+" and UNIT = "D").

!  All functions work for scalar or one dimensional array in input.

!  -1- ANGLE_DOMAIN function

!->function ANGLE_DOMAIN(ALPHA,PI,DOM,UNIT)

!     Converts longitudes in UNIT values under choisen DOMain.
!     The input (ALPHA) is a longitude (REAL) or a LOLA type structure ( or
!     array of them). The output has the same type than the input.       

!  -2- VAL_ functions

!->integer function VAL_LAT(LAT,NUM_ERR,PI,UNIT)

!     Test validity of LAT [-90.0,90.0] 
!     Return -1 or NUM_ERR if it's present in error case, 1 if it's ok.

!->integer function VAL_LON(LON,NUM_ERR,PI,DOM,UNIT)

!     Test validity of LON [-180.0,180.0[ or [0.0,360.0[ 
!     Return -1 or NUM_ERR if it's present in error case, 1 if it's ok.

!->integer function VAL_COORD(PT_COORD,NUM_ERR,PI,DOM,UNIT) 

!     Test validity of LAT [-90.0,90.0] and LON [-180.0,180.0[ or [0.0,360.0[
!     (depends the value of DOM) of a PT_COORD structure of type LOLA (in UNIT).
!     Return -1 or NUM_ERR if it's present in error case, 1 if it's ok.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author : Jean-Daniel GRIL , CNRM/GMAP/COOPE , Februry 08 2000

! Modified:
! In April 2001 by M. Janousek (A few modifs to port the deck to the model code)
! In November 2004 by JD Gril : more routines to manage angles
!                             : debug VAL_COORD_x
! 2005 by JD Gril : more functions for Mercator RT
! In June 2006 by JD Gril     : line too long (L607 > 132 col.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! ******************* Definition of parameters **********************************
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
  USE PARKIND1  ,ONLY : JPIM,    JPRB
!
  IMPLICIT NONE
!
! Include Kinds
! -------------
!
!* kindef: define default KIND macros
! --------------------------------------
!INTEGER, PARAMETER :: JPIM=4
!INTEGER, PARAMETER :: JPRB=8
!USE PARKIND1  ,ONLY : JPIM,    JPRB
! --------------------------------------
!
!
! ******************* Definition of type ****************************************
!
  TYPE LOLA
     SEQUENCE
     REAL(KIND=JPRB) :: LON, LAT
  END TYPE LOLA
!
! ******************* Definition of Interface ***********************************
!
  INTERFACE ANGLE_DOMAIN
     MODULE PROCEDURE ANGLE_DOMAIN_RS, ANGLE_DOMAIN_LOLAS, ANGLE_DOMAIN_RV, ANGLE_DOMAIN_LOLAV
  END INTERFACE
  INTERFACE VAL_LAT
     MODULE PROCEDURE VAL_LAT_S, VAL_LAT_V
  END INTERFACE
  INTERFACE VAL_LON
     MODULE PROCEDURE VAL_LON_S, VAL_LON_V
  END INTERFACE
  INTERFACE VAL_COORD
     MODULE PROCEDURE VAL_COORD_S, VAL_COORD_V
  END INTERFACE
  INTERFACE COSIN_TO_ANGLE
     MODULE PROCEDURE COSIN_TO_ANGLE_S, COSIN_TO_ANGLE_V
  END INTERFACE
  INTERFACE P_ASIN
     MODULE PROCEDURE P_ASIN_S, P_ASIN_V
  END INTERFACE
  INTERFACE P_ACOS
     MODULE PROCEDURE P_ACOS_S, P_ACOS_V
  END INTERFACE
  INTERFACE MINIMAX
     MODULE PROCEDURE MINIMAX_S, MINIMAX_V
  END INTERFACE
CONTAINS

! =================== FUNCTIONS =================================================

! ******************* Independants functions ************************************
!
! -------------------------------------------------------------------------------
  FUNCTION ANGLE_DOMAIN_RS(ALPHA,PI,DOM,UNIT) RESULT (BETA)
    REAL(KIND=JPRB), INTENT(IN)                           :: ALPHA
    CHARACTER (LEN=2), INTENT(IN), OPTIONAL               :: DOM
    CHARACTER (LEN=1), INTENT(IN), OPTIONAL               :: UNIT
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL                 :: PI
    REAL(KIND=JPRB)                                       :: BETA

    REAL(KIND=JPRB)   :: CVT, TPI, M
    CHARACTER (LEN=2) :: TDOM
    CHARACTER (LEN=1) :: TUNIT
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:ANGLE_DOMAIN_RS',0,ZHOOK_HANDLE)
    IF (PRESENT(PI)) THEN
       TPI = PI
    ELSE
       TPI = ASIN(1.0_JPRB)*2.0_JPRB
    ENDIF
    IF (PRESENT(DOM)) THEN
       IF ((DOM=='0+').OR.(DOM=='-+')) THEN
          TDOM = DOM
       ELSE
          TDOM = "-+"
       ENDIF
    ELSE
       TDOM = "-+"
    ENDIF
    IF (PRESENT(UNIT)) THEN
       IF ((UNIT=='R').OR.(UNIT=='D')) THEN
          TUNIT = UNIT
       ELSE
          TUNIT = "D"
       ENDIF
    ELSE
       TUNIT = "D"
    ENDIF

    IF (TUNIT=='R') THEN
       CVT = TPI
    ELSE
       CVT = 180.0_JPRB
    ENDIF

    IF (TDOM=='-+') THEN
       M = MOD(ALPHA,CVT)
       BETA = (M-CVT*MOD(REAL(INT(ALPHA/CVT),KIND=JPRB),2.0_JPRB))*SIGN(1.0_JPRB,ALPHA)*SIGN(1.0_JPRB,M)
    ELSE
       M = MOD(ALPHA,2.0_JPRB*CVT)
       BETA = M-2.0_JPRB*CVT*(SIGN(0.5_JPRB,ALPHA)-0.5_JPRB)
    ENDIF
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:ANGLE_DOMAIN_RS',1,ZHOOK_HANDLE)
  END FUNCTION ANGLE_DOMAIN_RS
! -------------------------------------------------------------------------------
   FUNCTION ANGLE_DOMAIN_LOLAS(ALPHA,PI,DOM,UNIT) RESULT (BETA)
    TYPE (LOLA), INTENT(IN)                                :: ALPHA
    CHARACTER (LEN=2), INTENT(IN), OPTIONAL                :: DOM
    CHARACTER (LEN=1), INTENT(IN), OPTIONAL                :: UNIT
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL                  :: PI
    TYPE (LOLA)                                            :: BETA

    REAL(KIND=JPRB)   :: TPI
    CHARACTER (LEN=2) :: TDOM
    CHARACTER (LEN=1) :: TUNIT
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:ANGLE_DOMAIN_LOLAS',0,ZHOOK_HANDLE)
    IF (PRESENT(PI)) THEN
       TPI = PI
    ELSE
       TPI = ASIN(1.0_JPRB)*2.0_JPRB
    ENDIF
    IF (PRESENT(DOM)) THEN
       IF ((DOM=='0+').OR.(DOM=='-+')) THEN
          TDOM = DOM
       ELSE
          TDOM = "-+"
       ENDIF
    ELSE
       TDOM = "-+"
    ENDIF
    IF (PRESENT(UNIT)) THEN
       IF ((UNIT=='R').OR.(UNIT=='D')) THEN
          TUNIT = UNIT
       ELSE
          TUNIT = "D"
       ENDIF
    ELSE
       TUNIT = "D"
    ENDIF

    BETA%LON = ANGLE_DOMAIN(ALPHA%LON,TPI,TDOM,TUNIT)
    BETA%LAT = ALPHA%LAT
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:ANGLE_DOMAIN_LOLAS',1,ZHOOK_HANDLE)
  END FUNCTION ANGLE_DOMAIN_LOLAS
! -------------------------------------------------------------------------------
  FUNCTION ANGLE_DOMAIN_RV(ALPHA,PI,DOM,UNIT) RESULT (BETA)
    REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)             :: ALPHA
    CHARACTER (LEN=2), INTENT(IN), OPTIONAL               :: DOM
    CHARACTER (LEN=1), INTENT(IN), OPTIONAL               :: UNIT
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL                 :: PI
    REAL(KIND=JPRB), DIMENSION(SIZE(ALPHA))               :: BETA

    REAL(KIND=JPRB)                         :: CVT, TPI
    REAL(KIND=JPRB), DIMENSION(SIZE(ALPHA)) :: Z_M
    CHARACTER (LEN=2)                       :: TDOM
    CHARACTER (LEN=1)                       :: TUNIT
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:ANGLE_DOMAIN_RV',0,ZHOOK_HANDLE)
    IF (PRESENT(PI)) THEN
       TPI = PI
    ELSE
       TPI = ASIN(1.0_JPRB)*2.0_JPRB
    ENDIF
    IF (PRESENT(DOM)) THEN
       IF ((DOM=='0+').OR.(DOM=='-+')) THEN
          TDOM = DOM
       ELSE
          TDOM = "-+"
       ENDIF
    ELSE
       TDOM = "-+"
    ENDIF
    IF (PRESENT(UNIT)) THEN
       IF ((UNIT=='R').OR.(UNIT=='D')) THEN
          TUNIT = UNIT
       ELSE
          TUNIT = "D"
       ENDIF
    ELSE
       TUNIT = "D"
    ENDIF

    IF (TUNIT=='R') THEN
       CVT = TPI
    ELSE
       CVT = 180.0_JPRB
    ENDIF

    IF (TDOM=='-+') THEN
       Z_M(:) = MOD(ALPHA(:),CVT)
       BETA = (Z_M(:)-CVT*MOD(REAL(INT(ALPHA(:)/CVT),KIND=JPRB),2.0_JPRB))*SIGN(1.0_JPRB,ALPHA(:))*SIGN(1.0_JPRB,Z_M(:))
    ELSE
       Z_M(:) = MOD(ALPHA(:),2.0_JPRB*CVT)
       BETA = Z_M(:)-2.0_JPRB*CVT*(SIGN(0.5_JPRB,ALPHA(:))-0.5_JPRB)
    ENDIF
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:ANGLE_DOMAIN_RV',1,ZHOOK_HANDLE)
  END FUNCTION ANGLE_DOMAIN_RV
! -------------------------------------------------------------------------------
  FUNCTION ANGLE_DOMAIN_LOLAV(YL_ALPHA,PI,DOM,UNIT) RESULT (YD_BETA)
    TYPE (LOLA), DIMENSION(:), INTENT(IN)               :: YL_ALPHA
    CHARACTER (LEN=2), INTENT(IN), OPTIONAL             :: DOM
    CHARACTER (LEN=1), INTENT(IN), OPTIONAL             :: UNIT
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL               :: PI
    TYPE (LOLA), DIMENSION(SIZE(YL_ALPHA))              :: YD_BETA

    REAL(KIND=JPRB)   :: TPI
    CHARACTER (LEN=2) :: TDOM
    CHARACTER (LEN=1) :: TUNIT
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:ANGLE_DOMAIN_LOLAV',0,ZHOOK_HANDLE)
    IF (PRESENT(PI)) THEN
       TPI = PI
    ELSE
       TPI = ASIN(1.0_JPRB)*2.0_JPRB
    ENDIF
    IF (PRESENT(DOM)) THEN
       IF ((DOM=='0+').OR.(DOM=='-+')) THEN
          TDOM = DOM
       ELSE
          TDOM = "-+"
       ENDIF
    ELSE
       TDOM = "-+"
    ENDIF
    IF (PRESENT(UNIT)) THEN
       IF ((UNIT=='R').OR.(UNIT=='D')) THEN
          TUNIT = UNIT
       ELSE
          TUNIT = "D"
       ENDIF
    ELSE
       TUNIT = "D"
    ENDIF

    YD_BETA(:)%LON = ANGLE_DOMAIN(YL_ALPHA(:)%LON,TPI,TDOM,TUNIT)
    YD_BETA(:)%LAT = YL_ALPHA(:)%LAT
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:ANGLE_DOMAIN_LOLAV',1,ZHOOK_HANDLE)
  END FUNCTION ANGLE_DOMAIN_LOLAV
! -------------------------------------------------------------------------------
  FUNCTION VAL_LAT_S(LAT,NUM_ERR,PI,UNIT) RESULT(ETAT)
    REAL(KIND=JPRB), INTENT(IN)                          :: LAT
    CHARACTER (LEN=1), INTENT(IN), OPTIONAL              :: UNIT
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL                :: PI
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL             :: NUM_ERR
    INTEGER(KIND=JPIM)                                   :: ETAT

    INTEGER(KIND=JPIM) :: TNE
    REAL(KIND=JPRB)    :: TPI, LATMXABS
    CHARACTER (LEN=1)  :: TUNIT
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_LAT_S',0,ZHOOK_HANDLE)
    IF (PRESENT(NUM_ERR))THEN
       TNE = NUM_ERR
    ELSE
       TNE = -1_JPIM
    ENDIF
    IF (PRESENT(PI)) THEN
       TPI = PI
    ELSE
       TPI = ASIN(1.0_JPRB)*2.0_JPRB
    ENDIF
    IF (PRESENT(UNIT)) THEN
       IF ((UNIT=='R').OR.(UNIT=='D')) THEN
          TUNIT = UNIT
       ELSE
          TUNIT = "D"
       ENDIF
    ELSE
       TUNIT = "D"
    ENDIF

    IF (TUNIT=='R') THEN
       LATMXABS = TPI/2.0_JPRB
    ELSE
       LATMXABS = 90.0_JPRB
    ENDIF

    IF (ABS(LAT) > LATMXABS) THEN
       ETAT = TNE
    ELSE
       ETAT = 1_JPIM
    ENDIF
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_LAT_S',1,ZHOOK_HANDLE)
  END FUNCTION VAL_LAT_S
! -------------------------------------------------------------------------------
  FUNCTION VAL_LAT_V(P_LAT,NUM_ERR,PI,UNIT) RESULT(ETAT)
    REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)                 :: P_LAT
    CHARACTER (LEN=1), INTENT(IN), OPTIONAL                   :: UNIT
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL                     :: PI
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL                  :: NUM_ERR
    INTEGER(KIND=JPIM)                                        :: ETAT

    INTEGER(KIND=JPIM) :: TNE
    REAL(KIND=JPRB)    :: TPI, Z_LATMXABS
    CHARACTER (LEN=1)  :: TUNIT
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_LAT_V',0,ZHOOK_HANDLE)
    IF (PRESENT(NUM_ERR))THEN
       TNE = NUM_ERR
    ELSE
       TNE = -1_JPIM
    ENDIF
    IF (PRESENT(PI)) THEN
       TPI = PI
    ELSE
       TPI = ASIN(1.0_JPRB)*2.0_JPRB
    ENDIF
    IF (PRESENT(UNIT)) THEN
       IF ((UNIT=='R').OR.(UNIT=='D')) THEN
          TUNIT = UNIT
       ELSE
          TUNIT = "D"
       ENDIF
    ELSE
       TUNIT = "D"
    ENDIF

    IF (TUNIT=='R') THEN
       Z_LATMXABS = TPI/2.0_JPRB
    ELSE
       Z_LATMXABS = 90.0_JPRB
    ENDIF

    IF (ANY(ABS(P_LAT(:)) > Z_LATMXABS)) THEN
       ETAT = TNE
    ELSE
       ETAT = 1_JPIM
    ENDIF
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_LAT_V',1,ZHOOK_HANDLE)
  END FUNCTION VAL_LAT_V
! -------------------------------------------------------------------------------
  FUNCTION VAL_LON_S(LON,NUM_ERR,PI,DOM,UNIT) RESULT(ETAT)
    REAL(KIND=JPRB), INTENT(IN)                                :: LON
    CHARACTER (LEN=2), INTENT(IN), OPTIONAL                    :: DOM
    CHARACTER (LEN=1), INTENT(IN), OPTIONAL                    :: UNIT
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL                      :: PI
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL                   :: NUM_ERR
    INTEGER(KIND=JPIM)                                         :: ETAT

    INTEGER(KIND=JPIM) :: TNE
    REAL(KIND=JPRB)    :: TPI, CVT, S, LONMIN, LONMAX
    CHARACTER (LEN=2)  :: TDOM
    CHARACTER (LEN=1)  :: TUNIT
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_LON_S',0,ZHOOK_HANDLE)
    IF (PRESENT(NUM_ERR))THEN
       TNE = NUM_ERR
    ELSE
       TNE = -1_JPIM
    ENDIF
    IF (PRESENT(PI)) THEN
       TPI = PI
    ELSE
       TPI = ASIN(1.0_JPRB)*2.0_JPRB
    ENDIF
    IF (PRESENT(DOM)) THEN
       IF ((DOM=='0+').OR.(DOM=='-+')) THEN
          TDOM = DOM
       ELSE
          TDOM = "-+"
       ENDIF
    ELSE
       TDOM = "-+"
    ENDIF
    IF (PRESENT(UNIT)) THEN
       IF ((UNIT=='R').OR.(UNIT=='D')) THEN
          TUNIT = UNIT
       ELSE
          TUNIT = "D"
       ENDIF
    ELSE
       TUNIT = "D"
    ENDIF

    IF (TUNIT=='R') THEN
       CVT = TPI
    ELSE
       CVT = 180.0_JPRB
    ENDIF
    IF (TDOM=='-+') THEN
       S = -1.0_JPRB
    ELSE
       S = 0.0_JPRB
    ENDIF
    LONMIN = S*CVT
    LONMAX =(2.0_JPRB +S)*CVT

    IF ((LON < LONMIN).OR.(LON >= LONMAX)) THEN
       ETAT = TNE
    ELSE
       ETAT = 1_JPIM
    ENDIF
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_LON_S',1,ZHOOK_HANDLE)
  END FUNCTION VAL_LON_S
! -------------------------------------------------------------------------------
  FUNCTION VAL_LON_V(LON,NUM_ERR,PI,DOM,UNIT) RESULT(ETAT)
    REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)                       :: LON
    CHARACTER (LEN=2), INTENT(IN), OPTIONAL                         :: DOM
    CHARACTER (LEN=1), INTENT(IN), OPTIONAL                         :: UNIT
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL                           :: PI
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL                        :: NUM_ERR
    INTEGER(KIND=JPIM)                                              :: ETAT
 
    INTEGER(KIND=JPIM) :: TNE
    REAL(KIND=JPRB)    :: TPI, Z_CVT, Z_S, Z_LONMIN, Z_LONMAX
    CHARACTER (LEN=2)  :: TDOM
    CHARACTER (LEN=1)  :: TUNIT
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_LON_V',0,ZHOOK_HANDLE)
    IF (PRESENT(NUM_ERR))THEN
       TNE = NUM_ERR
    ELSE
       TNE = -1_JPIM
    ENDIF
    IF (PRESENT(PI)) THEN
       TPI = PI
    ELSE
       TPI = ASIN(1.0_JPRB)*2.0_JPRB
    ENDIF
    IF (PRESENT(DOM)) THEN
       IF ((DOM=='0+').OR.(DOM=='-+')) THEN
          TDOM = DOM
       ELSE
          TDOM = "-+"
       ENDIF
    ELSE
       TDOM = "-+"
    ENDIF
    IF (PRESENT(UNIT)) THEN
       IF ((UNIT=='R').OR.(UNIT=='D')) THEN
          TUNIT = UNIT
       ELSE
          TUNIT = "D"
       ENDIF
    ELSE
       TUNIT = "D"
    ENDIF

    IF (TUNIT=='R') THEN
       Z_CVT = TPI
    ELSE
       Z_CVT = 180.0_JPRB
    ENDIF
    IF (TDOM=='-+') THEN
       Z_S = -1.0_JPRB
    ELSE
       Z_S = 0.0_JPRB
    ENDIF
    Z_LONMIN = Z_S*Z_CVT
    Z_LONMAX =(2.0_JPRB +Z_S)*Z_CVT

    IF ((ANY(LON(:) < Z_LONMIN)).OR.(ANY(LON(:) >= Z_LONMAX))) THEN
       ETAT = TNE
    ELSE
       ETAT = 1_JPIM
    ENDIF
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_LON_V',1,ZHOOK_HANDLE)
  END FUNCTION VAL_LON_V
! -------------------------------------------------------------------------------
  FUNCTION VAL_COORD_S(PT_COORD,NUM_ERR,PI,DOM,UNIT) RESULT(ETAT)
    TYPE (LOLA), INTENT(IN)                               :: PT_COORD
    CHARACTER (LEN=2), INTENT(IN), OPTIONAL               :: DOM
    CHARACTER (LEN=1), INTENT(IN), OPTIONAL               :: UNIT
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL                 :: PI
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL              :: NUM_ERR
    INTEGER(KIND=JPIM)                                    :: ETAT

    INTEGER(KIND=JPIM) :: TNE
    REAL(KIND=JPRB)    :: TPI
    CHARACTER (LEN=2)  :: TDOM
    CHARACTER (LEN=1)  :: TUNIT
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_COORD_S',0,ZHOOK_HANDLE)
    IF (PRESENT(NUM_ERR))THEN
       TNE = NUM_ERR
    ELSE
       TNE = -1_JPIM
    ENDIF
    IF (PRESENT(PI)) THEN
       TPI = PI
    ELSE
       TPI = ASIN(1.0_JPRB)*2.0_JPRB
    ENDIF
    IF (PRESENT(DOM)) THEN
       IF ((DOM=='0+').OR.(DOM=='-+')) THEN
          TDOM = DOM
       ELSE
          TDOM = "-+"
       ENDIF
    ELSE
       TDOM = "-+"
    ENDIF
    IF (PRESENT(UNIT)) THEN
       IF ((UNIT=='R').OR.(UNIT=='D')) THEN
          TUNIT = UNIT
       ELSE
          TUNIT = "D"
       ENDIF
    ELSE
       TUNIT = "D"
    ENDIF

    IF ((VAL_LON(PT_COORD%LON,TNE,TPI,TDOM,TUNIT) == 1_JPIM).AND.(VAL_LAT(PT_COORD%LAT,TNE,TPI,TUNIT) == 1_JPIM)) THEN
       ETAT = 1_JPIM
    ELSE
       ETAT = TNE
    ENDIF
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_COORD_S',1,ZHOOK_HANDLE)
  END FUNCTION VAL_COORD_S
! -------------------------------------------------------------------------------
  FUNCTION VAL_COORD_V(YD_PT_COORD,K_NUM_ERR,PI,CD_DOM,CD_UNIT) RESULT(ETAT)
    TYPE (LOLA), DIMENSION(:), INTENT(IN)                   :: YD_PT_COORD
    CHARACTER (LEN=2), INTENT(IN), OPTIONAL                 :: CD_DOM
    CHARACTER (LEN=1), INTENT(IN), OPTIONAL                 :: CD_UNIT
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL                   :: PI
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL                :: K_NUM_ERR
    INTEGER(KIND=JPIM)                                      :: ETAT
 
    INTEGER(KIND=JPIM) :: I_TNE
    CHARACTER (LEN=2)  :: CL_TDOM
    REAL(KIND=JPRB)    :: Z_TPI
    CHARACTER (LEN=1)  :: CL_TUNIT
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_COORD_V',0,ZHOOK_HANDLE)
    IF (PRESENT(K_NUM_ERR))THEN
       I_TNE = K_NUM_ERR
    ELSE
       I_TNE = -1_JPIM
    ENDIF
    IF (PRESENT(PI)) THEN
       Z_TPI = PI
    ELSE
       Z_TPI = ASIN(1.0_JPRB)*2.0_JPRB
    ENDIF
    IF (PRESENT(CD_DOM)) THEN
       IF ((CD_DOM=='0+').OR.(CD_DOM=='-+')) THEN
          CL_TDOM = CD_DOM
       ELSE
          CL_TDOM = "-+"
       ENDIF
    ELSE
       CL_TDOM = "-+"
    ENDIF
    IF (PRESENT(CD_UNIT)) THEN
       IF ((CD_UNIT=='R').OR.(CD_UNIT=='D')) THEN
          CL_TUNIT = CD_UNIT
       ELSE
          CL_TUNIT = "D"
       ENDIF
    ELSE
       CL_TUNIT = "D"
    ENDIF

    IF ((VAL_LON(YD_PT_COORD(:)%LON,I_TNE,Z_TPI,CL_TDOM,CL_TUNIT) == 1_JPIM).AND. &
         (VAL_LAT(YD_PT_COORD(:)%LAT,I_TNE,Z_TPI,CL_TUNIT) == 1_JPIM)) THEN 
       ETAT = 1_JPIM
    ELSE
       ETAT = I_TNE
    ENDIF
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:VAL_COORD_V',1,ZHOOK_HANDLE)
  END FUNCTION VAL_COORD_V
! -------------------------------------------------------------------------------

  FUNCTION LOLAR (COORD_DEG) RESULT (COORD_RAD)
    TYPE(LOLA), INTENT(IN)                      :: COORD_DEG
    TYPE(LOLA)                                  :: COORD_RAD

    REAL(KIND=JPRB) :: TPI,DTR
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:LOLAR',0,ZHOOK_HANDLE)
    TPI = ASIN(1.0_JPRB)*2.0_JPRB
    DTR = TPI/180.0_JPRB
    COORD_RAD%LON = COORD_DEG%LON*DTR
    COORD_RAD%LAT = COORD_DEG%LAT*DTR
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:LOLAR',1,ZHOOK_HANDLE)
  END FUNCTION LOLAR

! -------------------------------------------------------------------------------
  FUNCTION LOLAD (COORD_RAD) RESULT (COORD_DEG)
    TYPE(LOLA), INTENT(IN)                      :: COORD_RAD
    TYPE(LOLA)                                  :: COORD_DEG

    REAL(KIND=JPRB) :: TPI,RTD
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:LOLAD',0,ZHOOK_HANDLE)
    TPI = ASIN(1.0_JPRB)*2.0_JPRB
    RTD = 180.0_JPRB/TPI
    COORD_DEG%LON = COORD_RAD%LON*RTD
    COORD_DEG%LAT = COORD_RAD%LAT*RTD
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:LOLAD',1,ZHOOK_HANDLE)
  END FUNCTION LOLAD
  ! -------------------------------------------------------------------------------
  ! Function to compute Cosine,Sine to Angle
  ! -------------------------------------------------------------------------------
  FUNCTION COSIN_TO_ANGLE_S(COSINUS,SINUS) RESULT (ANGLE)
    REAL(KIND=JPRB), INTENT(IN)                  :: COSINUS,SINUS
    REAL(KIND=JPRB)                              :: ANGLE
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
 
    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:COSIN_TO_ANGLE_S',0,ZHOOK_HANDLE)
    ANGLE = P_ACOS(COSINUS)*SIGN(1.0_JPRB,SINUS)
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:COSIN_TO_ANGLE_S',1,ZHOOK_HANDLE)
  END FUNCTION COSIN_TO_ANGLE_S
  ! -------------------------------------------------------------------------------
  ! Function to compute Cosine,Sine to Angle
  ! -------------------------------------------------------------------------------
  FUNCTION COSIN_TO_ANGLE_V(COSINUS,SINUS) RESULT (ANGLE)
    REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)    :: COSINUS,SINUS
    REAL(KIND=JPRB), DIMENSION(SIZE(SINUS))      :: ANGLE
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
 

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:COSIN_TO_ANGLE_V',0,ZHOOK_HANDLE)
    ANGLE = P_ACOS(COSINUS)*SIGN(1.0_JPRB,SINUS)
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:COSIN_TO_ANGLE_V',1,ZHOOK_HANDLE)
  END FUNCTION COSIN_TO_ANGLE_V
  ! -------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------
  ! Function to compute Acos without error
  ! -------------------------------------------------------------------------------
  FUNCTION P_ACOS_S(COSINUS) RESULT (ANGLE)
    REAL(KIND=JPRB), INTENT(IN)                  :: COSINUS
    REAL(KIND=JPRB)                              :: ANGLE
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:P_ACOS_S',0,ZHOOK_HANDLE)
    ANGLE = ACOS(MINIMAX(COSINUS))
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:P_ACOS_S',1,ZHOOK_HANDLE)
  END FUNCTION P_ACOS_S
  ! -------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------
  ! Function to compute Acos without error
  ! -------------------------------------------------------------------------------
  FUNCTION P_ACOS_V(COSINUS) RESULT (ANGLE)
    REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)  :: COSINUS
    REAL(KIND=JPRB), DIMENSION(SIZE(COSINUS))   :: ANGLE
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:P_ACOS_V',0,ZHOOK_HANDLE)
    ANGLE = ACOS(MINIMAX(COSINUS))
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:P_ACOS_V',1,ZHOOK_HANDLE)
  END FUNCTION P_ACOS_V
  ! -------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------
  ! Function to compute Asin without error
  ! -------------------------------------------------------------------------------
  FUNCTION P_ASIN_S(SINUS) RESULT (ANGLE)
    REAL(KIND=JPRB), INTENT(IN)                  :: SINUS
    REAL(KIND=JPRB)                              :: ANGLE
    REAL(KIND=JPRB) :: ZHOOK_HANDLE


    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:P_ASIN_S',0,ZHOOK_HANDLE)
    ANGLE = ASIN(MINIMAX(SINUS))
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:P_ASIN_S',1,ZHOOK_HANDLE)
  END FUNCTION P_ASIN_S
  ! -------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------
  ! Function to compute Asin without error
  ! -------------------------------------------------------------------------------
  FUNCTION P_ASIN_V(SINUS) RESULT (ANGLE)
    REAL(KIND=JPRB), DIMENSION(:), INTENT(IN) :: SINUS
    REAL(KIND=JPRB), DIMENSION(SIZE(SINUS))   :: ANGLE
    REAL(KIND=JPRB) :: ZHOOK_HANDLE


    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:P_ASIN_V',0,ZHOOK_HANDLE)
    ANGLE = ASIN(MINIMAX(SINUS))
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:P_ASIN_V',1,ZHOOK_HANDLE)
  END FUNCTION P_ASIN_V
  ! -------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------
  ! Function MinMax
  ! -------------------------------------------------------------------------------
  FUNCTION MINIMAX_S(VAL,LIM) RESULT (VALO)
    REAL(KIND=JPRB), INTENT(IN)                      :: VAL
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL            :: LIM
    REAL(KIND=JPRB)                                  :: VALO

    REAL(KIND=JPRB) :: TLIM
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:MINIMAX_S',0,ZHOOK_HANDLE)
    IF (PRESENT(LIM)) THEN
        TLIM = LIM
    ELSE
        TLIM = 1.0_JPRB
    ENDIF
    VALO = MIN(TLIM,MAX(-1.0_JPRB*TLIM,VAL))
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:MINIMAX_S',1,ZHOOK_HANDLE)
  END FUNCTION MINIMAX_S
  ! -------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------
  ! Function MinMax
  ! -------------------------------------------------------------------------------
  FUNCTION MINIMAX_V(VAL,LIM) RESULT (VALO)
    REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)        :: VAL
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL            :: LIM
    REAL(KIND=JPRB), DIMENSION(SIZE(VAL))            :: VALO

    REAL(KIND=JPRB) :: TLIM
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:MINIMAX_V',0,ZHOOK_HANDLE)
    IF (PRESENT(LIM)) THEN
        TLIM = LIM
    ELSE
        TLIM = 1.0_JPRB
    ENDIF
    VALO = MIN(TLIM,MAX(-1.0_JPRB*TLIM,VAL))
  IF (LHOOK) CALL DR_HOOK('MODE_EGGANGLES:MINIMAX_V',1,ZHOOK_HANDLE)
  END FUNCTION MINIMAX_V

  ! -------------------------------------------------------------------------------
END MODULE MODE_EGGANGLES
