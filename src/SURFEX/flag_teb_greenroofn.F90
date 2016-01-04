!     #########
      SUBROUTINE FLAG_TEB_GREENROOF_n (TGRR, TGRO, TGRPE, T,  &
                                       KFLAG)
!     ##################################
!
!!****  *FLAG_TEB_GREENROOF_n* - routine to flag ISBA variables where green roofs are
!!                            not present
!!                         
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
!!      Original    10/2011
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_PROG_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_PARAM_TIME_t
USE MODD_TEB_n, ONLY : TEB_t
!
USE MODD_CO2V_PAR,       ONLY : XANFMINIT, XCONDCTMIN
!                                
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODI_FLAG_GR_SNOW
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
TYPE(TEB_VEG_PROG_t), INTENT(INOUT) :: TGRR
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_VEG_PARAM_TIME_t), INTENT(INOUT) :: TGRPE
TYPE(TEB_t), INTENT(INOUT) :: T
!
INTEGER, INTENT(IN) :: KFLAG ! 1 : to put physical values to run ISBA afterwards
!                            ! 2 : to flag with XUNDEF value for points without green roof
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL :: ZWR, ZTG, ZWG, ZRESA, ZANFM, ZDEF
INTEGER :: JL1, JL2 ! loop counter on layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('FLAG_TEB_GREENROOF_N',0,ZHOOK_HANDLE)
!
ZWR = XUNDEF
!
IF (KFLAG==1) THEN
  ZTG   = 300.
  ZWG   = 0.5
  ZRESA = 100.
  ZANFM = XANFMINIT
  ZDEF  = 0.
ELSEIF (KFLAG==2) THEN
  ZTG   = XUNDEF
  ZWG   = XUNDEF
  ZRESA = XUNDEF
  ZANFM = XUNDEF
  ZDEF  = XUNDEF
ENDIF
!
!-------------------------------------------------------------------------------
!     
  !
  DO JL1=1,TGRO%NGROUND_LAYER
    WHERE (T%CUR%XGREENROOF(:)==0.) 
      TGRR%CUR%XTG (:,JL1,1) = ZTG
      TGRR%CUR%XWG (:,JL1,1) = ZWG
      TGRR%CUR%XWGI(:,JL1,1) = ZDEF
    END WHERE
  END DO
  !
  WHERE (T%CUR%XGREENROOF(:)==0.) 
    TGRR%CUR%XWR  (:,1) = ZWR
    TGRR%CUR%XRESA(:,1) = ZRESA
  END WHERE
  !
  IF (TGRO%CPHOTO/='NON') THEN
    !
    WHERE (T%CUR%XGREENROOF(:)==0.)
      TGRR%CUR%XANFM (:,1) = ZANFM              
      TGRR%CUR%XAN   (:,1) = ZDEF
      TGRR%CUR%XANDAY(:,1) = ZDEF
      TGRR%CUR%XLE   (:,1) = ZDEF
    END WHERE
    !
  ELSE IF (TGRO%CPHOTO=='LAI' .OR. TGRO%CPHOTO=='LST' .OR. TGRO%CPHOTO=='NIT' .OR. TGRO%CPHOTO=='NCB') THEN
    !
    WHERE (T%CUR%XGREENROOF(:)==0.) TGRPE%CUR%XLAI(:,1) = ZDEF
    !
  ELSE IF (TGRO%CPHOTO=='AGS' .OR. TGRO%CPHOTO=='AST') THEN
    !
    DO JL1=1,SIZE(TGRR%CUR%XBIOMASS,2)
      WHERE (T%CUR%XGREENROOF(:)==0.)
        TGRR%CUR%XBIOMASS     (:,JL1,1) = ZDEF
        TGRR%CUR%XRESP_BIOMASS(:,JL1,1) = ZDEF
      END WHERE
    END DO
      !
  END IF
    !
!ENDIF
  !
!
!-------------------------------------------------------------------------------
!
!* Flag snow characteristics
!
 CALL FLAG_GR_SNOW(KFLAG,T%CUR%XGREENROOF(:)==0.,TGRR%CUR%TSNOW)
!
!
!* snow-free characteristics
!
IF (KFLAG==1) THEN
  WHERE (T%CUR%XGREENROOF==0.) TGRR%CUR%XSNOWFREE_ALB(:,1)      = 0.2
  WHERE (T%CUR%XGREENROOF==0.) TGRR%CUR%XSNOWFREE_ALB_VEG(:,1)  = 0.2
  WHERE (T%CUR%XGREENROOF==0.) TGRR%CUR%XSNOWFREE_ALB_SOIL(:,1) = 0.2
ELSEIF (KFLAG==2) THEN
  WHERE (T%CUR%XGREENROOF==0.) TGRR%CUR%XSNOWFREE_ALB(:,1)      = XUNDEF
  WHERE (T%CUR%XGREENROOF==0.) TGRR%CUR%XSNOWFREE_ALB_VEG(:,1)  = XUNDEF
  WHERE (T%CUR%XGREENROOF==0.) TGRR%CUR%XSNOWFREE_ALB_SOIL(:,1) = XUNDEF
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FLAG_TEB_GREENROOF_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FLAG_TEB_GREENROOF_n
