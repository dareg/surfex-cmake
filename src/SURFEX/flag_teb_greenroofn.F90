!     #########
      SUBROUTINE FLAG_TEB_GREENROOF_n (TGR, TGRO, TGRPE, T, TVG, &
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
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TEB_GREENROOF_PGD_EVOL_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
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
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_GREENROOF_PGD_EVOL_t), INTENT(INOUT) :: TGRPE
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
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
  DO JL1=1,TGRO%NLAYER_GR
    WHERE (T%XGREENROOF(:)==0.) 
      TGR%XTG (:,JL1) = ZTG
      TGR%XWG (:,JL1) = ZWG
      TGR%XWGI(:,JL1) = ZDEF
    END WHERE
  END DO
  !
  WHERE (T%XGREENROOF(:)==0.) 
    TGR%XWR  (:) = ZWR
    TGR%XRESA(:) = ZRESA
  END WHERE
  !
  IF (TVG%CPHOTO/='NON') THEN
    !
    WHERE (T%XGREENROOF(:)==0.)
      TGR%XANFM (:) = ZANFM              
      TGR%XAN   (:) = ZDEF
      TGR%XANDAY(:) = ZDEF
      TGR%XLE   (:) = ZDEF
    END WHERE
    !
  ELSE IF (TVG%CPHOTO=='LAI' .OR. TVG%CPHOTO=='LST' .OR. TVG%CPHOTO=='NIT' .OR. TVG%CPHOTO=='NCB') THEN
    !
    WHERE (T%XGREENROOF(:)==0.) TGRPE%XLAI(:) = ZDEF
    !
  ELSE IF (TVG%CPHOTO=='AGS' .OR. TVG%CPHOTO=='AST') THEN
    !
    DO JL1=1,SIZE(TGR%XBIOMASS,2)
      WHERE (T%XGREENROOF(:)==0.)
        TGR%XBIOMASS     (:,JL1) = ZDEF
        TGR%XRESP_BIOMASS(:,JL1) = ZDEF
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
 CALL FLAG_GR_SNOW(KFLAG,T%XGREENROOF(:)==0.,TGR%TSNOW)
!
!
!* snow-free characteristics
!
IF (KFLAG==1) THEN
  WHERE (T%XGREENROOF==0.) TGR%XSNOWFREE_ALB      = 0.2
  WHERE (T%XGREENROOF==0.) TGR%XSNOWFREE_ALB_VEG  = 0.2
  WHERE (T%XGREENROOF==0.) TGR%XSNOWFREE_ALB_SOIL = 0.2
ELSEIF (KFLAG==2) THEN
  WHERE (T%XGREENROOF==0.) TGR%XSNOWFREE_ALB      = XUNDEF
  WHERE (T%XGREENROOF==0.) TGR%XSNOWFREE_ALB_VEG  = XUNDEF
  WHERE (T%XGREENROOF==0.) TGR%XSNOWFREE_ALB_SOIL = XUNDEF
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FLAG_TEB_GREENROOF_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FLAG_TEB_GREENROOF_n
