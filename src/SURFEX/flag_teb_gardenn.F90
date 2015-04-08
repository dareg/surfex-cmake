!     #########
      SUBROUTINE FLAG_TEB_GARDEN_n(KFLAG)
!     ##################################
!
!!****  *FLAG_TEB_GARDEN_n* - routine to flag ISBA variables where gardens are
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
USE MODD_CO2V_PAR,       ONLY : XANFMINIT, XCONDCTMIN
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
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
INTEGER, INTENT(IN) :: KFLAG ! 1 : to put physical values to run ISBA afterwards
!                            ! 2 : to flag with XUNDEF value for points wihtout garden
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
IF (LHOOK) CALL DR_HOOK('FLAG_TEB_GARDEN_N',0,ZHOOK_HANDLE)
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
  DO JL1=1,TGDO%NGROUND_LAYER
    WHERE (T%XGARDEN(:)==0.) 
      TGD%XTG (:,JL1) = ZTG
      TGD%XWG (:,JL1) = ZWG
      TGD%XWGI(:,JL1) = ZDEF
    END WHERE
  END DO
  !
  WHERE (T%XGARDEN(:)==0.) 
    TGD%XWR  (:) = ZWR
    TGD%XRESA(:) = ZRESA
  END WHERE
  !
  IF (TVG%CPHOTO/='NON') THEN
    !
    WHERE (T%XGARDEN(:)==0.)
      TGD%XANFM (:) = ZANFM              
      TGD%XAN   (:) = ZDEF
      TGD%XANDAY(:) = ZDEF
      TGD%XLE   (:) = ZDEF
    END WHERE
    !
    IF (TVG%CPHOTO=='LAI' .OR. TVG%CPHOTO=='LST' .OR. TVG%CPHOTO=='NIT' .OR. TVG%CPHOTO=='NCB') THEN
      !
      WHERE (T%XGARDEN(:)==0.) TGDPE%XLAI(:) = ZDEF
      !
    ELSE IF (TVG%CPHOTO=='AGS' .OR. TVG%CPHOTO=='AST') THEN
      !
      DO JL1=1,SIZE(TGD%XBIOMASS,2)
        WHERE (T%XGARDEN(:)==0.)
          TGD%XBIOMASS     (:,JL1) = ZDEF
          TGD%XRESP_BIOMASS(:,JL1) = ZDEF
        END WHERE
      END DO
      !
    END IF
    !
  ENDIF
  !
!
!-------------------------------------------------------------------------------
!
!* Flag snow characteristics
!
 CALL FLAG_GR_SNOW(KFLAG,T%XGARDEN(:)==0.,TGD%TSNOW)
!
!
!* snow-free characteristics
!
IF (KFLAG==1) THEN
  WHERE (T%XGARDEN==0.) TGD%XSNOWFREE_ALB      = 0.2
  WHERE (T%XGARDEN==0.) TGD%XSNOWFREE_ALB_VEG  = 0.2
  WHERE (T%XGARDEN==0.) TGD%XSNOWFREE_ALB_SOIL = 0.2
ELSEIF (KFLAG==2) THEN
  WHERE (T%XGARDEN==0.) TGD%XSNOWFREE_ALB      = XUNDEF
  WHERE (T%XGARDEN==0.) TGD%XSNOWFREE_ALB_VEG  = XUNDEF
  WHERE (T%XGARDEN==0.) TGD%XSNOWFREE_ALB_SOIL = XUNDEF
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FLAG_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FLAG_TEB_GARDEN_n
