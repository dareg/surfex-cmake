!     #########
      SUBROUTINE FLAG_TEB_GARDEN_n (TGDR, TGDO, TGDPE, T,  &
                                    KFLAG)
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
!
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_PROG_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_TEB_VEG_PARAM_n, ONLY : TEB_VEG_PARAM_TIME_t
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
TYPE(TEB_VEG_PROG_t), INTENT(INOUT) :: TGDR
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_VEG_PARAM_TIME_t), INTENT(INOUT) :: TGDPE
TYPE(TEB_t), INTENT(INOUT) :: T
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
    WHERE (T%CUR%XGARDEN(:)==0.) 
      TGDR%CUR%XTG (:,JL1,1) = ZTG
      TGDR%CUR%XWG (:,JL1,1) = ZWG
      TGDR%CUR%XWGI(:,JL1,1) = ZDEF
    END WHERE
  END DO
  !
  WHERE (T%CUR%XGARDEN(:)==0.) 
    TGDR%CUR%XWR  (:,1) = ZWR
    TGDR%CUR%XRESA(:,1) = ZRESA
  END WHERE
  !
  IF (TGDO%CPHOTO/='NON') THEN
    !
    WHERE (T%CUR%XGARDEN(:)==0.)
      TGDR%CUR%XANFM (:,1) = ZANFM              
      TGDR%CUR%XAN   (:,1) = ZDEF
      TGDR%CUR%XANDAY(:,1) = ZDEF
      TGDR%CUR%XLE   (:,1) = ZDEF
    END WHERE
    !
    IF (TGDO%CPHOTO=='LAI' .OR. TGDO%CPHOTO=='LST' .OR. TGDO%CPHOTO=='NIT' .OR. TGDO%CPHOTO=='NCB') THEN
      !
      WHERE (T%CUR%XGARDEN(:)==0.) TGDPE%CUR%XLAI(:,1) = ZDEF
      !
    ELSE IF (TGDO%CPHOTO=='AGS' .OR. TGDO%CPHOTO=='AST') THEN
      !
      DO JL1=1,SIZE(TGDR%CUR%XBIOMASS,2)
        WHERE (T%CUR%XGARDEN(:)==0.)
          TGDR%CUR%XBIOMASS     (:,JL1,1) = ZDEF
          TGDR%CUR%XRESP_BIOMASS(:,JL1,1) = ZDEF
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
 CALL FLAG_GR_SNOW(KFLAG,T%CUR%XGARDEN(:)==0.,TGDR%CUR%TSNOW)
!
!
!* snow-free characteristics
!
IF (KFLAG==1) THEN
  WHERE (T%CUR%XGARDEN==0.) TGDR%CUR%XSNOWFREE_ALB(:,1)      = 0.2
  WHERE (T%CUR%XGARDEN==0.) TGDR%CUR%XSNOWFREE_ALB_VEG(:,1)  = 0.2
  WHERE (T%CUR%XGARDEN==0.) TGDR%CUR%XSNOWFREE_ALB_SOIL(:,1) = 0.2
ELSEIF (KFLAG==2) THEN
  WHERE (T%CUR%XGARDEN==0.) TGDR%CUR%XSNOWFREE_ALB(:,1)      = XUNDEF
  WHERE (T%CUR%XGARDEN==0.) TGDR%CUR%XSNOWFREE_ALB_VEG(:,1)  = XUNDEF
  WHERE (T%CUR%XGARDEN==0.) TGDR%CUR%XSNOWFREE_ALB_SOIL(:,1) = XUNDEF
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FLAG_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FLAG_TEB_GARDEN_n
