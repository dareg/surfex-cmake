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
!!	V. Masson   *Meteo France*	
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
USE MODD_TEB_n,          ONLY : XGARDEN
USE MODD_TEB_GARDEN_n,   ONLY : NGROUND_LAYER, NPATCH,              &
                                CPHOTO, CISBA, CRESPSL,             &
                                XTG, XWG, XWGI, XWR, XLAI, TSNOW,   &
                                XRESA, XANFM, XAN, XLE, XANDAY,     &
                                XBSLAI, XBIOMASS, XRESP_BIOMASS,    &
                                XLITTER, XSOILCARB, XLIGNIN_STRUC  
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
INTEGER :: JL1, JL2, JPATCH ! loop counter on layers
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
DO JPATCH=1,NPATCH
  !
  DO JL1=1,NGROUND_LAYER
    WHERE (XGARDEN(:)==0.) 
      XTG (:,JL1,JPATCH) = ZTG
      XWG (:,JL1,JPATCH) = ZWG
      XWGI(:,JL1,JPATCH) = ZDEF
    END WHERE
  END DO
  !
  WHERE (XGARDEN(:)==0.) 
    XWR  (:,JPATCH) = ZWR
    XRESA(:,JPATCH) = ZRESA
  END WHERE
  !
  IF (CPHOTO/='NON') THEN
    !
    WHERE (XGARDEN(:)==0.)
      XANFM (:,JPATCH) = ZANFM              
      XAN   (:,JPATCH) = ZDEF
      XANDAY(:,JPATCH) = ZDEF
      XLE   (:,JPATCH) = ZDEF
    END WHERE
    !
    IF (CPHOTO=='LAI' .OR. CPHOTO=='LST' .OR. CPHOTO=='NIT' .OR. CPHOTO=='NCB') THEN
      !
      WHERE (XGARDEN(:)==0.) XLAI(:,JPATCH) = ZDEF
      !
    ELSE IF (CPHOTO=='AGS' .OR. CPHOTO=='AST') THEN
      !
      DO JL1=1,SIZE(XBIOMASS,2)
        WHERE (XGARDEN(:)==0.)
          XBIOMASS     (:,JL1,JPATCH) = ZDEF
          XRESP_BIOMASS(:,JL1,JPATCH) = ZDEF
        END WHERE
      END DO
      !
    END IF
    !
    IF (CRESPSL=='CNT') THEN
      !
      DO JL2=1,SIZE(XLITTER,3)
        DO JL1=1,SIZE(XLITTER,2)
          WHERE(XGARDEN(:)==0.) XLITTER  (:,JL1,JL2,JPATCH) = ZDEF
        END DO
      END DO
      DO JL1=1,SIZE(XSOILCARB,2)
        WHERE(XGARDEN(:)==0.)   XSOILCARB    (:,JL1,JPATCH) = ZDEF
      END DO
      DO JL1=1,SIZE(XLIGNIN_STRUC,2)
        WHERE(XGARDEN(:)==0.)   XLIGNIN_STRUC(:,JL1,JPATCH) = ZDEF
      END DO
      !
    END IF
    !
  ENDIF
  !
END DO
!
!-------------------------------------------------------------------------------
!
!* Flag snow characteristics
!
CALL FLAG_GR_SNOW(KFLAG,XGARDEN(:)==0.,TSNOW)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FLAG_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FLAG_TEB_GARDEN_n
