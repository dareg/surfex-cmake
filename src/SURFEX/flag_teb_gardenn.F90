!     #########
      SUBROUTINE FLAG_TEB_GARDEN_n (VR, VO, PLAI, PGARDEN, KFLAG)
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
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
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
TYPE(ISBA_PROG_t), INTENT(INOUT) :: VR
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: VO
REAL, DIMENSION(:,:), INTENT(INOUT) :: PLAI
REAL, DIMENSION(:), INTENT(IN) :: PGARDEN
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
  DO JL1=1,VO%NGROUND_LAYER
    WHERE (PGARDEN(:)==0.) 
      VR%XTG (:,JL1,1) = ZTG
      VR%XWG (:,JL1,1) = ZWG
      VR%XWGI(:,JL1,1) = ZDEF
    END WHERE
  END DO
  !
  WHERE (PGARDEN(:)==0.) 
    VR%XWR  (:,1) = ZWR
    VR%XRESA(:,1) = ZRESA
  END WHERE
  !
  IF (VO%CPHOTO/='NON') THEN
    !
    WHERE (PGARDEN(:)==0.)
      VR%XANFM (:,1) = ZANFM              
      VR%XAN   (:,1) = ZDEF
      VR%XANDAY(:,1) = ZDEF
      VR%XLE   (:,1) = ZDEF
    END WHERE
    !
    IF (VO%CPHOTO=='LAI' .OR. VO%CPHOTO=='LST' .OR. VO%CPHOTO=='NIT' .OR. VO%CPHOTO=='NCB') THEN
      !
      WHERE (PGARDEN(:)==0.) PLAI(:,1) = ZDEF
      !
    ELSE IF (VO%CPHOTO=='AGS' .OR. VO%CPHOTO=='AST') THEN
      !
      DO JL1=1,SIZE(VR%XBIOMASS,2)
        WHERE (PGARDEN(:)==0.)
          VR%XBIOMASS     (:,JL1,1) = ZDEF
          VR%XRESP_BIOMASS(:,JL1,1) = ZDEF
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
 CALL FLAG_GR_SNOW(KFLAG,PGARDEN(:)==0.,VR%TSNOW)
!
!
!* snow-free characteristics
!
IF (KFLAG==1) THEN
  WHERE (PGARDEN==0.) VR%XSNOWFREE_ALB(:,1)      = 0.2
  WHERE (PGARDEN==0.) VR%XSNOWFREE_ALB_VEG(:,1)  = 0.2
  WHERE (PGARDEN==0.) VR%XSNOWFREE_ALB_SOIL(:,1) = 0.2
ELSEIF (KFLAG==2) THEN
  WHERE (PGARDEN==0.) VR%XSNOWFREE_ALB(:,1)      = XUNDEF
  WHERE (PGARDEN==0.) VR%XSNOWFREE_ALB_VEG(:,1)  = XUNDEF
  WHERE (PGARDEN==0.) VR%XSNOWFREE_ALB_SOIL(:,1) = XUNDEF
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FLAG_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FLAG_TEB_GARDEN_n
