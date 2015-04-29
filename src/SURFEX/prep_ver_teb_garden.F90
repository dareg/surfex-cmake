!     #########
SUBROUTINE PREP_VER_TEB_GARDEN (TGD, TGDO, TGDP, TOP, TVG)
!     #################################################################################
!
!!****  *PREP_VER_TEB_GARDEN* - change in ISBA fields due to altitude change
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified by B. Decharme  (01/2009), Optional Arpege deep soil temperature initialization
!!------------------------------------------------------------------
!

!
!
USE MODD_TEB_GARDEN_n, ONLY : TEB_GARDEN_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
!
USE MODD_ISBA_PAR,       ONLY : XWGMIN
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_PREP,           ONLY : XZS_LS, XT_CLIM_GRAD
USE MODD_CSTS,           ONLY : XTT, XDAY, XLMTT, XRHOLW
!
USE MODE_THERMOS
USE MODI_PREP_VER_SNOW
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!*      0.2    declarations of local variables
!
!
TYPE(TEB_GARDEN_t), INTENT(INOUT) :: TGD
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
!
INTEGER                         :: JL        ! loop counter on layers
INTEGER                         :: IWORK     ! Work integer
!
REAL, DIMENSION(:), ALLOCATABLE :: ZWGTOT    ! total water content
REAL, DIMENSION(:), ALLOCATABLE :: ZDW       ! variation of water in soil
REAL, DIMENSION(:), ALLOCATABLE :: ZZSFREEZE ! altitude where soil temperature equals XTT
INTEGER                         :: IDEEP_SOIL! layer corresponding to deep soil temperature
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWGI_CLIM_GRAD ! ice content vertical gradient
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTG_LS! temperature on initial orography
!
REAL                            :: ZGRADX = 5.E-4 ! slope of ice content gradient
REAL                            :: ZH0    = 5.E-1 ! constant used to define ice content gradient
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!*      1.0    Ice content climatologic gradient
!
IF (LHOOK) CALL DR_HOOK('PREP_VER_TEB_GARDEN',0,ZHOOK_HANDLE)
ALLOCATE(ZWGI_CLIM_GRAD (SIZE(TGD%XWG,1),SIZE(TGD%XWG,2)))
!
ZWGI_CLIM_GRAD(:,:) = ZGRADX * EXP( - TGDP%XDG(:,:) / ZH0 )
!-------------------------------------------------------------------------------------
!
!*      1.1    Temperature profile
!
ALLOCATE(ZTG_LS(SIZE(TGD%XTG,1),SIZE(TGD%XTG,2)))
ZTG_LS(:,:) = TGD%XTG(:,:)
!
DO JL=1,SIZE(TGD%XTG,2)
  WHERE(TGD%XTG(:,JL)/=XUNDEF) &
    TGD%XTG(:,JL) = TGD%XTG(:,JL) + XT_CLIM_GRAD  * (TOP%XZS - XZS_LS)  
END DO
!
!-------------------------------------------------------------------------------------
!
!*      1.2    Water and ice in the soil
!
ALLOCATE(ZZSFREEZE      (SIZE(TGD%XWG,1)))
ALLOCATE(ZWGTOT         (SIZE(TGD%XWG,1)))
ALLOCATE(ZDW            (SIZE(TGD%XWG,1)))
!
!* general case
!
IWORK=SIZE(TGD%XTG,2)
!
DO JL=1,IWORK
  !
  ZDW(:) = 0.
  ! altitude where deep soil freezes (diurnal surface response is not treated)
  ZZSFREEZE(:) = TOP%XZS + (XTT - TGD%XTG(:,JL)) / XT_CLIM_GRAD
  !
  WHERE(TGD%XTG(:,JL)/=XUNDEF) 
    !
    WHERE (ZTG_LS(:,JL) < XTT)
      !
      WHERE (TOP%XZS <= XZS_LS)
        !
        WHERE (TOP%XZS > ZZSFREEZE) 
          ZDW(:) = ZWGI_CLIM_GRAD(:,JL) * (TOP%XZS - XZS_LS)
        ELSEWHERE
          ZDW(:) = ZWGI_CLIM_GRAD(:,JL) * (ZZSFREEZE - XZS_LS) + ZGRADX * (TOP%XZS - ZZSFREEZE)
        ENDWHERE
        !
      ELSEWHERE
        !
        ZDW(:) = ZWGI_CLIM_GRAD(:,JL) * (TOP%XZS - XZS_LS)
        !
      ENDWHERE
      !
    ELSEWHERE
      !
      WHERE (TOP%XZS <= XZS_LS)
        !
        ZDW(:) = ZGRADX * (TOP%XZS - XZS_LS)
        !
      ELSEWHERE
        !
        ZDW(:) = ZWGI_CLIM_GRAD(:,JL) * (TOP%XZS - ZZSFREEZE)
        !
      END WHERE
      !
    END WHERE
    !
    ZWGTOT(:) = XUNDEF
    !
    WHERE(TGD%XWG(:,JL)/=XUNDEF)         
      ZWGTOT(:) = TGD%XWG(:,JL) + TGD%XWGI(:,JL)
    ENDWHERE        
    !
    WHERE(TGD%XWG(:,JL)/=XUNDEF)      
      TGD%XWGI(:,JL) = TGD%XWGI(:,JL) + ZDW(:)
      TGD%XWG (:,JL) = TGD%XWG (:,JL) - ZDW(:)
    ENDWHERE
    !
    WHERE (TGD%XWGI(:,JL) < 0..AND.TGD%XWGI(:,JL)/=XUNDEF) 
      TGD%XWGI(:,JL) = 0.
      TGD%XWG (:,JL) = ZWGTOT(:)
    END WHERE
    !
    WHERE (TGD%XWG(:,JL) < XWGMIN.AND.TGD%XWG(:,JL)/=XUNDEF)
      TGD%XWG (:,JL) = XWGMIN
      TGD%XWGI(:,JL) = ZWGTOT(:) - XWGMIN
    END WHERE
    !
    WHERE(TGD%XWGI(:,JL) > 0..AND.TGD%XWGI(:,JL)/=XUNDEF)
      TGD%XTG(:,JL) = MIN(XTT,TGD%XTG(:,JL))
    ELSEWHERE
      TGD%XTG(:,JL) = MAX(XTT,TGD%XTG(:,JL))
    ENDWHERE
    !
  ENDWHERE
  !
END DO
!
!* limits in force-restore case
!
IF (TVG%CISBA=='3-L') THEN 
  WHERE (TGD%XWGI(:,3) /= XUNDEF)
    TGD%XWG (:,3) = TGD%XWG(:,3)+TGD%XWGI(:,3)
    TGD%XWGI(:,3) = 0.
    TGD%XTG (:,3) = ZTG_LS(:,3)
  END WHERE
END IF
!
DEALLOCATE(ZZSFREEZE)
DEALLOCATE(ZWGI_CLIM_GRAD)
DEALLOCATE(ZWGTOT   )
DEALLOCATE(ZDW      )
!
!* masks where fields are not defined
WHERE (TGD%XTG(:,1:SIZE(TGD%XWG,2)) == XUNDEF)
  TGD%XWG (:,:) = XUNDEF
  TGD%XWGI(:,:) = XUNDEF
END WHERE
!
!-------------------------------------------------------------------------------------
!
!*      1.4    Snow variables
!
!* vertical shift
IF (TVG%CISBA=='DIF') THEN
  IDEEP_SOIL = TGDO%NGROUND_LAYER
ELSE
  IDEEP_SOIL = 2
END IF
 CALL PREP_VER_SNOW(TGD%TSNOW,XZS_LS,TOP%XZS,SPREAD(ZTG_LS,3,1),SPREAD(TGD%XTG,3,1),IDEEP_SOIL)
!
!-------------------------------------------------------------------------------------
!
!*      2.     Deallocation of large-scale orography
!
DEALLOCATE(ZTG_LS)
IF (LHOOK) CALL DR_HOOK('PREP_VER_TEB_GARDEN',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_VER_TEB_GARDEN
