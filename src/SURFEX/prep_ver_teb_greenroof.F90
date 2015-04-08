!     #########
SUBROUTINE PREP_VER_TEB_GREENROOF
!     #################################################################################
!
!!****  *PREP_VER_TEB_GREENROOF* - change in ISBA fields due to altitude change
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
!!     V. Masson  + A.Lemonsu & C.deMunck
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified by B. Decharme  (01/2009), Optional Arpege deep soil temperature initialization
!!------------------------------------------------------------------
!

!
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TGRP => TEB_GREENROOF_PGD
USE MODD_TEB_GREENROOF_n, ONLY : TGR => TEB_GREENROOF
USE MODD_ISBA_PAR,          ONLY : XWGMIN
USE MODD_SURF_PAR,          ONLY : XUNDEF
USE MODD_PREP,              ONLY : XZS_LS, XT_CLIM_GRAD
USE MODD_CSTS,              ONLY : XTT, XDAY, XLMTT, XRHOLW
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
IF (LHOOK) CALL DR_HOOK('PREP_VER_TEB_GREENROOF',0,ZHOOK_HANDLE)
ALLOCATE(ZWGI_CLIM_GRAD (SIZE(TGR%XWG,1),SIZE(TGR%XWG,2)))
!
ZWGI_CLIM_GRAD(:,:) = ZGRADX * EXP( - TGRP%XDG(:,:) / ZH0 )
!-------------------------------------------------------------------------------------
!
!*      1.1    Temperature profile
!
ALLOCATE(ZTG_LS(SIZE(TGR%XTG,1),SIZE(TGR%XTG,2)))
ZTG_LS(:,:) = TGR%XTG(:,:)
!
  DO JL=1,SIZE(TGR%XTG,2)
    WHERE(TGR%XTG(:,JL)/=XUNDEF) &
      TGR%XTG(:,JL) = TGR%XTG(:,JL) + XT_CLIM_GRAD  * (TOP%XZS - XZS_LS)  
  END DO
!
!-------------------------------------------------------------------------------------
!
!*      1.2    Water and ice in the soil
!
ALLOCATE(ZZSFREEZE      (SIZE(TGR%XWG,1)))
ALLOCATE(ZWGTOT         (SIZE(TGR%XWG,1)))
ALLOCATE(ZDW            (SIZE(TGR%XWG,1)))
!
!* general case
!
IWORK=SIZE(TGR%XTG,2)
!
  !
  DO JL=1,IWORK
    !
    ZDW(:) = 0.
    ! altitude where deep soil freezes (diurnal surface response is not treated)
    ZZSFREEZE(:) = TOP%XZS + (XTT - TGR%XTG(:,JL)) / XT_CLIM_GRAD
    !
    WHERE(TGR%XTG(:,JL)/=XUNDEF) 
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
      WHERE(TGR%XWG(:,JL)/=XUNDEF)         
        ZWGTOT(:) = TGR%XWG(:,JL) + TGR%XWGI(:,JL)
      ENDWHERE 
      !
      WHERE(TGR%XWG(:,JL)/=XUNDEF)      
        TGR%XWGI(:,JL) = TGR%XWGI(:,JL) + ZDW(:)
        TGR%XWG (:,JL) = TGR%XWG (:,JL) - ZDW(:)
      ENDWHERE      
      !
      WHERE (TGR%XWGI(:,JL) < 0..AND.TGR%XWGI(:,JL)/=XUNDEF) 
        TGR%XWGI(:,JL) = 0.
        TGR%XWG (:,JL) = ZWGTOT(:)
      END WHERE
      !
      WHERE (TGR%XWG(:,JL) < XWGMIN.AND.TGR%XWG(:,JL)/=XUNDEF)
        TGR%XWG (:,JL) = XWGMIN
        TGR%XWGI(:,JL) = ZWGTOT(:) - XWGMIN
      END WHERE
      !
      WHERE(TGR%XWGI(:,JL) > 0..AND.TGR%XWGI(:,JL)/=XUNDEF)
        TGR%XTG(:,JL) = MIN(XTT,TGR%XTG(:,JL))
      ELSEWHERE
        TGR%XTG(:,JL) = MAX(XTT,TGR%XTG(:,JL))
      ENDWHERE
      !
    ENDWHERE
    !
  END DO
  !
!
!
DEALLOCATE(ZZSFREEZE     )
DEALLOCATE(ZWGI_CLIM_GRAD)
DEALLOCATE(ZWGTOT        )
DEALLOCATE(ZDW           )
!
!* masks where fields are not defined
WHERE (TGR%XTG(:,1:SIZE(TGR%XWG,2)) == XUNDEF)
  TGR%XWG (:,:) = XUNDEF
  TGR%XWGI(:,:) = XUNDEF
END WHERE
!
!-------------------------------------------------------------------------------------
!
IDEEP_SOIL = TGRO%NLAYER_GR
 CALL PREP_VER_SNOW(TGR%TSNOW,XZS_LS,TOP%XZS,SPREAD(ZTG_LS,3,1),SPREAD(TGR%XTG,3,1),IDEEP_SOIL)
!
!-------------------------------------------------------------------------------------
!
!*      2.     Deallocation of large-scale orography
!
DEALLOCATE(ZTG_LS)
IF (LHOOK) CALL DR_HOOK('PREP_VER_TEB_GREENROOF',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_VER_TEB_GREENROOF
