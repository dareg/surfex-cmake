!     #########
SUBROUTINE PREP_VER_TEB_GARDEN (GDR, TGDO, PZS, PDG)
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
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
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
TYPE(ISBA_PROG_t), INTENT(INOUT) :: GDR
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: TGDO
REAL, DIMENSION(:), INTENT(IN) :: PZS
!
REAL, DIMENSION(:,:), INTENT(IN) :: PDG
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
ALLOCATE(ZWGI_CLIM_GRAD (SIZE(GDR%XWG,1),SIZE(GDR%XWG,2)))
!
ZWGI_CLIM_GRAD(:,:) = ZGRADX * EXP( - PDG(:,:) / ZH0 )
!-------------------------------------------------------------------------------------
!
!*      1.1    Temperature profile
!
ALLOCATE(ZTG_LS(SIZE(GDR%XTG,1),SIZE(GDR%XTG,2)))
ZTG_LS(:,:) = GDR%XTG(:,:,1)
!
DO JL=1,SIZE(GDR%XTG,2)
  WHERE(GDR%XTG(:,JL,1)/=XUNDEF) &
    GDR%XTG(:,JL,1) = GDR%XTG(:,JL,1) + XT_CLIM_GRAD  * (PZS - XZS_LS)  
END DO
!
!-------------------------------------------------------------------------------------
!
!*      1.2    Water and ice in the soil
!
ALLOCATE(ZZSFREEZE      (SIZE(GDR%XWG,1)))
ALLOCATE(ZWGTOT         (SIZE(GDR%XWG,1)))
ALLOCATE(ZDW            (SIZE(GDR%XWG,1)))
!
!* general case
!
IWORK=SIZE(GDR%XTG,2)
!
DO JL=1,IWORK
  !
  ZDW(:) = 0.
  ! altitude where deep soil freezes (diurnal surface response is not treated)
  ZZSFREEZE(:) = PZS + (XTT - GDR%XTG(:,JL,1)) / XT_CLIM_GRAD
  !
  WHERE(GDR%XTG(:,JL,1)/=XUNDEF) 
    !
    WHERE (ZTG_LS(:,JL) < XTT)
      !
      WHERE (PZS <= XZS_LS)
        !
        WHERE (PZS > ZZSFREEZE) 
          ZDW(:) = ZWGI_CLIM_GRAD(:,JL) * (PZS - XZS_LS)
        ELSEWHERE
          ZDW(:) = ZWGI_CLIM_GRAD(:,JL) * (ZZSFREEZE - XZS_LS) + ZGRADX * (PZS - ZZSFREEZE)
        ENDWHERE
        !
      ELSEWHERE
        !
        ZDW(:) = ZWGI_CLIM_GRAD(:,JL) * (PZS - XZS_LS)
        !
      ENDWHERE
      !
    ELSEWHERE
      !
      WHERE (PZS <= XZS_LS)
        !
        ZDW(:) = ZGRADX * (PZS - XZS_LS)
        !
      ELSEWHERE
        !
        ZDW(:) = ZWGI_CLIM_GRAD(:,JL) * (PZS - ZZSFREEZE)
        !
      END WHERE
      !
    END WHERE
    !
    ZWGTOT(:) = XUNDEF
    !
    WHERE(GDR%XWG(:,JL,1)/=XUNDEF)         
      ZWGTOT(:) = GDR%XWG(:,JL,1) + GDR%XWGI(:,JL,1)
    ENDWHERE        
    !
    WHERE(GDR%XWG(:,JL,1)/=XUNDEF)      
      GDR%XWGI(:,JL,1) = GDR%XWGI(:,JL,1) + ZDW(:)
      GDR%XWG (:,JL,1) = GDR%XWG (:,JL,1) - ZDW(:)
    ENDWHERE
    !
    WHERE (GDR%XWGI(:,JL,1) < 0..AND.GDR%XWGI(:,JL,1)/=XUNDEF) 
      GDR%XWGI(:,JL,1) = 0.
      GDR%XWG (:,JL,1) = ZWGTOT(:)
    END WHERE
    !
    WHERE (GDR%XWG(:,JL,1) < XWGMIN.AND.GDR%XWG(:,JL,1)/=XUNDEF)
      GDR%XWG (:,JL,1) = XWGMIN
      GDR%XWGI(:,JL,1) = ZWGTOT(:) - XWGMIN
    END WHERE
    !
    WHERE(GDR%XWGI(:,JL,1) > 0..AND.GDR%XWGI(:,JL,1)/=XUNDEF)
      GDR%XTG(:,JL,1) = MIN(XTT,GDR%XTG(:,JL,1))
    ELSEWHERE
      GDR%XTG(:,JL,1) = MAX(XTT,GDR%XTG(:,JL,1))
    ENDWHERE
    !
  ENDWHERE
  !
END DO
!
!* limits in force-restore case
!
IF (TGDO%CISBA=='3-L') THEN 
  WHERE (GDR%XWGI(:,3,1) /= XUNDEF)
    GDR%XWG (:,3,1) = GDR%XWG(:,3,1)+GDR%XWGI(:,3,1)
    GDR%XWGI(:,3,1) = 0.
    GDR%XTG (:,3,1) = ZTG_LS(:,3)
  END WHERE
END IF
!
DEALLOCATE(ZZSFREEZE)
DEALLOCATE(ZWGI_CLIM_GRAD)
DEALLOCATE(ZWGTOT   )
DEALLOCATE(ZDW      )
!
!* masks where fields are not defined
WHERE (GDR%XTG(:,1:SIZE(GDR%XWG,2),1) == XUNDEF)
  GDR%XWG (:,:,1) = XUNDEF
  GDR%XWGI(:,:,1) = XUNDEF
END WHERE
!
!-------------------------------------------------------------------------------------
!
!*      1.4    Snow variables
!
!* vertical shift
IF (TGDO%CISBA=='DIF') THEN
  IDEEP_SOIL = TGDO%NGROUND_LAYER
ELSE
  IDEEP_SOIL = 2
END IF
 CALL PREP_VER_SNOW(GDR%TSNOW,XZS_LS,PZS,SPREAD(ZTG_LS,3,1),GDR%XTG,IDEEP_SOIL)
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
