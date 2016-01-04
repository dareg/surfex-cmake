!     #########
SUBROUTINE PREP_VER_TEB_GARDEN (GDR, TGDO, TGDP, TOP, PDG)
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
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_PROG_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
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
TYPE(TEB_VEG_PROG_t), INTENT(INOUT) :: GDR
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(ISBA_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
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
ALLOCATE(ZWGI_CLIM_GRAD (SIZE(GDR%CUR%XWG,1),SIZE(GDR%CUR%XWG,2)))
!
ZWGI_CLIM_GRAD(:,:) = ZGRADX * EXP( - PDG(:,:) / ZH0 )
!-------------------------------------------------------------------------------------
!
!*      1.1    Temperature profile
!
ALLOCATE(ZTG_LS(SIZE(GDR%CUR%XTG,1),SIZE(GDR%CUR%XTG,2)))
ZTG_LS(:,:) = GDR%CUR%XTG(:,:,1)
!
DO JL=1,SIZE(GDR%CUR%XTG,2)
  WHERE(GDR%CUR%XTG(:,JL,1)/=XUNDEF) &
    GDR%CUR%XTG(:,JL,1) = GDR%CUR%XTG(:,JL,1) + XT_CLIM_GRAD  * (TOP%XZS - XZS_LS)  
END DO
!
!-------------------------------------------------------------------------------------
!
!*      1.2    Water and ice in the soil
!
ALLOCATE(ZZSFREEZE      (SIZE(GDR%CUR%XWG,1)))
ALLOCATE(ZWGTOT         (SIZE(GDR%CUR%XWG,1)))
ALLOCATE(ZDW            (SIZE(GDR%CUR%XWG,1)))
!
!* general case
!
IWORK=SIZE(GDR%CUR%XTG,2)
!
DO JL=1,IWORK
  !
  ZDW(:) = 0.
  ! altitude where deep soil freezes (diurnal surface response is not treated)
  ZZSFREEZE(:) = TOP%XZS + (XTT - GDR%CUR%XTG(:,JL,1)) / XT_CLIM_GRAD
  !
  WHERE(GDR%CUR%XTG(:,JL,1)/=XUNDEF) 
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
    WHERE(GDR%CUR%XWG(:,JL,1)/=XUNDEF)         
      ZWGTOT(:) = GDR%CUR%XWG(:,JL,1) + GDR%CUR%XWGI(:,JL,1)
    ENDWHERE        
    !
    WHERE(GDR%CUR%XWG(:,JL,1)/=XUNDEF)      
      GDR%CUR%XWGI(:,JL,1) = GDR%CUR%XWGI(:,JL,1) + ZDW(:)
      GDR%CUR%XWG (:,JL,1) = GDR%CUR%XWG (:,JL,1) - ZDW(:)
    ENDWHERE
    !
    WHERE (GDR%CUR%XWGI(:,JL,1) < 0..AND.GDR%CUR%XWGI(:,JL,1)/=XUNDEF) 
      GDR%CUR%XWGI(:,JL,1) = 0.
      GDR%CUR%XWG (:,JL,1) = ZWGTOT(:)
    END WHERE
    !
    WHERE (GDR%CUR%XWG(:,JL,1) < XWGMIN.AND.GDR%CUR%XWG(:,JL,1)/=XUNDEF)
      GDR%CUR%XWG (:,JL,1) = XWGMIN
      GDR%CUR%XWGI(:,JL,1) = ZWGTOT(:) - XWGMIN
    END WHERE
    !
    WHERE(GDR%CUR%XWGI(:,JL,1) > 0..AND.GDR%CUR%XWGI(:,JL,1)/=XUNDEF)
      GDR%CUR%XTG(:,JL,1) = MIN(XTT,GDR%CUR%XTG(:,JL,1))
    ELSEWHERE
      GDR%CUR%XTG(:,JL,1) = MAX(XTT,GDR%CUR%XTG(:,JL,1))
    ENDWHERE
    !
  ENDWHERE
  !
END DO
!
!* limits in force-restore case
!
IF (TGDO%CISBA=='3-L') THEN 
  WHERE (GDR%CUR%XWGI(:,3,1) /= XUNDEF)
    GDR%CUR%XWG (:,3,1) = GDR%CUR%XWG(:,3,1)+GDR%CUR%XWGI(:,3,1)
    GDR%CUR%XWGI(:,3,1) = 0.
    GDR%CUR%XTG (:,3,1) = ZTG_LS(:,3)
  END WHERE
END IF
!
DEALLOCATE(ZZSFREEZE)
DEALLOCATE(ZWGI_CLIM_GRAD)
DEALLOCATE(ZWGTOT   )
DEALLOCATE(ZDW      )
!
!* masks where fields are not defined
WHERE (GDR%CUR%XTG(:,1:SIZE(GDR%CUR%XWG,2),1) == XUNDEF)
  GDR%CUR%XWG (:,:,1) = XUNDEF
  GDR%CUR%XWGI(:,:,1) = XUNDEF
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
 CALL PREP_VER_SNOW(GDR%CUR%TSNOW,XZS_LS,TOP%XZS,SPREAD(ZTG_LS,3,1),GDR%CUR%XTG,IDEEP_SOIL)
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
