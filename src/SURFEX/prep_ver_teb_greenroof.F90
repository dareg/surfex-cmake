!     #########
SUBROUTINE PREP_VER_TEB_GREENROOF (GRR, GRRO, PZS, PDG)
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
!
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
!
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
!
TYPE(ISBA_PROG_t), INTENT(INOUT) :: GRR
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GRRO
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
IF (LHOOK) CALL DR_HOOK('PREP_VER_TEB_GREENROOF',0,ZHOOK_HANDLE)
ALLOCATE(ZWGI_CLIM_GRAD (SIZE(GRR%XWG,1),SIZE(GRR%XWG,2)))
!
ZWGI_CLIM_GRAD(:,:) = ZGRADX * EXP( - PDG(:,:) / ZH0 )
!-------------------------------------------------------------------------------------
!
!*      1.1    Temperature profile
!
ALLOCATE(ZTG_LS(SIZE(GRR%XTG,1),SIZE(GRR%XTG,2)))
ZTG_LS(:,:) = GRR%XTG(:,:,1)
!
  DO JL=1,SIZE(GRR%XTG,2)
    WHERE(GRR%XTG(:,JL,1)/=XUNDEF) &
      GRR%XTG(:,JL,1) = GRR%XTG(:,JL,1) + XT_CLIM_GRAD  * (PZS - XZS_LS)  
  END DO
!
!-------------------------------------------------------------------------------------
!
!*      1.2    Water and ice in the soil
!
ALLOCATE(ZZSFREEZE      (SIZE(GRR%XWG,1)))
ALLOCATE(ZWGTOT         (SIZE(GRR%XWG,1)))
ALLOCATE(ZDW            (SIZE(GRR%XWG,1)))
!
!* general case
!
IWORK=SIZE(GRR%XTG,2)
!
  !
  DO JL=1,IWORK
    !
    ZDW(:) = 0.
    ! altitude where deep soil freezes (diurnal surface response is not treated)
    ZZSFREEZE(:) = PZS + (XTT - GRR%XTG(:,JL,1)) / XT_CLIM_GRAD
    !
    WHERE(GRR%XTG(:,JL,1)/=XUNDEF) 
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
      WHERE(GRR%XWG(:,JL,1)/=XUNDEF)         
        ZWGTOT(:) = GRR%XWG(:,JL,1) + GRR%XWGI(:,JL,1)
      ENDWHERE 
      !
      WHERE(GRR%XWG(:,JL,1)/=XUNDEF)      
        GRR%XWGI(:,JL,1) = GRR%XWGI(:,JL,1) + ZDW(:)
        GRR%XWG (:,JL,1) = GRR%XWG (:,JL,1) - ZDW(:)
      ENDWHERE      
      !
      WHERE (GRR%XWGI(:,JL,1) < 0..AND.GRR%XWGI(:,JL,1)/=XUNDEF) 
        GRR%XWGI(:,JL,1) = 0.
        GRR%XWG (:,JL,1) = ZWGTOT(:)
      END WHERE
      !
      WHERE (GRR%XWG(:,JL,1) < XWGMIN.AND.GRR%XWG(:,JL,1)/=XUNDEF)
        GRR%XWG (:,JL,1) = XWGMIN
        GRR%XWGI(:,JL,1) = ZWGTOT(:) - XWGMIN
      END WHERE
      !
      WHERE(GRR%XWGI(:,JL,1) > 0..AND.GRR%XWGI(:,JL,1)/=XUNDEF)
        GRR%XTG(:,JL,1) = MIN(XTT,GRR%XTG(:,JL,1))
      ELSEWHERE
        GRR%XTG(:,JL,1) = MAX(XTT,GRR%XTG(:,JL,1))
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
WHERE (GRR%XTG(:,1:SIZE(GRR%XWG,2),1) == XUNDEF)
  GRR%XWG (:,:,1) = XUNDEF
  GRR%XWGI(:,:,1) = XUNDEF
END WHERE
!
!-------------------------------------------------------------------------------------
!
IDEEP_SOIL = GRRO%NGROUND_LAYER
 CALL PREP_VER_SNOW(GRR%TSNOW,XZS_LS,PZS,SPREAD(ZTG_LS,3,1),GRR%XTG,IDEEP_SOIL)
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
