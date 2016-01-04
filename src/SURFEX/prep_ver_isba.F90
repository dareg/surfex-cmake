!     #########
SUBROUTINE PREP_VER_ISBA (I)
!     #################################################################################
!
!!****  *PREP_VER_ISBA* - change in ISBA fields due to altitude change
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
!!      S. Riette   04/2010 Modification of XTG corrections after freezing
!!------------------------------------------------------------------
!

!
!
USE MODD_ISBA_n, ONLY : ISBA_t
!
USE MODD_ISBA_PAR,       ONLY : XWGMIN
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_PREP,           ONLY : XZS_LS, XT_CLIM_GRAD
USE MODD_PREP_ISBA,      ONLY : LSNOW_IDEAL
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
TYPE(ISBA_t), INTENT(INOUT) :: I
!
INTEGER                         :: JL        ! loop counter on layers
INTEGER                         :: JP        ! loop counter on patches
INTEGER                         :: IWORK     ! Work integer
!
REAL, DIMENSION(:), ALLOCATABLE :: ZWGTOT    ! total water content
REAL, DIMENSION(:), ALLOCATABLE :: ZDW       ! variation of water in soil
REAL, DIMENSION(:), ALLOCATABLE :: ZZSFREEZE ! altitude where soil temperature equals XTT
INTEGER                         :: IDEEP_SOIL! layer corresponding to deep soil temperature
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWGI_CLIM_GRAD ! ice content vertical gradient
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTG_LS! temperature on initial orography
!
REAL                            :: ZGRADX = 5.E-4 ! slope of ice content gradient
REAL                            :: ZH0    = 5.E-1 ! constant used to define ice content gradient
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!*      1.0    Ice content climatologic gradient
!
IF (LHOOK) CALL DR_HOOK('PREP_VER_ISBA',0,ZHOOK_HANDLE)
ALLOCATE(ZWGI_CLIM_GRAD (SIZE(I%R%XWG,1),SIZE(I%R%XWG,2),SIZE(I%R%XWG,3)))
!
ZWGI_CLIM_GRAD(:,:,:) = ZGRADX * EXP( - I%M%X%XDG(:,:,:) / ZH0 )
!-------------------------------------------------------------------------------------
!
!*      1.1    Temperature profile
!
ALLOCATE(ZTG_LS(SIZE(I%R%XTG,1),SIZE(I%R%XTG,2),SIZE(I%R%XTG,3)))
ZTG_LS(:,:,:) = I%R%XTG(:,:,:)
!
DO JP=1,SIZE(I%R%XTG,3)
  DO JL=1,SIZE(I%R%XTG,2)
    WHERE(I%R%XTG(:,JL,JP)/=XUNDEF) &
      I%R%XTG(:,JL,JP) = I%R%XTG(:,JL,JP) + XT_CLIM_GRAD  * (I%P%XZS - XZS_LS)  
  END DO
END DO
!
!-------------------------------------------------------------------------------------
!
!*      1.2    Water and ice in the soil
!
ALLOCATE(ZZSFREEZE      (SIZE(I%R%XWG,1)))
ALLOCATE(ZWGTOT         (SIZE(I%R%XWG,1)))
ALLOCATE(ZDW            (SIZE(I%R%XWG,1)))
!
!* general case
!
IF(I%O%LTEMP_ARP)THEN
  IWORK=SIZE(I%R%XWG,2)
ELSE
  IWORK=SIZE(I%R%XTG,2)
ENDIF
!
DO JP=1,SIZE(I%R%XWG,3)
  !
  DO JL=1,IWORK
    !
    ZDW(:) = 0.
    ! altitude where deep soil freezes (diurnal surface response is not treated)
    ZZSFREEZE(:) = I%P%XZS + (XTT - I%R%XTG(:,JL,JP)) / XT_CLIM_GRAD
    !
    WHERE(I%R%XTG(:,JL,JP)/=XUNDEF) 
      !
      WHERE (ZTG_LS(:,JL,JP) < XTT)
        !
        WHERE (I%P%XZS <= XZS_LS)
          !
          WHERE (I%P%XZS > ZZSFREEZE) 
            ZDW(:) = ZWGI_CLIM_GRAD(:,JL,JP) * (I%P%XZS - XZS_LS)
          ELSEWHERE
            ZDW(:) = ZWGI_CLIM_GRAD(:,JL,JP) * (ZZSFREEZE - XZS_LS) + ZGRADX * (I%P%XZS - ZZSFREEZE)
          ENDWHERE
          !
        ELSEWHERE
          !
          ZDW(:) = ZWGI_CLIM_GRAD(:,JL,JP) * (I%P%XZS - XZS_LS)
          !
        ENDWHERE
        !
      ELSEWHERE
        !
        WHERE (I%P%XZS <= XZS_LS)
          !
          ZDW(:) = ZGRADX * (I%P%XZS - XZS_LS)
          !
        ELSEWHERE
          !
          ZDW(:) = ZWGI_CLIM_GRAD(:,JL,JP) * (I%P%XZS - ZZSFREEZE)
          !
        END WHERE
        !
      END WHERE 
      !
      ZWGTOT(:) = XUNDEF
      !
      WHERE(I%R%XWG(:,JL,JP)/=XUNDEF)        
        ZWGTOT(:) = I%R%XWG(:,JL,JP) + I%R%XWGI(:,JL,JP)
      ENDWHERE
      !
      WHERE(I%R%XWG(:,JL,JP)/=XUNDEF)        
        I%R%XWGI(:,JL,JP) = I%R%XWGI(:,JL,JP) + ZDW(:)
        I%R%XWG (:,JL,JP) = I%R%XWG (:,JL,JP) - ZDW(:)
      ENDWHERE
      !
      WHERE (I%R%XWGI(:,JL,JP)<0.0.AND.I%R%XWGI(:,JL,JP)/=XUNDEF) 
        I%R%XWGI(:,JL,JP) = 0.
        I%R%XWG (:,JL,JP) = ZWGTOT(:)
      END WHERE
      !
      WHERE (I%R%XWG(:,JL,JP)<XWGMIN.AND.I%R%XWG(:,JL,JP)/=XUNDEF)
        I%R%XWG (:,JL,JP) = XWGMIN
        I%R%XWGI(:,JL,JP) = ZWGTOT(:) - XWGMIN
      END WHERE
      !
      WHERE(I%R%XWGI(:,JL,JP)>0.0.AND.I%R%XWGI(:,JL,JP)/=XUNDEF)
        I%R%XTG(:,JL,JP) = MIN(XTT,I%R%XTG(:,JL,JP))
      ELSEWHERE
        I%R%XTG(:,JL,JP) = MAX(XTT,I%R%XTG(:,JL,JP))
      ENDWHERE
      !
    END WHERE
    !
  END DO
  !
END DO
!
!* limits in force-restore case
!
IF (I%O%CISBA=='3-L') THEN 
  DO JP=1,SIZE(I%R%XWG,3)
     WHERE (I%R%XWGI(:,3,JP) /= XUNDEF)
       I%R%XWG (:,3,JP) = I%R%XWG(:,3,JP)+I%R%XWGI(:,3,JP)
       I%R%XWGI(:,3,JP) = 0.
       I%R%XTG (:,3,JP) = ZTG_LS(:,3,JP) + XT_CLIM_GRAD  * (I%P%XZS - XZS_LS)       
     END WHERE
     IF(I%O%LTEMP_ARP)THEN
        I%R%XTG (:,4:SIZE(I%R%XTG,2),JP) = ZTG_LS(:,4:SIZE(I%R%XTG,2),JP)
     ENDIF
  END DO
ELSEIF(I%O%CISBA=='2-L'.AND.I%O%LTEMP_ARP) THEN
  DO JP=1,SIZE(I%R%XWG,3)
     I%R%XTG (:,3:SIZE(I%R%XTG,2),JP) = ZTG_LS(:,3:SIZE(I%R%XTG,2),JP)
  END DO
END IF
!
DEALLOCATE(ZZSFREEZE)
DEALLOCATE(ZWGI_CLIM_GRAD)
DEALLOCATE(ZWGTOT   )
DEALLOCATE(ZDW      )
!
!* masks where fields are not defined
WHERE (I%R%XTG(:,1:SIZE(I%R%XWG,2),:) == XUNDEF)
  I%R%XWG (:,:,:) = XUNDEF
  I%R%XWGI(:,:,:) = XUNDEF
END WHERE
!
!-------------------------------------------------------------------------------------
!
!*      1.4    Snow variables
!
!* vertical shift
IF (.NOT.LSNOW_IDEAL) THEN
  IF (I%O%CISBA=='DIF') THEN
    IDEEP_SOIL = I%O%NGROUND_LAYER
  ELSE
    IDEEP_SOIL = 2
  END IF        
  CALL PREP_VER_SNOW(I%R%TSNOW,XZS_LS,I%P%XZS,ZTG_LS,I%R%XTG,IDEEP_SOIL)
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      2.     Deallocation of large-scale orography
!
DEALLOCATE(ZTG_LS)
IF (LHOOK) CALL DR_HOOK('PREP_VER_ISBA',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_VER_ISBA
