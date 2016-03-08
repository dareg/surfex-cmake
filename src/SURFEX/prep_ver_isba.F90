!     #########
SUBROUTINE PREP_VER_ISBA (IO, IR, PZS, PDG)
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
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
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
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
REAL, DIMENSION(:), INTENT(IN) :: PZS
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDG
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
ALLOCATE(ZWGI_CLIM_GRAD (SIZE(IR%XWG,1),SIZE(IR%XWG,2),SIZE(IR%XWG,3)))
!
ZWGI_CLIM_GRAD(:,:,:) = ZGRADX * EXP( - PDG(:,:,:) / ZH0 )
!-------------------------------------------------------------------------------------
!
!*      1.1    Temperature profile
!
ALLOCATE(ZTG_LS(SIZE(IR%XTG,1),SIZE(IR%XTG,2),SIZE(IR%XTG,3)))
ZTG_LS(:,:,:) = IR%XTG(:,:,:)
!
DO JP=1,SIZE(IR%XTG,3)
  DO JL=1,SIZE(IR%XTG,2)
    WHERE(IR%XTG(:,JL,JP)/=XUNDEF) &
      IR%XTG(:,JL,JP) = IR%XTG(:,JL,JP) + XT_CLIM_GRAD  * (PZS - XZS_LS)  
  END DO
END DO
!
!-------------------------------------------------------------------------------------
!
!*      1.2    Water and ice in the soil
!
ALLOCATE(ZZSFREEZE      (SIZE(IR%XWG,1)))
ALLOCATE(ZWGTOT         (SIZE(IR%XWG,1)))
ALLOCATE(ZDW            (SIZE(IR%XWG,1)))
!
!* general case
!
IF(IO%LTEMP_ARP)THEN
  IWORK=SIZE(IR%XWG,2)
ELSE
  IWORK=SIZE(IR%XTG,2)
ENDIF
!
DO JP=1,SIZE(IR%XWG,3)
  !
  DO JL=1,IWORK
    !
    ZDW(:) = 0.
    ! altitude where deep soil freezes (diurnal surface response is not treated)
    ZZSFREEZE(:) = PZS + (XTT - IR%XTG(:,JL,JP)) / XT_CLIM_GRAD
    !
    WHERE(IR%XTG(:,JL,JP)/=XUNDEF) 
      !
      WHERE (ZTG_LS(:,JL,JP) < XTT)
        !
        WHERE (PZS <= XZS_LS)
          !
          WHERE (PZS > ZZSFREEZE) 
            ZDW(:) = ZWGI_CLIM_GRAD(:,JL,JP) * (PZS - XZS_LS)
          ELSEWHERE
            ZDW(:) = ZWGI_CLIM_GRAD(:,JL,JP) * (ZZSFREEZE - XZS_LS) + ZGRADX * (PZS - ZZSFREEZE)
          ENDWHERE
          !
        ELSEWHERE
          !
          ZDW(:) = ZWGI_CLIM_GRAD(:,JL,JP) * (PZS - XZS_LS)
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
          ZDW(:) = ZWGI_CLIM_GRAD(:,JL,JP) * (PZS - ZZSFREEZE)
          !
        END WHERE
        !
      END WHERE 
      !
      ZWGTOT(:) = XUNDEF
      !
      WHERE(IR%XWG(:,JL,JP)/=XUNDEF)        
        ZWGTOT(:) = IR%XWG(:,JL,JP) + IR%XWGI(:,JL,JP)
      ENDWHERE
      !
      WHERE(IR%XWG(:,JL,JP)/=XUNDEF)        
        IR%XWGI(:,JL,JP) = IR%XWGI(:,JL,JP) + ZDW(:)
        IR%XWG (:,JL,JP) = IR%XWG (:,JL,JP) - ZDW(:)
      ENDWHERE
      !
      WHERE (IR%XWGI(:,JL,JP)<0.0.AND.IR%XWGI(:,JL,JP)/=XUNDEF) 
        IR%XWGI(:,JL,JP) = 0.
        IR%XWG (:,JL,JP) = ZWGTOT(:)
      END WHERE
      !
      WHERE (IR%XWG(:,JL,JP)<XWGMIN.AND.IR%XWG(:,JL,JP)/=XUNDEF)
        IR%XWG (:,JL,JP) = XWGMIN
        IR%XWGI(:,JL,JP) = ZWGTOT(:) - XWGMIN
      END WHERE
      !
      WHERE(IR%XWGI(:,JL,JP)>0.0.AND.IR%XWGI(:,JL,JP)/=XUNDEF)
        IR%XTG(:,JL,JP) = MIN(XTT,IR%XTG(:,JL,JP))
      ELSEWHERE
        IR%XTG(:,JL,JP) = MAX(XTT,IR%XTG(:,JL,JP))
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
IF (IO%CISBA=='3-L') THEN 
  DO JP=1,SIZE(IR%XWG,3)
     WHERE (IR%XWGI(:,3,JP) /= XUNDEF)
       IR%XWG (:,3,JP) = IR%XWG(:,3,JP)+IR%XWGI(:,3,JP)
       IR%XWGI(:,3,JP) = 0.
       IR%XTG (:,3,JP) = ZTG_LS(:,3,JP) + XT_CLIM_GRAD  * (PZS - XZS_LS)       
     END WHERE
     IF(IO%LTEMP_ARP)THEN
        IR%XTG (:,4:SIZE(IR%XTG,2),JP) = ZTG_LS(:,4:SIZE(IR%XTG,2),JP)
     ENDIF
  END DO
ELSEIF(IO%CISBA=='2-L'.AND.IO%LTEMP_ARP) THEN
  DO JP=1,SIZE(IR%XWG,3)
     IR%XTG (:,3:SIZE(IR%XTG,2),JP) = ZTG_LS(:,3:SIZE(IR%XTG,2),JP)
  END DO
END IF
!
DEALLOCATE(ZZSFREEZE)
DEALLOCATE(ZWGI_CLIM_GRAD)
DEALLOCATE(ZWGTOT   )
DEALLOCATE(ZDW      )
!
!* masks where fields are not defined
WHERE (IR%XTG(:,1:SIZE(IR%XWG,2),:) == XUNDEF)
  IR%XWG (:,:,:) = XUNDEF
  IR%XWGI(:,:,:) = XUNDEF
END WHERE
!
!-------------------------------------------------------------------------------------
!
!*      1.4    Snow variables
!
!* vertical shift
IF (.NOT.LSNOW_IDEAL) THEN
  IF (IO%CISBA=='DIF') THEN
    IDEEP_SOIL = IO%NGROUND_LAYER
  ELSE
    IDEEP_SOIL = 2
  END IF        
  CALL PREP_VER_SNOW(IR%TSNOW,XZS_LS,PZS,ZTG_LS,IR%XTG,IDEEP_SOIL)
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
