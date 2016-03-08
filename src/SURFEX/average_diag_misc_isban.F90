!     #############################
      SUBROUTINE AVERAGE_DIAG_MISC_ISBA_n (DGMI, DGMIP, IO, IP, IMX, PLAI, IR)
!     #############################
!
!
!!****  *AVERAGE_DIAG_MISC_ISBA_n*  
!!
!!    PURPOSE
!!    -------
!      Average the cumulated diagnostics from all ISBA tiles
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!      P. Le Moigne           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2004
!!      B. Decharme  2008    New diag Total albedo, Total SWI, & Flood
!!      B. Decharme 09/2009  New diag Total soil SWI
!!      B. Decharme  2012    Averaged LAI
!!      B. Decharme  2012    New diag for DIF:
!!                           F2 stress
!!                           Root zone swi, wg and wgi
!!                           swi, wg and wgi comparable to ISBA-FR-DG2 and DG3 layers
!!                           active layer thickness over permafrost
!!                           frozen layer thickness over non-permafrost
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t, DIAG_MISC_ISBA_PATCH_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODD_SURF_PAR,         ONLY : XUNDEF, NUNDEF
!
USE MODD_CSTS,             ONLY : XRHOLW
!
!
!
USE MODI_COMPUT_COLD_LAYERS_THICK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
TYPE(DIAG_MISC_ISBA_PATCH_t), INTENT(INOUT) :: DGMIP
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMX
REAL, DIMENSION(:,:), INTENT(INOUT) :: PLAI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
INTEGER                         :: JJ        ! grid-cell loop counter
INTEGER                         :: JPATCH    ! tile loop counter
INTEGER                         :: JLAYER    ! layer loop counter
REAL, DIMENSION(SIZE(IP%XPATCH,1)) :: ZSUMPATCH
REAL, DIMENSION(SIZE(IP%XPATCH,1)) :: ZSUMDG, ZSNOW, ZSUMFRD2, ZSUMFRD3, ZPONDF2
REAL, DIMENSION(SIZE(IP%XPATCH,1),SIZE(IP%XPATCH,2)) :: ZLAI
REAL                            :: ZWORK
INTEGER                         :: INI,INP,IDEPTH,IWORK
!
REAL, DIMENSION(SIZE(IMX%XDG,1),SIZE(IMX%XDG,2)) :: ZPOND, ZTG, ZDG
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!       0.     Initialization
!              --------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_MISC_ISBA_N',0,ZHOOK_HANDLE)
!
IF (.NOT.DGMI%LSURF_MISC_BUDGET) THEN
   IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_MISC_ISBA_N',1,ZHOOK_HANDLE)
   RETURN
ENDIF
!
INI=SIZE(IP%XPATCH,1)
INP=SIZE(IP%XPATCH,2)
!
ZSUMPATCH(:) = 0.0
DO JPATCH=1,INP
   DO JJ=1,INI
      ZSUMPATCH(JJ) = ZSUMPATCH(JJ) + IP%XPATCH(JJ,JPATCH)
   END DO
END DO
!
ZSUMFRD2(:)=0.0
ZSUMFRD3(:)=0.0
ZSUMDG  (:)=0.0
ZSNOW   (:)=0.0
ZPONDF2 (:)=0.0
!
WHERE(PLAI(:,:)/=XUNDEF)
      ZLAI(:,:)=PLAI(:,:)
ELSEWHERE
      ZLAI(:,:)=0.0
ENDWHERE
!
!-------------------------------------------------------------------------------
!
!       1.     Surface Miscellaneous terms
!              ---------------------------
!
DGMI%XHV  (:)   = 0.
DGMI%XPSNG(:)   = 0.
DGMI%XPSNV(:)   = 0.
DGMI%XPSN (:)   = 0.
DGMI%XSWI (:,:) = 0.
DGMI%XTSWI(:,:) = 0.
DGMI%XFSAT(:)   = 0.
DGMI%XFFG (:)   = 0.
DGMI%XFFV (:)   = 0.
DGMI%XFF  (:)   = 0.
DGMI%XTWSNOW(:) = 0.
DGMI%XTDSNOW(:) = 0.  
DGMI%XTTSNOW(:) = 0.
DGMI%XLAI   (:) = 0.
!   
DGMI%XSOIL_SWI  (:)  = 0.
DGMI%XSOIL_TSWI (:)  = 0.
DGMI%XSOIL_TWG  (:) = 0.
DGMI%XSOIL_TWGI (:) = 0.
DGMI%XSOIL_WG   (:) = 0.
DGMI%XSOIL_WGI  (:) = 0.
! 
IF(IO%CISBA=='DIF')THEN
!        
  DGMI%XALT   (:) = 0. 
  DGMI%XFLT   (:) = 0. 
!
ENDIF

IF(IO%CISBA=='DIF'.AND.DGMI%LSURF_MISC_DIF)THEN
!   
  DGMI%XFRD2_TSWI (:) = 0.
  DGMI%XFRD2_TWG  (:) = 0.
  DGMI%XFRD2_TWGI (:) = 0.
!   
  DGMI%XFRD3_TSWI (:) = 0.
  DGMI%XFRD3_TWG  (:) = 0.
  DGMI%XFRD3_TWGI (:) = 0.
!  
ENDIF
!
DO JPATCH=1,INP
!
!cdir nodep
  DO JJ=1,INI
!
    IF (ZSUMPATCH(JJ) > 0.) THEN
!
!     Halstead coefficient
      DGMI%XHV(JJ) = DGMI%XHV(JJ) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XHV(JJ)
!
!     Snow fractions
      DGMI%XPSNG(JJ) = DGMI%XPSNG(JJ) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XPSNG(JJ)
      DGMI%XPSNV(JJ) = DGMI%XPSNV(JJ) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XPSNV(JJ)
      DGMI%XPSN (JJ) = DGMI%XPSN (JJ) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XPSN (JJ)
!
!     Saturated fraction
      DGMI%XFSAT (JJ) = DGMI%XFSAT (JJ) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XFSAT (JJ)
!
!     Flood fractions
      DGMI%XFFG(JJ) = DGMI%XFFG(JJ) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XFFG(JJ)
      DGMI%XFFV(JJ) = DGMI%XFFV(JJ) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XFFV(JJ)
      DGMI%XFF (JJ) = DGMI%XFF (JJ) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XFF (JJ)
!
!     Total LAI
      DGMI%XLAI (JJ) = DGMI%XLAI(JJ)  + IP%XPATCH(JJ,JPATCH) * ZLAI (JJ,JPATCH)
!      
!     Snow total outputs
      DGMI%XTWSNOW(JJ) = DGMI%XTWSNOW(JJ) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XTWSNOW(JJ)
      DGMI%XTDSNOW(JJ) = DGMI%XTDSNOW(JJ) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XTDSNOW(JJ)
!      
      IF (DGMIP%AL(JPATCH)%XTWSNOW(JJ)>0.0) THEN
         DGMI%XTTSNOW(JJ) = DGMI%XTTSNOW(JJ) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XTTSNOW(JJ)
         ZSNOW      (JJ) = ZSNOW      (JJ) + IP%XPATCH(JJ,JPATCH)
      ENDIF
!
    ENDIF
!
  ENDDO
!
ENDDO
!
!-------------------------------------------------------------------------------
!
!       2.     Specific treatement following CISBA option
!              ------------------------------------------
!
!   Soil Wetness Index profile, Total Soil Wetness Index and 
!   Total Soil Water Content (Liquid+Solid) and Total Frozen Content
!
!---------------------------------------------
IF(IO%CISBA=='DIF')THEN ! DIF case
!---------------------------------------------
!
! Active and Frozen layers thickness
  ZTG(:,:)=0.0
  ZDG(:,:)=0.0
  DO JPATCH=1,INP
     DO JLAYER=1,IO%NGROUND_LAYER
        DO JJ=1,INI 
           ZTG(JJ,JLAYER) = ZTG(JJ,JLAYER) + IP%XPATCH(JJ,JPATCH) * IR%XTG(JJ,JLAYER,JPATCH)
           ZDG(JJ,JLAYER) = ZDG(JJ,JLAYER) + IP%XPATCH(JJ,JPATCH) * IMX%XDG(JJ,JLAYER,JPATCH)
        ENDDO
     ENDDO
  ENDDO
  CALL COMPUT_COLD_LAYERS_THICK(ZDG,ZTG,DGMI%XALT,DGMI%XFLT)
!    
  ZPOND(:,:)=0.0
  DO JPATCH=1,INP      
     IF(IP%NSIZE_NATURE_P(JPATCH) > 0 )THEN
       DO JLAYER = 1,IO%NGROUND_LAYER
!         cdir nodep 
          DO JJ=1,INI
             IDEPTH=IMX%NWG_LAYER(JJ,JPATCH)
             IF(JLAYER<=IDEPTH.AND.IDEPTH/=NUNDEF)THEN
               ZWORK=IP%XDZG(JJ,JLAYER,JPATCH)
               !Soil Wetness Index profile
               DGMI%XSWI (JJ,JLAYER) = DGMI%XSWI (JJ,JLAYER)+ZWORK*IP%XPATCH(JJ,JPATCH)*DGMIP%AL(JPATCH)%XSWI (JJ,JLAYER) 
               DGMI%XTSWI(JJ,JLAYER) = DGMI%XTSWI(JJ,JLAYER)+ZWORK*IP%XPATCH(JJ,JPATCH)*DGMIP%AL(JPATCH)%XTSWI(JJ,JLAYER)
               ZPOND     (JJ,JLAYER) = ZPOND     (JJ,JLAYER)+ZWORK*IP%XPATCH(JJ,JPATCH)
               !Total soil wetness index, total water and ice contents
               DGMI%XSOIL_SWI (JJ) = DGMI%XSOIL_SWI (JJ) + ZWORK * IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XSWI (JJ,JLAYER)
               DGMI%XSOIL_TSWI(JJ) = DGMI%XSOIL_TSWI(JJ) + ZWORK * IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XTSWI(JJ,JLAYER)
               ZSUMDG         (JJ) = ZSUMDG         (JJ) + ZWORK * IP%XPATCH(JJ,JPATCH)
               DGMI%XSOIL_TWG (JJ) = DGMI%XSOIL_TWG (JJ) + ZWORK * IP%XPATCH(JJ,JPATCH) * (IR%XWG(JJ,JLAYER,JPATCH) &
                                                         + IR%XWGI(JJ,JLAYER,JPATCH))
               DGMI%XSOIL_TWGI(JJ) = DGMI%XSOIL_TWGI(JJ) + ZWORK * IP%XPATCH(JJ,JPATCH) * IR%XWGI(JJ,JLAYER,JPATCH)
             ENDIF
          ENDDO
       ENDDO
     ENDIF
  ENDDO
!  
  WHERE(ZPOND(:,:)> 0.)
        DGMI%XSWI (:,:) = DGMI%XSWI (:,:) / ZPOND(:,:)
        DGMI%XTSWI(:,:) = DGMI%XTSWI(:,:) / ZPOND(:,:)
  ELSEWHERE
        DGMI%XSWI (:,:) = XUNDEF
        DGMI%XTSWI(:,:) = XUNDEF
  ENDWHERE
!
! ---------------------------------------------
  IF(DGMI%LSURF_MISC_DIF)THEN ! LSURF_MISC_DIF case
! ---------------------------------------------
!    
    DO JPATCH=1,INP
!  
      IF (IP%NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
!     
      DO JLAYER = 1,IO%NGROUND_LAYER
!       cdir nodep 
        DO JJ=1,INI
          IDEPTH=IMX%NWG_LAYER(JJ,JPATCH)
          IF(JLAYER<=IDEPTH.AND.IDEPTH/=NUNDEF)THEN
            !
            ! ISBA-FR-DG2 comparable soil wetness index, liquid water and ice contents
            ZWORK=MIN(IP%XDZG(JJ,JLAYER,JPATCH),MAX(0.0,IMX%XDG2(JJ,JPATCH)-IMX%XDG(JJ,JLAYER,JPATCH)&
                    +IP%XDZG(JJ,JLAYER,JPATCH)))
            DGMI%XFRD2_TSWI (JJ) = DGMI%XFRD2_TSWI (JJ) + ZWORK * IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XTSWI(JJ,JLAYER)
            DGMI%XFRD2_TWG  (JJ) = DGMI%XFRD2_TWG  (JJ) + ZWORK * IP%XPATCH(JJ,JPATCH) * IR%XWG  (JJ,JLAYER,JPATCH)
            DGMI%XFRD2_TWGI (JJ) = DGMI%XFRD2_TWGI (JJ) + ZWORK * IP%XPATCH(JJ,JPATCH) * IR%XWGI (JJ,JLAYER,JPATCH)
            ZSUMFRD2        (JJ) = ZSUMFRD2        (JJ) + ZWORK * IP%XPATCH(JJ,JPATCH)
            !
            ! ISBA-FR-DG3 comparable soil wetness index, liquid water and ice contents
            ZWORK=MIN(IP%XDZG(JJ,JLAYER,JPATCH),MAX(0.0,IMX%XDG(JJ,JLAYER,JPATCH)-IMX%XDG2(JJ,JPATCH)))
            DGMI%XFRD3_TSWI (JJ) = DGMI%XFRD3_TSWI (JJ) + ZWORK * IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XTSWI(JJ,JLAYER)
            DGMI%XFRD3_TWG  (JJ) = DGMI%XFRD3_TWG  (JJ) + ZWORK * IP%XPATCH(JJ,JPATCH) * IR%XWG  (JJ,JLAYER,JPATCH)
            DGMI%XFRD3_TWGI (JJ) = DGMI%XFRD3_TWGI (JJ) + ZWORK * IP%XPATCH(JJ,JPATCH) * IR%XWGI (JJ,JLAYER,JPATCH)
            ZSUMFRD3        (JJ) = ZSUMFRD3        (JJ) + ZWORK * IP%XPATCH(JJ,JPATCH)
            !
          ENDIF
        ENDDO
      ENDDO
!
    ENDDO
!    
    WHERE(ZSUMFRD2(:)>0.0) 
          DGMI%XFRD2_TSWI (:) = DGMI%XFRD2_TSWI (:) / ZSUMFRD2(:)
          DGMI%XFRD2_TWG  (:) = DGMI%XFRD2_TWG  (:) / ZSUMFRD2(:)
          DGMI%XFRD2_TWGI (:) = DGMI%XFRD2_TWGI (:) / ZSUMFRD2(:)          
    ELSEWHERE
          DGMI%XFRD2_TSWI (:) = XUNDEF
    ENDWHERE 
!    
    WHERE(ZSUMFRD3(:)>0.0) 
          DGMI%XFRD3_TSWI (:) = DGMI%XFRD3_TSWI (:) / ZSUMFRD3(:)
          DGMI%XFRD3_TWG  (:) = DGMI%XFRD3_TWG  (:) / ZSUMFRD3(:)
          DGMI%XFRD3_TWGI (:) = DGMI%XFRD3_TWGI (:) / ZSUMFRD3(:) 
    ELSEWHERE
          DGMI%XFRD3_TSWI (:) = XUNDEF
    ENDWHERE
!
! ---------------------------------------------
  ENDIF ! End LSURF_MISC_DIF case
! ---------------------------------------------
!
!---------------------------------------------
ELSE ! Force-restore case
!---------------------------------------------
! 
  DO JPATCH=1,INP
     DO JJ=1,INI     
        IF(ZSUMPATCH(JJ) > 0.)THEN
!
          DGMI%XSWI (JJ,1) = DGMI%XSWI (JJ,1) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XSWI (JJ,1)
          DGMI%XSWI (JJ,2) = DGMI%XSWI (JJ,2) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XSWI (JJ,2)
          DGMI%XTSWI(JJ,1) = DGMI%XTSWI(JJ,1) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XTSWI(JJ,1)
          DGMI%XTSWI(JJ,2) = DGMI%XTSWI(JJ,2) + IP%XPATCH(JJ,JPATCH) * DGMIP%AL(JPATCH)%XTSWI(JJ,2)
!
          DGMI%XSOIL_SWI (JJ) = DGMI%XSOIL_SWI (JJ) + &
                  IP%XPATCH(JJ,JPATCH) * IMX%XDG (JJ,2,JPATCH) * DGMIP%AL(JPATCH)%XSWI (JJ,2)
          DGMI%XSOIL_TSWI(JJ) = DGMI%XSOIL_TSWI(JJ) + IP%XPATCH(JJ,JPATCH) * &
                  IMX%XDG (JJ,2,JPATCH) * DGMIP%AL(JPATCH)%XTSWI(JJ,2)
          DGMI%XSOIL_TWG (JJ) = DGMI%XSOIL_TWG (JJ) + IP%XPATCH(JJ,JPATCH) * &
                  IMX%XDG (JJ,2,JPATCH) * (IR%XWG(JJ,2,JPATCH) + IR%XWGI(JJ,2,JPATCH))
          DGMI%XSOIL_TWGI(JJ) = DGMI%XSOIL_TWGI(JJ) + &
                  IP%XPATCH(JJ,JPATCH) * IMX%XDG (JJ,2,JPATCH) * IR%XWGI(JJ,2,JPATCH) 
! 
          ZSUMDG         (JJ) = ZSUMDG        (JJ) + IP%XPATCH(JJ,JPATCH) * IMX%XDG(JJ,IO%NGROUND_LAYER,JPATCH)        
!          
       ENDIF
     ENDDO
  ENDDO     
!     
  IF(IO%CISBA=='3-L')THEN
!          
    ZPOND(:,:)=0.0
    DO JPATCH=1,INP
       DO JJ=1,SIZE(IP%XPATCH,1)        
          IF(ZSUMPATCH(JJ) > 0.)THEN
!
            ZWORK=MAX(0.0,IMX%XDG(JJ,3,JPATCH)-IMX%XDG(JJ,2,JPATCH))
!
!           Remenber: no ice in the third layer of 3-L
            ZPOND          (JJ,3) = ZPOND          (JJ,3) + IP%XPATCH(JJ,JPATCH) * ZWORK
            DGMI%XSWI      (JJ,3) = DGMI%XSWI      (JJ,3) + IP%XPATCH(JJ,JPATCH) * ZWORK * DGMIP%AL(JPATCH)%XSWI (JJ,3)
            DGMI%XSOIL_SWI (JJ  ) = DGMI%XSOIL_SWI (JJ  ) + IP%XPATCH(JJ,JPATCH) * ZWORK * DGMIP%AL(JPATCH)%XSWI (JJ,3)  
            DGMI%XTSWI     (JJ,3) = DGMI%XTSWI     (JJ,3) + IP%XPATCH(JJ,JPATCH) * ZWORK * DGMIP%AL(JPATCH)%XTSWI(JJ,3)
            DGMI%XSOIL_TSWI(JJ  ) = DGMI%XSOIL_TSWI(JJ  ) + IP%XPATCH(JJ,JPATCH) * ZWORK * DGMIP%AL(JPATCH)%XTSWI(JJ,3)  
            DGMI%XSOIL_TWG (JJ  ) = DGMI%XSOIL_TWG (JJ  ) + IP%XPATCH(JJ,JPATCH) * ZWORK * IR%XWG  (JJ,3,JPATCH)  
!
          ENDIF
       ENDDO
    ENDDO
!
    WHERE(ZPOND(:,3)>0.0)
          DGMI%XSWI (:,3) = DGMI%XSWI (:,3) / ZPOND(:,3)
          DGMI%XTSWI(:,3) = DGMI%XTSWI(:,3) / ZPOND(:,3)
    ELSEWHERE
          DGMI%XSWI (:,3) = XUNDEF
          DGMI%XTSWI(:,3) = XUNDEF
    ENDWHERE
!
  ENDIF
  
!
!---------------------------------------------
ENDIF ! End ISBA soil scheme case   
!---------------------------------------------
!
!       3.     Final computation for grid-cell diag
!              ------------------------------------
!
!Total Soil Wetness Index and Soil Water Content (m3.m-3)
WHERE(ZSUMDG(:)>0.0)
      DGMI%XSOIL_SWI (:) = DGMI%XSOIL_SWI (:)/ZSUMDG(:)
      DGMI%XSOIL_TSWI(:) = DGMI%XSOIL_TSWI(:)/ZSUMDG(:)
      DGMI%XSOIL_WG  (:) = DGMI%XSOIL_TWG (:)/ZSUMDG(:)
      DGMI%XSOIL_WGI (:) = DGMI%XSOIL_TWGI(:)/ZSUMDG(:)
ENDWHERE
!       
!Total Soil Water Content (Liquid+Solid) and Total Frozen Content (kg/m2)
DGMI%XSOIL_TWG (:)= DGMI%XSOIL_TWG (:) * XRHOLW
DGMI%XSOIL_TWGI(:)= DGMI%XSOIL_TWGI(:) * XRHOLW
!
! Snow temperature  
WHERE(ZSNOW(:)>0.0)
      DGMI%XTTSNOW(:) = DGMI%XTTSNOW(:)/ZSNOW(:)
ELSEWHERE
      DGMI%XTTSNOW(:) = XUNDEF
ENDWHERE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_MISC_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE_DIAG_MISC_ISBA_n
