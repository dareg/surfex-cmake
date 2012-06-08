!     #############################
      SUBROUTINE AVERAGE_DIAG_MISC_ISBA_n
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
!!	P. Le Moigne           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2004
!!      B. Decharme  2008    New diag Total albedo, Total SWI, & Flood
!!      B. Decharme 09/2009  New diag Total soil SWI
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR,         ONLY : XUNDEF, NUNDEF
!
USE MODD_CSTS,             ONLY : XRHOLW
!
USE MODD_ISBA_n,           ONLY : CISBA, XPATCH, NGROUND_LAYER, &
                                  XDG, XDZG, XWG, XWGI,         &
                                  NSIZE_NATURE_P, NWG_LAYER,    &
                                  XDG2
!
USE MODD_DIAG_MISC_ISBA_n, ONLY : LSURF_MISC_BUDGET,      &
                                    XHV  , XAVG_HV   ,      &
                                    XSWI , XAVG_SWI  ,      &
                                    XTSWI, XAVG_TSWI ,      &
                                    XDPSNG, XAVG_PSNG ,     &
                                    XDPSNV, XAVG_PSNV ,     &
                                    XDPSN , XAVG_PSN  ,     &
                                    XALBT, XAVG_ALBT ,      &
                                    XDFFG , XAVG_FFG  ,     &
                                    XDFFV , XAVG_FFV  ,     &
                                    XDFF  , XAVG_FF   ,     &
                                    XDFSAT, XAVG_FSAT   ,   &
                                    XTWSNOW, XAVG_TWSNOW ,  &
                                    XTDSNOW, XAVG_TDSNOW ,  &
                                    XTTSNOW, XAVG_TTSNOW ,  &
                                    XSOIL_TSWI, XSOIL_TWG,  &
                                    XSOIL_TWGI, XSURF_TSWI, &
                                    XSURF_TWG, XSURF_TWGI,  &
                                    XROOT_TSWI, XROOT_TWG,  &
                                    XROOT_TWGI, XDEEP_TSWI, &
                                    XDEEP_TWG, XDEEP_TWGI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER                         :: JJ        ! grid-cell loop counter
INTEGER                         :: JPATCH    ! tile loop counter
INTEGER                         :: JLAYER    ! layer loop counter
REAL, DIMENSION(SIZE(XPATCH,1)) :: ZSUMPATCH
REAL, DIMENSION(SIZE(XPATCH,1)) :: ZSUMDG, ZSNOW, ZSUMSURF, ZSUMROOT, ZSUMDEEP
REAL                            :: ZWORK
INTEGER                         :: INI,INP,IDEPTH,IWORK
!
REAL, DIMENSION(SIZE(XPATCH,1)) :: ZPOND
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
IF (.NOT.LSURF_MISC_BUDGET) THEN
   IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_MISC_ISBA_N',1,ZHOOK_HANDLE)
   RETURN
ENDIF
!
INI=SIZE(XPATCH,1)
INP=SIZE(XPATCH,2)
!
ZSUMPATCH(:) = 0.0
DO JPATCH=1,INP
   DO JJ=1,INI
      ZSUMPATCH(JJ) = ZSUMPATCH(JJ) + XPATCH(JJ,JPATCH)
   END DO
END DO
!
ZSUMSURF(:)=0.0
ZSUMROOT(:)=0.0
ZSUMDEEP(:)=0.0
ZSUMDG  (:)=0.0
ZSNOW   (:)=0.0
!
!-------------------------------------------------------------------------------
!
!       1.     Surface Miscellaneous terms
!              ---------------------------
!
XAVG_HV  (:)   = 0.
XAVG_PSNG(:)   = 0.
XAVG_PSNV(:)   = 0.
XAVG_PSN (:)   = 0.
XAVG_ALBT(:)   = 0.
XAVG_SWI (:,:) = 0.
XAVG_TSWI(:,:) = 0.
XAVG_FSAT(:)   = 0.
XAVG_FFG (:)   = 0.
XAVG_FFV (:)   = 0.
XAVG_FF  (:)   = 0.
XAVG_TWSNOW(:) = 0.
XAVG_TDSNOW(:) = 0.  
XAVG_TTSNOW(:) = 0.
!   
XSOIL_TSWI(:)  = 0.
XSOIL_TWG  (:) = 0.
XSOIL_TWGI (:) = 0.
!   
IF(CISBA=='DIF')THEN
!        
  XSURF_TSWI (:) = 0.
  XSURF_TWG  (:) = 0.
  XSURF_TWGI (:) = 0.
!  
  XROOT_TSWI (:) = 0.
  XROOT_TWG  (:) = 0.
  XROOT_TWGI (:) = 0.
!   
  XDEEP_TSWI (:) = 0.
  XDEEP_TWG  (:) = 0.
  XDEEP_TWGI (:) = 0.
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
      XAVG_HV(JJ) = XAVG_HV(JJ) + XPATCH(JJ,JPATCH) * XHV(JJ,JPATCH)
!
!     Snow fractions
      XAVG_PSNG(JJ) = XAVG_PSNG(JJ) + XPATCH(JJ,JPATCH) * XDPSNG(JJ,JPATCH)
      XAVG_PSNV(JJ) = XAVG_PSNV(JJ) + XPATCH(JJ,JPATCH) * XDPSNV(JJ,JPATCH)
      XAVG_PSN (JJ) = XAVG_PSN (JJ) + XPATCH(JJ,JPATCH) * XDPSN (JJ,JPATCH)
!
!     Saturated fraction
      XAVG_FSAT (JJ) = XAVG_FSAT (JJ) + XPATCH(JJ,JPATCH) * XDFSAT (JJ,JPATCH)
!
!     Flood fractions
      XAVG_FFG(JJ) = XAVG_FFG(JJ) + XPATCH(JJ,JPATCH) * XDFFG(JJ,JPATCH)
      XAVG_FFV(JJ) = XAVG_FFV(JJ) + XPATCH(JJ,JPATCH) * XDFFV(JJ,JPATCH)
      XAVG_FF (JJ) = XAVG_FF (JJ) + XPATCH(JJ,JPATCH) * XDFF (JJ,JPATCH)
!
!     Total albedo
      XAVG_ALBT(JJ) = XAVG_ALBT(JJ) + XPATCH(JJ,JPATCH) * XALBT (JJ,JPATCH)
!      
!     Snow total outputs
      XAVG_TWSNOW(JJ) = XAVG_TWSNOW(JJ) + XPATCH(JJ,JPATCH) * XTWSNOW(JJ,JPATCH)
      XAVG_TDSNOW(JJ) = XAVG_TDSNOW(JJ) + XPATCH(JJ,JPATCH) * XTDSNOW(JJ,JPATCH)
!      
      IF (XTWSNOW(JJ,JPATCH)>0.0) THEN
         XAVG_TTSNOW(JJ) = XAVG_TTSNOW(JJ) + XPATCH(JJ,JPATCH) * XTTSNOW(JJ,JPATCH)
         ZSNOW      (JJ) = ZSNOW      (JJ) + XPATCH(JJ,JPATCH)
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
IF(CISBA=='DIF')THEN ! DIF case
!---------------------------------------------
!
! Ponderation for diag by layer
  ZPOND(:)=0.0
  DO JPATCH=1,INP
    IF (NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
!cdir nodep 
    DO JJ=1,INI
      ZPOND(JJ)=ZPOND(JJ)+XPATCH(JJ,JPATCH)
    ENDDO
  ENDDO
!
  DO JPATCH=1,INP
!  
    IF (NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
!     
!   Surface
    XSURF_TSWI(:) = XSURF_TSWI(:) + XPATCH(:,JPATCH) * XTSWI(:,1,JPATCH)                   * XDZG(:,1,JPATCH)
    XSURF_TWG (:) = XSURF_TWG (:) + XPATCH(:,JPATCH) * (XWG (:,1,JPATCH)+XWGI(:,1,JPATCH)) * XDZG(:,1,JPATCH)
    XSURF_TWGI(:) = XSURF_TWGI(:) + XPATCH(:,JPATCH) * XWGI (:,1,JPATCH)                   * XDZG(:,1,JPATCH)
!
    DO JLAYER = 1,NGROUND_LAYER
!cdir nodep 
      DO JJ=1,INI
        IDEPTH=NWG_LAYER(JJ,JPATCH)
        IF(JLAYER<=IDEPTH.AND.IDEPTH/=NUNDEF)THEN
          !      
          XAVG_SWI (JJ,JLAYER) = XAVG_SWI (JJ,JLAYER)+XPATCH(JJ,JPATCH)*XSWI (JJ,JLAYER,JPATCH) 
          XAVG_TSWI(JJ,JLAYER) = XAVG_TSWI(JJ,JLAYER)+XPATCH(JJ,JPATCH)*XTSWI(JJ,JLAYER,JPATCH)
          !
          ZWORK=XDZG(JJ,JLAYER,JPATCH)
          XSOIL_TSWI(JJ) = XSOIL_TSWI(JJ) + ZWORK * XPATCH(JJ,JPATCH) * XTSWI(JJ,JLAYER,JPATCH)
          XSOIL_TWG (JJ) = XSOIL_TWG (JJ) + ZWORK * XPATCH(JJ,JPATCH) * (XWG(JJ,JLAYER,JPATCH)+XWGI(JJ,JLAYER,JPATCH))
          XSOIL_TWGI(JJ) = XSOIL_TWGI(JJ) + ZWORK * XPATCH(JJ,JPATCH) * XWGI(JJ,JLAYER,JPATCH)   
          !
          ZWORK=MIN(XDZG(JJ,JLAYER,JPATCH),MAX(0.0,XDG2(JJ,JPATCH)-XDG(JJ,JLAYER,JPATCH)+XDZG(JJ,JLAYER,JPATCH)))
          XROOT_TSWI (JJ) = XROOT_TSWI (JJ) + ZWORK * XPATCH(JJ,JPATCH) * XTSWI(JJ,JLAYER,JPATCH)
          XROOT_TWG  (JJ) = XROOT_TWG  (JJ) + ZWORK * XPATCH(JJ,JPATCH) * (XWG(JJ,JLAYER,JPATCH)+XWGI(JJ,JLAYER,JPATCH))
          XROOT_TWGI (JJ) = XROOT_TWGI (JJ) + ZWORK * XPATCH(JJ,JPATCH) * XWGI(JJ,JLAYER,JPATCH)
          !
          ZWORK=MIN(XDZG(JJ,JLAYER,JPATCH),MAX(0.0,XDG(JJ,JLAYER,JPATCH)-XDG2(JJ,JPATCH)))
          XDEEP_TSWI (JJ) = XDEEP_TSWI (JJ) + ZWORK * XPATCH(JJ,JPATCH) * XTSWI(JJ,JLAYER,JPATCH)
          XDEEP_TWG  (JJ) = XDEEP_TWG  (JJ) + ZWORK * XPATCH(JJ,JPATCH) * (XWG(JJ,JLAYER,JPATCH)+XWGI(JJ,JLAYER,JPATCH))
          XDEEP_TWGI (JJ) = XDEEP_TWGI (JJ) + ZWORK * XPATCH(JJ,JPATCH) * XWGI(JJ,JLAYER,JPATCH)
          !
        ENDIF
      ENDDO
!
    XAVG_SWI (:,JLAYER) = XAVG_SWI (:,JLAYER) / ZPOND(:)
    XAVG_TSWI(:,JLAYER) = XAVG_TSWI(:,JLAYER) / ZPOND(:) 
!
    ENDDO
!
!cdir nodep
    DO JJ=1,INI      
      IDEPTH=NWG_LAYER(JJ,JPATCH)
      IF(IDEPTH/=NUNDEF)THEN
        ZSUMSURF (JJ) = ZSUMSURF (JJ) + XPATCH(JJ,JPATCH) * XDG(JJ,1,JPATCH)
        ZSUMROOT (JJ) = ZSUMROOT (JJ) + XPATCH(JJ,JPATCH) * XDG2(JJ,JPATCH)
        ZSUMDEEP (JJ) = ZSUMDEEP (JJ) + XPATCH(JJ,JPATCH) * (XDG(JJ,IDEPTH,JPATCH)-XDG2(JJ,JPATCH))
        ZSUMDG   (JJ) = ZSUMDG   (JJ) + XPATCH(JJ,JPATCH) * XDG(JJ,IDEPTH,JPATCH)
      ENDIF
    ENDDO
!
  ENDDO
!
  WHERE(ZSUMSURF(:)>0.0) XSURF_TSWI (:) = XSURF_TSWI (:) / ZSUMSURF(:)
  WHERE(ZSUMROOT(:)>0.0) XROOT_TSWI (:) = XROOT_TSWI (:) / ZSUMROOT(:)
  WHERE(ZSUMDEEP(:)>0.0) XDEEP_TSWI (:) = XDEEP_TSWI (:) / ZSUMDEEP(:)
!
  XSURF_TWG  (:) = XSURF_TWG  (:) * XRHOLW
  XSURF_TWGI (:) = XSURF_TWGI (:) * XRHOLW
  XROOT_TWG  (:) = XROOT_TWG  (:) * XRHOLW
  XROOT_TWGI (:) = XROOT_TWGI (:) * XRHOLW 
  XDEEP_TWG  (:) = XDEEP_TWG  (:) * XRHOLW
  XDEEP_TWGI (:) = XDEEP_TWGI (:) * XRHOLW 
!
!---------------------------------------------
ELSE ! Force-restore case
!---------------------------------------------
! 
  DO JPATCH=1,INP
!  
     IF (NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
!
!    Soil Wetness Index profile
     DO JLAYER = 1,NGROUND_LAYER
!cdir nodep
        DO JJ=1,INI
           XAVG_SWI (JJ,JLAYER) = XAVG_SWI (JJ,JLAYER) + XPATCH(JJ,JPATCH) * XSWI (JJ,JLAYER,JPATCH)
           XAVG_TSWI(JJ,JLAYER) = XAVG_TSWI(JJ,JLAYER) + XPATCH(JJ,JPATCH) * XTSWI(JJ,JLAYER,JPATCH)
        ENDDO 
      ENDDO
!
!cdir nodep
      DO JJ=1,INI
        XSOIL_TSWI(JJ) = XSOIL_TSWI(JJ) + XPATCH(JJ,JPATCH) * XTSWI(JJ,2,JPATCH) * XDG (JJ,2,JPATCH)
        XSOIL_TWG (JJ) = XSOIL_TWG (JJ) + XPATCH(JJ,JPATCH) * (XWG(JJ,2,JPATCH)+XWGI(JJ,2,JPATCH))* XDG (JJ,2,JPATCH)
        XSOIL_TWGI(JJ) = XSOIL_TWGI(JJ) + XPATCH(JJ,JPATCH) * XWGI(JJ,2,JPATCH)* XDG (JJ,2,JPATCH)
      ENDDO
      IF(CISBA=='3-L')THEN
!cdir nodep   
        DO JJ=1,SIZE(XPATCH,1)        
          XSOIL_TSWI(JJ) = XSOIL_TSWI(JJ) + XPATCH(JJ,JPATCH) * XTSWI(JJ,3,JPATCH)* (XDG(JJ,3,JPATCH)-XDG(JJ,2,JPATCH))  
          ! no ice in the third layer of 3-L
          XSOIL_TWG (JJ) = XSOIL_TWG (JJ) + XPATCH(JJ,JPATCH) * XWG(JJ,3,JPATCH)* (XDG(JJ,3,JPATCH)-XDG(JJ,2,JPATCH))  
        ENDDO
      ENDIF
!
!cdir nodep
      DO JJ=1,INI      
         ZSUMDG (JJ) = ZSUMDG (JJ) + XPATCH(JJ,JPATCH) * XDG(JJ,NGROUND_LAYER,JPATCH)
      ENDDO
!      
  ENDDO
!  
!---------------------------------------------
ENDIF   
!---------------------------------------------
!
!       3.     Final computation for grid-cell diag
!              ------------------------------------
!
!Total Soil Wetness Index
WHERE(ZSUMDG(:)>0.0)XSOIL_TSWI(:) = XSOIL_TSWI(:)/ZSUMDG(:)
        !Total Soil Water Content (Liquid+Solid) and Total Frozen Content (kg/m2)
XSOIL_TWG (:)= XSOIL_TWG (:) * XRHOLW
XSOIL_TWGI(:)= XSOIL_TWGI(:) * XRHOLW
!
! Snow temperature  
WHERE(ZSNOW(:)>0.0)
      XAVG_TTSNOW(:) = XAVG_TTSNOW(:)/ZSNOW(:)
ELSEWHERE
      XAVG_TTSNOW(:) = XUNDEF
ENDWHERE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_MISC_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE_DIAG_MISC_ISBA_n
