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
USE MODD_SURF_PAR,         ONLY : XUNDEF
!
USE MODD_CSTS,             ONLY : XRHOLW
!
USE MODD_ISBA_n,           ONLY : CISBA, XPATCH, NGROUND_LAYER, &
                                    XDG, XWG, XWGI  
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
                                    XTWSNOW, XAVG_TWSNOW ,  &
                                    XTDSNOW, XAVG_TDSNOW ,  &
                                    XTTSNOW, XAVG_TTSNOW ,  &
                                    XSOIL_TSWI, XSOIL_TWG,  &
                                    XSOIL_TWGI  
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER                         :: JPATCH    ! tile loop counter
INTEGER                         :: JLAYER    ! layer loop counter
REAL, DIMENSION(SIZE(XPATCH,1)) :: ZSUMPATCH
REAL, DIMENSION(SIZE(XPATCH,1)) :: ZSUMDG, ZSNOW
INTEGER                         :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!       0.     Initialization
!              --------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_MISC_ISBA_N',0,ZHOOK_HANDLE)
IF (LSURF_MISC_BUDGET) THEN
!
  ZSUMPATCH(:) = 0.0
  DO JPATCH=1,SIZE(XPATCH,2)
    ZSUMPATCH(:) = ZSUMPATCH(:) + XPATCH(:,JPATCH)
  END DO
!
  ZSUMDG(:)=0.0
  ZSNOW (:)=0.0
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
  DO JPATCH=1,SIZE(XPATCH,2)
!
!cdir nodep
    DO JJ=1,SIZE(ZSUMPATCH)
!
      IF (ZSUMPATCH(JJ) > 0.) THEN
!
!     Halstead coefficient
!
        XAVG_HV(JJ) = XAVG_HV(JJ) + XPATCH(JJ,JPATCH) * XHV(JJ,JPATCH)
!
!     Snow fractions
!
        XAVG_PSNG(JJ) = XAVG_PSNG(JJ) + XPATCH(JJ,JPATCH) * XDPSNG(JJ,JPATCH)
        XAVG_PSNV(JJ) = XAVG_PSNV(JJ) + XPATCH(JJ,JPATCH) * XDPSNV(JJ,JPATCH)
        XAVG_PSN (JJ) = XAVG_PSN (JJ) + XPATCH(JJ,JPATCH) * XDPSN (JJ,JPATCH)
!
!     Flood fractions
!
        XAVG_FFG(JJ) = XAVG_FFG(JJ) + XPATCH(JJ,JPATCH) * XDFFG(JJ,JPATCH)
        XAVG_FFV(JJ) = XAVG_FFV(JJ) + XPATCH(JJ,JPATCH) * XDFFV(JJ,JPATCH)
        XAVG_FF (JJ) = XAVG_FF (JJ) + XPATCH(JJ,JPATCH) * XDFF (JJ,JPATCH)
!
!     Total albedo
!
        XAVG_ALBT(JJ) = XAVG_ALBT(JJ) + XPATCH(JJ,JPATCH) * XALBT (JJ,JPATCH)
!      
!     Snow total outputs
!
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
!   Soil Wetness Index profile
!
    DO JLAYER = 1, NGROUND_LAYER
!
!cdir nodep
      DO JJ=1,SIZE(ZSUMPATCH)
        IF (ZSUMPATCH(JJ) > 0.    ) THEN
          XAVG_SWI (JJ,JLAYER) = XAVG_SWI (JJ,JLAYER) + XPATCH(JJ,JPATCH) * XSWI (JJ,JLAYER,JPATCH)
          XAVG_TSWI(JJ,JLAYER) = XAVG_TSWI(JJ,JLAYER) + XPATCH(JJ,JPATCH) * XTSWI(JJ,JLAYER,JPATCH)
        ENDIF
      ENDDO
 
    ENDDO
!
!   Total Soil Wetness Index
!
!   Total Soil Water Content (Liquid+Solid) and Total Frozen Content
!
    IF(CISBA=='DIF')THEN!     
!cdir nodep
      DO JJ=1,SIZE(XPATCH,1)
        XSOIL_TSWI(JJ) = XSOIL_TSWI(JJ) + XPATCH(JJ,JPATCH) * XTSWI(JJ,1,JPATCH) * XDG (JJ,1,JPATCH)
        XSOIL_TWG (JJ) = XSOIL_TWG (JJ) + XPATCH(JJ,JPATCH) * (XWG (JJ,1,JPATCH)+XWGI(JJ,1,JPATCH)) * XDG (JJ,1,JPATCH)
        XSOIL_TWGI(JJ) = XSOIL_TWGI(JJ) + XPATCH(JJ,JPATCH) * XWGI(JJ,1,JPATCH)                    * XDG (JJ,1,JPATCH)
      ENDDO
      DO JLAYER = 2, NGROUND_LAYER   
!cdir nodep 
        DO JJ=1,SIZE(XPATCH,1)
          IF (ZSUMPATCH(JJ) > 0.) THEN 
            XSOIL_TSWI(JJ) = XSOIL_TSWI(JJ) + XPATCH(JJ,JPATCH) * XTSWI(JJ,JLAYER,JPATCH) &
                                              * (XDG (JJ,JLAYER,JPATCH)-XDG(JJ,JLAYER-1,JPATCH))  
            XSOIL_TWG (JJ) = XSOIL_TWG (JJ) + XPATCH(JJ,JPATCH) * (XWG(JJ,JLAYER,JPATCH)+  &
                       XWGI(JJ,JLAYER,JPATCH)) * (XDG (JJ,JLAYER,JPATCH)-XDG(JJ,JLAYER-1,JPATCH))  
            XSOIL_TWGI(JJ) = XSOIL_TWGI(JJ) + XPATCH(JJ,JPATCH) * XWGI(JJ,JLAYER,JPATCH)   &
                                              * (XDG (JJ,JLAYER,JPATCH)-XDG(JJ,JLAYER-1,JPATCH))    
          ENDIF
        ENDDO       
      ENDDO
    ELSE
!cdir nodep
      DO JJ=1,SIZE(XPATCH,1)
        XSOIL_TSWI(JJ) = XSOIL_TSWI(JJ) + XPATCH(JJ,JPATCH) * XTSWI(JJ,2,JPATCH) * XDG (JJ,2,JPATCH)
        XSOIL_TWG (JJ) = XSOIL_TWG (JJ) + XPATCH(JJ,JPATCH) * (XWG(JJ,2,JPATCH)+XWGI(JJ,2,JPATCH)) * XDG (JJ,2,JPATCH)
        XSOIL_TWGI(JJ) = XSOIL_TWGI(JJ) + XPATCH(JJ,JPATCH) * XWGI(JJ,2,JPATCH)                    * XDG (JJ,2,JPATCH)
      ENDDO
      IF(CISBA=='3-L')THEN
!cdir nodep   
        DO JJ=1,SIZE(XPATCH,1)        
          XSOIL_TSWI(JJ) = XSOIL_TSWI(JJ) + XPATCH(JJ,JPATCH) * XTSWI(JJ,3,JPATCH) * &
                  (XDG(JJ,3,JPATCH)-XDG(JJ,2,JPATCH))  
          ! no ice in the third layer of 3-L
          XSOIL_TWG (JJ) = XSOIL_TWG (JJ) + XPATCH(JJ,JPATCH) * (XWG(JJ,3,JPATCH)+XWGI(JJ,3,JPATCH)) * &
                  (XDG(JJ,3,JPATCH)-XDG(JJ,2,JPATCH))  
        ENDDO
      ENDIF
    ENDIF   
    
    WHERE (ZSUMPATCH(:) > 0.)
        ZSUMDG   (:) = ZSUMDG    (:) + XPATCH(:,JPATCH) * XDG(:,NGROUND_LAYER,JPATCH)
    END WHERE     

  END DO
!
! Snow temperature
!  
  WHERE(ZSNOW(:)>0.0)
        XAVG_TTSNOW(:) = XAVG_TTSNOW(:)/ZSNOW(:)
  ELSEWHERE
        XAVG_TTSNOW(:) = XUNDEF
  ENDWHERE
!
! Total Soil Wetness Index
!
  WHERE(ZSUMDG(:)>0.0)XSOIL_TSWI(:) = XSOIL_TSWI(:)/ZSUMDG(:)
!
! Total Soil Water Content (Liquid+Solid) and Total Frozen Content (kg/m2)
!
  XSOIL_TWG (:)= XSOIL_TWG (:) * XRHOLW
  XSOIL_TWGI(:)= XSOIL_TWGI(:) * XRHOLW
!
END IF
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_MISC_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE_DIAG_MISC_ISBA_n
