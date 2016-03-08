!     #########
SUBROUTINE DIAG_MISC_ISBA_n (DGMIP, DGMIK, IPK, INIK, IMXK, IRK, AGK, IO, OSURF_MISC_BUDGET, &
                             PTSTEP, HSNOW, OAGRIP, PTIME, KSIZE, KPATCH, KMASK     )  
!     ###############################################################################
!
!!****  *DIAG_MISC-ISBA_n * - additional diagnostics for ISBA
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
!!     P. Le Moigne 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2004
!!      Modified    10/2004 by P. Le Moigne: Halstead coefficient
!!      B. Decharme   2008    Do not limit the SWI to 1
!!                            Add total SWI
!!      S. Lafont    03/2009 : change unit of carbon output in kg/m2/s
!!      A.L. Gibelin 04/2009 : Add respiration diagnostics
!!      A.L. Gibelin 07/2009 : Suppress RDK and transform GPP as a diagnostic
!!        S. Lafont  01/2011 : accumulate carbon variable between 2 outputs
!!       B. Decharme 05/2012 : Carbon fluxes in diag_evap
!!       B. Decharme 05/2012 : Active and frozen layers thickness for dif
!!       B. Decharme 06/2013 : Snow temp for EBA scheme (XP_SNOWTEMP not allocated)
!!
!!------------------------------------------------------------------
!
!
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_PATCH_t, DIAG_MISC_ISBA_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
!
USE MODD_CSTS,       ONLY : XTT
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
!                                     
USE MODD_TYPE_SNOW
!
USE MODI_COMPUT_COLD_LAYERS_THICK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(DIAG_MISC_ISBA_PATCH_t), INTENT(INOUT) :: DGMIP
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMIK
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IPK
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INIK
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMXK
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IRK
TYPE(AGRI_t), INTENT(INOUT) :: AGK
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
!
LOGICAL, INTENT(IN) :: OSURF_MISC_BUDGET
REAL,               INTENT(IN)    :: PTSTEP        ! timestep for  accumulated values 
 CHARACTER(LEN=*), INTENT(IN)      :: HSNOW         ! snow scheme
LOGICAL, INTENT(IN)               :: OAGRIP
REAL,    INTENT(IN)               :: PTIME   ! current time since midnight
INTEGER, INTENT(IN)               :: KSIZE, KPATCH
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
!    
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(IRK%XPSN,1))    :: ZSNOWTEMP
REAL, DIMENSION(SIZE(IRK%TSNOW%WSNOW,1),SIZE(IRK%TSNOW%WSNOW,2)) :: ZWORK
REAL, DIMENSION(SIZE(IRK%TSNOW%WSNOW,1),SIZE(IRK%TSNOW%WSNOW,2)) :: ZWORKTEMP
!
REAL, DIMENSION(KSIZE) :: ZALT, ZFLT
!
LOGICAL :: GMASK
INTEGER :: JJ, JI, JK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_MISC_ISBA_N',0,ZHOOK_HANDLE)
!
IF (OSURF_MISC_BUDGET) THEN
  !
  DGMIK%XSWI (:,:)=XUNDEF
  DGMIK%XTSWI(:,:)=XUNDEF  
  DO JJ=1,SIZE(IRK%XWG,2)
    DO JI=1,SIZE(IRK%XWG,1)
      IF(IRK%XWG (JI,JJ,1)/=XUNDEF)THEN    
        DGMIK%XSWI (JI,JJ) = (IRK%XWG (JI,JJ,1) - IPK%XWWILT(JI,JJ)) / (IPK%XWFC(JI,JJ) - IPK%XWWILT(JI,JJ))
        DGMIK%XTSWI(JI,JJ) = (IRK%XWG (JI,JJ,1) - IPK%XWWILT(JI,JJ)) / (IPK%XWFC(JI,JJ) - IPK%XWWILT(JI,JJ))
      ENDIF
      IF(IRK%XWGI (JI,JJ,1)/=XUNDEF)THEN    
        DGMIK%XTSWI(JI,JJ) = DGMIK%XTSWI(JI,JJ) +  IRK%XWGI(JI,JJ,1) / (IPK%XWFC(JI,JJ) - IPK%XWWILT(JI,JJ))
      ENDIF
    ENDDO
  ENDDO
  !
  DO JK=1,SIZE(DGMIK%XSWI,2)
!cdir nodep
    DO JJ=1,KSIZE
      JI                      =  KMASK         (JJ)
      !
      DGMIP%AL(KPATCH)%XSWI     (JI,JK)  =  DGMIK%XSWI        (JJ,JK)
      DGMIP%AL(KPATCH)%XTSWI    (JI,JK)  =  DGMIK%XTSWI       (JJ,JK)
      !
    END DO
  ENDDO  
  !
  DO JI = 1,SIZE(IRK%TSNOW%WSNOW,2)
!cdir nodep 
    DO JJ = 1,SIZE(IRK%TSNOW%WSNOW,1)
      ZWORK(JJ,JI)  = IRK%TSNOW%WSNOW(JJ,JI,1) / IRK%TSNOW%RHO(JJ,JI,1)
    ENDDO
  ENDDO
  !
  DGMIK%XTWSNOW=0.
  DGMIK%XTDSNOW=0.
  ZSNOWTEMP=0.  
  !
  IF (HSNOW/='EBA')THEN
     ZWORKTEMP(:,:) = DGMIK%XSNOWTEMP(:,:)
  ELSE
     ZWORKTEMP(:,1) = MIN(IRK%XTG(:,1,1),XTT)
  ENDIF
  !
  DO JI = 1,SIZE(IRK%TSNOW%WSNOW,2)
!cdir nodep 
    DO JJ = 1,SIZE(IRK%TSNOW%WSNOW,1)
      DGMIK%XTWSNOW(JJ) = DGMIK%XTWSNOW(JJ) + IRK%TSNOW%WSNOW(JJ,JI,1)      
      DGMIK%XTDSNOW(JJ) = DGMIK%XTDSNOW(JJ) + ZWORK (JJ,JI)
      ZSNOWTEMP(JJ) = ZSNOWTEMP(JJ) + ZWORKTEMP(JJ,JI) * ZWORK(JJ,JI)
    ENDDO
  ENDDO
  !
  WHERE(DGMIK%XTDSNOW(:)>0.0)
        ZSNOWTEMP(:)=ZSNOWTEMP(:)/DGMIK%XTDSNOW(:)
  ELSEWHERE
        ZSNOWTEMP(:)=XUNDEF
  ENDWHERE
  !
!cdir nodep
  DO JJ=1,KSIZE
     JI                     =  KMASK       (JJ)
     !
     DGMIP%AL(KPATCH)%XPSNG  (JI) = IRK%XPSNG     (JJ,1)
     DGMIP%AL(KPATCH)%XPSNV  (JI) = IRK%XPSNV     (JJ,1)
     DGMIP%AL(KPATCH)%XPSN   (JI) = IRK%XPSN      (JJ,1)     
     DGMIP%AL(KPATCH)%XFF    (JI) = INIK%XFF      (JJ,1)
     DGMIP%AL(KPATCH)%XFFG   (JI) = INIK%XFFG     (JJ,1)
     DGMIP%AL(KPATCH)%XFFV   (JI) = INIK%XFFV     (JJ,1)   
     DGMIP%AL(KPATCH)%XFSAT  (JI) = INIK%XFSAT    (JJ)       
     DGMIP%AL(KPATCH)%XHV    (JI) = DGMIK%XHV     (JJ)     
     DGMIP%AL(KPATCH)%XTWSNOW(JI) = DGMIK%XTWSNOW (JJ)
     DGMIP%AL(KPATCH)%XTDSNOW(JI) = DGMIK%XTDSNOW (JJ)
     DGMIP%AL(KPATCH)%XTTSNOW(JI) = ZSNOWTEMP     (JJ)
     !
  END DO
!
  IF (HSNOW=='3-L' .OR. HSNOW=='CRO') THEN
     !
    DO JK=1,SIZE(DGMIK%XSNOWLIQ,2)
!cdir nodep
      DO JJ=1,KSIZE
        JI                      =  KMASK         (JJ)
        !
        DGMIP%AL(KPATCH)%XSNOWLIQ (JI,JK)  =  DGMIK%XSNOWLIQ    (JJ,JK)
        DGMIP%AL(KPATCH)%XSNOWTEMP(JI,JK)  =  DGMIK%XSNOWTEMP   (JJ,JK)
        !
      END DO
    ENDDO
     !
  ENDIF
!
! cosine of solar zenith angle 
!

  IF (IO%CPHOTO/='NON'.AND.IO%LTR_ML) THEN
       !
!cdir nodep
       DO JJ=1,KSIZE
         JI = KMASK(JJ)
         !
         DGMIP%AL(KPATCH)%XFAPAR      (JI) = DGMIK%XFAPAR      (JJ)
         DGMIP%AL(KPATCH)%XFAPIR      (JI) = DGMIK%XFAPIR      (JJ)
         DGMIP%AL(KPATCH)%XFAPAR_BS   (JI) = DGMIK%XFAPAR_BS   (JJ)
         DGMIP%AL(KPATCH)%XFAPIR_BS   (JI) = DGMIK%XFAPIR_BS   (JJ)
         !
       ENDDO
       !
       ! Mask where vegetation evolution is performed (just before solar midnight)
       GMASK = ( PTIME - PTSTEP < 0. ) .AND. ( PTIME >= 0. )
       IF (GMASK) THEN
!cdir nodep
         DO JJ=1,KSIZE
           JI = KMASK(JJ)
           !
           IF (IRK%XMUS(JJ,1).NE.0.) THEN
             DGMIP%AL(KPATCH)%XDFAPARC   (JI) = IRK%XFAPARC   (JJ,1) / IRK%XMUS(JJ,1) 
             DGMIP%AL(KPATCH)%XDFAPIRC   (JI) = IRK%XFAPIRC   (JJ,1) / IRK%XMUS(JJ,1)
             DGMIP%AL(KPATCH)%XDLAI_EFFC (JI) = IRK%XLAI_EFFC (JJ,1) / IRK%XMUS(JJ,1)
           ENDIF
           !
         ENDDO
!cdir nodep         
         DO JJ=1,KSIZE   
           IRK%XFAPARC(JJ,1)   = 0.
           IRK%XFAPIRC(JJ,1)   = 0.
           IRK%XLAI_EFFC(JJ,1) = 0.
           IRK%XMUS(JJ,1)      = 0.
         ENDDO
       ENDIF
       !
  ENDIF
  !
  IF(IO%CISBA=='DIF')THEN
    ZALT(:)=0.0
    ZFLT(:)=0.0
    CALL COMPUT_COLD_LAYERS_THICK(IMXK%XDG(:,:,1),IRK%XTG(:,:,1),ZALT,ZFLT)
    DO JJ=1,KSIZE
       JI              =  KMASK(JJ)
       DGMIP%AL(KPATCH)%XALT(JI) =  ZALT(JJ) 
       DGMIP%AL(KPATCH)%XFLT(JI) =  ZFLT(JJ)  
    ENDDO
  ENDIF
  !
END IF
!
IF (OAGRIP) THEN
  !
!cdir nodep
  DO JJ=1,KSIZE
     JI                     =  KMASK         (JJ)
     !
     DGMIP%AL(KPATCH)%XSEUIL   (JI)  =  AGK%XTHRESHOLDSPT (JJ,1)
     !
  END DO
!
END IF
IF (LHOOK) CALL DR_HOOK('DIAG_MISC_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_MISC_ISBA_n
