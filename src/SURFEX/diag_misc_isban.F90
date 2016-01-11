!     #########
SUBROUTINE DIAG_MISC_ISBA_n (DGMIP, PKD, OSURF_MISC_BUDGET, &
                             PTSTEP, HISBA, HPHOTO, HSNOW, OAGRIP, OTR_ML, &
                            PTIME, KSIZE, KPATCH, KMASK, PSEUIL,          &
                            PPSN, PPSNG, PPSNV, PFF, PFFG, PFFV,          &
                            PWG, PWGI, PWFC, PWWILT, PWSNOW, PRSNOW,      &
                            PFAPARC, PFAPIRC, PLAI_EFFC, PMUS, PFSAT,     &
                            PDG, PTG                                      )  
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
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_PATCH_t
USE MODD_PACK_DIAG_ISBA, ONLY : PACK_DIAG_ISBA_t
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
!
TYPE(DIAG_MISC_ISBA_PATCH_t), INTENT(INOUT) :: DGMIP
TYPE(PACK_DIAG_ISBA_t), INTENT(INOUT) :: PKD
!
LOGICAL, INTENT(IN) :: OSURF_MISC_BUDGET
REAL,               INTENT(IN)    :: PTSTEP        ! timestep for  accumulated values 
 CHARACTER(LEN=*), INTENT(IN)      :: HISBA         ! ISBA scheme
 CHARACTER(LEN=*), INTENT(IN)      :: HPHOTO        ! type of photosynthesis
 CHARACTER(LEN=*), INTENT(IN)      :: HSNOW         ! snow scheme
LOGICAL, INTENT(IN)               :: OAGRIP
LOGICAL, INTENT(IN)               :: OTR_ML
REAL,    INTENT(IN)               :: PTIME   ! current time since midnight
INTEGER, INTENT(IN)               :: KSIZE, KPATCH
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
REAL, DIMENSION(:), INTENT(IN)    :: PSEUIL
!    
!Snow/flood fraction at t
REAL, DIMENSION(:), INTENT(IN)    :: PPSN
REAL, DIMENSION(:), INTENT(IN)    :: PPSNG
REAL, DIMENSION(:), INTENT(IN)    :: PPSNV
REAL, DIMENSION(:), INTENT(IN)    :: PFF
REAL, DIMENSION(:), INTENT(IN)    :: PFFG
REAL, DIMENSION(:), INTENT(IN)    :: PFFV
!
REAL, DIMENSION(:,:),  INTENT(IN) :: PWG           ! soil water content profile (m3/m3)
REAL, DIMENSION(:,:),  INTENT(IN) :: PWGI          ! soil solid water content profile (m3/m3)
REAL, DIMENSION(:,:),  INTENT(IN) :: PWFC          ! field capacity profile (m3/m3)
REAL, DIMENSION(:,:),  INTENT(IN) :: PWWILT        ! wilting point  profile (m3/m3)
REAL, DIMENSION(:,:),  INTENT(IN) :: PWSNOW        ! snow reservoir (kg/m2)
REAL, DIMENSION(:,:),  INTENT(IN) :: PRSNOW        ! snow density (kg/m3)
!
REAL, DIMENSION(:,:),  INTENT(IN) :: PDG           ! soil layer depth
REAL, DIMENSION(:,:),  INTENT(IN) :: PTG           ! soil temperature
!
REAL, DIMENSION(:), INTENT(INOUT) :: PFAPARC
REAL, DIMENSION(:), INTENT(INOUT) :: PFAPIRC
REAL, DIMENSION(:), INTENT(INOUT) :: PLAI_EFFC
REAL, DIMENSION(:), INTENT(INOUT) :: PMUS
!
REAL, DIMENSION(:), INTENT(IN)    :: PFSAT
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PPSN))    :: ZSNOWTEMP
REAL, DIMENSION(SIZE(PWSNOW,1),SIZE(PWSNOW,2)) :: ZWORK
REAL, DIMENSION(SIZE(PWSNOW,1),SIZE(PWSNOW,2)) :: ZWORKTEMP
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
  PKD%DGMI%XSWI (:,:)=XUNDEF
  PKD%DGMI%XTSWI(:,:)=XUNDEF  
  DO JJ=1,SIZE(PWG,2)
    DO JI=1,SIZE(PWG,1)
      IF(PWG (JI,JJ)/=XUNDEF)THEN    
        PKD%DGMI%XSWI (JI,JJ) = (PWG (JI,JJ) - PWWILT(JI,JJ)) / (PWFC(JI,JJ) - PWWILT(JI,JJ))
        PKD%DGMI%XTSWI(JI,JJ) = (PWG (JI,JJ) - PWWILT(JI,JJ)) / (PWFC(JI,JJ) - PWWILT(JI,JJ))
      ENDIF
      IF(PWGI (JI,JJ)/=XUNDEF)THEN    
        PKD%DGMI%XTSWI(JI,JJ) = PKD%DGMI%XTSWI(JI,JJ) +  PWGI(JI,JJ) / (PWFC(JI,JJ) - PWWILT(JI,JJ))
      ENDIF
    ENDDO
  ENDDO
  !
  DO JK=1,SIZE(PKD%DGMI%XSWI,2)
!cdir nodep
    DO JJ=1,KSIZE
      JI                      =  KMASK         (JJ)
      !
      DGMIP%AL(KPATCH)%XSWI     (JI,JK)  =  PKD%DGMI%XSWI        (JJ,JK)
      DGMIP%AL(KPATCH)%XTSWI    (JI,JK)  =  PKD%DGMI%XTSWI       (JJ,JK)
      !
    END DO
  ENDDO  
  !
  DO JI = 1,SIZE(PWSNOW,2)
!cdir nodep 
    DO JJ = 1,SIZE(PWSNOW,1)
      ZWORK(JJ,JI)  = PWSNOW(JJ,JI) / PRSNOW(JJ,JI)
    ENDDO
  ENDDO
  !
  PKD%DGMI%XTWSNOW=0.
  PKD%DGMI%XTDSNOW=0.
  ZSNOWTEMP=0.  
  !
  IF (HSNOW/='EBA')THEN
     ZWORKTEMP(:,:) = PKD%DGMI%XSNOWTEMP(:,:)
  ELSE
     ZWORKTEMP(:,1) = MIN(PTG(:,1),XTT)
  ENDIF
  !
  DO JI = 1,SIZE(PWSNOW,2)
!cdir nodep 
    DO JJ = 1,SIZE(PWSNOW,1)
      PKD%DGMI%XTWSNOW(JJ) = PKD%DGMI%XTWSNOW(JJ) + PWSNOW(JJ,JI)      
      PKD%DGMI%XTDSNOW(JJ) = PKD%DGMI%XTDSNOW(JJ) + ZWORK (JJ,JI)
      ZSNOWTEMP(JJ) = ZSNOWTEMP(JJ) + ZWORKTEMP(JJ,JI) * ZWORK(JJ,JI)
    ENDDO
  ENDDO
  !
  WHERE(PKD%DGMI%XTDSNOW(:)>0.0)
        ZSNOWTEMP(:)=ZSNOWTEMP(:)/PKD%DGMI%XTDSNOW(:)
  ELSEWHERE
        ZSNOWTEMP(:)=XUNDEF
  ENDWHERE
  !
!cdir nodep
  DO JJ=1,KSIZE
     JI                     =  KMASK       (JJ)
     !
     DGMIP%AL(KPATCH)%XHV      (JI)  =  PKD%DGMI%XHV       (JJ)
     DGMIP%AL(KPATCH)%XPSNG   (JI)  =  PPSNG       (JJ)
     DGMIP%AL(KPATCH)%XPSNV   (JI)  =  PPSNV       (JJ)
     DGMIP%AL(KPATCH)%XPSN    (JI)  =  PPSN        (JJ)     
     DGMIP%AL(KPATCH)%XFF     (JI)  =  PFF         (JJ)
     DGMIP%AL(KPATCH)%XFFG    (JI)  =  PFFG        (JJ)
     DGMIP%AL(KPATCH)%XFFV    (JI)  =  PFFV        (JJ)     
     DGMIP%AL(KPATCH)%XTWSNOW  (JI)  =  PKD%DGMI%XTWSNOW   (JJ)
     DGMIP%AL(KPATCH)%XTDSNOW  (JI)  =  PKD%DGMI%XTDSNOW   (JJ)
     DGMIP%AL(KPATCH)%XTTSNOW  (JI)  =  ZSNOWTEMP   (JJ)
     DGMIP%AL(KPATCH)%XFSAT   (JI)  =  PFSAT       (JJ)     
     !
  END DO
!
  IF (HSNOW=='3-L' .OR. HSNOW=='CRO') THEN
     !
    DO JK=1,SIZE(PKD%DGMI%XSNOWLIQ,2)
!cdir nodep
      DO JJ=1,KSIZE
        JI                      =  KMASK         (JJ)
        !
        DGMIP%AL(KPATCH)%XSNOWLIQ (JI,JK)  =  PKD%DGMI%XSNOWLIQ    (JJ,JK)
        DGMIP%AL(KPATCH)%XSNOWTEMP(JI,JK)  =  PKD%DGMI%XSNOWTEMP   (JJ,JK)
        !
      END DO
    ENDDO
     !
  ENDIF
!
! cosine of solar zenith angle 
!

  IF (HPHOTO/='NON'.AND.OTR_ML) THEN
       !
!cdir nodep
       DO JJ=1,KSIZE
         JI = KMASK(JJ)
         !
         DGMIP%AL(KPATCH)%XFAPAR      (JI) = PKD%DGMI%XFAPAR      (JJ)
         DGMIP%AL(KPATCH)%XFAPIR      (JI) = PKD%DGMI%XFAPIR      (JJ)
         DGMIP%AL(KPATCH)%XFAPAR_BS   (JI) = PKD%DGMI%XFAPAR_BS   (JJ)
         DGMIP%AL(KPATCH)%XFAPIR_BS   (JI) = PKD%DGMI%XFAPIR_BS   (JJ)
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
           IF (PMUS(JJ).NE.0.) THEN
             DGMIP%AL(KPATCH)%XDFAPARC   (JI) = PFAPARC   (JJ) / PMUS(JJ) 
             DGMIP%AL(KPATCH)%XDFAPIRC   (JI) = PFAPIRC   (JJ) / PMUS(JJ)
             DGMIP%AL(KPATCH)%XDLAI_EFFC (JI) = PLAI_EFFC (JJ) / PMUS(JJ)
           ENDIF
           !
         ENDDO
!cdir nodep         
         DO JJ=1,KSIZE   
           PFAPARC(JJ)   = 0.
           PFAPIRC(JJ)   = 0.
           PLAI_EFFC(JJ) = 0.
           PMUS(JJ)      = 0.
         ENDDO
       ENDIF
       !
  ENDIF
  !
  IF(HISBA=='DIF')THEN
    ZALT(:)=0.0
    ZFLT(:)=0.0
    CALL COMPUT_COLD_LAYERS_THICK(PDG,PTG,ZALT,ZFLT)
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
     DGMIP%AL(KPATCH)%XSEUIL   (JI)  =  PSEUIL (JJ)
     !
  END DO
!
END IF
IF (LHOOK) CALL DR_HOOK('DIAG_MISC_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_MISC_ISBA_n
