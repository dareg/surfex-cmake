!     #########
      SUBROUTINE DIAG_CPL_ESM_ISBA (I, &
                                    PTSTEP,PCPL_DRAIN,PCPL_RUNOFF,PCPL_EFLOOD, &
                                     PCPL_PFLOOD,PCPL_IFLOOD,PCPL_ICEFLUX         )  
!     #####################################################################
!
!!****  *DIAG_CPL_ESM_ISBA*  
!!
!!    PURPOSE
!!    -------
!         
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      B. Decharme           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
!
USE MODD_ISBA_n, ONLY : ISBA_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(ISBA_t), INTENT(INOUT) :: I
!
REAL, INTENT(IN)                   :: PTSTEP
REAL, DIMENSION(:,:), INTENT(IN)   :: PCPL_DRAIN
REAL, DIMENSION(:,:), INTENT(IN)   :: PCPL_RUNOFF
REAL, DIMENSION(:,:), INTENT(IN)   :: PCPL_EFLOOD
REAL, DIMENSION(:,:), INTENT(IN)   :: PCPL_PFLOOD
REAL, DIMENSION(:,:), INTENT(IN)   :: PCPL_IFLOOD
REAL, DIMENSION(:,:), INTENT(IN)   :: PCPL_ICEFLUX
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PCPL_DRAIN,1),SIZE(PCPL_DRAIN,2)) :: ZCPL_DRAIN
REAL, DIMENSION(SIZE(PCPL_DRAIN,1),SIZE(PCPL_DRAIN,2)) :: ZCPL_RECHARGE
!
REAL, DIMENSION(SIZE(I%IP%XPATCH,1)) :: ZSUMPATCH
REAL, DIMENSION(SIZE(I%IP%XPATCH,1)) :: ZBUDGET
!
INTEGER :: INI,INP
INTEGER :: JI, JPATCH ! tile loop counter
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_CPL_ESM_ISBA',0,ZHOOK_HANDLE)
!
!* Initialization
!  --------------
!
INI=SIZE(I%IP%XPATCH,1)
INP=SIZE(I%IP%XPATCH,2)
!
ZSUMPATCH(:) = 0.0
DO JPATCH=1,INP
  DO JI=1,INI
     ZSUMPATCH(JI) = ZSUMPATCH(JI) + I%IP%XPATCH(JI,JPATCH)
  ENDDO
ENDDO
!
ZCPL_RECHARGE(:,:) = 0.0
!
IF(I%O%CISBA/='DIF')THEN
! prevent small negatives values with ISBA-FR
  ZCPL_DRAIN(:,:)=MAX(0.0,PCPL_DRAIN(:,:))
ELSE
  ZCPL_DRAIN(:,:)=PCPL_DRAIN(:,:)
ENDIF
!
!* groundwater case
!  ----------------
!
IF(I%O%LWTD)THEN
  DO JPATCH=1,INP
    DO JI=1,INI
      IF(I%P%XGW(JI)>0.0.AND.ZSUMPATCH(JI)>0.0)THEN
        ZCPL_RECHARGE(JI,JPATCH) = PCPL_DRAIN(JI,JPATCH)
        ZCPL_DRAIN   (JI,JPATCH) = 0.0
      ENDIF
    ENDDO
  ENDDO
ENDIF
!
!* update ISBA - RRM coupling variable (kg/m2)
!  -------------------------------------------
!
!kg/m²
DO JPATCH=1,INP
  DO JI=1,INI
!  
     IF(ZSUMPATCH(JI)>0.0)THEN
       I%I%XCPL_DRAIN (JI) = I%I%XCPL_DRAIN (JI) + PTSTEP * ZCPL_DRAIN (JI,JPATCH) * I%IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI) 
       I%I%XCPL_RUNOFF(JI) = I%I%XCPL_RUNOFF(JI) + PTSTEP * PCPL_RUNOFF(JI,JPATCH) * I%IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI) 
     ENDIF
!
     IF(I%O%LGLACIER.AND.ZSUMPATCH(JI)>0.0)THEN
        I%I%XCPL_ICEFLUX(JI) = I%I%XCPL_ICEFLUX(JI) + PTSTEP * PCPL_ICEFLUX(JI,JPATCH) * I%IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
     ENDIF
!
     IF(I%O%LWTD.AND.ZSUMPATCH(JI)>0.0)THEN
        I%I%XCPL_RECHARGE(JI) = I%I%XCPL_RECHARGE(JI) + PTSTEP * ZCPL_RECHARGE(JI,JPATCH) * I%IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
     ENDIF
!   
     IF(I%O%LFLOOD.AND.ZSUMPATCH(JI)>0.0)THEN
        I%I%XCPL_EFLOOD  (JI) = I%I%XCPL_EFLOOD  (JI) + PTSTEP * PCPL_EFLOOD  (JI,JPATCH)*I%IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
        I%I%XCPL_PFLOOD  (JI) = I%I%XCPL_PFLOOD  (JI) + PTSTEP * PCPL_PFLOOD  (JI,JPATCH)*I%IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
        I%I%XCPL_IFLOOD  (JI) = I%I%XCPL_IFLOOD  (JI) + PTSTEP * PCPL_IFLOOD  (JI,JPATCH)*I%IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
     ENDIF
!    
  ENDDO
ENDDO
!
!* update ISBA Floodplains variable for mass conservation (kg/m2)
!  --------------------------------------------------------------
!
IF(I%O%LFLOOD)THEN
  ZBUDGET(:)=I%I%XPIFLOOD(:)+(I%I%XCPL_PFLOOD(:)-I%I%XCPL_IFLOOD(:)-I%I%XCPL_EFLOOD(:))
  WHERE(ZBUDGET (:)<=0.0)
        I%I%XPIFLOOD(:)=0.0
        I%I%XFFLOOD (:)=0.0
  ENDWHERE
ENDIF
!
IF (LHOOK) CALL DR_HOOK('DIAG_CPL_ESM_ISBA',1,ZHOOK_HANDLE)
!
END SUBROUTINE DIAG_CPL_ESM_ISBA
