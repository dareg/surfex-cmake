!     #########
      SUBROUTINE DIAG_CPL_ESM_ISBA (IO, P, IP, INI, PTSTEP, PCPL_DRAIN, PCPL_RUNOFF, &
                                    PCPL_EFLOOD,PCPL_PFLOOD,PCPL_IFLOOD,PCPL_ICEFLUX  )  
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
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PGD_t), INTENT(INOUT) :: P
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
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
REAL, DIMENSION(SIZE(IP%XPATCH,1)) :: ZSUMPATCH
REAL, DIMENSION(SIZE(IP%XPATCH,1)) :: ZBUDGET
!
INTEGER :: INJ,JP
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
INJ=SIZE(IP%XPATCH,1)
JP=SIZE(IP%XPATCH,2)
!
ZSUMPATCH(:) = 0.0
DO JPATCH=1,JP
  DO JI=1,INJ
     ZSUMPATCH(JI) = ZSUMPATCH(JI) + IP%XPATCH(JI,JPATCH)
  ENDDO
ENDDO
!
ZCPL_RECHARGE(:,:) = 0.0
!
IF(IO%CISBA/='DIF')THEN
! prevent small negatives values with ISBA-FR
  ZCPL_DRAIN(:,:)=MAX(0.0,PCPL_DRAIN(:,:))
ELSE
  ZCPL_DRAIN(:,:)=PCPL_DRAIN(:,:)
ENDIF
!
!* groundwater case
!  ----------------
!
IF(IO%LWTD)THEN
  DO JPATCH=1,JP
    DO JI=1,INJ
      IF(P%XGW(JI)>0.0.AND.ZSUMPATCH(JI)>0.0)THEN
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
DO JPATCH=1,JP
  DO JI=1,INJ
!  
     IF(ZSUMPATCH(JI)>0.0)THEN
       INI%XCPL_DRAIN (JI) = INI%XCPL_DRAIN (JI) + PTSTEP * ZCPL_DRAIN (JI,JPATCH) * IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI) 
       INI%XCPL_RUNOFF(JI) = INI%XCPL_RUNOFF(JI) + PTSTEP * PCPL_RUNOFF(JI,JPATCH) * IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI) 
     ENDIF
!
     IF(IO%LGLACIER.AND.ZSUMPATCH(JI)>0.0)THEN
        INI%XCPL_ICEFLUX(JI) = INI%XCPL_ICEFLUX(JI) + PTSTEP * PCPL_ICEFLUX(JI,JPATCH) * IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
     ENDIF
!
     IF(IO%LWTD.AND.ZSUMPATCH(JI)>0.0)THEN
        INI%XCPL_RECHARGE(JI) = INI%XCPL_RECHARGE(JI) + PTSTEP * ZCPL_RECHARGE(JI,JPATCH) * IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
     ENDIF
!   
     IF(IO%LFLOOD.AND.ZSUMPATCH(JI)>0.0)THEN
        INI%XCPL_EFLOOD  (JI) = INI%XCPL_EFLOOD  (JI) + PTSTEP * PCPL_EFLOOD  (JI,JPATCH)*IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
        INI%XCPL_PFLOOD  (JI) = INI%XCPL_PFLOOD  (JI) + PTSTEP * PCPL_PFLOOD  (JI,JPATCH)*IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
        INI%XCPL_IFLOOD  (JI) = INI%XCPL_IFLOOD  (JI) + PTSTEP * PCPL_IFLOOD  (JI,JPATCH)*IP%XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
     ENDIF
!    
  ENDDO
ENDDO
!
!* update ISBA Floodplains variable for mass conservation (kg/m2)
!  --------------------------------------------------------------
!
IF(IO%LFLOOD)THEN
  ZBUDGET(:)=INI%XPIFLOOD(:)+(INI%XCPL_PFLOOD(:)-INI%XCPL_IFLOOD(:)-INI%XCPL_EFLOOD(:))
  WHERE(ZBUDGET (:)<=0.0)
        INI%XPIFLOOD(:)=0.0
        INI%XFFLOOD (:)=0.0
  ENDWHERE
ENDIF
!
IF (LHOOK) CALL DR_HOOK('DIAG_CPL_ESM_ISBA',1,ZHOOK_HANDLE)
!
END SUBROUTINE DIAG_CPL_ESM_ISBA
