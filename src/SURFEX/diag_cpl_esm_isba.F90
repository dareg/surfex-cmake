!     #########
      SUBROUTINE DIAG_CPL_ESM_ISBA(PTSTEP,PCPL_DRAIN,PCPL_RUNOFF,PCPL_EFLOOD, &
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
!!	B. Decharme           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,        ONLY : XLVTT, XRHOLW
USE MODD_ISBA_GRID_n, ONLY : XMESH_SIZE
USE MODD_ISBA_n,      ONLY : XPATCH, LFLOOD, LGLACIER, &
                             CISBA, LWTD, XGW,         &
                             XFFLOOD, XPIFLOOD,        &
                             XCPL_DRAIN, XCPL_RUNOFF,  &
                             XCPL_RECHARGE,            &
                             XCPL_ICEFLUX, XCPL_EFLOOD,&
                             XCPL_PFLOOD, XCPL_IFLOOD
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
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
REAL, DIMENSION(SIZE(XPATCH,1)) :: ZSUMPATCH
REAL, DIMENSION(SIZE(XPATCH,1)) :: ZBUDGET
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
INI=SIZE(XPATCH,1)
INP=SIZE(XPATCH,2)
!
ZSUMPATCH(:) = 0.0
DO JPATCH=1,INP
  DO JI=1,INI
     ZSUMPATCH(JI) = ZSUMPATCH(JI) + XPATCH(JI,JPATCH)
  ENDDO
ENDDO
!
ZCPL_RECHARGE(:,:) = 0.0
!
IF(CISBA/='DIF')THEN
! prevent small negatives values with ISBA-FR
  ZCPL_DRAIN(:,:)=MAX(0.0,PCPL_DRAIN(:,:))
ELSE
  ZCPL_DRAIN(:,:)=PCPL_DRAIN(:,:)
ENDIF
!
!* groundwater case
!  ----------------
!
IF(LWTD)THEN
  DO JPATCH=1,INP
    DO JI=1,INI
      IF(XGW(JI)>0.0.AND.ZSUMPATCH(JI)>0.0)THEN
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
       XCPL_DRAIN (JI) = XCPL_DRAIN (JI) + PTSTEP * ZCPL_DRAIN (JI,JPATCH) * XPATCH(JI,JPATCH)/ZSUMPATCH(JI) 
       XCPL_RUNOFF(JI) = XCPL_RUNOFF(JI) + PTSTEP * PCPL_RUNOFF(JI,JPATCH) * XPATCH(JI,JPATCH)/ZSUMPATCH(JI) 
     ENDIF
!
     IF(LGLACIER.AND.ZSUMPATCH(JI)>0.0)THEN
        XCPL_ICEFLUX(JI) = XCPL_ICEFLUX(JI) + PTSTEP * PCPL_ICEFLUX(JI,JPATCH) * XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
     ENDIF
!
     IF(LWTD.AND.ZSUMPATCH(JI)>0.0)THEN
        XCPL_RECHARGE(JI) = XCPL_RECHARGE(JI) + PTSTEP * ZCPL_RECHARGE(JI,JPATCH) * XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
     ENDIF
!   
     IF(LFLOOD.AND.ZSUMPATCH(JI)>0.0)THEN
        XCPL_EFLOOD  (JI) = XCPL_EFLOOD  (JI) + PTSTEP * PCPL_EFLOOD  (JI,JPATCH)*XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
        XCPL_PFLOOD  (JI) = XCPL_PFLOOD  (JI) + PTSTEP * PCPL_PFLOOD  (JI,JPATCH)*XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
        XCPL_IFLOOD  (JI) = XCPL_IFLOOD  (JI) + PTSTEP * PCPL_IFLOOD  (JI,JPATCH)*XPATCH(JI,JPATCH)/ZSUMPATCH(JI)
     ENDIF
!    
  ENDDO
ENDDO
!
!* update ISBA Floodplains variable for mass conservation (kg/m2)
!  --------------------------------------------------------------
!
IF(LFLOOD)THEN
  ZBUDGET(:)=XPIFLOOD(:)+(XCPL_PFLOOD(:)-XCPL_IFLOOD(:)-XCPL_EFLOOD(:))
  WHERE(ZBUDGET (:)<=0.0)
        XPIFLOOD(:)=0.0
        XFFLOOD (:)=0.0
  ENDWHERE
ENDIF
!
IF (LHOOK) CALL DR_HOOK('DIAG_CPL_ESM_ISBA',1,ZHOOK_HANDLE)
!
END SUBROUTINE DIAG_CPL_ESM_ISBA
