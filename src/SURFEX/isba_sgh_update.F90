!     #######################################################################
      SUBROUTINE ISBA_SGH_UPDATE (PMESH_SIZE, IO, IP, INP, INI, IR, PRAIN )
!     #######################################################################
!
!!****  *SGH_UPDATE*  
!!
!!    PURPOSE
!!    -------
!
!     Calculates the evolution of the fraction, mu, of the grid cell
!     reached by the rain, the Topmodel saturated fraction and the diagnostic
!     wetland fraction.
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
!!      07/2011 (B. Decharme) : Add fsat diag for dt92
!!      (B. Decharme) 04/2013 : DIF lateral subsurface drainage
!!
!-------------------------------------------------------------------------------
!
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODD_SGH_PAR,     ONLY : NDIMTAB, XMTOKM, XSTOHR, X001,   &
                             XMUREGP, XMUREGA
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(IN) :: PMESH_SIZE
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: INP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
REAL, DIMENSION(:), INTENT(IN)   :: PRAIN
!                                   PRAIN   = rain rate (kg/m2/s)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PRAIN))          :: ZDIST, ZBETA    ! IO%CRAIN = SGH
!                                        ZDIST  = the cell scale (in km)
!                                        ZBETA  = cell scale dependency parameter
!
REAL, DIMENSION(SIZE(PRAIN))          :: ZD_TOP, ZW_TOP, ZQTOP  ! IO%CRUNOFF = SGH
!                                        ZW_TOP = ative TOPMODEL-soil moisture at 't' (m3 m-3)
!                                        ZD_TOP = Topmodel active layer
!                                        ZQTOP  = Topmodel lateral sub-surface flow (-)
!
INTEGER, DIMENSION(SIZE(PRAIN))       :: IUP,IDOWN  ! IO%CRUNOFF = SGH
!                                        change in xsat (or fsat) index
!
INTEGER, DIMENSION(SIZE(PRAIN))       :: NMASK      ! indices correspondance between arrays
!
REAL, DIMENSION(SIZE(PRAIN))          :: ZWSAT_AVG, ZWWILT_AVG
!                                        Average soil properties content
!
REAL                                  :: ZW_UP, ZW_DOWN
REAL                                  :: ZF_UP, ZF_DOWN, ZSLOPEF
REAL                                  :: ZQ_UP, ZQ_DOWN, ZSLOPEQ
!
INTEGER                               :: INJ, JJ, JI, JPATCH, JTAB, ICOUNT, &
                                         JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_SGH_UPDATE',0,ZHOOK_HANDLE)
!
INJ=SIZE(PRAIN,1)
INI%XFSAT(:) = 0.0
!
!*   1.0 Spatial distribution of precipitation
!    ---------------------------------------------
!
IF(IO%CRAIN=='SGH')THEN
!
  WHERE(PRAIN(:)>0.0)
    INI%XMUF (:) =1.0
  ELSEWHERE
    INI%XMUF (:) =0.0
  ENDWHERE

!        
! calculate the cell scale (in km)
!
  ZDIST(:) = SQRT(PMESH_SIZE(:))/XMTOKM
!
  WHERE(ZDIST(:)>=15.0)
!
!       calculate beta for the mu calculation
!         
    ZBETA (:) = XMUREGA + XMUREGP * EXP(-X001*ZDIST(:))
!
!       calculate mu, precip is in mm/hr
!         
    INI%XMUF (:) = 1.0 - EXP(-ZBETA(:)*(PRAIN(:)*XSTOHR))
!
  ENDWHERE
!
ENDIF
!
!*   2.0 Computation of the saturated fraction given by TOPMODEL 
!    -----------------------------------------------------------
!
IF(IO%CRUNOFF=='SGH')THEN
!
! Calculation of the ative TOPMODEL-soil moisture at 't' (m)
! ---------------------------------------------------------------
!
  ZQTOP     (:) = 0.0
  ZW_TOP    (:) = 0.0
  ZD_TOP    (:) = 0.0
  ZWSAT_AVG (:) = 0.0
  ZWWILT_AVG(:) = 0.0
!
  IF(IO%CISBA=='DIF')THEN        
!
    DO JPATCH=1,IO%NPATCH
      IF (INP%NSIZE_NATURE_P(JPATCH)>0 )THEN
      DO JL=1,IO%NLAYER_DUN
        DO JJ=1,INJ
          ZD_TOP    (JJ) = ZD_TOP    (JJ) + INP%XPATCH(JJ,JPATCH)*INP%XSOILWGHT(JJ,JL,JPATCH)
          ZWSAT_AVG (JJ) = ZWSAT_AVG (JJ) + INP%XPATCH(JJ,JPATCH)*INP%XSOILWGHT(JJ,JL,JPATCH)*INP%XWSAT(JJ,JL)
          ZWWILT_AVG(JJ) = ZWWILT_AVG(JJ) + INP%XPATCH(JJ,JPATCH)*INP%XSOILWGHT(JJ,JL,JPATCH)*INP%XWD0 (JJ,JL)
          ZW_TOP    (JJ) = ZW_TOP    (JJ) + INP%XPATCH(JJ,JPATCH)*INP%XSOILWGHT(JJ,JL,JPATCH)*IR%XWG(JJ,JL,JPATCH)
        ENDDO
      ENDDO
      ENDIF
    ENDDO
!
    WHERE(ZD_TOP(:)>0.0)
         ZWSAT_AVG (:) = ZWSAT_AVG (:)/ZD_TOP(:)
         ZWWILT_AVG(:) = ZWWILT_AVG(:)/ZD_TOP(:)
         ZW_TOP    (:) = ZW_TOP    (:)/ZD_TOP(:)
    ENDWHERE
!
  ELSE
!    
    DO JPATCH=1,IO%NPATCH
      IF (INP%NSIZE_NATURE_P(JPATCH)>0 )THEN
        DO JJ=1,INJ
          ZD_TOP(JJ) = ZD_TOP(JJ)+INP%XRUNOFFD(JJ,JPATCH)*INP%XPATCH(JJ,JPATCH)
          ZW_TOP(JJ) = ZW_TOP(JJ)+INP%XRUNOFFD(JJ,JPATCH)*INP%XPATCH(JJ,JPATCH)*IR%XWG(JJ,2,JPATCH)
        ENDDO
      ENDIF
    ENDDO
!  
    WHERE(ZD_TOP(:)>0.0)
          ZW_TOP(:) = ZW_TOP(:) / ZD_TOP(:)
    ENDWHERE
!      
    ZWSAT_AVG (:) = INP%XWSAT(:,1)
    ZWWILT_AVG(:) = INP%XWD0 (:,1)
!
  ENDIF
!
! Find the boundary
! -----------------
!
  NMASK(:)=0
  ICOUNT=0
  DO JJ=1,INJ  
     IF((IP%XTI_MEAN(JJ)/=XUNDEF.AND.ZW_TOP(JJ)<ZWSAT_AVG(JJ).AND.ZW_TOP(JJ)>ZWWILT_AVG(JJ)))THEN     
       ICOUNT=ICOUNT+1
       NMASK(ICOUNT)=JJ       
     ENDIF
     IF(ZW_TOP(JJ)>=ZWSAT_AVG(JJ))THEN
        INI%XFSAT (JJ) = 1.0
     ENDIF
  ENDDO
!     
! compare wt_array and WT
! -----------------------
!
  DO JTAB=1,NDIMTAB
     DO JJ=1,ICOUNT
        JI = NMASK(JJ)    
        IF(INI%XTAB_WTOP(JI,JTAB)>ZW_TOP(JI))THEN
          IUP(JJ)=JTAB
          IDOWN(JJ)=JTAB+1
        ELSEIF(INI%XTAB_WTOP(JI,JTAB)==ZW_TOP(JI))THEN
          IUP(JJ)=JTAB
          IDOWN(JJ)=JTAB
        ENDIF
     ENDDO    
  ENDDO 
!    
! calculate fsat
! --------------
!     
  DO JJ=1,ICOUNT
!  
     JI = NMASK(JJ)
!     
!    new range
     ZF_UP   = INI%XTAB_FSAT(JI,IUP  (JJ))
     ZF_DOWN = INI%XTAB_FSAT(JI,IDOWN(JJ))
     ZQ_UP   = INI%XTAB_QTOP(JI,IUP  (JJ))
     ZQ_DOWN = INI%XTAB_QTOP(JI,IDOWN(JJ))     
     ZW_UP   = INI%XTAB_WTOP(JI,IUP  (JJ))
     ZW_DOWN = INI%XTAB_WTOP(JI,IDOWN(JJ))
!     
!    Calculate new FSAT
     ZSLOPEF = 0.0
     ZSLOPEQ = 0.0
     IF(IUP(JJ)/=IDOWN(JJ))THEN
       ZSLOPEF = (ZF_UP-ZF_DOWN)/(ZW_UP-ZW_DOWN)
       ZSLOPEQ = (ZQ_UP-ZQ_DOWN)/(ZW_UP-ZW_DOWN)
     ENDIF
!     
     INI%XFSAT(JI) = ZF_DOWN+(ZW_TOP(JI)-ZW_DOWN)*ZSLOPEF
     ZQTOP(JI) = ZQ_DOWN+(ZW_TOP(JI)-ZW_DOWN)*ZSLOPEQ
!     
  ENDDO
!
! Subsurface flow by layer (m/s)
! ------------------------------
!
  IF(IO%CISBA=='DIF')THEN        
!
    DO JPATCH=1,IO%NPATCH
      IF(INP%NSIZE_NATURE_P(JPATCH)>0)THEN
        DO JL=1,IO%NLAYER_DUN
           DO JJ=1,INJ
            INI%XTOPQS(JJ,JL,JPATCH)=INP%XKANISO(JJ,JL)*INP%XCONDSAT(JJ,1,JPATCH)*ZQTOP(JJ)*&
                                     INP%XSOILWGHT(JJ,JL,JPATCH)/INP%XRUNOFFD(JJ,JPATCH)
           ENDDO
        ENDDO
      ENDIF
    ENDDO
!
  ENDIF
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ISBA_SGH_UPDATE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------
!
END SUBROUTINE ISBA_SGH_UPDATE
