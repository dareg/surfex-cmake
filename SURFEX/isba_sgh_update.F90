!     ################################################################
      SUBROUTINE ISBA_SGH_UPDATE(HISBA,HRUNOFF,HRAIN,PRAIN,PMUF,PFSAT)
!     ################################################################
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
!!	B. Decharme           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_n,      ONLY : NGROUND_LAYER, NPATCH, XPATCH, XWG, XWWILT, &
                             XWSAT, XDG, XTAB_FSAT, XTAB_WTOP, XRUNOFFB, &
                             XTI_MEAN
!
USE MODD_ISBA_GRID_n, ONLY : XMESH_SIZE
!
USE MODD_SGH_PAR,     ONLY : NDIMTAB, XMTOKM, XSTOHR, X001,      &
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
CHARACTER(LEN=*), INTENT(IN)     :: HISBA  ! type of ISBA version:
!                                          ! '2-L' (default)
!                                          ! '3-L'
!                                          ! 'DIF'
!
CHARACTER(LEN=*), INTENT(IN)     :: HRUNOFF! surface runoff formulation
!                                          ! 'WSAT'
!                                          ! 'DT92'
!                                          ! 'SGH ' Topmodel
!                                                     
!
CHARACTER(LEN=*), INTENT(IN)     :: HRAIN  ! Rainfall spatial distribution
                                           ! 'DEF' = No rainfall spatial distribution
                                           ! 'SGH' = Rainfall exponential spatial distribution
! 
REAL, DIMENSION(:), INTENT(IN)   :: PRAIN
!                                   PRAIN   = rain rate (kg/m2/s)
!
REAL, DIMENSION(:), INTENT(OUT)  :: PMUF
!                                   PMUF = fraction of the grid cell reached by the precipitation
!
REAL, DIMENSION(:), INTENT(OUT)  :: PFSAT
!                                   PFSAT   = Topmodel satured fraction
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PRAIN))          :: ZDIST, ZBETA    ! HRAIN = SGH
!                                        ZDIST  = the cell scale (in km)
!                                        ZBETA  = cell scale dependency parameter
!
REAL, DIMENSION(SIZE(PRAIN))          :: ZW_TOP, ZD_TOP  ! HRUNOFF = SGH
!                                        ZW_TOP = ative TOPMODEL-soil moisture at 't' (m3 m-3)
!                                        ZD_TOP = Depth of the ative layer (m)
!
REAL, DIMENSION(SIZE(PRAIN))          :: ZW    ! HRUNOFF = SGH
!                                        ZW = used for the calculation of ZW_TOP
!
REAL, DIMENSION(NDIMTAB)              :: ZCOMP  ! HRUNOFF = SGH
!                                        ZCOMP = array for wt comparaison with wt_array
!
INTEGER, DIMENSION(1)                 :: IUP,IDOWN  ! HRUNOFF = SGH
!                                        change in xsat (or fsat) index
!
LOGICAL, DIMENSION(SIZE(PRAIN))       :: LSEARCH   ! HRUNOFF = SGH
!
REAL                                  :: ZF_UP, ZF_DOWN, ZW_UP, ZW_DOWN, ZSLOPE
!
INTEGER                               :: I, J    
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_SGH_UPDATE',0,ZHOOK_HANDLE)
!
!*   1.0 Spatial distribution of precipitation
!    ---------------------------------------------
!
IF(HRAIN=='SGH')THEN
!
  WHERE(PRAIN(:)>0.0)
        PMUF (:) =1.0
  ELSEWHERE
        PMUF (:) =0.0
  ENDWHERE

!        
! calculate the cell scale (in km)
!
  ZDIST(:) = SQRT(XMESH_SIZE(:))/XMTOKM
!
  WHERE(ZDIST(:)>=10.0)
!
!       calculate beta for the mu calculation
!         
        ZBETA (:) = XMUREGA + XMUREGP * EXP(-X001*ZDIST(:))
!
!       calculate mu, precip is in mm/hr
!         
        PMUF (:) = 1.0 - EXP(-ZBETA(:)*(PRAIN(:)*XSTOHR))
!
  ENDWHERE
!
ENDIF
!
!*   2.0 Computation of the saturated fraction given by TOPMODEL
!    -----------------------------------------------------------
!
IF(HRUNOFF=='SGH')THEN
!
  PFSAT (:) =0.0
!
! Calculation of the ative TOPMODEL-soil moisture at 't' (m3 m-3)
! ---------------------------------------------------------------
!
  ZD_TOP(:) = 0.0
  ZW_TOP(:) = 0.0
!
  IF(HISBA=='DIF')THEN
    !
    DO I=1,NPATCH
      ZD_TOP(:) = ZD_TOP(:)+XDG     (:,NGROUND_LAYER-1,I)*XPATCH(:,I)
      !
      ZW(:) = XWG(:,1,I)*XDG(:,1,I)
      DO J=2,NGROUND_LAYER-1
        ZW(:) = ZW(:) + XWG(:,J,I)*(XDG(:,J,I)-XDG(:,J-1,I))
      ENDDO
      ZW(:) = ZW(:)/XDG(:,NGROUND_LAYER-1,I)
      !
      ZW_TOP(:) = ZW_TOP(:)+ZW(:)*XPATCH(:,I)
    ENDDO
    !
  ELSE         
    !    
    DO I=1,NPATCH
      ZD_TOP(:) = ZD_TOP(:)+XDG     (:,2,I)*XPATCH(:,I)
      !
      ZW(:) = XWG(:,2,I)*XDG(:,2,I)
      ZW_TOP(:) = ZW_TOP(:)+ZW(:)*XPATCH(:,I)
    ENDDO
    !    
  ENDIF
!
  ZW_TOP(:) = ZW_TOP(:)/ZD_TOP(:)
!
! Find the boundary
! -----------------
!
  LSEARCH(:)=(XTI_MEAN(:)/=XUNDEF)
!
  WHERE(ZW_TOP(:)>=XWSAT(:,1))
        PFSAT  (:) = 1.0
        LSEARCH(:) = .FALSE.
  ENDWHERE
  WHERE(ZW_TOP(:)<=XWWILT(:,1))
        PFSAT  (:) = 0.0
        LSEARCH(:) = .FALSE.
  ENDWHERE
!          
! calculate fsat
! --------------
!
  DO I=1,SIZE(PRAIN)
!
     IF(.NOT.LSEARCH(I))CYCLE
!
!    compare wt_array and WT
!
     ZCOMP(:)=XTAB_WTOP(I,:)-ZW_TOP(I)
!
!    calculate array index
     IUP  (:) = MINLOC(ZCOMP(:),ZCOMP(:)>=0.0)
     IDOWN(:) = MAXLOC(ZCOMP(:),ZCOMP(:)<=0.0)
!     
!    new range
     ZF_UP   = XTAB_FSAT(I,IUP  (1))
     ZF_DOWN = XTAB_FSAT(I,IDOWN(1))
     ZW_UP   = XTAB_WTOP(I,IUP  (1))
     ZW_DOWN = XTAB_WTOP(I,IDOWN(1))
!     
!    Calculate new FSAT
     ZSLOPE = 0.0
     IF(IUP(1)/=IDOWN(1)) ZSLOPE = (ZF_UP-ZF_DOWN)/(ZW_UP-ZW_DOWN)
!     
     PFSAT(I) = ZF_DOWN+(ZW_TOP(I)-ZW_DOWN)*ZSLOPE
!     
  ENDDO 
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ISBA_SGH_UPDATE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------
!
END SUBROUTINE ISBA_SGH_UPDATE
