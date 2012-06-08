SUBROUTINE EXTRAPOL_FIELDS(HPROGRAM,KLUOUT)
!!
!!    PURPOSE
!!    -------
!!  parameters defined by cover need to be extrapolated if LDATA_VEGTYPE and NOT LDATA_"PARAM"
!!  all ten-day periods are calculated one time for all, then written in PGD.txt
!!
!!    METHOD
!!    ------ 
!!  these parameters are: LAI, HT, DG, ROOTFRAC, IRRIG, WATSUP
!!  Parameters are calculated as in ecoclimap, by vegtype, and then extrapolated
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    S. Faroux        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    16/11/10
!!
!!    DECLARATIONS
!!
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODD_ISBA_n,         ONLY : XCOVER, CISBA, CPHOTO 
USE MODD_DATA_COVER,     ONLY : XDATA_LAI, XDATA_H_TREE,                                &
                                XDATA_IRRIG, XDATA_WATSUP,                              &
                                XDATA_GARDEN, XDATA_NATURE,                             &
                                XDATA_ROOT_DEPTH, XDATA_GROUND_DEPTH,                   &
                                XDATA_ROOT_EXTINCTION, XDATA_ROOT_LIN
!                                
USE MODD_DATA_ISBA_n,    ONLY : XPAR_LAI, XPAR_H_TREE, XPAR_ROOT_DEPTH,           &
                                XPAR_GROUND_DEPTH, XPAR_IRRIG, XPAR_WATSUP,       &
                                LDATA_VEGTYPE, LDATA_LAI, LDATA_H_TREE, LDATA_DG, &
                                LDATA_IRRIG, LDATA_WATSUP, LDATA_ROOTFRAC,        &
                                LDATA_GROUND_DEPTH, LDATA_ROOT_DEPTH, LDATA_Z0
!                                
USE MODD_ISBA_n,         ONLY : CISBA, NGROUND_LAYER
!
USE MODI_AV_PGD
USE MODI_INI_VAR_FROM_VEGTYPE_DATA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=6),       INTENT(IN)    :: HPROGRAM  ! host program
INTEGER,                INTENT(IN)    :: KLUOUT
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
CHARACTER(LEN=3)  :: YTREE, YNAT, YVEG, YDIF
INTEGER :: JTIME
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EXTRAPOL_FIELDS',0,ZHOOK_HANDLE)
!
YNAT ='NAT'
YTREE='TRE'
YVEG ='VEG'
YDIF ='DVG'
!
!            2. Extrapolations for land use or user
!            --------------------------------------
!
!   LAI
!   ---
IF (.NOT.LDATA_LAI) THEN
!
  DO JTIME=1,36 
!    
!   ECOCLIMAP spatial distribution field
    CALL AV_PGD(XPAR_LAI(:,JTIME,:),XCOVER,XDATA_LAI(:,JTIME,:),YVEG,'ARI',KDECADE=JTIME)
!
!   Extrapolation toward new vegtype distribution field from updated land-use map or user 
    CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,KLUOUT,'LAI: leaf area index',XPAR_LAI(:,JTIME,:))
!    
  ENDDO
  LDATA_LAI=.TRUE.
!  
ENDIF
!
!   H_TREE
!   ------
IF (.NOT.LDATA_H_TREE .AND. (CPHOTO/='NON' .OR. .NOT.LDATA_Z0)) THEN
!  
! ECOCLIMAP spatial distribution field       
  CALL AV_PGD(XPAR_H_TREE,XCOVER,XDATA_H_TREE,YTREE,'ARI')
!
! Extrapolation toward new vegtype distribution field from updated land-use map or user  
  CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,KLUOUT,'H_TREE: height of trees',XPAR_H_TREE)
!  
  LDATA_H_TREE=.TRUE.
!  
ENDIF
!
!   DG
!   --
!
!ROOT_DEPTH is needed for DIF, 2-L, 3-L 
IF (.NOT.LDATA_DG .AND. (CISBA/='DIF' .OR. LDATA_ROOTFRAC) .AND. .NOT.LDATA_ROOT_DEPTH) THEN
  CALL AV_PGD (XPAR_ROOT_DEPTH(:,:),XCOVER,XDATA_ROOT_DEPTH(:,:),YNAT,'ARI')
  CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,KLUOUT,'ROOTDEPTH', XPAR_ROOT_DEPTH(:,:))
  LDATA_ROOT_DEPTH = .TRUE.
ENDIF
!
!GROUND_DEPTH is needed for DIF and 3-L
IF (.NOT.LDATA_DG .AND. CISBA/='2-L' .AND. .NOT.LDATA_GROUND_DEPTH) THEN
  CALL AV_PGD (XPAR_GROUND_DEPTH(:,:),XCOVER,XDATA_GROUND_DEPTH(:,:),YNAT,'ARI')
  CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,KLUOUT,'GROUNDDEPTH', XPAR_GROUND_DEPTH(:,:))
  LDATA_GROUND_DEPTH = .TRUE.
ENDIF
!
!
!  IRRIG
!  -----
IF (.NOT.LDATA_IRRIG) THEN
  DO JTIME=1,36
!   ECOCLIMAP spatial distribution field       
    CALL AV_PGD(XPAR_IRRIG(:,JTIME,:),XCOVER,XDATA_IRRIG,YVEG,'ARI',KDECADE=JTIME)
!   Extrapolation toward new vegtype distribution field from updated land-use map or user  
    CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,KLUOUT,'IRRIG  ', XPAR_IRRIG(:,JTIME,:))
  ENDDO
  LDATA_IRRIG=.TRUE.
ENDIF
!
!   WATSUP
!   ------
IF (.NOT.LDATA_WATSUP) THEN
  DO JTIME=1,36
!   ECOCLIMAP spatial distribution field       
    CALL AV_PGD(XPAR_WATSUP(:,JTIME,:),XCOVER,XDATA_WATSUP,YVEG,'ARI',KDECADE=JTIME)  
!   Extrapolation toward new vegtype distribution field from updated land-use map or user  
    CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,KLUOUT,'WATSUP  ', XPAR_WATSUP(:,JTIME,:))
  ENDDO
  LDATA_WATSUP=.TRUE.
ENDIF
!
IF (LHOOK) CALL DR_HOOK('EXTRAPOL_FIELDS',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE EXTRAPOL_FIELDS
