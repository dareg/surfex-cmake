!     #########
      SUBROUTINE PUT_SFX_LAND(KLUOUT,OCPL_WTD,OCPL_FLOOD, &
                              PWTD,PFWTD,PFFLOOD,PPIFLOOD )  
!     #####################################################
!
!!****  *PUT_SFX_LAND* - routine to put some land surface variables to surfex
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	B. Decharme      *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2009
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODD_SURF_ATM_n, ONLY : NR_NATURE, NSIZE_NATURE
!
USE MODD_ISBA_n,     ONLY : XFFLOOD, XPIFLOOD,&
                            XGW, XWTD, XFWTD
!
USE MODI_PACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,           INTENT(IN)  :: KLUOUT
LOGICAL,           INTENT(IN)  :: OCPL_WTD
LOGICAL,           INTENT(IN)  :: OCPL_FLOOD
!
REAL, DIMENSION(:), INTENT(IN) :: PWTD
REAL, DIMENSION(:), INTENT(IN) :: PFWTD
REAL, DIMENSION(:), INTENT(IN) :: PFFLOOD
REAL, DIMENSION(:), INTENT(IN) :: PPIFLOOD
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
CHARACTER(LEN=50)     :: YCOMMENT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_LAND',0,ZHOOK_HANDLE)
!
!*       1.0   Initialization
!              --------------
!
IF(NSIZE_NATURE==0)THEN
  IF (LHOOK) CALL DR_HOOK('PUT_SFX_LAND',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
!*       2.0   Put variable over nature
!              ------------------------
!
IF(OCPL_WTD)THEN
!    
  XWTD    (:) = XUNDEF
  XFWTD   (:) = XUNDEF
!
  YCOMMENT='water table depth'
  CALL PACK_SAME_RANK(NR_NATURE(:),PWTD(:),XWTD(:))
  CALL CHECK_LAND(YCOMMENT,XWTD)
!  
  YCOMMENT='fraction of water table rise'
  CALL PACK_SAME_RANK(NR_NATURE(:),PFWTD(:),XFWTD(:))
  CALL CHECK_LAND(YCOMMENT,XFWTD)
!
  XFWTD(:)=XFWTD(:)*XGW(:)
!
ENDIF
!
IF(OCPL_FLOOD)THEN
!
  XFFLOOD (:) = XUNDEF
  XPIFLOOD(:) = XUNDEF
!
  YCOMMENT='Flood fraction'
  CALL PACK_SAME_RANK(NR_NATURE(:),PFFLOOD(:),XFFLOOD(:))
  CALL CHECK_LAND(YCOMMENT,XFFLOOD)
!  
  YCOMMENT='Potential flood infiltration'
  CALL PACK_SAME_RANK(NR_NATURE(:),PPIFLOOD(:),XPIFLOOD(:))
  CALL CHECK_LAND(YCOMMENT,XPIFLOOD)
!
! No flood for very smal flooded area (<0.1% of grid-cell)
!
  WHERE(XFFLOOD (:)<0.001)
        XFFLOOD (:)=0.0
        XPIFLOOD(:)=0.0
  ENDWHERE
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_LAND',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE CHECK_LAND(HCOMMENT,PFIELD)
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
CHARACTER(LEN=*),   INTENT(IN) :: HCOMMENT
REAL, DIMENSION(:), INTENT(IN) :: PFIELD
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_LAND:CHECK_LAND',0,ZHOOK_HANDLE)
!
IF(ANY(PFIELD(:)>=XUNDEF))THEN
  WRITE(KLUOUT,*)'PUT_SFX_LAND: problem after get '//TRIM(HCOMMENT)//' from OASIS'
  WRITE(KLUOUT,*)'PUT_SFX_LAND: some points not defined = ',COUNT(PFIELD(:)>=XUNDEF)
  CALL ABOR1_SFX('PUT_SFX_LAND: problem after get '//TRIM(HCOMMENT)//' from OASIS')          
ENDIF
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_LAND:CHECK_LAND',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHECK_LAND
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PUT_SFX_LAND
