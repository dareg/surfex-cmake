!     #########
      SUBROUTINE PUT_SFX_LAND (I, U, &
                               KLUOUT,OCPL_WTD,OCPL_FLOOD, &
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
!!      B. Decharme      *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2009
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
!
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
!
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
INTEGER,           INTENT(IN)  :: KLUOUT
LOGICAL,           INTENT(IN)  :: OCPL_WTD
LOGICAL,           INTENT(IN)  :: OCPL_FLOOD
!
REAL, DIMENSION(:), INTENT(IN) :: PWTD     ! water table depth (negative below soil surface) (m)
REAL, DIMENSION(:), INTENT(IN) :: PFWTD    ! fraction of water table rise (-)
REAL, DIMENSION(:), INTENT(IN) :: PFFLOOD  ! fraction of flooded area (-)
REAL, DIMENSION(:), INTENT(IN) :: PPIFLOOD !(kg/m2)
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
IF(U%NSIZE_NATURE==0)THEN
  IF (LHOOK) CALL DR_HOOK('PUT_SFX_LAND',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
!*       2.0   Put variable over nature
!              ------------------------
!
IF(OCPL_WTD)THEN
!    
  I%IP%XWTD    (:) = XUNDEF
  I%IP%XFWTD   (:) = XUNDEF
!
  YCOMMENT='water table depth'
  CALL PACK_SAME_RANK(U%NR_NATURE(:),PWTD(:),I%IP%XWTD(:))
  CALL CHECK_LAND(YCOMMENT,I%IP%XWTD)
!  
  YCOMMENT='fraction of water table rise'
  CALL PACK_SAME_RANK(U%NR_NATURE(:),PFWTD(:),I%IP%XFWTD(:))
  CALL CHECK_LAND(YCOMMENT,I%IP%XFWTD)
!
  WHERE(I%P%XGW(:)==0.0)
        I%IP%XWTD    (:) = XUNDEF
        I%IP%XFWTD   (:) = 0.0
  ELSEWHERE
        I%IP%XFWTD(:)=I%IP%XFWTD(:)*I%P%XGW(:)
  ENDWHERE
!
ENDIF
!
IF(OCPL_FLOOD)THEN
!
  I%I%XFFLOOD (:) = XUNDEF
  I%I%XPIFLOOD(:) = XUNDEF
!
  YCOMMENT='Flood fraction'
  CALL PACK_SAME_RANK(U%NR_NATURE(:),PFFLOOD(:),I%I%XFFLOOD(:))
  CALL CHECK_LAND(YCOMMENT,I%I%XFFLOOD)
!  
  YCOMMENT='Potential flood infiltration'
  CALL PACK_SAME_RANK(U%NR_NATURE(:),PPIFLOOD(:),I%I%XPIFLOOD(:))
  CALL CHECK_LAND(YCOMMENT,I%I%XPIFLOOD)
!
! No flood for very smal flooded area (<0.1% of grid-cell)
!
  WHERE(I%I%XFFLOOD (:)<0.001)
        I%I%XFFLOOD (:)=0.0
        I%I%XPIFLOOD(:)=0.0
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
