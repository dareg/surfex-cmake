!     #########
      SUBROUTINE INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,ILUOUT,HNAME,PFIELD,PDEF)
!     ##############################################################
!
!!
!!    This routine defines the points where extrapolation is done and calls
!!      interpol_field, for every vegtype
!!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_DATA_ISBA_n, ONLY : XPAR_VEGTYPE
USE MODD_SURF_ATM_n,     ONLY : NSIZE_FULL
!
USE MODI_INTERPOL_FIELD
USE MODI_UNPACK_SAME_RANK
USE MODI_PACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_SURF_MASK_n
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=6), INTENT(IN)    :: HPROGRAM  ! host model
INTEGER, INTENT(IN) :: ILUOUT
CHARACTER(LEN=*), INTENT(IN) :: HNAME
REAL, DIMENSION(:), OPTIONAL, INTENT(IN) :: PDEF 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PFIELD
!

!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL, DIMENSION(:), ALLOCATABLE :: ZFIELD_TOT
INTEGER, DIMENSION(:), ALLOCATABLE :: IMASK  ! mask for packing from complete field to nature field
INTEGER, DIMENSION(:), ALLOCATABLE :: NSIZE, NSIZE_TOT
INTEGER               :: IDIM
INTEGER               :: JVEGTYPE  ! loop counter on vegtypes
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!

!-------------------------------------------------------------------------------
!
!*    1.      Initializations
!             ---------------
!
IF (LHOOK) CALL DR_HOOK('INI_VAR_FROM_VEGTYPE_DATA',0,ZHOOK_HANDLE)
!
IDIM=SIZE(PFIELD,1)
!
ALLOCATE(IMASK(IDIM))
ALLOCATE(NSIZE(IDIM))
ALLOCATE(NSIZE_TOT(NSIZE_FULL))
ALLOCATE(ZFIELD_TOT(NSIZE_FULL))
!
DO JVEGTYPE=1,SIZE(PFIELD,2)
  NSIZE(:)=0
  WHERE (PFIELD(:,JVEGTYPE).NE.XUNDEF) NSIZE(:)=1
  WHERE (XPAR_VEGTYPE(:,JVEGTYPE)==0.) NSIZE(:)=-1
  CALL GET_SURF_MASK_n('NATURE',IDIM,IMASK,NSIZE_FULL,ILUOUT)
  CALL UNPACK_SAME_RANK(IMASK,NSIZE,NSIZE_TOT,-1)
  CALL UNPACK_SAME_RANK(IMASK,PFIELD(:,JVEGTYPE),ZFIELD_TOT)
  CALL INTERPOL_FIELD(HPROGRAM,ILUOUT,NSIZE_TOT,ZFIELD_TOT,HNAME,PDEF(JVEGTYPE))
  CALL PACK_SAME_RANK(IMASK,ZFIELD_TOT,PFIELD(:,JVEGTYPE))  
ENDDO
!
DEALLOCATE(IMASK)
DEALLOCATE(NSIZE)
DEALLOCATE(NSIZE_TOT)
DEALLOCATE(ZFIELD_TOT)
!
IF (LHOOK) CALL DR_HOOK('INI_VAR_FROM_VEGTYPE_DATA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_VAR_FROM_VEGTYPE_DATA
