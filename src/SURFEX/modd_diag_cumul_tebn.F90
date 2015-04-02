!     ############################
      MODULE MODD_DIAG_CUMUL_TEB_n
!     ############################
!
!!****  *MODD_DIAG_CUMUL_TEB - declaration of cumulated surface parameters for TEB scheme
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      C de Munck   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       19/02/2013
!
!
!*       0.   DECLARATIONS
!             ------------
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

!TYPE DIAG_CUMUL_TEB_OPTIONS_t
!------------------------------------------------------------------------------
!
!  LOGICAL :: LTEB_CUM   ! flag for cumulated terms of teb scheme
!
!------------------------------------------------------------------------------
!END TYPE DIAG_CUMUL_TEB_OPTIONS_t
!
TYPE DIAG_CUMUL_TEB_t
!------------------------------------------------------------------------------
!* miscellaneous variables
!
  REAL, POINTER, DIMENSION(:)   :: XRUNOFFC_TOWN      ! cumulateda ggregated water runoff for town      (kg/m2 town/s)
  REAL, POINTER, DIMENSION(:)   :: XRUNOFFC_GARDEN    ! cumulated water runoff for green areas          (kg/m2 garden/s)
  REAL, POINTER, DIMENSION(:)   :: XRUNOFFC_ROAD      ! cumulated water runoff for roads                (kg/m2 road/s)
  REAL, POINTER, DIMENSION(:)   :: XRUNOFFC_ROOF      ! cumulated aggregated water runoff for roofs     (kg/m2 roof/s)
  REAL, POINTER, DIMENSION(:)   :: XRUNOFFC_STRLROOF  ! cumulated water runoff for structural roofs     (kg/m2 structural roof/s)
  REAL, POINTER, DIMENSION(:)   :: XRUNOFFC_GREENROOF ! cumulated water runoff for greenroof            (kg/m2 greenroof/s)
  REAL, POINTER, DIMENSION(:)   :: XDRAINC_GREENROOF  ! cumulated water vertical drainage for greenroof (kg/m2 greenroof/s)
  REAL, POINTER, DIMENSION(:)   :: XDRAINC_GARDEN     ! cumulated water vertical drainage for gardens   (kg/m2 garden/s)
  REAL, POINTER, DIMENSION(:)   :: XIRRIGC_GREENROOF  ! cumulated water supply from summer irrigation   (kg/m2 greenroof/s)
  REAL, POINTER, DIMENSION(:)   :: XIRRIGC_GARDEN     ! cumulated water supply from summer irrigation   (kg/m2 garden/s)
  REAL, POINTER, DIMENSION(:)   :: XIRRIGC_ROAD       ! cumulated water supply from summer irrigation   (kg/m2 road/s)
  !
  REAL, POINTER, DIMENSION(:)   :: XHVACC_COOL        ! cumulated en. consump. of the cooling system [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XHVACC_HEAT        ! cumulated en. consump. of the heating system [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XTHER_PROD_BLDC    ! cumulated en. product. of thermal      solar panels [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XPHOT_PROD_BLDC    ! cumulated en. product. of photovoltaic solar panels [W m-2(bld)]
!
!------------------------------------------------------------------------------
END TYPE DIAG_CUMUL_TEB_t

!TYPE(DIAG_CUMUL_TEB_OPTIONS_t), ALLOCATABLE, TARGET, SAVE :: DIAG_CUMUL_TEB_OPTIONS_MODEL(:)
!LOGICAL, POINTER :: LTEB_CUM=>NULL()

TYPE(DIAG_CUMUL_TEB_t),         ALLOCATABLE, TARGET, SAVE :: DIAG_CUMUL_TEB_MODEL(:,:)

TYPE(DIAG_CUMUL_TEB_t), POINTER :: DIAG_CUMUL_TEB => NULL()
!$OMP THREADPRIVATE(DIAG_CUMUL_TEB)

CONTAINS
!
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!
SUBROUTINE DIAG_CUMUL_TEB_GOTO_MODEL(KFROM, KTO, LKFROM, KFROM_PATCH, KTO_PATCH)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
INTEGER, INTENT(IN) :: KFROM_PATCH, KTO_PATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DIAG_CUMUL_TEB_N:DIAG_CUMUL_TEB_GOTO_MODEL',0,ZHOOK_HANDLE)

DIAG_CUMUL_TEB => DIAG_CUMUL_TEB_MODEL(KTO,KTO_PATCH)

IF (LHOOK) CALL DR_HOOK('MODD_DIAG_CUMUL_TEB_N:DIAG_CUMUL_TEB_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE DIAG_CUMUL_TEB_GOTO_MODEL
!
!------------------------------------------------------------------------------
!
SUBROUTINE DIAG_CUMUL_TEB_ALLOC(KMODEL,KPATCH)
INTEGER, INTENT(IN) :: KMODEL, KPATCH
INTEGER :: J,JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_CUMUL_TEB_N:DIAG_CUMUL_TEB_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DIAG_CUMUL_TEB_MODEL(KMODEL,KPATCH))
DO J=1,KMODEL
 DO JP=1,KPATCH
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XRUNOFFC_TOWN)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XRUNOFFC_GARDEN)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XRUNOFFC_ROAD)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XRUNOFFC_ROOF)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XRUNOFFC_STRLROOF)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XRUNOFFC_GREENROOF)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XDRAINC_GREENROOF)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XDRAINC_GARDEN)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XIRRIGC_GREENROOF)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XIRRIGC_GARDEN)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XIRRIGC_ROAD)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XHVACC_COOL)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XHVACC_HEAT)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XTHER_PROD_BLDC)
  NULLIFY(DIAG_CUMUL_TEB_MODEL(J,JP)%XPHOT_PROD_BLDC)
 ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_CUMUL_TEB_N:DIAG_CUMUL_TEB_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_CUMUL_TEB_ALLOC
!
!------------------------------------------------------------------------------
!
SUBROUTINE DIAG_CUMUL_TEB_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_CUMUL_TEB_N:DIAG_CUMUL_TEB_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DIAG_CUMUL_TEB_MODEL)) DEALLOCATE(DIAG_CUMUL_TEB_MODEL)
IF (ASSOCIATED(DIAG_CUMUL_TEB)) NULLIFY(DIAG_CUMUL_TEB)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_CUMUL_TEB_N:DIAG_CUMUL_TEB_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_CUMUL_TEB_DEALLO
!
!------------------------------------------------------------------------------

END MODULE MODD_DIAG_CUMUL_TEB_n
