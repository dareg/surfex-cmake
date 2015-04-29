!     #################################################################################
SUBROUTINE WRITE_DIAG_SURF_ATM_n(HPROGRAM,HWRITE)
!     #################################################################################
!
!!****  *WRITE_DIAG_SURF_ATM_n * - Chooses the surface schemes for diagnostics writing
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!------------------------------------------------------------------
!
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
!
USE MODD_CH_SNAP_n, ONLY : CHN => CH_EMIS_SNAP
!
USE MODD_CH_EMIS_FIELD_n, ONLY : CHE => CH_EMIS_FIELD
!
USE MODD_SURF_CONF,      ONLY : CPROGNAME
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_CH_SURF_n, ONLY : CHU => CH_SURF
USE MODD_SV_n, ONLY : SV => SV
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
!
USE MODI_WRITE_DIAG_NATURE_n 
USE MODI_WRITE_DIAG_SEA_n 
USE MODI_WRITE_DIAG_INLAND_WATER_n 
USE MODI_WRITE_DIAG_TOWN_n 
!
USE MODI_WRITE_DIAG_SEB_SURF_ATM_n
!
USE MODI_WRITE_DIAG_CH_AGGR_n
USE MODI_WRITE_DIAG_CH_SNAP_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),   INTENT(IN)  :: HWRITE    ! 'PGD' : only physiographic fields are written
!                                            ! 'ALL' : all fields are written
!
!
!*      0.2    declarations of local variables
!
 CHARACTER(LEN=100) :: YCOMMENT
INTEGER            :: IRESP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SURF_ATM_N',0,ZHOOK_HANDLE)
CPROGNAME = HPROGRAM
!
IF (U%NDIM_SEA    >0) CALL WRITE_DIAG_SEA_n         (HPROGRAM,HWRITE)
IF (U%NDIM_WATER  >0) CALL WRITE_DIAG_INLAND_WATER_n(HPROGRAM,HWRITE)
IF (U%NDIM_NATURE >0) CALL WRITE_DIAG_NATURE_n      (HPROGRAM,HWRITE)
IF (U%NDIM_TOWN   >0) CALL WRITE_DIAG_TOWN_n        (HPROGRAM,HWRITE)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Writing
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
IF (DGU%XDIAG_TSTEP==XUNDEF .OR. ABS(NINT(U%TTIME%TIME/DGU%XDIAG_TSTEP)*DGU%XDIAG_TSTEP-U%TTIME%TIME)<1.E-3 ) THEN
  !
  IF (DGU%LFRAC) THEN
    CALL INIT_IO_SURF_n(HPROGRAM,'FULL  ','SURF  ','WRITE')
    YCOMMENT = '(-)'
    CALL WRITE_SURF(HPROGRAM,'FRAC_SEA   ',U%XSEA,   IRESP,HCOMMENT=YCOMMENT)
    CALL WRITE_SURF(HPROGRAM,'FRAC_NATURE',U%XNATURE,IRESP,HCOMMENT=YCOMMENT)
    CALL WRITE_SURF(HPROGRAM,'FRAC_WATER ',U%XWATER, IRESP,HCOMMENT=YCOMMENT)
    CALL WRITE_SURF(HPROGRAM,'FRAC_TOWN  ',U%XTOWN,  IRESP,HCOMMENT=YCOMMENT)
    CALL END_IO_SURF_n(HPROGRAM)
  END IF
  !
  IF (HWRITE/='PGD'.AND.DGU%LDIAG_GRID) CALL WRITE_DIAG_SEB_SURF_ATM_n(DGU, UG, &
                                                                       HPROGRAM)
  !
  IF (CHU%LCH_EMIS .AND. SV%NBEQ>0 .AND. CHU%LCH_SURF_EMIS) THEN
    IF (CHU%CCH_EMIS=='AGGR') THEN 
      CALL WRITE_DIAG_CH_AGGR_n(CHE, &
                                HPROGRAM)
    ELSE IF (CHU%CCH_EMIS=='SNAP') THEN
      CALL WRITE_DIAG_CH_SNAP_n(CHN, &
                                HPROGRAM)
    END IF
  END IF
  !  
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SURF_ATM_N',1,ZHOOK_HANDLE)
!
!--------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_SURF_ATM_n
