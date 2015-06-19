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
USE MODD_BEM_n, ONLY : B => BEM
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
USE MODD_DIAG_CUMUL_TEB_n, ONLY : DGCT => DIAG_CUMUL_TEB
USE MODD_DIAG_MISC_TEB_n, ONLY : DGMT => DIAG_MISC_TEB
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
USE MODD_DIAG_UTCI_TEB_n, ONLY : DGUT => DIAG_UTCI_TEB
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_PANEL_n, ONLY : TPN => TEB_PANEL
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
!
USE MODD_CH_SEAFLUX_n, ONLY : CHS => CH_SEAFLUX
USE MODD_DIAG_OCEAN_n, ONLY : DGO => DIAG_OCEAN
USE MODD_DIAG_SEAFLUX_n, ONLY : DGS => DIAG_SEAFLUX
USE MODD_DIAG_SEAICE_n, ONLY : DGSI => DIAG_SEAICE
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
!
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DGEI => DIAG_EVAP_ISBA
USE MODD_DIAG_ISBA_n, ONLY : DGI => DIAG_ISBA
USE MODD_DIAG_MISC_ISBA_n, ONLY : DGMI => DIAG_MISC_ISBA
USE MODD_DST_n, ONLY : DST => DST
USE MODD_GR_BIOG_n, ONLY : GB => GR_BIOG
USE MODD_ISBA_n, ONLY : I => ISBA
!
USE MODD_CH_FLAKE_n, ONLY : CHF => CH_FLAKE
USE MODD_CH_WATFLUX_n, ONLY : CHW => CH_WATFLUX
USE MODD_DIAG_FLAKE_n, ONLY : DGF => DIAG_FLAKE
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DGMF => DIAG_MISC_FLAKE
USE MODD_DIAG_WATFLUX_n, ONLY : DGW => DIAG_WATFLUX
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
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
IF (U%NDIM_SEA    >0) CALL WRITE_DIAG_SEA_n(CHS, DGO, DGS, DGSI, DGU, S, U, &
                                                     HPROGRAM,HWRITE)
IF (U%NDIM_WATER  >0) CALL WRITE_DIAG_INLAND_WATER_n(CHF, CHW, DGF, DGMF, DGU, DGW, F, U, W, &
                                                     HPROGRAM,HWRITE)
IF (U%NDIM_NATURE >0) CALL WRITE_DIAG_NATURE_n(CHI, DGEI, DGI, DGMI, DGU, DST, GB, I, U, &
                                                     HPROGRAM,HWRITE)
IF (U%NDIM_TOWN   >0) CALL WRITE_DIAG_TOWN_n(B, BOP, CHT, DGCT, DGMT, DGMTO, DGT, DGUT, U, TGDPE, TGDP, T, TOP, TPN, TVG, &
                                                     HPROGRAM,HWRITE)
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
