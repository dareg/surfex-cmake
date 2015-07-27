!     #############################################################
      SUBROUTINE INIT_INLAND_WATER_n(HPROGRAM,HINIT,                        &
                                   KI,KSV,KSW,                                &
                                   HSV,PCO2,PRHOA,                            &
                                   PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                                   PEMIS,PTSRAD,PTSURF,                       &
                                   KYEAR, KMONTH,KDAY, PTIME,                 &
                                   HATMFILE,HATMFILETYPE,                     &
                                   HTEST                                      )  
!     #############################################################
!
!!****  *INIT_INLAND_WATER_n* - routine to initialize inland water
!!
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!!      B. Decharme  04/2013 new coupling variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_WATFLUX_GRID_n, ONLY : WG => WATFLUX_GRID
!
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_BLD_DESCRIPTION_n, ONLY : BDD => BLD_DESC
USE MODD_CH_EMIS_FIELD_n, ONLY : CHE => CH_EMIS_FIELD
USE MODD_CH_FLAKE_n, ONLY : CHF => CH_FLAKE
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
USE MODD_CH_SEAFLUX_n, ONLY : CHS => CH_SEAFLUX
USE MODD_CH_SNAP_n, ONLY : CHN => CH_EMIS_SNAP
USE MODD_CH_SURF_n, ONLY : CHU => CH_SURF
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
USE MODD_CH_WATFLUX_n, ONLY : CHW => CH_WATFLUX
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_DATA_SEAFLUX_n, ONLY : DTS => DATA_SEAFLUX
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
USE MODD_DATA_TSZ0_n, ONLY : DTZ => DATA_TSZ0
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DGEI => DIAG_EVAP_ISBA
USE MODD_DIAG_FLAKE_n, ONLY : DGF => DIAG_FLAKE
USE MODD_DIAG_ISBA_n, ONLY : DGI => DIAG_ISBA
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DGMF => DIAG_MISC_FLAKE
USE MODD_DIAG_MISC_ISBA_n, ONLY : DGMI => DIAG_MISC_ISBA
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
USE MODD_DIAG_OCEAN_n, ONLY : DGO => DIAG_OCEAN
USE MODD_DIAG_SEAFLUX_n, ONLY : DGS => DIAG_SEAFLUX
USE MODD_DIAG_SEAICE_n, ONLY : DGSI => DIAG_SEAICE
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
USE MODD_DIAG_UTCI_TEB_n, ONLY : DGUT => DIAG_UTCI_TEB
USE MODD_DIAG_WATFLUX_n, ONLY : DGW => DIAG_WATFLUX
USE MODD_FLAKE_GRID_n, ONLY : FG => FLAKE_GRID
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_FLAKE_SBL_n, ONLY : FSB => FLAKE_SBL
USE MODD_GR_BIOG_n, ONLY : GB => GR_BIOG
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
USE MODD_ISBA_CANOPY_n, ONLY : ICP => ISBA_CANOPY
USE MODD_ISBA_n, ONLY : I => ISBA
USE MODD_OCEAN_n, ONLY : O => OCEAN
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
USE MODD_SEAFLUX_SBL_n, ONLY : SSB => SEAFLUX_SBL
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_SV_n, ONLY : SV => SV
USE MODD_TEB_CANOPY_n, ONLY : TCP => TEB_CANOPY
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GREENROOF_n, ONLY : TGR => TEB_GREENROOF
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
USE MODD_WATFLUX_SBL_n, ONLY : WSB => WATFLUX_SBL
!
USE MODD_DIAG_IDEAL_n, ONLY : DGL => DIAG_IDEAL
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
USE MODD_CSTS,       ONLY : XTT
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_INIT_FLAKE_n
!
USE MODI_INIT_IDEAL_FLUX
!
USE MODI_INIT_WATFLUX_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
 CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),                 INTENT(IN)  :: HINIT     ! choice of fields to initialize
INTEGER,                          INTENT(IN)  :: KI        ! number of points
INTEGER,                          INTENT(IN)  :: KSV       ! number of scalars
INTEGER,                          INTENT(IN)  :: KSW       ! number of short-wave spectral bands
 CHARACTER(LEN=6), DIMENSION(KSV), INTENT(IN)  :: HSV       ! name of all scalar variables
REAL,             DIMENSION(KI),  INTENT(IN)  :: PCO2      ! CO2 concentration (kg/m3)
REAL,             DIMENSION(KI),  INTENT(IN)  :: PRHOA     ! air density
REAL,             DIMENSION(KI),  INTENT(IN)  :: PZENITH   ! solar zenithal angle
REAL,             DIMENSION(KI),  INTENT(IN)  :: PAZIM     ! solar azimuthal angle (rad from N, clock)
REAL,             DIMENSION(KSW), INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),  INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSRAD    ! radiative temperature
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
INTEGER,                          INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,                          INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,                          INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                             INTENT(IN)  :: PTIME     ! current time since
                                                          !  midnight (UTC, s)
!
 CHARACTER(LEN=28),                INTENT(IN)  :: HATMFILE    ! atmospheric file name
 CHARACTER(LEN=6),                 INTENT(IN)  :: HATMFILETYPE! atmospheric file type
 CHARACTER(LEN=2),                 INTENT(IN)  :: HTEST       ! must be equal to 'OK'
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
!*       2.     Selection of surface scheme
!               ---------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_INLAND_WATER_N',0,ZHOOK_HANDLE)
IF (U%CWATER=='NONE  ') THEN
  PDIR_ALB=0.
  PSCA_ALB=0.
  PEMIS   =1.
  PTSRAD  =XTT
  PTSURF  =XTT
ELSE IF (U%CWATER=='FLUX  ') THEN
  CALL INIT_IDEAL_FLUX(DGL, &
                       HPROGRAM,HINIT,KI,KSV,KSW,HSV,PCO2,PRHOA,   &
                         PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB,  &
                         PEMIS,PTSRAD,PTSURF,'OK'                    )  
ELSE IF (U%CWATER=='WATFLX') THEN
  CALL INIT_WATFLUX_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                 DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, &
                                 DGT, DGUT, DGW, F, FSB, GB, IOB, ICP, I, O, S, &
                                 SSB, UG, U, SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, &
                                 TVG, WG, W, WSB, &
                      HPROGRAM,HINIT,KI,KSV,KSW,HSV,PCO2,PRHOA,     &
                        PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB,    &
                        PEMIS,PTSRAD,PTSURF,                          &
                        KYEAR,KMONTH,KDAY,PTIME,HATMFILE,HATMFILETYPE,&
                        'OK'                                          )  
ELSE IF (U%CWATER=='FLAKE ') THEN
  CALL INIT_FLAKE_n(BOP, BDD, CHE, CHF, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, &
                         DTS, DTT, DTZ, DGEI, DGF, DGI, DGMF, DGMI, DGMTO, DGO, DGS, &
                         DGSI, DGU, DGT, DGUT, DGW, FG, F, FSB, GB, IOB, ICP, &
                         I, O, S, SSB, UG, U, SV, TCP, TGD, TGDO, TGR, &
                         TGRO, T, TOP, TVG, W, WSB, &
                    HPROGRAM,HINIT,KI,KSV,KSW,HSV,PCO2,PRHOA,       &
                        PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB,    &
                        PEMIS,PTSRAD,PTSURF,                          &
                        KYEAR,KMONTH,KDAY,PTIME,HATMFILE,HATMFILETYPE,&
                        'OK')          
END IF
IF (LHOOK) CALL DR_HOOK('INIT_INLAND_WATER_N',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_INLAND_WATER_n
