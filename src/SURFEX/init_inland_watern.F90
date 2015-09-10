!     #############################################################
      SUBROUTINE INIT_INLAND_WATER_n (CHW, DTCO, DGU, DGW, IOB, UG, U, WG, W, WSB, &
                                      CHF, DGF, DGMF, FG, F, FSB, DGL,        &
                                      HPROGRAM,HINIT,                         &
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
!
!
USE MODD_CH_WATFLUX_n, ONLY : CH_WATFLUX_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_DIAG_WATFLUX_n, ONLY : DIAG_WATFLUX_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_WATFLUX_GRID_n, ONLY : WATFLUX_GRID_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
USE MODD_WATFLUX_SBL_n, ONLY : WATFLUX_SBL_t
USE MODD_CH_FLAKE_n, ONLY : CH_FLAKE_t
USE MODD_DIAG_FLAKE_n, ONLY : DIAG_FLAKE_t
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DIAG_MISC_FLAKE_t
USE MODD_FLAKE_GRID_n, ONLY : FLAKE_GRID_t
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_FLAKE_SBL_n, ONLY : FLAKE_SBL_t
USE MODD_DIAG_IDEAL_n, ONLY : DIAG_IDEAL_t
!
!
USE MODD_CSTS,       ONLY : XTT
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
TYPE(CH_WATFLUX_t), INTENT(INOUT) :: CHW
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(DIAG_WATFLUX_t), INTENT(INOUT) :: DGW
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(WATFLUX_GRID_t), INTENT(INOUT) :: WG
TYPE(WATFLUX_t), INTENT(INOUT) :: W
TYPE(WATFLUX_SBL_t), INTENT(INOUT) :: WSB
TYPE(CH_FLAKE_t), INTENT(INOUT) :: CHF
TYPE(DIAG_FLAKE_t), INTENT(INOUT) :: DGF
TYPE(DIAG_MISC_FLAKE_t), INTENT(INOUT) :: DGMF
TYPE(FLAKE_GRID_t), INTENT(INOUT) :: FG
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(FLAKE_SBL_t), INTENT(INOUT) :: FSB
TYPE(DIAG_IDEAL_t), INTENT(INOUT) :: DGL
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
  CALL INIT_WATFLUX_n(CHW, DTCO, DGU, DGW, IOB, UG, U, WG, W, WSB, &
                      HPROGRAM,HINIT,KI,KSV,KSW,HSV,PCO2,PRHOA,     &
                        PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB,    &
                        PEMIS,PTSRAD,PTSURF,                          &
                        KYEAR,KMONTH,KDAY,PTIME,HATMFILE,HATMFILETYPE,&
                        'OK'                                          )  
ELSE IF (U%CWATER=='FLAKE ') THEN
  CALL INIT_FLAKE_n(CHF, DTCO, DGF, DGMF, DGU, FG, F, FSB, IOB, &
                    UG, U, &
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
