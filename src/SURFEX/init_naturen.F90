!     #############################################################
      SUBROUTINE INIT_NATURE_n (AG, CHI, DTCO, DTI, DTZ, DGEI, DGI, DGL, DGMI, &
                                DGU, DST, GB, IOB, ICP, IG, I, SLT, UG, U, SV, &
                                  HPROGRAM,HINIT,OLAND_USE,                    &
                                   KI,KSV,KSW,                                &
                                   HSV,PCO2,PRHOA,                            &
                                   PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                                   PEMIS,PTSRAD,PTSURF,                       &
                                   KYEAR, KMONTH,KDAY, PTIME,                 &
                                   HATMFILE,HATMFILETYPE,                     &
                                   HTEST                                      )  
!     #############################################################
!
!!****  *INIT_NATURE_n* - routine to choose initialization of vegetation scheme
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
!!      Original    27/09/96 
!!       V.Masson   18/08/97 call to fmread directly with dates and strings
!!       V.Masson   15/03/99 new PGD treatment with COVER types
!        F.Solmon  06/00   adaptation for patch approach
!!      B. Decharme  04/2013 new coupling variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_IDEAL_n, ONLY : DIAG_IDEAL_t
USE MODD_DIAG_ISBA_n, ONLY : DIAG_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_ISBA_CANOPY_n, ONLY : ISBA_CANOPY_t
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SLT_n, ONLY : SLT_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SV_n, ONLY : SV_t
!
!
USE MODD_CSTS,       ONLY : XTT
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_INIT_IDEAL_FLUX
!
USE MODI_INIT_ISBA_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
!
TYPE(AGRI_t), INTENT(INOUT) :: AG
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_IDEAL_t), INTENT(INOUT) :: DGL
TYPE(DIAG_ISBA_t), INTENT(INOUT) :: DGI
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(ISBA_CANOPY_t), INTENT(INOUT) :: ICP
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(SLT_t), INTENT(INOUT) :: SLT
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SV_t), INTENT(INOUT) :: SV
!
 CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),                 INTENT(IN)  :: HINIT     ! choice of fields to initialize
LOGICAL,                          INTENT(IN)  :: OLAND_USE ! 
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
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
!*       2.     Selection of surface scheme
!               ---------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_NATURE_N',0,ZHOOK_HANDLE)
IF (U%CNATURE=='NONE  ') THEN
  PDIR_ALB=0.
  PSCA_ALB=0.
  PEMIS   =1.
  PTSRAD  =XTT
  PTSURF  =XTT
ELSE IF (U%CNATURE=='FLUX  ') THEN
  CALL INIT_IDEAL_FLUX(DGL, &
                       HPROGRAM,HINIT,KI,KSV,KSW,HSV,PCO2,PRHOA,     &
                           PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB,  &
                           PEMIS,PTSRAD,PTSURF,'OK'                    )  
ELSE IF (U%CNATURE=='ISBA  ' .OR. U%CNATURE=='TSZ0') THEN
  CALL INIT_ISBA_n(AG, CHI, DTCO, DTI, DTZ, DGEI, DGI, DGMI, &
                   DGU, DST, GB, IOB, ICP, IG, I, SLT, UG, U, SV, &
                   HPROGRAM,HINIT,OLAND_USE,KI,KSV,KSW,HSV,PCO2,PRHOA,     &
                     PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB,    &
                     PEMIS,PTSRAD,PTSURF,                          &
                     KYEAR,KMONTH,KDAY,PTIME,HATMFILE,HATMFILETYPE,&
                     'OK'                                          )  
END IF
IF (LHOOK) CALL DR_HOOK('INIT_NATURE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_NATURE_n
