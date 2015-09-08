!#############################################################
SUBROUTINE INIT_SURF_LANDUSE_n (YSC, &
                                HPROGRAM,HINIT,OLAND_USE,                  &
                               KI,KSV,KSW,                                &
                               HSV,PCO2,PRHOA,                            &
                               PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                               PEMIS,PTSRAD,PTSURF,                       &
                               KYEAR, KMONTH,KDAY, PTIME,                 &
                               HATMFILE,HATMFILETYPE,                     &
                               HTEST                                      )  
!#############################################################
!
!!****  *INIT_SURF_LANDUSE_n* - routine to initialize LAND USE 
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
!!    S. Faroux    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      modified    06-13  B. Decharme  : New coupling variable
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_SURFEX_n, ONLY : SURFEX_t
!
USE YOMHOOK   ,   ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
!
!
USE MODI_INIT_IO_SURF_n
USE MODI_END_IO_SURF_n
!
USE MODI_GET_TYPE_DIM_n
USE MODI_READ_SURF
!
USE MODI_SET_VEGTYPES_FRACTIONS
USE MODI_COMPUTE_ISBA_PARAMETERS
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(SURFEX_t), INTENT(INOUT) :: YSC
!
 CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),                 INTENT(IN)  :: HINIT     ! choice of fields to initialize
LOGICAL,                          INTENT(IN)  :: OLAND_USE ! choice of doing land use or not 
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
!
INTEGER,                          INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,                          INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,                          INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                             INTENT(IN)  :: PTIME     ! current time since
                                                           !  midnight (UTC, s)
!
 CHARACTER(LEN=28),                INTENT(IN)  :: HATMFILE    ! atmospheric file name
 CHARACTER(LEN=6),                 INTENT(IN)  :: HATMFILETYPE! atmospheric file type
 CHARACTER(LEN=2),                 INTENT(IN)  :: HTEST       ! must be equal to 'OK'
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
REAL, DIMENSION(:,:),ALLOCATABLE  :: ZWORK      ! 2D array to write data in file
INTEGER           :: JLAYER
INTEGER           :: ILU          ! 1D physical dimension
INTEGER           :: IRESP          ! Error code after redding
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=4)  :: YLVL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_SURF_LANDUSE_N',0,ZHOOK_HANDLE)
!
IF (HTEST/='OK') THEN
   CALL ABOR1_SFX('INIT_SURF_LANDUSEN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
IF (.NOT. OLAND_USE)THEN
   IF (LHOOK) CALL DR_HOOK('INIT_SURF_LANDUSE_N',1,ZHOOK_HANDLE)
   RETURN
ENDIF
!
IF (YSC%I%CISBA=='DIF') THEN
   CALL ABOR1_SFX('INIT_SURF_LANDUSEN: LAND USE NOT IMPLEMENTED WITH DIF')
ENDIF
!
!-------------------------------------------------------------------------------
!
!* initialization for I/O
!
 CALL INIT_IO_SURF_n(YSC%BOP, YSC%BDD, YSC%CHE, YSC%CHI, YSC%CHS, YSC%CHN, YSC%CHU, YSC%CHT, YSC%CHW, &
                     YSC%DTCO, YSC%DTS, YSC%DTT, YSC%DTZ, YSC%DGEI, YSC%DGF, YSC%DGI, YSC%DGMI, YSC%DGMTO, &
                     YSC%DGO, YSC%DGS, YSC%DGSI, YSC%DGU, YSC%DGT, YSC%DGUT, YSC%DGW, YSC%F, YSC%FSB, &
                     YSC%GB, YSC%IOB, YSC%ICP, YSC%I, YSC%O, YSC%S, YSC%SSB, YSC%UG, YSC%U, YSC%SV, &
                     YSC%TCP, YSC%TGD, YSC%TGDO, YSC%TGR, YSC%TGRO, YSC%T, YSC%TOP, YSC%TVG, YSC%W, YSC%WSB, &
                        HPROGRAM,'NATURE','ISBA  ','READ ')
!
!* 1D physical dimension
!
 CALL GET_TYPE_DIM_n(YSC%DTCO, YSC%U, &
                     'NATURE',ILU)
ALLOCATE(ZWORK(ILU,YSC%I%NPATCH))
!
YSC%DTI%LDATA_MIXPAR = .TRUE.
IF (.NOT.ASSOCIATED(YSC%DTI%XPAR_VEGTYPE)) ALLOCATE(YSC%DTI%XPAR_VEGTYPE(ILU,NVEGTYPE))
IF (YSC%DTI%NTIME==0) YSC%DTI%NTIME = 36
IF (.NOT.ASSOCIATED(YSC%DTI%XPAR_LAI)) ALLOCATE(YSC%DTI%XPAR_LAI(ILU,YSC%DTI%NTIME,NVEGTYPE))
IF (.NOT.ASSOCIATED(YSC%DTI%XPAR_H_TREE)) ALLOCATE(YSC%DTI%XPAR_H_TREE(ILU,NVEGTYPE))
IF (.NOT.ASSOCIATED(YSC%DTI%XPAR_ROOT_DEPTH)) ALLOCATE(YSC%DTI%XPAR_ROOT_DEPTH(ILU,NVEGTYPE))
IF (.NOT.ASSOCIATED(YSC%DTI%XPAR_GROUND_DEPTH)) ALLOCATE(YSC%DTI%XPAR_GROUND_DEPTH(ILU,NVEGTYPE))
IF (.NOT.ASSOCIATED(YSC%DTI%XPAR_IRRIG)) ALLOCATE(YSC%DTI%XPAR_IRRIG(ILU,YSC%DTI%NTIME,NVEGTYPE))
IF (.NOT.ASSOCIATED(YSC%DTI%XPAR_WATSUP)) ALLOCATE(YSC%DTI%XPAR_WATSUP(ILU,YSC%DTI%NTIME,NVEGTYPE))
!
!* read old patch fraction
!       
ALLOCATE(YSC%I%XPATCH_OLD(ILU,YSC%I%NPATCH))       
YRECFM = 'PATCH'
 CALL READ_SURF(YSC%IOB, &
                HPROGRAM,YRECFM,YSC%I%XPATCH_OLD(:,:),IRESP)
!
!* read old soil layer thicknesses (m)
!
ALLOCATE(YSC%I%XDG_OLD(ILU,YSC%I%NGROUND_LAYER,YSC%I%NPATCH))
!
DO JLAYER=1,YSC%I%NGROUND_LAYER
  WRITE(YLVL,'(I4)') JLAYER
  YRECFM='OLD_DG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(YSC%IOB, &
                HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  YSC%I%XDG_OLD(:,JLAYER,:)=ZWORK
END DO
DEALLOCATE(ZWORK)
!
!* End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!* read new fraction of each vege type
! and then extrapolate parameters defined by cover
!       
 CALL SET_VEGTYPES_FRACTIONS(YSC%BOP, YSC%BDD, YSC%CHE, YSC%CHI, YSC%CHS, YSC%CHN, YSC%CHU, YSC%CHT, &
                             YSC%CHW, YSC%DTS, YSC%DTT, YSC%DTZ, YSC%DGEI, YSC%DGF, YSC%DGI, YSC%DGMI, &
                             YSC%DGMTO, YSC%DGO, YSC%DGS, YSC%DGSI, YSC%DGU, YSC%DGT, YSC%DGUT, &
                             YSC%DGW, YSC%F, YSC%FSB, YSC%GB, YSC%ICP, YSC%O, YSC%S, YSC%SSB, YSC%SV, &
                             YSC%TCP, YSC%TGD, YSC%TGDO, YSC%TGR, YSC%TGRO, YSC%T, YSC%TOP, YSC%TVG, &
                             YSC%W, YSC%WSB, YSC%DTCO, YSC%DTI, YSC%IOB, YSC%IG, YSC%I, YSC%UG, YSC%U, &
                             HPROGRAM)
!
!* re-initialize ISBA with new parameters
!       
 CALL COMPUTE_ISBA_PARAMETERS(YSC%AG, YSC%BOP, YSC%BDD, YSC%CHE, YSC%CHI, YSC%CHS, YSC%CHN, YSC%CHU, &
                              YSC%CHT, YSC%CHW, YSC%DTCO, YSC%DTI, YSC%DTS, YSC%DTT, YSC%DTZ, YSC%DGEI, &
                              YSC%DGF, YSC%DGI, YSC%DGMI, YSC%DGMTO, YSC%DGO, YSC%DGS, YSC%DGSI, &
                              YSC%DGU, YSC%DGT, YSC%DGUT, YSC%DGW, YSC%DST, YSC%F, YSC%FSB, YSC%GB, &
                              YSC%IOB, YSC%ICP, YSC%IG, YSC%I, YSC%O, YSC%S, YSC%SSB, YSC%SLT, YSC%UG, &
                              YSC%U, YSC%SV, YSC%TCP, YSC%TGD, YSC%TGDO, YSC%TGR, YSC%TGRO, YSC%T, &
                              YSC%TOP, YSC%TVG, YSC%W, YSC%WSB, &
                              HPROGRAM,HINIT,OLAND_USE,                  &
                             ILU,KSV,KSW,                                &
                             HSV,PCO2,PRHOA,                            &
                             PZENITH,PSW_BANDS,PDIR_ALB,PSCA_ALB,       &
                             PEMIS,PTSRAD,PTSURF,                       &
                             HTEST                                      )
!-------------------------------------------------------------------------------
!                       
IF (LHOOK) CALL DR_HOOK('INIT_SURF_LANDUSE_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_SURF_LANDUSE_n                           
