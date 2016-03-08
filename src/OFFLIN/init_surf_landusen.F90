!#############################################################
SUBROUTINE INIT_SURF_LANDUSE_n (DTCO, DUO, U, UG, IM, SV, SLT, DST, &
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
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_DIAG_n, ONLY : DIAG_OPTIONS_t
USE MODD_SURFEX_n, ONLY : ISBA_MODEL_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_SLT_n, ONLY : SLT_t
USE MODD_SV_n, ONLY : SV_t
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM
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
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(DATA_COVER_t) :: DTCO
TYPE(DIAG_OPTIONS_t) :: DUO
TYPE(SURF_ATM_t) :: U
TYPE(SURF_ATM_GRID_t) :: UG
TYPE(ISBA_MODEL_t) :: IM
TYPE(SV_t), INTENT(INOUT) :: SV
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(SLT_t), INTENT(INOUT) :: SLT
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
INTEGER           :: JLAYER, INFOMPI
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
IF (IM%I%O%CISBA=='DIF') THEN
   CALL ABOR1_SFX('INIT_SURF_LANDUSEN: LAND USE NOT IMPLEMENTED WITH DIF')
ENDIF
!
!-------------------------------------------------------------------------------
!
CALL MPI_BCAST(UG%NGRID_FULL_PAR,KIND(UG%NGRID_FULL_PAR)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
IF (NRANK/=NPIO) ALLOCATE(UG%XGRID_FULL_PAR(UG%NGRID_FULL_PAR))
 CALL MPI_BCAST(UG%XGRID_FULL_PAR,&
                SIZE(UG%XGRID_FULL_PAR)*KIND(UG%XGRID_FULL_PAR)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
!   
!* initialization for I/O
!
CALL INIT_IO_SURF_n(DTCO, U, HPROGRAM,'NATURE','ISBA  ','READ ')
!
!* 1D physical dimension
!
 CALL GET_TYPE_DIM_n(DTCO, U, 'NATURE',ILU)
ALLOCATE(ZWORK(ILU,IM%I%O%NPATCH))
!
IM%DTI%LDATA_MIXPAR = .TRUE.
IF (.NOT.ASSOCIATED(IM%DTI%XPAR_VEGTYPE)) ALLOCATE(IM%DTI%XPAR_VEGTYPE(ILU,NVEGTYPE))
IF (IM%DTI%NTIME==0) IM%DTI%NTIME = 36
IF (.NOT.ASSOCIATED(IM%DTI%XPAR_LAI)) ALLOCATE(IM%DTI%XPAR_LAI(ILU,IM%DTI%NTIME,NVEGTYPE))
IF (.NOT.ASSOCIATED(IM%DTI%XPAR_H_TREE)) ALLOCATE(IM%DTI%XPAR_H_TREE(ILU,NVEGTYPE))
IF (.NOT.ASSOCIATED(IM%DTI%XPAR_ROOT_DEPTH)) ALLOCATE(IM%DTI%XPAR_ROOT_DEPTH(ILU,NVEGTYPE))
IF (.NOT.ASSOCIATED(IM%DTI%XPAR_GROUND_DEPTH)) ALLOCATE(IM%DTI%XPAR_GROUND_DEPTH(ILU,NVEGTYPE))
IF (.NOT.ASSOCIATED(IM%DTI%XPAR_IRRIG)) ALLOCATE(IM%DTI%XPAR_IRRIG(ILU,IM%DTI%NTIME,NVEGTYPE))
IF (.NOT.ASSOCIATED(IM%DTI%XPAR_WATSUP)) ALLOCATE(IM%DTI%XPAR_WATSUP(ILU,IM%DTI%NTIME,NVEGTYPE))
!
!* read old patch fraction
!       
ALLOCATE(IM%I%IP%XPATCH_OLD(ILU,IM%I%O%NPATCH))       
YRECFM = 'PATCH'
 CALL READ_SURF(HPROGRAM,YRECFM,IM%I%IP%XPATCH_OLD(:,:),IRESP)
!
!* read old soil layer thicknesses (m)
!
ALLOCATE(IM%I%M%X%XDG_OLD(ILU,IM%I%O%NGROUND_LAYER,IM%I%O%NPATCH))
!
DO JLAYER=1,IM%I%O%NGROUND_LAYER
  WRITE(YLVL,'(I4)') JLAYER
  YRECFM='OLD_DG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  IM%I%M%X%XDG_OLD(:,JLAYER,:)=ZWORK
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
 CALL SET_VEGTYPES_FRACTIONS(DTCO, IM%DTI, IM%IG%NDIM, IM%I%O, IM%I%P, UG, U, HPROGRAM)
!
!* re-initialize ISBA with new parameters
!       ²
 CALL COMPUTE_ISBA_PARAMETERS(DTCO, DUO%LREAD_BUDGETC, UG, U, &
                              IM%AG, IM%CHI, IM%DTI, IM%DGI, IM%GB, &
                              IM%ICP, IM%ISS, IM%IG, IM%I, DST, &
                              SLT, SV, HPROGRAM, HINIT, OLAND_USE, &
                              ILU,KSV,KSW, HSV,PCO2,PRHOA,  &
                             PZENITH,PSW_BANDS,PDIR_ALB,PSCA_ALB,      &
                             PEMIS,PTSRAD,PTSURF,HTEST                 )
!-------------------------------------------------------------------------------
!                       
IF (LHOOK) CALL DR_HOOK('INIT_SURF_LANDUSE_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_SURF_LANDUSE_n                           
