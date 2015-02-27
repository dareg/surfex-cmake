!#############################################################
SUBROUTINE INIT_SURF_LANDUSE_n(HPROGRAM,HINIT,OLAND_USE,                  &
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
USE YOMHOOK   ,   ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
!
USE MODD_ISBA_n,  ONLY : XPATCH_OLD,XDG_OLD,CISBA,NGROUND_LAYER, NPATCH
USE MODD_DATA_ISBA_n,    ONLY : XPAR_VEGTYPE, XPAR_LAI, XPAR_H_TREE, XPAR_ROOT_DEPTH, &
                                XPAR_GROUND_DEPTH, XPAR_IRRIG, XPAR_WATSUP, &
                                LDATA_MIXPAR, NTIME
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
IF (CISBA=='DIF') THEN
   CALL ABOR1_SFX('INIT_SURF_LANDUSEN: LAND USE NOT IMPLEMENTED WITH DIF')
ENDIF
!
!-------------------------------------------------------------------------------
!
!* initialization for I/O
!
 CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','ISBA  ','READ ')
!
!* 1D physical dimension
!
 CALL GET_TYPE_DIM_n('NATURE',ILU)
ALLOCATE(ZWORK(ILU,NPATCH))
!
LDATA_MIXPAR = .TRUE.
IF (.NOT.ASSOCIATED(XPAR_VEGTYPE)) ALLOCATE(XPAR_VEGTYPE(ILU,NVEGTYPE))
IF (NTIME==0) NTIME = 36
IF (.NOT.ASSOCIATED(XPAR_LAI)) ALLOCATE(XPAR_LAI(ILU,NTIME,NVEGTYPE))
IF (.NOT.ASSOCIATED(XPAR_H_TREE)) ALLOCATE(XPAR_H_TREE(ILU,NVEGTYPE))
IF (.NOT.ASSOCIATED(XPAR_ROOT_DEPTH)) ALLOCATE(XPAR_ROOT_DEPTH(ILU,NVEGTYPE))
IF (.NOT.ASSOCIATED(XPAR_GROUND_DEPTH)) ALLOCATE(XPAR_GROUND_DEPTH(ILU,NVEGTYPE))
IF (.NOT.ASSOCIATED(XPAR_IRRIG)) ALLOCATE(XPAR_IRRIG(ILU,NTIME,NVEGTYPE))
IF (.NOT.ASSOCIATED(XPAR_WATSUP)) ALLOCATE(XPAR_WATSUP(ILU,NTIME,NVEGTYPE))
!
!* read old patch fraction
!       
ALLOCATE(XPATCH_OLD(ILU,NPATCH))       
YRECFM = 'PATCH'
 CALL READ_SURF(HPROGRAM,YRECFM,XPATCH_OLD(:,:),IRESP)
!
!* read old soil layer thicknesses (m)
!
ALLOCATE(XDG_OLD(ILU,NGROUND_LAYER,NPATCH))
!
DO JLAYER=1,NGROUND_LAYER
  WRITE(YLVL,'(I4)') JLAYER
  YRECFM='OLD_DG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  XDG_OLD(:,JLAYER,:)=ZWORK
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
 CALL SET_VEGTYPES_FRACTIONS(HPROGRAM)
!
!* re-initialize ISBA with new parameters
!       
 CALL COMPUTE_ISBA_PARAMETERS(HPROGRAM,HINIT,OLAND_USE,                  &
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
