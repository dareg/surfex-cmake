!############################################################################
SUBROUTINE INIT_COUPLING_SURF_TRIP_n(HPROGRAM, KI, KSW, OFLOOD_ISBA, HGRID, &
                                       KDIMTAB, PZENITH, PSW_BANDS, PEMIS,    &
                                       PTSRAD, PDIR_ALB, PSCA_ALB             )  
!############################################################################
!
!!****  *INIT_COUPLING_SURF_TRIP_n* - routine to test and initialyse the SURFACE-TRIP coupling
!!
!!    PURPOSE
!!    -------
!!
!     For instance, no interpolation between the ISBA and the TRIP grids that
!     must be the same (LAT LON REG 1° or 0.5°). In the near future some
!     interpolation will be added.
!
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
!!	B. Decharme   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2008
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_COUPLING_VAR_SFX_TRIP
!
USE MODI_GET_LUOUT
USE MODI_GET_GRID_CONF_ISBA_n
USE MODI_GET_GRID_CONF_TRIP_n
USE MODI_GET_CONF_TRIP_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_UPDATE_ESM_SURF_ATM_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM    ! program calling surf. schemes
CHARACTER(LEN=10),                INTENT(IN)  :: HGRID       
LOGICAL,                          INTENT(IN)  :: OFLOOD_ISBA
INTEGER,                          INTENT(IN)  :: KDIMTAB
INTEGER,                          INTENT(IN)  :: KI         ! Surfex grid dimension
INTEGER,                          INTENT(IN)  :: KSW        ! Number of spectral bands
!
REAL,             DIMENSION(KI),  INTENT(IN)  :: PZENITH   ! solar zenithal angle
REAL,             DIMENSION(KSW), INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),  INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSRAD    ! radiative temperature
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
! 1D TRIP Dimension
!
REAL, DIMENSION(:), ALLOCATABLE      :: ZTRIP_FFLOOD    ! Flooded fraction for ISBA (-)
REAL, DIMENSION(:), ALLOCATABLE      :: ZTRIP_PIFLOOD   ! Flood potential infiltration for ISBA (kg/m2/s)
!
! 1D SFX Dimension
!
REAL, DIMENSION(:), ALLOCATABLE      :: ZSFX_FFLOOD    ! Flooded fraction for ISBA (-)
REAL, DIMENSION(:), ALLOCATABLE      :: ZSFX_PIFLOOD   ! Flood potential infiltration for ISBA (kg/m2/s)
!
! Other
!
LOGICAL :: LFLOOD_TRIP
!
REAL    :: ZLONMIN_ISBA
REAL    :: ZLONMIN_TRIP
REAL    :: ZLONMAX_ISBA
REAL    :: ZLONMAX_TRIP
REAL    :: ZLATMIN_ISBA
REAL    :: ZLATMIN_TRIP
REAL    :: ZLATMAX_ISBA
REAL    :: ZLATMAX_TRIP
REAL    :: ZRES_ISBA
REAL    :: ZRES_TRIP
REAL    :: ZT_COUPLING
!
INTEGER :: ICOMP_TRIP
INTEGER :: ICOMP_ISBA
!
INTEGER :: ILON_ISBA
INTEGER :: ILON_TRIP
INTEGER :: ILAT_ISBA
INTEGER :: ILAT_TRIP
!
INTEGER :: ILUOUT, JL, INI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_COUPLING_SURF_TRIP_N',0,ZHOOK_HANDLE)
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!*       1.0   Test ISBA - TRIP option compatibility
!              -------------------------------------
!
CALL GET_CONF_TRIP_n(LFLOOD_TRIP,ZT_COUPLING)
!
ICOMP_TRIP=0
IF(LFLOOD_TRIP)ICOMP_TRIP=1
!
ICOMP_ISBA=0
IF(OFLOOD_ISBA)ICOMP_ISBA=1
!
IF(ICOMP_TRIP/=ICOMP_ISBA)THEN
  WRITE(ILUOUT,*)'INIT_COUPLING_SURF_TRIPN: IN TRIP LFLOODT = ',LFLOOD_TRIP
  WRITE(ILUOUT,*)'INIT_COUPLING_SURF_TRIPN: IN ISBA LFLOOD  = ',OFLOOD_ISBA
  CALL ABOR1_SFX('INIT_COUPLING_SURF_TRIPN: ERROR OF (DE)ACTIVATION OF FLOODING SCHEME')
ENDIF
!
!*       2.0   Test ISBA - TRIP grid compatibility
!              -------------------------------------
!
IF(HGRID/="LONLAT REG")THEN
  CALL ABOR1_SFX('INIT_COUPLING_SURF_TRIPN: ISBA-TRIP REQUIRE LONLAT REG GRID')
ENDIF
!
CALL GET_GRID_CONF_ISBA_n(ZLONMIN_ISBA,ZLONMAX_ISBA,ZLATMIN_ISBA,ZLATMAX_ISBA,ILON_ISBA,ILAT_ISBA,JL)
ZRES_ISBA = (ZLONMAX_ISBA - ZLONMIN_ISBA) / ILON_ISBA
!
CALL GET_GRID_CONF_TRIP_n(ZLONMIN_TRIP,ZLONMAX_TRIP,ZLATMIN_TRIP,ZLATMAX_TRIP,ZRES_TRIP,ILON_TRIP,ILAT_TRIP)
!
IF(ZRES_TRIP/=ZRES_ISBA)THEN
  WRITE(ILUOUT,*)'INIT_COUPLING_SURF_TRIPN: TRIP RESOLUTION = ',ZRES_TRIP
  WRITE(ILUOUT,*)'INIT_COUPLING_SURF_TRIPN: ISBA RESOLUTION = ',ZRES_ISBA
  CALL ABOR1_SFX('INIT_COUPLING_SURF_TRIPN: ISBA AND TRIP REQUIRE SAME RESOLUTION GRID')
ENDIF
!
IF(ZLONMIN_ISBA/=ZLONMIN_TRIP.OR.ZLONMAX_ISBA/=ZLONMAX_TRIP)THEN
   CALL ABOR1_SFX('INIT_COUPLING_SURF_TRIPN: WRONG CONFIGURATION FOR LONGITUDE')
ENDIF
!
IF(ZLATMIN_ISBA/=ZLATMIN_TRIP.OR.ZLATMAX_ISBA/=ZLATMAX_TRIP)THEN
  CALL ABOR1_SFX('INIT_COUPLING_SURF_TRIPN: WRONG CONFIGURATION FOR LATITUDE')
ENDIF        
!
!
!
!  
!*       3.0   Initialyse ISBA - TRIP variable
!              -------------------------------
!
IF(OFLOOD_ISBA.AND.LFLOOD_TRIP)THEN
!
  INI=ILON_TRIP*ILAT_TRIP
!
  ALLOCATE(ZTRIP_FFLOOD   (INI))
  ALLOCATE(ZTRIP_PIFLOOD  (INI))
!  
  ALLOCATE(ZSFX_FFLOOD   (KI))
  ALLOCATE(ZSFX_PIFLOOD  (KI))
!  
! in TRIP dimension INI
  CALL GET_COUPLING_VAR_TRIP_n(ZTRIP_FFLOOD(:),ZTRIP_PIFLOOD(:))
!  
!   Interpolation from TRIP grid to SFX grid
    IF(KI==INI)THEN
      ZSFX_FFLOOD (:) = ZTRIP_FFLOOD (:)
      ZSFX_PIFLOOD(:) = ZTRIP_PIFLOOD(:)
    ELSE
      CALL ABOR1_SFX('COUPLING_SURF_TRIPN: TRIP and SFX are not on the same grid')    
    ENDIF
!
    WHERE(ZSFX_FFLOOD (:)<0.01)
          ZSFX_FFLOOD (:)=0.0
          ZSFX_PIFLOOD(:)=0.0
    ENDWHERE
!    
  CALL PUT_COUPLING_VAR_SFX_n(ZSFX_FFLOOD,ZSFX_PIFLOOD,PTSTEP=ZT_COUPLING)
! 
  DEALLOCATE(ZTRIP_FFLOOD )
  DEALLOCATE(ZTRIP_PIFLOOD)
  DEALLOCATE(ZSFX_FFLOOD  )
  DEALLOCATE(ZSFX_PIFLOOD )
!
  CALL UPDATE_ESM_SURF_ATM_n(HPROGRAM, KI, KSW, PZENITH, PSW_BANDS,   &
                               PTSRAD, PDIR_ALB, PSCA_ALB, PEMIS        )  
!
ENDIF
IF (LHOOK) CALL DR_HOOK('INIT_COUPLING_SURF_TRIP_N',1,ZHOOK_HANDLE)
!   
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_COUPLING_SURF_TRIP_n      
