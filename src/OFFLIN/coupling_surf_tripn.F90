!######
SUBROUTINE COUPLING_SURF_TRIP_n (HPROGRAM,KI,KSW,ORESTART,KYEAR,KMONTH,KTRIP, &
                                   PTIME,PDURATION,PZENITH,PSW_BANDS,PEMIS,     &
                                   PTSRAD,PDIR_ALB,PSCA_ALB                     )  
!###################################################################
!
!!****  *COUPLING_SURF_TRIP_n*  
!!
!!    PURPOSE
!!    -------
!!   
!!    Driver for the coupling between SURFEX and TRIP
!!      
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!	B. Decharme     
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/02/05 
!!      Modif.      28/05/08 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TRIP_GRID_n, ONLY : GMASK
USE MODE_COUPLING_VAR_SFX_TRIP
!
USE MODI_GET_LUOUT
USE MODI_GET_CONF_ISBA_n
USE MODI_GET_CONF_TRIP_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_TRIP_INTERFACE
USE MODI_UPDATE_ESM_SURF_ATM_n
!
USE MODI_ABOR1_SFX
USE MODI_GET_TRIP_SIZE_n
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=6), INTENT(IN)         :: HPROGRAM ! program calling surf. schemes
INTEGER,          INTENT(IN)         :: KI       ! Surfex grid dimension
INTEGER,          INTENT(IN)         :: KSW      ! Number of spectral bands
!
LOGICAL, INTENT(IN)                  :: ORESTART ! write restart file
!
INTEGER, INTENT(IN)                  :: KYEAR    ! current year (UTC)
INTEGER, INTENT(IN)                  :: KMONTH   ! current month (UTC)
INTEGER, INTENT(INOUT)               :: KTRIP    ! number of Trip timestep counter
!
REAL,    INTENT(IN)                  :: PTIME    ! current time since start of the run (s)
REAL,    INTENT(IN)                  :: PDURATION! duration of run                     (s)
!
REAL,             DIMENSION(KI),  INTENT(IN)  :: PZENITH   ! solar zenithal angle
REAL,             DIMENSION(KSW), INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),  INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSRAD    ! radiative temperature
!
!*      0.2    declarations of local variables
!
LOGICAL :: LTRIP
LOGICAL :: LFLOOD
!
REAL    :: ZTSTEP_COUPLING
!
INTEGER :: ILON  !TRIP lon dimension
INTEGER :: ILAT  !TRIP lat dimension
INTEGER :: INI   !TRIP total dimension
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: ILUOUT
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SURF_TRIP_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
! * 1. Get ISBA and TRIP  configuration
!      
 CALL GET_CONF_ISBA_n(LTRIP,LFLOOD)
!
IF (.NOT.LTRIP .AND. LHOOK) CALL DR_HOOK('COUPLING_SURF_TRIP_N',1,ZHOOK_HANDLE)
IF (.NOT.LTRIP) RETURN
!
 CALL GET_CONF_TRIP_n(PTSTEP_COUPLING=ZTSTEP_COUPLING)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
IF(LTRIP.AND.MOD(PTIME,ZTSTEP_COUPLING) == 0)THEN
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
  KTRIP = KTRIP + 1
!
  CALL GET_TRIP_SIZE_n(ILON,ILAT)
!
  INI = ILON*ILAT
!
  CALL COUPLING_SURF_TRIP_DIM(ILON,ILAT)
!
!-------------------------------------------------------------------------------
!
! * 2. Put TRIP variables into SURFEX if flooding scheme
!
  IF(LFLOOD) CALL COUPLING_SURF_TRIP_FLOOD(INI)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SURF_TRIP_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!
SUBROUTINE COUPLING_SURF_TRIP_DIM(KLON,KLAT)
!
INTEGER, INTENT(IN) :: KLON
INTEGER, INTENT(IN) :: KLAT
!
! 2D LAT/LON Dimension
!
REAL, DIMENSION(KLON,KLAT)    :: Z2D_DRAIN       ! Drainage for TRIP (kg/m2)      
REAL, DIMENSION(KLON,KLAT)    :: Z2D_RUNOFF      ! Runoff for TRIP (kg/m2)   
REAL, DIMENSION(KLON,KLAT)    :: Z2D_SRC_FLOOD   ! Flood source budget for TRIP (kg/m2)
!
! 1D TRIP Dimension
!
REAL, DIMENSION(KLON*KLAT)    :: ZTRIP_DRAIN     ! Cumulative drainage for TRIP (kg)
REAL, DIMENSION(KLON*KLAT)    :: ZTRIP_RUNOFF    ! Cumulative Runoff for TRIP (kg)
REAL, DIMENSION(KLON*KLAT)    :: ZTRIP_SRC_FLOOD ! Cumulative flood budget for TRIP (kg)
!
! 1D SFX Dimension
!
REAL, DIMENSION(KI)    :: ZSFX_DRAIN     ! Cumulative drainage from ISBA (kg)
REAL, DIMENSION(KI)    :: ZSFX_RUNOFF    ! Cumulative Runoff from ISBA (kg)
REAL, DIMENSION(KI)    :: ZSFX_SRC_FLOOD ! Cumulative flood budget from ISBA (kg)
!
INTEGER :: I,J,ICOUNT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SURF_TRIP_N:COUPLING_SURF_TRIP_DIM',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------  
!
! * 1. Initialyse local variables
!
! 2D LAT/LON Dimension
!
Z2D_DRAIN    (:,:) = 0.0
Z2D_RUNOFF   (:,:) = 0.0
Z2D_SRC_FLOOD(:,:) = 0.0
!
! 1D TRIP Dimension
! 
ZTRIP_DRAIN    (:) = 0.0
ZTRIP_RUNOFF   (:) = 0.0
ZTRIP_SRC_FLOOD(:) = 0.0
! 
ZSFX_DRAIN    (:) = 0.0
ZSFX_RUNOFF   (:) = 0.0
ZSFX_SRC_FLOOD(:) = 0.0
!
!-------------------------------------------------------------------------------
!
! * 2. Get SURFEX variables for TRIP in kg
!
 CALL GET_COUPLING_VAR_SFX_n(LFLOOD,ZSFX_DRAIN(:),ZSFX_RUNOFF(:),ZSFX_SRC_FLOOD(:))
!
!-------------------------------------------------------------------------------
!
! * 3. Interpolation from SFX grid to TRIP grid (kg)
!
IF(KI==INI)THEN
  ZTRIP_RUNOFF   (:)=ZSFX_RUNOFF   (:)
  ZTRIP_DRAIN    (:)=ZSFX_DRAIN    (:)
  IF(LFLOOD)THEN
    ZTRIP_SRC_FLOOD(:)=ZSFX_SRC_FLOOD(:)   
  ENDIF
ELSE
  CALL ABOR1_SFX('COUPLING_SURF_TRIPN: TRIP and SFX are not on the same grid')
ENDIF
!
! 2d lat/lon TRIP grid array
!
ICOUNT=0
DO J=1,KLAT
  DO I=1,KLON
    ICOUNT=ICOUNT+1
      IF(.NOT.GMASK(I,J))CYCLE
      Z2D_RUNOFF   (I,J)= ZTRIP_RUNOFF   (ICOUNT)
      Z2D_DRAIN    (I,J)= ZTRIP_DRAIN    (ICOUNT)
   ENDDO
ENDDO
!
IF(LFLOOD)THEN
  ICOUNT=0
  DO J=1,KLAT
    DO I=1,KLON
      ICOUNT=ICOUNT+1
      IF(.NOT.GMASK(I,J))CYCLE
      Z2D_SRC_FLOOD(I,J)= ZTRIP_SRC_FLOOD(ICOUNT)
    ENDDO
  ENDDO
ENDIF
!
! Mass conservation
!
!
!-------------------------------------------------------------------------------
!
! * 5. Call Trip coupling
!
 CALL TRIP_INTERFACE(ILUOUT,ILON,ILAT,ORESTART,KYEAR,KMONTH,KTRIP,PDURATION,&
                      Z2D_RUNOFF(:,:),Z2D_DRAIN(:,:),Z2D_SRC_FLOOD(:,:)      )
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SURF_TRIP_N:COUPLING_SURF_TRIP_DIM',1,ZHOOK_HANDLE)
!
END SUBROUTINE COUPLING_SURF_TRIP_DIM
!
!
SUBROUTINE COUPLING_SURF_TRIP_FLOOD(KNI)
!
INTEGER, INTENT(IN) :: KNI
!
REAL, DIMENSION(KNI) :: ZTRIP_FFLOOD       ! Flooded fraction from TRIP (-)
REAL, DIMENSION(KNI) :: ZTRIP_PIFLOOD     ! Flood potential infiltration from TRIP (kg)
!
REAL, DIMENSION(KI)    :: ZSFX_FFLOOD    ! Flooded fraction for ISBA (-)
REAL, DIMENSION(KI)    :: ZSFX_PIFLOOD   ! Flood potential infiltration for ISBA (kg)
!    
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SURF_TRIP_N:COUPLING_SURF_TRIP_FLOOD',0,ZHOOK_HANDLE)
!
ZTRIP_PIFLOOD      (:) = 0.0
ZTRIP_FFLOOD       (:) = 0.0
!
ZSFX_FFLOOD   (:) = 0.0
ZSFX_PIFLOOD  (:) = 0.0
!   
!   TRIP dimension INI in kg
 CALL GET_COUPLING_VAR_TRIP_n(ZTRIP_FFLOOD(:),ZTRIP_PIFLOOD(:))
!
!   Interpolation from TRIP grid to SFX grid in kg
IF(KI==KNI)THEN
  ZSFX_FFLOOD (:) = ZTRIP_FFLOOD (:)
  ZSFX_PIFLOOD(:) = ZTRIP_PIFLOOD(:)
ELSE
  CALL ABOR1_SFX('COUPLING_SURF_TRIPN: TRIP and SFX are not on the same grid')
ENDIF
!
!   Mass conservation
!   
!
!   Put into SFX
!
WHERE(ZSFX_FFLOOD (:)<0.01)
  ZSFX_FFLOOD (:)=0.0
  ZSFX_PIFLOOD(:)=0.0
ENDWHERE
!
 CALL PUT_COUPLING_VAR_SFX_n(ZSFX_FFLOOD,ZSFX_PIFLOOD)
!
!-------------------------------------------------------------------------------
!
! * 8. Update radiative properties with flooding
!
 CALL UPDATE_ESM_SURF_ATM_n(HPROGRAM, KI, KSW, PZENITH, PSW_BANDS,   &
                           PTSRAD, PDIR_ALB, PSCA_ALB, PEMIS        )
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SURF_TRIP_N:COUPLING_SURF_TRIP_FLOOD',1,ZHOOK_HANDLE)
!
END SUBROUTINE COUPLING_SURF_TRIP_FLOOD
!
!-------------------------------------------------------------------------------
END SUBROUTINE COUPLING_SURF_TRIP_n
