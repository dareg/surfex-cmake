!     #########
SUBROUTINE PREP_TEB_GRIB(HPROGRAM,HSURF,HFILE,KLUOUT,PFIELD)
!     #################################################################################
!
!!****  *PREP_TEB_GRIB* - prepares TEB field from operational GRIB
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

!
USE MODD_TYPE_DATE_SURF
!
USE MODI_PREP_GRIB_GRID
USE MODE_READ_GRIB
USE MODI_INTERP_GRID
!
USE MODD_PREP,       ONLY : CINGRID_TYPE, CINTERP_TYPE
USE MODD_GRID_GRIB,  ONLY : CGRIB_FILE, NNI
USE MODD_PREP_TEB,   ONLY : XGRID_ROAD, XGRID_WALL, XGRID_ROOF, XTI_BLD, XTI_ROAD
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=7),   INTENT(IN)  :: HSURF     ! type of field
CHARACTER(LEN=28),  INTENT(IN)  :: HFILE     ! name of file
INTEGER,            INTENT(IN)  :: KLUOUT    ! logical unit of output listing
REAL,DIMENSION(:,:), POINTER    :: PFIELD    ! field to interpolate horizontally
!
!*      0.2    declarations of local variables
!
TYPE (DATE_TIME)                :: TZTIME_GRIB    ! current date and time
CHARACTER(LEN=6)                :: YINMODEL ! model from which GRIB file originates
REAL, DIMENSION(:)  , POINTER   :: ZMASK => NULL()          ! Land mask
REAL, DIMENSION(:),   POINTER   :: ZFIELD1D => NULL() ! 1D field read
REAL, DIMENSION(:,:), POINTER   :: ZFIELD => NULL()   ! field read
REAL, DIMENSION(:,:), POINTER   :: ZD => NULL()             ! depth of field in the soil
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!*      1.     Reading of grid
!              ---------------
!
IF (LHOOK) CALL DR_HOOK('PREP_TEB_GRIB',0,ZHOOK_HANDLE)
!
IF (TRIM(HFILE).NE.CGRIB_FILE) CGRIB_FILE=""
!
CALL PREP_GRIB_GRID(HFILE,KLUOUT,YINMODEL,CINGRID_TYPE,TZTIME_GRIB)
!
CALL READ_GRIB_LAND_MASK(HFILE,KLUOUT,YINMODEL,ZMASK)
!
!---------------------------------------------------------------------------------------
SELECT CASE(HSURF)
!---------------------------------------------------------------------------------------
!
!*     2.      Orography
!              ---------
!
  CASE('ZS     ')
    SELECT CASE (YINMODEL)
      CASE ('ECMWF ','ARPEGE','ALADIN','MOCAGE')
        CALL READ_GRIB_ZS_LAND(HFILE,KLUOUT,YINMODEL,ZMASK,ZFIELD1D)
        ALLOCATE(PFIELD(SIZE(ZFIELD1D),1))
        PFIELD(:,1) = ZFIELD1D(:)
        DEALLOCATE(ZFIELD1D)
    END SELECT
!
!*      3.     Profile of temperatures in roads
!              --------------------------------
!
  CASE('T_ROAD')
     !* reading of the profile and its depth definition
     SELECT CASE(YINMODEL)
       CASE('ECMWF ')
         CALL READ_GRIB_TG_ECMWF(HFILE,KLUOUT,YINMODEL,ZMASK,ZFIELD,ZD)
       CASE('ARPEGE','ALADIN','MOCAGE')
         CALL READ_GRIB_TG_METEO_FRANCE(HFILE,KLUOUT,YINMODEL,ZMASK,ZFIELD,ZD)
     END SELECT
     !* if deep road temperature is prescribed
     IF (XTI_ROAD/=XUNDEF) THEN
       ZFIELD(:,2:) = XTI_ROAD 
     END IF
     CALL TEB_PROFILE_GRIB(XGRID_ROAD)

!*      4.     Profile of temperatures in walls
!              --------------------------------

  CASE('T_WALL')
     CALL READ_GRIB_T_TEB(HFILE,KLUOUT,YINMODEL,XTI_BLD,ZMASK,ZFIELD,ZD)
     CALL TEB_PROFILE_GRIB(XGRID_WALL)


!*      5.     Profile of temperatures in roofs
!              --------------------------------

  CASE('T_ROOF')    
     CALL READ_GRIB_T_TEB(HFILE,KLUOUT,YINMODEL,XTI_BLD,ZMASK,ZFIELD,ZD)
     CALL TEB_PROFILE_GRIB(XGRID_ROOF)

!
!*      6.     Canyon air temperature
!              ----------------------
!
  CASE('T_CAN  ')
    SELECT CASE (YINMODEL)
      CASE ('ECMWF ','ARPEGE','ALADIN','MOCAGE')
        CALL READ_GRIB_T2(HFILE,KLUOUT,YINMODEL,ZMASK,ZFIELD1D)
        ALLOCATE(PFIELD(SIZE(ZFIELD1D),1))
        PFIELD(:,1) = ZFIELD1D(:)
        DEALLOCATE(ZFIELD1D)
    END SELECT
!
!*      7.      Canyon air humidity
!               -------------------
!
  CASE('Q_CAN  ')
    SELECT CASE (YINMODEL)
      CASE ('ECMWF ','ARPEGE','ALADIN','MOCAGE')
        ALLOCATE(PFIELD(NNI,1))
        PFIELD(:,1) = 0.
    END SELECT

!
!*      9.     Deep road temperature
!              ---------------------

  CASE('TI_ROAD')    
     IF (XTI_ROAD==XUNDEF) THEN
       CALL READ_GRIB_T2(HFILE,KLUOUT,YINMODEL,ZMASK,ZFIELD1D)
       ALLOCATE(PFIELD(SIZE(ZFIELD1D),1))
       PFIELD(:,1) = ZFIELD1D(:)
       DEALLOCATE(ZFIELD1D)
     ELSE
       ALLOCATE(PFIELD(NNI,1))
       PFIELD = XTI_ROAD
     END IF


!*      9.     Building temperature
!              --------------------

  CASE('TI_BLD ')    
     ALLOCATE(PFIELD(NNI,1))
     PFIELD = XTI_BLD

!*     10.     Other quantities (water reservoirs)
!              ----------------

  CASE DEFAULT
    ALLOCATE(PFIELD(NNI,1))
    PFIELD = 0.

END SELECT
!
DEALLOCATE(ZMASK)
!
!*      4.     Interpolation method
!              --------------------
!
CINTERP_TYPE='HORIBL'
!
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_TEB_GRIB',1,ZHOOK_HANDLE)
CONTAINS
!
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
SUBROUTINE TEB_PROFILE_GRIB(PGRID)
!-------------------------------------------------------------------------------------
!
REAL, DIMENSION(:),   INTENT(IN)  :: PGRID  ! destination grid
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!-------------------------------------------------------------------------------------
!
!* interpolation on fine vertical grid
IF (LHOOK) CALL DR_HOOK('TEB_PROFILE_GRIB',0,ZHOOK_HANDLE)
ALLOCATE(PFIELD(SIZE(ZFIELD,1),SIZE(PGRID)))
CALL INTERP_GRID(ZD,ZFIELD,PGRID,PFIELD)
!
!* end
DEALLOCATE(ZFIELD)
DEALLOCATE(ZD)
IF (LHOOK) CALL DR_HOOK('TEB_PROFILE_GRIB',1,ZHOOK_HANDLE)

END SUBROUTINE TEB_PROFILE_GRIB
!
!-------------------------------------------------------------------------------------
END SUBROUTINE PREP_TEB_GRIB
