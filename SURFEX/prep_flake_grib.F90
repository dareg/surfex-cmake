!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE PREP_FLAKE_GRIB(HPROGRAM,HSURF,HFILE,KLUOUT,PFIELD)
!     #################################################################################
!
!!****  *PREP_FLAKE_GRIB* - prepares FLAKE field from operational GRIB
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
!!     S. Malardel
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!------------------------------------------------------------------
!salgado - Is equal to prep_watflux_grib.f90
!        - The other variables are not read from grib files
!        - may be changed in future, namely if AROME will use FLAKE
!
USE MODE_READ_GRIB
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE MODD_GRID_GRIB,  ONLY : CGRIB_FILE
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
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
 CHARACTER(LEN=6)              :: YINMODEL ! model from which GRIB file originates
REAL, DIMENSION(:)  , POINTER   :: ZMASK => NULL()          ! Land mask
REAL, DIMENSION(:), POINTER :: ZFIELD => NULL()   ! field read
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!*      1.     Reading of grid
!              ---------------
!
IF (LHOOK) CALL DR_HOOK('PREP_FLAKE_GRIB',0,ZHOOK_HANDLE)
!
IF (TRIM(HFILE).NE.CGRIB_FILE) CGRIB_FILE=""
!
 CALL READ_GRIB_LAND_MASK(HFILE,KLUOUT,YINMODEL,ZMASK)
!
!*      2.     Reading of field
!              ----------------
!
!--------------------
SELECT CASE(HSURF)
!--------------------
!
!* 1.  Orography
!      ---------
!
  CASE('ZS     ')
    SELECT CASE (YINMODEL)
      CASE ('ECMWF ','ARPEGE','ALADIN','MOCAGE','HIRLAM','RACMO ')
        CALL READ_GRIB_ZS_LAND(HFILE,KLUOUT,YINMODEL,ZMASK,ZFIELD)
        ALLOCATE(PFIELD(SIZE(ZFIELD),1))
        PFIELD(:,1) = ZFIELD(:)
        DEALLOCATE(ZFIELD)
      CASE DEFAULT
        CALL ABOR1_SFX('PREP_FLAKE_GRIB (ZS): UNSUPPORTED MODEL ' // TRIM(YINMODEL))
    END SELECT

!
!* 3.  Temperature profiles
!      --------------------
!
  CASE('TS     ')
    SELECT CASE (YINMODEL)
      CASE ('ECMWF ','ARPEGE','ALADIN','MOCAGE','HIRLAM','RACMO ')
        CALL READ_GRIB_TSWATER(HFILE,KLUOUT,YINMODEL,ZMASK,ZFIELD)
        ALLOCATE(PFIELD(SIZE(ZFIELD),1))
        PFIELD(:,1) = ZFIELD(:)
        DEALLOCATE(ZFIELD)
      CASE DEFAULT
        CALL ABOR1_SFX('PREP_FLAKE_GRIB (TS): UNSUPPORTED MODEL ' // TRIM(YINMODEL))
    END SELECT

END SELECT
!
DEALLOCATE(ZMASK)
!
!*      4.     Interpolation method
!              --------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_FLAKE_GRIB',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------------
END SUBROUTINE PREP_FLAKE_GRIB
