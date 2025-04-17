!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE PREP_SNOW_GRIB(HPROGRAM,HSURF,HFILE,KLUOUT,KLAYER,PFIELD)
!     #################################################################################
!
!!****  *PREP_SNOW_GRIB* - prepares snow field from operational GRIB
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
!!      Original     01/2004
!!      C. Ardilouze 06/2013 read snow albedo and density (for Erai-land)
!!------------------------------------------------------------------
!
!
USE MODI_ABOR1_SFX  
USE MODE_READ_GRIB
!
USE MODD_TYPE_DATE_SURF
!
USE MODD_PREP_SNOW,      ONLY : NGRID_LEVEL, XGRID_SNOW, LSNOW_FRAC_ECMWF
USE MODD_PREP_ISBA,      ONLY : XRM_WM_ECMWF
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_GRID_GRIB,      ONLY : CGRIB_FILE, NNI, CINMODEL
USE MODD_SNOW_PAR,       ONLY : XANSMIN, XANSMAX, XRHOSMAX
USE MODD_CSTS,           ONLY : XTT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN)    :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=10),   INTENT(IN)   :: HSURF     ! type of field
CHARACTER(LEN=28),  INTENT(IN)    :: HFILE     ! name of file
INTEGER,            INTENT(IN)    :: KLUOUT    ! logical unit of output listing
INTEGER,            INTENT(IN)    :: KLAYER    ! Number of layer of output snow scheme
REAL,DIMENSION(:,:,:), POINTER    :: PFIELD    ! field to interpolate horizontally
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(:)  , POINTER   :: ZMASK => NULL()          ! Land mask
REAL, DIMENSION(:)  , POINTER   :: ZLFRAC => NULL()         ! Fraction of land
REAL, DIMENSION(:),   POINTER   :: ZFIELD1D => NULL()       ! field read
REAL, DIMENSION(:),   POINTER   :: ZHEAT => NULL()          ! heat in snow
REAL, DIMENSION(:),   POINTER   :: ZRHO => NULL()          ! density of snow
INTEGER                         :: JVEGTYPE       ! loop counter on vegtypes
INTEGER                         :: JLAYER         ! loop on snow fine grid
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!*      1.     Reading of grid
!              ---------------
!
IF (LHOOK) CALL DR_HOOK('PREP_SNOW_GRIB',0,ZHOOK_HANDLE)
!
IF (TRIM(HFILE).NE.CGRIB_FILE) CGRIB_FILE=""
!
 IF (CINMODEL.EQ.'ECMWF') THEN 
   CALL READ_GRIB_LAND_MASK(HFILE,KLUOUT,CINMODEL,ZMASK,PLFRAC=ZLFRAC)
 ELSE
   CALL READ_GRIB_LAND_MASK(HFILE,KLUOUT,CINMODEL,ZMASK)
 END IF
 IF (SIZE(ZMASK).NE.NNI) CALL ABOR1_SFX("PREP_SNOW_GRIB: size of LSM differs from NNI")
 IF (CINMODEL.EQ.'ECMWF'.AND.LSNOW_FRAC_ECMWF.AND.XRM_WM_ECMWF.EQ.0.) &
      CALL ABOR1_SFX("PREP_SNOW_GRIB: with LSNOW_FRAC_ECMWF, XRM_WM_ECMWF should be more than 0.0")
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading of the physical field for urban areas
!              ---------------------------------------------
!
IF (HSURF(7:8)=='RO') THEN
  ! 
  SELECT CASE(HSURF(1:3))
    CASE('DEP')
      ALLOCATE(PFIELD(NNI,1,1))
    CASE('ALB','WWW')
      ALLOCATE(PFIELD(NNI,1,1))
    CASE('HEA','RHO')
      ALLOCATE(PFIELD(NNI,1,1))
  END SELECT
  !
  PFIELD(:,:,:) = 0.
!
!-------------------------------------------------------------------------------------
!
!*      3.     Reading of the physical field for vegetated areas
!              -------------------------------------------------
!
ELSE
!
  SELECT CASE(HSURF(1:3))
!
!*      3.1    Total snow content (kg/m2)
!
  CASE('WWW')
     IF (CINMODEL.EQ.'ECMWF'.AND.LSNOW_FRAC_ECMWF) THEN
       CALL READ_GRIB_SNOW_VEG_AND_DEPTH(HFILE,KLUOUT,CINMODEL,ZMASK,PLFRAC=ZLFRAC,PSNV=ZFIELD1D) 
     ELSE
       CALL READ_GRIB_SNOW_VEG_AND_DEPTH(HFILE,KLUOUT,CINMODEL,ZMASK,PSNV=ZFIELD1D)
     END IF     
     !
     ALLOCATE(PFIELD(SIZE(ZFIELD1D),1,1))
     PFIELD(:,1,1)=ZFIELD1D(:)
     DEALLOCATE(ZFIELD1D)
!
!
!*      3.2    Total snow depth (m)
!
  CASE('DEP')
     IF (CINMODEL.EQ.'ECMWF'.AND.LSNOW_FRAC_ECMWF) THEN
       CALL READ_GRIB_SNOW_VEG_AND_DEPTH(HFILE,KLUOUT,CINMODEL,ZMASK,PLFRAC=ZLFRAC,PSNVD=ZFIELD1D)
     ELSE
       CALL READ_GRIB_SNOW_VEG_AND_DEPTH(HFILE,KLUOUT,CINMODEL,ZMASK,PSNVD=ZFIELD1D)       
     END IF
     !
     ALLOCATE(PFIELD(SIZE(ZFIELD1D),1,1))
     PFIELD(:,1,1)=ZFIELD1D(:)
     DEALLOCATE(ZFIELD1D)
!
!
!*      3.3    Profile of heat in the snow
!
  CASE('HEA')
     !* read temperature
     IF (CINMODEL.EQ.'ECMWF'.AND.LSNOW_FRAC_ECMWF) THEN
        CALL READ_GRIB_TS(HFILE,KLUOUT,CINMODEL,ZMASK,ZFIELD1D,PLFRAC=ZLFRAC)
     ELSE
        CALL READ_GRIB_TS(HFILE,KLUOUT,CINMODEL,ZMASK,ZFIELD1D)
     END IF
     WHERE (ZFIELD1D/=XUNDEF) ZFIELD1D(:) = MIN(ZFIELD1D,XTT)
     !
     ALLOCATE(PFIELD(SIZE(ZFIELD1D),1,1))
     PFIELD(:,1,1)=ZFIELD1D(:)
     DEALLOCATE(ZFIELD1D)
!
!*      3.4    Albedo
!
  CASE('ALB')    
    IF (CINMODEL.EQ.'ECMWF'.AND.LSNOW_FRAC_ECMWF) THEN
      CALL READ_GRIB_SNOW_ALB(HFILE,KLUOUT,CINMODEL,ZMASK,ZFIELD1D,PLFRAC=ZLFRAC)
    ELSE
      CALL READ_GRIB_SNOW_ALB(HFILE,KLUOUT,CINMODEL,ZMASK,ZFIELD1D)
    END IF
    ALLOCATE(PFIELD(SIZE(ZFIELD1D),1,1))
    PFIELD(:,1,1)=ZFIELD1D(:)
    DEALLOCATE(ZFIELD1D)
!
!*      3.5    Density
!
  CASE('RHO')    
    IF (CINMODEL.EQ.'ECMWF'.AND.LSNOW_FRAC_ECMWF) THEN
      CALL READ_GRIB_SNOW_DEN(HFILE,KLUOUT,CINMODEL,ZMASK,ZFIELD1D,PLFRAC=ZLFRAC)
    ELSE
      CALL READ_GRIB_SNOW_DEN(HFILE,KLUOUT,CINMODEL,ZMASK,ZFIELD1D)
    END IF
    ALLOCATE(PFIELD(SIZE(ZFIELD1D),1,1))
    PFIELD(:,1,1)=ZFIELD1D(:)
    DEALLOCATE(ZFIELD1D)
!
!*      3.6    SG1: initial grain is partially rounded
!
  CASE('SG1')
    ALLOCATE(PFIELD(NNI,1,1))
    IF (CINMODEL.EQ.'ECMWF'.AND.LSNOW_FRAC_ECMWF) THEN
       WHERE(ZLFRAC(:).LT.XRM_WM_ECMWF)
          PFIELD(:,1,1) = XUNDEF
       ELSEWHERE
          PFIELD(:,1,1) = -20
       END WHERE
    ELSE
       WHERE(ZMASK(:).NE.1.)
          PFIELD(:,1,1) = XUNDEF
       ELSEWHERE
          PFIELD(:,1,1) = -20
       END WHERE
    END IF
!
!*      3.7    SG2: initial grain is partially rounded
!
  CASE('SG2')
    ALLOCATE(PFIELD(NNI,1,1))
    IF (CINMODEL.EQ.'ECMWF'.AND.LSNOW_FRAC_ECMWF) THEN
       WHERE(ZLFRAC(:).LT.XRM_WM_ECMWF)
          PFIELD(:,1,1) = XUNDEF
       ELSEWHERE
          PFIELD(:,1,1) = 80
       END WHERE
    ELSE
       WHERE(ZMASK(:).NE.1.)
          PFIELD(:,1,1) = XUNDEF
       ELSEWHERE
          PFIELD(:,1,1) = 80
       END WHERE
    END IF
!
!*      3.8    AGE: snow is 3-days old
!
  CASE('AGE')
    ALLOCATE(PFIELD(NNI,1,1))
    IF (CINMODEL.EQ.'ECMWF'.AND.LSNOW_FRAC_ECMWF) THEN
       WHERE(ZLFRAC(:).LT.XRM_WM_ECMWF)
          PFIELD(:,1,1) = XUNDEF
       ELSEWHERE
          PFIELD(:,1,1) = 3
       END WHERE
    ELSE
       WHERE(ZMASK(:).NE.1.)
          PFIELD(:,1,1) = XUNDEF
       ELSEWHERE
          PFIELD(:,1,1) = 3
       END WHERE
    END IF
!
!*      3.9    HIS: 0 by default
!
  CASE('HIS')
    ALLOCATE(PFIELD(NNI,1,1))
    IF (CINMODEL.EQ.'ECMWF'.AND.LSNOW_FRAC_ECMWF) THEN
       WHERE(ZLFRAC(:).LT.XRM_WM_ECMWF)
          PFIELD(:,1,1) = XUNDEF
       ELSEWHERE
          PFIELD(:,1,1) = 0
       END WHERE
    ELSE
       WHERE(ZMASK(:).NE.1.)
          PFIELD(:,1,1) = XUNDEF
       ELSEWHERE
          PFIELD(:,1,1) = 0
       END WHERE
    END IF
!
  END SELECT
!
END IF
!
DEALLOCATE(ZMASK)
IF (CINMODEL.EQ.'ECMWF') THEN
   DEALLOCATE(ZLFRAC)
END IF
!
!-------------------------------------------------------------------------------------
!
!*      4.     Interpolation method
!              --------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_SNOW_GRIB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
END SUBROUTINE PREP_SNOW_GRIB
