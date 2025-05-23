!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE PREP_WATFLUX_EXTERN (GCP,HPROGRAM,HSURF,HFILE,HFILETYPE,HFILEPGD,HFILEPGDTYPE,KLUOUT,PFIELD)
!     #################################################################################
!
!
USE MODD_GRID_CONF_PROJ_n, ONLY : GRID_CONF_PROJ_t
!
USE MODD_PREP,       ONLY : CINGRID_TYPE, CINTERP_TYPE
USE MODD_SURFEX_MPI, ONLY : NRANK,NPIO
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE MODI_PREP_GRID_EXTERN
USE MODI_READ_SURF
USE MODI_OPEN_AUX_IO_SURF
USE MODI_CLOSE_AUX_IO_SURF
USE MODI_ABOR1_SFX
USE MODI_GET_LUOUT
!
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(GRID_CONF_PROJ_t),INTENT(INOUT) :: GCP
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=7),   INTENT(IN)  :: HSURF     ! type of field
 CHARACTER(LEN=28),  INTENT(IN)  :: HFILE     ! name of file
 CHARACTER(LEN=6),   INTENT(IN)  :: HFILETYPE ! type of input file
 CHARACTER(LEN=28),  INTENT(IN)  :: HFILEPGD     ! name of file
 CHARACTER(LEN=6),   INTENT(IN)  :: HFILEPGDTYPE ! type of input file
INTEGER,            INTENT(IN)  :: KLUOUT    ! logical unit of output listing
REAL,DIMENSION(:,:), POINTER    :: PFIELD    ! field to interpolate horizontally
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(:), ALLOCATABLE :: ZMASK
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
INTEGER           :: IRESP          ! reading return code
INTEGER           :: ILUOUT
INTEGER           :: IDIM_WATER
!
INTEGER           :: IVERSION       ! total 1D dimension
INTEGER           :: INI            ! total 1D dimension
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!*      1.     Preparation of IO for reading in the file
!              -----------------------------------------
!
!* Note that all points are read, even those without physical meaning.
!  These points will not be used during the horizontal interpolation step.
!  Their value must be defined as XUNDEF.
!
IF (LHOOK) CALL DR_HOOK('PREP_WATFLUX_EXTERN',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading of grid
!              ---------------
!
CALL OPEN_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE,'FULL  ')
CALL PREP_GRID_EXTERN(GCP,HFILEPGDTYPE,KLUOUT,CINGRID_TYPE,CINTERP_TYPE,INI)
!
CALL READ_SURF(HFILEPGDTYPE,'DIM_WATER',IDIM_WATER,IRESP,HDIR='-')
!
YRECFM='VERSION'
 CALL READ_SURF(HFILEPGDTYPE,YRECFM,IVERSION,IRESP)
!
IF (NRANK/=NPIO) INI = 0
!
ALLOCATE(ZMASK(INI))
IF (IVERSION>=7) THEN
  YRECFM='FRAC_WATER'
  CALL READ_SURF(HFILEPGDTYPE,YRECFM,ZMASK,IRESP,HDIR='A')       
ELSE
  ZMASK(:) = 1.
ENDIF  
!
IF (IDIM_WATER==0) THEN
  CALL GET_LUOUT(HPROGRAM,ILUOUT)
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) 'No inland water data available in input file ',HFILE
  WRITE(ILUOUT,*) 'Please change your input file '
  WRITE(ILUOUT,*) '             or '
  WRITE(ILUOUT,*) 'specify inland water temperature XTS_WATER_UNIF'
  CALL ABOR1_SFX('PREP_WATFLUX_EXTERN: No inland water data available in input file')
END IF
!
IF (NRANK/=NPIO) INI=0
!
!---------------------------------------------------------------------------------------
SELECT CASE(HSURF)
!---------------------------------------------------------------------------------------
!
!*     3.      Orography
!              ---------
!
  CASE('ZS     ')
    ALLOCATE(PFIELD(INI,1))
    YRECFM='ZS'
    CALL READ_SURF(HFILEPGDTYPE,YRECFM,PFIELD(:,1),IRESP,HDIR='E')
    CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)
!
!*      4.  Sea surface temperature
!           -----------------------
!
  CASE('TSWATER')
    ALLOCATE(PFIELD(INI,1))
    YRECFM='TS_WATER'
    CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)
    CALL OPEN_AUX_IO_SURF(HFILE,HFILETYPE,'WATER ')
    CALL READ_SURF(HFILETYPE,YRECFM,PFIELD(:,1),IRESP,HDIR='E')
    CALL CLOSE_AUX_IO_SURF(HFILE,HFILETYPE)
    WHERE (ZMASK(:)==0.) PFIELD(:,1) = XUNDEF    
!
!---------------------------------------------------------------------------------------
END SELECT
!-------------------------------------------------------------------------------------
!
DEALLOCATE(ZMASK)
!
!*      6.     End of IO
!              ---------
!
IF (LHOOK) CALL DR_HOOK('PREP_WATFLUX_EXTERN',1,ZHOOK_HANDLE)
!
!---------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_WATFLUX_EXTERN
