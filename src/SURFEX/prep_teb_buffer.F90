!     #########
SUBROUTINE PREP_TEB_BUFFER(HPROGRAM,HSURF,KLUOUT,PFIELD)
!     #################################################################################
!
!!****  *PREP_TEB_BUFFER* - prepares TEB field from operational BUFFER
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
!!      Original    03/2005
!!------------------------------------------------------------------
!

!
USE MODD_TYPE_DATE_SURF
!
USE MODI_PREP_BUFFER_GRID
USE MODE_READ_BUFFER
USE MODI_INTERP_GRID
!
USE MODD_PREP,       ONLY : CINTERP_TYPE
USE MODD_GRID_BUFFER,  ONLY : NNI
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
INTEGER,            INTENT(IN)  :: KLUOUT    ! logical unit of output listing
REAL,DIMENSION(:,:), POINTER    :: PFIELD    ! field to interpolate horizontally
!
!*      0.2    declarations of local variables
!
TYPE (DATE_TIME)                :: TZTIME_BUF    ! current date and time
CHARACTER(LEN=6)                :: YINMODEL ! model from which BUFFER originates
REAL, DIMENSION(:),   POINTER   :: ZFIELD1D ! 1D field read
REAL, DIMENSION(:,:), POINTER   :: ZFIELD   ! field read
REAL, DIMENSION(:,:), POINTER   :: ZD             ! depth of field in the soil
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!*      1.     Reading of grid
!              ---------------
!
IF (LHOOK) CALL DR_HOOK('PREP_TEB_BUFFER',0,ZHOOK_HANDLE)
CALL PREP_BUFFER_GRID(KLUOUT,YINMODEL,TZTIME_BUF)
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
      CASE ('ALADIN')
        CALL READ_BUFFER_ZS(KLUOUT,YINMODEL,ZFIELD1D)
        ALLOCATE(PFIELD(NNI,1))
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
       CASE('ALADIN')
         CALL READ_BUFFER_TG(KLUOUT,YINMODEL,ZFIELD,ZD)
     END SELECT
     !* if deep road temperature is prescribed
     IF (XTI_ROAD/=XUNDEF) THEN
       ZFIELD(:,2:) = XTI_ROAD 
     END IF
     CALL TEB_PROFILE_BUFFER(XGRID_ROAD)

!*      4.     Profile of temperatures in walls
!              --------------------------------

  CASE('T_WALL')
     CALL READ_BUFFER_T_TEB(KLUOUT,YINMODEL,XTI_BLD,ZFIELD,ZD)
     CALL TEB_PROFILE_BUFFER(XGRID_WALL)


!*      5.     Profile of temperatures in roofs
!              --------------------------------

  CASE('T_ROOF')    
     CALL READ_BUFFER_T_TEB(KLUOUT,YINMODEL,XTI_BLD,ZFIELD,ZD)
     CALL TEB_PROFILE_BUFFER(XGRID_ROOF)

!
!*      6.     Canyon air temperature
!              ----------------------
!
  CASE('T_CAN  ')
    SELECT CASE (YINMODEL)
      CASE ('ALADIN')
        CALL READ_BUFFER_T2(KLUOUT,YINMODEL,ZFIELD1D)
        ALLOCATE(PFIELD(NNI,1))
        PFIELD(:,1) = ZFIELD1D(:)
        DEALLOCATE(ZFIELD1D)
    END SELECT
!
!*      7.      Canyon air humidity
!               -------------------
!
  CASE('Q_CAN  ')
    SELECT CASE (YINMODEL)
      CASE ('ALADIN')
        ALLOCATE(PFIELD(NNI,1))
        PFIELD(:,1) = 0.
    END SELECT

!
!*      9.     Deep road temperature
!              ---------------------

  CASE('TI_ROAD')    
     IF (XTI_ROAD==XUNDEF) THEN
       CALL READ_BUFFER_T2(KLUOUT,YINMODEL,ZFIELD1D)
       ALLOCATE(PFIELD(NNI,1))
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
!*      4.     Interpolation method
!              --------------------
!
CINTERP_TYPE='BUFFER'
!
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_TEB_BUFFER',1,ZHOOK_HANDLE)
CONTAINS
!
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
SUBROUTINE TEB_PROFILE_BUFFER(PGRID)
!-------------------------------------------------------------------------------------
!
REAL, DIMENSION(:),   INTENT(IN)  :: PGRID  ! destination grid
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!* interpolation on fine vertical grid
IF (LHOOK) CALL DR_HOOK('TEB_PROFILE_BUFFER',0,ZHOOK_HANDLE)
ALLOCATE(PFIELD(SIZE(ZFIELD,1),SIZE(PGRID)))
CALL INTERP_GRID(ZD,ZFIELD,PGRID,PFIELD)
!
!* end
DEALLOCATE(ZFIELD)
DEALLOCATE(ZD)
IF (LHOOK) CALL DR_HOOK('TEB_PROFILE_BUFFER',1,ZHOOK_HANDLE)

END SUBROUTINE TEB_PROFILE_BUFFER
!
!-------------------------------------------------------------------------------------
END SUBROUTINE PREP_TEB_BUFFER
