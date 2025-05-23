!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE PREP_WATFLUX_BUFFER(HPROGRAM,HSURF,KLUOUT,PFIELD)
!     #################################################################################
!
!!****  *PREP_WATFLUX_BUFFER* - prepares WATFLUX field from operational BUFFER
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
USE MODE_READ_BUFFER
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE MODI_PREP_BUFFER_GRID
!
USE MODD_PREP,       ONLY : CINTERP_TYPE
USE MODD_GRID_BUFFER,  ONLY : NNI
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
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
 CHARACTER(LEN=6)              :: YINMODEL ! model from which BUFFER file originates
REAL, DIMENSION(:), POINTER :: ZFIELD   ! field read
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!*      1.     Reading of grid
!              ---------------
!
IF (LHOOK) CALL DR_HOOK('PREP_WATFLUX_BUFFER',0,ZHOOK_HANDLE)
 CALL PREP_BUFFER_GRID(KLUOUT,YINMODEL,TZTIME_BUF)

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
      CASE ('ALADIN')
        CALL READ_BUFFER_ZS(KLUOUT,YINMODEL,ZFIELD)
        ALLOCATE(PFIELD(NNI,1))
        PFIELD(:,1) = ZFIELD(:)
        DEALLOCATE(ZFIELD)
    END SELECT

!
!* 3.  Temperature profiles
!      --------------------
!
  CASE('TSWATER')
    SELECT CASE (YINMODEL)
      CASE ('ALADIN')
        CALL READ_BUFFER_T2(KLUOUT,YINMODEL,ZFIELD)
        ALLOCATE(PFIELD(NNI,1))
        PFIELD(:,1) = ZFIELD(:)
        DEALLOCATE(ZFIELD)
    END SELECT

END SELECT
!
!*      4.     Interpolation method
!              --------------------
!
CINTERP_TYPE='BUFFER'
IF (LHOOK) CALL DR_HOOK('PREP_WATFLUX_BUFFER',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------------
END SUBROUTINE PREP_WATFLUX_BUFFER
