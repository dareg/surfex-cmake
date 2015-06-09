!     ###########################################################
      SUBROUTINE ZOOM_PGD_SURF_ATM(HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE)
!     ###########################################################

!!
!!    PURPOSE
!!    -------
!!   This program prepares the physiographic data fields.
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson                   Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     13/10/03
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
!
USE MODD_SURF_ATM_SSO_n, ONLY : USS => SURF_ATM_SSO
!
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
USE MODD_CH_SURF_n, ONLY : CHU => CH_SURF
USE MODD_DUMMY_SURF_FIELDS_n, ONLY : DUU => DUMMY_SURF_FIELDS
!
USE MODI_INI_CSTS
USE MODI_READ_NAM_WRITE_COVER_TEX
USE MODI_PGD_GRID
USE MODI_OPEN_AUX_IO_SURF
USE MODI_CLOSE_AUX_IO_SURF
USE MODI_READ_SURF
USE MODI_ZOOM_PGD_COVER
USE MODI_ZOOM_PGD_OROGRAPHY
USE MODI_INIT_READ_DATA_COVER
USE MODI_INI_DATA_COVER
USE MODI_SURF_VERSION
USE MODI_ZOOM_PGD_INLAND_WATER
USE MODI_ZOOM_PGD_NATURE
USE MODI_ZOOM_PGD_SEA
USE MODI_ZOOM_PGD_TOWN
USE MODI_READ_COVER_GARDEN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
 CHARACTER(LEN=6),     INTENT(IN)  :: HPROGRAM    ! program calling
 CHARACTER(LEN=28),    INTENT(IN)  :: HINIFILE    ! input atmospheric file name
 CHARACTER(LEN=6),     INTENT(IN)  :: HINIFILETYPE! input atmospheric file type
 CHARACTER(LEN=28),    INTENT(IN)  :: HFILE       ! output atmospheric file name
 CHARACTER(LEN=6),     INTENT(IN)  :: HFILETYPE   ! output atmospheric file type
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: IRESP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
!
!*    1.      Set default constant values 
!             ---------------------------
!
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_SURF_ATM',0,ZHOOK_HANDLE)
 CALL SURF_VERSION
!
 CALL INI_CSTS
!
 CALL READ_NAM_WRITE_COVER_TEX(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*    2.      Initialisation of output grid and schemes
!             -----------------------------------------
!
 CALL PGD_GRID(UG, U, &
               HPROGRAM,HFILE,HFILETYPE,.TRUE.,UG%CGRID,UG%NGRID_PAR,UG%XGRID_PAR)
!
 CALL OPEN_AUX_IO_SURF(HINIFILE,HINIFILETYPE,'FULL  ')
 CALL READ_SURF(HINIFILETYPE,'SEA',   U%CSEA,   IRESP)
 CALL READ_SURF(HINIFILETYPE,'NATURE',U%CNATURE,IRESP)
 CALL READ_SURF(HINIFILETYPE,'WATER', U%CWATER, IRESP)
 CALL READ_SURF(HINIFILETYPE,'TOWN',  U%CTOWN,  IRESP)
 CALL READ_COVER_GARDEN(HINIFILETYPE,U%LGARDEN)
 CALL INIT_READ_DATA_COVER(HPROGRAM)
 CALL INI_DATA_COVER(DTCO, U)
 CALL CLOSE_AUX_IO_SURF(HINIFILE,HINIFILETYPE)
!
!-------------------------------------------------------------------------------
!
!*    3.      surface cover
!             -------------
!
 CALL ZOOM_PGD_COVER(HPROGRAM,HINIFILE,HINIFILETYPE,U%LECOCLIMAP)
!
!-------------------------------------------------------------------------------
!
!*    4.      Orography
!             ---------
!
 CALL ZOOM_PGD_OROGRAPHY(UG, U, USS, &
                         HPROGRAM,U%XSEA,U%XWATER,HINIFILE,HINIFILETYPE)
!
!_______________________________________________________________________________
!
!*    5.      Additionnal fields for nature scheme
!             ------------------------------------
!
IF (U%NDIM_NATURE>0)                                 &
  CALL ZOOM_PGD_NATURE(HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE,U%LECOCLIMAP)  
!_______________________________________________________________________________
!
!*    6.      Additionnal fields for town scheme
!             ----------------------------------
!
IF (U%NDIM_TOWN>0)                                 &
  CALL ZOOM_PGD_TOWN(HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE,U%LECOCLIMAP,U%LGARDEN)  
!_______________________________________________________________________________
!
!*    7.      Additionnal fields for inland water scheme
!             ------------------------------------------
!
IF (U%NDIM_WATER>0)                                 &
  CALL ZOOM_PGD_INLAND_WATER(HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE,U%LECOCLIMAP)  
!_______________________________________________________________________________
!
!*    8.      Additionnal fields for sea scheme
!             ---------------------------------
!
IF (U%NDIM_SEA>0)                                 &
  CALL ZOOM_PGD_SEA(HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE)  
!
!_______________________________________________________________________________
!
!*    9.      Dummy fields
!             ------------
!
DUU%NDUMMY_NBR = 0
!_______________________________________________________________________________
!
!*   10.      Chemical Emission fields
!             ------------------------
!
CHU%LCH_EMIS = .FALSE.
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_SURF_ATM',1,ZHOOK_HANDLE)
!_______________________________________________________________________________
!
END SUBROUTINE ZOOM_PGD_SURF_ATM
