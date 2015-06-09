!     ###########################################################
      SUBROUTINE PGD_SURF_ATM(HPROGRAM,HFILE,HFILETYPE,OZS)
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
!!      A. Lemonsu      05/2009         Ajout de la clef LGARDEN pour TEB
!!      J. Escobar      11/2013         Add USE MODI_READ_NAM_PGD_CHEMISTRY
!!      B. Decharme     02/2014         Add LRM_RIVER
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
!
USE MODD_SURF_CONF,       ONLY : CPROGNAME
USE MODD_PGD_GRID,        ONLY : LLATLONMASK
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_CH_SURF_n, ONLY : CHU => CH_SURF
!
USE MODI_GET_LUOUT
USE MODI_READ_PGD_ARRANGE_COVER
USE MODI_READ_PGD_COVER_GARDEN
USE MODI_INI_DATA_COVER
USE MODI_READ_PGD_SCHEMES
USE MODI_READ_NAM_PGD_CHEMISTRY
USE MODI_READ_NAM_WRITE_COVER_TEX
USE MODI_WRITE_COVER_TEX_START
USE MODI_WRITE_COVER_TEX_COVER
USE MODI_LATLON_GRID
USE MODI_PUT_PGD_GRID
USE MODI_LATLONMASK
USE MODI_PGD_FRAC
USE MODI_PGD_COVER
USE MODI_PGD_OROGRAPHY
USE MODI_PGD_NATURE
USE MODI_PGD_TOWN
USE MODI_PGD_INLAND_WATER
USE MODI_PGD_SEA
USE MODI_PGD_DUMMY
USE MODI_PGD_CHEMISTRY
USE MODI_PGD_CHEMISTRY_SNAP
USE MODI_WRITE_COVER_TEX_END
USE MODI_INIT_READ_DATA_COVER
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
CHARACTER(LEN=6),     INTENT(IN)  :: HPROGRAM ! program calling
CHARACTER(LEN=28),    INTENT(IN)  :: HFILE    ! atmospheric file name
CHARACTER(LEN=6),     INTENT(IN)  :: HFILETYPE! atmospheric file type
LOGICAL,              INTENT(IN)  :: OZS      ! .true. if orography is imposed by atm. model
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
LOGICAL :: LRM_RIVER   !delete inland river coverage. Default is false
!
INTEGER :: ILUOUT ! logical unit of output listing file
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PGD_SURF_ATM',0,ZHOOK_HANDLE)
!
LRM_RIVER = .FALSE.
!
CPROGNAME=HPROGRAM
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!*    1.      Set default constant values 
!             ---------------------------
!
 CALL READ_PGD_ARRANGE_COVER(HPROGRAM,U%LWATER_TO_NATURE,U%LTOWN_TO_ROCK)
!
 CALL READ_PGD_COVER_GARDEN(HPROGRAM,U%LGARDEN)
!
 CALL INIT_READ_DATA_COVER(HPROGRAM)
!
 CALL INI_DATA_COVER(DTCO, U)
!
!*    1.2     surface schemes
 CALL READ_PGD_SCHEMES(HPROGRAM,U%CNATURE,U%CSEA,U%CTOWN,U%CWATER)
!
!*    1.3     prints all parameters in a Latex file
 CALL READ_NAM_WRITE_COVER_TEX(HPROGRAM)
!
 CALL WRITE_COVER_TEX_START(HPROGRAM)
 CALL WRITE_COVER_TEX_COVER
!-------------------------------------------------------------------------------
!
!*    2.      Grid
!             ----
!
ALLOCATE(UG%XLAT(U%NSIZE_FULL))
ALLOCATE(UG%XLON(U%NSIZE_FULL))
ALLOCATE(UG%XMESH_SIZE(U%NSIZE_FULL))
ALLOCATE(UG%XJPDIR(U%NSIZE_FULL))
 CALL LATLON_GRID(UG%CGRID,UG%NGRID_PAR,U%NSIZE_FULL,ILUOUT,UG%XGRID_PAR,UG%XLAT,UG%XLON,UG%XMESH_SIZE,UG%XJPDIR)
!
!
!*    2.3     Stores the grid in the module MODD_PGD_GRID
!
 CALL PUT_PGD_GRID(UG%CGRID,U%NSIZE_FULL,UG%NGRID_PAR,UG%XGRID_PAR)
!
!*    2.4     mask to limit the number of input data to read
 CALL LATLONMASK      (UG%CGRID,UG%NGRID_PAR,UG%XGRID_PAR,LLATLONMASK)
!
!-------------------------------------------------------------------------------
!
!*    3.      surface cover
!             -------------
!
 CALL PGD_FRAC(HPROGRAM,U%LECOCLIMAP)
IF (U%LECOCLIMAP) CALL PGD_COVER(HPROGRAM,LRM_RIVER)
!
!-------------------------------------------------------------------------------
!
!*    4.      Orography
!             ---------
!
 CALL PGD_OROGRAPHY(HPROGRAM,U%XSEA,U%XWATER,HFILE,HFILETYPE,OZS)
!
!_______________________________________________________________________________
!
!*    5.      Additionnal fields for nature scheme
!             ------------------------------------
!
IF (U%NDIM_NATURE>0) CALL PGD_NATURE(HPROGRAM,U%LECOCLIMAP)  
!_______________________________________________________________________________
!
!*    6.      Additionnal fields for town scheme
!             ----------------------------------
!
IF (U%NDIM_TOWN>0) CALL PGD_TOWN(HPROGRAM,U%LECOCLIMAP,U%LGARDEN)  
!_______________________________________________________________________________
!
!*    7.      Additionnal fields for inland water scheme
!             ------------------------------------------
!
IF (U%NDIM_WATER>0) CALL PGD_INLAND_WATER(HPROGRAM,U%LECOCLIMAP,LRM_RIVER)   
!_______________________________________________________________________________
!
!*    8.      Additionnal fields for sea scheme
!             ---------------------------------
!
IF (U%NDIM_SEA>0) CALL PGD_SEA(HPROGRAM)  
!
!_______________________________________________________________________________
!
!*    9.      Dummy fields
!             ------------
!
 CALL PGD_DUMMY(HPROGRAM)
!_______________________________________________________________________________
!
!*   10.      Chemical Emission fields
!             ------------------------
!
 CALL READ_NAM_PGD_CHEMISTRY(HPROGRAM,CHU%CCH_EMIS)
IF (CHU%CCH_EMIS=='SNAP') THEN
  CALL PGD_CHEMISTRY_SNAP(HPROGRAM,CHU%LCH_EMIS)
ELSE IF (CHU%CCH_EMIS=='AGGR') THEN
  CALL PGD_CHEMISTRY(HPROGRAM,CHU%LCH_EMIS)
ENDIF
!_______________________________________________________________________________
!
!*   11.     Writing in cover latex file
!            ---------------------------
!
 CALL WRITE_COVER_TEX_END(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('PGD_SURF_ATM',1,ZHOOK_HANDLE)
!_______________________________________________________________________________
!
END SUBROUTINE PGD_SURF_ATM
