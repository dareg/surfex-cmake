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
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
USE MODD_DATA_ISBA_n, ONLY : DTI => DATA_ISBA
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
USE MODD_ISBA_n, ONLY : I => ISBA
!
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_BLD_DESCRIPTION_n, ONLY : BDD => BLD_DESC
USE MODD_DATA_BEM_n, ONLY : DTB => DATA_BEM
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_GRID_n, ONLY : TG => TEB_GRID
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
!
USE MODD_DATA_SEAFLUX_n, ONLY : DTS => DATA_SEAFLUX
USE MODD_SEAFLUX_GRID_n, ONLY : SG => SEAFLUX_GRID
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
!
USE MODD_FLAKE_GRID_n, ONLY : FG => FLAKE_GRID
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_WATFLUX_GRID_n, ONLY : WG => WATFLUX_GRID
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
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
 CALL PGD_GRID(IOB, &
               UG, U, &
               HPROGRAM,HFILE,HFILETYPE,.TRUE.,UG%CGRID,UG%NGRID_PAR,UG%XGRID_PAR)
!
 CALL OPEN_AUX_IO_SURF(IOB, &
                       HINIFILE,HINIFILETYPE,'FULL  ')
 CALL READ_SURF(IOB, &
                HINIFILETYPE,'SEA',   U%CSEA,   IRESP)
 CALL READ_SURF(IOB, &
                HINIFILETYPE,'NATURE',U%CNATURE,IRESP)
 CALL READ_SURF(IOB, &
                HINIFILETYPE,'WATER', U%CWATER, IRESP)
 CALL READ_SURF(IOB, &
                HINIFILETYPE,'TOWN',  U%CTOWN,  IRESP)
 CALL READ_COVER_GARDEN(IOB, &
                        HINIFILETYPE,U%LGARDEN)
 CALL INIT_READ_DATA_COVER(HPROGRAM)
 CALL INI_DATA_COVER(DTCO, U)
 CALL CLOSE_AUX_IO_SURF(HINIFILE,HINIFILETYPE)
!
!-------------------------------------------------------------------------------
!
!*    3.      surface cover
!             -------------
!
 CALL ZOOM_PGD_COVER(DTCO, IOB, UG, U, &
                     HPROGRAM,HINIFILE,HINIFILETYPE,U%LECOCLIMAP)
!
!-------------------------------------------------------------------------------
!
!*    4.      Orography
!             ---------
!
 CALL ZOOM_PGD_OROGRAPHY(DTCO, IOB, &
                         UG, U, USS, &
                         HPROGRAM,U%XSEA,U%XWATER,HINIFILE,HINIFILETYPE)
!
!_______________________________________________________________________________
!
!*    5.      Additionnal fields for nature scheme
!             ------------------------------------
!
IF (U%NDIM_NATURE>0)                                 &
  CALL ZOOM_PGD_NATURE(CHI, DTCO, DTI, IOB, IG, I, UG, U, USS, &
                       HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE,U%LECOCLIMAP)  
!_______________________________________________________________________________
!
!*    6.      Additionnal fields for town scheme
!             ----------------------------------
!
IF (U%NDIM_TOWN>0)                                 &
  CALL ZOOM_PGD_TOWN(BOP, BDD, DTB, DTCO, DTT, IOB, UG, U, TGDO, TGDP, TG, &
                                TOP, TVG, &
                     HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE,U%LECOCLIMAP,U%LGARDEN)  
!_______________________________________________________________________________
!
!*    7.      Additionnal fields for inland water scheme
!             ------------------------------------------
!
IF (U%NDIM_WATER>0)                                 &
  CALL ZOOM_PGD_INLAND_WATER(DTCO, FG, F, UG, U, USS, WG, W, &
                             HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE,U%LECOCLIMAP)  
!_______________________________________________________________________________
!
!*    8.      Additionnal fields for sea scheme
!             ---------------------------------
!
IF (U%NDIM_SEA>0)                                 &
  CALL ZOOM_PGD_SEA(DTCO, DTS, IOB, SG, S, UG, U, &
                    HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE)  
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
