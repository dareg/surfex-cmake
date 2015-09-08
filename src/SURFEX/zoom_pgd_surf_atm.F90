!     ###########################################################
      SUBROUTINE ZOOM_PGD_SURF_ATM (YSC, &
                                    HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE)
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
!
!
USE MODD_SURFEX_n, ONLY : SURFEX_t
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
!
TYPE(SURFEX_t), INTENT(INOUT) :: YSC
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
 CALL PGD_GRID(YSC%IOB, &
               YSC%UG, YSC%U, &
               HPROGRAM,HFILE,HFILETYPE,.TRUE.,YSC%UG%CGRID,YSC%UG%NGRID_PAR,YSC%UG%XGRID_PAR)
!
 CALL OPEN_AUX_IO_SURF(YSC%IOB, &
                       HINIFILE,HINIFILETYPE,'FULL  ')
 CALL READ_SURF(YSC%IOB, &
                HINIFILETYPE,'SEA',   YSC%U%CSEA,   IRESP)
 CALL READ_SURF(YSC%IOB, &
                HINIFILETYPE,'NATURE',YSC%U%CNATURE,IRESP)
 CALL READ_SURF(YSC%IOB, &
                HINIFILETYPE,'WATER', YSC%U%CWATER, IRESP)
 CALL READ_SURF(YSC%IOB, &
                HINIFILETYPE,'TOWN',  YSC%U%CTOWN,  IRESP)
 CALL READ_COVER_GARDEN(YSC%IOB, &
                        HINIFILETYPE,YSC%U%LGARDEN)
 CALL INIT_READ_DATA_COVER(HPROGRAM)
 CALL INI_DATA_COVER(YSC%DTCO, YSC%U)
 CALL CLOSE_AUX_IO_SURF(HINIFILE,HINIFILETYPE)
!
!-------------------------------------------------------------------------------
!
!*    3.      surface cover
!             -------------
!
 CALL ZOOM_PGD_COVER(YSC%DTCO, YSC%IOB, YSC%UG, YSC%U, &
                     HPROGRAM,HINIFILE,HINIFILETYPE,YSC%U%LECOCLIMAP)
!
!-------------------------------------------------------------------------------
!
!*    4.      Orography
!             ---------
!
 CALL ZOOM_PGD_OROGRAPHY(YSC%DTCO, YSC%IOB, &
                         YSC%UG, YSC%U, YSC%USS, &
                         HPROGRAM,YSC%U%XSEA,YSC%U%XWATER,HINIFILE,HINIFILETYPE)
!
!_______________________________________________________________________________
!
!*    5.      Additionnal fields for nature scheme
!             ------------------------------------
!
IF (YSC%U%NDIM_NATURE>0)                                 &
  CALL ZOOM_PGD_NATURE(YSC%CHI, YSC%DTCO, YSC%DTI, YSC%IOB, YSC%IG, YSC%I, YSC%UG, YSC%U, YSC%USS, &
                       HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE,YSC%U%LECOCLIMAP)  
!_______________________________________________________________________________
!
!*    6.      Additionnal fields for town scheme
!             ----------------------------------
!
IF (YSC%U%NDIM_TOWN>0)                                 &
  CALL ZOOM_PGD_TOWN(YSC%B, YSC%DGCT, YSC%DGMT, YSC%T, YSC%TGD, YSC%TGDPE, YSC%TGR, YSC%TGRPE, &
                     YSC%BOP, YSC%BDD, YSC%DTB, YSC%DTCO, YSC%DTT, YSC%IOB, YSC%UG, YSC%U, &
                     YSC%TGDO, YSC%TGDP, YSC%TG, YSC%TOP, YSC%TVG, &
                     HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE,YSC%U%LECOCLIMAP,YSC%U%LGARDEN)  
!_______________________________________________________________________________
!
!*    7.      Additionnal fields for inland water scheme
!             ------------------------------------------
!
IF (YSC%U%NDIM_WATER>0)                                 &
  CALL ZOOM_PGD_INLAND_WATER(YSC%DTCO, YSC%FG, YSC%F, YSC%UG, YSC%U, YSC%USS, YSC%WG, YSC%W, &
                             HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE,YSC%U%LECOCLIMAP)  
!_______________________________________________________________________________
!
!*    8.      Additionnal fields for sea scheme
!             ---------------------------------
!
IF (YSC%U%NDIM_SEA>0)                                 &
  CALL ZOOM_PGD_SEA(YSC%DTCO, YSC%DTS, YSC%IOB, YSC%SG, YSC%S, YSC%UG, YSC%U, &
                    HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE)  
!
!_______________________________________________________________________________
!
!*    9.      Dummy fields
!             ------------
!
YSC%DUU%NDUMMY_NBR = 0
!_______________________________________________________________________________
!
!*   10.      Chemical Emission fields
!             ------------------------
!
YSC%CHU%LCH_EMIS = .FALSE.
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_SURF_ATM',1,ZHOOK_HANDLE)
!_______________________________________________________________________________
!
END SUBROUTINE ZOOM_PGD_SURF_ATM
