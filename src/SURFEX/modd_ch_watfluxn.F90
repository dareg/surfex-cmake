!     #####################
      MODULE MODD_CH_WATFLUX_n
!     ######################
!
!!
!!    PURPOSE
!!    -------
!     
!   
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!
!!    AUTHOR
!!    ------
!!  P. Tulet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!  16/07/03 (P. Tulet)  restructured for externalization
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE CH_WATFLUX_t
!
  CHARACTER(LEN=6)              :: CCH_DRY_DEP              ! deposition scheme
  REAL, DIMENSION(:,:), POINTER :: XDEP                     ! final dry deposition  
                                                            ! velocity  for lakes
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSV            ! name of the scalar var.                                                            
  INTEGER    :: NSV_CHSBEG, NSV_CHSEND                      ! chemical begin and ending
                                                            ! index of the HSV/CSV array
  INTEGER    :: NBEQ                                        ! number of chemical species
                                                            ! in the surface scheme

  INTEGER    :: NSV_DSTBEG, NSV_DSTEND                      ! dust begin and ending
  INTEGER    :: NSV_SLTBEG, NSV_SLTEND                      ! sea salt begin and ending
  INTEGER    :: NSV_AERBEG, NSV_AEREND                      ! aerosol begin and ending
  INTEGER    :: NDSTEQ                                      ! number of dust species
  INTEGER    :: NSLTEQ                                      ! number of sea salt species
  INTEGER    :: NAEREQ                                      ! number of aerosol species

  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES      ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES
!
END TYPE CH_WATFLUX_t
!
CONTAINS
!
SUBROUTINE CH_WATFLUX_INIT(YCH_WATFLUX)
TYPE(CH_WATFLUX_t), INTENT(INOUT) :: YCH_WATFLUX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_WATFLUX_N:CH_WATFLUX_INIT",0,ZHOOK_HANDLE)
NULLIFY(YCH_WATFLUX%XDEP)
NULLIFY(YCH_WATFLUX%CSV)
NULLIFY(YCH_WATFLUX%CCH_NAMES)
NULLIFY(YCH_WATFLUX%CAER_NAMES)
NULLIFY(YCH_WATFLUX%CDSTNAMES)
NULLIFY(YCH_WATFLUX%CSLTNAMES)
YCH_WATFLUX%CCH_DRY_DEP=' '
YCH_WATFLUX%NSV_CHSBEG=0
YCH_WATFLUX%NSV_CHSEND=0
YCH_WATFLUX%NSV_DSTBEG=0
YCH_WATFLUX%NSV_DSTEND=0
YCH_WATFLUX%NBEQ=0
YCH_WATFLUX%NSV_SLTBEG=0
YCH_WATFLUX%NSV_SLTEND=0
YCH_WATFLUX%NSV_AERBEG=0
YCH_WATFLUX%NSV_AEREND=0
YCH_WATFLUX%NDSTEQ=0
YCH_WATFLUX%NSLTEQ=0
YCH_WATFLUX%NAEREQ=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_WATFLUX_N:CH_WATFLUX_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_WATFLUX_INIT


END MODULE MODD_CH_WATFLUX_n
