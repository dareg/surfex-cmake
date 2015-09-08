!     #####################
      MODULE MODD_CH_SEAFLUX_n
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

TYPE CH_SEAFLUX_t
!
  CHARACTER(LEN=6)              :: CCH_DRY_DEP              !  deposition scheme
  REAL, DIMENSION(:,:), POINTER :: XDEP                     ! final dry deposition  
                                                            ! velocity  for sea
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSV            ! name of the scalar var.                                                            
  INTEGER    :: NSV_CHSBEG, NSV_CHSEND                      ! chemical begin and ending
                                                            ! index of the HSV/CSV array
  INTEGER    :: NSV_DSTBEG, NSV_DSTEND                      ! dust begin and ending
  INTEGER    :: NSV_SLTBEG, NSV_SLTEND                      ! sea salt begin and ending
  INTEGER    :: NSV_AERBEG, NSV_AEREND                      ! aerosol begin and ending
  INTEGER    :: NBEQ                                        ! number of chemical species
  INTEGER    :: NDSTEQ                                      ! number of dust species
  INTEGER    :: NSLTEQ                                      ! number of sea salt species
  INTEGER    :: NAEREQ                                      ! number of aerosol species
                                                            ! in the surface scheme
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES      ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES                                                            
!
END TYPE CH_SEAFLUX_t



CONTAINS

!




SUBROUTINE CH_SEAFLUX_INIT(YCH_SEAFLUX)
TYPE(CH_SEAFLUX_t), INTENT(INOUT) :: YCH_SEAFLUX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_SEAFLUX_N:CH_SEAFLUX_INIT",0,ZHOOK_HANDLE)
NULLIFY(YCH_SEAFLUX%XDEP)
NULLIFY(YCH_SEAFLUX%CSV)
NULLIFY(YCH_SEAFLUX%CCH_NAMES)
NULLIFY(YCH_SEAFLUX%CAER_NAMES)
NULLIFY(YCH_SEAFLUX%CDSTNAMES)
NULLIFY(YCH_SEAFLUX%CSLTNAMES)
YCH_SEAFLUX%CCH_DRY_DEP=' '
YCH_SEAFLUX%NBEQ=0
YCH_SEAFLUX%NDSTEQ=0
YCH_SEAFLUX%NSLTEQ=0
YCH_SEAFLUX%NAEREQ=0
YCH_SEAFLUX%NSV_CHSBEG=0
YCH_SEAFLUX%NSV_CHSEND=0
YCH_SEAFLUX%NSV_DSTBEG=0
YCH_SEAFLUX%NSV_DSTEND=0
YCH_SEAFLUX%NSV_SLTBEG=0
YCH_SEAFLUX%NSV_SLTEND=0
YCH_SEAFLUX%NSV_AERBEG=0
YCH_SEAFLUX%NSV_AEREND=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_SEAFLUX_N:CH_SEAFLUX_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_SEAFLUX_INIT


END MODULE MODD_CH_SEAFLUX_n
