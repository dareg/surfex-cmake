!     #####################
      MODULE MODD_CH_FLAKE_n
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
!!      Modified    04/2013, P. Le Moigne: FLake chemistry
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

TYPE CH_FLAKE_t
!
  CHARACTER(LEN=6)              :: CCH_DRY_DEP              ! deposition scheme
  REAL, DIMENSION(:,:), POINTER :: XDEP                     ! final dry deposition  
                                                            ! velocity  for lakes
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSV            ! name of the scalar var                                                            
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

                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES      ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES

!
END TYPE CH_FLAKE_t



CONTAINS

!





SUBROUTINE CH_FLAKE_INIT(YCH_FLAKE)
TYPE(CH_FLAKE_t), INTENT(INOUT) :: YCH_FLAKE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_FLAKE_N:CH_FLAKE_INIT",0,ZHOOK_HANDLE)
NULLIFY(YCH_FLAKE%XDEP)
NULLIFY(YCH_FLAKE%CSV)
NULLIFY(YCH_FLAKE%CCH_NAMES)
NULLIFY(YCH_FLAKE%CAER_NAMES)
NULLIFY(YCH_FLAKE%CDSTNAMES)
NULLIFY(YCH_FLAKE%CSLTNAMES)
YCH_FLAKE%CCH_DRY_DEP=' '
YCH_FLAKE%NSV_CHSBEG=0
YCH_FLAKE%NSV_CHSEND=0
YCH_FLAKE%NSV_DSTBEG=0
YCH_FLAKE%NSV_DSTEND=0
YCH_FLAKE%NBEQ=0
YCH_FLAKE%NSV_SLTBEG=0
YCH_FLAKE%NSV_SLTEND=0
YCH_FLAKE%NSV_AERBEG=0
YCH_FLAKE%NSV_AEREND=0
YCH_FLAKE%NDSTEQ=0
YCH_FLAKE%NSLTEQ=0
YCH_FLAKE%NAEREQ=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_FLAKE_N:CH_FLAKE_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_FLAKE_INIT


END MODULE MODD_CH_FLAKE_n
