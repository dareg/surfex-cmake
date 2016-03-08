!     ####################
      MODULE MODD_SSO_n
!     ######################
!
!!****  *MODD_SSO - declaration of surface parameters related to orography
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE SSO_t
!
!-----------------------------------------------------------------------------------------------------
!
! Type of roughness
!
 CHARACTER(LEN=4) :: CROUGH     ! type of orographic roughness
!                              ! 'NONE'
                               ! 'Z01D'
                               ! 'Z04D'
                               ! 'BE04'

!-----------------------------------------------------------------------------------------------------
!
! Subgrid orography parameters
!
  REAL, DIMENSION(:), POINTER :: XAOSIP,XAOSIM,XAOSJP,XAOSJM
! directional A/S quantities in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
! They are used in soil routines to compute effective roughness length
!
  REAL, DIMENSION(:), POINTER :: XHO2IP,XHO2IM,XHO2JP,XHO2JM
! directional h/2 quantities in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
! They are used in soil routines to compute effective roughness length
!
  REAL, DIMENSION(:,:), POINTER :: XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM
! directional total roughness lenghts in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
!
  REAL, DIMENSION(:), POINTER   :: XZ0EFFJPDIR    ! heading of J direction (deg from N clockwise)

  REAL, DIMENSION(:), POINTER   :: XZ0REL         ! relief roughness length                 (m)
!
  REAL, DIMENSION(:), POINTER   :: XSSO_SLOPE         ! slope of S.S.O.
  REAL, DIMENSION(:), POINTER   :: XSSO_ANIS          ! anisotropy of S.S.O.
  REAL, DIMENSION(:), POINTER   :: XSSO_DIR           ! direction of S.S.O. (deg from N clockwise) 
  REAL, DIMENSION(:), POINTER   :: XSSO_STDEV         ! S.S.O. standard deviation           (m)
!
!
  REAL, DIMENSION(:), POINTER   :: XAVG_ZS        ! averaged orography                      (m)
  REAL, DIMENSION(:), POINTER   :: XSIL_ZS        ! silhouette orography                    (m)
  REAL, DIMENSION(:), POINTER   :: XMAX_ZS        ! maximum subgrid orography               (m)
  REAL, DIMENSION(:), POINTER   :: XMIN_ZS        ! minimum subgrid orography               (m)
! Zo threshold
  REAL   :: XFRACZ0                                ! Z0=Min(Z0, Href/XFRACZ0)
  REAL   :: XCOEFBE                                ! Beljaars coefficient         
!-----------------------------------------------------------------------------------------------------
!
!


END TYPE SSO_t



CONTAINS

!




SUBROUTINE SSO_INIT(YSSO)
TYPE(SSO_t), INTENT(INOUT) :: YSSO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SSO_N:SSO_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YSSO%XAOSIP)
  NULLIFY(YSSO%XAOSIM)
  NULLIFY(YSSO%XAOSJP)
  NULLIFY(YSSO%XAOSJM)
  NULLIFY(YSSO%XHO2IP)
  NULLIFY(YSSO%XHO2IM)
  NULLIFY(YSSO%XHO2JP)
  NULLIFY(YSSO%XHO2JM)
  NULLIFY(YSSO%XZ0EFFIP)
  NULLIFY(YSSO%XZ0EFFIM)
  NULLIFY(YSSO%XZ0EFFJP)
  NULLIFY(YSSO%XZ0EFFJM)
  NULLIFY(YSSO%XZ0EFFJPDIR)  
  NULLIFY(YSSO%XZ0REL)
  NULLIFY(YSSO%XSSO_SLOPE)
  NULLIFY(YSSO%XSSO_ANIS)
  NULLIFY(YSSO%XSSO_DIR)
  NULLIFY(YSSO%XSSO_STDEV)
  NULLIFY(YSSO%XAVG_ZS)
  NULLIFY(YSSO%XSIL_ZS)
  NULLIFY(YSSO%XMAX_ZS)
  NULLIFY(YSSO%XMIN_ZS)
YSSO%CROUGH=' '
YSSO%XFRACZ0=2.
YSSO%XCOEFBE=2.
IF (LHOOK) CALL DR_HOOK("MODD_SSO_N:SSO_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE SSO_INIT


END MODULE MODD_SSO_n
