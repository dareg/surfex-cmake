!     ####################
      MODULE MODD_SSO1_n
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

TYPE SSO_1P_t
!
!-----------------------------------------------------------------------------------------------------
!
! Type of roughness
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
  REAL, DIMENSION(:), POINTER :: XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM
! directional total roughness lenghts in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
!
  REAL, DIMENSION(:), POINTER   :: XZ0REL         ! relief roughness length                 (m)
!
  REAL, DIMENSION(:), POINTER   :: XSSO_SLOPE         ! slope of S.S.O.
!   
!-----------------------------------------------------------------------------------------------------
!
END TYPE SSO_1P_t
!
TYPE SSO_NP_t
!
TYPE(SSO_1P_t), ALLOCATABLE :: AL(:) 
!
END TYPE SSO_NP_t
!
CONTAINS
!
SUBROUTINE SSO_1P_INIT(ISS)
TYPE(SSO_1P_t), INTENT(INOUT) :: ISS
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SSO1_N:SSO_1P_INIT",0,ZHOOK_HANDLE)
!
  NULLIFY(ISS%XAOSIP)
  NULLIFY(ISS%XAOSIM)
  NULLIFY(ISS%XAOSJP)
  NULLIFY(ISS%XAOSJM)
  NULLIFY(ISS%XHO2IP)
  NULLIFY(ISS%XHO2IM)
  NULLIFY(ISS%XHO2JP)
  NULLIFY(ISS%XHO2JM)
  NULLIFY(ISS%XZ0EFFIP)
  NULLIFY(ISS%XZ0EFFIM)
  NULLIFY(ISS%XZ0EFFJP)
  NULLIFY(ISS%XZ0EFFJM) 
  NULLIFY(ISS%XZ0REL)
  NULLIFY(ISS%XSSO_SLOPE)
!
IF (LHOOK) CALL DR_HOOK("MODD_SSO1_N:SSO_1P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE SSO_1P_INIT
!
SUBROUTINE SSO_NP_INIT(ISS,KPATCH)
TYPE(SSO_NP_t), INTENT(INOUT) :: ISS
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SSO1_N:SSO_NP_INIT",0,ZHOOK_HANDLE)
!
ALLOCATE(ISS%AL(KPATCH))
DO JP=1,KPATCH
  CALL SSO_1P_INIT(ISS%AL(JP))
ENDDO
!
IF (LHOOK) CALL DR_HOOK("MODD_SSO1_N:SSO_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE SSO_NP_INIT
!
END MODULE MODD_SSO1_n
