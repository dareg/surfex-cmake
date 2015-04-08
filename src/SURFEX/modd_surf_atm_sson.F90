!     ####################
      MODULE MODD_SURF_ATM_SSO_n
!     ######################
!
!!****  *MODD_SURF_ATM_SSO - declaration of surface parameters related to orography
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

TYPE SURF_ATM_SSO_t
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
  REAL, DIMENSION(:), POINTER :: XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM
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


END TYPE SURF_ATM_SSO_t

TYPE(SURF_ATM_SSO_t), ALLOCATABLE, TARGET, SAVE :: SURF_ATM_SSO_MODEL(:)

TYPE(SURF_ATM_SSO_t), POINTER :: SURF_ATM_SSO => NULL()
!$OMP THREADPRIVATE(SURF_ATM_SSO)

CONTAINS

SUBROUTINE SURF_ATM_SSO_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_SURF_ATM_SSO_N:SURF_ATM_SSO_GOTO_MODEL',0,ZHOOK_HANDLE)

SURF_ATM_SSO => SURF_ATM_SSO_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_SURF_ATM_SSO_N:SURF_ATM_SSO_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE SURF_ATM_SSO_GOTO_MODEL

SUBROUTINE SURF_ATM_SSO_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_SSO_N:SURF_ATM_SSO_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(SURF_ATM_SSO_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XAOSIP)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XAOSIM)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XAOSJP)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XAOSJM)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XHO2IP)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XHO2IM)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XHO2JP)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XHO2JM)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XZ0EFFIP)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XZ0EFFIM)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XZ0EFFJP)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XZ0EFFJM)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XZ0EFFJPDIR)  
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XZ0REL)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XSSO_SLOPE)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XSSO_ANIS)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XSSO_DIR)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XSSO_STDEV)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XAVG_ZS)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XSIL_ZS)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XMAX_ZS)
  NULLIFY(SURF_ATM_SSO_MODEL(J)%XMIN_ZS)
ENDDO
SURF_ATM_SSO_MODEL(:)%CROUGH=' '
SURF_ATM_SSO_MODEL(:)%XFRACZ0=2.
SURF_ATM_SSO_MODEL(:)%XCOEFBE=2.
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_SSO_N:SURF_ATM_SSO_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE SURF_ATM_SSO_ALLOC

SUBROUTINE SURF_ATM_SSO_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_SSO_N:SURF_ATM_SSO_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(SURF_ATM_SSO_MODEL)) DEALLOCATE(SURF_ATM_SSO_MODEL)
IF (ASSOCIATED(SURF_ATM_SSO)) NULLIFY(SURF_ATM_SSO)
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_SSO_N:SURF_ATM_SSO_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE SURF_ATM_SSO_DEALLO

END MODULE MODD_SURF_ATM_SSO_n
