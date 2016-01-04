!##################
MODULE MODD_ISBA_PGD_n
!##################
!
!!****  *MODD_ISBA - declaration of packed surface parameters for ISBA scheme
!!
!!    PURPOSE
!!    -------
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
!!      A. Boone   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/09/02
!!      A.L. Gibelin    04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin    04/2009 : TAU_WOOD for NCB option 
!!      A.L. Gibelin    05/2009 : Add carbon spinup
!!      A.L. Gibelin    06/2009 : Soil carbon variables for CNT option
!!      A.L. Gibelin    07/2009 : Suppress RDK and transform GPP as a diagnostic
!!      A.L. Gibelin    07/2009 : Suppress PPST and PPSTF as outputs
!!      P. Samuelsson   02/2012 : MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PGD_t
! 
! General surface parameters:
!
REAL, POINTER, DIMENSION(:)   :: XZS               ! relief                                  (m)
REAL, POINTER, DIMENSION(:,:) :: XCOVER            ! fraction of each ecosystem              (-)
LOGICAL, POINTER, DIMENSION(:):: LCOVER            ! GCOVER(i)=T --> ith cover field is not 0.
!
! Topmodel statistics
!
REAL, POINTER, DIMENSION(:)      :: XTI_MIN,XTI_MAX,XTI_MEAN,XTI_STD,XTI_SKEW
!
REAL, POINTER, DIMENSION(:,:)    :: XSAND          ! sand fraction                           (-)
REAL, POINTER, DIMENSION(:,:)    :: XCLAY          ! clay fraction                           (-)
REAL, POINTER, DIMENSION(:,:)    :: XSOC           ! soil organic carbon content             (kg/m2)
REAL, POINTER, DIMENSION(:)      :: XPERM          ! permafrost distribution                 (-)
REAL, POINTER, DIMENSION(:)      :: XGW            ! groundwater distribution                (-)
REAL, POINTER, DIMENSION(:)      :: XPH            ! soil pH
REAL, POINTER, DIMENSION(:)      :: XFERT          ! soil fertilisation rate (kgN/ha/h)
REAL, POINTER, DIMENSION(:)      :: XRUNOFFB       ! sub-grid dt92 surface runoff slope parameter (-)  
REAL, POINTER, DIMENSION(:)      :: XWDRAIN        ! continuous drainage parameter           (-)
!
! Subgrid orography parameters
!
! directional A/S quantities in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
! They are used in soil routines to compute effective roughness length
!
REAL, DIMENSION(:), POINTER :: XAOSIP,XAOSIM,XAOSJP,XAOSJM
!
! directional h/2 quantities in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
! They are used in soil routines to compute effective roughness length
!
REAL, DIMENSION(:), POINTER :: XHO2IP,XHO2IM,XHO2JP,XHO2JM
!
REAL, DIMENSION(:), POINTER   :: XZ0EFFJPDIR    ! heading of J direction (deg from N clockwise)
REAL, DIMENSION(:), POINTER   :: XSSO_SLOPE     ! slope of S.S.O.                         (-)
REAL, DIMENSION(:), POINTER   :: XSSO_STDEV     ! relief  standard deviation              (m)
!
END TYPE ISBA_PGD_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS

SUBROUTINE ISBA_PGD_INIT(YISBA_PGD)
TYPE(ISBA_PGD_t), INTENT(INOUT) :: YISBA_PGD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PGD_N:ISBA_PGD_INIT",0,ZHOOK_HANDLE)
!
NULLIFY(YISBA_PGD%XZS)
NULLIFY(YISBA_PGD%XCOVER)
NULLIFY(YISBA_PGD%LCOVER)
!
NULLIFY(YISBA_PGD%XTI_MIN)
NULLIFY(YISBA_PGD%XTI_MAX)
NULLIFY(YISBA_PGD%XTI_MEAN)
NULLIFY(YISBA_PGD%XTI_STD)
NULLIFY(YISBA_PGD%XTI_SKEW)
!
NULLIFY(YISBA_PGD%XSAND)
NULLIFY(YISBA_PGD%XCLAY)
NULLIFY(YISBA_PGD%XSOC)
NULLIFY(YISBA_PGD%XPERM)
NULLIFY(YISBA_PGD%XGW)
NULLIFY(YISBA_PGD%XPH)
NULLIFY(YISBA_PGD%XFERT)
NULLIFY(YISBA_PGD%XRUNOFFB)
NULLIFY(YISBA_PGD%XWDRAIN)
!
NULLIFY(YISBA_PGD%XAOSIP)
NULLIFY(YISBA_PGD%XAOSIM)
NULLIFY(YISBA_PGD%XAOSJP)
NULLIFY(YISBA_PGD%XAOSJM)
NULLIFY(YISBA_PGD%XHO2IP)
NULLIFY(YISBA_PGD%XHO2IM)
NULLIFY(YISBA_PGD%XHO2JP)
NULLIFY(YISBA_PGD%XHO2JM)
NULLIFY(YISBA_PGD%XZ0EFFJPDIR)
NULLIFY(YISBA_PGD%XSSO_SLOPE)
NULLIFY(YISBA_PGD%XSSO_STDEV)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PGD_N:ISBA_PGD_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PGD_INIT

END MODULE MODD_ISBA_PGD_n
