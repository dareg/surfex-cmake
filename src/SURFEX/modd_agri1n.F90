!     ##################
      MODULE MODD_AGRI1_n
!     ##################
!
!!****  *MODD_AGRI_n - declaration of SEEDING date for summer crops 
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
!!      P. LE MOIGNE   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/2006
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
TYPE AGRI_1P_t
!
INTEGER, POINTER, DIMENSION (:)   :: NIRRINUM       
                                        ! Stage for Irrigation (4 stages)
!
LOGICAL, POINTER,DIMENSION(:)     :: LIRRIGATE 
                                        ! True if irrigation performed
!
LOGICAL, POINTER,DIMENSION(:)     :: LIRRIDAY 
                                        ! True if irrigation occurs during present day
!                                          
REAL, POINTER, DIMENSION(:)       :: XTHRESHOLDSPT 
                                        ! Spatialized threshold
!
END TYPE AGRI_1P_t
!-------------------------------------------------------------------------------
!
TYPE AGRI_NP_t
!
TYPE(AGRI_1P_t), ALLOCATABLE :: AL(:) 
!
END TYPE AGRI_NP_t
!
CONTAINS
!
SUBROUTINE AGRI_1P_INIT(AG)
TYPE(AGRI_1P_t), INTENT(INOUT) :: AG
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_AGRI1_N:AGRI_1P_INIT",0,ZHOOK_HANDLE)
!
  NULLIFY(AG%NIRRINUM)
  NULLIFY(AG%LIRRIGATE)
  NULLIFY(AG%LIRRIDAY)
  NULLIFY(AG%XTHRESHOLDSPT)
!
IF (LHOOK) CALL DR_HOOK("MODD_AGRI1_N:AGRI_1P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE AGRI_1P_INIT
!
SUBROUTINE AGRI_NP_INIT(AG,KPATCH)
TYPE(AGRI_NP_t), INTENT(INOUT) :: AG
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_AGRI1_N:AGRI_NP_INIT",0,ZHOOK_HANDLE)
!
ALLOCATE(AG%AL(KPATCH))
!
DO JP=1,KPATCH
  CALL AGRI_1P_INIT(AG%AL(JP))
ENDDO
!
IF (LHOOK) CALL DR_HOOK("MODD_AGRI1_N:AGRI_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE AGRI_NP_INIT

!-------------------------------------------------------------------------------
!
END MODULE MODD_AGRI1_n
