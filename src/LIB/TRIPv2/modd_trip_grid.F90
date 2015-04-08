!##################
MODULE MODD_TRIP_GRID
!##################
!
!!****  *MODD_TRIP_GRID - declaration of grid for TRIP scheme
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       05/2008
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE TRIP_GRID_t
!-------------------------------------------------------------------------------
!
INTEGER, POINTER, DIMENSION(:,:)  :: NGRCN       ! Flow direction (1->8)
INTEGER, POINTER, DIMENSION(:,:)  :: NSEQ        ! River sequence
INTEGER                           :: NSEQMAX     ! maximum down flow
INTEGER, POINTER, DIMENSION(:,:)  :: NNEXTX      ! returns x and y point 
INTEGER, POINTER, DIMENSION(:,:)  :: NNEXTY      ! of destination grid:
!                                                        8 1 2
!                                                        7   3
!                                                        6 5 4
INTEGER, POINTER, DIMENSION(:,:)  :: NBASID      ! basin number id
INTEGER                           :: NBASMIN     ! minimum basin number
INTEGER                           :: NBASMAX     ! maximum basin number
!-------------------------------------------------------------------------------
!
REAL, POINTER,  DIMENSION(:)     :: XTRIP_GRID ! lits of parameters used to define the grid
!
!-------------------------------------------------------------------------------
!
REAL, POINTER, DIMENSION(:,:)    :: XAREA      ! 2d grid area [m*m]
REAL, POINTER, DIMENSION(:,:)    :: XLEN       ! distance between grids       [m]
!
!-------------------------------------------------------------------------------
!
LOGICAL, POINTER, DIMENSION(:,:) :: GMASK      !Logical Mask for TRIP grid
LOGICAL, POINTER, DIMENSION(:,:) :: GMASK_VEL  !Logical Mask for variable velocity scheme
LOGICAL, POINTER, DIMENSION(:,:) :: GMASK_GW   !Logical Mask for Groundwater grid
LOGICAL, POINTER, DIMENSION(:,:) :: GMASK_FLD  !Logical Mask for floodplain scheme
LOGICAL, POINTER, DIMENSION(:,:) :: GMASK_GRE  !Logical Mask for Greenland grid
LOGICAL, POINTER, DIMENSION(:,:) :: GMASK_ANT  !Logical Mask for Antartactic grid
!
!-------------------------------------------------------------------------------
!
END TYPE TRIP_GRID_t
!
TYPE(TRIP_GRID_t), ALLOCATABLE, TARGET, SAVE :: TRIP_GRID_MODEL(:)

TYPE(TRIP_GRID_t), POINTER :: TRIP_GRID => NULL()
!$OMP THREADPRIVATE(TRIP_GRID)

CONTAINS

SUBROUTINE TRIP_GRID_GOTO_MODEL(KTO)
INTEGER, INTENT(IN) :: KTO
REAL(KIND=JPRB)     :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TRIP_GRID:TRIP_GRID_GOTO_MODEL',0,ZHOOK_HANDLE)

TRIP_GRID => TRIP_GRID_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_TRIP_GRID:TRIP_GRID_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE TRIP_GRID_GOTO_MODEL

SUBROUTINE TRIP_GRID_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TRIP_GRID:TRIP_GRID_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TRIP_GRID_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(TRIP_GRID_MODEL(J)%NGRCN)
  NULLIFY(TRIP_GRID_MODEL(J)%NSEQ)
  NULLIFY(TRIP_GRID_MODEL(J)%NNEXTX)
  NULLIFY(TRIP_GRID_MODEL(J)%NNEXTY)
  NULLIFY(TRIP_GRID_MODEL(J)%NBASID)
  NULLIFY(TRIP_GRID_MODEL(J)%XTRIP_GRID)
  NULLIFY(TRIP_GRID_MODEL(J)%XAREA)
  NULLIFY(TRIP_GRID_MODEL(J)%XLEN)
  NULLIFY(TRIP_GRID_MODEL(J)%GMASK)
  NULLIFY(TRIP_GRID_MODEL(J)%GMASK_VEL)
  NULLIFY(TRIP_GRID_MODEL(J)%GMASK_GW)
  NULLIFY(TRIP_GRID_MODEL(J)%GMASK_FLD)
  NULLIFY(TRIP_GRID_MODEL(J)%GMASK_GRE)
  NULLIFY(TRIP_GRID_MODEL(J)%GMASK_ANT)
ENDDO
TRIP_GRID_MODEL(:)%NSEQMAX=0
IF (LHOOK) CALL DR_HOOK("MODD_TRIP_GRID:TRIP_GRID_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TRIP_GRID_ALLOC

SUBROUTINE TRIP_GRID_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TRIP_GRID:TRIP_GRID_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TRIP_GRID_MODEL)) DEALLOCATE(TRIP_GRID_MODEL)
IF (ASSOCIATED(TRIP_GRID)) NULLIFY(TRIP_GRID)
IF (LHOOK) CALL DR_HOOK("MODD_TRIP_GRID:TRIP_GRID_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TRIP_GRID_DEALLO

END MODULE MODD_TRIP_GRID
