!     #########
SUBROUTINE PREP_ISBA_UNIF(KLUOUT,HSURF,PFIELD)
!     #################################################################################
!
!!****  *PREP_ISBA_UNIF* - prepares ISBA field from prescribed values
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      P. Samuelsson  02/2012  MEB
!!------------------------------------------------------------------
!
!
USE MODD_PREP,           ONLY : CINTERP_TYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_PREP_ISBA,      ONLY : XHUG_SURF, XHUG_ROOT, XHUG_DEEP,    &
                                  XTG_SURF, XTG_ROOT, XTG_DEEP,     &
                                  XWR_DEF, XWRV_DEF, XWRVN_DEF,     &
                                  XQC_DEF,                          &
                                  XHUGI_SURF, XHUGI_ROOT, XHUGI_DEEP  
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER,            INTENT(IN)  :: KLUOUT    ! output listing logical unit
 CHARACTER(LEN=7),   INTENT(IN)  :: HSURF     ! type of field
REAL, POINTER, DIMENSION(:,:,:) :: PFIELD    ! field to interpolate horizontally
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_ISBA_UNIF',0,ZHOOK_HANDLE)
SELECT CASE(HSURF)
!
!*      3.0    Orography
!
  CASE('ZS     ')
    ALLOCATE(PFIELD(1,1,1))
    PFIELD = 0.
   
!
!*      3.1    Profile of soil relative humidity
!
  CASE('WG     ')
    ALLOCATE(PFIELD(1,3,1))
    PFIELD(:,1,1) = XHUG_SURF
    PFIELD(:,2,1) = XHUG_ROOT
    PFIELD(:,3,1) = XHUG_DEEP

!*      3.2    Profile of soil humidity for ice

  CASE('WGI    ')
    ALLOCATE(PFIELD(1,3,1))
    PFIELD(:,1,1) = XHUGI_SURF
    PFIELD(:,2,1) = XHUGI_ROOT
    PFIELD(:,3,1) = XHUGI_DEEP

!*      3.3    Profile of temperatures

  CASE('TG     ')
    ALLOCATE(PFIELD(1,3,1))
    PFIELD(:,1,1) = XTG_SURF
    PFIELD(:,2,1) = XTG_ROOT
    PFIELD(:,3,1) = XTG_DEEP

!*      3.4    Other quantities

  CASE('WR     ')
    ALLOCATE(PFIELD(1,1,1))
    PFIELD = XWR_DEF

  CASE('WRV    ')
    ALLOCATE(PFIELD(1,1,NVEGTYPE))
    PFIELD = XWRV_DEF

  CASE('WRVN   ')
    ALLOCATE(PFIELD(1,1,NVEGTYPE))
    PFIELD = XWRVN_DEF

  CASE('TV     ')
    ALLOCATE(PFIELD(1,1,NVEGTYPE))
    PFIELD = XTG_SURF

  CASE('TC     ')
    ALLOCATE(PFIELD(1,1,NVEGTYPE))
    PFIELD = XTG_SURF

  CASE('QC     ')
    ALLOCATE(PFIELD(1,1,NVEGTYPE))
    PFIELD = XQC_DEF

  CASE('LAI    ')
    ALLOCATE(PFIELD(1,1,1))
    PFIELD = XUNDEF

  CASE('ICE_STO')
    ALLOCATE(PFIELD(1,1,1))
    PFIELD = 0.0

END SELECT
!
!*      4.     Interpolation method
!              --------------------
!
CINTERP_TYPE='UNIF  '
IF (LHOOK) CALL DR_HOOK('PREP_ISBA_UNIF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
END SUBROUTINE PREP_ISBA_UNIF
