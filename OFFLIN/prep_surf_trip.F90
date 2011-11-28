!     #################################################################################
SUBROUTINE PREP_SURF_TRIP(HPROGRAM)
!     #################################################################################
!
!!****  *PREP_SURF_TRIP* - driver for TRIP fields preparation
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
!!     B. Decharme 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/2008
!!------------------------------------------------------------------
!
!
USE MODI_GET_CONF_ISBA_n
USE MODI_PREP_TRIP
USE MODI_PREP_COUPLING_SURF_TRIP
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_GET_CONF_ISBA_n
!
USE MODI_GET_TYPE_DIM_n
!
USE MODI_SURF_VERSION
!
USE MODI_GET_LUOUT
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM    ! program calling surf. schemes
!
!
!*      0.2    declarations of local variables
!
CHARACTER(LEN=10) :: YGRID
!
LOGICAL :: LTRIP,LFLOOD
!
INTEGER :: ILU, ILUOUT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PREP_SURF_TRIP',0,ZHOOK_HANDLE)
CALL SURF_VERSION
!-------------------------------------------------------------------------------------
!
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
! * 1. Get ISBA configuration
!      
CALL GET_CONF_ISBA_n(LTRIP,LFLOOD,YGRID)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! TRIP parameters configuration:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF(.NOT.LTRIP)THEN
!        
  IF(LFLOOD)THEN
    WRITE(ILUOUT,*)'Error : In NAM_SGH_ISBAn, LFLOOD = True but LTRIP = False'
    CALL ABOR1_SFX('PREP_SURF_TRIP: LFLOOD=T BUT LTRIP=F')
  ENDIF
!
ELSE
!
  CALL GET_TYPE_DIM_n('FULL  ',ILU) 
!
  CALL PREP_COUPLING_SURF_TRIP(HPROGRAM,ILU,LFLOOD,YGRID)
!
  CALL PREP_TRIP(HPROGRAM)
!
ENDIF
IF (LHOOK) CALL DR_HOOK('PREP_SURF_TRIP',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_SURF_TRIP
