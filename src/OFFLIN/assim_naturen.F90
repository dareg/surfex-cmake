!     ###############################################################################
SUBROUTINE ASSIM_NATURE_n(HPROGRAM,KI,                                    &
                          PCON_RAIN, PSTRAT_RAIN, PCON_SNOW, PSTRAT_SNOW, &
                          PCLOUDS,   PLSM,        PEVAPTR,   PEVAP,       & 
                          PSWEC,     PTSC,                                &
                          PTS,       PT2M,        PHU2M,     PSWE,        &
                          HTEST )

!     ###############################################################################
!
!!****  *ASSIM_NATURE_n * - Chooses the surface assimilation schemes for natural continental parts  
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     T. Aspelien
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2012
!!--------------------------------------------------------------------
!
USE MODD_SURF_ATM_n, ONLY : CNATURE
USE YOMHOOK,         ONLY : LHOOK,   DR_HOOK
USE PARKIND1,        ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_ASSIM_ISBA_n
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),    INTENT(IN) :: HPROGRAM  ! program calling surf. schemes
INTEGER,             INTENT(IN) :: KI
REAL, DIMENSION(:), INTENT(IN) :: PCON_RAIN
REAL, DIMENSION(:), INTENT(IN) :: PSTRAT_RAIN
REAL, DIMENSION(:), INTENT(IN) :: PCON_SNOW
REAL, DIMENSION(:), INTENT(IN) :: PSTRAT_SNOW
REAL, DIMENSION(:), INTENT(IN) :: PCLOUDS
REAL, DIMENSION(:), INTENT(IN) :: PLSM
REAL, DIMENSION(:), INTENT(IN) :: PEVAPTR
REAL, DIMENSION(:), INTENT(IN) :: PEVAP
REAL, DIMENSION(:), INTENT(IN) :: PSWEC
REAL, DIMENSION(:), INTENT(IN) :: PTSC
REAL, DIMENSION(:), INTENT(IN) :: PTS
REAL, DIMENSION(:), INTENT(IN) :: PT2M
REAL, DIMENSION(:), INTENT(IN) :: PHU2M
REAL, DIMENSION(:), INTENT(IN) :: PSWE
CHARACTER(LEN=2),    INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
REAL(KIND=JPRB)                :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_N',0,ZHOOK_HANDLE)

IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_NATURE_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF

IF (CNATURE=='ISBA  ') THEN

  CALL ASSIM_ISBA_n(HPROGRAM,KI,                                    &
                    PCON_RAIN, PSTRAT_RAIN, PCON_SNOW, PSTRAT_SNOW, &
                    PCLOUDS,   PLSM,        PEVAPTR,   PEVAP,       &
                    PSWEC,     PTSC,                                &
                    PTS,       PT2M,        PHU2M,     PSWE,        &
                    HTEST )
ELSE
  WRITE(*,*) 'No assimilation done for scheme: ',TRIM(CNATURE)
END IF

IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_NATURE_n
