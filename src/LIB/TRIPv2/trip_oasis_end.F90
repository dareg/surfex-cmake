!#########
SUBROUTINE TRIP_OASIS_END
!########################
!
!!****  *TRIP_OASIS_END* - end coupling TRIP - OASIS
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
#ifdef TRIPOASIS
USE MOD_OASIS
#endif
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER                    :: IERR   ! Error info
!
!-------------------------------------------------------------------------------
#ifdef TRIPOASIS
!-------------------------------------------------------------------------------
!
CALL OASIS_TERMINATE(IERR)
IF (IERR/=OASIS_OK) THEN
   WRITE(*,'(A)'   )'Error OASIS terminate'
   WRITE(*,'(A,I4)')'Return code from oasis_terminate : ',IERR
   CALL ABORT
   STOP
ENDIF
!
!-------------------------------------------------------------------------------
#endif
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIP_OASIS_END
