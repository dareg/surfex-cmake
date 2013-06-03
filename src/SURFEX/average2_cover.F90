!     #########################
      SUBROUTINE AVERAGE2_COVER(HPROGRAM)
!     #########################
!
!!**** *AVERAGE2_COVER* computes the cover fractions
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson         Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    10/12/97
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_PGDWORK,   ONLY : NSIZE, XSUMCOVER
USE MODD_SURF_ATM_n, ONLY : XCOVER, LCOVER
!
USE MODD_PGD_GRID,       ONLY : CGRID
!
USE MODI_SUM_ON_ALL_PROCS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
 CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM      ! Type of program
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
REAL, DIMENSION(:), ALLOCATABLE :: ZUNITY
!
INTEGER :: JCOVER, ICPT ! loop counter on cover classes
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
!*    1.     Average values
!            --------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGE2_COVER',0,ZHOOK_HANDLE)
ALLOCATE(ZUNITY(SIZE(NSIZE)))
ZUNITY (:) = 0.
!
DO JCOVER=1,SIZE(XSUMCOVER,2)
  ICPT = SUM_ON_ALL_PROCS(HPROGRAM,CGRID,XSUMCOVER(:,JCOVER)/=0., 'COV')
  IF (ICPT>0) LCOVER(JCOVER) = .TRUE.
ENDDO
!
ALLOCATE(XCOVER(SIZE(NSIZE),COUNT(LCOVER)))
!
ICPT = 0
DO JCOVER=1,SIZE(XSUMCOVER,2)
  IF (LCOVER(JCOVER)) THEN
    ICPT = ICPT+1
    WHERE (NSIZE(:)/=0)
      XCOVER(:,ICPT)=XSUMCOVER(:,JCOVER) /NSIZE(:)
      ZUNITY(:)=ZUNITY(:) + XCOVER(:,ICPT)
    ENDWHERE
  ENDIF
END DO
!
DO JCOVER=1,SIZE(XCOVER,2)
  WHERE (NSIZE(:) /=0 )
    XCOVER(:,JCOVER)=XCOVER(:,JCOVER) / ZUNITY(:)
  END WHERE
END DO
!
!-------------------------------------------------------------------------------
DEALLOCATE(ZUNITY)
IF (LHOOK) CALL DR_HOOK('AVERAGE2_COVER',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE2_COVER
