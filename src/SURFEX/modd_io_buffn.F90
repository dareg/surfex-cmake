!     ##################
      MODULE MODD_IO_BUFF_n
!     ##################
!
!!****  *MODD_IO_IO_BUFF - 
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
!!      V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!    
!
!*       0.   DECLARATIONS
!
!!      
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE IO_BUFF_t

 CHARACTER(LEN=12), DIMENSION(3000) :: CREC   ! list of records already read/written
INTEGER                            :: NREC   ! number of records read/written

END TYPE IO_BUFF_t
!
!CONTAINS
!
!!
!
!
!!
!
!SUBROUTINE IO_BUFF_INIT(Y!SUBROUTINE IO_BUFF)
!TYPE(SUBROUTINE IO_BUFF_t), INTENT(INOUT) :: YSUBROUTINE IO_BUFF
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!IF (LHOOK) CALL DR_HOOK("MODD_IO_BUFF_N:IO_BUFF_INIT",0,ZHOOK_HANDLE)
!YIO_BUFF%NREC=0
!IF (LHOOK) CALL DR_HOOK("MODD_IO_BUFF_N:IO_BUFF_INIT",1,ZHOOK_HANDLE)
!END SUBROUTINE IO_BUFF_INIT
!
!
END MODULE MODD_IO_BUFF_n
