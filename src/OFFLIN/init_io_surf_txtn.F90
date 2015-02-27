!     #########
      SUBROUTINE INIT_IO_SURF_TXT_n(HMASK,HACTION)
!     ######################
!
!!****  *INIT_IO_SURF_TXT_n* Keep in memory the output files
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      P. Le Moigne 04/2004: distinguish in and out file name
!!      P. Le Moigne 04/2006: special HACTION='GTMSK' to initialize
!!                            a mask different of 'FULL ' in order 
!!                            to read dimensions only.
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE MODD_IO_SURF_TXT, ONLY : NMASK, NFULL, CMASK
!
USE MODI_GET_LUOUT
USE MODI_GET_DIM_FULL_n
USE MODI_GET_SIZE_FULL_n
USE MODI_GET_TYPE_DIM_n
USE MODI_INIT_IO_SURF_MASK_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HMASK    
 CHARACTER(LEN=5),  INTENT(IN)  :: HACTION    
!
INTEGER                        :: ILU,IRET, IL, ILUOUT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_TXT_N',0,ZHOOK_HANDLE)
!
 CALL GET_LUOUT('TEXTE ',ILUOUT)
!
 CALL GET_DIM_FULL_n(NFULL)
!
 CALL GET_SIZE_FULL_n('TEXTE ',NFULL,ILU)
!
IL = ILU
 CALL GET_TYPE_DIM_n(HMASK,IL)
 CALL INIT_IO_SURF_MASK_n(HMASK, IL, ILUOUT, ILU, NMASK)
!
CMASK = HMASK
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_TXT_N',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
!
END SUBROUTINE INIT_IO_SURF_TXT_n
