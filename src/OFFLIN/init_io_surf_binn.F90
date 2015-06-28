!     #########
      SUBROUTINE INIT_IO_SURF_BIN_n(HMASK,HACTION)
!     ######################
!
!!****  *INIT_IO_SURF_BIN_n* Keep in memory the output files
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
!!      A. Lemonsu  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_IO_SURF_BIN, ONLY : NMASK, NFULL, CMASK
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
!
IMPLICIT NONE
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HMASK    
 CHARACTER(LEN=5),  INTENT(IN)  :: HACTION    
!
INTEGER                        :: ILU, IRET, IL, ILUOUT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_BIN_N',0,ZHOOK_HANDLE)
!
 CALL GET_LUOUT('BINARY',ILUOUT)
!
 CALL GET_DIM_FULL_n(U, &
                     NFULL)
!
 CALL GET_SIZE_FULL_n(U, &
                      'BINARY ',NFULL,ILU)
!
IL = ILU
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     HMASK,IL)
 CALL INIT_IO_SURF_MASK_n(DTCO, U, &
                          HMASK, IL, ILUOUT, ILU, NMASK)
!
CMASK = HMASK
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_BIN_N',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
END SUBROUTINE INIT_IO_SURF_BIN_n
