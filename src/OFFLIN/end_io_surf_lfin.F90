!     #########
      SUBROUTINE END_IO_SURF_LFI_n(HPROGRAM)
!     #######################################################
!
!!****  *END_IO_SURF_LFI_n* - routine to close IO files
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
!!      S.Malardel   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2003 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_IO_SURF_LFI, ONLY : CLUOUT_LFI, CFILE_LFI, NFULL, CMASK, CFILEOUT_LFI, CFILEIN_LFI, &
                             NMASK
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
INTEGER :: IRET ! error code
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('END_IO_SURF_LFI_N',0,ZHOOK_HANDLE)
!
NFULL = 0
!
CMASK = '      '
!
NMASK=>NULL()
!
IF (CFILE_LFI==CFILEOUT_LFI .AND. NRANK==NPIO .OR. CFILE_LFI==CFILEIN_LFI) THEN
  CALL FMCLOS(CFILE_LFI,'KEEP',CLUOUT_LFI,IRET)
ENDIF
!
CFILE_LFI = '                              '
!
IF (LHOOK) CALL DR_HOOK('END_IO_SURF_LFI_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE END_IO_SURF_LFI_n
