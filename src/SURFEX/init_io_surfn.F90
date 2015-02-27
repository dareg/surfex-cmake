!     #########
      SUBROUTINE INIT_IO_SURF_n(HPROGRAM,HMASK,HSCHEME,HACTION)
!     #######################################################
!
!
!!****  *INIT_IO_SURF_n* - chooses the routine to initialize IO
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
!!      Modified    04/2004 by P. LeMoigne: add HACTION if ASCII mode selected
!!      Modified    01/2009 by B. Decjharme: add HPROGRAM if FA mode selected
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
#ifdef SFX_ASC
USE MODI_INIT_IO_SURF_ASC_n
#endif
#ifdef SFX_BIN
USE MODI_INIT_IO_SURF_BIN_n
#endif
#ifdef SFX_FA
USE MODI_INIT_IO_SURF_FA_n
#endif
#ifdef SFX_LFI
USE MODI_INIT_IO_SURF_LFI_n
#endif
#ifdef SFX_NC
USE MODI_INIT_IO_SURF_NC_n
#endif
#ifdef SFX_OL
USE MODI_INIT_IO_SURF_OL_n
#endif
#ifdef SFX_TXT
USE MODI_INIT_IO_SURF_TXT_n
#endif
#ifdef SFX_MNH
USE MODI_MNHINIT_IO_SURF_n
#endif
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
 CHARACTER(LEN=6),  INTENT(IN)  :: HMASK
 CHARACTER(LEN=6),  INTENT(IN)  :: HSCHEME  ! scheme used
 CHARACTER(LEN=5),  INTENT(IN)  :: HACTION  ! action performed ('READ ','WRITE')
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_N',0,ZHOOK_HANDLE)
IF (HPROGRAM=='MESONH') THEN
#ifdef SFX_MNH
  CALL MNHINIT_IO_SURF_n(HPROGRAM,HMASK,HACTION)
#endif
END IF
!
IF (HPROGRAM=='OFFLIN' ) THEN
#ifdef SFX_OL
  CALL INIT_IO_SURF_OL_n(HPROGRAM,HMASK,HSCHEME,HACTION)
#endif
ENDIF
!
IF (HPROGRAM=='ASCII ' ) THEN
#ifdef SFX_ASC
  CALL INIT_IO_SURF_ASC_n(HMASK,HACTION)
#endif
ENDIF
!
IF (HPROGRAM=='TEXTE ' ) THEN
#ifdef SFX_TXT
  CALL INIT_IO_SURF_TXT_n(HMASK,HACTION)
#endif
ENDIF
!
IF (HPROGRAM=='BINARY' ) THEN
#ifdef SFX_BIN
  CALL INIT_IO_SURF_BIN_n(HMASK,HACTION)
#endif
ENDIF
!
IF (HPROGRAM=='AROME ' ) THEN
#ifdef SFX_ARO
  CALL AROINIT_IO_SURF_n(HPROGRAM,HMASK,HACTION)
#endif
ENDIF
!
IF (HPROGRAM=='FA    ' ) THEN
#ifdef SFX_FA
  CALL INIT_IO_SURF_FA_n(HPROGRAM,HMASK,HACTION)
#endif
ENDIF
!
IF (HPROGRAM=='LFI   ' ) THEN
#ifdef SFX_LFI
  CALL INIT_IO_SURF_LFI_n(HMASK,HACTION)
#endif
ENDIF
!
IF (HPROGRAM=='NC    ' ) THEN
#ifdef SFX_NC
  CALL INIT_IO_SURF_NC_n(HMASK,HACTION)
#endif
ENDIF
!
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_IO_SURF_n
