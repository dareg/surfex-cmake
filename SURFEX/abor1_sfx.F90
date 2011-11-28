!     #############################################################
      SUBROUTINE ABOR1_SFX(YTEXT)
!     #############################################################
!
!!****  *ABOR1_SFX* - abor1 subroutine
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
!!	P. Le Moigne   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2008 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_GET_LUOUT
USE MODI_CLOSE_FILE
USE MODD_SURF_CONF, ONLY : CPROGNAME
!      
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=*),  INTENT(IN)  :: YTEXT
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
CHARACTER(LEN=6)  :: YPROGRAM      
INTEGER           :: ILUOUT         ! logical unit of output file      
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!* get output listing file logical unit
!
IF (LHOOK) CALL DR_HOOK('ABOR1_SFX',0,ZHOOK_HANDLE)
YPROGRAM = CPROGNAME
!      
CALL GET_LUOUT(YPROGRAM,ILUOUT)
!
IF (YPROGRAM=='ASCII ' .OR. YPROGRAM=='TEXTE ' .OR. YPROGRAM=='BINARY') THEN
   WRITE(*,*)YTEXT
   WRITE(*,*)'---------------------------------------------------------------------------'
   WRITE(*,*) 'MORE DETAILS ABOUT THE CRASH IN THE OUTPUT LISTING: SEE THE FILE NAMED    '
   WRITE(*,*) 'LISTING_[NAME OF THE RUNNING .EXE: PGD, PREP, OFFLINE].txt              '
   WRITE(*,*)'---------------------------------------------------------------------------'
ENDIF
!
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '--------------------   FATAL ERROR in SURFEX  -----------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '-                                                                         -'
WRITE(ILUOUT,*)YTEXT
WRITE(ILUOUT,*) '-                                                                         -'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
CALL CLOSE_FILE(YPROGRAM,ILUOUT)
!
CALL ABORT
STOP
IF (LHOOK) CALL DR_HOOK('ABOR1_SFX',1,ZHOOK_HANDLE)
!
END SUBROUTINE ABOR1_SFX
