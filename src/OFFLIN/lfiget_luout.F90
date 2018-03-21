!     #########
      SUBROUTINE LFIGET_LUOUT(HPROGRAM,KLUOUT)
!     #######################################################
!
!!****  *LFIGET_LUOUT* - get output listing logical unit
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling GROUND
INTEGER,           INTENT(OUT) :: KLUOUT   ! Logical unit of output listing
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
                                    ! at the open of the file in LFI  routines 
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LFIGET_LUOUT',0,ZHOOK_HANDLE)
KLUOUT=10  ! This value is not used by FMATTR and FMOPEN routines, for which
IF (LHOOK) CALL DR_HOOK('LFIGET_LUOUT',1,ZHOOK_HANDLE)
!          ! logical unit start at number 11
!-------------------------------------------------------------------------------
!
END SUBROUTINE LFIGET_LUOUT
