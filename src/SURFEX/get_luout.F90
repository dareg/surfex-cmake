!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE GET_LUOUT(HPROGRAM,KLUOUT)
!     #######################################################
!
!!****  *GET_LUOUT* - routine to get output listing logical unit
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
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
#ifdef SFX_LFI
USE MODI_LFIGET_LUOUT
#endif
#ifdef SFX_MNH
USE MODI_MNHGET_LUOUT
#endif
USE MODD_SURFEX_HOST
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling GROUND
INTEGER,           INTENT(OUT) :: KLUOUT   ! Logical unit of output listing
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_LUOUT',0,ZHOOK_HANDLE)
IF (HPROGRAM=='MESONH') THEN
#ifdef SFX_MNH
  CALL MNHGET_LUOUT(HPROGRAM,KLUOUT)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
  CALL YRSURFEX_HOST%GET_LUOUT(HPROGRAM,KLUOUT)
ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef SFX_LFI
  CALL LFIGET_LUOUT(HPROGRAM,KLUOUT)
#endif
ELSE
  KLUOUT = 10
END IF
IF (LHOOK) CALL DR_HOOK('GET_LUOUT',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_LUOUT

