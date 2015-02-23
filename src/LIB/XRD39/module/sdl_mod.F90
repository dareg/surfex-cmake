MODULE SDL_MOD

!RJ: stripped version of SDL_MOD for aborts in FA/LFI only, full in SDL_MODULE

USE PARKIND1  ,ONLY : JPIM  ,JPRB

IMPLICIT NONE

SAVE

PRIVATE

PUBLIC :: SDL_SRLABORT

CONTAINS

SUBROUTINE SDL_SRLABORT

! Purpose :
! -------
!   To abort in serial environment

!RJ: simple abort here is more portable, gives backtrace, but is extension
!RJ CALL EC_RAISE(SIGABRT)
CALL ABORT()
STOP 'SDL_SRLABORT'

END SUBROUTINE SDL_SRLABORT

END MODULE SDL_MOD
