C Jan-2011 P. Marguinaud Interface to thread-safe FA
      SUBROUTINE FAIXLA
      USE FA_MOD, ONLY : FA_COM_DEFAULT,
     S                   FA_COM_DEFAULT_INIT,
     S                   NEW_FA_DEFAULT
      IMPLICIT NONE
      IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()
      CALL FAIXLA_MT (FA_COM_DEFAULT)
      END SUBROUTINE

