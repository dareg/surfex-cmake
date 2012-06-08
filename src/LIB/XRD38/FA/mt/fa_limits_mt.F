C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FA_LIMITS_MT(FA, KPXPAH,KPXIND,KPXGEO,KPXNIV)
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK

      ! -----
      ! Purpose:
      !   Returns limits necessary for proper dimensioning of FA fields
      !   (arrays entering subroutine FACIES). It was created because
      !   of two reasons:
      !
      !   1) Subroutine FALIMU does not provide all needed limits.
      !
      !   2) Direct use of parameters JPX* from include "facomp.h" would
      !   be more convenient (no need to define interface arrays as
      !   allocatable). But since this include contains CPP directives,
      !   its use in external utilities might give wrong numbers with
      !   improper compile options. Moreover, include is written
      !   in fixed form and cannot be inserted into free form source.
      !
      ! Created:
      !   25-Jun-2010, J. Masek
      !
      ! Modified:
      !
      ! -----

      TYPE(FA_COM) :: FA
      INTEGER :: KPXPAH,KPXIND,KPXGEO,KPXNIV

      KPXPAH=FA%JPXPAH
      KPXIND=FA%JPXIND
      KPXGEO=FA%JPXGEO
      KPXNIV=FA%JPXNIV

      END SUBROUTINE 

