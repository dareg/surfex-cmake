      SUBROUTINE DECFP_MF (PFVAL,KEXP,KMANT)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C
C********************************************************************
C*
C*    NAME      : DECFP
C*
C*    FUNCTION  : DECODE GRIB CODE REPRESENTATION TO
C*                FLOATING POINT NUMBER.
C*
C*    INPUT     : KEXP  - 8 BIT SIGNED EXPONENT
C*                KMANT - 24 BIT MANTISSA
C*
C*    OUTPUT    : PFVAL - FLOATING POINT NUMBER
C*                KEXP  - UNCHANGED
C*                KMANT - UNCHANGED
C*
C*    JOHN HENNESSY , ECMWF , 15TH APRIL 1985
C*
C*    Modified by Jean CLOCHARD, February 1990, to comply with "DOCTOR",
C*    and to get a more straight-forward code.
C*
C********************************************************************
C
#include "precision.h"
C
      INTEGER KEXP, KMANT
C
      REAL (KIND=JPDBLR) PFVAL
C
      REAL (KIND=JPDBLR) ZC16
C
      SAVE ZC16
      DATA ZC16 / 16. /
C
C      Direct decoding, depending on sign bit, included in KEXP.
C
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('DECFP_MF',0,ZHOOK_HANDLE)
      IF (KEXP.LT.128) THEN
        PFVAL = KMANT*(ZC16**(KEXP-70))
      ELSE
        PFVAL = -KMANT*(ZC16**(KEXP-198))
      ENDIF
C
C
      IF (LHOOK) CALL DR_HOOK('DECFP_MF',1,ZHOOK_HANDLE)
      END
