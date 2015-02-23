      SUBROUTINE CONFI (PFVAL,KEXP,KMANT,PNFVAL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C
C
C
C
C
C********************************************************************
C*
C*    NAME      : CONFI
C*
C*    FUNCTION  : CONVERT FLOATING POINT NUMBER FROM MACHINE
C*                REPRESENTATION TO GRIB REPRESENTATION,
C*                THE RESULT NUMBER NOT EXCEEDING THE INPUT NUMBER.
C*    INPUT     : PFVAL  - FLOATING POINT NUMBER TO BE CONVERTED.
C*
C*    OUTPUT    : KEXP   - 8 BIT SIGNED EXPONENT
C*                KMANT  - 24 BIT MANTISSA
C*                PFVAL  - UNCHANGED
C*                PNFVAL - "EXACT" VALUE REPRESENTED BY (KEXP,KMANT).
C*
C*    Jean CLOCHARD , French DMN, January 1990.
C*    Nearly rewritten from CONFP subroutine from John HENNESSY, ECMWF.
C*
C*    Overflowing values are truncated, with a message on the listing.
C*
C********************************************************************
C
#include "precision.h"
C
      INTEGER KEXP, KMANT
C
      REAL (KIND=JPDBLR) PFVAL, PNFVAL
C
      INTEGER IAUXIL, IEXP
C
      REAL (KIND=JPDBLR) ZEPS, ZREF, ZC16, ZMANT
C
      LOGICAL LLPOSI
C
      INTRINSIC LOG, ABS, MIN, MAX, INT
C
      SAVE ZEPS, ZC16
      DATA ZEPS / 1.E-12 /
      DATA ZC16 / 16. /
C
C
C
C           Elimination of sign.
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CONFI',0,ZHOOK_HANDLE)
      LLPOSI=PFVAL.GE.0.0
      ZREF = ABS (PFVAL)
C
C
C
C
C     EXPONENT
C
      IF (ZREF.EQ.0.0) THEN
        KEXP = 0
      ELSE
        IEXP = INT ( ( LOG(ZREF)/LOG(16.0) + 65. ) +ZEPS )
        KEXP= MAX (0, MIN (127,IEXP))
      ENDIF
C
C
C
C
C
C
C     MANTISSA
C
  201 CONTINUE
C
      ZMANT = ZREF/(ZC16**(KEXP-70))
C
      IF (LLPOSI) THEN
        KMANT=INT (ZMANT)
      ELSE
C
C        Special case for negative values... because the "INT" function
C     is not equivalent to the "integer part" mathematical function
C     for this range of values.
C
        IAUXIL=2+INT (ZMANT)
        KMANT=IAUXIL - INT ( -ZMANT + REAL (IAUXIL,JPDBLR) )
      ENDIF
C
      IF (KMANT.GE.2**24) THEN
C
        IF (KEXP.LT.127) THEN
C
C     Some rounding error ocurred in the computation of KEXP, and could
C     not be compensated by ZEPS. Incrementation of KEXP, and new value
C     of KMANT computed.
C
          KEXP=KEXP+1
          GOTO 201
        ELSE
          PRINT *,
     S'*/*/* OVERFLOW OF GRIB FLOATING-POINT REPRESENTATION WITH ',PFVAL
          KMANT=2**24-1
        ENDIF
C
      ENDIF
C
C     COMPUTE "EXACT" VALUE REPRESENTED, AND ADD SIGN BIT TO EXPONENT.
C
      IF (LLPOSI) THEN
        PNFVAL =  KMANT*(ZC16**(KEXP-70))
      ELSE
        PNFVAL = -KMANT*(ZC16**(KEXP-70))
        KEXP=KEXP+128
      ENDIF
C
C
      IF (LHOOK) CALL DR_HOOK('CONFI',1,ZHOOK_HANDLE)
      END
