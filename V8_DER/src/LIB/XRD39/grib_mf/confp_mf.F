      SUBROUTINE CONFP_MF (PFVAL,KEXP,KMANT)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C
C
C
C
C
C********************************************************************
C*
C*    NAME      : CONFP
C*
C*    FUNCTION  : CONVERT FLOATING POINT NUMBER FROM MACHINE
C*                REPRESENTATION TO GRIB REPRESENTATION.
C*
C*    INPUT     : PFVAL - FLOATING POINT NUMBER TO BE CONVERTED.
C*
C*    OUTPUT    : KEXP  - 8 BIT SIGNED EXPONENT
C*                KMANT - 24 BIT MANTISSA
C*                PFVAL - UNCHANGED.
C*
C*    JOHN HENNESSY , ECMWF , 15TH APRIL 1985
C*
C*    Modified by Jean CLOCHARD, February 1990, to comply with "DOCTOR",
C*    and to get a better precision on the truncation made:
C*       according to FM92 GRIB specification and recommandation,
C*    replacement of "1.0/ALOG(2.0)" by "1.+ZEPS", where ZEPS
C*    is a small positive value avoiding rounding errors to lead to
C*    a mantissa greater than (or equal to) 2**24.
C*
C*    (recommandation: ZEPS=1.E-12 for machines with word length of 60
C*                     bits or more, ZEPS=1.E-8 for 32 bits word length)
C*
C*       Tested with 10**5 pseudo-random values through the RANF
C*    function on Cyber 960 under NOS/VE 1.4.2, the modification gives
C*    an enhancement of more than 10, both on standard deviation and
C*    maximum error. Similar results are obtained when changing test
C*    interval from (0,1) to (-256,256).
C*
C*    Overflowing values are truncated, with a message on the listing.
C*
C********************************************************************
C
#include "precision.h"
C
      INTEGER KEXP, KMANT
C
      REAL (KIND=JPDBLR) PFVAL
C
      INTEGER IEXP
C
      REAL (KIND=JPDBLR) ZEPS, ZREF, ZC16
C
      INTRINSIC ABS, LOG, MAX, MIN
C
      SAVE ZEPS, ZC16
      DATA ZEPS / 1.E-12 /
      DATA ZC16 / 16. /
C
C
C
C           Elimination of sign.
C
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CONFP_MF',0,ZHOOK_HANDLE)
      ZREF = ABS (PFVAL)
C
C
C
C     EXPONENT
C
      IF (ZREF.EQ.0.0) THEN
        KEXP = 0
      ELSE
        IEXP = INT ( ( LOG(ZREF)/LOG(ZC16) + 65. ) +ZEPS )
        KEXP= MAX (0, MIN (127,IEXP))
      ENDIF
C
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
      KMANT = NINT ( ZREF/(ZC16**(KEXP-70)) )
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
C     ADD SIGN BIT TO EXPONENT.
C
      IF (PFVAL.LT.0.0) KEXP = KEXP + 128
C
C
C
C
      IF (LHOOK) CALL DR_HOOK('CONFP_MF',1,ZHOOK_HANDLE)
      END
