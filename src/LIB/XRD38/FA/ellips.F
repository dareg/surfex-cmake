C Jan-2011 P. Marguinaud Interface to thread-safe FA
      SUBROUTINE ELLIPS(KSMAX,KMSMAX,KNTMP,KMTMP)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C
C ***ELLIPS*** - General routine for computing elliptic truncation
C
C    Purpose.
C    --------
C       Computation of zonal and meridional limit wavenumbers within the ellipse
C    Interface:
C    ----------
C                   *CALL* *ELLIPS *
C
C        Explicit arguments :
C        --------------------
C
C        Implicit arguments :
C        --------------------
C
C
C     Method.
C     -------
C        See documentation
C
C     Externals.   NONE.
C     ----------
C
C     Reference.
C     ----------
C        ARPEGE/ALADIN documentation
C
C     Author.
C     -------
C        G. Radnoti LACE 97/04/04
C
C     Modifications.
C-------------------------------------------------------------
C        J.Vivoda, 99/05/19  treating NSMAX=0 and NMSMAX=0
C
#include "precision.h"
C
      INTEGER KSMAX, KMSMAX
      INTEGER KNTMP(0:KMSMAX),KMTMP(0:KSMAX)
C
      INTEGER JM, JN
C
      REAL (KIND=JPDBLR) ZEPS, ZKN, ZKM, ZAUXIL
C
       REAL(KIND=JPRB) :: ZHOOK_HANDLE
       IF (LHOOK) CALL DR_HOOK('ELLIPS',0,ZHOOK_HANDLE)
       ZEPS=1.E-10
       ZAUXIL=0.
C
C 1. Computing meridional limit wavenumbers along zonal wavenumbers
C
 100  CONTINUE
      DO 110 JM=1,KMSMAX-1
      ZKN = REAL(KSMAX,JPDBLR)/REAL(KMSMAX,JPDBLR)*
     1 SQRT(MAX(ZAUXIL,REAL(KMSMAX**2-JM**2,JPDBLR)))
        KNTMP(JM)=INT(ZKN+ZEPS)
 110  CONTINUE

      IF( KMSMAX.EQ.0 )THEN
         KNTMP(0)=KSMAX
      ELSE
         KNTMP(0)=KSMAX
         KNTMP(KMSMAX)=0
      ENDIF
C
C 2. Computing zonal limit wavenumbers along meridional wavenumbers
C             
 200  CONTINUE
      DO 210 JN=1,KSMAX-1
       ZKM = REAL(KMSMAX,JPDBLR)/REAL(KSMAX,JPDBLR)*
     1 SQRT(MAX(ZAUXIL,REAL(KSMAX**2-JN**2,JPDBLR)))
        KMTMP(JN)=INT(ZKM+ZEPS)
 210  CONTINUE

      IF( KSMAX.EQ.0 )THEN
         KMTMP(0)=KMSMAX
      ELSE
         KMTMP(0)=KMSMAX
         KMTMP(KSMAX)=0
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('ELLIPS',1,ZHOOK_HANDLE)
      END      

