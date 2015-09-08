      SUBROUTINE PRTBIN_MF (KIN,KNBIT,KOUT,KERR)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C
C**** PRTBIN - Binary to decimal conversion.
C
C     Purpose.
C     --------
C
C           Produces a decimal number with ones and zeroes
C           corresponding to the ones and zeroes of the input
C           binary number.
C           eg input number 1011 binary, output number 1011 decimal.
C
C**   Interface.
C     ----------
C
C           CALL PRTBIN_MF (KIN,KNBIT,KOUT,KERR)
C
C           Integer    K.
C           Real       P.
C           Logical    O.
C           Character  H.
C
C               Input Parameters.
C               -----------------
C
C               KIN   - Integer variable containing binary number.
C
C               KNBIT - Number of bits in binary number.
C
C               Output Parameters.
C               -----------------
C
C               KOUT  - Integer variable containing decimal value
C                       with ones and zeroes corresponding to those of
C                       the input binary number.
C
C               KERR  - 0, If no error.
C                       1, Number of bits in binary number exceeds
C                          maximum allowed.
C
C     Method.
C     -------
C
C           Masking expression used is not ANSI standard.
C
C     Externals.
C     ----------
C
C           None.
C
C     Reference.
C     ----------
C
C           None.
C
C     Comments.
C     ---------
C
C           Routine contains sections 0, 1 and section 9.
C
C     Author.
C     -------
C
C           John Hennessy     ECMWF  October 1985
C
C     Modifications.
C     --------------
C
C           John Hennessy     ECMWF  March 1991
C           Made to conform to current programming standards.
C
C     ---------------------------------------------------------------
C
C
      INTEGER KIN, KNBIT, KOUT, KERR
C
      INTEGER J101, IK, IMASC, ITEMP
C
C
C
C
C
C
C
C
C*    Section 0. Definition of variables. Check on parameters.
C     -----------------------------------------------------------------
C
C     Check length of binary number.
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('PRTBIN_MF',0,ZHOOK_HANDLE)
      IF (KNBIT.LT.0.OR.KNBIT.GT.15)
     C   THEN
             KERR = 1
             WRITE (*,9000) KNBIT
             GO TO 900
         ELSE
             KERR = 0
         ENDIF
C
C     -----------------------------------------------------------------
C
C
C
C
C
C
C
C
C
C
C*    Section 1. Generate required number.
C     -----------------------------------------------------------------
C
  100 CONTINUE
C
      KOUT = 0
C
      DO 101 J101=1,KNBIT
         IK    = J101 - 1
         IMASC = 2**IK
         ITEMP = IAND(KIN,IMASC)
         IF (ITEMP.NE.0) KOUT = KOUT + 10**IK
  101 CONTINUE
C
C     -----------------------------------------------------------------
C
C
C
C
C
C
C
C
C
C
C*    Section 9. Format statements. Return to calling routine.
C     -----------------------------------------------------------------
C
  900 CONTINUE
C
 9000 FORMAT (1H ,'PRTBIN : Binary number too long - ',I3,' bits.')
C
      IF (LHOOK) CALL DR_HOOK('PRTBIN_MF',1,ZHOOK_HANDLE)
      END
