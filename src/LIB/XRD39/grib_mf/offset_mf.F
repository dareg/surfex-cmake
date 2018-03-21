      SUBROUTINE OFFSET_MF  (KOFF,KVAL,KWORD,KBYTE,KNBIT,KLEN,KERR)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C
C
C
C
C
C********************************************************************
C*
C*    NAME      : OFFSET
C*
C*    FUNCTION  : CALCULATES THE WORD AND BIT OFFSET OF THE START
C*                OF THE NEXT BIT-FIELD IN AN ARRAY OF BINARY
C*                DATA, FOLLOWING THE INSERTION/EXTRACTION OF ONE
C*                OR MORE CONSECUTIVE BIT-FIELDS.
C*
C*    INPUT     : KOFF  - BIT OFFSET AT WHICH LAST INSERTION OR
C*                        EXTRACTION STARTED.
C*                KVAL - NUMBER OF FIELDS INSERTED OR EXTRACTED IN
C*                        LAST OPERATION.
C*                KWORD - WORD NUMBER OF ARRAY AT WHICH LAST OPERATION
C*                        STARTED.
C*                KBYTE - LENGTH, IN BITS, OF LAST FIELD(S) INSERTED
C*                        OR EXTRACTED.
C*                KNBIT - NUMBER OF BITS IN COMPUTER WORD.
C*                KLEN  - NUMBER OF WORDS IN ARRAY.
C*
C*    OUTPUT    : KOFF  - BIT OFFSET AT WHICH NEXT INSERTION OR
C*                        EXTRACTION STARTS.
C*                KVAL - UNCHANGED.
C*                KWORD - WORD NUMBER OF ARRAY AT WHICH NEXT OPERATION
C*                        STARTS.
C*                KBYTE - UNCHANGED.
C*                KNBIT - UNCHANGED.
C*                KLEN  - UNCHANGED.
C*
C*                KERR  - EQUAL 0 IF NO ERROR.
C*                        EQUAL -2 IF NEXT OPERATION IS OUTSIDE ARRAY
C*                        BOUNDS.
C*
C*    GENERAL   : OFFSET CALLS -----
C*
C*    AUTHOR    : J.HENNESSY  15.4.85
C*
C*    MODIFIED  : J.HENNESSY  28.11.85
C*    Modified by Jean CLOCHARD, February 1990, to comply with "DOCTOR".
C*
C********************************************************************
C
      INTEGER KOFF, KVAL, KWORD, KBYTE, KNBIT, KLEN, KERR
C
      INTEGER IBITL, INTER, IOFFS
C
C
C     CLEAR ERROR INDICATOR
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('OFFSET_MF',0,ZHOOK_HANDLE)
      KERR = 0
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
C     CALCULATE NEXT WORD AND BIT POSITIONS.
C
      IBITL = KVAL * KBYTE
      INTER = IBITL / KNBIT
      IOFFS = IBITL - INTER * KNBIT
      KWORD = KWORD + INTER
      KOFF  = KOFF  + IOFFS
C
      IF (KOFF.GE.KNBIT)
     C   THEN
             KWORD = KWORD + 1
             KOFF  = KOFF  - KNBIT
         END IF
C
C
C
C
C
C
C
C
C     CHECK THAT NEXT WORD TO BE ACCESSED LIES WITHIN THE ARRAY BOUNDS.
C
      IF (KWORD.GT.KLEN)
     C   THEN
             KERR = -2
             WRITE (*,9001) KWORD,KLEN
 9001        FORMAT (1H ,'WORD ',I8,' IS OUTSIDE ARRAY BOUNDS ',I8)
         ENDIF
C
C
      IF (LHOOK) CALL DR_HOOK('OFFSET_MF',1,ZHOOK_HANDLE)
      END
