      SUBROUTINE DECOGA (PFDATA,KLENF,KBITS,KNBIT,KB1PAR,
     C                   KB2PAR,PVERT,KLENV,KGRIB,KLENG,KWORD,
     C              KJLENV,KJLENF,KCPACK,KSCALP,KERR,PMIN,PMAX,LDARPE)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C
C
C
C
C
C
C********************************************************************
C*
C*     NAME     : DECOGA
C*
C*     FUNCTION : DECODE WMO GRIB CODED DATA.
C*
C*     INPUT    : KNBIT   - NUMBER OF BITS IN COMPUTER WORD
C*
C*                KGRIB   - INTEGER ARRAY CONTAINING DATA IN GRIB CODE.
C*                KLENG   - LENGTH OF ARRAY KGRIB.
C*
C*                PVERT   - REAL ARRAY TO RECEIVE VERTICAL COORDINATE
C*                          PARAMETERS.
C*                KLENV   - LENGTH OF ARRAY PVERT.
C*
C*                PFDATA  - REAL ARRAY TO RECEIVE DECODED DATA VALUES.
C*                KLENF   - LENGTH OF ARRAY PFDATA.
C*
C*                          *****  IF THIS IS 1 ONLY
C*                          THE PRODUCT AND GRID DEFINITION BLOCKS
C*                          ARE DECODED. *****
C*
C*                KB1PAR  - INTEGER ARRAY (DIMENSION 19) TO RECEIVE
C*                          PRODUCT DEFINITION INFORMATION.
C*
C*                KB2PAR  - INTEGER ARRAY (DIMENSION 11) TO RECEIVE GRID
C*                          DEFINITION INFORMATION.
C*
C*                PMIN    - MINIMUM VALUE OF THE FIELD DATA.
C*
C*                PMAX    - MAXIMUM VALUE OF THE FIELD DATA.
C*
C*                  *****   These 2 last values have only to be supplied
C*                          if next argument is .TRUE. .   *****
C*
C*                LDARPE  - .TRUE., modifications for ARPEGE coding
C*                                  have been included when coding data;
C*                          .FALSE., no such modifications.
C*
C*     OUTPUT   : PARAMETERS FROM BLOCK 1 , PRODUCT DEFINITION BLOCK.
C*                --------------------------------------------------
C*                KB1PAR  - INTEGER ARRAY (DIMENSION 19) OF PRODUCT
C*                          DEFINITION INFORMATION.
C*
C*                WORD      CONTENTS
C*                ----      --------
C*                  1       ORIGINATING CENTRE IDENTIFIER.
C*                  2       MODEL IDENTIFICATION.
C*                  3       GRID DEFINITION.
C*                  4       FLAG ( CODE TABLE 1)
C*                  5       PARAMETER IDENTIFIER ( CODE TABLE 2 ).
C*                  6       TYPE OF LEVEL (CODE TABLE 3).
C*                 7-8      VALUE(S) OF LEVELS (CODE TABLE 3).
C*                  9       YEAR OF DATA
C*                 10       MONTH OF DATA
C*                 11       DAY OF DATA
C*                 12       HOUR OF DATA
C*                 13       MINUTE OF DATA
C*                 14       TIME UNIT (CODE TABLE 4).
C*                 15       TIME RANGE ONE
C*                 16       TIME RANGE TWO
C*                 17       TIME RANGE FLAG (CODE TABLE 5).
C*                 18       NUMBER AVERAGED OR ACCUMULATED.
C*                 19       NUMBER MISSING FROM AVERAGES/ACCUMULATIONS.
C*
C*                          VALUE(S) OF LEVEL CAN OCCUPY 2 WORDS..
C*                          FOR A LAYER THE FIRST WORD DEFINES THE TOP
C*                          AND THE SECOND THE BOTTOM OF THE LAYER.
C*                          FOR A SINGLE LEVEL, ONLY THE FIRST WORD IS
C*                          USED.
C*
C*                PARAMETERS FROM BLOCK 2 , GRID DEFINITION BLOCK.
C*                -----------------------------------------------
C*                KB2PAR - INTEGER ARRAY (DIMENSION 17) CONTAINING GRID
C*                         DEFINITION INFORMATION.
C*                         USE VARIES WITH DATA REPRESENTATION TYPE.
C*
C*                WORD      LAT/LONG GRID
C*                ----      -------------
C*                 1        DATA REPRESENTATION TYPE (CODE TABLE 6).
C*                 2        NO. OF POINTS ALONG A LATITUDE
C*                 3        NO. OF POINTS ALONG A MERIDIAN.
C*                 4        LATITUDE OF ORIGIN (SOUTH - IVE).
C*                 5        LONGITUDE OF ORIGIN (WEST - IVE).
C*                 6        RESOLUTION FLAG. (CODE TABLE 7)
C*                 7        LATITUDE OF EXTREME POINT (SOUTH - IVE).
C*                 8        LONGITUDE OF EXTREME POINT (WEST - IVE).
C*                 9        LATITUDE INCREMENT
C*                10        LONGITUDE INCREMENT
C*                11        SCANNING MODE FLAGS (CODE TABLE 8)
C*
C*                WORD      GAUSSIAN GRID
C*                ----      -------------
C*                1-9       AS FOR LAT/LONGITUDE GRID.
C*                10        THE NUMBER OF LATITUDE LINES BETWEEN A POLE
C*                          AND THE EQUATOR.
C*                11        SCANNING MODE FLAGS (CODE TABLE 8)
C*
C*                WORD      SPHERICAL HARMONICS
C*                ----      -------------------
C*                  1       DATA REPRESENTATION TYPE (CODE TABLE 6)
C*                  2       J - PENTAGONAL RESOLUTION PARAMETER
C*                  3       K - PENTAGONAL RESOLUTION PARAMETER
C*                  4       M - PENTAGONAL RESOLUTION PARAMETER
C*                  5       REPRESENTATION TYPE (CODE TABLE 9)
C*                  6       REPRESENTATION MODE (CODE TABLE 10)
C*                7-11      NOT USED
C*
C*
C*                USE FOR OTHER DATA REPRESENTATION TYPES IS
C*
C*                PVERT  - REAL ARRAY OF VERTICAL COORDINATE PARAMETERS
C*                KJLENV - NUMBER OF VALUES IN THIS ARRAY.
C*
C*                USED FOR HYBRID LEVELS ONLY.
C*
C*                PARAMETERS FROM BLOCK 3 , BIT-MAP DEFINITION BLOCK.
C*                --------------------------------------------------
C*
C*                TO BE DEFINED LATER.
C*
C*                PARAMETERS FROM BLOCK 4 , BINARY DATA BLOCK.
C*                -------------------------------------------
C*                PFDATA  - ARRAY OF FLOATING POINT VALUES.
C*                KJLENF  - NUMBER OF VALUES IN THIS ARRAY.
C*                          KJLENF IS NEGATIVE IF MISSING DATA CODED
C*
C*                KBITS   - NUMBER OF BITS FOR CODED DATA VALUES.
C*
C*                KWORD   - NUMBER OF WORDS DECODED.
C*
C*                KCPACK  - KCPACK GREATER THAN ZERO INDICATES COMPLEX
C*                          PACKING, WITH A *TRIANGULAR* SUB-TRUNCATION
C*                          "UNPACKED" OF ORDER KCPACK.
C
C*                KSCALP  - FOR COMPLEX PACKING OF SPECTRAL COEFFICIENTS
C*                          ONLY, LAPLACIAN SCALING FACTOR.
C*
C*                 *****   THE 2 NEXT VALUES ARE SUPPLIED BY THE ROUTINE
C*                         ONLY IF ARGUMENT "LDARPE" IS .FALSE. .  *****
C*
C*                PMIN    - MINIMUM VALUE OF THE FIELD DATA, OR A SLIGHT
C*                          UNDER-APPROXIMATION OF THE "TRUE" MINIMUM OF
C*                          DECODED FIELD DATA.
C*
C*                PMAX    - OVER-ESTIMATION OF THE MAXIMUM VALUE OF THE
C*                          FIELD DATA (COMPUTED WITH PMIN, THE SCALE
C*                          FACTOR "ISCALE" OF GRIB, AND KBITS).
C*
C*
C*                KERR    - ERROR INDICATOR.
C*
C*                          0, NO ERROR.
C*
C*                         -1, NUMBER OF BITS PER DATA VALUE EXCEEDS
C*                             WORD LENGTH, OR MAXIMUM ALLOWED.
C*                         -2, INPUT ARRAY HOLDS FEWER CODED VALUES
C*                             THAN EXPECTED.
C*                         -3, OUTPUT ARRAY TOO SMALL TO HOLD DECODED
C*                             DATA VALUES.
C*                         -4, CURRENTLY UNDEFINED OPTION SELECTED.
C*
C*                          1, NO BLOCK 5 FOUND.
C                           2, ERROR IN CONVERTING OLD TO NEW FORM OF
C                              BLOCK 1 OR 2.
C*                          3, NO BLOCK 0 FOUND.
C*
C*    GENERAL  : DECOGA CALLS GBYTE
C*                            GSBYTE
C*                            OFFSET
C*                            DECFP
C*
C*
C*    MODIFIED : J. HENNESSY 17.06.87
C*               CONVERT OLD FORMAT BLOCK 1 AND 2 TO NEW FORMAT.
C*               ( DISTINCTION BETWEEN ANALYSIS AND INITIALISED
C*               NOT POSSIBLE ).
C*
C*    Modifications by Mats HAMRUD, ECMWF, to handle complex unpacking
C*    of spectral harmonics data, 1988.
C*
C*    Modifications by Jean CLOCHARD, French DMN, February 1990:
C*
C*    -  to remove automatic arrays introduced by the last modification;
C*    -  to comply with the "DOCTOR" norm;
C*    -  to allow modifications of the unpacking for ARPEGE files use,
C*       in an optional way.
C*
C********************************************************************
C
C
C
C
C
C
C
C
C

      USE SDL_MOD   , ONLY : SDL_SRLABORT

#include "precision.h"
C
      INTEGER KLENF, KBITS, KNBIT, KLENV, KLENG, KWORD, KJLENV, KJLENF
      INTEGER KCPACK, KSCALP, KERR
C
      INTEGER KB1PAR(19),KB2PAR(17)
      INTEGER (KIND=JPDBLE) KGRIB(KLENG)
C
      REAL (KIND=JPDBLR) PMIN, PMAX
C
      REAL (KIND=JPDBLR) PFDATA(KLENF), PVERT(KLENV)
C
      LOGICAL LDARPE
C
      INTEGER IMAX, IOFF, IBYTE, INVAL, ISNEW, J, ITEMP, INC, ILBLK
      INTEGER ILNIL, IEXP, IMANT, ILEN, IPW, IPB, ILBIN, IFLAG, IREP
      INTEGER ISCALX, ISCALE, IMISS, ISPDA, IPREMC, ISSUIV, IL, ILBVAL
      INTEGER (KIND=JPDBLE) ILEXP, ILMANT, ILFLAG, ILSCALX, ILLNIL
      INTEGER (KIND=JPDBLE) ILLBLK, ILLBIN, ILBITS, ILSPDA, ILSCALP
      INTEGER (KIND=JPDBLE) ILCPACK
C
C
      INTEGER IBLOCK(24), ILAT(2)
      INTEGER (KIND=JPDBLE) ILBLOCK(24), ILB2PAR(17), ILLAT(2)
C
      REAL (KIND=JPDBLR) ZSCALE
C
C
C     SET MAXIMUM NUMBER OF BITS PER DATA FIELD.
C
      SAVE IMAX
      DATA IMAX /60/
C
C     CLEAR ERROR INDICATOR.
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('DECOGA',0,ZHOOK_HANDLE)
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
C
C
C
C
C
C
C
C********************************************************************
C*
C*    BLOCK 0 - INDICATOR BLOCK.
C*
C********************************************************************
C
C
C     EXTRACT 4 OCTETS CONTAINING ASCII G R I B .
C
      KWORD = 1
      IOFF  = 0
      IBYTE = 8
      INVAL  = 4
C
      CALL GSBYTE_MF (KGRIB(1),ILBLOCK(1),IOFF,IBYTE,0,INVAL,
     S             KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
      IBLOCK=ILBLOCK
C
      IF (KERR.NE.0) THEN
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ELSEIF (IBLOCK(1).NE.71.AND.IBLOCK(2).NE.82.AND.
     S        IBLOCK(3).NE.73.AND.IBLOCK(4).NE.66) THEN
        KERR=3
        WRITE (UNIT=*,FMT=*) 'NO ''GRIB'' GROUP (BLOCK 0) FOUND'
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C
C
C
C
C
C
C********************************************************************
C*
C*    BLOCK 1 - PRODUCT DEFINITION BLOCK.
C*
C********************************************************************
C
C
C     EXTRACT FIRST 4 OCTETS OF BLOCK 1.
C
      INVAL  = 4
C
      CALL GSBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,IBYTE,0,INVAL,
     S             KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
      IBLOCK=ILBLOCK
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C     IF THESE OCTETS DO NOT CONTAIN 0, 0, 24 AND 0 RESPECTIVELY
C     THEN THE DATA IS IN THE OLD VERSION OF THE CODE.
C
      IF (IBLOCK(1).EQ.0.AND.IBLOCK(2).EQ.0.AND.IBLOCK(3).EQ.24.
     C    AND.IBLOCK(4).EQ.0)
     C      THEN
                ISNEW = 1
                INVAL = 20
            ELSE
                ISNEW = 0
                INVAL = 16
            ENDIF
C
C     EXTRACT NEXT 16 OR 20 OCTETS OF BLOCK 1 - NUMBER DEPENDS
C     ON VERSION OF CODE.
C
      CALL GSBYTE_MF (KGRIB(KWORD),ILBLOCK(5),IOFF,IBYTE,0,INVAL,
     S             KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
      IBLOCK=ILBLOCK
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C     TRANSFER PRODUCT DEFINITION INFORMATION TO OUTPUT ARRAY.
C     EXCEPT FOR NUMBER AVERAGED/ACCUMULATED AND MISSING FROM AVER/ACCUM
C
      DO 100 J=1,17
         KB1PAR(J) = IBLOCK(J+4*ISNEW)
  100 CONTINUE
C
C
C
C
C     NUMBER AVERAGED OCCUPIES 2 OCTETS.
C
      KB1PAR(18) = IBLOCK(18+4*ISNEW) * 256 + IBLOCK(19+4*ISNEW)
C
C     NUMBER MISSING FROM AVERAGES/ACCUMULATIONS.
C
      KB1PAR(19)=IBLOCK(20+4*ISNEW)
C
C
C***
C*    DESCRIPTION OF LEVEL OR LAYER ( CODE TABLE 3 ).
C***
C
C     CERTAIN LEVEL TYPES REQUIRE THAT THE DESCRIPTION OCCUPY
C     BOTH 8 BIT FIELDS. PUT SINGLE VALUE IN FIRST LEVEL WORD.
C
      ITEMP = 32 + (ISNEW*224)
C
      IF( (KB1PAR(6).EQ. 20).OR.
     X    (KB1PAR(6).EQ.100).OR.
     X    (KB1PAR(6).EQ.103).OR.
     X    (KB1PAR(6).EQ.105).OR.
     X    (KB1PAR(6).EQ.107).OR.
     X    (KB1PAR(6).EQ.109).OR.
     X    (KB1PAR(6).EQ.111).OR.
     X    (KB1PAR(6).EQ.113).OR.
     X    (KB1PAR(6).EQ.115).OR.
     X    (KB1PAR(6).EQ.117).OR.
     X    (KB1PAR(6).EQ.125).OR.
     X    (KB1PAR(6).EQ.127).OR.
     X    (KB1PAR(6).EQ.160).OR.
     X    (KB1PAR(6).EQ.210) ) THEN
C
              KB1PAR(7) = KB1PAR(7) * ITEMP + KB1PAR(8)
              KB1PAR(8) = 0
          ENDIF
C
C
C
C
C
C
C
C
C
C***
C*    TIME RANGE.
C***
C
C     ONE TIME RANGE CAN OCCUPY TWO OCTETS.
C
      IF (KB1PAR(17).EQ.10)
     C   THEN
             KB1PAR(15) = KB1PAR(15) * ITEMP + KB1PAR(16)
             KB1PAR(16) = 0
         ENDIF
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
C
C********************************************************************
C*
C*    BLOCK 2 - GRID DESCRIPTION BLOCK.
C*
C********************************************************************
C
C***
C*    CHECK FLAG INDICATING IF BLOCK 2 IS INCLUDED IN CODE.
C***
C
C     INC IS SET TO 1 , IF BLOCK 2 IS INCLUDED.
C
      INC = 0
C
C
      IF (ISNEW.EQ.0)
     C   THEN
C            OLD VERSION OF CODE
C
C            BLOCKS INCLUDED    BINARY VALUE      DECIMAL VALUE
C
C                NONE             00000000               0
C                 2               00000001               1
C                 3               00000010               2
C             2 AND 3             00000011               3
C
             IF (KB1PAR(4).EQ.1.OR.KB1PAR(4).EQ.3) INC = 1
         ELSE
C            NEW VERSION
C
C            BLOCKS INCLUDED    BINARY VALUE      DECIMAL VALUE
C
C                NONE             00000000               0
C                 2               10000000             128
C                 3               01000000              64
C             2 AND 3             11000000             192
C
             IF (KB1PAR(4).EQ.128.OR.KB1PAR(4).EQ.192) INC = 1
         ENDIF
C
      IF (INC.EQ.1)
     C   THEN
C
C
C***
C*           LENGTH OF GRID DESCRIPTION BLOCK.
C***
C
             CALL GBYTE_MF (KGRIB(KWORD),ILLBLK,IOFF,24)
             ILBLK=ILLBLK
             CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C
C
C***
C*           NUMBER OF UNUSED BITS AT END OF BLOCK.
C***
C
             CALL GBYTE_MF (KGRIB(KWORD),ILLNIL,IOFF,8)
             ILNIL=ILLNIL
             CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C
C
C
C***
C*           NEXT OCTET IS RESERVED.
C***
C
             CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C
C
C
C
C
C
C***
C*           DATA REPRESENTATION TYPE.
C***
C
             CALL GBYTE_MF (KGRIB(KWORD),ILB2PAR(1),IOFF,8)
             KB2PAR=ILB2PAR
             CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C            LAT/LONGITUDE GRID, GAUSSIAN GRID AND SPHERICAL HARMONICS
C            ARE THE ONLY DATA REPRESENTATIONS HANDLED.
C
             IF (KB2PAR(1).NE.0.AND.KB2PAR(1).NE.4.AND.KB2PAR(1).NE.50
     C           .AND.KB2PAR(1).NE.80)
     C          THEN
                     WRITE (*,*)'GRID DESCRIPTION BLOCK NOT YET DEFINED'
                     KERR = -4
                     IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                     RETURN
                ENDIF
C
C
C
C
C***
C*           LAT/LONG OR GAUSSIAN GRID.
C***
C
             IF (KB2PAR(1).EQ.0.OR.KB2PAR(1).EQ.4)
     C         THEN
C
C                  NUMBER OF LAT/LONG POINTS.
C
                   CALL GSBYTE_MF (KGRIB(KWORD),ILB2PAR(2),IOFF,16,0,2,
     S                          KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
                   KB2PAR=ILB2PAR
                   IF (KERR.NE.0)  THEN 
                     IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                     RETURN
                   ENDIF
C
C                  LAT/LONG OF ORIGIN.
C
                   CALL GSBYTE_MF (KGRIB(KWORD),ILLAT(1),IOFF,24,0,2,
     S                          KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
                   ILAT=ILLAT
                   IF (KERR.NE.0)  THEN 
                     IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                     RETURN
                   ENDIF
C
C                  IF SIGN BIT SET TO 1 , VALUES ARE NEGATIVE.
C
                  DO 200 J=1,2
C
                  IF (ILAT(J).LE.2**23) THEN
                    KB2PAR(J+3) = ILAT(J)
                  ELSE
                    KB2PAR(J+3) = 2**23 - ILAT(J)
                  ENDIF
C
  200             CONTINUE
C
C
C                  RESOLUTION FLAG.
C
                   CALL GBYTE_MF (KGRIB(KWORD),ILB2PAR(6),IOFF,8)
                   KB2PAR=ILB2PAR
                   CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
                   IF (KERR.NE.0)  THEN 
                     IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                     RETURN
                   ENDIF
C
C
C                  LAT/LONG OF EXTREME POINTS.
C
                   CALL GSBYTE_MF (KGRIB(KWORD),ILLAT(1),IOFF,24,0,2,
     S                          KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
                   ILAT=ILLAT
                   IF (KERR.NE.0)  THEN 
                     IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                     RETURN
                   ENDIF
C
C                  IF SIGN BIT SET TO 1 , VALUES ARE NEGATIVE.
C
                  DO 300 J=1,2
C
                  IF (ILAT(J).LE.2**23) THEN
                    KB2PAR(J+6) = ILAT(J)
                  ELSE
                    KB2PAR(J+6) = 2**23 - ILAT(J)
                  ENDIF
C
  300             CONTINUE
C
C                 DIRECTION INCREMENTS / NUMBER OF LATITUDE LINES.
C
                  CALL GSBYTE_MF (KGRIB(KWORD),ILB2PAR(9),IOFF,16,0,2,
     S                         KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
                  KB2PAR=ILB2PAR
                  IF (KERR.NE.0)  THEN 
                    IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                    RETURN
                  ENDIF
C
C                 SCANNING MODE FLAGS.
C
                  CALL GBYTE_MF (KGRIB(KWORD),ILB2PAR(11),IOFF,8)
                  KB2PAR=ILB2PAR
                  CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
                  IF (KERR.NE.0)  THEN 
                    IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                    RETURN
                  ENDIF
C
C                 4 RESERVED OCTETS.
C
                  CALL OFFSET_MF (IOFF,4,KWORD,8,KNBIT,KLENG,KERR)
                  IF (KERR.NE.0)  THEN 
                    IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                    RETURN
                  ENDIF
C
              ENDIF
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
C
C
C
C
C***
C*            SPHERICAL HARMONIC DATA.
C***
C
              IF (KB2PAR(1).EQ.50.OR.KB2PAR(1).EQ.80)
     C           THEN
C
C                    PENTAGONAL RESOLUTION PARAMETERS.
C
                     CALL GSBYTE_MF (KGRIB(KWORD),ILB2PAR(2),IOFF,16,0,
     S                            3,KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
                     KB2PAR=ILB2PAR
                     IF (KERR.NE.0)  THEN 
                       IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                       RETURN
                     ENDIF
C
C                    REPRESENTATION TYPE AND MODE.
C
                     CALL GSBYTE_MF (KGRIB(KWORD),ILB2PAR(5),IOFF,8,0,2,
     S                            KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
                     KB2PAR=ILB2PAR
                     IF (KERR.NE.0)  THEN 
                       IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                       RETURN
                     ENDIF
C
C                    18 RESERVED OCTETS.
C
                     CALL OFFSET_MF (IOFF,18,KWORD,8,KNBIT,KLENG,KERR)
                     IF (KERR.NE.0)  THEN 
                       IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                       RETURN
                     ENDIF
C
                 ENDIF
C
C**
C*    STRETCHED AND ROTATED SHPERICAL HARMONICS.
C**
         IF(KB2PAR(1).EQ.80) THEN
C
C*    SOUTHERN POLE OF ROTATED GRID
C
                   CALL GSBYTE_MF (KGRIB(KWORD),ILLAT(1),IOFF,24,0,2,
     S                          KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
                   ILAT=ILLAT
                   IF (KERR.NE.0)  THEN 
                     IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                     RETURN
                   ENDIF
C
C                  IF SIGN BIT SET TO 1 , VALUES ARE NEGATIVE.
C
                  DO 410 J=1,2
C
                  IF (ILAT(J).LE.2**23) THEN
                    KB2PAR(J+11) = ILAT(J)
                  ELSE
                    KB2PAR(J+11) = 2**23 - ILAT(J)
                  ENDIF
C
  410             CONTINUE
C
C*    ANGLE OF ROTATION
C
      CALL GBYTE_MF (KGRIB(KWORD),ILEXP,IOFF,8)
      IEXP=ILEXP
      CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      CALL GBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,24)
      IMANT=ILMANT
      CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
             CALL DECFP_MF (KB2PAR(14),IEXP,IMANT)
C
C*    POLE OF STRETCHING
C
                   CALL GSBYTE_MF (KGRIB(KWORD),ILLAT(1),IOFF,24,0,2,
     S                          KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
                   ILAT=ILLAT
                   IF (KERR.NE.0)  THEN 
                     IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                     RETURN
                   ENDIF
C
C                  IF SIGN BIT SET TO 1 , VALUES ARE NEGATIVE.
C
                  DO 420 J=1,2
C
                  IF (ILAT(J).LE.2**23) THEN
                    KB2PAR(J+14) = ILAT(J)
                  ELSE
                    KB2PAR(J+14) = 2**23 - ILAT(J)
                  ENDIF
C
  420             CONTINUE
C
C*    STRETCHING FACTOR
C
      CALL GBYTE_MF (KGRIB(KWORD),ILEXP,IOFF,8)
      IEXP=ILEXP
      CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      CALL GBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,24)
      IMANT=ILMANT
      CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CALL DECFP_MF (KB2PAR(17),IEXP,IMANT)
C
      ENDIF
C
C
C
C            LENGTH IS 32 OCTETS FOR LAT/LONG, GAUSSIAN AND SPHERICAL
C            HARMONICS, 52 OCTETS FOR STRETCHED AND ROTATED SPHERICAL
C            HARMONICS.FOR ANY DATA  ON HYBRID LEVELS THE
C            VERTICAL COORDINATES ARE ADDED.
C            GET NUMBER OF VERTICAL COORDINATE PARAMETERS, IF ANY.
C
         IF(KB2PAR(1).EQ.80) THEN
            ILEN=52
         ELSE
            ILEN=32
         ENDIF
C
             KJLENV = (ILBLK - ILEN) / 4
C
             IF(KLENV.LT.KJLENV) THEN
                    KERR = -3
                    WRITE (*,9001) KJLENV,KLENV
 9001               FORMAT (1H ,'NUMBER OF VERTICAL COORDINATES - ',I4,
     C                          ', ARRAY SIZE IS - ',I4)
                    IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                    RETURN
                ENDIF
C
C
C
C***
C*               VERTICAL COORDINATE PARAMETERS FOR HYBRID LEVELS.
C***
C
                 IF (KJLENV.NE.0) THEN
C
                    DO 400 J=1,KJLENV
                    CALL GBYTE_MF (KGRIB(KWORD),ILEXP,IOFF,8)
                    IEXP=ILEXP
                    CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
                    IF (KERR.NE.0)  THEN 
                      IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                      RETURN
                    ENDIF
                    CALL GBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,24)
                    IMANT=ILMANT
                    CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
                    IF (KERR.NE.0)  THEN 
                      IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                      RETURN
                    ENDIF
                    CALL DECFP_MF (PVERT(J),IEXP,IMANT)
  400               CONTINUE
C
                 ENDIF
C
         ENDIF
C
C
C
C     RETURN IF ONLY PRODUCT AND GRID DEFINITION BLOCKS REQUIRED.
C
      IF (KLENF.EQ.1)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C
C
C
C
C
C
C
C********************************************************************
C*
C*    BLOCK 3 - BIT MAP BLOCK.
C*
C********************************************************************
C
C     INC IS SET TO 1 , IF BLOCK 3 IS INCLUDED.
C
      INC = 0
C
C
      IF (ISNEW.EQ.0)
     C   THEN
C            OLD VERSION OF CODE
C
C            BLOCKS INCLUDED    BINARY VALUE      DECIMAL VALUE
C
C                NONE             00000000               0
C                 2               00000001               1
C                 3               00000010               2
C             2 AND 3             00000011               3
C
             IF (KB1PAR(4).EQ.2.OR.KB1PAR(4).EQ.3) INC = 1
         ELSE
C            NEW VERSION
C
C            BLOCKS INCLUDED    BINARY VALUE      DECIMAL VALUE
C
C                NONE             00000000               0
C                 2               10000000             128
C                 3               01000000              64
C             2 AND 3             11000000             192
C
             IF (KB1PAR(4).EQ.64.OR.KB1PAR(4).EQ.192) INC = 1
         ENDIF
C
      IF (INC.EQ.1)
     C     THEN
               WRITE (*,*)'BIT MAP BLOCK NOT YET DEFINED'
               KERR = -4
               IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
               RETURN
           ENDIF
C
C
C
C
C**********************************************************************
C*
C*    CONVERT VALUES IN BLOCKS 1 AND 2 TO NEW CODE VALUES.
C*
C**********************************************************************
C
C
      IF (ISNEW.EQ.0)
     C   THEN
C
C            CHANGE CODE FOR FLAG INDICATING INCLUSION OR OMISSION
C            OF BLOCKS 2 AND 3.
C
C            BINARY VALUE OLD CODE = 00000001
C            BINARY VALUE NEW CODE = 10000000 (128 DECIMAL)
C
             IF (KB1PAR(4).EQ.1)
     C          THEN
                    KB1PAR(4) = 128
                ELSE
                    WRITE (*,*) 'DECOGA : BLOCK FLAG ERROR.'
                    KERR = 2
                    IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                    RETURN
                ENDIF
C
C            CHANGE CODE FOR TIME UNIT.
C
             IF (KB1PAR(14).EQ.40)
     C          THEN
                    KB1PAR(14) = 1
                ELSE
                    WRITE (*,*) 'DECOGA : TIME UNIT ERROR.'
                    KERR = 2
                    IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                    RETURN
                ENDIF
C
C            CONVERT FLAGS FOR LAT/LONG AND GAUSSIAN GRID DATA.
C
             IF (KB2PAR(1).EQ.0.OR.KB2PAR(1).EQ.4)
     C          THEN
C
C                   CONVERT SCANNING MODE FLAG.
C
                    IF (KB2PAR(11).EQ.1)
     C                 THEN
                           KB2PAR(11) = 0
                       ELSE
                           WRITE (*,*) 'DECOGA : SCAN MODE FLAG ERROR.'
                           KERR = 2
                           IF (LHOOK) 
     C                       CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                           RETURN
                       ENDIF
C
C                   CONVERT RESOLUTION FLAG.
C
                    IF (KB2PAR(6).EQ.3)
     C                 THEN
                           KB2PAR(6) = 128
                       ELSE
                           WRITE (*,*) 'DECOGA : RESOLUTION FLAG ERROR.'
                           KERR = 2
                           IF (LHOOK) 
     C                       CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                           RETURN
                       ENDIF
                ENDIF
         ENDIF
C
C
C
C
C
C
C********************************************************************
C*
C*    BLOCK 4 - BINARY DATA BLOCK.
C*
C********************************************************************
C
C
C
C***
C*    GET LENGTH OF BINARY DATA BLOCK.
C***
C
C
      IPW = KWORD
      IPB = IOFF
C
      CALL GBYTE_MF (KGRIB(KWORD),ILLBIN,IOFF,24)
      ILBIN=ILLBIN
      CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C
C
C
C***
C*    4 BIT FLAG / 4 BIT COUNT OF UNUSED BITS AT END OF BLOCK OCTET.
C***
C
      CALL GBYTE_MF (KGRIB(KWORD),ILFLAG,IOFF,8)
      IFLAG=ILFLAG
      CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      IF (ISNEW.EQ.0)
     C    THEN
C             OLD VERSION OF CODE
C
C             0000---- GRID POINT
C             0001---- SPHERICAL HARMONICS
C
              IREP = IFLAG / 16
              ILNIL = IFLAG - IREP * 16
C
          ELSE
C             NEW VERSION OF CODE
C
C             0000---- GRID POINT
C             1000---- SPHERICAL HARMONICS
C             1100----       -''-           COMPLEX PACKING
C
              IREP = IFLAG / 128
              ILNIL = MOD(IFLAG,16)
C
          ENDIF
C
C
C
C
C
C
C
C
C***
C*    GET SCALE FACTOR
C***
C
      CALL GBYTE_MF (KGRIB(KWORD),ILSCALX,IOFF,16)
      ISCALX=ILSCALX
      CALL OFFSET_MF (IOFF,1,KWORD,16,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      IF (ISCALX.LE.2**15) THEN
        ISCALE = ISCALX
      ELSE
        ISCALE = 2**15 - ISCALX
      ENDIF
C
C
C
C
C***
C*    GET REFERENCE VALUE (PMIN) IN GRIB FORMAT (IEXP,IMANT)
C***
C
C
      CALL GBYTE_MF (KGRIB(KWORD),ILEXP,IOFF,8)
      IEXP=ILEXP
      CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      CALL GBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,24)
      IMANT=ILMANT
      CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C
C     CHECK FOR MISSING DATA INDICATORS.
C
      IMISS = 0
      IF (ISCALX.EQ.65535.AND.IEXP.EQ.255.AND.IMANT.EQ.16777215)
     C   IMISS = 1
C
C
C
C***
C*    GET NUMBER OF BITS IN EACH DATA VALUE.
C***
C
      CALL GBYTE_MF (KGRIB(KWORD),ILBITS,IOFF,8)
      KBITS=ILBITS
      CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C
C     CHECK NUMBER OF BITS PER DATA FIELD.
C
      IF (KBITS.GT.KNBIT.OR.KBITS.GT.IMAX)
     C    THEN
              KERR = -1
              WRITE (*,9002) KBITS,KNBIT,IMAX
 9002         FORMAT (1H ,'NUMBER OF BITS PER DATA VALUE, ',I3,
     C                'EXCEEDS WORD LENGTH, ',I3,' OR MAXIMUM ',
     C                ' PERMITTED VALUE, ',I3)
              CALL SDL_SRLABORT
C
      ELSEIF (LDARPE) THEN
C
C     Minimum value given as input argument... it has not been decoded
C     before as stored in GRIB field (should be zero).
C        But the scale factor may be computed, now.
C
        ZSCALE= ( PMAX - PMIN ) / FLOAT ( 2**KBITS-1 )
C
      ELSEIF (IMISS.EQ.0) THEN
C
C
C     CONVERT REFERENCE VALUE AND SCALE FACTOR.
C
C
        CALL DECFP_MF (PMIN,IEXP,IMANT)
        ZSCALE = 2.**ISCALE
        PMAX=PMIN+FLOAT (2**(KBITS+1)-3) * 2.** (ISCALE-1)
C
      ENDIF
C
C**
C*    IF COMPLEX PACKING GET ADDITIONAL PARAMETERS
C**
      IF(IFLAG.GE.192) THEN
C
C     GET START OCTET OF PACKED DATA
C
      CALL GBYTE_MF (KGRIB(KWORD),ILSPDA,IOFF,16)
      ISPDA=ILSPDA
      CALL OFFSET_MF (IOFF,1,KWORD,16,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      IPREMC=1+(8*ISPDA-1)/KNBIT
C
      IF (IPREMC.GT.KLENF) THEN
        KERR=-3
              WRITE (*,9004) KLENF,ISPDA
 9004         FORMAT (' OUTPUT ARRAY LENGTH (',I7,
     S                ' WORDS) TOO SHORT, START OCTET OF PACKED DATA =',
     S                I9)
              CALL SDL_SRLABORT
      ENDIF
C
C     GET LAPLACIAN SCALING FACTOR
C
      CALL GBYTE_MF (KGRIB(KWORD),ILSCALP,IOFF,16)
      KSCALP=ILSCALP
      CALL OFFSET_MF (IOFF,1,KWORD,16,KNBIT,KLENG,KERR)
C
      IF (KERR.NE.0) THEN
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ELSEIF (KSCALP.GT.2**15) THEN
        KSCALP=2**15-KSCALP
      ENDIF
C
C     GET TRUNCATION OF UNPACKED COEFFS.
C
      DO 430 J=1,3
      CALL GBYTE_MF (KGRIB(KWORD),ILCPACK,IOFF,8)
      KCPACK=ILCPACK
      CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
 430  CONTINUE
C
C     COMPUTE NUMBER OF UNPACKED DATA ITEMS
C
      IREP=(KCPACK+1)**2
      ENDIF
C
C***
C*    IF DATA IS IN SPHERICAL HARMONIC FORM, NEXT OCTETS
C*    CONTAIN UNPACKED COEFFICIENT(S).
C***
C
      IF (IREP.NE.0)
     C    THEN
C
C            EVENTUALLY,
C            GET UNPACKED COEFFICIENTS IN GRIB FORMAT AND
C            CONVERT TO FLOATING POINT.
C
        IF (IMISS.EQ.1) THEN
C
          DO 440 J=1,IREP
          PFDATA(J)=0.
  440     CONTINUE
C
        ELSEIF (LDARPE) THEN
C
C       The "unpacked" values should then contain zeroes,
C     and are skipped: the corresponding values of PFDATA array
C     are then not provided. A positioning in KGRIB is made.
C
          CALL OFFSET_MF (IOFF,IREP,KWORD,32,KNBIT,KLENG,KERR)
          IF (KERR.NE.0)  THEN 
            IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
C
        ELSE
C
         DO 450 J=1,IREP
             CALL GBYTE_MF (KGRIB(KWORD),ILEXP,IOFF,8)
             IEXP=ILEXP
             CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
             CALL GBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,24)
             IMANT=ILMANT
             CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
                    CALL DECFP_MF (PFDATA(J),IEXP,IMANT)
 450     CONTINUE
C
        ENDIF
C
        IF (IFLAG.GE.192) THEN
C
C            Proper positioning on start octet of packed data, because
C        it may be not the first octet following the "unpacked" data.
C
          ISSUIV=19+IREP*4
C
          IF (ISSUIV.LT.ISPDA) THEN
            CALL OFFSET_MF (IOFF,ISPDA-ISSUIV,KWORD,8,KNBIT,KLENG,KERR)
            IF (KERR.NE.0)  THEN 
              IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
              RETURN
            ENDIF
          ENDIF
C
        ENDIF
C
      ENDIF
C
C
C
C
C*****
C*    DECODE DATA VALUES TO FLOATING POINT AND STORE IN PFDATA.
C*****
C
C     FIRST CALCULATE THE NUMBER OF DATA VALUES.
C
      KJLENF = ILBIN - 11 - IREP * 4
      IF(IFLAG.GE.192) KJLENF=KJLENF-7-(ISPDA-ISSUIV)
      KJLENF = (KJLENF * 8 - ILNIL) / KBITS
C
C     CHECK LENGTH OF OUTPUT ARRAY.
C
      IF (KJLENF+IREP.GT.KLENF)
     C   THEN
             KERR = -3
             WRITE (*,9003) KJLENF,KLENF
 9003        FORMAT (1H ,'NUMBER OF VALUES TO BE DECODED IS - ',I7,
     C                    ', ARRAY SIZE - ',I7)
             CALL SDL_SRLABORT
         ENDIF
C
      IF (IMISS.EQ.0.AND.ZSCALE.GT.0.) THEN
        CALL GSBYTE_MF (KGRIB(KWORD),PFDATA(IREP+1),IOFF,KBITS,0,KJLENF,
     S               KNBIT,'D',KLENG,KERR,KWORD,.TRUE.)
        IF (KERR.NE.0)  THEN 
          IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
          RETURN
        ENDIF
        CALL UNPAGB (PFDATA(IREP+1),PFDATA(IREP+1),PMIN,PMAX,KBITS,
     S               ZSCALE,KJLENF,LDARPE)
      ELSE
C
        IF (ZSCALE.LE.0.) THEN
C
          DO 460 J= IREP+1,IREP+KJLENF
          PFDATA(J) = PMIN
  460     CONTINUE
C
        ELSE
C
          DO 470 J= IREP+1,IREP+KJLENF
          PFDATA(J) = 0.
  470     CONTINUE
C
        ENDIF
C
        CALL OFFSET_MF (IOFF,KJLENF,KWORD,KBITS,KNBIT,KLENG,KERR)
        IF (KERR.NE.0)  THEN 
          IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
          RETURN
        ENDIF
      ENDIF
C
C
C     INCLUDE UNPACKED COEFFICIENT IN COUNT, IF SPHERICAL HARMONIC DATA.
C
      KJLENF = IREP + KJLENF
C
C
C
C
C***
C*    SKIP ANY ZERO FILL AT END OF DATA BLOCK.
C***
C
      IL = (KWORD-IPW) * KNBIT + IOFF -IPB
      IL = ILBIN * 8 - IL
C
      IF (IL.NE.0) THEN
              CALL OFFSET_MF (IOFF,1,KWORD,IL,KNBIT,KLENG,KERR)
              IF (KERR.NE.0)  THEN 
                IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
                RETURN
              ENDIF
      ENDIF
C
C
C
C
C********************************************************************
C*
C*    BLOCK 5 - END BLOCK.
C*
C********************************************************************
C
C     CHECK 7 7 7 7 AT END OF CODED DATA.
C
C
      INVAL=4
      ILBVAL=8
      CALL GSBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,ILBVAL,0,INVAL,
     S             KNBIT,'D',KLENG,KERR,KWORD,.FALSE.)
      IBLOCK=ILBLOCK
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C     55 IS THE VALUE OF ASCII CHARACTER '7'
C
      DO 600 J=1,INVAL
         IF (IBLOCK(J).NE.55) KERR = 1
  600 CONTINUE
C
      IF (KERR.EQ.1) WRITE (*,*) ' NO 7777 GROUP FOUND '
C
C
C
C     SET NUMBER OF VALUES NEGATIVE, IF MISSING DATA VALUES.
C
      IF (IMISS.EQ.1) KJLENF = - KJLENF
C
      IF (LHOOK) CALL DR_HOOK('DECOGA',1,ZHOOK_HANDLE)
      END
