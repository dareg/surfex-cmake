      SUBROUTINE CODEGA (PFDATA,KLENF,KBITS,KNBIT,KB1PAR,
     C                   KB2PAR,PVERT,KLENV,KGRIB,KLENG,KWORD,
     C                   KROUND,KCPACK,KSCALP,KERR,PMIN,PMAX,LDARPE)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C
C
C********************************************************************
C*
C*     NAME     : CODEGA
C*
C*     FUNCTION : CODE ARRAY OF FLOATING POINT VALUES
C*                IN WMO GRIB CODE.
C*
C*     INPUT    : PARAMETERS FOR BLOCK 1 , PRODUCT DEFINITION BLOCK.
C*                --------------------------------------------------
C*                KB1PAR - INTEGER ARRAY (DIMENSION 19) OF PRODUCT
C*                         DEFINITION INFORMATION.
C*
C*                WORD      CONTENTS
C*                ----      --------
C*                  1       ORIGINATING CENTRE IDENTIFIER.
C*                  2       MODEL IDENTIFICATION.
C*                  3       GRID DEFINITION.
C*                  4       FLAG ( CODE TABLE 1 ).
C*                  5       PARAMETER IDENTIFIER (CODE TABLE 2).
C*                  6       TYPE OF LEVEL (CODE TABLE 3).
C*                 7-8      VALUE(S) OF LEVEL (CODE TABLE 3).
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
C*                          VALUE(S) OF LEVEL CAN OCCUPY 2 WORDS.
C*                          FOR A LAYER THE FIRST WORD DEFINES THE TOP
C*                          AND THE SECOND THE BOTTOM OF THE LAYER.
C*                          FOR A SINGLE LEVEL, ONLY THE FIRST WORD IS
C*                          USED.
C*
C*                PARAMETERS FOR BLOCK 2 , GRID DEFINITION BLOCK.
C*                -----------------------------------------------
C*                KB2PAR - INTEGER ARRAY (DIMENSION 17) CONTAINING GRID
C*                         DEFINITION INFORMATION.
C*                         USE VARIES WITH REPRESENTATION TYPE.
C*
C*                WORD      LAT/LONG GRID
C*                ----      -------------
C*                 1        DATA REPRESENTATION TYPE (CODE TABLE 6)
C*                 2        NO. OF POINTS ALONG A LATITUDE
C*                 3        NO. OF POINTS ALONG A MERIDIAN
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
C*                TO BE DEFINED LATER.
C*
C*                PVERT  - REAL ARRAY OF VERTICAL COORDINATE PARAMETERS
C*                KLENV  - LENGTH OF THIS ARRAY (NUMBER OF COORDINATE
C*                         PARAMETERS).
C*
C*                USED FOR HYBRID LEVELS ONLY.
C*
C*
C*                PARAMETERS FOR BLOCK 3 , BIT-MAP DEFINITION BLOCK.
C*                --------------------------------------------------
C*
C*                TO BE DEFINED LATER.
C*
C*                PARAMETERS FOR BLOCK 4 , BINARY DATA BLOCK.
C*                -------------------------------------------
C*                PFDATA  - ARRAY OF FLOATING POINT VALUES.
C*
C*                          ****
C*                              VALUES IN THIS ARRAY ARE OVERWRITTEN
C*                              BY THE SUBROUTINE.
C*                                                              ****
C*
C*                KLENF   - LENGTH OF THIS ARRAY (NUMBER OF DATA
C*                          VALUES TO BE CODED)
C*                KBITS   - NUMBER OF BITS FOR CODED DATA VALUES.
C*
C*                KCPACK  - KCPACK GREATER THAN ZERO INDICATES COMPLEX
C*                          PACKING, WITH A *TRIANGULAR* SUB-TRUNCATION
C*                          "UNPACKED" OF ORDER KCPACK.
C
C*                KSCALP  - FOR COMPLEX PACKING OF SPECTRAL COEFFICIENTS
C*                          ONLY, LAPLACIAN SCALING FACTOR.
C*
C*                          TO INDICATE MISSING DATA VALUES, KLENF IS
C*                          MADE NEGATIVE AND ALL DATA VALUES SET TO 0.
C*
C*                OTHER PARAMETERS
C*                ----------------
C*
C*                KNBIT    - NUMBER OF BITS IN COMPUTER WORD.
C*
C*                KGRIB   - ARRAY TO RECEIVE DATA IN GRIB CODE.
C*                KLENG   - LENGTH OF THIS ARRAY.
C*
C*                KROUND  - 0 , NO ROUNDING.
C*                          1 , ROUND TO NEAREST 120 OCTETS.
C*
C*                LDARPE  - .TRUE., modifications for ARPEGE coding.
C*                          .FALSE., no such modifications.
C*
C*    OUTPUT    : KGRIB   - DATA CODED IN GRIB CODE.
C*                KWORD   - NUMBER OF WORDS OCCUPIED BY CODED DATA.
C*
C*                KERR    - ERROR INDICATOR.
C*
C*                           0, NO ERROR.
C*
C*                          -1, NUMBER OF BITS PER DATA VALUE EXCEEDS
C*                              WORD LENGTH, OR MAXIMUM ALLOWED.
C*                          -2, OUTPUT ARRAY TOO SMALL TO HOLD CODED
C*                              DATA VALUES.
C*                          -3, CURRENTLY UNDEFINED OPTION SELECTED.
C*                          -4, MISSING DATA FIELD CONTAINS NON-ZERO.
C*                          -5, INVALID COMPLEX PACKING MODE.
C*                          -6, INVALID LAPLACIAN SCALING FACTOR.
C*
C*                           1, INVALID ORIGINATING CENTRE.
C*                           2, INVALID MODEL IDENTIFICATION.
C*                           3, INVALID GRID DEFINITION.
C*                           4, INVALID BLOCK INDICATOR FLAG
C*                           5, INVALID PARAMETER IDENTIFIER.
C*                           6, INVALID LEVEL TYPE INDICATOR.
C*                           7, ERROR IN DESCRIPTION OF LEVELS.
C*                           8, DATE/TIME ERROR.
C*                           9, INVALID TIME UNIT OR RANGE.
C*                          10, INVALID DATA REPRESENTATION TYPE.
C*                          11, INVALID EXTREME POINT/INCREMENT FLAG.
C*                          12, INVALID SCANNING MODE FLAGS.
C*                          13, INVALID COMPLEX PACKING MODE.
C*
C*                PMIN    - MINIMUM VALUE OF THE FIELD DATA.
C*
C*                PMAX    - MAXIMUM VALUE OF THE FIELD DATA.
C*
C*    JOHN HENNESSY , ECMWF , NOVEMBER 1985
C*
C*    Modifications by Mats HAMRUD, ECMWF, to handle complex packing of
C*    spectral harmonics data, 1988.
C*
C*    Modifications by Jean CLOCHARD, French DMN, January 1990:
C*
C*    -  to remove automatic arrays introduced by the last modification
C*     ( *** BUT INPUT ARRAY PFDATA IS OVERWRITTEN AGAIN *** );
C*    -  to comply with the "DOCTOR" norm;
C*    -  to use (standard case) as reference value the highest 32-bit
C*       value not exceeding the minimum of the field, this reference
C*       value being used to compute departures from this value;
C*    -  to allow modifications of the packing for ARPEGE files use,
C*       in an optional way.
C*
C********************************************************************
C
C
C
C
C
C
#include "precision.h"
C
      INTEGER KLENF, KBITS, KNBIT, KLENG, KWORD, KROUND, KCPACK, KSCALP
      INTEGER KERR, KLENV
C
      REAL (KIND=JPDBLR) PFDATA(*), PVERT(KLENV)
C
      REAL (KIND=JPDBLR) PMIN, PMAX
C
      LOGICAL LDARPE
C
      INTEGER (KIND=JPDBLE) KGRIB(KLENG)
      INTEGER KB1PAR(19), KB2PAR(17)
C
      INTEGER IMAX, ILENF, IMISS, J, IBYTE, INVAL, IOFF, ITEMP, IERR
      INTEGER IERY, IERD, IERM, IERH, IERN, I, IEXP, IMANT, IPW, IPB
      INTEGER IREP, IFLAG, ICPACK, ILEN, ISCALE, ISCALX, ISTPA, IRESTE
      INTEGER ISCALP, IAUXIL, ILENFM, ILBIN, IL, ILNIL, IBITS, INUMBI
      INTEGER (KIND=JPDBLE) IL8, ILEXP, ILMANT, ILFLAG, ILSCALX
      INTEGER (KIND=JPDBLE) ILAUXIL, ILLBIN, ILBITS, ILSTPA, ILSCALP
C
      INTEGER IBLOCK(24)
      INTEGER ILAT(2)
      INTEGER (KIND=JPDBLE) ILBLOCK(24), ILB2PAR(17), ILLAT(2)
C
      REAL (KIND=JPDBLR) ZCOEFF, ZEPSIL, ZAUXIL, ZSCALE, ZS, ZREFER
      REAL (KIND=JPDBLR) ZAUXI2
C
C     SET MAXIMUM NUMBER OF BITS PER DATA FIELD.
C
      SAVE IMAX
      DATA IMAX /60/
C
C     CLEAR ERROR INDICATOR.
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CODEGA',0,ZHOOK_HANDLE)
      KERR = 0
C
C
C
C
C
C********************************************************************
C*
C*    CHECK FOR MISSING DATA INDICATORS.
C*
C********************************************************************
C
C
      IF (KLENF.LT.0)
     C   THEN
             ILENF = - KLENF
             IMISS = 1
C
             DO 50 J=1,ILENF
C
                IF (PFDATA(J).NE.0.) THEN
                       KERR = -4
                       WRITE (*,9012)
 9012                  FORMAT (TR1,'NON-ZERO VALUE IN MISSING DATA',
     C                          ' FIELD.')
                       IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
                       RETURN
                   ENDIF
   50        CONTINUE
C
         ELSE
             ILENF = KLENF
             IMISS = 0
         ENDIF
C
C
C
C
C********************************************************************
C*
C*    CHECK NUMBER OF BITS PER DATA FIELD.
C*
C********************************************************************
C
C
      IF (KBITS.GT.KNBIT.OR.KBITS.GT.IMAX)
     C    THEN
              KERR = -1
              WRITE (*,9000) KBITS,KNBIT,IMAX
 9000         FORMAT (TR1,'NUMBER OF BITS PER DATA VALUE, ',I3,
     C                'EXCEEDS WORD LENGTH, ',I3,' OR MAXIMUM ',
     C                ' PERMITTED VALUE, ',I3)
              IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C
C***
C*    THE LETTERS G R I B ARE INSERTED IN THE
C*    FIRST 4 ELEMENTS OF 'IBLOCK'.
C***
C
      IBLOCK(1) = 71
      IBLOCK(2) = 82
      IBLOCK(3) = 73
      IBLOCK(4) = 66
C
      KWORD = 1
      IBYTE = 8
      INVAL  = 4
      IOFF  = 0
C
C     PACK 4 8-BIT FIELDS IN BLOCK 0 OF CODED ARRAY.
C
      ILBLOCK=IBLOCK
      CALL GSBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,IBYTE,0,INVAL,
     S             KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C*    BLOCK 1 - PRODUCT DEFINITION BLOCK.
C*
C********************************************************************
C
C
C
C***
C*    LENGTH OF BLOCK 1 (IN OCTETS) IN FIRST 3 OCTETS.
C***
C
      IBLOCK(1) = 0
      IBLOCK(2) = 0
      IBLOCK(3) = 24
C
C
C
C
C***
C*    FOURTH OCTET RESERVED AND SET TO 0.
C***
C
      IBLOCK(4) = 0
C
C
C
C
C
C
C***
C*    IDENTIFICATION OF CENTRE
C***
      IBLOCK(5) = KB1PAR(1)
      IF (KB1PAR(1).LT.1.OR.KB1PAR(1).GT.98)
     C   THEN
             KERR = 1
             WRITE (*,9001) KB1PAR(1)
 9001        FORMAT (TR1,'INVALID ORIGINATING CENTRE ',I3)
             IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C*    MODEL IDENTIFICATION
C***
C
      IBLOCK(6) = KB1PAR(2)
      IF (KB1PAR(2).LT.1.OR.KB1PAR(2).GT.255)
     C   THEN
             KERR = 2
             WRITE (*,9002) KB1PAR(2)
 9002        FORMAT (TR1,'INVALID MODEL IDENTIFICATION ',I4)
             IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C***
C*    GRID DEFINITION
C***
C
      IBLOCK(7) = KB1PAR(3)
      IF (KB1PAR(3).LT.1.OR.KB1PAR(3).GT.255)
     C   THEN
             KERR = 3
             WRITE (*,9003) KB1PAR(3)
 9003        FORMAT (TR1,'INVALID GRID IDENTIFICATION ',I4)
             IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C***
C*    OPTIONAL BLOCK INDICATOR FLAGS ( CODE TABLE 1 ).
C***
C
C     BLOCKS INCLUDED    BINARY VALUE      DECIMAL VALUE
C
C         NONE             00000000               0
C          2               10000000             128
C          3               01000000              64
C         2 AND 3          11000000             192
C
      IBLOCK(8) = KB1PAR(4)
C
C     SHIFT 2-BIT FLAG FIELD RIGHT BY 6 BITS TO FACILITATE CHECKING
C
      ITEMP = KB1PAR(4) / 64
C
      IF (ITEMP.LT.0.OR.ITEMP.GT.3)    KERR = 4
C
C     CHECK IF GRID SPECIFIED IS UNCATALOGUED.
C
      IF (KB1PAR(3).EQ.255.AND.KB1PAR(4).EQ.0)  KERR = 4
      IF (KB1PAR(3).EQ.255.AND.KB1PAR(4).EQ.64) KERR = 4
C
      IF (KERR.NE.0) THEN
           CALL PRTBIN_MF (KB1PAR(4),8,ITEMP,IERR)
           WRITE (*,9004) ITEMP
 9004      FORMAT (TR1,'INVALID BLOCK INDICATOR FLAG ',I8.8)
           IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
           RETURN
      ENDIF
C
C
C
C
C***
C*    PARAMETER IDENTIFIER ( CODE TABLE 2 ).
C***
C
      IBLOCK(9) = KB1PAR(5)
      IF (KB1PAR(5).LT.0.OR.KB1PAR(5).GT.255)
     C   THEN
             KERR = 5
             WRITE (*,9005) KB1PAR(5)
 9005        FORMAT (TR1,'INVALID PARAMETER ',I4)
             IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
             RETURN
         ENDIF
C
C
C
C
C***
C*    INDICATOR OF LEVEL TYPE ( CODE TABLE 3 ).
C***
C
      IBLOCK(10) = KB1PAR(6)
C
C*    CODE FIGURES FOR LEVEL TYPE RANGE FROM 0 - 210, WITH 200
C     BEING USED FOR PSEUDO-LEVELS.
C
      IF (KB1PAR(6).LT.0.OR.KB1PAR(6).GT.210) KERR = 6
      IF (KERR.EQ.6)
     C   THEN
             WRITE (*,9006) KB1PAR(6)
 9006        FORMAT (TR1,'INVALID LEVEL TYPE ',I4)
             IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C*    DESCRIPTION OF LEVEL OR LAYER ( CODE TABLE 3 ).
C***
C
C
C
C     FOR CERTAIN SPECIAL LEVELS NO FURTHER DESCRIPTION IS
C     NECESSARY, AND THE FIELDS ARE SET TO 0.
C
      IF (KB1PAR(6).LT.100.OR.KB1PAR(6).EQ.102)
     C    THEN
              KB1PAR(7) = 0
              KB1PAR(8) = 0
          ENDIF
C
      IBLOCK(11) = KB1PAR(7)
      IBLOCK(12) = KB1PAR(8)
C
C     CERTAIN LEVEL TYPES REQUIRE THAT THE DESCRIPTION OCCUPY
C     BOTH 8 BIT FIELDS.
C
      IF( (KB1PAR(6).NE. 20).AND.
     X    (KB1PAR(6).NE.100).AND.
     X    (KB1PAR(6).NE.103).AND.
     X    (KB1PAR(6).NE.105).AND.
     X    (KB1PAR(6).NE.107).AND.
     X    (KB1PAR(6).NE.109).AND.
     X    (KB1PAR(6).NE.111).AND.
     X    (KB1PAR(6).NE.113).AND.
     X    (KB1PAR(6).NE.115).AND.
     X    (KB1PAR(6).NE.117).AND.
     X    (KB1PAR(6).NE.125).AND.
     X    (KB1PAR(6).NE.127).AND.
     X    (KB1PAR(6).NE.160).AND.
     X    (KB1PAR(6).NE.210) ) THEN
C
              IF (KB1PAR(7).GT.255.OR.KB1PAR(8).GT.255) KERR = 7
          ELSE
C
C             8 LOW ORDER BITS IN IBLOCK(12)
C             HIGH ORDER BITS IN IBLOCK(11)
C
              IBLOCK(12) = IBLOCK(11)
              IBLOCK(11) = IBLOCK(11) / 256
              IBLOCK(12) = IBLOCK(12) - IBLOCK(11) * 256
C
C             CHECK MAXIMUM POSSIBLE IN 16 BITS.
C
              IF (KB1PAR(7).GT.65535) KERR = 7
          ENDIF
C
      IF (KERR.NE.0)
     C    THEN
              WRITE (*,9007) KB1PAR(7),KB1PAR(8)
 9007         FORMAT (TR1,'LEVEL DESCRIPTION ERROR ',I8,3X,I8)
              IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C*    DATE AND TIME. CHECK VALID RANGES AND MISSING DATA VALUES.
C***
C
      IBLOCK(13) = KB1PAR(9)
      IERY = 0
      IF (KB1PAR(9).LT.0.OR.KB1PAR(9).GT.99) IERY = 1
      IF (KB1PAR(9).EQ.255) IERY = 0
      IBLOCK(14) = KB1PAR(10)
      IERM = 0
      IF (KB1PAR(10).LT.1.OR.KB1PAR(10).GT.12) IERM = 1
      IF (KB1PAR(10).EQ.255) IERM = 0
      IBLOCK(15) = KB1PAR(11)
      IERD = 0
      IF (KB1PAR(11).LT.1.OR.KB1PAR(11).GT.31) IERD = 1
      IF (KB1PAR(11).EQ.255) IERD = 0
      IBLOCK(16) = KB1PAR(12)
      IERH = 0
      IF (KB1PAR(12).LT.0.OR.KB1PAR(12).GT.23) IERH = 1
      IF (KB1PAR(12).EQ.255) IERH = 0
      IBLOCK(17) = KB1PAR(13)
      IERN = 0
      IF (KB1PAR(13).LT.0.OR.KB1PAR(13).GT.59) IERN = 1
      IF (KB1PAR(13).EQ.255) IERN = 0
C
      KERR = IERY + IERM + IERD + IERH + IERN
C
      IF (KERR.NE.0)
     C    THEN
              KERR = 8
              WRITE (*,9008) KB1PAR(9),KB1PAR(10),KB1PAR(11),KB1PAR(12),
     C                       KB1PAR(13)
 9008         FORMAT (TR1,'INVALID DATE/TIME ',3I2,' / ',2I2)
              IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C***
C*    UNIT OF TIME - ( CODE TABLE 4 ).
C*    TIME RANGE(S) AND TIME RANGE INDICATOR ( CODE TABLE 5 ).
C***
C
C     TIME UNIT.
C
      IBLOCK(18) = KB1PAR(14)
      IF (KB1PAR(14).LT.0.OR.KB1PAR(14).GT.7) KERR = 9
C
C*    TIME UNIT CODES ARE IN THE RANGE 0 - 7, WITH 254 USED FOR SECONDS.
C
      IF (KB1PAR(14).EQ.254) KERR = 0
C
C     THERE CAN BE ONE OR TWO TIME RANGES, EACH IN 1 OCTECT OR
C     ONE TIME RANGE IN 2 OCTECTS.
C
      IF (KB1PAR(17).EQ.0.OR.KB1PAR(17).EQ.1
     C .OR.KB1PAR(17).EQ.10)   KB1PAR(16)  = 0
C
C     CROSS CHECK RANGE VALUES AND FLAG VALUES.
C
      IF (KB1PAR(15).GT.65535.OR.KB1PAR(15).LT.0)   KERR = 9
      IF (KB1PAR(15).GT.255.AND.KB1PAR(17).NE.10)   KERR = 9
      IF (KB1PAR(16).GT.255.OR.KB1PAR(16).LT.0)     KERR = 9
      IF (KB1PAR(17).LT.0.OR.KB1PAR(17).GT.10)      KERR = 9
      IF (KERR.NE.0)
     C    THEN
              WRITE (*,9009) KB1PAR(14),KB1PAR(15),KB1PAR(16),KB1PAR(17)
 9009         FORMAT (TR1,'TIME UNIT/TIME 1/TIME 2/INDICATOR ERROR - ',
     C                I8,2X,'/',I8,2X,'/',I8,2X,'/',I8)
              IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
              RETURN
          ENDIF
C
      IBLOCK(19) = KB1PAR(15)
      IBLOCK(20) = KB1PAR(16)
C
C     ONE TIME RANGE OCCUPYING BOTH OCTETS.
C
      IF (KB1PAR(17).EQ.10)
     C   THEN
C            8 LOW ORDER BITS IN IBLOCK(20)
C            HIGH ORDER BITS IN IBLOCK(19)
C
             IBLOCK(20) = IBLOCK(19)
             IBLOCK(19) = IBLOCK(19) / 256
             IBLOCK(20) = IBLOCK(20) - IBLOCK(19) * 256
         ENDIF
C
C     TIME RANGE INDICATOR
C
      IBLOCK(21) = KB1PAR(17)
C
C
C
C
C
C
C
C***
C*    NUMBER AVERAGED OR ACCUMULATED.
C***
C
      IBLOCK(22) = 0
      IBLOCK(23) = 0
C
C
C     CHECK CONSISTENCY INDICATOR AND NUMBER FIELD.
C
      IF (KB1PAR(17).EQ.3.AND.KB1PAR(18).EQ.0)
     C   THEN
             KERR = 12
             WRITE (*,9013) KB1PAR(17),KB1PAR(18)
 9013        FORMAT (TR1,'INDICATOR/NUMBER AVERAGED ERROR - ',
     C                I8,2X,'/',I8)
             IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
             RETURN
         ENDIF
C
C     VALUE IN TWO OCTETS.
C
      IBLOCK(22) = KB1PAR(18)
C
C     8 LOW ORDER BITS IN IBLOCK(23)
C     HIGH ORDER BITS IN IBLOCK(22)
C
      IBLOCK(23) = IBLOCK(22)
      IBLOCK(22) = IBLOCK(22) / 256
      IBLOCK(23) = IBLOCK(23) - IBLOCK(22) * 256
C
C
C
C
C
C
C***
C*    NUMBER MISSING FROM AVERAGES/ACCUMULATIONS.
C***
C
      IBLOCK(24) = KB1PAR(19)
C
C
C
C
C
C
C
C***
C*
C*    PACK 24 8-BIT FIELDS OF BLOCK 1 IN CODED ARRAY.
C*
C***
C
      IBYTE = 8
      INVAL = 24
C
      ILBLOCK=IBLOCK
      CALL GSBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,IBYTE,0,INVAL,
     S             KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C
      IF (KB1PAR(4).LT.128) GOTO 333
C
C
C***
C*    SET ARRAY IBLOCK TO 0 . IT IS USED TO ZERO FILL RESERVED OCTETS.
C***
C
      DO 100 J=1,24
         IBLOCK(J) = 0
  100 CONTINUE
C
C
C
C
C***
C*    LAT/LONGITUDE GRID, GAUSSIAN GRID AND SPHERICAL HARMONICS
C*    ARE THE ONLY DATA REPRESENTATIONS HANDLED.
C***
C
      IF (KB1PAR(4).EQ.128.OR.KB1PAR(4).EQ.192)
     C   THEN
             IF (KB2PAR(1).NE.0.AND.KB2PAR(1).NE.4.AND.KB2PAR(1).NE.50
     C           .AND.KB2PAR(1).NE.80)
     C          THEN
                     WRITE (*,*)'GRID DESCRIPTION BLOCK NOT YET DEFINED'
                     KERR = -3
                     IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
                     RETURN
                ENDIF
C
C
C
C
C***
C*       LENGTH OF GRID DESCRIPTION BLOCK.
C***
C
C        LENGTH IS 32 OCTETS FOR LAT/LONG, GAUSSIAN AND SPHERICAL
C        HARMONICS . FOR ANY DATA  ON HYBRID LEVELS THE
C        VERTICAL COORDINATES ARE ADDED.
C
         IF(KB2PAR(1).EQ.80) THEN
            I=52
         ELSE
            I = 32
         ENDIF
         IF (KB1PAR(6).GT.108) I = I + KLENV * 4
C
         IL8=I
         CALL SBYTE_MF (KGRIB(KWORD),IL8,IOFF,24)
         CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
         IF (KERR.NE.0)  THEN 
           IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
           RETURN
         ENDIF
C
C
C
C
C***
C*       NUMBER OF UNUSED BITS AT END OF BLOCK.
C*       CURRENT DEFINITION IMPLIES 0.
C***
C
         ILBLOCK=IBLOCK
         CALL SBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,8)
         CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
         IF (KERR.NE.0)  THEN 
           IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
           RETURN
         ENDIF
C
C
C
C
C
C***
C*       NEXT OCTET IS RESERVED.
C***
C
         ILBLOCK=IBLOCK
         CALL SBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,8)
         CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
         IF (KERR.NE.0)  THEN 
           IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
           RETURN
         ENDIF
C
C
C
C
C
C
C***
C*       DATA REPRESENTATION TYPE.
C***
C
         ILB2PAR=KB2PAR
         CALL SBYTE_MF (KGRIB(KWORD),ILB2PAR(1),IOFF,8)
         CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
         IF (KERR.NE.0)  THEN 
           IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C*       LAT/LONG OR GAUSSIAN GRID.
C***
C
         IF (KB2PAR(1).EQ.0.OR.KB2PAR(1).EQ.4)
     C    THEN
C
C            NUMBER OF LAT/LONG POINTS.
C
             ILB2PAR=KB2PAR
             CALL GSBYTE_MF (KGRIB(KWORD),ILB2PAR(2),IOFF,16,0,2,
     S                    KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C            LAT/LONG OF ORIGIN.
C            SIGN BIT SET TO 1 IF VALUES ARE NEGATIVE.
C
             DO 200 J=1,2
C
             IF (KB2PAR(J+3).GE.0) THEN
                ILAT(J) = KB2PAR(J+3)
             ELSE
                ILAT(J) = 2**23 - KB2PAR(J+3)
             ENDIF
C
  200        CONTINUE
C
             ILBLOCK=IBLOCK
             CALL GSBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,24,0,2,
     S                    KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C            RESOLUTION FLAG.
C
C            INCREMENTS       BINARY CODE     DECIMAL VALUE
C
C             NOT GIVEN        00000000              0
C              GIVEN           10000000            128
C
C
C            SHIFT 1 BIT FLAG FIELD TO LOW ORDER BIT
C
             ITEMP = KB2PAR(6) / 128
C
             IF (ITEMP.LE.0.OR.ITEMP.GT.1)
     C          THEN
                    CALL PRTBIN_MF (KB2PAR(6),8,ITEMP,IERR)
                    WRITE (*,9011) ITEMP
 9011               FORMAT (TR1,'INVALID RESOLUTION FLAG ',I8.8)
                    KERR = 11
                    IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
                    RETURN
                ENDIF
C
C            SET TO ALL 1-BITS UNUSED INCREMENT FIELDS.
C
             IF (KB2PAR(6).EQ.0)
     C          THEN
C                   1111111111111111 BINARY = 65535 DECIMAL
                    KB2PAR(9)  = 65535
                    KB2PAR(10) = 65535
                ENDIF
C
             ILB2PAR=KB2PAR
             CALL SBYTE_MF (KGRIB(KWORD),ILB2PAR(6),IOFF,8)
             CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C            LAT/LONG OF EXTREME POINT.
C
C            SIGN BIT SET TO 1 IF VALUES ARE NEGATIVE.
C
             DO 300 J=1,2
C
             IF (KB2PAR(J+6).GE.0) THEN
                ILAT(J) = KB2PAR(J+6)
             ELSE
                ILAT(J) = 2**23 - KB2PAR(J+6)
             ENDIF
C
  300        CONTINUE
C
             ILLAT=ILAT
             CALL GSBYTE_MF (KGRIB(KWORD),ILLAT(1),IOFF,24,0,2,
     S                    KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C
C
C            DIRECTION INCREMENTS / NUMBER OF LATITUDE LINES.
C
             ILB2PAR=KB2PAR
             CALL GSBYTE_MF (KGRIB(KWORD),ILB2PAR(9),IOFF,16,0,2,
     S                    KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C            SCANNING MODE FLAGS.
C
C            VALID VALUES       VALUE / 32
C             BINARY             DECIMAL
C
C              00000000               0
C              10000000               4
C              01000000               2
C              11000000               6
C              00100000               1
C              10100000               5
C              01100000               3
C              11100000               7
C
C            SHIFT 3 BIT FLAG FIELD TO LOW ORDER BITS.
C
             ITEMP = KB2PAR(11) / 32
C
             IF (ITEMP.LT.0.OR.ITEMP.GT.7)
     C          THEN
                    CALL PRTBIN_MF (KB2PAR(11),8,ITEMP,IERR)
                    WRITE (*,9014) ITEMP
 9014               FORMAT (TR1,'INVALID SCANNING MODE FLAGS ',I8.8)
                    KERR = 12
                    IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
                    RETURN
                ENDIF
C
             ILB2PAR=KB2PAR
             CALL SBYTE_MF (KGRIB(KWORD),ILB2PAR(11),IOFF,8)
             CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C            4 RESERVED OCTETS.
C
             ILBLOCK=IBLOCK
             CALL GSBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,8,0,4,
     S                    KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C***
C*       SPHERICAL HARMONIC DATA.
C***
C
         IF (KB2PAR(1).EQ.50.OR.KB2PAR(1).EQ.80)
     C    THEN
C
C            PENTAGONAL RESOLUTION PARAMETERS.
C
             ILB2PAR=KB2PAR
             CALL GSBYTE_MF (KGRIB(KWORD),ILB2PAR(2),IOFF,16,0,3,
     S                    KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C            REPRESENTATION TYPE AND MODE.
C
             ILB2PAR=KB2PAR
             CALL GSBYTE_MF (KGRIB(KWORD),ILB2PAR(5),IOFF,8,0,2,
     S                    KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C            18 RESERVED OCTETS.
C
             ILBLOCK=IBLOCK
             CALL GSBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,8,0,18,
     S                    KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
          ENDIF
C
C
C***
C*    ADD PARAMETERS NEEDED TO DESCRIBE STRETCHED AND ROTATED
C*    SPHERICAL HARMONICS DATA
C***
C
      IF(KB2PAR(1).EQ.80) THEN
C
C*    POLE OF ROTATION
C
             DO 310 J=1,2
C
             IF (KB2PAR(J+11).GE.0) THEN
                ILAT(J) = KB2PAR(J+11)
             ELSE
                ILAT(J) = 2**23 - KB2PAR(J+11)
             ENDIF
C
  310        CONTINUE
C
             ILLAT=ILAT
             CALL GSBYTE_MF (KGRIB(KWORD),ILLAT(1),IOFF,24,0,2,
     S                    KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C*    ANGLE OF ROTATION
             ILB2PAR=KB2PAR
             CALL CONFP_MF (ILB2PAR(14),IEXP,IMANT)
             ILEXP=IEXP
             CALL SBYTE_MF (KGRIB(KWORD),ILEXP,IOFF,8)
             CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
             ILMANT=IMANT
             CALL SBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,24)
             CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C*    POLE OF STRETCHING
C
             DO 320 J=1,2
C
             IF (KB2PAR(J+14).GE.0) THEN
                ILAT(J) = KB2PAR(J+14)
             ELSE
                ILAT(J) = 2**23 - KB2PAR(J+14)
             ENDIF
C
  320        CONTINUE
             ILLAT=ILAT
             CALL GSBYTE_MF (KGRIB(KWORD),ILLAT(1),IOFF,24,0,2,
     S                    KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
C
C*    STRETCHING  FACTOR
             ILB2PAR=KB2PAR
             CALL CONFP_MF (ILB2PAR(17),IEXP,IMANT)
             ILEXP=IEXP
             CALL SBYTE_MF (KGRIB(KWORD),IEXP,IOFF,8)
             CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
             ILMANT=IMANT
             CALL SBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,24)
             CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
         ENDIF
C
C
C
C***
C*       ADD VERTICAL COORDINATE PARAMETERS FOR HYBRID LEVELS.
C***
C
         IF (KB1PAR(6).GT.108) THEN
C
             DO 400 J=1,KLENV
                 CALL CONFP_MF (PVERT(J),IEXP,IMANT)
                 ILEXP=IEXP
                 CALL SBYTE_MF (KGRIB(KWORD),ILEXP,IOFF,8)
                 CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
                 IF (KERR.NE.0)  THEN 
                   IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
                   RETURN
                 ENDIF
                 ILMANT=IMANT
                 CALL SBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,24)
                 CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
                 IF (KERR.NE.0)  THEN 
                   IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
                   RETURN
                 ENDIF
  400        CONTINUE
C
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
  333 CONTINUE
C
C
C********************************************************************
C*
C*    BLOCK 3 (IF REQUIRED) - BIT MAP BLOCK.
C*
C********************************************************************
C
      IF (KB1PAR(4).EQ.64.OR.KB1PAR(4).EQ.192)
     C     THEN
               WRITE (*,*)'BIT MAP BLOCK NOT YET DEFINED'
               KERR = -3
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C*    BLOCK 4 - BINARY DATA BLOCK.
C*
C********************************************************************
C
C
C***
C*    RETAIN POINTERS TO WORD AND BIT POSITION OF BINARY DATA
C*    BLOCK LENGTH FIELD. ENTER LENGTH AS 0.
C***
C
      IPW = KWORD
      IPB = IOFF
C
      ILBLOCK=IBLOCK
      CALL SBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,24)
      CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C
C
C
C***
C*    4 BIT FLAG / 4 BIT COUNT OF UNUSED BITS AT END OF BLOCK OCTET.
C***
C
C     FLAG IS 1000 FOR SPHERICAL HARMONICS, 0000 FOR LAT/LONG
C     OR GAUSSIAN GRID, 1100 FOR SPHERICAL HARM. COMPLEX PACKING
C
      IREP = 0
      IFLAG=0
      ICPACK=KCPACK
C
      IF (KB2PAR(1).EQ.50.OR.KB2PAR(1).EQ.80) THEN
         IFLAG=128
         IREP=1
C
         IF (ICPACK.LT.0) THEN
           WRITE (*,*) 'CODEGA : COMPLEX PACKING CODE ERROR'
           KERR = -5
           IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
           RETURN
         ELSEIF (IABS (KSCALP).GE.2**15) THEN
           WRITE (*,*) 'CODEGA : LAPLACIAN SCALING FACTOR ERROR'
           KERR = -6
           IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
           RETURN
         ENDIF
C
         IF(ICPACK.NE.0.AND.KB2PAR(6).EQ.2) THEN
            IFLAG=IFLAG+64
            IREP=(ICPACK+1)**2
         ELSE
            ICPACK=0
         ENDIF
C
      ENDIF
C
      ILFLAG=IFLAG
      CALL SBYTE_MF (KGRIB(KWORD),ILFLAG,IOFF,8)
      CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C
C
C
C***
C*    FIND MAXIMUM AND MINIMUM VALUES IN DATA ARRAY. FOR
C*    DATA IN SPHERICAL HARMONIC FORM THE FIRST WORD CONTAINS
C*    THE REAL (0,0)COEFFICIENT, WHICH IS TREATED SEPARATELY.
C*    FOR COMPLEX PACKING AND SPECTRAL DATA REPRESENTATION MODE=2
C*    THE FIRST (ICPACK+1)**2 COEFFICENTS ARE NOT PACKED
C***
C
C
      ILEN = ILENF - IREP
      CALL MXMN_MF (PFDATA(IREP+1),ILEN,PMAX,PMIN)
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
C***
C*    COMPUTE REFERENCE VALUE AND SCALE FACTOR.
C***
C
      IF (LDARPE) THEN
C
        ZCOEFF = PMAX - PMIN
        ZEPSIL=1.E-290
C
        IF ( ZCOEFF .LE. ZEPSIL ) THEN
          ZAUXIL=MIN ( ABS (PMIN), ABS (PMAX) )
          IF ( ZAUXIL .LE. ZEPSIL ) ZAUXIL=0.
C
          PMAX=SIGN (ZAUXIL,PMAX)
          PMIN=PMAX
          ZSCALE=0.
        ELSE
          ZSCALE=REAL (2**KBITS-1) / ZCOEFF
        ENDIF
C
C             Scale factor and reference value forced to zero.
C
        ISCALE=0
        ZREFER=PMIN
        IEXP=0
        IMANT=0
C
      ELSEIF (IMISS.EQ.1) THEN
C
C           Scaling factor,
C           EXPONENT AND MANTISSA SET TO ALL 1-BITS FOR MISSING DATA.
C
        ISCALE = 65535
        IEXP   = 255
        IMANT  = 16777215
        ZREFER = PMIN
C
      ELSE
C
C     CONVERT MINIMUM VALUE (PMIN) TO GRIB FORMAT (IEXP,IMANT).
C
C     Its image decoded back is being used as reference value to compute
C     (in PACKGB) normalized integer field values that will be
C     chained into the binary string.
C
        CALL CONFI (PMIN,IEXP,IMANT,ZREFER)
C
        ZS = (PMAX-ZREFER)/(2**(KBITS+1)-1)
        ZAUXIL=1.
        ZAUXI2=2.
C
C     CONFI may return ZREFER > PMAX if the range
C     of PMAX-PMIN is smaller than the GRIB accuracy
C
        IF (ZS.GT.0.) THEN
           ZS = LOG(ZS)/LOG(ZAUXI2) + ZAUXI2
        ELSE
           ZS=0.
        ENDIF
        ISCALE = MIN(INT(ZS),INT(ZS+SIGN(ZAUXIL,ZS)))
C
C      Absolute value of ISCALE is limited, to avoid problems due to
C      fields constant except on a few points.
C
        ISCALE = MAX(-99,MIN(99,ISCALE))
        ZSCALE = ZAUXI2** (-ISCALE)
C
      ENDIF
C
C
C
C
C     SET SIGN BIT (BIT 16) AND PUT SCALE FACTOR INTO REMAINING
C     15 BITS OF CODED 16 BIT FIELD.
C
      IF (ISCALE.GE.0) THEN
        ISCALX = ISCALE
      ELSE
        ISCALX= 2**15 - ISCALE
      ENDIF
C
      ILSCALX=ISCALX
      CALL SBYTE_MF (KGRIB(KWORD),ILSCALX,IOFF,16)
      CALL OFFSET_MF (IOFF,1,KWORD,16,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C
C
C
C
      IF (IEXP.EQ.0.AND.IMANT.EQ.0) THEN
        ILMANT=IMANT
        CALL SBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,32)
        CALL OFFSET_MF (IOFF,1,KWORD,32,KNBIT,KLENG,KERR)
      ELSE
        ILEXP=IEXP
        CALL SBYTE_MF (KGRIB(KWORD),ILEXP,IOFF,8)
        CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
        IF (KERR.NE.0)  THEN 
          IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
          RETURN
        ENDIF
        ILMANT=IMANT
        CALL SBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,24)
        CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
      ENDIF
C
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C
C
C
C***
C*    NUMBER OF BITS IN EACH DATA VALUE.
C***
C
      ILBITS=KBITS
      CALL SBYTE_MF (KGRIB(KWORD),ILBITS,IOFF,8)
      CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C
C
C
C***
C*    IF COMPLEX PACKING DESCRIBE PACKING
C***
C
      IF (IREP.NE.0.AND.ICPACK.NE.0) THEN
C
C     START OF PACKED DATA
C
      ISTPA=19+IREP*4
C
C Removing of the padding (06 sept 2000, not correct and useless)
C
C      IRESTE=MOD (IOFF+40+IREP*32,KNBIT)
CC
C      IF (IRESTE.NE.0) THEN
CC
C        IF ( MOD (KNBIT,KBITS).EQ.0 .AND. MOD (IRESTE,8).EQ.0 .AND.
C     S       MOD (KBITS,IRESTE).NE.0 ) THEN
CC
CC        Packed data will start on a word boundary,
CC     and each packed data will be into a single word.
CC
C          ISTPA=ISTPA+IRESTE/8
C        ENDIF
CC
C      ENDIF
C
      ILSTPA=ISTPA
      CALL SBYTE_MF (KGRIB(KWORD),ILSTPA,IOFF,16)
      CALL OFFSET_MF (IOFF,1,KWORD,16,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C     SCALING FACTOR.
C
      IF (KSCALP.GE.0) THEN
        ISCALP=KSCALP
      ELSE
        ISCALP=2**15-KSCALP
      ENDIF
C
      ILSCALP=ISCALP
      CALL SBYTE_MF (KGRIB(KWORD),ILSCALP,IOFF,16)
      CALL OFFSET_MF (IOFF,1,KWORD,16,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C     TRUNCATION FOR UNPACKED PART OF SPECTRAL DATA, CURRENTLY
C     ONLY TRIANGULAR TRUNCATION SUPPORTED
C
      IAUXIL=ICPACK*(1+2**8*(1+2**8))
      ILAUXIL=IAUXIL
      CALL SBYTE_MF (KGRIB(KWORD),ILAUXIL,IOFF,24)
      CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      ENDIF
C***
C*    IF SPHERICAL HARMONICS DATA, NEXT 4 OCTETS ARE DIFFERENT FORMAT.
C***
C
      IF (IREP.NE.0)
     C    THEN
C
C            STORE IREP COEFFICIENTS IN FLOATING
C            POINT FORM.
C
        IF (LDARPE) THEN
C
C                THE "UNPACKED" PART IS FILLED WITH ZEROES.
C
          IEXP=0
          IMANT=0
C
          DO 510 J=1,IREP
          ILMANT=IMANT
          CALL SBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,32)
          CALL OFFSET_MF (IOFF,1,KWORD,32,KNBIT,KLENG,KERR)
          IF (KERR.NE.0)  THEN 
            IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
  510     CONTINUE
C
        ELSE
C
C                STANDARD CASE.
C
             DO 520 J=1,IREP
             CALL CONFP_MF (PFDATA(J),IEXP,IMANT)
C
             ILEXP=IEXP
             CALL SBYTE_MF (KGRIB(KWORD),ILEXP,IOFF,8)
             CALL OFFSET_MF (IOFF,1,KWORD,8,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
             ILMANT=IMANT
             CALL SBYTE_MF (KGRIB(KWORD),ILMANT,IOFF,24)
             CALL OFFSET_MF (IOFF,1,KWORD,24,KNBIT,KLENG,KERR)
             IF (KERR.NE.0)  THEN 
               IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
               RETURN
             ENDIF
 520         CONTINUE
C
          ENDIF
C
        ENDIF
C
C
C
C
C***
C*    SCALE AND STORE DATA VALUES.
C***
C
C     DO 600 J = IREP+1 , ILENF
C        IPDATA    = NINT ( (PFDATA(J)-ZREFER) * ZSCALE )
C        PFDATA(J)=OR (IPDATA,0)
C 600 CONTINUE
C
      ILENFM = ILENF - IREP
      CALL PACKGB (PFDATA(IREP+1),PFDATA(IREP+1),ZREFER,ZSCALE,ILENFM)
      CALL GSBYTE_MF (KGRIB(KWORD),PFDATA(IREP+1),IOFF,KBITS,0,ILENFM,
     S             KNBIT,'C',KLENG,KERR,KWORD,.TRUE.)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
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
C
C
C
C***
C*    ENTER LENGTH OF BINARY DATA BLOCK, HAVING ENSURED THAT
C*    THE LENGTH IS AN EVEN NUMBER OF OCTETS.
C***
C
C     LENGTH OF BINARY DATA BLOCK IN BITS.
C
      ILBIN = (KWORD-IPW) * KNBIT + IOFF - IPB
C
      IL = MOD (ILBIN,16)
C
C     FILL UNUSED PORTION OF LAST 2 OCTETS WITH BINARY ZEROES.
C
      ILNIL = 0
C
      IF (IL.NE.0) THEN
              ILNIL = 16 - IL
              ILBLOCK=IBLOCK
              CALL SBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,ILNIL)
              CALL OFFSET_MF (IOFF,1,KWORD,ILNIL,KNBIT,KLENG,KERR)
              IF (KERR.NE.0)  THEN 
                IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
                RETURN
              ENDIF
          ENDIF
C
      ILBIN = (KWORD-IPW) * KNBIT + IOFF - IPB
C
C     ENTER LENGTH - IN OCTETS - OF DATA BLOCK.
C
      ILBIN = ILBIN / 8
      ILLBIN=ILBIN
      CALL SBYTE_MF (KGRIB(IPW),ILLBIN,IPB,24)
      CALL OFFSET_MF (IPB,1,IPW,24,KNBIT,KLENG,KERR)
C
C     ENTER NUMBER OF UNUSED BITS IN FLAG/BIT COUNT FIELD.
C
      IFLAG = IFLAG + ILNIL
      ILFLAG=IFLAG
      CALL SBYTE_MF (KGRIB(IPW),ILFLAG,IPB,8)
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
C
C
C
C********************************************************************
C*
C*    BLOCK 5 - END BLOCK.
C*
C********************************************************************
C
C
C***
C*    ADD 7 7 7 7 TO CODED DATA.
C***
C
      IBLOCK(1) = 55
      IBLOCK(2) = 55
      IBLOCK(3) = 55
      IBLOCK(4) = 55
C
      ILBLOCK=IBLOCK
      CALL GSBYTE_MF (KGRIB(KWORD),ILBLOCK(1),IOFF,8,0,4,
     S             KNBIT,'C',KLENG,KERR,KWORD,.FALSE.)
      IF (KERR.NE.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C
C***
C*    SET ANY UNUSED PART OF LAST WORD TO BINARY ZEROES.
C***
C
      IF (IOFF.NE.KNBIT)
     C   THEN
             IBITS = KNBIT - IOFF
             ILBLOCK=IBLOCK
             CALL SBYTE_MF (KGRIB(KWORD),ILBLOCK(5),IOFF,IBITS)
         ENDIF
C
C
C
C***
C*    ROUND TO 120 OCTETS, IF REQUIRED.
C***
C
      IF (KROUND.EQ.1)
     C   THEN
             INUMBI = KWORD * KNBIT
             I = INUMBI / 960
             I = I * 960
             I = INUMBI - I
             IF (I.NE.0) I = (960 - I) / KNBIT
C
             DO 700 J=KWORD+1,KWORD+I
                KGRIB(J) = 0
  700        CONTINUE
C
             KWORD = KWORD + I
         ENDIF
C
      IF (LHOOK) CALL DR_HOOK('CODEGA',1,ZHOOK_HANDLE)
      END

