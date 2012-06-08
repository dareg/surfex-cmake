C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIENG_MT (LFI, KNUMER, KNIMES, KCODE, LDFATA, 
     S                      CDMESS, CDNSPR,CDACTI )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        THIS SUBROUTINE PRINTS STANDARD MESSAGES FROM LFI INDEXED-
C     FILE SOFTWARE, ABORTING PROGRAM IF REQUIRED.
C        Messages related to "debugging mode" are directly printed
C     by concerned subroutines.
C
C        This subroutine is the english version, translated from the
C     the original (french) one, and is always called through the
C     "hat" routine LFIEMS.
C     For french version see subroutine LFIEFR.
C     ( Pour la version francaise, voir LFIEFR )
C**
C  DUMMY ARGUMENTS : KNUMER ==> Logical Unit concerned, if any;
C  ( all INPUT ones )           ( if LFI%JPNIL ==> no Logical Unit )
C                    KNIMES ==> Level (0,1,2) of Message;
C                    KCODE  ==> Response code of action concerned;
C                    LDFATA ==> True if Abort of program is required;
C                    CDMESS ==> If KNIMES#0, Message to print;
C                    CDNSPR ==> Subroutine name which calls LFIEMS;
C                    CDACTI ==> Name of FORTRAN input/output action
C                               if KCODE >0, else... it depends !
C*
C     !----------------------------------------------------------------!
C     !   TABLE OF POSSIBLE VALUES FOR RESPONSE CODES OF LFI SOFTWARE  !
C     !----------------------------------------------------------------!
C
C-----------------------------------------------------------------------
C      0 ==> No error has been detected, everything is OK.
C-----------------------------------------------------------------------
Cpositive==> It is the (absolute value of) FORTRAN response code
C  value     from an OPEN, READ, WRITE, CLOSE or INQUIRE instruction;
C            see vendor's reference manual for exact meaning.
C-----------------------------------------------------------------------
C     -1 ==> Logical Unit currently not opened for the software.
C-----------------------------------------------------------------------
C     -2 ==> "LEVEL" value outside [0-2] range .
C-----------------------------------------------------------------------
C     -3 ==> Bad lock option (internal subroutine "LFIVER") .
C-----------------------------------------------------------------------
C     -4 ==> Explicit change for Multi-Tasking mode, but almost one unit
C            is currently open-problems may arise (subroutine "LFIINI").
C-----------------------------------------------------------------------
C     -5 ==> Logical Unit is currently opened (LFIOUV, LFIAFM, LFISFM) .
C-----------------------------------------------------------------------
C     -6 ==> Not enough space within tables to open requested Unit.
C            (LFIOUV)
C-----------------------------------------------------------------------
C     -7 ==> Invalid "STATUS" for FORTRAN instruction "OPEN" (LFIOUV) .
C-----------------------------------------------------------------------
C     -8 ==> Incompatible values given for "LDNOMM" and "CDSTTO":
C            a file which "STATUS" is 'OLD' or 'NEW' must have a name .
C            (LFIOUV) (THIS REPONSE CODE HAS CURRENTLY NO MORE SENSE)
C-----------------------------------------------------------------------
C     -9 ==> Incompatibility between "STATUS" 'NEW' or 'OLD' and (respe-
C            ctively) file existence or non-existence (LFIOUV) .
C-----------------------------------------------------------------------
C    -10 ==> The file is not a LFI one, or may not be treated through
C            this configuration or version of the software (LFIOUV) .
C-----------------------------------------------------------------------
C    -11 ==> File not closed after a modification (LFIOUV): this
C            error is not fatal if "LDERFA" is .FALSE., but in such a
C            case file integrity and data coherence are not guaranteed.
C            Note that once a file has got such problem, this response
C            code will stay even after a subsequent modification.
C-----------------------------------------------------------------------
C    -12 ==> File has a "STATUS" 'OLD' but an error occurred when
C            reading the first physical record of file (LFIOUV) .
C-----------------------------------------------------------------------
C    -13 ==> File is already open for another LFI logical unit.
C            (LFIOUV)
C-----------------------------------------------------------------------
C    -14 ==> Incorrect value for INTEGER argument (generally negative) .
C-----------------------------------------------------------------------
C    -15 ==> Incorrect CHARACTER argument (too long, for instance).
C-----------------------------------------------------------------------
C    -16 ==> Incoherence in Tables, File, internal calls, software.
C            THIS ERROR MAY NEVER BE FILTERED. ALWAYS FATAL.
C-----------------------------------------------------------------------
C    -17 ==> Too many logical records to store an extra one (LFIECR) .
C            (note that logical records consist of user-readable data
C             records, but also of holes cataloged in index... which are
C             created when existing records may not be rewritten in
C             place, or when records are suppressed; such holes may be
C             "re-cycled")
C-----------------------------------------------------------------------
C    -18 ==> A logical record name formed only with SPACES is invalid.
C            (for internal use of LFI software, holes in index are
C             described by a blank record name)
C-----------------------------------------------------------------------
C    -19 ==> File opened with "STATUS" set to 'SCRATCH', so may not be
C            kept at CLOSE time: 'KEEP' is illicit for "CDSTTC" (LFIFER)
C            if this error is not fatal, then a FORTRAN "CLOSE" without
C            "STATUS" parameter is performed, in the same manner as if
C            "CDSTTC" is neither 'KEEP' nor 'DELETE'.
C-----------------------------------------------------------------------
C    -20 ==> No logical record with such name found within logical unit.
C            (LFILEC, LFIREN, LFISUP)
C-----------------------------------------------------------------------
C    -21 ==> Requested logical record is LONGER (has more data) in file;
C            if this error is not fatal, then a PARTIAL read is
C            performed, at requested length.
C            (LFILAP, LFILAS, LFILEC)
C-----------------------------------------------------------------------
C    -22 ==> Requested logical record SHORTER (has less data) in file;
C            even if this error is not fatal, NO READING OF DATA OCCURS.
C            (LFILAP, LFILAS, LFILEC) .
C-----------------------------------------------------------------------
C    -23 ==> No or no more "NEXT" record to read (LFILAS) .
C-----------------------------------------------------------------------
C    -24 ==> The character variable given as actual output argument is
C            TOO SHORT to store the record NAME, even when suppressing
C            any spaces at the end of the name.
C            (LFICAP, LFICAS, LFILAP, LFILAS)
C-----------------------------------------------------------------------
C    -25 ==> The new name of the logical record is (already) used for
C            another logical record within the file (LFIREN).
C-----------------------------------------------------------------------
C    -26 ==> No or no more "PREVIOUS" logical record to read (LFILAP).
C-----------------------------------------------------------------------
C    -27 ==> Insufficient CONTIGUOUS space within tables to treat the
C            "multiple" file requested (LFIOUV) .
C-----------------------------------------------------------------------
C    -28 ==> Multiply factor (of elementary physical record length) too
C            big for the current configuration of the software.
C            (LFIOUV, LFIAFM, LFIFMD)
C-----------------------------------------------------------------------
C    -29 ==> Not enough space within tables to store the multiply factor
C            to be associated to logical unit (LFIAFM) .
C-----------------------------------------------------------------------
C    -30 ==> Logical unit number invalid for FORTRAN.
C-----------------------------------------------------------------------
C    -31 ==> Logical unit has no multiply factor predefined.
C            (LFISFM)
C-----------------------------------------------------------------------
C
#ifndef f77
      USE SDL_MOD   , ONLY : SDL_SRLABORT
#endif
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KNUMER, KNIMES, KCODE, ILDMES, ILBLAN, INLNOM, INUMER
      INTEGER IDECBL, IPOSBL, ILACTI, ILACT2, ILNSPR, ILMESU, IJL, J, IJ
      INTEGER INBALO, ILMESA, INLIGN, IDECAL
C
      LOGICAL LDFATA
C
      CHARACTER  CDNSPR*(*), CLJOLI*6, CDMESS*(*), CLMESA*80, CDACTI*(*)
C
#include "lficom2.h"
#include "lficom_mt.h"
C**
C     1.  -  INITIALISATIONS.
C-----------------------------------------------------------------------
C
C        Search for "useful" length of argument CDACTI.
C        (i.e. not taking into account any blank characters at the end)
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIENG_MT',0,ZHOOK_HANDLE)
      IDECBL=0
C
  101 CONTINUE
      IPOSBL=IDECBL+INDEX (CDACTI(IDECBL+1:),' ')
C
      IF (IPOSBL.LE.IDECBL) THEN
        ILACTI=LEN (CDACTI)
      ELSEIF (CDACTI(IPOSBL:).EQ.' ') THEN
        ILACTI=IPOSBL-1
      ELSE
        IDECBL=IPOSBL
        GOTO 101
      ENDIF
C
      ILACT2=MIN0 (ILACTI,LFI%JPNCPN)
      ILACTI=MIN0 (ILACT2,8)
      ILNSPR=MIN0 (LEN (CDNSPR),LFI%JPLSPX)
C
C        Prefix (and possible suffix) for the message(s).
C
      IF (LDFATA) THEN
        CLJOLI=' *****'
      ELSEIF (KNIMES.EQ.0.OR.KCODE.NE.0) THEN
        CLJOLI=' */*/*'
      ELSE
        CLJOLI=' /////'
      ENDIF
C
      IF (KNIMES.NE.0) THEN
C**
C     2.  -  PRINTS MESSAGE PREPARED BY SUBROUTINE WHICH CALLED LFIEMS.
C-----------------------------------------------------------------------
C
        ILMESU=MIN0 (LEN (CLMESS)-
     S         LEN (CLJOLI)-ILNSPR-4,LEN (CDMESS))
        CLMESS=CLJOLI//' '//CDNSPR(1:ILNSPR)//' - '
     S             //CDMESS(1:ILMESU)
        WRITE (UNIT=LFI%NULOUT,FMT='(A)') CLMESS
      ENDIF
C
      IF (KNIMES.EQ.0.OR.LDFATA) THEN
C**
C     3.  -  CONSTITUTION OF "AD HOC" MESSAGE, DEPENDING OF *KCODE*.
C-----------------------------------------------------------------------
C
C     Before, check if logical unit considered corresponds to a logical
C     unit currently opened for LFI software (or not).
C
        IF (KNUMER.EQ.LFI%JPNIL) THEN
          IJL=0
        ELSE
C
          DO 301 J=1,LFI%NBFIOU
          IJL=LFI%NUMIND(J)
          IF (KNUMER.EQ.LFI%NUMERO(IJL)) GOTO 302
  301     CONTINUE
C
          IJL=0
        ENDIF
C
  302   CONTINUE
C
        IF (KCODE.GT.0) THEN
C
          IF ((CDACTI.EQ.'READ'.OR.CDACTI.EQ.'WRITE')
     S        .AND.LFI%NUMAPH(IJL).GT.0) THEN
            WRITE (UNIT=CLMESS,FMT='(''ERROR "'',A,''"'',I7,
     S             '',UNIT'',I3,'',REC.NUM'',I6,'',*'',I6,
     S             '' WORDS'')') CDACTI(1:ILACTI),KCODE,KNUMER,
     S                            LFI%NUMAPH(IJL),
     S                            LFI%JPLARD*LFI%MFACTM(IJL)
          ELSE
            WRITE (UNIT=CLMESS,
     S             FMT='(''FORTRAN "'',A,''" ERROR, CODE=''
     S             ,I7,'', UNIT='',I3)') CDACTI(1:ILACTI),KCODE,KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-1) THEN
          WRITE (UNIT=CLMESS,FMT='(''LOGICAL UNIT'',I3,
     S           '' NOT OPENED FOR LFI SOFTWARE'')') KNUMER
C
        ELSEIF (KCODE.EQ.-2) THEN
C
          IF (KNUMER.EQ.LFI%JPNIL) THEN
            CLMESS='ACTUAL VALUE FOR LEVEL "KNIVAU" '//
     S                 'OUTSIDE [0-2] RANGE'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S   '(''MESSAGE LEVEL OUTSIDE [0-2] RANGE, UNIT'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-3) THEN
          ILDMES=MIN0 (8,LEN (CDMESS))
          CLMESS='UNKNOWN ACTION '''//CDMESS(1:ILDMES)
     S           //''' ON LOCKS'
C
        ELSEIF (KCODE.EQ.-4) THEN
          CLMESS='EXPL.CHANGE OF MULTI-TASKING MODE '//
     S               'WITH UNIT(S) OPENED'
C
        ELSEIF (KCODE.EQ.-5) THEN
          WRITE (UNIT=CLMESS,FMT='(''LOGICAL UNIT'',I3,
     S          '' ALREADY OPENED FOR LFI - AND SHOULD NOT.'')') KNUMER
C
        ELSEIF (KCODE.EQ.-6) THEN
          WRITE (UNIT=CLMESS,FMT='(I3,'' ENTRIES,'',
     S    '' NOT ENOUGH PLACE WITHIN TABLES FOR UNIT'',I3)')
     S    LFI%JPNXFI,KNUMER
C
        ELSEIF (KCODE.EQ.-7) THEN
          WRITE (UNIT=CLMESS,FMT='(''FORTRAN STATUS'''''',A,
     S           '''''' UNKNOWN, UNIT'',I3)') CDACTI(1:ILACTI),KNUMER
C
        ELSEIF (KCODE.EQ.-8) THEN
          WRITE (UNIT=CLMESS,FMT='(''UNIT'',I3,'' OF STATUS ''''''
     S,A,'''''' MUST HAVE AN EXPLICIT NAME'')') KNUMER,CDACTI(1:ILACTI)
C
        ELSEIF (KCODE.EQ.-9) THEN
C
        IF (CDACTI.EQ.'OLD') THEN
          WRITE (UNIT=CLMESS,FMT=
     S'(''STATUS ''''OLD'''' BUT FILE DOES NOT EXIST, UNIT'',
     S      I3)') KNUMER
        ELSE
          ILBLAN=INDEX (CDACTI(1:ILACTI),' ')
          IF (ILBLAN.GT.1) ILACTI=ILBLAN-1
          WRITE (UNIT=CLMESS,FMT=
     S'(''STATUS '''''',A,'''''' BUT FILE ALREADY EXISTS, UNIT'',
     S  I3)') CDACTI(1:ILACTI),KNUMER
        ENDIF
C
        ELSEIF (KCODE.EQ.-10) THEN
          WRITE (UNIT=CLMESS,FMT='(''INCOMPATIBILITY'',
     S           '' FILE / SOFTWARE, UNIT'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-11) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''UNIT'',I3,'' NOT CLOSED AFTER '',
     S           ''ITS LAST MODIFICATION'')') KNUMER
C
        ELSEIF (KCODE.EQ.-12) THEN
          WRITE (UNIT=CLMESS,FMT='(''UNIT'',I3,
     S  '' OF STATUS ''''OLD'''' - READ OF FIRST RECORD FAILED'')')
     S       KNUMER
C
        ELSEIF (KCODE.EQ.-13) THEN
          INLNOM=1
          INUMER=LFI%JPNIL
C
          DO 131 J=1,LFI%NBFIOU
          IJ=LFI%NUMIND(J)
C
          IF (CDACTI.EQ.LFI%CNOMFI(IJ)) THEN
            INUMER=LFI%NUMERO(IJ)
            INLNOM=MIN0 (LFI%NLNOMF(IJ),LEN (CLMESS)-3)
            GOTO 132
          ENDIF
C
  131     CONTINUE
C
  132     CONTINUE
          CLMESS=' '''//CDACTI(1:INLNOM)//''''
          WRITE (UNIT=LFI%NULOUT,FMT='(A)') CLMESS
          WRITE (UNIT=CLMESS,FMT='(''UNIT'',I3,'' - FILE '',
     S           ''ALREADY OPEN WITH UNIT'',I3)') KNUMER,INUMER
C
        ELSEIF (KCODE.EQ.-14) THEN
C
          IF (CDNSPR.EQ.'LFIECR'.OR.CDNSPR.EQ.'LFILEC'.OR.
     S        CDNSPR.EQ.'LFILAS'.OR.CDNSPR.EQ.'LFILAP') THEN
            WRITE (UNIT=CLMESS,FMT=
     S   '(''INCORRECT RECORD LENGTH, UNIT'',I3)') KNUMER
          ELSEIF (KNUMER.EQ.LFI%JPNIL) THEN
            CLMESS='INCORRECT ENTRY IN *LFI%NUMERO* TABLE'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S   '(''INCORRECT INTEGER TYPE ARGUMENT, UNIT'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-15) THEN
          WRITE (UNIT=CLMESS,FMT='(''RECORD NAME INCORRECT OR '',
     S           ''TOO LONG, UNIT'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-16) THEN
          WRITE (UNIT=CLMESS,FMT='(''INCOHERENCE (TABLES, FILE, '',
     S           ''INTERNAL CALLS, SOFTWARE), UNIT'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-17) THEN
C
          IF (IJL.NE.0) THEN
            INBALO=LFI%MDES1D(IXM(LFI%JPNALO,IJL))
          ELSE
            INBALO=LFI%JPNIL
          ENDIF
C
          WRITE (UNIT=CLMESS,
     S           FMT='(I6,'' RECORDS, INDEX FULL, UNIT'',
     S           I3)') INBALO,KNUMER
C
        ELSEIF (KCODE.EQ.-18) THEN
          WRITE (UNIT=CLMESS,FMT='(''BLANK RECORD NAME IS INVALID'',
     S           '', UNIT'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-19) THEN
          WRITE (UNIT=CLMESS,FMT='(''UNIT'',I3,
     S           '' IS ''''SCRATCH'''', SO MAY NOT BE KEPT'')') KNUMER
C
        ELSEIF (KCODE.EQ.-20) THEN
          WRITE (UNIT=CLMESS,FMT='(''RECORD "'',A,
     S           ''" NOT FOUND, UNIT'',I3)') CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-21) THEN
          WRITE (UNIT=CLMESS,FMT='(''RECORD "'',A,
     S    ''" *LONGER* THAN REQUESTED, UNIT'',I3)')
     S      CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-22) THEN
          WRITE (UNIT=CLMESS,FMT='(''RECORD "'',A,
     S    ''" *SHORTER* THAN REQUESTED-UNIT'',I3)')
     S      CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-23) THEN
          WRITE (UNIT=CLMESS,FMT='(''NO/NO MORE NEXT RECORD'',
     S    '' TO READ, UNIT'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-24) THEN
          WRITE (UNIT=CLMESS,FMT='(''CHARAC. VARIABLE TOO SHORT '',
     S    ''FOR "'',A,''", UNIT'',I3)')
     S      CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-25) THEN
          WRITE (UNIT=CLMESS,FMT='(''NEW RECORD NAME: "'',A,
     S    ''" ALREADY USED, UNIT'',I3)')
     S      CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-26) THEN
          WRITE (UNIT=CLMESS,FMT='(''NO/NO MORE PREVIOUS RECORD '',
     S    '' TO READ, UNIT'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-27) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''INSUFFICIENT CONTIGUOUS SPACE WI'',
     S    ''THIN TABLES, UNIT'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-28) THEN
C
          IF (KNUMER.EQ.LFI%JPNIL) THEN
            WRITE (UNIT=CLMESS,
     S             FMT='(''NEW DEFAULT MULTIPLY FACTOR EX'',
     S      ''CEEDS MAXIMUM ('',I3,'')'')') LFI%JPFACX
          ELSE
            WRITE (UNIT=CLMESS,FMT='(''SPECIFIED MULTIPLY FACTOR '',
     S      ''EXCEEDS MAXIMUM ('',I3,''), UNIT'',I3)') LFI%JPFACX,KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-29) THEN
          WRITE (UNIT=CLMESS,FMT='(I3,'' ENTRIES,'',
     S    '' NO MORE PLACE FOR MULTIPLY FACTOR, UNIT'',I3)')
     S    LFI%JPXUFM,KNUMER
C
        ELSEIF (KCODE.EQ.-30) THEN
          WRITE (UNIT=CLMESS,FMT='(''INVALID FORTRAN LOGICAL UNIT'',
     S           '' NUMBER:'',I8)') KNUMER
C
        ELSEIF (KCODE.EQ.-31) THEN
          WRITE (UNIT=CLMESS,FMT='(''LOGICAL UNIT NUMBER'',I3,
     S       '' HAS NO PREDEFINED MULTIPLY FACTOR'')') KNUMER
C
C                  For unexpected error codes...
C
        ELSEIF (KNUMER.EQ.LFI%JPNIL) THEN
          WRITE (UNIT=CLMESS,FMT='(''*UNKNOWN* GLOBAL ERROR CODE'',
     S                             I6)') KCODE
        ELSE
          WRITE (UNIT=CLMESS,FMT='(''*UNKNOWN* ERROR CODE'',I6,
     S           '' ON LOGICAL UNIT'',I3)') KCODE,KNUMER
        ENDIF
C
        ILMESA=LEN (CLMESA)
        ILMESU=ILMESA-1-2*LEN (CLJOLI)-ILNSPR-4
        CLMESA=CLJOLI//' '//CDNSPR(1:ILNSPR)//' - '//
     S         CLMESS(1:ILMESU)//CLJOLI
        WRITE (UNIT=LFI%NULOUT,FMT='(A)') CLMESA
C
C           If logical unit corresponds to a LFI logical unit
C     already opened, its name is printed.
C
        IF (IJL.NE.0) THEN
C
          IF (LFI%NLNOMF(IJL).LE.LFI%JPLFTX) THEN
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') CLJOLI
     S             //' NAME - APPEARENT BUT'
     S             //' COMPLETE - OF LFI LOGICAL UNIT CONCERNED:'
          ELSE
            WRITE (UNIT=CLMESS,FMT='(A,
     S             '' NAME - APPEARENT, AND TRUNCATED BY'',I4,
     S       '' CARACTERES - OF LFI LOGICAL UNIT CONCERNED:'')')
     S      CLJOLI,LFI%NLNOMF(IJL)-LFI%JPLFTX
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') CLMESS
          ENDIF
C
          INLIGN=(LFI%NLNOMF(IJL)-1)/LFI%JPLFIX
          IDECAL=0
C
          DO 801 J=1,INLIGN
          WRITE (UNIT=LFI%NULOUT,FMT='(A)')
     S           LFI%CNOMFI(IJL)(IDECAL+1:IDECAL+LFI%JPLFIX)//'...'
          IDECAL=IDECAL+LFI%JPLFIX
  801     CONTINUE
C
          IF (LFI%NLNOMF(IJL).LE.LFI%JPLFTX) THEN
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') 
     S             LFI%CNOMFI(IJL)(IDECAL+1:LFI%NLNOMF(IJL))
          ELSE
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') 
     S             LFI%CNOMFI(IJL)(IDECAL+1:LFI%JPLFTX)
     S             //'...'
          ENDIF
C
          IF (LFI%CNOMSY(IJL).NE.LFI%CNOMFI(IJL)) THEN
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') CLJOLI//
     S ' *SYSTEM* NAME (APPEARENT) OF LFI LOGICAL UNIT CONCERNED:'
            INLIGN=(LFI%NLNOMS(IJL)-1)/LFI%JPLFIX
            IDECAL=0
C
            DO 802 J=1,INLIGN
            WRITE (UNIT=LFI%NULOUT,FMT='(A)')
     S             LFI%CNOMSY(IJL)(IDECAL+1:IDECAL+LFI%JPLFIX)//'...'
            IDECAL=IDECAL+LFI%JPLFIX
  802       CONTINUE
C
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') 
     S             LFI%CNOMSY(IJL)(IDECAL+1:LFI%NLNOMS(IJL))
          ENDIF
C
        ENDIF
C
        WRITE (UNIT=LFI%NULOUT,FMT='(A)') CLMESA
        IF (LDFATA.AND.KCODE.NE.0) THEN
C
C            Aborts program.
C
          CALL SDL_SRLABORT
        ENDIF
C
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIENG_MT',1,ZHOOK_HANDLE)
      END

