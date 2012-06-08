C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFICFG_MT (LFI)
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C         Sous-programme imprimant les "PARAMETER" de base definissant
C     la configuration dans laquelle tourne le logiciel de fichiers
C     indexes LFI.
C**
C         Ce sous-programme n'a pas d'arguments.
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER IREP, INUMER, INIMES
#include "lficom_mt.h"
C
C**
C     1.  -  INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFICFG_MT',0,ZHOOK_HANDLE)
      IF (LFI%LFICFG_LLPREA) THEN
        CALL LFIINI_MT (LFI, 2)
        LFI%LFICFG_LLPREA=.FALSE.
      ENDIF
C
C         Envoi d'une banniere.
C
      IREP=0
      INUMER=LFI%JPNIL
      CLNSPR='LFICFG'
      WRITE (UNIT=LFI%NULOUT,FMT='(///)')
C
      IF (LFI%LFRANC) THEN
        CLMESS='Configuration du logiciel LFI'//
     S             ' ("PARAMETER" de base):'
      ELSE
        CLMESS='Configuration of LFI software (basic "PARAMETER"):'
      ENDIF
C
      INIMES=2
      LLFATA=.FALSE.
      CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CLACTI)
      WRITE (UNIT=LFI%NULOUT,FMT='(/)')
C**
C     2.  -  IMPRESSION DES "PARAMETER" DE BASE DU LOGICIEL.
C-----------------------------------------------------------------------
C
      IF (LFI%LFRANC) THEN
        WRITE (UNIT=LFI%NULOUT,FMT=9005)LFI%JPNBIM
        WRITE (UNIT=LFI%NULOUT,FMT=9007)LFI%JPNBIC
        WRITE (UNIT=LFI%NULOUT,FMT=9010)LFI%JPNCMO
        WRITE (UNIT=LFI%NULOUT,FMT=9020)LFI%JPNCPN
        WRITE (UNIT=LFI%NULOUT,FMT=9030)LFI%JPLARD,LFI%JPLARC
        WRITE (UNIT=LFI%NULOUT,FMT=9035)LFI%JPFACX
        WRITE (UNIT=LFI%NULOUT,FMT=9037)LFI%JPXUFM
        WRITE (UNIT=LFI%NULOUT,FMT=9040)LFI%JPRECL
        WRITE (UNIT=LFI%NULOUT,FMT=9050)LFI%JPNXFI
        WRITE (UNIT=LFI%NULOUT,FMT=9060)LFI%JPNPIA
        WRITE (UNIT=LFI%NULOUT,FMT=9070)LFI%JPNXPI
        WRITE (UNIT=LFI%NULOUT,FMT=9080)LFI%JPNPIS
        WRITE (UNIT=LFI%NULOUT,FMT=9090)LFI%JPNAPP
        WRITE (UNIT=LFI%NULOUT,FMT=9100)LFI%JPLDOC
        WRITE (UNIT=LFI%NULOUT,FMT=9110)LFI%JPNPDF
        WRITE (UNIT=LFI%NULOUT,FMT=9120)LFI%JPNXPR
        WRITE (UNIT=LFI%NULOUT,FMT=9130)LFI%JPNIL
        WRITE (UNIT=LFI%NULOUT,FMT=9135)LFI%JPLFTX
        WRITE (UNIT=LFI%NULOUT,FMT=9137)LFI%JPLFIX
      ELSE
        WRITE (UNIT=LFI%NULOUT,FMT=9006)LFI%JPNBIM
        WRITE (UNIT=LFI%NULOUT,FMT=9008)LFI%JPNBIC
        WRITE (UNIT=LFI%NULOUT,FMT=9011)LFI%JPNCMO
        WRITE (UNIT=LFI%NULOUT,FMT=9021)LFI%JPNCPN
        WRITE (UNIT=LFI%NULOUT,FMT=9031)LFI%JPLARD,LFI%JPLARC
        WRITE (UNIT=LFI%NULOUT,FMT=9036)LFI%JPFACX
        WRITE (UNIT=LFI%NULOUT,FMT=9038)LFI%JPXUFM
        WRITE (UNIT=LFI%NULOUT,FMT=9041)LFI%JPRECL
        WRITE (UNIT=LFI%NULOUT,FMT=9051)LFI%JPNXFI
        WRITE (UNIT=LFI%NULOUT,FMT=9061)LFI%JPNPIA
        WRITE (UNIT=LFI%NULOUT,FMT=9071)LFI%JPNXPI
        WRITE (UNIT=LFI%NULOUT,FMT=9081)LFI%JPNPIS
        WRITE (UNIT=LFI%NULOUT,FMT=9091)LFI%JPNAPP
        WRITE (UNIT=LFI%NULOUT,FMT=9101)LFI%JPLDOC
        WRITE (UNIT=LFI%NULOUT,FMT=9111)LFI%JPNPDF
        WRITE (UNIT=LFI%NULOUT,FMT=9121)LFI%JPNXPR
        WRITE (UNIT=LFI%NULOUT,FMT=9131)LFI%JPNIL
        WRITE (UNIT=LFI%NULOUT,FMT=9136)LFI%JPLFTX
        WRITE (UNIT=LFI%NULOUT,FMT=9138)LFI%JPLFIX
      ENDIF
C
C         Envoi d'un message terminal.
C
      WRITE (UNIT=LFI%NULOUT,FMT='(/)')
C
      IF (LFI%LFRANC) THEN
        CLMESS='Fin d''impression de la '//
     S             'Configuration du logiciel LFI'
      ELSE
        CLMESS='End of dump of LFI software configuration'
      ENDIF
C
      CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,LLFATA,
     S                CLMESS,CLNSPR,CLACTI)
      WRITE (UNIT=LFI%NULOUT,FMT='(///)')
C
      IF (LHOOK) CALL DR_HOOK('LFICFG_MT',1,ZHOOK_HANDLE)
      RETURN
C
 9005 FORMAT(' Nombre de Bits par mot machine.......................',
     S       I5)
C
 9006 FORMAT(' Number of Bits per computer word.....................',
     S       I5)
C
 9007 FORMAT(' Nombre de Bits par caractere machine.................',
     S       I5)
C
 9008 FORMAT(' Number of Bits per computer character................',
     S       I5)
C
 9010 FORMAT('    ... d''ou Nombre de caracteres par mot.............',
     S       I5)
C
 9011 FORMAT('    ... so Number of characters per word is...........',
     S       I5)
C
 9020 FORMAT(' Nombre maximum de Caracteres par NOM d''ARTICLE.......',
     S       I5)
C
 9021 FORMAT(' Maximum number of Characters per RECORD NAME.........',
     S       I5)
C
 9030 FORMAT(' Longueur d''Article PHYSIQUE Elementaire des Fichiers.',
     S       I5,' Mots (',I5,' Caracteres )')
C
 9031 FORMAT(' Elementary PHYSICAL Record Length of Files...........',
     S       I5,' Words (',I5,' Characters )')
C
 9035 FORMAT(' Facteur Multiplicatif Maximal de cette Longueur......',
     S       I5)
C
 9036 FORMAT(' Maximum Multiply Factor of this elementary Length....',
     S       I5)
C
 9037 FORMAT(' Nombre maximum d''associations explicites',/,
     S       '   entre Unites Logiques et Facteurs Multiplicatifs...',
     S       I5)
C
 9038 FORMAT(' Maximum number of explicit associations',/,
     S       '         between Logical Units and Multiply Factors...',
     S       I5)
C
 9040 FORMAT(' "RECL" Elementaire pour l''"OPEN" des Fichiers........',
     S       I5)
C
 9041 FORMAT(' Elementary "RECL" parameter for "OPEN" of Files......',
     S       I5)
C
 9050 FORMAT(' Nombre maximum de Fichiers ouverts en meme temps',/,
     S       '        (si leur "facteur multiplicatif" vaut 1 ).....',
     S       I5)
C
 9051 FORMAT(' Maximum number of Files open at the same time',/,
     S       '               (if their "multiply factor" is 1 ).....',
     S       I5)
C
 9060 FORMAT(' Nombre de *PAIRES* de "PAGES d''INDEX"',/,
     S       '        (en memoire) *PREALLOUEES* par Fichier........',
     S       I5)
C
 9061 FORMAT(' Number of *PREALLOCATED PAIRS* of "INDEX PAGES"',/,
     S       '        (in software commons) per File................',
     S       I5)
C
 9070 FORMAT(' Nombre TOTAL de *PAIRES* de "PAGES d''INDEX"',/,
     S       '        (en memoire) ALLOUABLES.......................',
     S       I5)
C
 9071 FORMAT(' TOTAL number of *PAIRS* of "INDEX PAGES"',/,
     S       '        (in software commons) allocatable.............',
     S       I5)
C
 9080 FORMAT('    ... d''ou Nombre de P.P.I. non preallouees.........',
     S       I5)
C
 9081 FORMAT('    ... so Number of P.I.P. not preallocated is.......',
     S       I5)
C
 9090 FORMAT(' Nombre Maxi. utile de NOMS d''ARTICLES',/,
     S       '        par PAGE ou ARTICLE d''INDEX Elementaire.......',
     S       I5)
C
 9091 FORMAT(' Maximum number of usable RECORD NAMES',/,
     S       '        per Elementary INDEX PAGE or RECORD...........',
     S       I5)
C
 9100 FORMAT(' Longueur de la Partie DOCUMENTAIRE du 1er Article....',
     S       I5,' Mots')
C
 9101 FORMAT(' Length of DOCUMENTARY part in 1st physical record....',
     S       I5,' Words')
C
 9110 FORMAT(' Nombre de PAGES de DONNEES par Fichier...............',
     S       I5)
C
 9111 FORMAT(' Number of DATA PAGES per FILE........................',
     S       I5)
C
 9120 FORMAT(' Nombre Maximum de PAIRES d''ARTICLES d''INDEX',/,
     S       '        Reservables a la Creation d''un Fichier........',
     S       I5)
C
 9121 FORMAT(' Maximum number of PAIRS of INDEX RECORDS',/,
     S       '        that may be RESERVED at File Creation.........',
     S       I5)
C
 9130 FORMAT(' Code "VALEUR ABSENTE"', /,
     S       '      pour certaines Tables D''ENTIERS.................',
     S       I5)
C
 9131 FORMAT(' "MISSING VALUE" code', /,
     S       '      for some INTEGER-type Tables....................',
     S       I5)
C
 9135 FORMAT(' Longueur Maximale des Noms de Fichiers',/,
     S       '                 traitable proprement.................',
     S       I5,' Caracteres')
C
 9136 FORMAT(' Maximum Length of File Names',/,
     S       '                    carefully handled.................',
     S       I5,' Characters')
C
 9137 FORMAT(' Longueur Maximale imprimable des Noms de Fichiers....',
     S       I5,' Caracteres')
C
 9138 FORMAT(' Maximum printable Length of File Names (most cases)..',
     S       I5,' Characters')
C
      END

