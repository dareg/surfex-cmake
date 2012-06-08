C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAISC1_MT (FA,  KREP, KRANG )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Ce sous-programme initialise un tableau "reference" de
C      l'en-tete GRIB, section 1.
C      (routine appelee une seule fois pour un fichier donne)
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C*
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANG
C
      INTEGER I, INUMER, INIMES, IRANGC
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     0.  -  INITIALISATIONS ET CONTROLES
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAISC1_MT',0,ZHOOK_HANDLE)
      KREP=0
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
      IRANGC=FA%NUCADR(KRANG)
C**
C     1.  -  INIT. DU TAB. FA%NSEC1 REPRESENTANT LA SECTION 1 DE GRIBEX
C-----------------------------------------------------------------------
C
C  2: identification of centre
C
C (defaut=85 pour Toulouse; pour en changer, utiliser FAREGU)
      FA%NSEC1(2,KRANG) = 85
C  3: generating process identification number, alloc by the orig. centre
      IF (FA%LIMLAM(IRANGC)) THEN
C Il s'agit du modele Aladin
        FA%NSEC1(3,KRANG) = 177
      ELSE
C Il s'agit du modele Arpege
        IF (FA%MADATE(7,KRANG).GT.0) THEN
C     prevision
          FA%NSEC1(3,KRANG) = 211
        ELSE
C     analyse
          FA%NSEC1(3,KRANG) = 201
C     analyse initialisee -> prevision
          IF (FA%MADATE(9,KRANG).EQ.1) FA%NSEC1(3,KRANG) = 211
        ENDIF
      ENDIF
C  4: grid definition
C     =255 for a non-catalogued grid (description follows in KSEC2)
      FA%NSEC1(4,KRANG) = 255
C  5: flag showing whether sections 2 and 3 are present
C     128 --> Section 2 is included, Section 3 is omitted (no bitmap)
      FA%NSEC1(5,KRANG) = 128
C  6 a 9: to be initialized later (specific to each field)
      FA%NSEC1(6:9,KRANG) = 0
C 10 a 21: valeurs deduites de FA%MADATE(1:11,KRANG)
C
C rappel :   2000 -->  an 100 (=FA%NSEC1(10))  siecle 20 (=FA%NSEC1(21))
C            2001 -->  an   1  "        "   siecle 21  "        "
C
      FA%NSEC1(10,KRANG) = 1 + MOD(FA%MADATE(1,KRANG) - 1 , 100)
      FA%NSEC1(21,KRANG) = 1 + (FA%MADATE(1,KRANG) - 1) / 100
      DO I=1,10
        FA%NSEC1(10+I,KRANG) = FA%MADATE(1+I,KRANG)
      ENDDO
C FA%NSEC1(18,KRANG)=10 signifie un codage sur 2 octets de l'echeance
C Ce n'est pas le cas, donc on revient a 0 pour GRIBEX
      IF (FA%NSEC1(18,KRANG)==10) FA%NSEC1(18,KRANG)=0
C FA%MADATE(10,KRANG) peut contenir l'echeance precedente dans le cas
C d'un calcul sur une periode (min, max par exemple): c'est une
C convention dans FA (depuis fin 2000) qui est incompatible avec
C GRIBEX. On retire donc cette valeur ici dans FA%NSEC1, sachant
C qu'elle sera utilisee plus tard.
      FA%NSEC1(19,KRANG) = 0
C
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      LLFATA=LLMOER (KREP,KRANG)
C
      IF (FA%LFAMOP.OR.LLFATA) THEN
        INIMES=2
        CLNSPR='FAISC1'
        INUMER=FA%JPNIIL
C
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4)')
     S     KREP, KRANG
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR, CLNSPR,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAISC1_MT',1,ZHOOK_HANDLE)
      END

