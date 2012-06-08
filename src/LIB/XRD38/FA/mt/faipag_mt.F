C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAIPAG_MT (FA,  KREP,   KRANG,  CDPREF, KNIVAU, CDSUFF,
     S                      KNIPAR )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C      Initialisation de quelques descripteurs de l'entete Gribex
C      (section 1) relatifs au parametre, a partir du nom FA du champ.
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
C                KNIVAU (Entree) ==> Niveau vertical eventuel;
C                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
C    ( Tableau ) KNIPAR (Sortie) ==> quelques descripteurs de la section 1 de
C                                    GRIBEX (KNIPAR(1)  =KSEC1(1),
C                                            KNIPAR(2:5)=KSEC1(6:9) ) et un
C                                    indicateur de type de champ (KNIPAR(6)=
C                                    KSEC1(18)):0->RAS; 2->min/max; 4->cumul
C
C     Original  : 06 juillet 2004, Denis Paradis DSI/DEV
C     --------
C
C     Modifications
C     -------------
C       R. El Ouaraini: 03-Oct-2006, enlever le commentaire de l'initialisation de KNIPAR(3) pour
C                  les types de niveaux :  hauteur, iso-tourb potent, isentrope et modele.
C
C*
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANG, KNIVAU, KNIPAR(6)
C
      CHARACTER (LEN=*) CDPREF, CDSUFF
C
      INTEGER INUMER, INIMES, J
C
      INTRINSIC LEN_TRIM
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     0.  -  INITIALISATIONS PREALABLES
C-----------------------------------------------------------------------
C
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAIPAG_MT',0,ZHOOK_HANDLE)
      KREP=0
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C
C  DEFAUTS:
C
C Numero de version de la table de code parametres
      KNIPAR(1)=1
C Indicateur de parametre (255=> valeur manquante)
      KNIPAR(2)=255
C Indicateur de type de niveau (1=> surface)
      KNIPAR(3)=1
C Niveau 1, Niveau 2 et type de champs
      KNIPAR(4:6)=0
C**
C     1.  -  UTILISATION DE LA TABLE DE CORRESPONDANCE
C-----------------------------------------------------------------------
C
C
      JMEM = 0
      DO J = 1,FA%NBPARC
        IF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.
     S      FA%CIPREF(J)(1:LEN_TRIM(FA%CIPREF(J))) .AND.
     S      CDSUFF(1:LEN_TRIM(CDSUFF)).EQ.
     S      FA%CISUFF(J)(1:LEN_TRIM(FA%CISUFF(J)))) THEN
          JMEM = J
          EXIT
        ENDIF
      ENDDO
      IF (FA%LFAMOP.AND.JMEM.EQ.0) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         'FAIPAG: WARNING, pas de reference GRIB pour ',
     S         CDPREF(1:LEN_TRIM(CDPREF))//CDSUFF(1:LEN_TRIM(CDSUFF))
        WRITE (UNIT=FA%NULOUT,FMT=*)'       Les defauts seront utilises'
        GOTO 1001
      ENDIF
      IF (JMEM.NE.0) KNIPAR(1:6) = FA%NCODPA(JMEM,1:6)
C**
C     2.  -  INITIALISATION DU NIVEAU (AUTRE QUE 0)
C--------------------------------------------------
C
C     2.1 -  Champs sur un niveau isobare
C
      IF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.'P') THEN
C  La pression est sur 5 chiffres: on la ramene a l'hPa
C     et on recree le niveau 1000 hPa
C
C  Si KNIVAU < 100, la pression fait moins d'un hPa et
C  on utilise une extension du GRIB introduite par le CEP:
C       KSEC1(7) = 210 (au lieu de 100)
C    et KSEC1(8) = pression en Pa
C
        IF (KNIVAU .GE. 100) THEN
          KNIPAR(3)=100
          KNIPAR(4)=KNIVAU/100
        ELSEIF (KNIVAU==0) THEN
          KNIPAR(3)=100
          KNIPAR(4)=1000
        ELSE
          KNIPAR(3)=210
          KNIPAR(4)=KNIVAU
        ENDIF
C
C     2.2 -  Champs sur un niveau hauteur
C
      ELSEIF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.'H') THEN
       KNIPAR(3)=105
        KNIPAR(4)=KNIVAU
C
C     2.3 -  Champs sur un niveau iso-tourbillon-potentiel
C
C            ( unite SI = K m2 s-1 kg-1 = 10+6 PVU
C              mais l'unite retenu est le milliPVU: 10-9 SI)
      ELSEIF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.'V') THEN
       KNIPAR(3)=117
C KNIVAU est exprime en 1/10 PVU
        KNIPAR(4)=KNIVAU*100
        IF (KNIVAU==0) KNIPAR(4)=1000
C
C     2.4 -  Champs sur un niveau isentrope
C
      ELSEIF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.'T') THEN
       KNIPAR(3)=113
        KNIPAR(4)=KNIVAU
C
C     2.5 -  Champs sur un niveau modele
C
      ELSEIF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.'S') THEN
       KNIPAR(3)=109
        KNIPAR(4)=KNIVAU
      ENDIF
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
        CLNSPR='FAIPAG'
        INUMER=FA%JPNIIL
C
        WRITE (UNIT=FA%NULOUT,FMT=*)'FAIPAG: KNIPAR(1:6) = ',KNIPAR(1:6)
        WRITE (UNIT=FA%NULOUT,FMT=*)
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4,
     S       '', CDPREF='''''',A,'''''', KNIVAU='',I6,
     S       '', CDSUFF='''''',A,'''''''')')
     S     KREP, KRANG, CDPREF(1:LEN_TRIM(CDPREF)), KNIVAU,
     S     CDSUFF(1:LEN_TRIM(CDSUFF))
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CDPREF,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAIPAG_MT',1,ZHOOK_HANDLE)
      END

