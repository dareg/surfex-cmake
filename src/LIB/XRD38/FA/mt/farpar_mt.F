C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FARPAR_MT (FA,  KREP, CDPREF, CDSUFF, KCODPA, KNUM )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme de reglage de la correspondance "nom d'article FA"
C              <-> "descripteurs GRIB du parametre+niveau"
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C    ( Tableau ) CDPREF (Entree) ==> Prefixe pour les KNUM noms d'article;
C    ( Tableau ) CDSUFF (Entree) ==> Suffixe pour les KNUM noms d'article;
C    ( Tableau ) KCODPA (Entree) ==> 6 descripteurs GRIB pour chacun
C                                    des KNUM parametres:
C      KCODPA(J,1) = KSEC1(1)  version de la table parametres
C      KCODPA(J,2) = KSEC1(6)  indicateur du parametre
C      KCODPA(J,3) = KSEC1(7)  indicateur du type de niveau
C      KCODPA(J,4) = KSEC1(8)  niveau
C      KCODPA(J,5) = KSEC1(9)  2ieme nv si couche, sinon 0
C      KCODPA(J,6) = KSEC1(18) indicateur du type de champ
C                              (0 sf si min/max:2 ou si cumul:4)
C
C                KNUM   (Entree  ==> Nombre de parametres a regler
C                          et        (dimension de CDPREF, CDSUFF et KCODPA)
C                        Sortie) ==> Nb de nouveaux parametres pouvant encore
C                                    etre definis lors d'un appel ulterieur.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUM
      INTEGER KCODPA(KNUM,6)
C
      CHARACTER (LEN=*) CDPREF(KNUM), CDSUFF(KNUM)
C
      INTEGER J, JJ, INUMER, INIMES, JMEM
#include "facom_mt.h"
C
      INTRINSIC LEN_TRIM
C
C**
C     0.  -  CONTROLES ET INITIALISATIONS PREALABLES
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FARPAR_MT',0,ZHOOK_HANDLE)
      IF (KNUM.LT.1) THEN
        KREP=-129
        IF (FA%LFAMOP) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FARPAR: Nb de parametres ',KNUM,' incorrect'
        ENDIF
        GOTO 1001
      ENDIF
      DO J = 1,KNUM
        IF ( LEN_TRIM(CDPREF(J)).LE.0 .OR.
     S       LEN_TRIM(CDPREF(J)).GT.FA%JPXPRF ) THEN
          KREP=-129
          IF (FA%LFAMOP) THEN
            WRITE (UNIT=FA%NULOUT,FMT=*)
     S              'FARPAR: Longueur du prefixe ',CDPREF(J),
     S              ' incorrecte : ',LEN_TRIM(CDPREF(J))
          ENDIF
          GOTO 1001
        ENDIF
        IF ( LEN_TRIM(CDSUFF(J)).LE.0 .OR.
     S       LEN_TRIM(CDSUFF(J)).GT.FA%JPXNOM-LEN_TRIM(CDPREF(J)) ) THEN
          KREP=-129
          IF (FA%LFAMOP) THEN
            WRITE (UNIT=FA%NULOUT,FMT=*)
     S              'FARPAR: Longueur du suffixe ',CDSUFF(J),
     S              ' incorrecte : ',LEN_TRIM(CDSUFF(J))
          ENDIF
          GOTO 1001
        ENDIF
        DO JJ = 1,3
          IF (KCODPA(J,JJ).LT.1 .OR. KCODPA(J,JJ).GT.255) THEN
            KREP=-129
            IF (FA%LFAMOP) THEN
              WRITE (UNIT=FA%NULOUT,FMT=*)
     S                'FARPAR: descripteur GRIB num ',JJ,
     S                ' pour le parametre num ',J,' ( ',
     S                CDPREF(J)//CDSUFF(J),' ) incorrect : ',
     S                KCODPA(J,JJ)
            ENDIF
            GOTO 1001
          ENDIF
        ENDDO
        IF (KCODPA(J,6).LT.0 .OR. KCODPA(J,6).GT.124) THEN
          KREP=-129
          IF (FA%LFAMOP) THEN
            WRITE (UNIT=FA%NULOUT,FMT=*)
     S              'FARPAR: descripteur GRIB, KSEC1(18),',
     S              ' pour le parametre num ',J,' ( ',
     S              CDPREF(J)//CDSUFF(J),' ) incorrect : ',
     S              KCODPA(J,6)
          ENDIF
          GOTO 1001
        ENDIF
        IF (KCODPA(J,4).LT.0) THEN
          KREP=-129
          IF (FA%LFAMOP) THEN
            WRITE (UNIT=FA%NULOUT,FMT=*)
     S              'FARPAR: descripteur GRIB, KSEC1(8),',
     S              ' pour le parametre num ',J,' ( ',
     S              CDPREF(J)//CDSUFF(J),' ) incorrect : ',
     S              KCODPA(J,4)
          ENDIF
          GOTO 1001
        ENDIF
      ENDDO
C
C**
C     2.  -  Prise en compte des nouvelles correspondances
C---------------------------------------------------------
C
C
      DO J = 1,KNUM
C Recherche prealable de l'eventuelle existence de la definition
C de ce parametre (il faudra alors l'ecraser).
        JMEM = 0
        DO JJ = 1,FA%NBPARC
          IF (CDPREF(J )(1:LEN_TRIM(CDPREF(J ))).EQ.
     S        FA%CIPREF(JJ)(1:LEN_TRIM(FA%CIPREF(JJ))) .AND.
     S        CDSUFF(J )(1:LEN_TRIM(CDSUFF(J ))).EQ.
     S        FA%CISUFF(JJ)(1:LEN_TRIM(FA%CISUFF(JJ)))) THEN
            JMEM = JJ
            EXIT
          ENDIF
        ENDDO
        IF (JMEM==0) THEN
          IF (FA%NBPARC+1.GT.FA%JPXPAR) THEN
            KREP=-129
            IF (FA%LFAMOP) THEN
              WRITE (UNIT=FA%NULOUT,FMT=*)
     S                'FARPAR: le nb lim de parametres definissables ',
     S                FA%JPXPAR
              WRITE (UNIT=FA%NULOUT,FMT=*)
     S                '        est depasse avec la definition de ',
     S                KNUM,' nouveaux parametres'
            ENDIF
            GOTO 1001
          ENDIF
          JMEM = FA%NBPARC + 1
          FA%NBPARC = FA%NBPARC + 1
        ENDIF
        FA%CIPREF(JMEM)     = 
     S         CDPREF(J)(1:LEN_TRIM(CDPREF(J)))
        FA%CISUFF(JMEM)     = 
     S         CDSUFF(J)(1:LEN_TRIM(CDSUFF(J)))
        FA%NCODPA(JMEM,1:6) = KCODPA(J,1:6)
        IF (FA%LFAMOP) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FARPAR: Prise en compte de ',CDPREF(J)//CDSUFF(J)
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           '        associe a KSEC1(1,6:9 et 18) = ',
     S           FA%NCODPA(JMEM,1:6)
        ENDIF
      ENDDO
      KNUM = FA%JPXPAR - FA%NBPARC
C
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
C
      IF (FA%LFAMOP) THEN
        INIMES=2
        CLNSPR='FARPAR'
        INUMER=FA%JPNIIL
C
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4)') KREP
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLNSPR,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FARPAR_MT',1,ZHOOK_HANDLE)
      END

