C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAREGI_MT (FA,  CDCLEF, KVAL, KOPT )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Ce sous-programme controle (lecture/ecriture) les options
C      de compression de GRIBEX implicites.
C      (REGlage des options Implicites de codage de gribex)
C**
C    Arguments : CDCLEF (Entree) ==> Mot clef precisant l'action a faire;
C                KVAL   (Sortie  ==> Valeur lue ou ecrite;
C                      ou Entree)
C                KOPT   (Entree) ==> Flag: 0->lecture 1->ecriture;
C*
C     Signification des divers elements du tableau FA%NCODGRI
C
C     FA%NCODGRI(1) = type de codage (0->option HOPER='C'; 1->option HOPER='K')
C     FA%NCODGRI(2) = KSEC4(6), indicateur de la presence de flags
C                               additionnels (0->non; 16->oui)
C     FA%NCODGRI(3) = KSEC4(7)
C     FA%NCODGRI(4) = KSEC4(9), indicateur de la presence de bitmaps
C                             secondaires (0->non; 32->oui)
C     FA%NCODGRI(5) = KSEC4(10), indicateur pour le nb de bits des
C                        groupes de pts de grille (0->const.; 16->different)
C     FA%NCODGRI(6) = KSEC4(11), nb de bits pour les groupes de pts de
C                      grille, quand il est constant.
C                        Si negatif, le logiciel calcule un nb optimal a partir
C                        de -KSEC4(11).
C     FA%NCODGRI(7) = KSEC4(12), indicateur pour les extensions generales de la
C                        compression (0->non; 8->oui)
C     FA%NCODGRI(8) = KSEC4(13), indicateur pour le rearrangement boustrophedo
C                        nique (0->non; 4->oui)
C     FA%NCODGRI(9) = KSEC4(14) (valeurs possibles: -1, 0 et 2)
C     FA%NCODGRI(10) = KSEC4(15) (valeurs possibles: -1, 0 et 1), sert avec
C                         KSEC4(14) a definir la technique de la difference
C                         spatiale.Si l'un des 2 est negatif, l'ordre de
C                         differentiation est estime dynamiquement, sinon
C                         l'ordre vaut KSEC4(14)+KSEC4(15)
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KNUMER, KVAL, KOPT
C
      CHARACTER*5 CDCLEF
C
      INTEGER INUMER, INIMES, IREP, INBITSMAX
C
      LOGICAL LLVERG
#include "facom2.h"
#include "facom_mt.h"
C
C
C
C**
C     0.  -  INITIALISATIONS ET ALLOCATIONS PREALABLES
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAREGI_MT',0,ZHOOK_HANDLE)
      IREP=0
      IF (FA%FAREGI_LLPREA) THEN
C
C          A la premiere utilisation, appel au sous-programme "FARINE".
C
        CALL FARINE_MT(FA, 2)
        FA%FAREGI_LLPREA=.FALSE.
      ENDIF
      INBITSMAX=MAX(FA%NBIPDG,FA%NBICSP)
C
C**
C     1.  -  SPECIFICATION D'UN NOUVEAU CODAGE
C-----------------------------------------------------------------------
C
      IF (KOPT==1) THEN
C
C   Pas de compression (sous-tronc et puissance laplacien
C   ajoutees systematiquement + tard pour les coeff spectraux)
C
        IF (CDCLEF=='BASIC'.OR.CDCLEF=='basic') THEN
          FA%NCODGRI(1)=0
          FA%NCODGRI(2)=0
          FA%NCODGRI(3)=0
          FA%NCODGRI(4)=0
          FA%NCODGRI(5)=0
          FA%NCODGRI(6)=0
C
C   Comme le "BASIC" avec une compression
C   ligne a ligne pour les points de grille
C
        ELSEIF (CDCLEF=='PACK1'.OR.CDCLEF=='pack1') THEN
          FA%NCODGRI(1)=0
          FA%NCODGRI(2)=16
          FA%NCODGRI(3)=0
          FA%NCODGRI(4)=0
          FA%NCODGRI(5)=16
          FA%NCODGRI(6)=0
C
C   Comme le "BASIC" avec une compression
C   pour les points de grille ou le nb de bits est le meme
C   dans chaque groupe de points de grille
C
        ELSEIF (CDCLEF=='PACK2'.OR.CDCLEF=='pack2') THEN
          FA%NCODGRI(1)=0
          FA%NCODGRI(2)=16
          FA%NCODGRI(3)=0
          FA%NCODGRI(4)=32
          FA%NCODGRI(5)=0
C Un nb de bits optimal sera recherche par le logiciel
          FA%NCODGRI(6)=-99
C
C   Comme le "BASIC" avec une compression
C   "OMM general" pour les points de grille
C
        ELSEIF (CDCLEF=='PACK3'.OR.CDCLEF=='pack3') THEN
          FA%NCODGRI(1)=0
          FA%NCODGRI(2)=16
          FA%NCODGRI(3)=0
          FA%NCODGRI(4)=32
          FA%NCODGRI(5)=16
          FA%NCODGRI(6)=0
C
C   Compression "aggressive": le logiciel va
C   tenter la compression ligne a ligne et comparer la
C   longueur de message obtenue avec celle en l'absence de
C   compression et au final retenir la meilleure solution.
C
        ELSEIF (CDCLEF=='APAC1'.OR.CDCLEF=='apac1') THEN
          FA%NCODGRI(1)=1
          FA%NCODGRI(2)=16
          FA%NCODGRI(3)=0
          FA%NCODGRI(4)=0
          FA%NCODGRI(5)=16
          FA%NCODGRI(6)=0
C
C   Compression "aggressive": le logiciel va
C   suivre la demarche "APAC1" en testant en plus le nb de bits
C   constant par groupe de pts de grille et retenir la meilleure
C   compression
C
        ELSEIF (CDCLEF=='APAC2'.OR.CDCLEF=='apac2') THEN
          FA%NCODGRI(1)=1
          FA%NCODGRI(2)=16
          FA%NCODGRI(3)=0
          FA%NCODGRI(4)=0
          FA%NCODGRI(5)=0
C Un nb de bits optimal sera recherche par le logiciel
          FA%NCODGRI(6)=-99
C
C   Compression "aggressive": le logiciel va
C   suivre la demarche "APAC1" en testant en plus la compression
C   "general OMM" et retenir la meilleure compression.
C
        ELSEIF (CDCLEF=='APAC3'.OR.CDCLEF=='apac3') THEN
          FA%NCODGRI(1)=1
          FA%NCODGRI(2)=16
          FA%NCODGRI(3)=0
          FA%NCODGRI(4)=32
          FA%NCODGRI(5)=16
          FA%NCODGRI(6)=0
C
C   Compression "aggressive": le logiciel va
C   suivre la demarche "APAC3" en testant en plus le nb de bits
C   constant par groupe de pts de grille et retenir la meilleure
C   compression.
C
        ELSEIF (CDCLEF=='APAC4'.OR.CDCLEF=='apac4') THEN
          FA%NCODGRI(1)=1
          FA%NCODGRI(2)=16
          FA%NCODGRI(3)=0
          FA%NCODGRI(4)=32
          FA%NCODGRI(5)=0
C Un nb de bits optimal sera recherche par le logiciel
          FA%NCODGRI(6)=-99
C
C   Specification du nb de bits a utiliser dans le cadre
C   de la compression avec nb de bits constant
C   par groupe de pts de grille
C   
        ELSEIF (CDCLEF=='WIDPA'.OR.CDCLEF=='widpa') THEN
          IF (KVAL.LT.1-INBITSMAX.OR.KVAL.GT.INBITSMAX-1) THEN
            IREP=-97
            GOTO 1001
          ENDIF
          FA%NCODGRI(6)=KVAL
C
C   Demande supplementaire de l'extension generale de la compression
C   (si KVAL=1, sinon c'est le retrait de cette option)
C
        ELSEIF (CDCLEF=='GEXTE'.OR.CDCLEF=='gexte') THEN
          IF (KVAL.EQ.1) THEN
            CALL GRSX2O(1)
            FA%NCODGRI(7)=8
          ELSE
            FA%NCODGRI(7)=0
          ENDIF
C
C   Demande supplementaire du rearrangement boustrophedonique dans la
C   compression (si KVAL=1, sinon c'est le retrait de cette option)
C
        ELSEIF (CDCLEF=='BOUST'.OR.CDCLEF=='boust') THEN
          IF (KVAL.EQ.1) THEN
            CALL GRSX2O(1)
            FA%NCODGRI(8)=4
          ELSE
            FA%NCODGRI(8)=0
          ENDIF
C
C   Demande supplementaire de la difference spatiale dans la compression.
C   KVAL donne l'ordre de differentiation:
C   -1-> calcul dynamique par GRIBEX; 1 a 3->ordre; 0->desactiv; autre->err
C
        ELSEIF (CDCLEF=='DIFFE'.OR.CDCLEF=='diffe') THEN
          IF (KVAL.EQ.-1) THEN
            CALL GRSX2O(1)
            CALL GRSN2O(1)
            FA%NCODGRI( 9)=0
            FA%NCODGRI(10)=-1
          ELSEIF (KVAL.EQ.1) THEN
            CALL GRSX2O(1)
            CALL GRSN2O(1)
            FA%NCODGRI( 9)=0
            FA%NCODGRI(10)=1
          ELSEIF (KVAL.EQ.2) THEN
            CALL GRSX2O(1)
            CALL GRSN2O(1)
            FA%NCODGRI( 9)=2
            FA%NCODGRI(10)=0
          ELSEIF (KVAL.EQ.3) THEN
            CALL GRSX2O(1)
            CALL GRSN2O(1)
            FA%NCODGRI( 9)=2
            FA%NCODGRI(10)=1
          ELSEIF (KVAL.EQ.0) THEN
            CALL GRSN2O(0)
            FA%NCODGRI( 9)=0
            FA%NCODGRI(10)=0
          ELSE
            IREP=-125
            GOTO 1001
          ENDIF
        ENDIF
C**
C     2.  -  DEMANDE D'INFORMATION
C-----------------------------------------------------------------------
C
      ELSEIF (KOPT==0) THEN
C
C   Obtention des mots-clef disponibles
C
        IF (CDCLEF=='CLEFS'.OR.CDCLEF=='clefs') THEN
          KVAL=0
          WRITE (UNIT=FA%NULOUT,FMT=*)
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'Mots clef disponibles pour FAREGI:'
          WRITE (UNIT=FA%NULOUT,FMT=*)
          WRITE (UNIT=FA%NULOUT,FMT=*)'BASIC: pas de compression'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'PACK1: BASIC avec une compression'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       ligne a ligne pour les pts de grille'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'PACK2: BASIC avec une compression'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       ou le nb de bits est cst pour les groupes'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'PACK3: BASIC avec une compression'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       OMM general pour les points de grille'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'APAC1: compression agressive'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       BASIC et PACK1 sont testes'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'APAC2: compression agressive'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       BASIC, PACK1 et PACK2 sont testes'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'APAC3: compression agressive'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       BASIC, PACK1 et PACK3 sont testes'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'APAC4: compression agressive'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       BASIC, PACK1, PACK2 et PACK3 testes'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'WIDPA: lecture/ecriture du nb de bits'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       a utiliser pour les groupes de points'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       de grille dans le cas PACK2'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'GEXTE: les extensions generales de la compression'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       sont activees (KVAL=1) ou desactiv. (KVAL=0)'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'BOUST: le rearrangement boustrophedonique est'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       active (KVAL=1) ou desactive (KVAL=0)'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'DIFFE: la differenciation spatiale est'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       activee (KVAL=ordre de differ. (1 a 3)'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       ou -1 (calcul dyn)) ou desactivee (0)'
          WRITE (UNIT=FA%NULOUT,FMT=*)
C
C   Lecture du nb de bits a utiliser dans le cadre
C   de la compression avec nb de bits constant
C   par groupe de pts de grille
C   
        ELSEIF (CDCLEF=='WIDPA'.OR.CDCLEF=='widpa') THEN
          KVAL=FA%NCODGRI(6)
C
C   Lecture de la presence ou non de la compression
C   "general extended"
C
        ELSEIF (CDCLEF=='GEXTE'.OR.CDCLEF=='gexte') THEN
          KVAL = FA%NCODGRI(7) / 8
C
C   Lecture de la presence ou non du rearrangement
C   boustrophedonique.
C
        ELSEIF (CDCLEF=='BOUST'.OR.CDCLEF=='boust') THEN
          KVAL = FA%NCODGRI(8) / 4
C
C   Lecture de la presence ou non de la differentiation spatiale
C
        ELSEIF (CDCLEF=='DIFFE'.OR.CDCLEF=='diffe') THEN
          KVAL=FA%NCODGRI( 9)+FA%NCODGRI(10)
        ELSE
          IREP=-125
          GOTO 1001
        ENDIF
C**
C     3.  -  OPTION INCONNUE
C-----------------------------------------------------------------------
C
      ELSE
        IREP=-125
        GOTO 1001
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      LLFATA=LLMOER (IREP,0)
C
      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=FA%NIMSGA
      ENDIF
C
      IF (INIMES.EQ.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('FAREGI_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FAREGI'
      INUMER=FA%JPNIIL
C
      WRITE (UNIT=CLMESS,FMT='(''IREP='',I2,
     S        '', CDCLEF='''''',A,'''''', KVAL='',I12,
     S        '', KOPT='',I4)')
     S        IREP,CDCLEF,KVAL,KOPT
      INUMER=FA%JPNIIL
      CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CLNSPR,.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FAREGI_MT',1,ZHOOK_HANDLE)
      END

