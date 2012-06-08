C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAREGU_MT (FA,  KNUMER, CDCLEF, KVAL, KOPT )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Ce sous-programme controle (lecture/ecriture) les options
C      de compression de GRIBEX, pour chacune des unites logiques et
C      certains descripteurs GRIB communs a l'unite logique.
C      (REGLAGE des options de codage de gribex pour une Unite)
C**
C    Arguments : KNUMER (Entree) ==> Numero de l'unite logique;
C                CDCLEF (Entree) ==> Mot clef precisant l'action a faire;
C                KVAL   (Sortie  ==> Valeur lue ou a ecrire;
C                      ou Entree)
C                KOPT   (Entree) ==> Flag: 0->lecture 1->ecriture;
C*
C     Signification des divers elements du tableau FA%NCOGRIF
C
C     FA%NCOGRIF(1,IRANG) = type de codage (0->option HOPER='C'
C                                        1->option HOPER='K')
C     FA%NCOGRIF(2,IRANG) = KSEC4(6), indicateur de la presence de flags
C                        additionnels (0->non; 16->oui)
C     FA%NCOGRIF(3,IRANG) = KSEC4(7)
C     FA%NCOGRIF(4,IRANG) = KSEC4(9), indicateur de la presence de bitmaps
C                        secondaires (0->non; 32->oui)
C     FA%NCOGRIF(5,IRANG) = KSEC4(10), indicateur pour le nb de bits des
C                        groupes de pts de grille (0->const.; 16->different)
C     FA%NCOGRIF(6,IRANG) = KSEC4(11), nb de bits pour les groupes de pts de grille
C                        quand il est constant.
C                        Si negatif, le logiciel calcule un nb optimal a partir
C                        de -KSEC4(11).
C     FA%NCOGRIF(7,IRANG) = KSEC4(12), indicateur pour les extensions generales de
C                        la compression (0->non; 8->oui)
C     FA%NCOGRIF(8,IRANG) = KSEC4(13), indicateur pour le rearrangement boustrophedo
C                        nique (0->non; 4->oui)
C     FA%NCOGRIF(9,IRANG) = KSEC4(14) (valeurs possibles: -1, 0 et 2)
C     FA%NCOGRIF(10,IRANG) = KSEC4(15) (valeurs possibles: -1, 0 et 1), sert avec
C                         KSEC4(14) a definir la technique de la difference
C                         spatiale. Si l'un des 2 est negatif, l'ordre de
C                         differentiation est estime dynamiquement, sinon
C                         l'ordre = KSEC4(14)+KSEC4(15)
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KNUMER, KVAL, KOPT
C
      CHARACTER*5 CDCLEF
C
      INTEGER IRANG, INIMES, IREP, INBITSMAX
C
      LOGICAL LLVERG
C
#include "facom2.h"
#include "facom_mt.h"
C
C**
C     0.  -  INITIALISATIONS ET ALLOCATIONS PREALABLES
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAREGU_MT',0,ZHOOK_HANDLE)
      IREP=0
C
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ENDIF
C
C Appel prealable a FAISC1 pour initialiser FA%NSEC1(2:21,IRANG).
C On le fait ici plutot que dans FAINIG pour ne pas ecraser
C les eventuelles modifs apportees par IDCEN et/ou IDMOD
C
      IF (FA%LISEC1(IRANG)) THEN
        CALL FAISC1_MT(FA, IREP,IRANG)
        IF (IREP.NE.0) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FAREGU: ERROR ',IREP,' dans appel a FAISC1 !!'
          GOTO 1001
        ENDIF
        FA%LISEC1(IRANG)=.FALSE.
      ENDIF
C
      INBITSMAX=MAX(FA%NBFPDG(IRANG),FA%NBFCSP(IRANG))
C**
C     1.  -  SPECIFICATION D'UN NOUVEAU CODAGE
C-----------------------------------------------------------------------
C
      IF (KOPT==1) THEN
C
C   Pas de compression (sous-tronc et puissance laplacien ajoutees
C   systematiquement + tard pour les coeff spectraux)
C
        IF (CDCLEF=='BASIC'.OR.CDCLEF=='basic') THEN
          FA%NCOGRIF(1,IRANG)=0
          FA%NCOGRIF(2,IRANG)=0
          FA%NCOGRIF(3,IRANG)=0
          FA%NCOGRIF(4,IRANG)=0
          FA%NCOGRIF(5,IRANG)=0
          FA%NCOGRIF(6,IRANG)=0
C
C   Comme le "BASIC" avec une compression
C   ligne a ligne pour les points de grille
C
        ELSEIF (CDCLEF=='PACK1'.OR.CDCLEF=='pack1') THEN
          FA%NCOGRIF(1,IRANG)=0
          FA%NCOGRIF(2,IRANG)=16
          FA%NCOGRIF(3,IRANG)=0
          FA%NCOGRIF(4,IRANG)=0
          FA%NCOGRIF(5,IRANG)=16
          FA%NCOGRIF(6,IRANG)=0
C
C   Comme le "BASIC" avec une compression
C   pour les points de grille ou le nb de bits est le meme
C   dans chaque groupe de points de grille
C
        ELSEIF (CDCLEF=='PACK2'.OR.CDCLEF=='pack2') THEN
          FA%NCOGRIF(1,IRANG)=0
          FA%NCOGRIF(2,IRANG)=16
          FA%NCOGRIF(3,IRANG)=0
          FA%NCOGRIF(4,IRANG)=32
          FA%NCOGRIF(5,IRANG)=0
C Un nb de bits optimal sera recherche par le logiciel
          FA%NCOGRIF(6,IRANG)=-99
C
C   Comme le "BASIC" avec une compression general OMM
C   pour les points de grille
C
        ELSEIF (CDCLEF=='PACK3'.OR.CDCLEF=='pack3') THEN
          FA%NCOGRIF(1,IRANG)=0
          FA%NCOGRIF(2,IRANG)=16
          FA%NCOGRIF(3,IRANG)=0
          FA%NCOGRIF(4,IRANG)=32
          FA%NCOGRIF(5,IRANG)=16
          FA%NCOGRIF(6,IRANG)=0
C
C   Compression "aggressive": le logiciel va
C   tenter la compression ligne a ligne puis l'absence
C   de compression et retenir la meilleure methode
C
        ELSEIF (CDCLEF=='APAC1'.OR.CDCLEF=='apac1') THEN
          FA%NCOGRIF(1,IRANG)=1
          FA%NCOGRIF(2,IRANG)=16
          FA%NCOGRIF(3,IRANG)=0
          FA%NCOGRIF(4,IRANG)=0
          FA%NCOGRIF(5,IRANG)=16
          FA%NCOGRIF(6,IRANG)=0
C
C   Compression "aggressive": le logiciel va
C   tenter la compression type "APAC1" puis celle avec le nb de bits
C   constant par groupe de pts de grille et retenir la meilleure
C
        ELSEIF (CDCLEF=='APAC2'.OR.CDCLEF=='apac2') THEN
          FA%NCOGRIF(1,IRANG)=1
          FA%NCOGRIF(2,IRANG)=16
          FA%NCOGRIF(3,IRANG)=0
          FA%NCOGRIF(4,IRANG)=0
          FA%NCOGRIF(5,IRANG)=0
C Un nb de bits optimal sera recherche par le logiciel
          FA%NCOGRIF(6,IRANG)=-99
C
C   Compression "aggressive": le logiciel va tenter
C   la compression type "APAC1" puis la compression generale
C   OMM et retenir la meilleure
C
        ELSEIF (CDCLEF=='APAC3'.OR.CDCLEF=='apac3') THEN
          FA%NCOGRIF(1,IRANG)=1
          FA%NCOGRIF(2,IRANG)=16
          FA%NCOGRIF(3,IRANG)=0
          FA%NCOGRIF(4,IRANG)=32
          FA%NCOGRIF(5,IRANG)=16
          FA%NCOGRIF(6,IRANG)=0
C
C   Compression "aggressive": le logiciel va
C   tenter la compression type "APAC3" puis celle avec le nb de bits
C   constant par groupe de pts de grille et retenir la meilleure
C
        ELSEIF (CDCLEF=='APAC4'.OR.CDCLEF=='apac4') THEN
          FA%NCOGRIF(1,IRANG)=1
          FA%NCOGRIF(2,IRANG)=16
          FA%NCOGRIF(3,IRANG)=0
          FA%NCOGRIF(4,IRANG)=32
          FA%NCOGRIF(5,IRANG)=0
C Un nb de bits optimal sera recherche par le logiciel
          FA%NCOGRIF(6,IRANG)=-99
C
C   Specification du nb de bits a utiliser dans le cadre
C   de la compression avec nb de bits constant
C   par groupe de pts de grille
C   
        ELSEIF (CDCLEF=='WIDPA'.OR.CDCLEF=='widpa') THEN
          IF (KVAL.LT.1-INBITSMAX.OR.KVAL.GT.INBITSMAX-1) THEN
            IREP=-97
            WRITE (UNIT=FA%NULOUT,FMT='(A)')'Dans FAREGU, action WIDPA:'
            WRITE (UNIT=FA%NULOUT,FMT='(A57,I8)')
     S   '!!  ERREUR  !! Valeur incorrecte, non prise en compte: ',KVAL
            GOTO 1001
          ENDIF
          FA%NCOGRIF(6,IRANG)=KVAL
C
C   Demande supplementaire de la compression avec extension
C   generale a la norme OMM (si KVAL=1, sinon c'est le retrait de l'option)
C
        ELSEIF (CDCLEF=='GEXTE'.OR.CDCLEF=='gexte') THEN
          IF (KVAL.EQ.1) THEN
            CALL GRSX2O(1)
            FA%NCOGRIF(7,IRANG)=8
          ELSE
            FA%NCOGRIF(7,IRANG)=0
          ENDIF
C
C   Demande supplementaire du rearrangement boustrophedonique dans la compression
C   (si KVAL=1, sinon c'est le retrait de cette option)
C
        ELSEIF (CDCLEF=='BOUST'.OR.CDCLEF=='boust') THEN
          IF (KVAL.EQ.1) THEN
            CALL GRSX2O(1)
            FA%NCOGRIF(8,IRANG)=4
          ELSE
            FA%NCOGRIF(8,IRANG)=0
          ENDIF
C
C   Demande supplementaire de la difference spatiale dans la
C   compression. KVAL donne l'ordre de differentiation
C   (-1-> calcul dynamique par GRIBEX; 1 a 3->ordre; 0->desactiv; autre->err)
C
        ELSEIF (CDCLEF=='DIFFE'.OR.CDCLEF=='diffe') THEN
          IF (KVAL.EQ.-1) THEN
            CALL GRSX2O(1)
            CALL GRSN2O(1)
            FA%NCOGRIF( 9,IRANG)=0
            FA%NCOGRIF(10,IRANG)=-1
          ELSEIF (KVAL.EQ.1) THEN
            CALL GRSX2O(1)
            CALL GRSN2O(1)
            FA%NCOGRIF( 9,IRANG)=0
            FA%NCOGRIF(10,IRANG)=1
          ELSEIF (KVAL.EQ.2) THEN
            CALL GRSX2O(1)
            CALL GRSN2O(1)
            FA%NCOGRIF( 9,IRANG)=2
            FA%NCOGRIF(10,IRANG)=0
          ELSEIF (KVAL.EQ.3) THEN
            CALL GRSX2O(1)
            CALL GRSN2O(1)
            FA%NCOGRIF( 9,IRANG)=2
            FA%NCOGRIF(10,IRANG)=1
          ELSEIF (KVAL.EQ.0) THEN
            CALL GRSN2O(0)
            FA%NCOGRIF( 9,IRANG)=0
            FA%NCOGRIF(10,IRANG)=0
          ELSE
            IREP=-125
            WRITE (UNIT=FA%NULOUT,FMT='(A)') 
     S             'Dans FAREGU, action DIFFE:'
            WRITE (UNIT=FA%NULOUT,FMT='(A57,I8)')
     S   '!!  ERREUR  !! Valeur incorrecte, non prise en compte: ',KVAL
            GOTO 1001
          ENDIF
C
C   Specification de l'identificateur du centre meteo (defaut=85 pour
C   Toulouse; pour Reading, il vaut 98). Sera utilise pour initialiser
C   KSEC1(2), le 2ieme elt de la section 1 de GRIBEX
C
        ELSEIF (CDCLEF=='IDCEN'.OR.CDCLEF=='idcen') THEN
          IF (KVAL.LT.7.OR.KVAL.GT.99) THEN
            IREP=-125
            WRITE (UNIT=FA%NULOUT,FMT='(A)')  
     S             'Dans FAREGU, action IDCEN:'
            WRITE (UNIT=FA%NULOUT,FMT='(A57,I8)')
     S   '!!  ERREUR  !! Valeur incorrecte, non prise en compte: ',KVAL
            GOTO 1001
          ENDIF
          FA%NSEC1(2,IRANG)=KVAL
C
C   Specification de l'identificateur de modele.
C   FAISC1 initialise automatiquement a
C      177 pour ALADIN
C      211 pour les previsions ARPEGE
C      201 pour les analyses ARPEGE
C   Sera utilise pour initialiser KSEC1(3).
C
        ELSEIF (CDCLEF=='IDMOD'.OR.CDCLEF=='idmod') THEN
          IF (KVAL.LT.0.OR.KVAL.GT.255) THEN
            IREP=-125
            WRITE (UNIT=FA%NULOUT,FMT='(A)')  
     S             'Dans FAREGU, action IDMOD:'
            WRITE (UNIT=FA%NULOUT,FMT='(A57,I8)')
     S   '!!  ERREUR  !! Valeur incorrecte, non prise en compte: ',KVAL
            GOTO 1001
          ENDIF
          FA%NSEC1(3,IRANG)=KVAL
        ELSE
          IREP=-125
          WRITE (UNIT=FA%NULOUT,FMT='(A)') 
     S     '!!  ERREUR  !! Dans FAREGU, action inconnue: '//CDCLEF
          GOTO 1001
        ENDIF
C**
C     2.  -  DEMANDE D'INFORMATION
C-----------------------------------------------------------------------
C
      ELSEIF (KOPT==0) THEN
C
C   Obtention des mots-clef disponibles
C
        IF (CDCLEF=='CLEFS'.OR. CDCLEF=='clefs' .OR.
     S      CDCLEF=='HELP' .OR. CDCLEF=='help') THEN
          KVAL=0
          WRITE (UNIT=FA%NULOUT,FMT=*)
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'Mots clef disponibles pour FAREGU:'
          WRITE (UNIT=FA%NULOUT,FMT=*)
          WRITE (UNIT=FA%NULOUT,FMT=*)'BASIC: pas de compression'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'PACK1: BASIC avec une compression'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       ligne a ligne pour les pts de grille'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'PACK2: BASIC avec une compression avec'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       nb de bits cst pour les groupes'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'PACK3: BASIC avec une compression generale'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       OMM pour les points de grille'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'APAC1: compression agressive:'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       BASIC et PACK1 sont testes'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'APAC2: compression agressive:'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       BASIC, PACK1 et PACK2 sont testes'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'APAC3: compression agressive:'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       BASIC, PACK1 et PACK3 sont testes'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'APAC4: compression agressive:'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       BASIC, PACK1, PACK2 et PACK3 testes'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'WIDPA: lecture/ecriture du nb de bits'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       a utiliser pour les groupes de points'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       de grille dans le cas PACK2'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'GEXTE: la compression avec extensions generales'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       activees (KVAL=1) ou desactivees (KVAL=0)'
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
     S    'IDCEN: lect/ecriture de l''identificateur du'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       centre meteo'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    'IDMOD: lect/ecriture de l''identificateur du'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S    '       modele'
          WRITE (UNIT=FA%NULOUT,FMT=*)
C
C   Lecture du nb de bits a utiliser dans le cadre
C   de la compression avec nb de bits
C   constant par groupe de pts de grille
C   
        ELSEIF (CDCLEF=='WIDPA'.OR.CDCLEF=='widpa') THEN
          KVAL=FA%NCOGRIF(6,IRANG)
C
C   Lecture de la presence ou non de la compression
C   "general extended"
C
        ELSEIF (CDCLEF=='GEXTE'.OR.CDCLEF=='gexte') THEN
          KVAL = FA%NCOGRIF(7,IRANG)/8
C
C   Lecture de la presence ou non du rearrangement
C   boustrophedonique.
C
        ELSEIF (CDCLEF=='BOUST'.OR.CDCLEF=='boust') THEN
          KVAL = FA%NCOGRIF(8,IRANG)/4
C
C   Lecture de la presence ou non de la differentiation spatiale
C
        ELSEIF (CDCLEF=='DIFFE'.OR.CDCLEF=='diffe') THEN
          KVAL=FA%NCOGRIF( 9,IRANG)+FA%NCOGRIF(10,IRANG)
C
C   Lecture de l'identificateur du centre meteo (defaut=85 pour
C   Toulouse; pour Reading, il vaut 98). Sera utilise pour initialiser
C   KSEC1(2), le 2ieme elt de la section 1 de GRIBEX
C
        ELSEIF (CDCLEF=='IDCEN'.OR.CDCLEF=='idcen') THEN
          KVAL=FA%NSEC1(2,IRANG)
C
C   Lecture de l'identificateur du modele
C
        ELSEIF (CDCLEF=='IDMOD'.OR.CDCLEF=='idmod') THEN
          KVAL=FA%NSEC1(3,IRANG)
        ELSE
          IREP=-125
          WRITE (UNIT=FA%NULOUT,FMT='(A)') 
     S      '!!  ERREUR  !! Dans FAREGU, action inconnue: '//CDCLEF
          GOTO 1001
        ENDIF
C**
C     3.  -  OPTION INCONNUE
C-----------------------------------------------------------------------
C
      ELSE
        IREP=-125
        WRITE (UNIT=FA%NULOUT,FMT='(A57,I8)')
     S    '!!  ERREUR  !! Dans FAREGU, option inconnue: KOPT= ',KOPT
        GOTO 1001
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      LLFATA=LLMOER (IREP,IRANG)
C
      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=IXNVMS(IRANG)
      ENDIF
C
      IF (INIMES.EQ.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('FAREGU_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FAREGU'
C
      WRITE (UNIT=CLMESS,FMT='(''IREP='',I4,'', KNUMER='',I3,
     S         '', CDCLEF='''''',A,'''''', KVAL='',I12,
     S        '', KOPT='',I4)')
     S              IREP,KNUMER,CDCLEF,KVAL,KOPT
      CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CLNSPR,.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FAREGU_MT',1,ZHOOK_HANDLE)
      END

