C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FANOUV_MT (FA,  KREP, KNUMER, LDNOMM, CDNOMF, CDSTTU, 
     S                      LDERFA,
     S                     LDIMST, KNIMES, KNBARP, KNBARI, CDNOMC )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme d'initialisation SANS OUVERTURE d'une unite logique 
C        "Fichier ARPEGE". Il s'agit d'un fichier indexe, 
C        traite par le logiciel LFI.
C
C      FANOUV est derive de FAITOU, mais ne fait pas l'appel a la couche LFI
C             pour l'ouverture reelle.
C
C      utilise pour la seule compression des donnees par des processeurs
C      qui ne font pas d'ecriture effective sur disque
C
C**
C     ARGUMENTS : Ce sont les memes que pour "LFIOUV", avec CDNOMC comme
C                 argument supplementaire.
C
C                 KREP   (Sortie) ==> Code-reponse du sous-programme;
C                 KNUMER (Entree) ==> Numero de l'unite logique;
C                 LDNOMM (Entree) ==> Vrai si l'unite logique doit etre
C                                     associee a un NOM de Fichier EXP-
C                                     LICITE lors de l'"OPEN" FORTRAN;
C                 CDNOMF (Entree) ==> Nom de fichier explicite, si
C                                     *LDNOMM* est VRAI - Meme si ce
C                                     n'est pas le cas, ce *DOIT* ETRE
C                                     UN OBJET DE TYPE "CHARACTER" .
C                 CDSTTU (Entree) ==> "STATUS" pour l'"OPEN" FORTRAN
C                                     ('OLD','NEW','UNKNOWN','SCRATCH')
C                                     par defaut, mettre 'UNKNOWN';
C                 LDERFA (Entree) ==> Option d'erreur fatale;
C                 LDIMST (Entree) ==> Option impression de Statistiques
C                                     au moment de la fermeture;
C                 KNIMES (Entree) ==> Niveau de la Messagerie (0,1 ou 2)
C                                     ( 0==>Rien, 2==>Tout )
C                 KNBARP (Entree) ==> Nombre d'articles logiques prevus,
C                                     ce qui n'est utilise que lors de
C                                     la Creation du fichier,
C                                     et qui n'empeche quand meme pas
C                                     d'avoir plus d'articles logiques;
C                 KNBARI (Sortie) ==> Nombre d'articles logiques de don-
C                                     nees sur le fichier, initialement.
C                                     (zero si creation)
C                 CDNOMC (Entree) ==> Nom du CADRE associe au fichier.
C*
C     N.B. :  Pour un fichier en mode creation, ce cadre doit avoir ete
C          defini au prealable (via le sous-programme FACADE, ou par
C          l'ouverture d'un fichier preexistant).
C             Pour un fichier ARPEGE preexistant, le cadre est lu sur le
C          fichier; s'il etait deja defini auparavant, il y a controle
C          de coherence entre les deux versions du cadre.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER, KNIMES, KNBARP, KNBARI
C
      CHARACTER CPNOMD*(*)
      PARAMETER ( CPNOMD='%%%%% FICHIER SANS NOM %%%%%' )
C
      INTEGER IRANG, INUMER, IRANMS, IREPOU, ILNOMC, ILOMIN, IREP, J
      INTEGER INBARP, IRANER, IRANGC
      INTEGER ILONGA, ITRONC, ILACTI, INIMES, INXLON
      INTEGER ITYPTR, IPHASE, IGARDE, IPOSEX, IPUILA
C
      INTEGER IDATEF (FA%JPLDAT)
      INTEGER (KIND=JPDBLE) ILDIMEN(FA%JPCADI),
     S                      ILRDPOL(FA%JPXPAH+FA%JPXIND)
      INTEGER (KIND=JPDBLE) ILPNVER, ILDATEF(FA%JPLDAT)
C
      LOGICAL LDNOMM, LDERFA, LDIMST, LLNOUF, LLNOUC, LLRLFI
      LOGICAL LLMODC, LLREDF, LLMODA
C
      CHARACTER CDNOMF*(*), CDSTTU*(*), CDNOMC*(*)
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES DIVERS
C-----------------------------------------------------------------------
C
C     Controle sommaire sur les arguments...
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FANOUV_MT',0,ZHOOK_HANDLE)
      IRANG=0
      IRANER=0
      IRANMS=0
      IREPOU=FA%JPNIIL
      LLRLFI=.FALSE.
      ILNOMC=LEN (CDNOMC)
      ILOMIN=MIN ( LEN (CDNOMF), LEN (CDSTTU), ILNOMC)
C
C        L'appel ci-dessous est legerement anticipe, de maniere a
C     initialiser les variables globales du logiciel s'il s'agit
C     du premier appel a un sous-programme de ce logiciel.
C
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (ILOMIN.LE.0) THEN
        IREP=-65
        GOTO 1001
      ELSEIF (IRANG.NE.0) THEN
C
C            Controle de non-ouverture prealable (au sens du logiciel)
C
        IREP=-55
        IRANMS=IRANG
        GOTO 1001
      ENDIF
C
C             Verrouillage global, si necessaire.
C
C        A-t-on deja atteint le nombre limite de fichiers ARPEGE
C     ouverts simultanement ? Si non, on cherche un emplacement libre
C     dans la table FA%NULOGI (logiquement, il devrait en exister un)
C
      IF (FA%NFIOUV.GE.FA%JPNXFA) THEN
        IREP=-56
        GOTO 1001
      ELSE
C
        DO 101 J=1,FA%JPNXFA
C
        IF (FA%NULOGI(J).EQ.FA%JPNIIL) THEN
          IRANG=J
          GOTO 102
        ENDIF
C
  101   CONTINUE
C
        IREP=-66
        GOTO 1001
C
  102   CONTINUE
C
      ENDIF
C
C**
C     2.  -  CONTROLES SPECIFIQUES AU LOGICIEL DE FICHIERS ARPEGE.
C-----------------------------------------------------------------------
C
      LLNOUF=KNBARI.EQ.0
      CALL FANUCA_MT (FA, CDNOMC,IRANGC,.FALSE.)
      LLNOUC=IRANGC.EQ.0
C
      IF (LLNOUF) THEN
C
        IF (LLNOUC) THEN
          IREP=-57
          GOTO 1001
        ELSE
C
C         Fichier en mode creation et cadre predefini... OK a ce niveau.
C
C       On ecrit les articles definissant le cadre sur le fichier,
C     ainsi qu'un article ayant pour nom l'identificateur "par defaut",
C     (en fait, le nom du cadre) de maniere a ce que cet article soit
C     sequentiellement celui qui suit le dernier article du cadre.
C
          ILNOMC=FA%NLCCAD(IRANGC)
C
        ENDIF
C
      ENDIF
C*
C        Controle de la Date fichier, et stockage dans FA%MADATE.
C
!  IDATEF arbitraire pour contenter FACOND
      IDATEF(1) = 1993
      IDATEF(2) = 9
      IDATEF(3) = 2
      IDATEF(4) = 0
      IDATEF(5) = 0
      IDATEF(6) = 1
      IDATEF(7) = 0
      IDATEF(8) = 0
      IDATEF(9) = 1
      IDATEF(10) = 0
      IDATEF(11) = 0

      DO J=1,FA%JPLDAT
        FA%MADATE(J,IRANG)=IDATEF(J)
      ENDDO

C**
C     3.  -  ON MET A JOUR LES TABLES RELATIVES AUX FICHIERS.
C-----------------------------------------------------------------------
C
  300 CONTINUE
!
      IREPOU=0
!
C
      FA%NFIOUV=FA%NFIOUV+1
      FA%NULIND(FA%NFIOUV)=IRANG
      FA%NULOGI(IRANG)=KNUMER
      FA%NUCADR(IRANG)=IRANGC
C
      FA%LNOMME(IRANG)=LDNOMM
      FA%NIVOMS(IRANG)=KNIMES
      FA%LERRFA(IRANG)=LDERFA
      FA%LCREAF(IRANG)=.FALSE.
      FA%NBFPDG(IRANG)=FA%NBIPDG
      FA%NBFCSP(IRANG)=FA%NBICSP
      FA%NPUFLA(IRANG)=FA%NPUILA
      FA%NMFDPL(IRANG)=FA%NMIDPL
      FA%NFGRIB(IRANG)=FA%NIGRIB
      FA%CIDENT(IRANG)=CDNOMC
      ITRONC=FA%MTRONC(IRANGC)
      ITYPTR=FA%NTYPTR(IRANGC)
C
      IF (ITYPTR.LT.0) THEN
        FA%NSTROF(IRANG)=MIN (FA%NSTROI,ITRONC-1,-ITYPTR-1)
      ELSE
        FA%NSTROF(IRANG)=MIN (FA%NSTROI,ITRONC-1)
      ENDIF
C
C Appel a FAINOC pour interpreter les eventuels defauts
C de -1 pris par FA%NBFPDG, FA%NBFCSP, FA%NSTROF et FA%NPUFLA en
C IRANG-ieme position.
C
      CALL FAINOC_MT (FA,  IRANG )
C
      IRANER=IRANG
      IRANMS=IRANG
      IPUILA=FA%NPUFLA(IRANG)
C
      FA%NCOGRIF(:,IRANG)=FA%NCODGRI(:)
      FA%NRASHO(IRANG) = 0
      FA%NRASVE(IRANG) = 0
C
C   L'initialisation de FLAP1Dx sera faite dans FACSIM
C
      FA%LIFLAP(IRANG)=.TRUE.
C
C        On incremente le nombre de fichiers attaches au cadre specifie.
C
      FA%NULCAD(IRANGC)=FA%NULCAD(IRANGC)+1
      IREP=IREPOU
      GOTO 1001
C**
C     9.  - CI-DESSOUS, ETIQUETTES DE BRANCHEMENT EN CAS D'ERREUR E/S.
C-----------------------------------------------------------------------
C
  901 CONTINUE
      CLACTI='INQUIRE'
C
C      AU CAS OU, ON FORCE LE CODE-REPONSE ENTREE/SORTIE A ETRE POSITIF.
C
      IREP=IABS (IREP)
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "LFIEMS" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANER)
C
      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=IXNVMS (IRANMS)
      ENDIF
C
      IF (.NOT.LLFATA.AND.INIMES.EQ.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('FANOUV_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FANOUV'
C
      IF (INIMES.EQ.2) THEN
C
        IF (ILNOMC.GT.0) THEN
          ILACTI=MIN (LEN (CLACTI),ILNOMC)
          CLACTI(1:ILACTI)=CDNOMC(1:ILNOMC)
        ELSE
          ILACTI=8
          CLACTI=FA%CHAINC(:ILACTI)
        ENDIF
C
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S     '', LDNOMM= '',L1,'', CDSTTU='''''',A7,'''''', LDERFA= '',L1,
     S     '',  LDIMST= '',L1,
     S         '', KNIMES='',I2,'', KNBARP='',I6,'' KNBARI='',I6)')
     S   KREP,KNUMER,LDNOMM,CDSTTU,LDERFA,LDIMST,KNIMES,KNBARP,KNBARI
        CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI(1:ILACTI),LLRLFI)
        CLMESS='CDNOMC='''//CLACTI(1:ILACTI)//''''
        CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S               CLNSPR,
     S               CLACTI(1:ILACTI),LLRLFI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FANOUV_MT',1,ZHOOK_HANDLE)
      END

