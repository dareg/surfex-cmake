C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAITOU_MT (FA,  KREP, KNUMER, LDNOMM, CDNOMF, CDSTTU, 
     S                      LDERFA, LDIMST, KNIMES, KNBARP, KNBARI, 
     S                      CDNOMC )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme d'OUVERTURE d'une unite logique "Fichier ARPEGE"
C     Il s'agit d'un fichier indexe, traite par le logiciel LFI.
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
      INTEGER INBARP, IRANER, IRANGC, INPAHE, INLATI, ISULEI, INPIND
      INTEGER INPGEO, INIVER, ILONGA, ITRONC, ILACTI, INIMES, INXLON
      INTEGER ITYPTR, IPHASE, IGARDE, IPOSEX, IPUILA
C
      INTEGER IDIMEN (FA%JPCADI), IRDPOL (FA%JPXPAH+FA%JPXIND)
      INTEGER IDATEF (FA%JPLDAT)
      INTEGER (KIND=JPDBLE) ILDIMEN(FA%JPCADI),
     S                      ILRDPOL(FA%JPXPAH+FA%JPXIND)
      INTEGER (KIND=JPDBLE) ILPNVER, ILDATEF(FA%JPLDAT)
C
      REAL (KIND=JPDBLR) ZCHMID (FA%JPCAFS), ZSINLA (FA%JPXGEO)
      REAL (KIND=JPDBLR) ZHYBRI (0:(1+FA%JPXNIV)*2)
C
      LOGICAL LDNOMM, LDERFA, LDIMST, LLVERG, LLNOUF, LLNOUC, LLRLFI
      LOGICAL LLMODC, LLREDF, LLMODA, LLMLAM
C
      CHARACTER CDNOMF*(*), CDSTTU*(*), CDNOMC*(*)
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES DIVERS, ET OUVERTURE DU FICHIER AU SENS "LFI".
C-----------------------------------------------------------------------
C
C     Controle sommaire sur les arguments...le reste est "sous-traite"
C     au sous-programme LFIOUV.
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAITOU_MT',0,ZHOOK_HANDLE)
      IRANG=0
      IRANER=0
      IRANMS=0
      IREPOU=FA%JPNIIL
      LLRLFI=.FALSE.
      LLVERG=.FALSE.
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
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
      LLVERG=FA%LFAMUL
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
C              Ouverture du fichier au sens du logiciel LFI.
C     (on ajoute au nombre d'articles prevus par l'utilisateur les
C      articles constituant le cadre, la date et l'identificateur)
C
      INBARP=KNBARP+7
      CALL LFIOUV_MT (FA%LFI, IREPOU,KNUMER,LDNOMM,CDNOMF,CDSTTU,
     S             LDERFA,LDIMST,
     S             KNIMES,INBARP,KNBARI)
C
      IF (IREPOU.NE.0.AND.IREPOU.NE.-11) THEN
        IREP=IREPOU
        LLRLFI=.TRUE.
        GOTO 1001
      ENDIF
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
          LLMLAM=FA%LIMLAM(IRANGC)
C          
          IDIMEN(1)=FA%MTRONC(IRANGC)
          INLATI=FA%NLATIT(IRANGC)
          IF (.NOT.LLMLAM) THEN
             INPAHE=(1+INLATI)/2
          ELSE
             INPAHE=8
             ISULEI=FA%NOZPAR(1,IRANGC)
C
             INPIND=2*ISULEI+4
             INPGEO=18
          ENDIF                   
          IDIMEN(2)=INLATI
          IDIMEN(3)=FA%NXLOPA(IRANGC)
          INIVER=FA%NNIVER(IRANGC)
          IDIMEN(4)=INIVER
          IDIMEN(5)=FA%NTYPTR(IRANGC)
          ZCHMID(1)=FA%SSLAPO(IRANGC)
          ZCHMID(2)=FA%SCLOPO(IRANGC)
          ZCHMID(3)=FA%SSLOPO(IRANGC)
          ZCHMID(4)=FA%SCODIL(IRANGC)
          ZHYBRI(0)=FA%SPREFE(IRANGC)
          ILNOMC=FA%NLCCAD(IRANGC)
          CLNOMA=CDNOMC
C
          IF (.NOT.LLMLAM) THEN
C
             DO 201 J=1,INPAHE
             IRDPOL(J)=FA%NLOPAR(J,IRANGC)
             IRDPOL(INPAHE+J)=FA%NOZPAR(J,IRANGC)
             ZSINLA(J)=FA%SINLAT(J,IRANGC)
  201        CONTINUE
C
          ELSE
             DO 211 J=1,INPGEO
                ZSINLA(J)=FA%SINLAT(J,IRANGC)
  211        CONTINUE
             DO 213 J=1,INPAHE
                IRDPOL(J)=FA%NLOPAR(J,IRANGC)
  213        CONTINUE
             DO 215 J=1,INPIND
                IRDPOL(INPAHE+J)=FA%NOZPAR(J,IRANGC)
  215        CONTINUE
C
          ENDIF
C
          DO 202 J=0,INIVER
          ZHYBRI(J+1)=FA%SFOHYB(1,J,IRANGC)
          ZHYBRI(J+2+INIVER)=FA%SFOHYB(2,J,IRANGC)
  202     CONTINUE
C
          LLRLFI=.TRUE.
          ILDIMEN=IDIMEN
          CALL LFIECR_MT(FA%LFI,IREP,KNUMER,FA%CPCADI,ILDIMEN,FA%JPCADI)
          IDIMEN=ILDIMEN
          IF (IREP.NE.0) GOTO 1001
C
          CALL LFIECR_MT(FA%LFI,IREP,KNUMER,FA%CPCAFS,ZCHMID,FA%JPCAFS)
          IF (IREP.NE.0) GOTO 1001
C
         IF (.NOT.LLMLAM) THEN
C
             ILONGA=INPAHE*2
             ILRDPOL=IRDPOL
             CALL LFIECR_MT(FA%LFI,IREP,KNUMER,FA%CPCARP,ILRDPOL,ILONGA)
             IRDPOL=ILRDPOL
             IF (IREP.NE.0) GOTO 1001
C
             ILONGA=INPAHE
             CALL LFIECR_MT(FA%LFI,IREP,KNUMER,FA%CPCASL,ZSINLA,ILONGA)
             IF (IREP.NE.0) GOTO 1001
C
          ELSE
C
             ILONGA=INPAHE+INPIND
             ILRDPOL=IRDPOL
             CALL LFIECR_MT(FA%LFI,IREP,KNUMER,FA%CPCARP,ILRDPOL,ILONGA)
             IRDPOL=ILRDPOL
             IF (IREP.NE.0) GOTO 1001
C
             ILONGA=INPGEO
             CALL LFIECR_MT(FA%LFI,IREP,KNUMER,FA%CPCASL,ZSINLA,ILONGA)
             IF (IREP.NE.0) GOTO 1001
C
          ENDIF
C           
          ILONGA=1+(1+INIVER)*2
          CALL LFIECR_MT (FA%LFI, IREP,KNUMER,FA%CPCACH,ZHYBRI,ILONGA)
          IF (IREP.NE.0) GOTO 1001
C
          ILPNVER=FA%JPNVER
          CALL LFIECR_MT (FA%LFI, IREP,KNUMER,
     S                 CLNOMA(1:ILNOMC),ILPNVER,1)
          IF (IREP.NE.0) GOTO 1001
C
          LLRLFI=.FALSE.
          GOTO 300
        ENDIF
C
      ENDIF
C*
C     2.1 - Fichier preexistant...lecture et controle du Cadre "Fichier"
C-----------------------------------------------------------------------
C
      CALL LFINFO_MT (FA%LFI, IREP,KNUMER,FA%CPCADI,ILONGA,IPOSEX)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ELSEIF (ILONGA.EQ.0) THEN
        IREP=-60
        GOTO 1001
      ELSEIF (ILONGA.NE.FA%JPCADI) THEN
        IREP=-61
        GOTO 1001
      ENDIF
C
      ILDIMEN=IDIMEN
      CALL LFILEC_MT (FA%LFI, IREP,KNUMER,FA%CPCADI,ILDIMEN,FA%JPCADI)
      IDIMEN=ILDIMEN
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ENDIF
C
      CALL LFINFO_MT (FA%LFI, IREP,KNUMER,FA%CPCAFS,ILONGA,IPOSEX)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ELSEIF (ILONGA.EQ.0) THEN
        IREP=-60
        GOTO 1001
      ELSEIF (ILONGA.NE.FA%JPCAFS) THEN
        IREP=-61
        GOTO 1001
      ENDIF
C
      CALL LFILEC_MT (FA%LFI, IREP,KNUMER,FA%CPCAFS,ZCHMID,FA%JPCAFS)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ENDIF
C
C        Coherence des dimensions par rapport aux valeurs "licites",
C     que l'on doit faire avant de poursuivre les lectures.
C
      IF(IDIMEN(5).LE.0) LLMLAM = .TRUE.
      ITRONC=IDIMEN(1)
      INLATI=IDIMEN(2)
      INPAHE=(1+INLATI)/2
      INXLON=IDIMEN(3)
      INIVER=IDIMEN(4)
      ITYPTR=IDIMEN(5)
      IPHASE=1
      IGARDE=1
      CALL FACADI_MT(FA, IREP,CDNOMC,ITYPTR,ZCHMID(1),ZCHMID(2),
     S               ZCHMID(3),ZCHMID(4),ITRONC,INLATI,INXLON,IRDPOL(1),
     S               IRDPOL(FA%JPXPAH+1),ZSINLA,
     S               INIVER,ZHYBRI(0),ZHYBRI(1),ZHYBRI(FA%JPXNIV+2),
     S               LLMODC,LLREDF,IPHASE,IRANGC,ILNOMC,IGARDE)
      IF (IREP.NE.0) GOTO 1001
C
      CALL LFINFO_MT (FA%LFI, IREP,KNUMER,FA%CPCARP,ILONGA,IPOSEX)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ELSEIF (ILONGA.EQ.0) THEN
        IREP=-60
        GOTO 1001
      ELSEIF (ILONGA.NE.INPAHE*2) THEN
        IF (.NOT.LLMLAM) THEN
           IREP=-61
           GOTO 1001
        ENDIF    
      ENDIF
C
      ILRDPOL=IRDPOL
      CALL LFILEC_MT (FA%LFI, IREP,KNUMER,FA%CPCARP,ILRDPOL,ILONGA)
      IRDPOL=ILRDPOL
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ENDIF
C
      CALL LFINFO_MT (FA%LFI, IREP,KNUMER,FA%CPCASL,ILONGA,IPOSEX)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ELSEIF (ILONGA.EQ.0) THEN
        IREP=-60
        GOTO 1001
      ELSEIF (ILONGA.NE.INPAHE) THEN
        IF (.NOT.LLMLAM) THEN
           IREP=-61
           GOTO 1001
        ENDIF      
      ENDIF
C
      CALL LFILEC_MT (FA%LFI, IREP,KNUMER,FA%CPCASL,ZSINLA,ILONGA)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ENDIF
C
      CALL LFINFO_MT (FA%LFI, IREP,KNUMER,FA%CPCACH,ILONGA,IPOSEX)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ELSEIF (ILONGA.EQ.0) THEN
        IREP=-60
        GOTO 1001
      ELSEIF (ILONGA.NE.1+(1+INIVER)*2) THEN
        IF (.NOT.LLMLAM) THEN
           IREP=-61
           GOTO 1001
        ENDIF      
      ENDIF
C
      CALL LFILEC_MT (FA%LFI, IREP,KNUMER,FA%CPCACH,ZHYBRI,ILONGA)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ENDIF
C
C        Tests complementaires sur les valeurs lues.
C
      IPHASE=2
      CALL FACADI_MT(FA, IREP,CDNOMC,ITYPTR,ZCHMID(1),ZCHMID(2),
     S               ZCHMID(3),ZCHMID(4),ITRONC,INLATI,INXLON,IRDPOL(1),
     S               IRDPOL(INPAHE+1),ZSINLA,
     S               INIVER,ZHYBRI(0),ZHYBRI(1),ZHYBRI(INIVER+2),
     S               LLMODC,LLREDF,IPHASE,IRANGC,ILNOMC,IGARDE)
      IF (IREP.NE.0) GOTO 1001
C*
C     2.2 - Fichier preexistant...l'identificateur du fichier est le
C           premier article suivant les articles du cadre.
C-----------------------------------------------------------------------
C
      CALL LFICAS_MT (FA%LFI, IREP,KNUMER,CLNOMA,ILONGA,
     S                IPOSEX,.FALSE.)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ELSEIF (ILONGA.EQ.0) THEN
        IREP=-110
        GOTO 1001
      ENDIF
C
C*
C     2.3 - Fichier preexistant...lecture et controle de l'article DATE.
C-----------------------------------------------------------------------
C
      CALL LFINFO_MT (FA%LFI, IREP,KNUMER,FA%CPDATE,ILONGA,IPOSEX)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ELSEIF (ILONGA.EQ.0) THEN
        IREP=-62
        GOTO 1001
      ELSEIF (ILONGA.NE.FA%JPLDAT) THEN
        IREP=-63
        GOTO 1001
      ENDIF
C
      ILDATEF=IDATEF
      CALL LFILEC_MT (FA%LFI, IREP,KNUMER,FA%CPDATE,ILDATEF,FA%JPLDAT)
      IDATEF=ILDATEF
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ENDIF
C
C        La ligne ci-dessous evite a FANDAI de croire, eventuellement,
C     a une redefinition de date.
C
      FA%LCREAF(IRANG)=.TRUE.
C
C        Controle de la Date fichier, et stockage dans FA%MADATE.
C
      CALL FANDAI_MT (FA, IREP,IRANG,IDATEF,LLMODA)
      IF (IREP.NE.0) GOTO 1001
C
C         Definition du Cadre proprement dite.
C
      IPHASE=3
      CALL FACADI_MT(FA, IREP,CDNOMC,ITYPTR,ZCHMID(1),ZCHMID(2),
     S               ZCHMID(3),ZCHMID(4),ITRONC,INLATI,INXLON,IRDPOL(1),
     S               IRDPOL(INPAHE+1),ZSINLA,
     S               INIVER,ZHYBRI(0),ZHYBRI(1),ZHYBRI(INIVER+2),
     S               LLMODC,LLREDF,IPHASE,IRANGC,ILNOMC,IGARDE)
      IF (IREP.NE.0) GOTO 1001
C**
C     3.  -  ON MET A JOUR LES TABLES RELATIVES AUX FICHIERS.
C-----------------------------------------------------------------------
C
  300 CONTINUE
C
      FA%NFIOUV=FA%NFIOUV+1
      FA%NULIND(FA%NFIOUV)=IRANG
      FA%NULOGI(IRANG)=KNUMER
      FA%NUCADR(IRANG)=IRANGC
C
      FA%LNOMME(IRANG)=LDNOMM
      FA%NIVOMS(IRANG)=KNIMES
      FA%LERRFA(IRANG)=LDERFA
      FA%LCREAF(IRANG)=LLNOUF
      FA%NBFPDG(IRANG)=FA%NBIPDG
      FA%NBFCSP(IRANG)=FA%NBICSP
      FA%NPUFLA(IRANG)=FA%NPUILA
      FA%NMFDPL(IRANG)=FA%NMIDPL
      FA%NFGRIB(IRANG)=FA%NIGRIB
      FA%CIDENT(IRANG)=CLNOMA
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
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'ASGN')
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
C           Deverrouillage global eventuel.
C
      IF (LLVERG) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'OFF')
C
      IF (.NOT.LLFATA.AND.INIMES.EQ.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('FAITOU_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FAITOU'
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
     S               CLNSPR,CLACTI(1:ILACTI),LLRLFI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAITOU_MT',1,ZHOOK_HANDLE)
      END

