C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACADE_MT (FA,  CDNOMC, KTYPTR, PSLAPO, PCLOPO, PSLOPO,
     S                    PCODIL, KTRONC, KNLATI, KNXLON, KNLOPA,
     S                    KNOZPA, PSINLA, KNIVER, PREFER, PAHYBR,
     S                    PBHYBR, LDGARD )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme servant a DEfinir un CADre, voire a le
C     redefinir.
C**
C        Arguments : CDNOMC ==> Nom symbolique du cadre;
C  (tous d'Entree)   KTYPTR ==> Type de transformation horizontale;
C                    PSLAPO ==> Sinus de la latitude du pole d'interet;
C                    PCLOPO ==> Cosinus " " longitude "   "       "   ;
C                    PSLOPO ==> Sinus   " " longitude "   "       "   ;
C                    PCODIL ==> Coefficient de dilatation;
C                    KTRONC ==> Troncature;
C                    KNLATI ==> Nombre de latitudes (de pole a pole);
C                    KNXLON ==> Nombre maxi de longitudes par parallele;
C         (Tableau)  KNLOPA ==> Nombre de longitudes par parallele;
C                               (du pole nord vers l'equateur seulement)
C         (Tableau)  KNOZPA ==> Nombre d'onde zonal maxi par parallele;
C                               (du pole nord vers l'equateur seulement)
C         (Tableau)  PSINLA ==> Sinus des latitudes de l'hemisphere nord
C                               (du pole nord vers l'equateur seulement)
C                    KNIVER ==> Nombre de niveaux verticaux;
C                    PREFER ==> Pression de reference (facteur multipli-
C                               catif de la premiere fonction de la
C                               coordonnee hybride)
C         (Tableau)  PAHYBR ==> Valeurs de la fonction "A" de la coordo-
C                               nnee hybride AUX LIMITES DE COUCHES;
C         (Tableau)  PBHYBR ==> Valeurs de la fonction "B" de la coordo-
C                               nnee hybride AUX LIMITES DE COUCHES;
C                    LDGARD ==> Vrai si le cadre doit etre conserve meme
C                               apres la fermeture du dernier fichier
C                               qui s'y rattache.
C*
C        La "redefinition" d'un cadre est possible a l'une de ces
C     conditions:
C
C     - le cadre a ete defini, mais n'a aucun fichier qui s'y rattache;
C     - le cadre defini a au moins un fichier qui s'y rattache, et les
C       nouveaux parametres de definition sont identiques a ceux deja
C       definis.
C
C        Toute "redefinition" de cadre donne lieu a une messagerie
C     de niveau 1, donc non masquee par defaut.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KTYPTR, KTRONC, KNLATI, KNXLON, KNIVER
C
      INTEGER KNLOPA ((1+KNLATI)/2), KNOZPA ((1+KNLATI)/2)
C
      REAL (KIND=JPDBLR) PSLAPO, PCLOPO, PSLOPO, PCODIL, PREFER
C
      REAL (KIND=JPDBLR) PSINLA ((1+KNLATI)/2), PAHYBR (0:KNIVER)
      REAL (KIND=JPDBLR) PBHYBR (0:KNIVER)
C
      CHARACTER CDNOMC*(*)
C
      LOGICAL LDGARD
C
      INTEGER IPHASE, IGARDE, IREP, IRANGC, ILNOMC, INIMES, INUMER
C
      LOGICAL LLREDF, LLMODC
C
C
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  SI PREMIERE UTILISATION, APPEL AU SOUS-PROGRAMME "FARINE".
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACADE_MT',0,ZHOOK_HANDLE)
      IF (FA%FACADE_LLPREA) THEN
        CALL FARINE_MT (FA, 2)
        FA%FACADE_LLPREA=.FALSE.
      ENDIF
C**
C     2.  -  LE TRAVAIL EST SOUS-TRAITE AU SOUS-PROGRAMME "FACADI".
C-----------------------------------------------------------------------
C
      IPHASE=0
C
      IF (LDGARD) THEN
        IGARDE=2
      ELSE
        IGARDE=0
      ENDIF
C
C             Verrouillage global prealable, si necessaire.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
      CALL FACADI_MT (FA, IREP,CDNOMC,KTYPTR,PSLAPO,PCLOPO,
     S             PSLOPO,PCODIL,
     S             KTRONC,KNLATI,KNXLON,KNLOPA,KNOZPA,PSINLA,KNIVER,
     S             PREFER,PAHYBR,PBHYBR,LLMODC,LLREDF,IPHASE,IRANGC,
     S             ILNOMC,IGARDE)
      ILNOMC=MIN (ILNOMC,FA%NCPCAD)
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE sous-programme "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
C
C          Deverrouillage global eventuel.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'OFF')
C
      LLFATA=LLMOER(IREP,0)
C
      IF (LLFATA) THEN
        INIMES=2
      ELSEIF (FA%NIMSGA.EQ.0) THEN
        INIMES=0
      ELSEIF (LLMODC) THEN
        INIMES=1
        WRITE (UNIT=CLMESS,FMT=
     S  '(''PARAMETRES NUMERIQUES DU CADRE '''''',A,'''''' MODIFIES '',
     S    '' - CONSERVATION A LA FERMETURE DU DERNIER FICHIER= '',L1)')
     S   CDNOMC(1:ILNOMC),LDGARD
      ELSEIF (LLREDF) THEN
        INIMES=1
        WRITE (UNIT=CLMESS,FMT='(''CADRE '''''',A,
     S '''''' REDEFINI - MEMES PARAMETRES NUMERIQUES - '',
     S '' CONSERVATION A LA FERMETURE DU DERNIER FICHIER= '',L1)')
     S   CDNOMC(1:ILNOMC),LDGARD
      ELSEIF (FA%NIMSGA.EQ.2) THEN
        INIMES=2
      ELSE
        INIMES=0
      ENDIF
C
      IF (INIMES.EQ.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('FACADE_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FACADE'
      INUMER=FA%JPNIIL
C
      IF (INIMES.EQ.1.AND.FA%NIMSGA.EQ.2) THEN
C
C        Cas ou il faut en fait 2 messages.
C
        CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI,.FALSE.)
        INIMES=2
      ENDIF
C
      IF (INIMES.EQ.2) THEN
C
        IF (IREP.EQ.-65.AND.ILNOMC.EQ.1) THEN
          ILNOMC=8
          CLACTI(1:ILNOMC)=FA%CHAINC(:ILNOMC)
        ELSE
          ILNOMC=MIN (LEN (CLACTI),ILNOMC,FA%NCPCAD)
          CLACTI(1:ILNOMC)=CDNOMC(1:ILNOMC)
        ENDIF
C
        WRITE (UNIT=CLMESS,
     S         FMT='(''ARGUMENTS SIMPLES= '''''',A,'''''',''
     S  ,I2,4('','',F7.4),3('','',I6),'','',I5,'','',F11.4,'', '',L1)')
     S       CLACTI(1:ILNOMC),KTYPTR,PSLAPO,PCLOPO,PSLOPO,PCODIL,
     S       KTRONC,KNLATI,KNXLON,KNIVER,PREFER,LDGARD
      ENDIF
C
      CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CLACTI(1:ILNOMC),.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FACADE_MT',1,ZHOOK_HANDLE)
      END

