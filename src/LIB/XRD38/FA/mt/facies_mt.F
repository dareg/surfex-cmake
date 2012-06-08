C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACIES_MT (FA,  CDNOMC, KTYPTR, PSLAPO, PCLOPO, PSLOPO,
     S                    PCODIL, KTRONC, KNLATI, KNXLON, KNLOPA,
     S                    KNOZPA, PSINLA, KNIVER, PREFER, PAHYBR,
     S                    PBHYBR, LDGARD )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme servant a obtenir le contenu d'un Cadre.
C     ( FACIES... par analogie avec FADIES, avec "C" pour Cadre... )
C**
C        Arguments : CDNOMC ==> Nom symbolique du cadre;
C  (tous de Sortie,  KTYPTR ==> Type de transformation horizontale
C   sauf CDNOMC)     PSLAPO ==> Sinus de la latitude du pole d'interet;
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
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KTYPTR, KTRONC, KNLATI, KNXLON, KNIVER
C
      INTEGER KNLOPA (FA%JPXPAH), KNOZPA (FA%JPXIND)
C
      REAL (KIND=JPDBLR) PSLAPO, PCLOPO, PSLOPO, PCODIL, PREFER
      REAL (KIND=JPDBLR) PSINLA (FA%JPXGEO), PAHYBR (0:FA%JPXNIV)
      REAL (KIND=JPDBLR) PBHYBR (0:FA%JPXNIV)
C
      CHARACTER CDNOMC*(*)
C
      LOGICAL LDGARD
C
      INTEGER IPHASE, IREP, IRANGC, ILNOMC, INIMES, INUMER, ILCDNO, J
      INTEGER INPAHE, ISULEI, INPIND, INPGEO
C
      LOGICAL LLVERG, LLMLAM
#include "facom_mt.h"
C
C
C**
C     0.  -  SI PREMIERE UTILISATION, APPEL AU SOUS-PROGRAMME "FARINE".
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACIES_MT',0,ZHOOK_HANDLE)
      IF (FA%FACIES_LLPREA) THEN
        CALL FARINE_MT (FA, 2)
        FA%FACIES_LLPREA=.FALSE.
      ENDIF
C**
C     1.  -  CONTROLE DE L'ARGUMENT "CDNOMC".
C-----------------------------------------------------------------------
C
      LLVERG=.FALSE.
      ILCDNO=LEN (CDNOMC)
      ILNOMC=1
C
      IF (ILCDNO.LE.0) THEN
        IREP=-65
        GOTO 1001
      ELSEIF (CDNOMC.EQ.' ') THEN
        IREP=-68
        GOTO 1001
      ENDIF
C
      DO 101 J=ILCDNO,1,-1
C
      IF (CDNOMC(J:J).NE.' ') THEN
        ILNOMC=J
        GOTO 102
      ENDIF
C
  101 CONTINUE
C
  102 CONTINUE
C
      IF (ILNOMC.GT.FA%NCPCAD) THEN
        IREP=-65
        GOTO 1001
      ENDIF
C**
C     2.  -  RECHERCHE DU CADRE DANS LES TABLES.
C-----------------------------------------------------------------------
C
C             Verrouillage global prealable, si necessaire.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
      LLVERG=FA%LFAMUL
C
      CALL FANUCA_MT (FA, CDNOMC,IRANGC,.FALSE.)
C
      IF (IRANGC.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ENDIF
C**
C     3.  -  TRANSFERT DES TABLES DU LOGICIEL DANS LES ARGUMENTS.
C-----------------------------------------------------------------------
C
      LLMLAM=FA%LIMLAM(IRANGC)
C
      KTYPTR=FA%NTYPTR(IRANGC)
      KTRONC=FA%MTRONC(IRANGC)
      KNIVER=FA%NNIVER(IRANGC)
      KNLATI=FA%NLATIT(IRANGC)
C
      IF (.NOT.LLMLAM) THEN
         INPAHE=(1+KNLATI)/2
      ELSE
         ISULEI=FA%NOZPAR(1,IRANGC)
         INPAHE=8
         INPIND=2*ISULEI+4
         INPGEO=18
      ENDIF
C
      KNXLON=FA%NXLOPA(IRANGC)
      PSLAPO=FA%SSLAPO(IRANGC)
      PCLOPO=FA%SCLOPO(IRANGC)
      PSLOPO=FA%SSLOPO(IRANGC)
      PCODIL=FA%SCODIL(IRANGC)
      PREFER=FA%SPREFE(IRANGC)
      LDGARD=FA%NGARDE(IRANGC).EQ.2.OR.
     S       (FA%NGARDE(IRANGC).EQ.1.AND.FA%LIGARD)
C
      IF (.NOT.LLMLAM) THEN
C
         DO 301 J=1,INPAHE
         KNLOPA(J)=FA%NLOPAR(J,IRANGC)
         KNOZPA(J)=FA%NOZPAR(J,IRANGC)
         PSINLA(J)=FA%SINLAT(J,IRANGC)
  301    CONTINUE
C
      ELSE
C
         DO 311 J=1,INPAHE
            KNLOPA(J)=FA%NLOPAR(J,IRANGC)
  311    CONTINUE
         DO 313 J=1,INPIND
            KNOZPA(J)=FA%NOZPAR(J,IRANGC)
  313    CONTINUE
         DO 315 J=1,INPGEO
            PSINLA(J)=FA%SINLAT(J,IRANGC)
  315    CONTINUE
C
      ENDIF
C
      DO 302 J=0,KNIVER
      PAHYBR(J)=FA%SFOHYB(1,J,IRANGC)
      PBHYBR(J)=FA%SFOHYB(2,J,IRANGC)
  302 CONTINUE
C
      IREP=0
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE sous-programme "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
C
C          Deverrouillage global eventuel.
C
      IF (LLVERG) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'OFF')
C
      LLFATA=IREP.NE.0.AND.FA%NRFAGA.NE.2
C
      IF (LLFATA.OR.FA%NIMSGA.EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('FACIES_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FACIES'
C
      IF (IREP.EQ.-65.AND.ILCDNO.LE.0) THEN
        ILNOMC=8
        CLACTI(1:ILNOMC)=FA%CHAINC(:ILNOMC)
      ELSE
        ILNOMC=MIN (LEN (CLACTI),ILNOMC)
        CLACTI(1:ILNOMC)=CDNOMC(1:ILNOMC)
      ENDIF
C
      ILNOMC=MIN (ILNOMC,FA%NCPCAD)
      WRITE (UNIT=CLMESS,
     S       FMT='(''ARGUMENTS SIMPLES= '''''',A,'''''','',
     S   I2,4('','',F7.4),3('','',I6),'','',I5,'','',F11.4,'', '',L1)')
     S     CLACTI(1:ILNOMC),KTYPTR,PSLAPO,PCLOPO,PSLOPO,PCODIL,
     S     KTRONC,KNLATI,KNXLON,KNIVER,PREFER,LDGARD
      INUMER=FA%JPNIIL
      CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CLACTI(1:ILNOMC),.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FACIES_MT',1,ZHOOK_HANDLE)
      END

