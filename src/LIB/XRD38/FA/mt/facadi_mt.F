C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACADI_MT (FA,  KREP, CDNOMC, KTYPTR, PSLAPO, PCLOPO, 
     S                    PSLOPO,
     S                    PCODIL, KTRONC, KNLATI, KNXLON, KNLOPA,
     S                    KNOZPA, PSINLA, KNIVER, PREFER, PAHYBR,
     S                    PBHYBR, LDMODC, LDREDF, KPHASE, KRANGC,
     S                    KLNOMC, KGARDE )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme A USAGE INTERNE AU LOGICIEL. Fait la plupart
C     des controles en vue de Definir un CADre, voire le redefinir.
C        En mode multi-taches, il doit y avoir verrouillage global
C     de la zone d'appel au sous-programme.
C**
C        Arguments : KREP   ==> Code-reponse du sous-programme;
C                    CDNOMC ==> Nom symbolique du cadre;
C  (tous d'Entree)   KTYPTR ==> Type de transformation horizontale;
C   sauf KRANGC      PSLAPO ==> Sinus de la latitude du pole d'interet;
C     et KLNOMC)     PCLOPO ==> Cosinus " " longitude "   "       "   ;
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
C                    LDMODC ==> Vrai s'il y a modification d'un cadre
C                               deja defini au prealable;
C                    LDREDF ==> Vrai s'il y a redefinition d'un cadre
C                               au sens large du terme (avec ou sans
C                               modification).
C                    KPHASE ==> Indique quelle(s) phase(s) du sous-prog.
C                               on doit executer:
C                               0 ==> Toutes,
C                               1 ==> Controle des variables simples,
C                               2 ==> Controle des tableaux,
C                               3 ==> Definition du cadre seule.
C        (Sortie)    KRANGC ==> Rang du cadre dans les tables.
C (Sortie si phase 1,KLNOMC ==> Longueur en caracteres du nom de cadre.
C  Entree sinon ! )
C                    KGARDE ==> Option de conservation du cadre
C                               apres la fermeture du dernier fichier
C                               qui s'y rattache. A noter que lors dans
C                               le cas d'une definition dynamique de
C                               cadre (appel par FAITOU, avec KGARDE=1),
C                               une redefinition de cadre n'est toleree
C                               qu'a l'identique.
C
C     N.B. :    En mode multi-taches, si l'on appelle le sous-programme
C            avec KPHASE=0 ou KPHASE=3, on doit verrouiller dans le
C            programme appelant l'appel au sous-programme.
C               Par ailleurs, LDMODC et LDREDF ne sont definis que si
C            KPHASE=0 ou KPHASE=3.
C*
C        La "redefinition" d'un cadre est possible a l'une de ces
C     conditions:
C
C     - le cadre a ete defini, mais n'a aucun fichier qui s'y rattache;
C     - le cadre defini a au moins un fichier qui s'y rattache, et les
C       nouveaux parametres de definition sont identiques a ceux deja
C       definis (a l'exception de l'option de conservation).
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KTYPTR, KTRONC, KNLATI, KNXLON, KNIVER, KREP, KPHASE
      INTEGER KRANGC, KLNOMC, KGARDE
C
      INTEGER KNLOPA (FA%JPXPAH), KNOZPA (FA%JPXIND)
C
      REAL (KIND=JPDBLR) PSLAPO, PCLOPO, PSLOPO, PCODIL, PREFER
C
      REAL (KIND=JPDBLR) PSINLA (FA%JPXGEO), PAHYBR (0:KNIVER)
      REAL (KIND=JPDBLR) PBHYBR (0:KNIVER)
      REAL (KIND=JPDBLR),PARAMETER ::  ZEPS=1.E-15_JPDBLR
C
      CHARACTER CDNOMC*(*)
C
      LOGICAL LDREDF, LDMODC
C
      INTEGER INEXPL, INGEOM, INPAHE, ILCDNO, J, IPREC, ICOMPT, IMSMAX
      INTEGER ISFLAM, JL, IK, INIMES, INUMER, ILNOMC
C
      INTEGER IESN0 (0:FA%JPXTRO), ISUHIG (0:FA%JPXTRO)
      INTEGER IKNTMP(0:FA%JPXTRO), IKMTMP(0:FA%JPXTRO), 
     S        ICPL4N(0:FA%JPXTRO)
C
      REAL (KIND=JPDBLR) ZMIN, ZPMIN, ZPMAX, ZPMINP, ZPMAXP, ZMSMAX
C
      LOGICAL LLMLAM
#include "facom_mt.h"
C**
C     0.  -  AIGUILLAGE EN FONCTION DE *KPHASE*.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACADI_MT',0,ZHOOK_HANDLE)
      KREP=0
      LDREDF=.FALSE.
      LDMODC=.FALSE.
C
      IF (KTYPTR .LE. 0 ) THEN
         LLMLAM = .TRUE.
         INEXPL=8
         INGEOM=18
      ELSE
         LLMLAM = .FALSE.
      ENDIF
C
      INPAHE=(1+KNLATI)/2
C
      IF (KPHASE.EQ.2) THEN
        GOTO 200
      ELSEIF (KPHASE.EQ.3) THEN
        GOTO 300
      ELSEIF (KPHASE.LT.0.OR.KPHASE.GT.3) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C**
C     1.  -  CONTROLE DES VARIABLES SIMPLES (SYNTAXE ET COHERENCE).
C            (sauf pression de reference)
C-----------------------------------------------------------------------
C
      ILCDNO=LEN (CDNOMC)
      KLNOMC=1
C
      IF (ILCDNO.LE.0) THEN
        KREP=-65
        GOTO 1001
      ELSEIF (CDNOMC.EQ.' ') THEN
        KREP=-68
        GOTO 1001
      ELSEIF (KGARDE.LT.0.OR.KGARDE.GT.2) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C
      DO 101 J=ILCDNO,1,-1
C
      IF (CDNOMC(J:J).NE.' ') THEN
        KLNOMC=J
        GOTO 102
      ENDIF
C
  101 CONTINUE
C
  102 CONTINUE
C
      IF (KLNOMC.GT.FA%NCPCAD) THEN
        KREP=-65
        GOTO 1001
      ENDIF
C
      IF (KTRONC.LE.0.OR.KTRONC.GT.FA%NXTRON) THEN
        KREP=-70
        GOTO 1001
      ELSEIF (KNLATI.LE.0.OR.KNLATI.GT.FA%NXLATI) THEN
        KREP=-71
        GOTO 1001
      ELSEIF (KNIVER.LE.0.OR.KNIVER.GT.FA%NXNIVV) THEN
        KREP=-72
        GOTO 1001
      ELSEIF (KNXLON.LE.0.OR.KNXLON.GT.FA%NXLONG) THEN
        KREP=-83
        GOTO 1001
      ENDIF

      IF (LLMLAM) THEN
        IF (-2*KTYPTR+1.GT.KNXLON) THEN
          KREP=-115
          GOTO 1001
        ELSEIF (2*KTRONC+1.GT.KNLATI) THEN
          KREP=-116
          GOTO 1001
        ENDIF
      ELSE
        IF (PCODIL.LT.1.) THEN
          KREP=-73
          GOTO 1001
        ELSEIF (KTYPTR.LE.0.OR.KTYPTR.GT.FA%NTYPTX) THEN
          KREP=-109
          GOTO 1001
        ELSEIF (MAX(ABS(PSLAPO),ABS(PCLOPO),ABS(PSLOPO)).GT.1.) THEN
          KREP=-100
          GOTO 1001
        ELSEIF (ABS (1.-(PCLOPO**2+PSLOPO**2)).GT.1.E-5) THEN
          KREP=-101
          GOTO 1001
        ELSEIF (2*KTRONC+1.GT.KNXLON) THEN
          KREP=-84
          GOTO 1001
        ELSEIF (2*KTRONC+1.GT.4*(KNLATI/2)) THEN
C
C       Le test ci-dessus est "dur" car il fait l'hypothese que,
C     dans le cas ou KNLATI est impair, la grille comporte les poles.
C
          KREP=-79
          GOTO 1001
        ENDIF
      ENDIF
C
      IF (KPHASE.EQ.1) THEN
        GOTO 1001
      ENDIF
C**
C     2.  -  CONTROLE DES TABLEAUX (SYNTAXE ET COHERENCE).
C            (et de la pression de reference)
C-----------------------------------------------------------------------
C
  200 CONTINUE
C
C
      IF (PREFER.LT.0..OR.PREFER.GT.REAL (10*FA%MPRESX)) THEN
        KREP=-108
        GOTO 1001
      ENDIF
C
C     No Mount Everest test
C
      IF (.FALSE.) THEN
      DO 202 J=0,KNIVER
      IPREC=MAX (0,J-1)
      ZMIN=MIN (PAHYBR(J),PBHYBR(J))
      ZPMIN=PREFER*PAHYBR(J)+FA%SPSMIN*PBHYBR(J)
      ZPMAX=PREFER*PAHYBR(J)+FA%SPSMAX*PBHYBR(J)
      ZPMINP=PREFER*PAHYBR(IPREC)+FA%SPSMIN*PBHYBR(IPREC)
      ZPMAXP=PREFER*PAHYBR(IPREC)+FA%SPSMAX*PBHYBR(IPREC)
C
      IF (ZMIN.LT.0..OR.PBHYBR(J).GT.1.) THEN
        KREP=-80
        GOTO 1001
      ELSEIF (J.NE.0.AND.(PBHYBR(J).LT.PBHYBR(IPREC).OR.
     S                 ZPMIN.LE.ZPMINP.OR.ZPMAX.LE.ZPMAXP)) THEN
        KREP=-81
        GOTO 1001
      ENDIF
C
  202 CONTINUE
      ENDIF ! No Mount Everest test
C
      IF (.NOT.LLMLAM) THEN
C
         DO 201 J=1,INPAHE
         IPREC=MAX (1,J-1)
C
         IF (KNLOPA(J).LE.0.OR.KNLOPA(J).GT.KNXLON) THEN
           KREP=-74
           GOTO 1001
         ELSEIF (KNLOPA(J).LT.KNLOPA(IPREC)) THEN
           KREP=-75
           GOTO 1001
         ELSEIF (KNOZPA(J).LT.0.OR.KNOZPA(J).GT.KTRONC) THEN
           KREP=-76
           GOTO 1001
         ELSEIF (KNOZPA(J).LT.KNOZPA(IPREC)) THEN
           KREP=-77
           GOTO 1001
         ELSEIF ((2*KNOZPA(J)+1).GT.KNLOPA(J)) THEN
           KREP=-78
           GOTO 1001
         ELSEIF (ABS (PSINLA(J)).GT.1.) THEN
           KREP=-102
           GOTO 1001
         ELSEIF (PSINLA(J).GE.PSINLA(IPREC).AND.J.NE.1) THEN
           KREP=-103
           GOTO 1001
         ENDIF
C
  201    CONTINUE
C
      ELSE
C
C        *****  ERROR HANDLING FOR LAM CASE
C
         IF (ABS(KNLOPA(2)).GT.1) THEN
           KREP=-117
           GOTO 1001
         ELSEIF (KNLOPA(3).LE.0.OR.KNLOPA(3).GT.KNXLON) THEN
           KREP=-118
           GOTO 1001
         ELSEIF (KNLOPA(4).LT.KNLOPA(3).OR.KNLOPA(4).GT.KNXLON) THEN
           KREP=-119
           GOTO 1001
         ELSEIF (KNLOPA(5).LE.0.OR.KNLOPA(5).GT.KNLATI) THEN
           KREP=-120
           GOTO 1001
         ELSEIF (KNLOPA(6).LE.KNLOPA(5).OR.KNLOPA(6).GT.KNLATI) THEN
           KREP=-121
           GOTO 1001
         ELSEIF (2*KNLOPA(7).GT.(KNLOPA(4)-KNLOPA(3))) THEN
           KREP=-122
           GOTO 1001
         ELSEIF (2*KNLOPA(8).GT.(KNLOPA(6)-KNLOPA(5))) THEN
           KREP=-123
           GOTO 1001
         ENDIF
C
      ENDIF
C
      IF (KPHASE.EQ.2) GOTO 1001
C**
C     3.  -  CONTROLES LIES A LA DEFINITION DU CADRE PROPREMENT DITE.
C-----------------------------------------------------------------------
C
  300 CONTINUE
C
C        Le nom de cadre specifie est-il deja defini ?
C
      CALL FANUCA_MT (FA, CDNOMC,KRANGC,.FALSE.)
      LDREDF=KRANGC.NE.0
      IF (LDREDF) GOTO 500
C
C        En arrivant ici, il s'agit donc d'un nouveau cadre.
C
      IF (FA%NCADEF.GE.FA%JPNXCA) THEN
C
C        Trop de cadres deja definis pour en stocker un de plus.
C
        KREP=-56
        GOTO 1001
      ENDIF
C
C       Recherche d'un emplacement disponible dans les tables de cadres,
C     lequel devrait en bonne logique exister...
C
      DO 302 J=1,FA%JPNXCA
C
      IF (FA%CNOMCA(J).EQ.' ') THEN
        KRANGC=J
        GOTO 303
      ENDIF
C
  302 CONTINUE
C
      KREP=-66
      GOTO 1001
C
  303 CONTINUE
C
C           Nouveau cadre, mise a jour des tables partagees de cadres.
C
      FA%NCADEF=FA%NCADEF+1
      FA%NCAIND(FA%NCADEF)=KRANGC
      FA%CNOMCA(KRANGC)=CDNOMC
      FA%NLCCAD(KRANGC)=KLNOMC
C**
C     4.  -  STOCKAGE DES PARAMETRES DU CADRE (NOUVEAU, OU REDEFINI).
C-----------------------------------------------------------------------
C
  400 CONTINUE
C
      FA%NULCAD(KRANGC)=0
      FA%NTYPTR(KRANGC)=KTYPTR
      FA%MTRONC(KRANGC)=KTRONC
      FA%NNIVER(KRANGC)=KNIVER
      FA%NLATIT(KRANGC)=KNLATI
      FA%NXLOPA(KRANGC)=KNXLON
      FA%SSLAPO(KRANGC)=PSLAPO
      FA%SCLOPO(KRANGC)=PCLOPO
      FA%SSLOPO(KRANGC)=PSLOPO
      FA%SCODIL(KRANGC)=PCODIL
      FA%SPREFE(KRANGC)=PREFER
C
      FA%LIMLAM(KRANGC)=LLMLAM
      FA%NSFLAM(KRANGC)=0
C
      IF (.NOT.LDREDF.OR.KGARDE.NE.1) FA%NGARDE(KRANGC)=KGARDE
C
      IF (.NOT.LLMLAM) THEN
         ICOMPT=0
C
         DO 401 J=1,INPAHE
         ICOMPT=ICOMPT+KNLOPA(J)
         FA%NLOPAR(J,KRANGC)=KNLOPA(J)
         FA%NOZPAR(J,KRANGC)=KNOZPA(J)
         FA%SINLAT(J,KRANGC)=PSINLA(J)
  401    CONTINUE
C
         IF (KNLATI.EQ.2*INPAHE) THEN
           FA%NVAPDG(KRANGC)=ICOMPT*2
         ELSE
           FA%NVAPDG(KRANGC)=ICOMPT*2-KNLOPA(INPAHE)
         ENDIF
C
      ELSE
C *****  CALCULATION OF KNOZPA(), THEN ALSO SETTING OF FACOM1-TABLES  *****
C
         ZMSMAX=REAL(-KTYPTR)
         IMSMAX = -KTYPTR
         ISFLAM = 0
         CALL ELLIPS(KTRONC,IMSMAX,IKNTMP,IKMTMP)
CDP      CALL ELLIPS(IMSMAX,KTRONC,IKNTMP,IKMTMP)
C
C Initialisation de FA%NOMPAR (du module FAMODU)
C
         FA%NOMPAR(2,KRANGC) = 0
         DO JL=0,IMSMAX
           FA%NOMPAR(2*JL+3,KRANGC) = FA%NOMPAR(2*JL+2,KRANGC) + 1
           FA%NOMPAR(2*JL+4,KRANGC) = FA%NOMPAR(2*JL+3,KRANGC)
     S                             + 4*(IKNTMP(JL)+1) -1
         ENDDO
         FA%NOMPAR(1,KRANGC) = KTRONC
         FA%NOMPAR(2,KRANGC) = IMSMAX
C
         DO 717 JL=0,KTRONC
            IK=IKMTMP(JL)
CDP         IK=IKNTMP(JL)
            ICPL4N(JL)=4*(IK+1)
            ISFLAM = ISFLAM + 4*(IK+1)
  717    CONTINUE
C
         IESN0(0)=1
C
         DO 771 J=1,KTRONC
            IESN0(J)=IESN0(J-1)+ICPL4N(J-1)
  771    CONTINUE
C
C -----  NOW SETTING OF TABLES  -----
         DO 707 J=1,INEXPL
            FA%NLOPAR(J,KRANGC)=KNLOPA(J)
  707    CONTINUE
         DO 770 J=1,INGEOM
            FA%SINLAT(J,KRANGC)=PSINLA(J)
  770    CONTINUE
         FA%NOZPAR(1,KRANGC)=KTRONC
         FA%NOZPAR(2,KRANGC)=IMSMAX
C
         DO 700 J=0,KTRONC
            FA%NOZPAR(2*J+3,KRANGC)=IESN0(J)
            FA%NOZPAR(2*J+4,KRANGC)=IESN0(J)+ICPL4N(J)-1
  700    CONTINUE
  
         IF (FA%NOZPAR(2*KTRONC+4,KRANGC).NE.
     S       FA%NOMPAR(2*IMSMAX+4,KRANGC))
     S   THEN
           KREP=-127
           GOTO 1001
         ENDIF
C
         FA%NSFLAM(KRANGC)=ISFLAM
C
C *****  DETERMINATION OF FA%NVAPDG()  *****
C
         FA%NVAPDG(KRANGC)=KNLATI*KNXLON
C
      ENDIF
C
      DO 402 J=0,KNIVER
      FA%SFOHYB(1,J,KRANGC)=PAHYBR(J)
      FA%SFOHYB(2,J,KRANGC)=PBHYBR(J)
  402 CONTINUE
C
      GOTO 1001
C**
C     5.  -  TENTATIVE DE REDEFINITION D'UN CADRE. CONTROLES AD HOC.
C-----------------------------------------------------------------------
C
  500 CONTINUE
C
      IF (FA%MTRONC(KRANGC).NE.KTRONC.OR.FA%NNIVER(KRANGC).NE.KNIVER.OR.
     S    FA%NLATIT(KRANGC).NE.KNLATI.OR.FA%NXLOPA(KRANGC).NE.KNXLON.OR.
     S    (ABS(FA%SSLAPO(KRANGC)-PSLAPO)>ZEPS) .OR.
     S    (ABS(FA%SCLOPO(KRANGC)-PCLOPO)>ZEPS) .OR.
     S    (ABS(FA%SSLOPO(KRANGC)-PSLOPO)>ZEPS) .OR.
     S    (ABS(FA%SCODIL(KRANGC)-PCODIL)>ZEPS) .OR.      
     S    FA%NTYPTR(KRANGC).NE.KTYPTR.OR.
     S    (ABS(FA%SPREFE(KRANGC)-PREFER)>ZEPS)) GOTO 505
C
      IF (.NOT.LLMLAM) THEN
         DO 5010 J=1,INPAHE
         IF (FA%NLOPAR(J,KRANGC).NE.KNLOPA(J).OR.
     S       FA%NOZPAR(J,KRANGC).NE.KNOZPA(J).OR.
     S       (ABS(FA%SINLAT(J,KRANGC)-PSINLA(J))>ZEPS)) GOTO 505
 5010    CONTINUE
      ELSE
         DO 5011 J=1,8
         IF (FA%NLOPAR(J,KRANGC).NE.KNLOPA(J)) GOTO 505
 5011    CONTINUE
         DO 5012 J=1,18
         IF (ABS(FA%SINLAT(J,KRANGC)-PSINLA(J))>ZEPS) GOTO 505
 5012    CONTINUE
      ENDIF
C
      DO 502 J=0,KNIVER
      IF ((ABS(FA%SFOHYB(1,J,KRANGC)-PAHYBR(J))>ZEPS).OR.
     S    (ABS(FA%SFOHYB(2,J,KRANGC)-PBHYBR(J))>ZEPS)) GOTO 505
  502 CONTINUE
C
C        Si on arrive ici, il y a redefinition a l'identique,
C     du moins pour les parametres numeriques.
C        L'option de conservation du cadre peut, elle, etre modifiee
C     dans le cas d'une definition non dynamique.
C
      IF (KGARDE.NE.1) FA%NGARDE(KRANGC)=KGARDE
      GOTO 1001
C
  505 CONTINUE
      LDMODC=.TRUE.
C
C        Il y a donc redefinition avec changement de parametre(s),
C     ce qui n'est possible que s'il n'y a pas de fichier rattache,
C     et s'il ne s'agit pas d'une definition dynamique de cadre
C     (appel par FAITOU avec KGARDE=1).
C
      IF (KGARDE.EQ.1) THEN
        KREP=-58
      ELSEIF (FA%NULCAD(KRANGC).NE.0) THEN
        KREP=-59
      ELSE
        GOTO 400
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE EVENTUELLE,
C            VIA LE sous-programme "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
C
      LLFATA=KREP.NE.0.AND.FA%NRFAGA.NE.2
C
      IF (FA%LFAMOP.OR.LLFATA) THEN
        INIMES=2
        CLNSPR='FACADI'
        INUMER=FA%JPNIIL
C
        IF (KREP.EQ.-65.AND.ILCDNO.LE.0) THEN
          ILNOMC=8
          CLACTI(1:ILNOMC)=FA%CHAINC(:ILNOMC)
        ELSE
          ILNOMC=MIN (KLNOMC,FA%NCPCAD,LEN (CLACTI))
          CLACTI(1:ILNOMC)=CDNOMC(1:ILNOMC)
        ENDIF
C
        WRITE (UNIT=CLMESS,FMT='(''ARGUM.SIMPLES='',I4,'','''''',A,
     S         '''''''',4('','',F7.4),4('','',I4),'','',F10.3,
     S         2('','',L1),2('','',I2),'','',I3,'','',I1)')
     S  KREP,CLACTI(1:ILNOMC),PSLAPO,PCLOPO,PSLOPO,PCODIL,
     S  KTRONC,KNLATI,KNXLON,KNIVER,PREFER,LDMODC,LDREDF,KPHASE,
     S  KRANGC,KLNOMC,KGARDE
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI(1:ILNOMC),.FALSE.)
      ELSEIF (KTRONC.LE.FA%NSTROI.AND.(KPHASE.EQ.0.OR.KPHASE.EQ.1)) THEN
        INIMES=1
        CLNSPR='FACADI'
        INUMER=FA%JPNIIL
        ILNOMC=MIN (KLNOMC,FA%NCPCAD)
        WRITE (UNIT=CLMESS,
     S         FMT='(''TRONCATURE ('',I2,'') INFERIEURE '',
     S ''OU EGALE A LA SOUS-TRONCATURE "NON COMPACTEE" IMPLICITE ('',I2,
     S ''), CADRE '''''',A,'''''''')') KTRONC,FA%NSTROI,CDNOMC(1:ILNOMC)
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FACADI_MT',1,ZHOOK_HANDLE)
      END

