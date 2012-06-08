C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FANION_MT (FA,  KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, 
     S                      LDEXIS,
     S                      LDCOSP, KNGRIB, KNBITS, KSTRON, KPUILA )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme renseignant sur l'EXISTENCE et les CARACTERISTI-
C     QUES eventuelles d'un Article de type CHAMP dans un Fichier ARPEGE
C       ( LDEXIS est le "fanion" leve si l'article existe )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
C                KNIVAU (Entree) ==> Niveau vertical eventuel;
C                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
C                LDEXIS (Sortie) ==> Vrai si l'article de type CHAMP
C                                    existe bien dans le Fichier;
C                LDCOSP (Sortie) ==> Vrai si le champ est represente
C                                    par des coefficients spectraux;
C                KNGRIB (Sortie) ==> Niveau de codage GRIB;
C                KNBITS (Sortie) ==> Nombre de bits de codage eventuel;
C                KSTRON (Sortie) ==> Sous-troncature non codee   " -le;
C                KPUILA (Sortie) ==> Puissance de laplacien eventuelle.
C
C        KNBITS n'a de sens que si l'article existe et a ete code;
C     de meme pour KSTRON et KPUILA, qui ne sont applicables qu'a un
C     champ represente en coefficients spectraux.
C        Les arguments de sortie n'ayant pas de sens sont mis a
C     0 pour les entiers, .FALSE. pour les logiques.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER, KNIVAU, KNGRIB, KNBITS, KSTRON, KPUILA
C
      INTEGER IREP, ILPRFU, ILSUFU, ILNOMU, ILONGA, IRANG, INIMES
      INTEGER ILPREF, ILSUFF, IPOFIN, IPOSEX, IRANGC, ILCHAM
C
      INTEGER (KIND=JPDBLE) IVALCO (5)
      INTEGER  IB1PAR (3)
C
      LOGICAL LLVERF, LLRLFI, LDCOSP, LDEXIS, LLTEMP, LLNOMU, LLMLAM
C
      CHARACTER CDPREF*(*), CDSUFF*(*)
      CHARACTER CLPREF*(FA%JPXNOM), CLSUFF*(FA%JPXSUF)
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FANION_MT',0,ZHOOK_HANDLE)
      IREP=0
      LLVERF=.FALSE.
      LLTEMP=.FALSE.
      LLRLFI=.FALSE.
      LLNOMU=.FALSE.
      ILPRFU=LEN (CDPREF)
      ILSUFU=LEN (CDSUFF)
      LDEXIS=.FALSE.
      LDCOSP=.FALSE.
      KNGRIB=0
      KNBITS=0
      KSTRON=0
      KPUILA=0
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ENDIF
C
C         Verrouillage eventuel du fichier.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'ON')
      LLVERF=FA%LFAMUL
C
      IF (FA%LCREAF(IRANG)) GOTO 1001
C**
C     2.  -  FABRICATION DU NOM D'ARTICLE VIA LE SOUS-PROGRAMME "FANFAR"
C            ( controles de CDPREF, KNIVAU, CDSUFF inclus )
C-----------------------------------------------------------------------
C
      CALL FANFAR_MT (FA, IREP,IRANG,CDPREF,KNIVAU,CDSUFF,CLNOMA,
     S                IB1PAR,ILPRFU,ILSUFU,ILNOMU)
      IF (IREP.NE.0) GOTO 1001
      LLNOMU=.TRUE.
C**
C     3.  -  RECHERCHE DE L'ARTICLE SUR LE FICHIER, LECTURE PARTIELLE.
C-----------------------------------------------------------------------
C
      CALL LFINFO_MT (FA%LFI, IREP,KNUMER,CLNOMA(1:ILNOMU),
     S             ILONGA,IPOSEX)
      LLRLFI=IREP.NE.0
      IF (LLRLFI.OR.ILONGA.EQ.0) GOTO 1001
      LDEXIS=.TRUE.
C
      IF (ILONGA.GT.FA%JPXCHA+2) THEN
        IREP=-90
        GOTO 1001
      ENDIF
C
      IF (FA%LERRFA(IRANG)) THEN
C
C        Le fichier est gere en mode "toute erreur est fatale".
C     Ce mode etant normalement couple au mode correspondant du logiciel
C     LFI, on va temporairement annuler l'option LFI afin de pouvoir
C     faire une lecture partielle de l'entete de l'article Champ.
C
        CALL LFIERF_MT (FA%LFI, IREP,KNUMER,.FALSE.)
        LLRLFI=IREP.NE.0
        IF (LLRLFI) GOTO 1001
        LLTEMP=.TRUE.
      ENDIF
C
      CALL LFILEC_MT (FA%LFI, IREP,KNUMER,CLNOMA(1:ILNOMU),IVALCO,5)
C
      IF (IREP.EQ.0) THEN
        IREP=-93
        GOTO 1001
      ELSEIF (IREP.NE.-21) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ELSEIF (IVALCO(1).LT.-1 .OR. IVALCO(1).GT.3 .OR.
     S        IVALCO(2).LT.0  .OR. IVALCO(2).GT.1 .OR.
     S  (IVALCO(1).GT.0 .AND. IVALCO(2).EQ.1 .AND. IVALCO(4).LT.0)) THEN
        IREP=-91
        GOTO 1001
      ELSE
        IREP=0
        KNGRIB=IVALCO(1)
        LDCOSP=IVALCO(2).EQ.1
      ENDIF
C
      IRANGC=FA%NUCADR(IRANG)
      LLMLAM=FA%LIMLAM(IRANGC)
C
      IF (LDCOSP) THEN
        IF (LLMLAM) THEN
           ILCHAM=FA%NSFLAM(IRANGC)
        ELSE
           IF (KNGRIB.EQ.3 .OR. KNGRIB.EQ.-1) THEN
             ILCHAM=(1+FA%MTRONC(IRANGC))*(2+FA%MTRONC(IRANGC))
           ELSE
             ILCHAM=(1+FA%MTRONC(IRANGC))**2
           ENDIF
        ENDIF
      ELSE
        ILCHAM=FA%NVAPDG(IRANGC)
      ENDIF
C
      IF (KNGRIB.EQ.-1 .OR. KNGRIB.EQ.0) THEN
C
C          Cas ou il n'y a aucun codage... controle longueur d'article
C
        IF (ILONGA.LT.(ILCHAM+2)) THEN
          IREP=-93
          GOTO 1001
        ELSEIF (ILONGA.GT.(ILCHAM+2)) THEN
          IREP=-94
          IF (LLMOER(IREP,IRANG)) GOTO 1001
        ENDIF
C
      ELSE
C
C        Cas avec codage GRIB (standard ou non).
C
        KNBITS=IVALCO(3)
C
        IF (LDCOSP) THEN
          KSTRON=IVALCO(4)
          KPUILA=IVALCO(5)
C
          IF (KNGRIB.EQ.2.AND.ILONGA.LT.(5+(1+KSTRON)**2)) THEN
            IREP=-93
            GOTO 1001
          ENDIF
C
        ENDIF
C
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
C
      IF (LLTEMP) THEN
C
C         On remet le fichier en mode "toute erreur fatale" au niveau
C     du logiciel LFI.
C
        CALL LFIERF_MT (FA%LFI, IREP,KNUMER,.TRUE.)
        LLRLFI=IREP.NE.0
      ENDIF
C
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANG)
C
C        Deverrouillage eventuel du fichier.
C
      IF (LLVERF) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'OFF')
C
      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=IXNVMS(IRANG)
      ENDIF
C
      IF (.NOT.LLFATA.AND.INIMES.NE.2)  THEN 
        IF (LHOOK) CALL DR_HOOK('FANION_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FANION'
C
      IF (ILPRFU.GE.1) THEN
        ILPREF=MIN (ILPRFU,LEN (CLPREF))
        CLPREF(1:ILPREF)=CDPREF(1:ILPREF)
      ELSE
        ILPREF=8
        CLPREF(1:ILPREF)=FA%CHAINC(:ILPREF)
      ENDIF
C
      IF (ILSUFU.GE.1) THEN
        ILSUFF=MIN (ILSUFU,LEN (CLSUFF))
        CLSUFF(1:ILSUFF)=CDSUFF(1:ILSUFF)
      ELSE
        ILSUFF=8
        CLSUFF(1:ILSUFF)=FA%CHAINC(:ILSUFF)
      ENDIF
C
      IF (.NOT.LLNOMU) THEN
        ILNOMU=MIN (ILPREF,FA%NCPCAD)
        CLNOMA(1:ILNOMU)=CLPREF(1:ILPREF)
      ENDIF
C
      WRITE (UNIT=CLMESS,
     S       FMT='(''ARGUMENTS:'',I4,'','',I3,'','''''',A,
     S       '''''','',I6,'','''''',A,'''''', LDEXIS= '',L1,
     S       '', LDCOSP= '',L1,'', KNGRIB='',I2,'', KNBITS='',I3,
     S       '',KSTRON='',I3,'',KPUILA='',I6)')
     S   KREP,KNUMER,CLPREF(1:ILPREF),KNIVAU,CLSUFF(1:ILSUFF),LDEXIS,
     S   LDCOSP,KNGRIB,KNBITS,KSTRON,KPUILA
      CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CLNOMA(1:ILNOMU),LLRLFI)
C
      IF (LHOOK) CALL DR_HOOK('FANION_MT',1,ZHOOK_HANDLE)
      END

