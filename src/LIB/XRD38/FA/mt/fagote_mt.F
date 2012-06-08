C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAGOTE_MT (FA,  KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP, 
     S                      KSTRON, KPUILA, KDMOPL )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'ajuster, pour un fichier ARPEGE
C     ouvert, les options liees au codage GRIB des champs.
C     CES OPTIONS NE SONT UTILISEES QUE POUR (RE)ECRIRE DES CHAMPS
C     codes en GRIB.
C       ( Grib, Options Techniques Effectives )
C**
C     Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                 KNUMER (Entree) ==> Numero d'Unite Logique concernee;
C                 KNGRIB (Entree) ==> Niveau de codage GRIB (-1,0,1,2,3);
C                 KNBPDG (Entree) ==> Nombre de bits par valeur point-
C                                     de-grille;
C                 KNBCSP (Entree) ==> Nombre de bits par partie reelle/
C                                     imaginaire de coeff. spectral;
C                 KSTRON (Entree) ==> Sous-troncature non compactee;
C                 KPUILA (Entree) ==> Puissance de laplacien;
C                 KDMOPL (Entree) ==> Degre de modulation de KPUILA.
C
C     Remarque:  KSTRON egal a -1 est accepte; dans ce cas on indexera
C                (pour chaque champ spectral ecrit) la sous-troncature
C                effective sur la troncature et sur le nombre de bits
C                par valeur compactee.
C                 
C     MODIF : 30/03/2007 JM AUDOIN  FA%LFAMOP Pour limiter Impression
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA
      INTEGER KDMOPL
C
      INTEGER IMINIM, IREP, IRANGC, ITRONC, J, INIMES, IRANG, ITYPTR
C
      LOGICAL LLVERF, LLMLAM
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAGOTE_MT',0,ZHOOK_HANDLE)
      LLVERF=.FALSE.
      IMINIM=MIN (2+KNGRIB,2+KNBPDG,2+KNBCSP,2+KSTRON,1+KDMOPL)
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ELSEIF (IMINIM.LE.0) THEN
        IREP=-64
        GOTO 1001
      ELSEIF (KNBPDG*KNBCSP.EQ.0 .AND. KNGRIB.GT.0) THEN
        IREP=-124
        GOTO 1001
      ELSEIF (KNGRIB.GT.3) THEN
        IREP=-96
        GOTO 1001
      ELSEIF (MAX (KNBPDG,KNBCSP).GT.FA%NBIMAX) THEN
        IREP=-97
        GOTO 1001
      ELSEIF (IABS (KPUILA).GT.2**15-1) THEN
        IREP=-98
        GOTO 1001
      ENDIF
C
C         Verrouillage eventuel du fichier.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'ON')
      LLVERF=FA%LFAMUL
C
      IRANGC=FA%NUCADR(IRANG)
      ITRONC=FA%MTRONC(IRANGC)
      ITYPTR=FA%NTYPTR(IRANGC)
      LLMLAM=FA%LIMLAM(IRANGC)
C
      IF (KSTRON.GE.ITRONC) THEN
        IREP=-99
        GOTO 1001
      ELSEIF (ITYPTR.LT.0.AND.KSTRON.GE.(-ITYPTR)) THEN
        IREP=-99
        GOTO 1001
      ENDIF
C**
C     2.  -  STOCKAGE DES NOUVEAUX PARAMETRES.
C-----------------------------------------------------------------------
C
      IF (KPUILA.NE.FA%NPUFLA(IRANG)) THEN
        FA%NPUFLA(IRANG)=KPUILA
        FA%LIFLAP(IRANG)=.TRUE.
      ENDIF
      IF (KSTRON.NE.FA%NSTROF(IRANG)) THEN
        FA%NSTROF(IRANG)=KSTRON
        FA%LISC2F(IRANG)=.TRUE.
      ENDIF
C
      IF (FA%LFAMOP.AND.(FA%NFGRIB(IRANG).EQ.3
     S   .OR.FA%NFGRIB(IRANG).EQ.-1)
     S  .AND.(KNGRIB.LT.3.AND.KNGRIB.GT.-1))           THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S        'FAGOTE: WARNING!! Les champs spectraux NE devront',
     S        ' PAS etre ranges comme dans le modele (rangt horiz.)',
     S        ' pour l''unite logique ',KNUMER
        WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
      ENDIF
      IF (FA%LFAMOP.AND.(FA%NFGRIB(IRANG).LT.3
     S   .AND.FA%NFGRIB(IRANG).GT.-1)
     S  .AND.(KNGRIB.EQ.3.OR.KNGRIB.EQ.-1))           THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S        'FAGOTE: WARNING!! Les champs spectraux devront',
     S        ' etre ranges comme dans le modele (rangt verti.) pour',
     S        ' l''unite logique ',KNUMER
        WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
      ENDIF
      FA%NFGRIB(IRANG)=KNGRIB
      FA%NBFPDG(IRANG)=KNBPDG
      FA%NBFCSP(IRANG)=KNBCSP
      FA%NMFDPL(IRANG)=KDMOPL
      IREP=0
C
C Appel a FAINOC pour interpreter les eventuels defauts
C de -1 pris par FA%NBFPDG, FA%NBFCSP, FA%NSTROF et FA%NPUFLA en
C IRANG-ieme position.
C
      CALL FAINOC_MT (FA,  IRANG )
C
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
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
        IF (LHOOK) CALL DR_HOOK('FAGOTE_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FAGOTE'
C
C***** FAZZZZ - KREP=iiii, KNUMER=iii, KNGRIB=ii, KNBPDG=iii,   *****
C*****          KNBCSP=iii, KSTRON=ii, KPUILA=iii, KDMOPL=iii  *****
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', KNGRIB='',I2,'', KNBPDG='',I3,'', KNBCSP='',I3,
     S       '', KSTRON='',I2,'', KPUILA='',I3,'', KDMOPL='',I3)')
     S   KREP,KNUMER,KNGRIB,KNBPDG,KNBCSP,KSTRON,KPUILA,KDMOPL
      CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR, CLACTI,.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FAGOTE_MT',1,ZHOOK_HANDLE)
      END

