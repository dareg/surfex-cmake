C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAGIOT_MT (FA,  KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA,
     S                    KDMOPL )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet de modifier les options implicites
C     liees au codage GRIB des champs.
C     CES OPTIONS NE SONT UTILISEES QUE POUR (RE)ECRIRE DES CHAMPS
C     codes en GRIB, et les nouvelles valeurs implicites ne serviront
C     que LORS d'une OUVERTURE de FICHIER ULTERIEURE.
C       ( Grib, Implicites Options Techniques )
C**
C     Arguments : KNGRIB (Entree) ==> Niveau de codage GRIB (-1,0,1,2,3);
C                 KNBPDG (Entree) ==> Nombre de bits par valeur point-
C                                     de-grille;
C                 KNBCSP (Entree) ==> Nombre de bits par partie reelle/
C                                     imaginaire de coeff. spectral;
C                 KSTRON (Entree) ==> Sous-troncature non compactee;
C                 KPUILA (Entree) ==> Puissance de laplacien;
C                 KDMOPL (Entree) ==> Degre de modulation de KPUILA.
C
C     N.B.: Il doit y avoir coherence vis-a-vis des cadres deja definis
C           et vis-a-vis des limites usagers.
C           ( ce qui en pratique, ne concerne que KSTRON )
C
C     Remarque:  KSTRON egal a -1 est accepte, et dans ce cas
C                on indexera (a chaque ouverture de fichier) la sous-
C                troncature effective sur la troncature.
C
C      MODIF 30/03/2007 JM AUDOIN FA%LFAMOP pour limiter IMPRESSION
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA, KDMOPL
C
      INTEGER IMINIM, IREP, INIMES, INUMER, J, IRANGC
C
      LOGICAL LLVERG
C
C
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAGIOT_MT',0,ZHOOK_HANDLE)
      IF (FA%FAGIOT_LLPREA) THEN
C
C          A la premiere utilisation, appel au sous-programme "FARINE".
C
        CALL FARINE_MT(FA, 2)
        FA%FAGIOT_LLPREA=.FALSE.
      ENDIF
C
      LLVERG=.FALSE.
      IMINIM=MIN (2+KNGRIB,2+KNBPDG,2+KNBCSP,2+KSTRON,1+KDMOPL)
C
      IF (IMINIM.LE.0) THEN
        IREP=-64
        GOTO 1001
      ELSEIF (KNBPDG*KNBCSP.EQ.0) THEN
        IREP=-64
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
C         Verrouillage global eventuel.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
      LLVERG=FA%LFAMUL
C
      IF (KSTRON.GE.FA%NXTRON) THEN
        IREP=-113
        GOTO 1001
      ENDIF
C
C        Coherence de "KSTRON" vis-a-vis de la troncature des cadres
C     deja definis.
C
      DO 101 J=1,FA%NCADEF
      IRANGC=FA%NCAIND(J)
C
      IF (KSTRON.GE.FA%MTRONC(IRANGC)) THEN
        IREP=-99
        GOTO 1001
      ENDIF
C
  101 CONTINUE
C**
C     2.  -  STOCKAGE DES NOUVEAUX PARAMETRES.
C-----------------------------------------------------------------------
C
      IF (FA%LFAMOP.AND.(FA%NIGRIB.EQ.-1.OR.FA%NIGRIB.EQ.3).AND.
     S    (KNGRIB.GT.-1.AND.KNGRIB.LT.3))     THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S          'FAGIOT: CAUTION!! Les champs spectraux ARPEGE ne',
     S          ' devront pas etre ranges comme dans le modele ARPEGE'
        WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
      ENDIF
      IF (FA%LFAMOP.AND.(KNGRIB.EQ.-1.OR.KNGRIB.EQ.3).AND.
     S    (FA%NIGRIB.GT.-1.AND.FA%NIGRIB.LT.3))     THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S          'FAGIOT: CAUTION!! Les champs spectraux ARPEGE devront',
     S          ' etre ranges comme dans le modele ARPEGE'
        WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
      ENDIF
      FA%NIGRIB=KNGRIB
      FA%NBIPDG=KNBPDG
      FA%NBICSP=KNBCSP
      FA%NSTROI=KSTRON
      FA%NPUILA=KPUILA
      FA%NMIDPL=KDMOPL
      IREP=0
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      LLFATA=LLMOER (IREP,0)
C
C        Deverrouillage global eventuel.
C
      IF (LLVERG) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'OFF')
C
      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=FA%NIMSGA
      ENDIF
C
      IF (.NOT.LLFATA.AND.INIMES.NE.2)  THEN 
        IF (LHOOK) CALL DR_HOOK('FAGIOT_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FAGIOT'
C
      WRITE (UNIT=CLMESS,FMT='(''KNGRIB='',I2,'', KNBPDG='',I3,
     S       '', KNBCSP='',I3,'', KSTRON='',I2,'', KPUILA='',I3,
     S       '', KDMOPL='',I3)')
     S   KNGRIB,KNBPDG,KNBCSP,KSTRON,KPUILA,KDMOPL
      INUMER=FA%JPNIIL
      CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CLACTI,.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FAGIOT_MT',1,ZHOOK_HANDLE)
      END

