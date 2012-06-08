C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAVEUR_MT (FA,  KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP, 
     S                    KSTRON, KPUILA, KDMOPL )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'obtenir, pour un fichier ARPEGE
C     ouvert, les options courantes liees au codage GRIB des champs.
C     CES OPTIONS NE SONT UTILISEES QUE POUR (RE)ECRIRE DES CHAMPS
C     codes en GRIB.
C       ( Visualisation (?) options Effectives pour l'UtilisateuR )
C**
C     Arguments : KREP   ==> Code-reponse du sous-programme;
C     (tous de    KNUMER ==> Numero d'Unite Logique concernee;
C      SORTIE)    KNGRIB ==> Niveau de codage GRIB (-1,0,1,2,3);
C                 KNBPDG ==> Nombre de bits par valeur point-de-grille;
C                 KNBCSP ==> Nombre de bits par partie reelle/imaginaire
C                            de coefficient spectral;
C                 KSTRON ==> Sous-troncature non compactee;
C                 KPUILA ==> Puissance de laplacien;
C                 KDMOPL ==> Degre de modulation de KPUILA.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA
      INTEGER KDMOPL
C
      INTEGER IREP, IRANG, INIMES
C
      LOGICAL LLVERF
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAVEUR_MT',0,ZHOOK_HANDLE)
      LLVERF=.FALSE.
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
C**
C     2.  -  RECOPIE DES VALEURS EN COMMON DANS LES ARGUMENTS.
C-----------------------------------------------------------------------
C
      KNGRIB=FA%NFGRIB(IRANG)
      KNBPDG=FA%NBFPDG(IRANG)
      KNBCSP=FA%NBFCSP(IRANG)
      KSTRON=FA%NSTROF(IRANG)
      KPUILA=FA%NPUFLA(IRANG)
      KDMOPL=FA%NMFDPL(IRANG)
      IREP=0
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
        IF (LHOOK) CALL DR_HOOK('FAVEUR_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FAVEUR'
C
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', KNGRIB='',I2,'', KNBPDG='',I3,'', KNBCSP='',I3,
     S       '', KSTRON='',I2,'', KPUILA='',I3,'', KDMOPL='',I3)')
     S   KREP,KNUMER,KNGRIB,KNBPDG,KNBCSP,KSTRON,KPUILA,KDMOPL
      CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CLACTI,.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FAVEUR_MT',1,ZHOOK_HANDLE)
      END

