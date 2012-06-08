C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAVORI_MT (FA,  KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA,
     S                    KDMOPL )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet de connaitre les options implicites
C     courantes liees au codage GRIB des champs.
C     CES OPTIONS NE SONT UTILISEES QUE POUR (RE)ECRIRE DES CHAMPS
C     codes en GRIB, et les valeurs implicites ne servent
C     que LORS d'une OUVERTURE de FICHIER.
C       ( Visualisation (?) Options R (??) Implicites )
C**
C     Arguments : KNGRIB ==> Niveau de codage GRIB (-1,0,1,2,3);
C     (tous de    KNBPDG ==> Nombre de bits par valeur point-de-grille;
C      SORTIE)    KNBCSP ==> Nombre de bits par partie reelle/imaginaire
C                            de coeff. spectral;
C                 KSTRON ==> Sous-troncature non compactee;
C                 KPUILA ==> Puissance de laplacien;
C                 KDMOPL ==> Degre de modulation de KPUILA.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA, KDMOPL
C
      INTEGER IREP, INIMES, INUMER
C
      LOGICAL LLVERG
C
C
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAVORI_MT',0,ZHOOK_HANDLE)
      IF (FA%FAVORI_LLPREA) THEN
C
C          A la premiere utilisation, appel au sous-programme "FARINE".
C
        CALL FARINE_MT(FA, 2)
        FA%FAVORI_LLPREA=.FALSE.
      ENDIF
C
C         Verrouillage global eventuel.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
      LLVERG=FA%LFAMUL
C**
C     2.  -  RECOPIE DES VALEURS EN COMMON DANS LES ARGUMENTS.
C-----------------------------------------------------------------------
C
      KNGRIB=FA%NIGRIB
      KNBPDG=FA%NBIPDG
      KNBCSP=FA%NBICSP
      KSTRON=FA%NSTROI
      KPUILA=FA%NPUILA
      KDMOPL=FA%NMIDPL
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE VIA LE SOUS-PROGRAMME "FAIPAR"
C-----------------------------------------------------------------------
C
 1001 CONTINUE
C
C        Deverrouillage global eventuel.
C
      IF (LLVERG) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'OFF')
C
      IF (FA%NIMSGA.EQ.2) THEN
        IREP=0
        INIMES=2
        CLNSPR='FAVORI'
        WRITE (UNIT=CLMESS,FMT='(''KNGRIB='',I2,'', KNBPDG='',I3,
     S         '', KNBCSP='',I3,'', KSTRON='',I3,'', KPUILA='',I3,
     S         '', KDMOPL='',I3)')
     S    KNGRIB,KNBPDG,KNBCSP,KSTRON,KPUILA,KDMOPL
        INUMER=FA%JPNIIL
        CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAVORI_MT',1,ZHOOK_HANDLE)
      END

