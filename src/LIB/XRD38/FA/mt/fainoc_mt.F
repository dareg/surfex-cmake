C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAINOC_MT (FA,  KRANG )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'INterpreter, pour un fichier ARPEGE
C     ouvert, les Options par defaut (-1) du Codage GRIB des champs:
C     FA%NBFPDG(KRANG), FA%NBFCSP(KRANG), FA%NSTROF(KRANG), FA%NPUFLA(KRANG).
C     Cette routine doit etre appelee par FAITOU ou FANOUV ou FAGOTE
C     pour ne pas laisser le defaut -1 lors du decodage ou du codage
C     GRIB.
C
C**
C     Arguments : KRANG  (Entree) ==> Rang de l'unite logique;
C
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KRANG
C
      INTEGER IRANGC, ITRONC, INBITS, ITYPTR, IAUXIL, IREP, INIMES
      INTEGER INUMER
C
      LOGICAL LLVERF, LLMLAM
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  INITIALISATIONS PREALABLES.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAINOC_MT',0,ZHOOK_HANDLE)
      LLVERF=.FALSE.
C
C         Verrouillage eventuel du fichier.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(KRANG),'ON')
      LLVERF=FA%LFAMUL
C
      IRANGC=FA%NUCADR(KRANG)
      LLMLAM=FA%LIMLAM(IRANGC)
      ITRONC=FA%MTRONC(IRANGC)
      ITYPTR=FA%NTYPTR(IRANGC)
C
C**
C     2.  -  INTERPRETATION DES OPTIONS PAR DEFAUT.
C-----------------------------------------------------------------------
C
C On distingue le cas ARPEGE du cas ALADIN (LLMLAM=.T.).
C
C
C Evaluation du nombre de bits par valeur point-de-grille
C
      IF (FA%NBFPDG(KRANG).LT.0) THEN
        IF (LLMLAM) THEN
          FA%NBFPDG(KRANG)=16
        ELSE
          FA%NBFPDG(KRANG)=16
        ENDIF
      ENDIF
C
C Evaluation du nombre de bits par partie reelle/imagin. de coeff. spectral
C
      IF (FA%NBFCSP(KRANG).LT.0) THEN
        IF (LLMLAM) THEN
          FA%NBFCSP(KRANG)=18
        ELSE
          FA%NBFCSP(KRANG)=16
        ENDIF
      ENDIF
C
C Evaluation de la sous-troncature non compactee
C
      IF (FA%NSTROF(KRANG).LT.0) THEN
        INBITS=FA%NBFCSP(KRANG)
        IF (LLMLAM) THEN
          IAUXIL=MAX ( ITRONC, -ITYPTR )
          IAUXIL=MAX ( 10, ((1+IAUXIL)*25)/(10*INBITS), (1+IAUXIL)/10 )
          IAUXIL=MIN ( IAUXIL, ITRONC-1, -ITYPTR-1 )
          FA%NSTROF(KRANG)=IAUXIL
        ELSE
          IAUXIL=MAX ( 10, 480/INBITS-10, (1+ITRONC)/10 )
          IAUXIL=MIN ( IAUXIL, ITRONC-1 )
          FA%NSTROF(KRANG)=IAUXIL
        ENDIF
      ENDIF
C
C Evaluation de la puissance de laplacien
C
      IF (FA%NPUFLA(KRANG).LT.0) THEN
        IF (LLMLAM) THEN
          FA%NPUFLA(KRANG)=2
        ELSE
          FA%NPUFLA(KRANG)=1
        ENDIF
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
C
C        Deverrouillage eventuel du fichier.
C
      IF (LLVERF) CALL LFIVER_MT (FA%LFI, FA%VRFICH(KRANG),'OFF')
C
      IF (IXNVMS(KRANG).NE.2)  THEN 
        IF (LHOOK) CALL DR_HOOK('FAINOC_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      INUMER=FA%JPNIIL
      INIMES=IXNVMS (KRANG)
      IREP=0
      CLNSPR='FAINOC'
      WRITE (UNIT=CLMESS,FMT='(''KRANG='',I4)') KRANG
      CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,.FALSE.,CLMESS,
     S                CLNSPR,' ',.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FAINOC_MT',1,ZHOOK_HANDLE)
      END

