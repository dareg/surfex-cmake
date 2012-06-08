C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAIS2F_MT (FA,  KREP, KRANG )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Ce sous-programme initialise un tableau "reference" de
C      l'en-tete GRIB, section 2.
C      (routine appelee une seule fois pour un fichier Aladin donne)
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C*
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANG
C
      INTEGER IRANGC, JM, JMAX, ILOW, IADD, INUMER, INIMES
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     0.  -  INITIALISATIONS ET CONTROLES
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAIS2F_MT',0,ZHOOK_HANDLE)
      KREP=0
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
      IRANGC=FA%NUCADR(KRANG)
C**
C     1.  -  INITIALISATION DU TABLEAU FA%NSC2ALF
C-----------------------------------------------------------------------
C
C  Les valeurs de ce tableau representent les nb de pts
C le long de chaque "parallele" (ici, le nb de coeff spectraux
C pour un meme m (nb d'onde zonal), excepte le triangle et les axes non
C compactes). Il s'agit en effet de deguiser un champ spectral 
C Aladin en champ pts de grille (grille lat-lon) pour profiter
C du compactage, voire de la compression, GRIBEX.
C Le rangt des CSP est fait verticalement (par colonne de m=cst)
C et pour chaque couple (m,n) correspond 4 CSP.
C
      JMAX = (FA%NOZPAR(6,IRANGC)-FA%NOZPAR(5,IRANGC)+1)/4 -1
      DO JM=1,JMAX
        ILOW=2+2*JM+1
        IADD=4* MAX(FA%NSTROF(KRANG)+1-JM,1)
C
        FA%NSC2ALF(JM,KRANG)=FA%NOMPAR(ILOW+1,IRANGC)-
     S                        (FA%NOMPAR(ILOW,IRANGC)+IADD)+1
      ENDDO
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      LLFATA=LLMOER (KREP,KRANG)
C
      IF (FA%LFAMOP.OR.LLFATA) THEN
        INIMES=2
        CLNSPR='FAIS2F'
        INUMER=FA%JPNIIL
C
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4)')
     S     KREP, KRANG
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLNSPR,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAIS2F_MT',1,ZHOOK_HANDLE)
      END

