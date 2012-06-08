C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAIFLA_MT (FA, KRANG)
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme est charge des Initialisations des
C     tableaux FLAp1d., utilises pour aplatir le spectre des champs
C     d'un fichier avant le compactage (coefficients spectraux seulement).
C**
C
C
C**
C     ARGUMENTS :    KRANG  (Entree) ==> Rang de l'unite logique
C
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KRANG
C
      INTEGER J, IRANGC, IPUILA, ITRONC
#include "facom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAIFLA_MT',0,ZHOOK_HANDLE)
      IF (FA%FAIFLA_LLPREA) THEN
C
C
C
C         C'EST LE PREMIER APPEL au sous-programme,
C         on commence par les allocations des differents tableaux
C
        ALLOCATE (FA%FLAP1D (FA%JPXTRO,FA%JPNXFA))
        ALLOCATE (FA%FLAP1DA (FA%JPXTRO*FA%JPXTRO,FA%JPNXFA))
        FA%FAIFLA_LLPREA=.FALSE.
      ENDIF
C
      IPUILA = FA%NPUFLA(KRANG)
      IRANGC = FA%NUCADR(KRANG)
      ITRONC = FA%MTRONC(IRANGC)
      IF (IPUILA.GT.0) THEN
C
        DO J=1,ITRONC
          FA%FLAP1D(J,KRANG)=FA%XLAP1D(J,0)**IPUILA
        ENDDO
        DO J=1,FA%JPXTRO*FA%JPXTRO
          FA%FLAP1DA(J,KRANG)=FA%XLAP1DA(J,0)**IPUILA
        ENDDO
C
      ELSEIF (IPUILA.LT.0) THEN
C
        DO J=1,ITRONC
          FA%FLAP1D(J,KRANG)=FA%XLAP1D(J,1)**(-IPUILA)
        ENDDO
        DO J=1,FA%JPXTRO*FA%JPXTRO
          FA%FLAP1DA(J,KRANG)=FA%XLAP1DA(J,1)**(-IPUILA)
        ENDDO
C
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAIFLA_MT',1,ZHOOK_HANDLE)
      END

