C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAIXLA_MT (FA)
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme est charge des Initialisations des
C     tableaux XLAp.d., utilises pour aplatir le spectre des champs
C     avant le compactage (coefficients spectraux seulement).
C**
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER J, JN, IDEBUT, IFIN, INDM, INDN, JJPUIS
#include "facom_mt.h"
C
C
C        On commence par les allocations des differents tableaux
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAIXLA_MT',0,ZHOOK_HANDLE)
      ALLOCATE (FA%XLAP1D (FA%JPXTRO,0:1))
      ALLOCATE (FA%XLAP1DA (FA%JPXTRO*FA%JPXTRO,0:1))
      ALLOCATE (FA%XLAP2D (2:FA%JPXCSP,FA%JPUILA,0:1))
      ALLOCATE (FA%XLAP2DA (FA%JPXTRO*FA%JPXTRO,FA%JPUILA,0:1))
C
      DO JN=1,FA%JPXTRO
        FA%XLAP1D(JN,0)=REAL (JN*(JN+1))
        FA%XLAP1D(JN,1)=1./FA%XLAP1D(JN,0)
      ENDDO
C
      DO JN=1,FA%JPXTRO*FA%JPXTRO
        INDN=1+(JN-1)/FA%JPXTRO
        INDM=JN-(INDN-1)*FA%JPXTRO
        FA%XLAP1DA(JN,0)=REAL (INDN**2+INDM**2)
        FA%XLAP1DA(JN,1)=1./FA%XLAP1DA(JN,0)
      ENDDO
C
      DO JN=1,FA%JPXTRO
        IDEBUT=JN**2+1
        IFIN=(1+JN)**2
C
        DO JJPUIS=1,FA%JPUILA
C
        DO J=IDEBUT,IFIN
          FA%XLAP2D(J,JJPUIS,0)=FA%XLAP1D(JN,0)**JJPUIS
          FA%XLAP2D(J,JJPUIS,1)=1./FA%XLAP2D(J,JJPUIS,0)
        ENDDO
C        
        ENDDO
C
      ENDDO
C
      DO JJPUIS=1,FA%JPUILA
C
      DO JN=1,FA%JPXTRO**2 
        FA%XLAP2DA(JN,JJPUIS,0)=FA%XLAP1DA(JN,0)**JJPUIS
        FA%XLAP2DA(JN,JJPUIS,1)=1./FA%XLAP2DA(JN,JJPUIS,0)
      ENDDO
C
      ENDDO
C
C
      IF (LHOOK) CALL DR_HOOK('FAIXLA_MT',1,ZHOOK_HANDLE)
      END

