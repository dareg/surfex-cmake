      SUBROUTINE PACKGB ( PFDATA, KPACKD, PREFER, PSCALE, KLENG )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C
C********************************************************************
C*
C*    NAME      : PACKGB
C*
C*    FUNCTION  : COMPUTES INDIVIDUAL "PACKED" VALUES (FIELD FOR GRIB),
C*                THE RESULT CONSISTS OF ONE DATA READY TO PACK
C*                WITHIN A BIT STRING PER COMPUTER WORD.
C*                   This subroutine has been designed to avoid explicit
C*                mixed use of REAL and INTEGER type values within the
C*                dummy-argument array PFDATA of CODEGA, this explicit
C*                use leading to non-standard code. The following code
C*                enables use of the same actual argument for the 2
C*                dummy-argument arrays.
C*
C*    INPUT     : PFDATA - FLOATING-POINT VALUES TO BE PACKED.
C*                PREFER - REFERENCE VALUE OF THE FIELD: SHOULD BE THE
C*                         MINIMUM VALUE, OR AN "UNDER-APPROXIMATION"
C*                         OF THE MINIMUM VALUE).
C*                PSCALE - SCALING FACTOR.
C*                KLENG  - NUMBER OF VALUES TO TREAT.
C*
C*    OUTPUT    : KPACKD - (POSITIVE) INTEGER VALUES "READY TO PACK"
C*
C*    AUTHOR    : J.CLOCHARD, FRENCH WEATHER SERVICE, 01/03/90.
C*
C********************************************************************
C*
#include "precision.h"
C
C     JP_STRIDE= pas permettant la correspondance entre les elements
C                d'un tableau de reels (KIND=JPDBLR): PFDATA(J) et
C                les elements d'un tableau d'entiers KPACKD(JP_STRIDE*J)
C                defini comme un tableau d'entiers representes sur
C                autant de bits que les reels.
C
      INTEGER, PARAMETER :: JP_STRIDE = JPDBLR / JP_SIMPLE_ENTIER
!
! If integers are on 32 bits, don't be afraid by the number of bits of the 
! real argument which is real and on 64 bits ... it's a trick : we pack 
! PFDATA to KPACKD in the same area, but with half bit and we are kind 
! enough to have a stride of 2 to access KPACKD ... ! (see also unpagb.F)
!
      INTEGER KLENG
C
      INTEGER KPACKD (JP_STRIDE*KLENG)
C
      REAL (KIND=JPDBLR) PREFER, PSCALE
C
      REAL (KIND=JPDBLR) PFDATA (KLENG)
C
      INTEGER J, II
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
C**
C     1.  -  STRAIGHT-FORWARD COMPUTING.
C
      IF (LHOOK) CALL DR_HOOK('PACKGB',0,ZHOOK_HANDLE)
!$OMP PARALLEL DO PRIVATE(J,II) SCHEDULE(STATIC,4096)
      DO 101 J=1,KLENG
#if defined(LITTLE_ENDIAN) || defined(LITTLE)
        II=JP_STRIDE*J -1
#else
        II=JP_STRIDE*J
#endif
        KPACKD(II)=NINT ( ( PFDATA(J) - PREFER ) * PSCALE )
  101 CONTINUE
!$OMP END PARALLEL DO
C
      IF (LHOOK) CALL DR_HOOK('PACKGB',1,ZHOOK_HANDLE)
      END
