      SUBROUTINE UNPAGB ( KPDATA, PFDATA, PMIN, PMAX, KBITS, PSCALE,
     S                    KLENG, LDARPE )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C
C********************************************************************
C*
C*    NAME      : UNPAGB
C*
C*    FUNCTION  : COMPUTES INDIVIDUAL "UNPACKED" VALUES (FIELD FROM GRIB
C*                ), THE INPUT CONSISTING OF ONE DATA JUST UNPACKED
C*                FROM A BIT STRING PER COMPUTER WORD.
C*                   This subroutine has been designed to avoid explicit
C*                mixed use of REAL and INTEGER type values within the
C*                dummy-argument array PFDATA of DECOGA, this explicit
C*                use leading to non-standard code. The following code
C*                enables use of the same actual argument for the 2
C*                dummy-argument arrays.
C*
C*    INPUT     : KPDATA - (POSITIVE) INTEGER VALUES "JUST UNPACKED"
C*                PMIN   - MINIMUM VALUE, OR AN "UNDER-APPROXIMATION"
C*                         OF THE MINIMUM VALUE).
C*                PMAX   - MAXIMUM VALUE, OR A "OVER-APPROXIMATION"
C*                         OF THE MAXIMUM VALUE).
C*                KBITS  - NUMBER OF BITS PER CODED VALUE.
C*                PSCALE - SCALE FACTOR TO APPLY.
C*                KLENG  - NUMBER OF VALUES TO BE TREATED.
C*                LDARPE  - .TRUE., modifications for ARPEGE coding
C*                                  have been included when coding data;
C*                          .FALSE., no such modifications.
C*
C*    PMAX and KBITS are used only if LDARPE is .TRUE. .
C*
C*
C*    OUTPUT    : PFDATA - FLOATING-POINT VALUES.
C*
C*    AUTHOR    : J.CLOCHARD, FRENCH WEATHER SERVICE, 01/03/90.
C*
C********************************************************************
C*
#include "precision.h"
C
C     JP_STRIDE= pas permettant la correspondance entre les elements
C                d'un tableau de reels (KIND=JPDBLR): PFDATA(J) et
C                les elements d'un tableau d'entiers KPDATA(JP_STRIDE*J)
C                defini comme un tableau d'entiers representes sur
C                autant de bits que les reels.
C
      INTEGER, PARAMETER :: JP_STRIDE = JPDBLR / JP_SIMPLE_ENTIER
!
! If integers are on 32 bits, don't be afraid by the number of bits of the 
! real argument which is real and on 64 bits ... it's a trick : we unpack 
! KPDATA to PFDATA in the same area, but we are kind enought to have a
! stride of 2 to access KPDATA ... ! (see also packgb.F)
!
      INTEGER KLENG, KBITS
C
      INTEGER KPDATA (JP_STRIDE*KLENG)
C
      REAL (KIND=JPDBLR) PMIN, PMAX, PSCALE
C
      REAL (KIND=JPDBLR) PFDATA (KLENG)
C
      INTEGER J, II, IAUXI1, IAUXI2
C
      LOGICAL LDARPE
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('UNPAGB',0,ZHOOK_HANDLE)
      IF (LDARPE) THEN
C**
C     1.  -  DIRECT COMPUTING WITH 2 CASES, WHICH ENABLES PERFECT
C            RESPECT OF MINIMUM AND MAXIMUM PROVIDED THAT PMIN AND PMAX
C            ARE THESE VALUES.
C
        IAUXI1=2**(KBITS-1)
        IAUXI2=2*IAUXI1-1
C
C            Here, PSCALE is (PMAX-PMIN)/(FLOAT(IAUXI2).
C
!$OMP PARALLEL DO PRIVATE(J,II) SCHEDULE(STATIC,4096)
        DO 101 J=KLENG,1,-1
C
#if defined(LITTLE_ENDIAN) || defined(LITTLE)
          II=JP_STRIDE*J -1
#else
          II=JP_STRIDE*J
#endif
          IF (KPDATA(II).LT.IAUXI1) THEN
            PFDATA(J)=PMIN+PSCALE*REAL (KPDATA(II),JPDBLR)
          ELSE
            PFDATA(J)=PMAX-PSCALE*REAL (IAUXI2-KPDATA(II),JPDBLR)
          ENDIF
C
101     CONTINUE
!$OMP END PARALLEL DO
C
      ELSE
C**
C     2.  -  DIRECT COMPUTING, WHICH ENABLES PERFECT RESPECT
C            OF MINIMUM PROVIDED THAT PMIN IS THIS VALUE.
C       (in standard GRIB, there is no estimation of the field maximum,
C        just a "gross" over-approximation can be given)
C
C
!$OMP PARALLEL DO PRIVATE(J,II) SCHEDULE(STATIC,4096)
        DO 201 J=KLENG,1,-1
#if defined(LITTLE_ENDIAN) || defined(LITTLE)
          II=JP_STRIDE*J -1
#else
          II=JP_STRIDE*J
#endif
          PFDATA(J)=PMIN+PSCALE*REAL (KPDATA(II),JPDBLR)
201     CONTINUE
!$OMP END PARALLEL DO
C
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('UNPAGB',1,ZHOOK_HANDLE)
      END
