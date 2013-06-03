! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
SUBROUTINE FADECOGA                                           &
&           (PFDATA, KLENF, KBITS, KNBIT, KB1PAR, KB2PAR,      &
&           PVERT, KLENV, KGRIB, KLENG, KWORD, KJLENV, KJLENF, &
&           KCPACK, KSCALP, KERR, PMIN, PMAX, LDARPE)
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KLENF                                  ! IN   
REAL (KIND=JPDBLR)     PFDATA     (KLENF)                     !   OUT
INTEGER (KIND=JPLIKB)  KBITS                                  !   OUT
INTEGER (KIND=JPLIKB)  KNBIT                                  ! IN   
INTEGER (KIND=JPLIKB)  KB1PAR     (19)                        !   OUT
INTEGER (KIND=JPLIKB)  KB2PAR     (17)                        !   OUT
INTEGER (KIND=JPLIKB)  KLENV                                  ! IN   
REAL (KIND=JPDBLR)     PVERT      (KLENV)                     ! IN   
INTEGER (KIND=JPLIKB)  KLENG                                  ! IN   
INTEGER (KIND=JPDBLE)  KGRIB      (KLENG)                     ! IN   
INTEGER (KIND=JPLIKB)  KWORD                                  !   OUT
INTEGER (KIND=JPLIKB)  KJLENV                                 !   OUT
INTEGER (KIND=JPLIKB)  KJLENF                                 !   OUT
INTEGER (KIND=JPLIKB)  KCPACK                                 !   OUT
INTEGER (KIND=JPLIKB)  KSCALP                                 !   OUT
INTEGER (KIND=JPLIKB)  KERR                                   !   OUT
REAL (KIND=JPDBLR)     PMIN                                   !   OUT
REAL (KIND=JPDBLR)     PMAX                                   !   OUT
LOGICAL                LDARPE                                 ! IN   
! Local integers
INTEGER (KIND=JPLIKM)  ILENF                                  ! IN   
INTEGER (KIND=JPLIKM)  IBITS                                  !   OUT
INTEGER (KIND=JPLIKM)  INBIT                                  ! IN   
INTEGER (KIND=JPLIKM)  IB1PAR     (19)                        !   OUT
INTEGER (KIND=JPLIKM)  IB2PAR     (17)                        !   OUT
INTEGER (KIND=JPLIKM)  ILENV                                  ! IN   
INTEGER (KIND=JPLIKM)  ILENG                                  ! IN   
INTEGER (KIND=JPLIKM)  IWORD                                  !   OUT
INTEGER (KIND=JPLIKM)  IJLENV                                 !   OUT
INTEGER (KIND=JPLIKM)  IJLENF                                 !   OUT
INTEGER (KIND=JPLIKM)  ICPACK                                 !   OUT
INTEGER (KIND=JPLIKM)  ISCALP                                 !   OUT
INTEGER (KIND=JPLIKM)  IERR                                   !   OUT
! Convert arguments

ILENF      = INT (     KLENF, JPLIKM)
INBIT      = INT (     KNBIT, JPLIKM)
ILENV      = INT (     KLENV, JPLIKM)
ILENG      = INT (     KLENG, JPLIKM)

CALL DECOGA                                                   &
&           (PFDATA, ILENF, IBITS, INBIT, IB1PAR, IB2PAR,      &
&           PVERT, ILENV, KGRIB, ILENG, IWORD, IJLENV, IJLENF, &
&           ICPACK, ISCALP, IERR, PMIN, PMAX, LDARPE)

KBITS      = INT (     IBITS, JPLIKB)
KB1PAR     = INT (    IB1PAR, JPLIKB)
KB2PAR     = INT (    IB2PAR, JPLIKB)
KWORD      = INT (     IWORD, JPLIKB)
KJLENV     = INT (    IJLENV, JPLIKB)
KJLENF     = INT (    IJLENF, JPLIKB)
KCPACK     = INT (    ICPACK, JPLIKB)
KSCALP     = INT (    ISCALP, JPLIKB)
KERR       = INT (      IERR, JPLIKB)

END SUBROUTINE
