! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FA_LIMITS_MT_BB                              &
&                       (FA, KPXPAH,KPXIND,KPXGEO,KPXNIV)
USE FA_MOD, ONLY : FA_COM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPLIKB
IMPLICIT NONE

! -----
! Purpose:
!   Returns limits necessary for proper dimensioning of FA fields
!   (arrays entering subroutine FACIES). It was created because
!   of two reasons:
!
!   1) Subroutine FALIMU does not provide all needed limits.
!
!   2) Direct use of parameters JPX* from include "facomp.h" would
!   be more convenient (no need to define interface arrays as
!   allocatable). But since this include contains CPP directives,
!   its use in external utilities might give wrong numbers with
!   improper compile options. Moreover, include is written
!   in fixed form and cannot be inserted into free form source.
!
! Created:
!   25-Jun-2010, J. Masek
!
! Modified:
!
! -----

TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) :: KPXPAH,KPXIND,KPXGEO,KPXNIV

KPXPAH=FA%JPXPAH
KPXIND=FA%JPXIND
KPXGEO=FA%JPXGEO
KPXNIV=FA%JPXNIV

END SUBROUTINE 



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FA_LIMITS_BB                    &
&           (KPXPAH, KPXIND, KPXGEO, KPXNIV)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPLIKB
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KPXPAH                                 !   OUT
INTEGER (KIND=JPLIKB)  KPXIND                                 !   OUT
INTEGER (KIND=JPLIKB)  KPXGEO                                 !   OUT
INTEGER (KIND=JPLIKB)  KPXNIV                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FA_LIMITS_MT_BB                           &
&           (FA, KPXPAH, KPXIND, KPXGEO, KPXNIV)

END SUBROUTINE

SUBROUTINE FA_LIMITS_MM                    &
&           (KPXPAH, KPXIND, KPXGEO, KPXNIV)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KPXPAH                                 !   OUT
INTEGER (KIND=JPLIKM)  KPXIND                                 !   OUT
INTEGER (KIND=JPLIKM)  KPXGEO                                 !   OUT
INTEGER (KIND=JPLIKM)  KPXNIV                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FA_LIMITS_MT_MM                           &
&           (FA, KPXPAH, KPXIND, KPXGEO, KPXNIV)

END SUBROUTINE

SUBROUTINE FA_LIMITS_MT_MM                     &
&           (FA, KPXPAH, KPXIND, KPXGEO, KPXNIV)
USE FA_MOD, ONLY : FA_COM,              &
&                   FA_COM_DEFAULT_INIT, &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPLIKM, JPLIKB
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KPXPAH                                 !   OUT
INTEGER (KIND=JPLIKM)  KPXIND                                 !   OUT
INTEGER (KIND=JPLIKM)  KPXGEO                                 !   OUT
INTEGER (KIND=JPLIKM)  KPXNIV                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IPXPAH                                 !   OUT
INTEGER (KIND=JPLIKB)  IPXIND                                 !   OUT
INTEGER (KIND=JPLIKB)  IPXGEO                                 !   OUT
INTEGER (KIND=JPLIKB)  IPXNIV                                 !   OUT
! Convert arguments


CALL FA_LIMITS_MT_BB                           &
&           (FA, IPXPAH, IPXIND, IPXGEO, IPXNIV)

KPXPAH     = INT (    IPXPAH, JPLIKM)
KPXIND     = INT (    IPXIND, JPLIKM)
KPXGEO     = INT (    IPXGEO, JPLIKM)
KPXNIV     = INT (    IPXNIV, JPLIKM)

END SUBROUTINE

!INTF KPXPAH          OUT 
!INTF KPXIND          OUT 
!INTF KPXGEO          OUT 
!INTF KPXNIV          OUT 
