! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FAIFLA_MT_BB          &
&                     (FA, KRANG)
USE FA_MOD, ONLY : FA_COM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!        Ce sous-programme est charge des Initialisations des
!     tableaux FLAp1d., utilises pour aplatir le spectre des champs
!     d'un fichier avant le compactage (coefficients spectraux seulement).
!**
!
!
!**
!     ARGUMENTS :    KRANG  (Entree) ==> Rang de l'unite logique
!
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KRANG
!
INTEGER (KIND=JPLIKB) J, IRANGC, IPUILA, ITRONC

!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FAIFLA_MT',0,ZHOOK_HANDLE)
!
IPUILA = FA%FICHIER(KRANG)%NPUFLA
IRANGC = FA%FICHIER(KRANG)%NUCADR
ITRONC = FA%CADRE(IRANGC)%MTRONC
!
IF (.NOT. ASSOCIATED (FA%FICHIER(KRANG)%FLAP1D)) &
& ALLOCATE (FA%FICHIER(KRANG)%FLAP1D (ITRONC))
IF (.NOT. ASSOCIATED (FA%FICHIER(KRANG)%FLAP1DA)) &
& ALLOCATE (FA%FICHIER(KRANG)%FLAP1DA (FA%JPXTRO*FA%JPXTRO))

!
IF (IPUILA.GT.0) THEN
!
  DO J=1,ITRONC
    FA%FICHIER(KRANG)%FLAP1D(J)=FA%XLAP1D(J,0)**IPUILA
  ENDDO
  DO J=1,FA%JPXTRO*FA%JPXTRO
    FA%FICHIER(KRANG)%FLAP1DA(J)=FA%XLAP1DA(J,0)**IPUILA
  ENDDO
!
ELSEIF (IPUILA.LT.0) THEN
!
  DO J=1,ITRONC
    FA%FICHIER(KRANG)%FLAP1D(J)=FA%XLAP1D(J,1)**(-IPUILA)
  ENDDO
  DO J=1,FA%JPXTRO*FA%JPXTRO
    FA%FICHIER(KRANG)%FLAP1DA(J)=FA%XLAP1DA(J,1)**(-IPUILA)
  ENDDO
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FAIFLA_MT',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FAIFLA_BB          &
&           (KRANG)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KRANG                                  ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAIFLA_MT_BB          &
&           (FA, KRANG)

END SUBROUTINE

SUBROUTINE FAIFLA_MM          &
&           (KRANG)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAIFLA_MT_MM          &
&           (FA, KRANG)

END SUBROUTINE

SUBROUTINE FAIFLA_MT_MM          &
&           (FA, KRANG)
USE FA_MOD, ONLY : FA_COM,              &
&                   FA_COM_DEFAULT_INIT, &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
! Local integers
INTEGER (KIND=JPLIKB)  IRANG                                  ! IN   
! Convert arguments

IRANG      = INT (     KRANG, JPLIKB)

CALL FAIFLA_MT_BB          &
&           (FA, IRANG)


END SUBROUTINE

!INTF KRANG         IN    
