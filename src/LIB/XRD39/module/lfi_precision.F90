MODULE LFI_PRECISION

#include "lfisuffix.h"

!
!----- DESCRIPTION DES "PARAMETER" DU LOGICIEL DE FICHIERS INDEXES -----
!
!     JPDBLE= PRECISION UTILISE POUR LES ENTIERS:
!     * SI JPDBLE=8 les INTEGER (KIND=JPDBLE) seront a priori en 64 BITS
!     * SI JPDBLE=4 les INTEGER (KIND=JPDBLE) seront a priori en 32 BITS
!
!     JPDBLR= PRECISION UTILISE POUR LES FLOTTANTS (REELS):
!     * SI JPDBLR=8 les REAL (KIND=JPDBLR) seront a priori en 64 BITS
!     * SI JPDBLR=4 les REAL (KIND=JPDBLR) seront a priori en 32 BITS
!
!     (les conventions peuvent dependre de la plate-forme consideree)
!
!     JP_SIMPLE_ENTIER= sous-type entier permettant de representer
!                       l'intervalle +/-(10**9 - 1)
!
      INTEGER, PARAMETER :: JPDBLE = 8, JPDBLR = 8
      INTEGER, PARAMETER :: JP_SIMPLE_ENTIER = SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: JPLIKB = _JPLIKB_, JPLIKM = _JPLIKM_

END MODULE LFI_PRECISION
