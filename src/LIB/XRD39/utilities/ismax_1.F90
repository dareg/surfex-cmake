! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
      FUNCTION ISMAX_1_BB &
     &        ( KNELEM, PCHAMP) RESULT (IISMAX_1)
      USE LFI_PRECISION, ONLY : JPLIKB, JPDBLR
      IMPLICIT NONE
 
!    Search for position of maximum value in array PCHAMP.
!    Simplified version (good for vector computers) assuming stride=1
!     Original  F. Vana - ONPP/CHMI - 23-Mar-2010

      INTEGER (KIND=JPLIKB) KNELEM, IISMAX_1
      REAL (KIND=JPDBLR) PCHAMP (*)

      IF (KNELEM <= 0) THEN
        IISMAX_1=0
      ELSE
        IISMAX_1=INT (MAXLOC(ARRAY=PCHAMP(1:KNELEM) ,DIM=1), JPLIKB)
      ENDIF
      RETURN
      END
