      INTEGER FUNCTION ISMAX_1 ( KNELEM, PCHAMP)

#include "precision.h"
 
c    Search for position of maximum value in array PCHAMP.
c    Simplified version (good for vector computers) assuming stride=1
c     Original  F. Vana - ONPP/CHMI - 23-Mar-2010

      INTEGER KNELEM
      REAL (KIND=JPDBLR) PCHAMP (*)

      IF (KNELEM <= 0) THEN
        ISMAX_1=0
      ELSE
        ISMAX_1=MAXLOC(ARRAY=PCHAMP(1:KNELEM) ,DIM=1)
      ENDIF
      RETURN
      END
