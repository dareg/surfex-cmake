! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"


#undef JLIK
#undef _ELLIPS_
#define JLIK JPLIKB
#define _ELLIPS_  ELLIPS_BB
#include "ellips.h"

#undef JLIK
#undef _ELLIPS_
#define JLIK JPLIKM
#define _ELLIPS_  ELLIPS_MM
#include "ellips.h"
