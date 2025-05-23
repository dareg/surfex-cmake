!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE SURFEX_DEALLO(YDSURFEX)
!
USE MODD_SURFEX_n, ONLY : SURFEX_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
USE MODI_DEALLOC_SURF_ATM_n
!
IMPLICIT NONE
!
TYPE (SURFEX_t), INTENT(INOUT) :: YDSURFEX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("SURFEX_DEALLO",0,ZHOOK_HANDLE)
!
CALL DEALLOC_SURF_ATM_n(YDSURFEX)
!
IF (LHOOK) CALL DR_HOOK("SURFEX_DEALLO",1,ZHOOK_HANDLE)
!
END SUBROUTINE SURFEX_DEALLO
