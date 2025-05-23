!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE ERROR_WRITE_SURF_TXT(HREC,KRESP)
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
 CHARACTER(LEN=*), INTENT(IN) :: HREC
INTEGER,          INTENT(OUT):: KRESP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ERROR_WRITE_SURF_TXT',0,ZHOOK_HANDLE)
KRESP = -1
!
WRITE(*,*) ' '
WRITE(*,*) 'WARNING'
WRITE(*,*) '-------'
WRITE(*,*) ' '
WRITE(*,*) 'error when writing article', HREC
WRITE(*,*) "default value may be used; who knows?"
WRITE(*,*) ' '
!
IF (LHOOK) CALL DR_HOOK('ERROR_WRITE_SURF_TXT',1,ZHOOK_HANDLE)
!
END SUBROUTINE ERROR_WRITE_SURF_TXT
