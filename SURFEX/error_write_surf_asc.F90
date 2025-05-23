!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE ERROR_WRITE_SURF_ASC(HREC,KRESP)
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
USE MODD_IO_SURF_ASC, ONLY : NLUOUT
!
IMPLICIT NONE
!
 CHARACTER(LEN=*), INTENT(IN) :: HREC
INTEGER,          INTENT(OUT):: KRESP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('ERROR_WRITE_SURF_ASC',0,ZHOOK_HANDLE)
KRESP = -1
!
WRITE(NLUOUT,*) ' '
WRITE(NLUOUT,*) 'WARNING'
WRITE(NLUOUT,*) '-------'
WRITE(NLUOUT,*) ' '
WRITE(NLUOUT,*) 'error when writing article', HREC
WRITE(NLUOUT,*) "default value may be used; who knows?"
WRITE(NLUOUT,*) ' '
IF (LHOOK) CALL DR_HOOK('ERROR_WRITE_SURF_ASC',1,ZHOOK_HANDLE)
!
END SUBROUTINE ERROR_WRITE_SURF_ASC
