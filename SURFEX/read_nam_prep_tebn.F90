!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE READ_NAM_PREP_TEB_n(HPROGRAM)
!     #######################################################
!
!---------------------------------------
!
USE MODN_PREP_TEB_SNOW
USE MODN_PREP_TEB
!
USE MODD_SURF_PAR, ONLY : XUNDEF, NUNDEF
!
USE MODI_DEFAULT_PREP_TEB
USE MODI_READ_PREP_TEB_SNOW
!
USE MODI_TEST_NAM_VAR_SURF
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODE_POS_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
!
INTEGER :: ILUNAM         ! logical unit of namelist file
INTEGER :: ILUOUT
LOGICAL :: GFOUND         ! Return code when searching namelist
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!---------------------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_NAM_PREP_TEB_N',0,ZHOOK_HANDLE)
NYEAR=NUNDEF
NMONTH=NUNDEF
NDAY=NUNDEF
XTIME=XUNDEF
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!
 CALL DEFAULT_PREP_TEB
!
 CALL OPEN_NAMELIST(HPROGRAM,ILUNAM)
 CALL POSNAM(ILUNAM,'NAM_PREP_TEB',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_PREP_TEB)
 CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
!
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CROAD_DIR',CROAD_DIR,'UNIF','ORIE')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CWALL_OPT',CWALL_OPT,'UNIF','TWO ')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CTYPE',   CTYPE,   '      ','GRIB  ','MESONH','ASCII ','LFI   ','FA    ')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CTYPEPGD',CTYPEPGD,'      ','GRIB  ','MESONH','ASCII ','LFI   ','FA    ')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CTYPE_WS',CTYPE_WS,'      ','GRIB  ','MESONH','ASCII ','LFI   ')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CTYPE_TS',CTYPE_TS,'      ','GRIB  ','MESONH','ASCII ','LFI   ')
!
 CALL READ_PREP_TEB_SNOW(HPROGRAM,CSNOW_ROOF,NSNOW_ROOF,CSNOW_ROAD,NSNOW_ROAD)
IF (LHOOK) CALL DR_HOOK('READ_NAM_PREP_TEB_N',1,ZHOOK_HANDLE)
!
!------------------------------------
!
END SUBROUTINE READ_NAM_PREP_TEB_n
