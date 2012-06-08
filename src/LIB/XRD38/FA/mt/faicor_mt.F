C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAICOR_MT (FA)
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C      Initialisation des tableaux permettant la correspondance
C      entre nom du champ dans FA et ses descripteurs GRIB.
C
C      **  FA%NBPARC ne doit pas depasser FA%JPXPAR !! **
C**
C
      TYPE(FA_COM) :: FA
      INTEGER J, INB, INBPAR
#include "facom_mt.h"
C**
C     1.  -  TRAITEMENT DES CHAMPS MONO NIVEAU
C-----------------------------------------------------------------------
C
C     1.1 -  Champs "profond"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAICOR_MT',0,ZHOOK_HANDLE)
      INBPAR = 0
      FA%NBPARC = INBPAR + 7
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'PROF'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'TEMPERATURE'
      FA%CISUFF(INBPAR+2 ) = 'RESERV.GLACE'
      FA%CISUFF(INBPAR+3 ) = 'RESERV.EAU'
      FA%CISUFF(INBPAR+4 ) = 'RUISSELLEMEN'
      FA%CISUFF(INBPAR+5 ) = 'FLUX DE GEL'
      FA%CISUFF(INBPAR+6 ) = 'PROP.RMAX.EA'
      FA%CISUFF(INBPAR+7 ) = 'XRUISSELLEME'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,  11, 111,  10,   0, 0 /)
      FA%NCODPA(INBPAR+2 ,1:6) = (/   1, 152, 112,   0, 250, 0 /)
      FA%NCODPA(INBPAR+3 ,1:6) = (/   1,  86, 112,   0, 250, 0 /)
      FA%NCODPA(INBPAR+4 ,1:6) = (/ 149,   2, 112,  10, 300, 4 /)
      FA%NCODPA(INBPAR+5 ,1:6) = (/ 149,  73, 112,  10, 300, 4 /)
      FA%NCODPA(INBPAR+6 ,1:6) = (/ 149,  77, 112,  10, 300, 0 /)
      FA%NCODPA(INBPAR+7 ,1:6) = (/ 149,  85, 112,  10, 300, 0 /)
      
C
C     1.2 -  Champs de surface
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 83
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'SURF'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'PRESSION'
      FA%CISUFF(INBPAR+2 ) = 'GEOPOTENTIEL'
      FA%CISUFF(INBPAR+3 ) = 'TEMPERATURE'
      FA%CISUFF(INBPAR+4 ) = 'FLU.MSUBL.NE'
      FA%CISUFF(INBPAR+5 ) = 'FLU.MEVAP.EA'
      FA%CISUFF(INBPAR+6 ) = 'PREC.EAU.GEC'
      FA%CISUFF(INBPAR+7 ) = 'PREC.EAU.CON'
      FA%CISUFF(INBPAR+8 ) = 'RESERV.NEIGE'
      FA%CISUFF(INBPAR+9 ) = 'NEBUL.TOTALE'
      FA%CISUFF(INBPAR+10) = 'NEBUL.CONVEC'
      FA%CISUFF(INBPAR+11) = 'NEBUL.BASSE'
      FA%CISUFF(INBPAR+12) = 'NEBUL.MOYENN'
      FA%CISUFF(INBPAR+13) = 'NEBUL.HAUTE'
      FA%CISUFF(INBPAR+14) = 'PREC.NEI.CON'
      FA%CISUFF(INBPAR+15) = 'PREC.NEI.GEC'
      FA%CISUFF(INBPAR+16) = 'IND.TERREMER'
      FA%CISUFF(INBPAR+17) = 'Z0.FOIS.G'
      FA%CISUFF(INBPAR+18) = 'ALBEDO'
      FA%CISUFF(INBPAR+19) = 'PROP.VEGETAT'
      FA%CISUFF(INBPAR+20) = 'FLU.RAY.SOLA'
      FA%CISUFF(INBPAR+21) = 'FLU.RAY.THER'
      FA%CISUFF(INBPAR+22) = 'FLU.LAT.MEVA'
      FA%CISUFF(INBPAR+23) = 'FLU.LAT.MSUB'
      FA%CISUFF(INBPAR+24) = 'FLU.CHA.SENS'
      FA%CISUFF(INBPAR+25) = 'ETA.GEOPOTEN'
      FA%CISUFF(INBPAR+26) = 'ETP.GEOPOTEN'
      FA%CISUFF(INBPAR+27) = 'TENS.DMOG.ZO'
      FA%CISUFF(INBPAR+28) = 'TENS.TURB.ZO'
      FA%CISUFF(INBPAR+29) = 'TENS.DMOG.ME'
      FA%CISUFF(INBPAR+30) = 'TENS.TURB.ME'
      FA%CISUFF(INBPAR+31) = 'RESERV.GLACE'
      FA%CISUFF(INBPAR+32) = 'RESERV.EAU'
      FA%CISUFF(INBPAR+33) = 'RAYT.LUNE.DE'
      FA%CISUFF(INBPAR+34) = 'CAPE.MOD.XFU'
      FA%CISUFF(INBPAR+35) = 'TOT.WAT.VAPO'
C Debut des champs pour Mocage :
CDP: champ TCO3 remplace par O3_TOT le 18 oct 2004 par coherence
C avec la BDAP qui venait de proceder au changement.
      FA%CISUFF(INBPAR+36) = 'O3_TOT'
      FA%CISUFF(INBPAR+37) = 'CH4_TOT'
      FA%CISUFF(INBPAR+38) = 'NO_TOT'
      FA%CISUFF(INBPAR+39) = 'NO2_TOT'
      FA%CISUFF(INBPAR+40) = 'CO_TOT'
C Fin des champs pour Mocage
      FA%CISUFF(INBPAR+41) = 'ALBEDO NEIGE'
      FA%CISUFF(INBPAR+42) = 'ALBEDO HISTO'
      FA%CISUFF(INBPAR+43) = 'ALBEDO.SOLNU'
      FA%CISUFF(INBPAR+44) = 'ALBEDO.VEG'
      FA%CISUFF(INBPAR+45) = 'EMISSIVITE'
      FA%CISUFF(INBPAR+46) = 'ET.GEOPOTENT'
      FA%CISUFF(INBPAR+47) = 'VAR.GEOP.ANI'
      FA%CISUFF(INBPAR+48) = 'VAR.GEOP.DIR'
      FA%CISUFF(INBPAR+49) = 'IND.VEG.DOMI'
      FA%CISUFF(INBPAR+50) = 'RESI.STO.MIN'
      FA%CISUFF(INBPAR+51) = 'PROP.ARGILE'
      FA%CISUFF(INBPAR+52) = 'PROP.SABLE'
      FA%CISUFF(INBPAR+53) = 'EPAIS.SOL'
      FA%CISUFF(INBPAR+54) = 'IND.FOLIAIRE'
      FA%CISUFF(INBPAR+55) = 'RES.EVAPOTRA'
      FA%CISUFF(INBPAR+56) = 'GZ0.THERM'
      FA%CISUFF(INBPAR+57) = 'RESERV.INTER'
      FA%CISUFF(INBPAR+58) = 'FL.Q TURBUL'
      FA%CISUFF(INBPAR+59) = 'FL.CT TURBUL'
      FA%CISUFF(INBPAR+60) = 'FONTE NEIGE'
      FA%CISUFF(INBPAR+61) = 'CHAL. DS SOL'
      FA%CISUFF(INBPAR+62) = 'EAU DANS SOL'
      FA%CISUFF(INBPAR+63) = 'RUISSELLEMEN'
      FA%CISUFF(INBPAR+64) = 'EVAPOTRANSPI'
      FA%CISUFF(INBPAR+65) = 'TRANSPIRATIO'
      FA%CISUFF(INBPAR+66) = 'RUISS. INTER'
      FA%CISUFF(INBPAR+67) = 'FLUX DE GEL'
      FA%CISUFF(INBPAR+68) = 'RAYT DIR SUR'
      FA%CISUFF(INBPAR+69) = 'RAYT DIFF DE'
      FA%CISUFF(INBPAR+70) = 'RAYT THER DE'
      FA%CISUFF(INBPAR+71) = 'PROP.RMAX.EA'
      FA%CISUFF(INBPAR+72) = 'DENSIT.NEIGE'
      FA%CISUFF(INBPAR+73) = 'CFU.Q.TURBUL'
      FA%CISUFF(INBPAR+74) = 'CFU.CT.TURBU'
      FA%CISUFF(INBPAR+75) = 'TIME.PREC.TO'
      FA%CISUFF(INBPAR+76) = 'XFONTE NEIGE'
      FA%CISUFF(INBPAR+77) = 'XCHAL. DS SO'
      FA%CISUFF(INBPAR+78) = 'XEAU DANS SO'
      FA%CISUFF(INBPAR+79) = 'XFLU.MEVAP.E'
      FA%CISUFF(INBPAR+80) = 'XRUISSELLEME'
      FA%CISUFF(INBPAR+81) = 'XEVAPOTRANSP'
      FA%CISUFF(INBPAR+82) = 'XTRANSPIRATI'
      FA%CISUFF(INBPAR+83) = 'XRUISS. INTE'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,   1,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+2 ,1:6) = (/   1,   8,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+3 ,1:6) = (/   1,  11,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+4 ,1:6) = (/ 149,  45,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+5 ,1:6) = (/ 149,  46,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+6 ,1:6) = (/   1,  62,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+7 ,1:6) = (/   1,  63,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+8 ,1:6) = (/   1,  65,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+9 ,1:6) = (/   1,  71,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+10,1:6) = (/   1,  72,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+11,1:6) = (/   1,  73,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+12,1:6) = (/   1,  74,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+13,1:6) = (/   1,  75,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+14,1:6) = (/   1,  78,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+15,1:6) = (/   1,  79,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+16,1:6) = (/   1,  81,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+17,1:6) = (/ 149,  47,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+18,1:6) = (/   1,  84,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+19,1:6) = (/   1,  87,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+20,1:6) = (/ 128, 111,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+21,1:6) = (/ 128, 112,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+22,1:6) = (/ 149,  52,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+23,1:6) = (/ 149,  53,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+24,1:6) = (/ 128, 122,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+25,1:6) = (/ 149,  54,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+26,1:6) = (/ 149,  55,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+27,1:6) = (/ 149,  56,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+28,1:6) = (/ 149,  57,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+29,1:6) = (/ 149,  58,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+30,1:6) = (/ 149,  59,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+31,1:6) = (/   1, 152,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+32,1:6) = (/   1,  86,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+33,1:6) = (/   1, 158,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+34,1:6) = (/   1, 160,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+35,1:6) = (/   1, 167,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+36,1:6) = (/   1,  10,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+37,1:6) = (/ 159,  21,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+38,1:6) = (/ 159,  22,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+39,1:6) = (/ 159,  23,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+40,1:6) = (/ 159,  24,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+41,1:6) = (/ 149,   5,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+42,1:6) = (/ 149,  48,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+43,1:6) = (/ 149,  49,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+44,1:6) = (/ 149,  50,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+45,1:6) = (/ 149,  60,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+46,1:6) = (/ 149,  93,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+47,1:6) = (/ 128, 161,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+48,1:6) = (/ 128, 162,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+49,1:6) = (/ 149,  61,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+50,1:6) = (/ 149,  62,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+51,1:6) = (/ 149,  63,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+52,1:6) = (/ 149,  64,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+53,1:6) = (/ 149,  51,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+54,1:6) = (/ 149,  65,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+55,1:6) = (/ 149,  66,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+56,1:6) = (/ 149,  67,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+57,1:6) = (/ 149,  68,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+58,1:6) = (/ 149,  69,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+59,1:6) = (/ 149,  70,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+60,1:6) = (/   1,  99,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+61,1:6) = (/ 149,   8,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+62,1:6) = (/   1, 150,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+63,1:6) = (/ 149,   2,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+64,1:6) = (/   1, 168,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+65,1:6) = (/ 149,  71,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+66,1:6) = (/ 149,  72,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+67,1:6) = (/ 149,  73,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+68,1:6) = (/ 149,  74,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+69,1:6) = (/ 149,  75,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+70,1:6) = (/ 149,  76,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+71,1:6) = (/ 149,  77,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+72,1:6) = (/ 149,   3,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+73,1:6) = (/ 149,  78,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+74,1:6) = (/ 149,  79,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+75,1:6) = (/ 149,  80,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+76,1:6) = (/ 149,  81,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+77,1:6) = (/ 149,  82,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+78,1:6) = (/ 149,  83,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+79,1:6) = (/ 149,  84,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+80,1:6) = (/ 149,  85,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+81,1:6) = (/ 149,  86,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+82,1:6) = (/ 149,  87,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+83,1:6) = (/ 149,  88,   1,   0,   0, 0 /)
C
C  Geopotentiel de surface en spectral (relief)
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 1
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'SPECSURF'
      ENDDO
      FA%CISUFF(INBPAR+1) = 'GEOPOTEN'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,   8,   1,   0,   0, 0 /)
C
C     1.3 -  Champs au sommet
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 2
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'SOMM'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'FLU.RAY.SOLA'
      FA%CISUFF(INBPAR+2 ) = 'FLU.RAY.THER'
      FA%NCODPA(INBPAR+1 ,1:6) = (/ 128, 113,   8,   0,   0, 4 /)
      FA%NCODPA(INBPAR+2 ,1:6) = (/ 128, 114,   8,   0,   0, 4 /)
C
C     1.4 -  Champs de la Couche Limite de Surface
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 11
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'CLS'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'TEMPERATURE'
      FA%CISUFF(INBPAR+2 ) = 'MAXI.TEMPERAT'
      FA%CISUFF(INBPAR+3 ) = 'MINI.TEMPERAT'
      FA%CISUFF(INBPAR+4 ) = 'VENT.ZONAL'
      FA%CISUFF(INBPAR+5 ) = 'VENT.MERIDIEN'
      FA%CISUFF(INBPAR+6 ) = 'HUMI.SPECIFIQ'
      FA%CISUFF(INBPAR+7 ) = 'HUMI.RELATIVE'
      FA%CISUFF(INBPAR+8 ) = 'MAXI.HUMI.REL'
      FA%CISUFF(INBPAR+9 ) = 'MINI.HUMI.REL'
      FA%CISUFF(INBPAR+10) = 'U.RAF.MOD.XFU'
      FA%CISUFF(INBPAR+11) = 'V.RAF.MOD.XFU'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,  11, 105,   2,   0, 0 /)
      FA%NCODPA(INBPAR+2 ,1:6) = (/   1,  15, 105,   2,   0, 2 /)
      FA%NCODPA(INBPAR+3 ,1:6) = (/   1,  16, 105,   2,   0, 2 /)
      FA%NCODPA(INBPAR+4 ,1:6) = (/   1,  33, 105,  10,   0, 0 /)
      FA%NCODPA(INBPAR+5 ,1:6) = (/   1,  34, 105,  10,   0, 0 /)
      FA%NCODPA(INBPAR+6 ,1:6) = (/   1,  51, 105,   2,   0, 0 /)
      FA%NCODPA(INBPAR+7 ,1:6) = (/   1,  52, 105,   2,   0, 0 /)
      FA%NCODPA(INBPAR+8 ,1:6) = (/ 149,  90, 105,   2,   0, 2 /)
      FA%NCODPA(INBPAR+9 ,1:6) = (/ 149,  91, 105,   2,   0, 2 /)
      FA%NCODPA(INBPAR+10,1:6) = (/   1, 163, 105,  10,   0, 0 /)
      FA%NCODPA(INBPAR+11,1:6) = (/   1, 164, 105,  10,   0, 0 /)
C
C     1.5 -  Champs au niveau tropopause
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 2
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'ICAO'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'TROP.PRESSUR'
      FA%CISUFF(INBPAR+2 ) = 'TROP.TEMPERA'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,   1,   7,   0,   0, 0 /)
      FA%NCODPA(INBPAR+2 ,1:6) = (/   1,  11,   7,   0,   0, 0 /)
C
C     1.6 -  Champs au niveau jet
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 3
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'JET'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'PRESSURE'
      FA%CISUFF(INBPAR+2 ) = 'VENT_ZONAL'
      FA%CISUFF(INBPAR+3 ) = 'VENT_MERIDIEN'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,   1,   6,   0,   0, 0 /)
      FA%NCODPA(INBPAR+2 ,1:6) = (/   1,  33,   6,   0,   0, 0 /)
      FA%NCODPA(INBPAR+3 ,1:6) = (/   1,  34,   6,   0,   0, 0 /)
C
C     1.7 -  Champs pour la Couche Limite Planetaire
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 2
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'CLP'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'MHAUT.MOD.XFU'
      FA%CISUFF(INBPAR+2 ) = 'MOCON.MOD.XFU'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1, 165,   1,   0,   0, 0 /)
      FA%NCODPA(INBPAR+2 ,1:6) = (/   1, 166,   1,   0,   0, 0 /)
C
C     1.8 -  Champs au niveau mer
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 1
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'MSL'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'PRESSURE'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,   2, 102,   0,   0, 0 /)
C
C     1.9 - Relief
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 1
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'INT'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'SURFGEOPOTENT'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,   8,   1,   0,   0, 0 /)
C
C     1.10- ATMO
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 4
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'ATMO'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'NEBUL.CONVEC'
      FA%CISUFF(INBPAR+2 ) = 'NEBUL.BASSE'
      FA%CISUFF(INBPAR+3 ) = 'NEBUL.MOYENN'
      FA%CISUFF(INBPAR+4 ) = 'NEBUL.HAUTE'
      FA%NCODPA(INBPAR+1 ,1:6) = (/ 128,  72,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+2 ,1:6) = (/ 128,  73,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+3 ,1:6) = (/ 128,  74,   1,   0,   0, 4 /)
      FA%NCODPA(INBPAR+4 ,1:6) = (/ 128,  75,   1,   0,   0, 4 /)
C
C     1.11- TOP
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 1
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'TOP'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'RAYT DIR SOM'
      FA%NCODPA(INBPAR+1 ,1:6) = (/ 149, 174,   8,   0,   0, 4 /)
C**
C     2.  -  TRAITEMENT DES CHAMPS SUR NIVEAUX
C-----------------------------------------------------------------------
C
C     2.1 -  Champs sur un niveau isobare
C
C  ATTENTION: comme les champs HAUTEUR sont tres proches
C             des champs ISOBARE, un copier-coller sera fait
C             pour les INB champs suivants:
C             Revoir la methode si un champ PRESSION ne suivait
C             pas cette regle !!
      INB = 32
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + INB
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'P'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'POT_VORTIC'
      FA%CISUFF(INBPAR+2 ) = 'GEOPOTENTI'
      FA%CISUFF(INBPAR+3 ) = 'TEMPERATUR'
      FA%CISUFF(INBPAR+4 ) = 'TEMPE_POTE'
      FA%CISUFF(INBPAR+5 ) = 'THETA_PRIM'
      FA%CISUFF(INBPAR+6 ) = 'VENT_ZONAL'
      FA%CISUFF(INBPAR+7 ) = 'VENT_MERID'
      FA%CISUFF(INBPAR+8 ) = 'VITESSE_VE'
      FA%CISUFF(INBPAR+9 ) = 'ABS_VORTIC'
      FA%CISUFF(INBPAR+10) = 'DIVERGENCE'
      FA%CISUFF(INBPAR+11) = 'HUMI.SPECI'
      FA%CISUFF(INBPAR+12) = 'HUMI_RELAT'
C Debut des champs pour Mocage :
      FA%CISUFF(INBPAR+13) = 'O3'
      FA%CISUFF(INBPAR+14) = 'NO2'
      FA%CISUFF(INBPAR+15) = 'NO'
      FA%CISUFF(INBPAR+16) = 'PAN'
      FA%CISUFF(INBPAR+17) = 'HNO3'
      FA%CISUFF(INBPAR+18) = 'SO2'
      FA%CISUFF(INBPAR+19) = 'CO'
      FA%CISUFF(INBPAR+20) = 'CH4'
      FA%CISUFF(INBPAR+21) = 'C2H6'
      FA%CISUFF(INBPAR+22) = 'NC4H10'
      FA%CISUFF(INBPAR+23) = 'C2H4'
      FA%CISUFF(INBPAR+24) = 'C3H6'
      FA%CISUFF(INBPAR+25) = 'OXYL'
      FA%CISUFF(INBPAR+26) = 'C5H8'
      FA%CISUFF(INBPAR+27) = 'APINEN'
      FA%CISUFF(INBPAR+28) = 'HCHO'
      FA%CISUFF(INBPAR+29) = 'CH3CHO'
      FA%CISUFF(INBPAR+30) = 'GLYOX'
      FA%CISUFF(INBPAR+31) = 'MGLYOX'
      FA%CISUFF(INBPAR+32) = 'CH3COE'
C Fin des champs pour Mocage
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,   4, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+2 ,1:6) = (/   1,   6, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+3 ,1:6) = (/   1,  11, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+4 ,1:6) = (/   1,  13, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+5 ,1:6) = (/   1,  15, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+6 ,1:6) = (/   1,  33, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+7 ,1:6) = (/   1,  34, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+8 ,1:6) = (/   1,  39, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+9 ,1:6) = (/   1,  41, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+10,1:6) = (/   1,  44, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+11,1:6) = (/   1,  51, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+12,1:6) = (/   1,  52, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+13,1:6) = (/ 159,   1, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+14,1:6) = (/ 159,   2, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+15,1:6) = (/ 159,   3, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+16,1:6) = (/ 159,   4, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+17,1:6) = (/ 159,   5, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+18,1:6) = (/ 159,   6, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+19,1:6) = (/ 159,   7, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+20,1:6) = (/ 159,   8, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+21,1:6) = (/ 159,   9, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+22,1:6) = (/ 159,  10, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+23,1:6) = (/ 159,  11, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+24,1:6) = (/ 159,  12, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+25,1:6) = (/ 159,  13, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+26,1:6) = (/ 159,  14, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+27,1:6) = (/ 159,  15, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+28,1:6) = (/ 159,  16, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+29,1:6) = (/ 159,  17, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+30,1:6) = (/ 159,  18, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+31,1:6) = (/ 159,  19, 100,   0,   0, 0 /)
      FA%NCODPA(INBPAR+32,1:6) = (/ 159,  20, 100,   0,   0, 0 /)
C
C     2.2 -  Champs sur un niveau hauteur
C
C   ATTENTION: on reprend les INB definitions du niveau ISOBARE!!
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + INB
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'H'
        FA%CISUFF(J) = FA%CISUFF(J-INB)
        FA%NCODPA(J,1:2) = FA%NCODPA(J-INB,1:2)
        FA%NCODPA(J,3) = 105
        FA%NCODPA(J,4:6) = FA%NCODPA(J-INB,4:6)
      ENDDO
C
C     2.3 -  Champs sur un niveau iso-tourbillon-potentiel
C
C  ATTENTION: comme les champs ISENTROPE sont tres proches
C             des champs ISO-PV, un copier-coller sera fait
C             pour les INB champs suivants:
C             Revoir la methode si un champ ISO-PV ne suivait
C             pas cette regle !!
      INB = 10
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + INB
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'V'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'GEOPOTENTIEL'
      FA%CISUFF(INBPAR+2 ) = 'TEMPERATURE'
      FA%CISUFF(INBPAR+3 ) = 'TEMPE_POTENT'
      FA%CISUFF(INBPAR+4 ) = 'VENT_ZONAL'
      FA%CISUFF(INBPAR+5 ) = 'VENT_MERIDIE'
      FA%CISUFF(INBPAR+6 ) = 'VITESSE_VE'
      FA%CISUFF(INBPAR+7 ) = 'ABS_VORTICIT'
      FA%CISUFF(INBPAR+8 ) = 'DIVERGENCE'
      FA%CISUFF(INBPAR+9 ) = 'HUMI.SPECI'
      FA%CISUFF(INBPAR+10) = 'HUMI_RELAT'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,   6, 117,   0,   0, 0 /)
      FA%NCODPA(INBPAR+2 ,1:6) = (/   1,  11, 117,   0,   0, 0 /)
      FA%NCODPA(INBPAR+3 ,1:6) = (/   1,  13, 117,   0,   0, 0 /)
      FA%NCODPA(INBPAR+4 ,1:6) = (/   1,  33, 117,   0,   0, 0 /)
      FA%NCODPA(INBPAR+5 ,1:6) = (/   1,  34, 117,   0,   0, 0 /)
      FA%NCODPA(INBPAR+6 ,1:6) = (/   1,  39, 117,   0,   0, 0 /)
      FA%NCODPA(INBPAR+7 ,1:6) = (/   1,  41, 117,   0,   0, 0 /)
      FA%NCODPA(INBPAR+8 ,1:6) = (/   1,  44, 117,   0,   0, 0 /)
      FA%NCODPA(INBPAR+9 ,1:6) = (/   1,  51, 117,   0,   0, 0 /)
      FA%NCODPA(INBPAR+10,1:6) = (/   1,  52, 117,   0,   0, 0 /)

C
C     2.4 -  Champs sur un niveau isentrope
C
C   ATTENTION: on reprend les INB definitions du niveau ISO-PV!!
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + INB
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'T'
        FA%CISUFF(J) = FA%CISUFF(J-INB)
        FA%NCODPA(J,1:2) = FA%NCODPA(J-INB,1:2)
        FA%NCODPA(J,3) = 113
        FA%NCODPA(J,4:6) = FA%NCODPA(J-INB,4:6)
      ENDDO
C Champs suppl differents de ceux d'ISO-PV
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 1
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'T'
      ENDDO
      FA%CISUFF(INBPAR+1) = 'POT_VORTICIT'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,   4, 113,   0,   0, 0 /)
C
C     2.5 -  Champs sur un niveau modele
C
C
      INBPAR = FA%NBPARC
      FA%NBPARC = INBPAR + 12
      IF (FA%NBPARC.GT.FA%JPXPAR) GOTO 1001
      DO J=INBPAR+1,FA%NBPARC
        FA%CIPREF(J) = 'S'
      ENDDO
      FA%CISUFF(INBPAR+1 ) = 'TEMPERATURE'
      FA%CISUFF(INBPAR+2 ) = 'WIND.U.PHYS'
      FA%CISUFF(INBPAR+3 ) = 'VENT_ZONAL'
      FA%CISUFF(INBPAR+4 ) = 'WIND.V.PHYS'
      FA%CISUFF(INBPAR+5 ) = 'VENT_MERIDIE'
      FA%CISUFF(INBPAR+6 ) = 'FONC.COURANT'
      FA%CISUFF(INBPAR+7 ) = 'POT.VITESSE'
      FA%CISUFF(INBPAR+8 ) = 'VORTICITY'
      FA%CISUFF(INBPAR+9 ) = 'DIVERGENCE'
      FA%CISUFF(INBPAR+10) = 'HUMI.SPECIFI'
      FA%CISUFF(INBPAR+11) = 'VERTIC.DIVER'
      FA%CISUFF(INBPAR+12) = 'PRESS.DEPART'
      FA%NCODPA(INBPAR+1 ,1:6) = (/   1,  11, 109,   0,   0, 0 /)
      FA%NCODPA(INBPAR+2 ,1:6) = (/   1,  33, 109,   0,   0, 0 /)
      FA%NCODPA(INBPAR+3 ,1:6) = (/   1,  33, 109,   0,   0, 0 /)
      FA%NCODPA(INBPAR+4 ,1:6) = (/   1,  34, 109,   0,   0, 0 /)
      FA%NCODPA(INBPAR+5 ,1:6) = (/   1,  34, 109,   0,   0, 0 /)
      FA%NCODPA(INBPAR+6 ,1:6) = (/   1,  35, 109,   0,   0, 0 /)
      FA%NCODPA(INBPAR+7 ,1:6) = (/   1,  36, 109,   0,   0, 0 /)
      FA%NCODPA(INBPAR+8 ,1:6) = (/   1,  43, 109,   0,   0, 0 /)
      FA%NCODPA(INBPAR+9 ,1:6) = (/   1,  44, 109,   0,   0, 0 /)
      FA%NCODPA(INBPAR+10,1:6) = (/   1,  51, 109,   0,   0, 0 /)
      FA%NCODPA(INBPAR+11,1:6) = (/ 149,  91, 109,   0,   0, 0 /)
      FA%NCODPA(INBPAR+12,1:6) = (/ 149,  92, 109,   0,   0, 0 /)
C
C**
C    10.  -  PHASE TERMINALE
C-----------------------------------------------------------------------
C
      IF (LHOOK) CALL DR_HOOK('FAICOR_MT',1,ZHOOK_HANDLE)
      RETURN
C
 1001 CONTINUE
      WRITE (UNIT=FA%NULOUT,FMT=*)
     S            'FAICOR: ERROR, ',FA%NBPARC,' fields are initialized'
      WRITE (UNIT=FA%NULOUT,FMT=*)
     S            '        but only ',FA%JPXPAR,' are possible.'
      WRITE (UNIT=FA%NULOUT,FMT=*)'        Change FA%JPXPAR in famodu.F'
      STOP
C
      END

