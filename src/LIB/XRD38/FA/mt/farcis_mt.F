C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FARCIS_MT (FA,  KREP, KRANG, PCHAMP, KSTRON, KPUILA )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C     elimination de la "puissance de laplacien" d'un champ en coeffi-
C     cients spectraux issu d'un codage GRIB, de maniere a restituer
C     le champ "d'origine" (a la precision du codage pres) .
C     ( Reconstitution des CoeffIcients Spectraux )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C    ( Tableau ) PCHAMP (Entree ET Sortie) ==> Champ en coef. spectraux;
C                KSTRON (Entree) ==> Sous-troncature non compactee;
C                KPUILA (Entree) ==> Puissance de laplacien utilisee.
C
C      ( Les 2 derniers parametres sont ceux qui ont ete effectivement
C        utilises lors de l'ecriture du champ )
C*
C       En mode multi-taches, il doit y avoir verrouillage du fichier
C     concerne avant l'appel au sous-programme.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANG, KSTRON, KPUILA
C
      REAL (KIND=JPDBLR) PCHAMP (FA%JPXCSP)
C
      INTEGER IRANGC, ITRONC, INUMER, IDIMNC, ILCHAM, IMTRONC, IPUISX, J
      INTEGER INDICE, JN, JM, INDLAP, IMLIM, IOFF, IM, JIND, IPUIS2
      INTEGER IRAPOR, IPUISR, INIMES, IDEB, IFIN
C
      LOGICAL LLMLAM
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FARCIS_MT',0,ZHOOK_HANDLE)
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C
      IF (FA%LIXLAP) THEN
        CALL FAIXLA_MT (FA)
        FA%LIXLAP=.FALSE.
      ENDIF
C
      IRANGC=FA%NUCADR(KRANG)
      ITRONC=FA%MTRONC(IRANGC)
      LLMLAM=FA%LIMLAM(IRANGC)
C
      IF (LLMLAM) IMTRONC=FA%NOZPAR(2,IRANGC)
      IF (ITRONC.LE.KSTRON) THEN
        KREP=-88
        GOTO 1001
      ELSEIF (LLMLAM.AND.IMTRONC.LE.KSTRON) THEN
        KREP=-88
        GOTO 1001
      ELSEIF (LLMLAM.AND.(IMTRONC.GT.3*ITRONC.
     S    OR.ITRONC.GT.3*IMTRONC)) THEN
C Il s'agit d'un garde-fou, modifiable (ne pas oublier FAPULA et FACSIM)
        KREP=-114
        GOTO 1001
      ELSE        
        KREP=0
      ENDIF
C
      IDIMNC=(1+KSTRON)**2
      IF (LLMLAM) THEN
        ILCHAM=FA%NSFLAM(IRANGC)
      ELSE  
        ILCHAM=(1+ITRONC)**2
      ENDIF  
C**
C     2.  -  RECONSTITUTION DU CHAMP "D'ORIGINE", DEBARRASSE DE LA
C            PUISSANCE DE LAPLACIEN QUI N'AFFECTE QUE LA PARTIE HORS
C            SOUS-TRONCATURE NON COMPACTEE.
C-----------------------------------------------------------------------
C
C        On essaie d'eviter l'exponentiation, en preferant multiplier
C     que diviser.
C
      IF (KPUILA.NE.0) THEN
C
        IPUISX=IABS (KPUILA)
C
        IF (KPUILA.GT.0) THEN
          INDICE=1
        ELSE
          INDICE=0
        ENDIF
C
        IF (IPUISX.LE.FA%JPUILA) THEN
C
          IF (LLMLAM) THEN
!$OMP PARALLEL DO PRIVATE(JN,IMLIM,IDEB,IFIN,JIND,IOFF,IM,INDLAP)
!$OMP&  IF(FA%LOPENMP)
            DO 2010 JN=1,ITRONC
            IMLIM=KSTRON-JN
            IDEB=MAX(FA%NOZPAR(2*JN+3,IRANGC)+4*(1+IMLIM),
     S             FA%NOZPAR(2*JN+3,IRANGC)+4)
            IFIN=FA%NOZPAR(2*JN+4,IRANGC)
            DO 2010 JIND=IDEB,IFIN
            IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
            IM=IOFF/4
            INDLAP=((JN-1)*FA%JPXTRO)+IM
            PCHAMP(JIND)=PCHAMP(JIND)*FA%XLAP2DA(INDLAP,IPUISX,INDICE)
 2010       CONTINUE         
!$OMP END PARALLEL DO
          ELSE
            DO 201 J=IDIMNC+1,ILCHAM
            PCHAMP(J)=PCHAMP(J)*FA%XLAP2D(J,IPUISX,INDICE)
  201       CONTINUE
          ENDIF
        ELSEIF (IPUISX.LE.2*FA%JPUILA) THEN
          IPUIS2=IPUISX/2
C
          IF (IPUISX.EQ.2*IPUIS2) THEN
C
            IF (LLMLAM) THEN
!$OMP PARALLEL DO PRIVATE(JN,IMLIM,IDEB,IFIN,JIND,IOFF,IM,INDLAP)
!$OMP&  IF(FA%LOPENMP)
              DO 2020 JN=1,ITRONC
              IMLIM=KSTRON-JN
              IDEB=MAX(FA%NOZPAR(2*JN+3,IRANGC)+4*(1+IMLIM),
     S               FA%NOZPAR(2*JN+3,IRANGC)+4)
              IFIN=FA%NOZPAR(2*JN+4,IRANGC)
              DO 2020 JIND=IDEB,IFIN
              IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
              IM=IOFF/4
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              PCHAMP(JIND)=PCHAMP(JIND)
     S          *( FA%XLAP2DA(INDLAP,IPUIS2,INDICE)**2 )
 2020         CONTINUE               
!$OMP END PARALLEL DO
            ELSE
              DO 202 J=IDIMNC+1,ILCHAM
              PCHAMP(J)=PCHAMP(J)*( FA%XLAP2D(J,IPUIS2,INDICE)**2 )
  202         CONTINUE
            ENDIF
C
          ELSE
C
            IF (LLMLAM) THEN
!$OMP PARALLEL DO PRIVATE(JN,IMLIM,IDEB,IFIN,JIND,IOFF,IM,INDLAP) 
!$OMP&  IF(FA%LOPENMP)
              DO 2030 JN=1,ITRONC
              IMLIM=KSTRON-JN
              IDEB=MAX(FA%NOZPAR(2*JN+3,IRANGC)+4*(1+IMLIM),
     S               FA%NOZPAR(2*JN+3,IRANGC)+4)
              IFIN=FA%NOZPAR(2*JN+4,IRANGC)
              DO 2030 JIND=IDEB,IFIN
              IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
              IM=IOFF/4
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              PCHAMP(JIND)=PCHAMP(JIND)
     S          *( FA%XLAP2DA(INDLAP,FA%JPUILA,INDICE)
     S            *FA%XLAP2DA(INDLAP,IPUISX-FA%JPUILA,INDICE) )
 2030         CONTINUE              
!$OMP END PARALLEL DO
            ELSE
              DO 203 J=IDIMNC+1,ILCHAM
              PCHAMP(J)=PCHAMP(J)*( FA%XLAP2D(J,FA%JPUILA,INDICE)
     S                          *FA%XLAP2D(J,IPUISX-FA%JPUILA,INDICE) )
  203         CONTINUE
            ENDIF
C
          ENDIF
C
        ELSE
          IRAPOR=1+(IPUISX-1)/FA%JPUILA
          IPUISR=IPUISX/IRAPOR
C
          IF (IPUISX.EQ.IRAPOR*IPUISR) THEN
C
            IF (LLMLAM) THEN
!$OMP PARALLEL DO PRIVATE(JN,IMLIM,IDEB,IFIN,JIND,IOFF,IM,INDLAP)
!$OMP&  IF(FA%LOPENMP)
              DO 2040 JN=1,ITRONC
              IMLIM=KSTRON-JN
              IDEB=MAX(FA%NOZPAR(2*JN+3,IRANGC)+4*(1+IMLIM),
     S               FA%NOZPAR(2*JN+3,IRANGC)+4)
              IFIN=FA%NOZPAR(2*JN+4,IRANGC)
              DO 2040 JIND=IDEB,IFIN
              IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
              IM=IOFF/4
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              PCHAMP(JIND)=PCHAMP(JIND)
     S           *( FA%XLAP2DA(INDLAP,IPUISR,INDICE)**IRAPOR )
 2040         CONTINUE                
!$OMP END PARALLEL DO
            ELSE
              DO 204 J=IDIMNC+1,ILCHAM
              PCHAMP(J)=PCHAMP(J)*( FA%XLAP2D(J,IPUISR,INDICE)**IRAPOR )
  204         CONTINUE
            ENDIF
C
          ELSE
C
            IF (LLMLAM) THEN
!$OMP PARALLEL DO PRIVATE(JN,IMLIM,IDEB,IFIN,JIND,IOFF,IM,INDLAP)
!$OMP&  IF(FA%LOPENMP)
              DO 2050 JN=1,ITRONC
              IMLIM=KSTRON-JN
              IDEB=MAX(FA%NOZPAR(2*JN+3,IRANGC)+4*(1+IMLIM),
     S               FA%NOZPAR(2*JN+3,IRANGC)+4)
              IFIN=FA%NOZPAR(2*JN+4,IRANGC)
              DO 2050 JIND=IDEB,IFIN
              IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
              IM=IOFF/4
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              PCHAMP(JIND)=PCHAMP(JIND)
     S          *( FA%XLAP2DA(INDLAP,FA%JPUILA,INDICE)**(IRAPOR-1)
     S         *FA%XLAP2DA(INDLAP,IPUISX-FA%JPUILA*(IRAPOR-1),INDICE) )
 2050         CONTINUE             
!$OMP END PARALLEL DO
            ELSE
              DO 205 J=IDIMNC+1,ILCHAM
              PCHAMP(J)=PCHAMP(J)* 
     S                 (FA%XLAP2D(J,FA%JPUILA,INDICE)**(IRAPOR-1)
     S                *FA%XLAP2D(J,IPUISX-FA%JPUILA*(IRAPOR-1),INDICE) )
  205         CONTINUE
            ENDIF
C
          ENDIF
C
        ENDIF
C
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE EVENTUELLE,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      LLFATA=LLMOER (KREP,KRANG)
C
      IF (FA%LFAMOP.OR.LLFATA) THEN
        INIMES=2
        CLNSPR='FARCIS'
        INUMER=FA%JPNIIL
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4,
     S    '', PCHAMP(1)='',G12.5,'', KSTRON='',I4,'', KPUILA='',I3)')
     S     KREP,KRANG,PCHAMP(1),KSTRON,KPUILA
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S               CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FARCIS_MT',1,ZHOOK_HANDLE)
      END

