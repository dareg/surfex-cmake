C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FANFAR_MT (FA,  KREP, KRANG, CDPREF, KNIVAU, CDSUFF, 
     S                      CDNOMA, KB1PAR, KLPRFU, KLSUFU, KLNOMU )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C     fabrication d'un nom article devant contenir un champ horizontal.
C              (Nom a Fabriquer pour un ARticle)
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
C                KNIVAU (Entree) ==> Niveau vertical eventuel;
C                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
C                CDNOMA (Sortie) ==> Nom d'article;
C    ( Tableau ) KB1PAR (Sortie) ==> 3 elements du "Bloc 1" des
C                                    interfaces GRIB concernes par la
C                                    coordonnee et le niveau verticaux;
C                KLPRFU (Sortie) ==> Longueur utile du prefixe;
C                KLSUFU (Sortie) ==> Longueur utile du suffixe;
C                KLNOMU (Sortie) ==> Longueur utile du nom d'article.
C*
C       En mode multi-taches, il doit y avoir verrouillage du fichier
C     concerne avant l'appel au sous-programme.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANG, KNIVAU, KLPRFU, KLSUFU, KLNOMU
      INTEGER KB1PAR (3)
C
      CHARACTER CDPREF*(*), CDSUFF*(*), CDNOMA*(*), CLAUXI*(FA%JPXNOM)
C
      INTEGER ILPREF, ILSUFF, ILCDNO, ILPRFU, ILSUFU, J, INCHIF, INIMES
      INTEGER INUMER, ILACTI, ILNOMA, ILAUXI, ITYNIV, INIVAU
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FANFAR_MT',0,ZHOOK_HANDLE)
      ILPREF=LEN (CDPREF)
      ILSUFF=LEN (CDSUFF)
      ILCDNO=LEN (CDNOMA)
      ILPRFU=MAX (0,ILPREF)
      ILSUFU=MAX (0,ILSUFF)
      KLNOMU=MAX (0,ILCDNO)
C
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C
      IF (ILCDNO.LE.0.OR.ILCDNO.GT.FA%NCPCAD) THEN
        KREP=-66
        GOTO 1001
      ELSEIF (MIN (ILPREF,ILSUFF).LE.0) THEN
        KREP=-65
        GOTO 1001
      ELSEIF (CDPREF.EQ.' '.OR.CDSUFF.EQ.' ') THEN
        KREP=-86
        GOTO 1001
      ENDIF
C
C       Decompte du nombre de caracteres utiles dans prefixe et suffixe.
C
      DO 103 J=ILPREF,1,-1
C
      IF (CDPREF(J:J).NE.' ') THEN
        ILPRFU=J
        GOTO 104
      ENDIF
C
  103 CONTINUE
C
      KREP=-66
      GOTO 1001
C
  104 CONTINUE
C
      IF (ILPRFU.GT.FA%JPXPRF) THEN
        KREP=-87
        GOTO 1001
      ENDIF
C
      DO 105 J=ILSUFF,1,-1
C
      IF (CDSUFF(J:J).NE.' ') THEN
        ILSUFU=J
        GOTO 106
      ENDIF
C
  105 CONTINUE
C
      KREP=-66
      GOTO 1001
C
  106 CONTINUE
C**
C     2.  -  DIFFERENTS CAS, SELON LE PREFIXE.
C-----------------------------------------------------------------------
C
      DO 201 J=1,FA%JPTNIV
C
      IF (CDPREF.EQ.FA%CTNPRF(J)) THEN
        ITYNIV=J
        GOTO 202
      ENDIF
C
  201 CONTINUE
C
      ITYNIV=0
C
  202 CONTINUE
C
      INCHIF=FA%NIVDSC(0,ITYNIV)
C
      IF (INCHIF.EQ.0) THEN
C
        INIVAU=0
C
      ELSEIF (KNIVAU.LT.FA%NIVDSC(1,ITYNIV).OR.
     S        KNIVAU.GT.FA%NIVDSC(2,ITYNIV)) THEN
C
        KREP=-64
        GOTO 1001
C
      ELSEIF (CDPREF.EQ.'S'.AND.KNIVAU.GT.
     S        FA%NNIVER(FA%NUCADR(KRANG))) THEN
C
        KREP=-64
        GOTO 1001
C
      ELSEIF (CDPREF.EQ.'L'.AND.KNIVAU.GT.
     S        FA%NNIVER(FA%NUCADR(KRANG))) THEN
C
        KREP=-64
        GOTO 1001
C
      ELSE
C
        INIVAU=KNIVAU
C
      ENDIF
C
      KB1PAR(1)=FA%NIVDSC(3,ITYNIV)
      KB1PAR(2)=FA%NIVDSC(4,ITYNIV)
      KB1PAR(3)=INIVAU
C
      ILSUFU=MIN (ILCDNO-ILPRFU-INCHIF,ILSUFU)
      KLNOMU=ILPRFU+INCHIF+ILSUFU
C
      IF (INCHIF.NE.0) THEN
        WRITE (UNIT=CLNOMA,FMT='(I8.8)') KNIVAU
        CDNOMA=CDPREF(1:ILPRFU)//CLNOMA(9-INCHIF:8)//CDSUFF(1:ILSUFU)
      ELSE
        CDNOMA=CDPREF(1:ILPRFU)//CDSUFF(1:ILSUFU)
      ENDIF
C
      IF (CDNOMA.EQ.FA%CPCACH.OR.CDNOMA.EQ.FA%CPCADI.OR.CDNOMA.EQ.
     S    FA%CPCAFS.OR.
     S    CDNOMA.EQ.FA%CPCARP.OR.CDNOMA.EQ.FA%CPDATE.OR.
     S    CDNOMA.EQ.FA%CIDENT(KRANG)) THEN
        KREP=-111
        GOTO 1001
      ENDIF
C
      KREP=0
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE EVENTUELLE,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      KLPRFU=ILPRFU
      KLSUFU=ILSUFU
      LLFATA=LLMOER (KREP,KRANG)
C
      IF (FA%LFAMOP.OR.LLFATA) THEN
        INIMES=2
        CLNSPR='FANFAR'
        INUMER=FA%JPNIIL
C
        IF (ILPRFU.GE.1) THEN
          ILACTI=MIN (ILPRFU,FA%NCPCAD)
          CLACTI(1:ILACTI)=CDPREF(1:ILACTI)
        ELSE
          ILACTI=8
          CLACTI(1:ILACTI)=FA%CHAINC(:ILACTI)
        ENDIF
C
        IF (ILSUFU.GE.1) THEN
          ILNOMA=MIN (ILSUFU,FA%NCPCAD)
          CLNOMA(1:ILNOMA)=CDSUFF(1:ILNOMA)
        ELSE
          ILNOMA=8
          CLNOMA(1:ILNOMA)=FA%CHAINC(:ILNOMA)
        ENDIF
C
        IF (KLNOMU.GE.1) THEN
          ILAUXI=MIN (KLNOMU,LEN(CLAUXI))
          CLAUXI(1:ILAUXI)=CDNOMA(1:ILAUXI)
        ELSE
          ILAUXI=8
          CLAUXI(1:ILAUXI)=FA%CHAINC(:ILAUXI)
        ENDIF
C
        WRITE (UNIT=CLMESS,
     S         FMT='(''ARGUMENTS='',2(I4,'', ''),'''''''',A,
     S         '''''','',I6,'', '''''',A,'''''',  '''''',A,'''''''',
     S         2('','',I4),'','',I6,3('','',I3))')
     S     KREP,KRANG,CLACTI(1:ILACTI),KNIVAU,CLNOMA(1:ILNOMA),
     S     CLAUXI(1:ILAUXI),KB1PAR,KLPRFU,KLSUFU,KLNOMU
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR, CLACTI(1:ILACTI),.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FANFAR_MT',1,ZHOOK_HANDLE)
      END






