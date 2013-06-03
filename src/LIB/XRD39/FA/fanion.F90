! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FANION_MT_BB                                          &
&                     (FA,  KREP, KNUMER, CDPREF, KNIVAU, CDSUFF,  &
&                      LDEXIS,                                     &
&                      LDCOSP, KNGRIB, KNBITS, KSTRON, KPUILA )
USE FA_MOD, ONLY : FA_COM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!      Sous-programme renseignant sur l'EXISTENCE et les CARACTERISTI-
!     QUES eventuelles d'un Article de type CHAMP dans un Fichier ARPEGE
!       ( LDEXIS est le "fanion" leve si l'article existe )
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KNUMER (Entree) ==> Numero de l'unite logique;
!                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
!                KNIVAU (Entree) ==> Niveau vertical eventuel;
!                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
!                LDEXIS (Sortie) ==> Vrai si l'article de type CHAMP
!                                    existe bien dans le Fichier;
!                LDCOSP (Sortie) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux;
!                KNGRIB (Sortie) ==> Niveau de codage GRIB;
!                KNBITS (Sortie) ==> Nombre de bits de codage eventuel;
!                KSTRON (Sortie) ==> Sous-troncature non codee   " -le;
!                KPUILA (Sortie) ==> Puissance de laplacien eventuelle.
!
!        KNBITS n'a de sens que si l'article existe et a ete code;
!     de meme pour KSTRON et KPUILA, qui ne sont applicables qu'a un
!     champ represente en coefficients spectraux.
!        Les arguments de sortie n'ayant pas de sens sont mis a
!     0 pour les entiers, .FALSE. pour les logiques.
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KNUMER, KNIVAU, KNGRIB
INTEGER (KIND=JPLIKB) KNBITS, KSTRON, KPUILA
!
INTEGER (KIND=JPLIKB) IREP, ILPRFU, ILSUFU, ILNOMU
INTEGER (KIND=JPLIKB) ILONGA, IRANG, INIMES
INTEGER (KIND=JPLIKB) ILPREF, ILSUFF, IPOSEX
INTEGER (KIND=JPLIKB) IRANGC, ILCHAM
!
INTEGER (KIND=JPDBLE) IVALCO (5)
INTEGER (KIND=JPLIKB)  IB1PAR (3)
!
LOGICAL LLVERF, LLRLFI, LDCOSP, LDEXIS, LLTEMP, LLNOMU, LLMLAM
!
CHARACTER CDPREF*(*), CDSUFF*(*)
CHARACTER CLPREF*(FA%JPXNOM), CLSUFF*(FA%JPXSUF)
!
CHARACTER(LEN=FA%JPXNOM) CLNOMA
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FANION_MT',0,ZHOOK_HANDLE)
IREP=0
LLVERF=.FALSE.
LLTEMP=.FALSE.
LLRLFI=.FALSE.
LLNOMU=.FALSE.
ILPRFU=INT (LEN (CDPREF), JPLIKB)
ILSUFU=INT (LEN (CDSUFF), JPLIKB)
LDEXIS=.FALSE.
LDCOSP=.FALSE.
KNGRIB=0
KNBITS=0
KSTRON=0
KPUILA=0
CALL FANUMU_MT_BB                &
&               (FA, KNUMER,IRANG)
!
IF (IRANG.EQ.0) THEN
  IREP=-51
  GOTO 1001
ENDIF
!
!         Verrouillage eventuel du fichier.
!
IF (FA%LFAMUL) CALL LFIVER_MT_BB                             &
&                              (FA%LFI, FA%FICHIER(IRANG)%VRFICH,'ON')
LLVERF=FA%LFAMUL
!
IF (FA%FICHIER(IRANG)%LCREAF) GOTO 1001
!**
!     2.  -  FABRICATION DU NOM D'ARTICLE VIA LE SOUS-PROGRAMME "FANFAR"
!            ( controles de CDPREF, KNIVAU, CDSUFF inclus )
!-----------------------------------------------------------------------
!
CALL FANFAR_MT_BB                                          &
&               (FA, IREP,IRANG,CDPREF,KNIVAU,CDSUFF,CLNOMA, &
&                IB1PAR,ILPRFU,ILSUFU,ILNOMU)
IF (IREP.NE.0) GOTO 1001
LLNOMU=.TRUE.
!**
!     3.  -  RECHERCHE DE L'ARTICLE SUR LE FICHIER, LECTURE PARTIELLE.
!-----------------------------------------------------------------------
!
CALL LFINFO_MT_BB                                    &
&               (FA%LFI, IREP,KNUMER,CLNOMA(1:ILNOMU), &
&             ILONGA,IPOSEX)
LLRLFI=IREP.NE.0
IF (LLRLFI.OR.ILONGA.EQ.0) GOTO 1001
LDEXIS=.TRUE.
!
IF (ILONGA.GT.FA%JPXCHA+2) THEN
  IREP=-90
  GOTO 1001
ENDIF
!
IF (FA%FICHIER(IRANG)%LERRFA) THEN
!
!        Le fichier est gere en mode "toute erreur est fatale".
!     Ce mode etant normalement couple au mode correspondant du logiciel
!     LFI, on va temporairement annuler l'option LFI afin de pouvoir
!     faire une lecture partielle de l'entete de l'article Champ.
!
  CALL LFIERF_MT_BB                           &
&                 (FA%LFI, IREP,KNUMER,.FALSE.)
  LLRLFI=IREP.NE.0
  IF (LLRLFI) GOTO 1001
  LLTEMP=.TRUE.
ENDIF
!
CALL LFILEC_MT_BB                                    &
&               (FA%LFI, IREP,KNUMER,CLNOMA(1:ILNOMU), &
&               IVALCO,5_JPLIKB )
!
IF (IREP.EQ.0) THEN
  IREP=-93
  GOTO 1001
ELSEIF (IREP.NE.-21) THEN
  LLRLFI=.TRUE.
  GOTO 1001
ELSEIF (IVALCO(1).LT.-1 .OR. IVALCO(1).GT.3 .OR.                   &
&        IVALCO(2).LT.0  .OR. IVALCO(2).GT.1 .OR.                   &
&  (IVALCO(1).GT.0 .AND. IVALCO(2).EQ.1 .AND. IVALCO(4).LT.0)) THEN
  IREP=-91
  GOTO 1001
ELSE
  IREP=0
  KNGRIB=IVALCO(1)
  LDCOSP=IVALCO(2).EQ.1
ENDIF
!
IRANGC=FA%FICHIER(IRANG)%NUCADR
LLMLAM=FA%CADRE(IRANGC)%LIMLAM
!
IF (LDCOSP) THEN
  IF (LLMLAM) THEN
     ILCHAM=FA%CADRE(IRANGC)%NSFLAM
  ELSE
     IF (KNGRIB.EQ.3 .OR. KNGRIB.EQ.-1) THEN
       ILCHAM=(1+FA%CADRE(IRANGC)%MTRONC)*(2+FA%CADRE(IRANGC)%MTRONC)
     ELSE
       ILCHAM=(1+FA%CADRE(IRANGC)%MTRONC)**2
     ENDIF
  ENDIF
ELSE
  ILCHAM=FA%CADRE(IRANGC)%NVAPDG
ENDIF
!
IF (KNGRIB.EQ.-1 .OR. KNGRIB.EQ.0) THEN
!
!          Cas ou il n'y a aucun codage... controle longueur d'article
!
  IF (ILONGA.LT.(ILCHAM+2)) THEN
    IREP=-93
    GOTO 1001
  ELSEIF (ILONGA.GT.(ILCHAM+2)) THEN
    IREP=-94
    IF (LLMOER(IREP,IRANG)) GOTO 1001
  ENDIF
!
ELSE
!
!        Cas avec codage GRIB (standard ou non).
!
  KNBITS=IVALCO(3)
!
  IF (LDCOSP) THEN
    KSTRON=IVALCO(4)
    KPUILA=IVALCO(5)
!
    IF (KNGRIB.EQ.2.AND.ILONGA.LT.(5+(1+KSTRON)**2)) THEN
      IREP=-93
      GOTO 1001
    ENDIF
!
  ENDIF
!
ENDIF
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
!
IF (LLTEMP) THEN
!
!         On remet le fichier en mode "toute erreur fatale" au niveau
!     du logiciel LFI.
!
  CALL LFIERF_MT_BB                          &
&                 (FA%LFI, IREP,KNUMER,.TRUE.)
  LLRLFI=IREP.NE.0
ENDIF
!
KREP=IREP
LLFATA=LLMOER (IREP,IRANG)
!
!        Deverrouillage eventuel du fichier.
!
IF (LLVERF) CALL LFIVER_MT_BB                              &
&                           (FA%LFI, FA%FICHIER(IRANG)%VRFICH,'OFF')
!
IF (LLFATA) THEN
  INIMES=2
ELSE
  INIMES=IXNVMS(IRANG)
ENDIF
!
IF (.NOT.LLFATA.AND.INIMES.NE.2)  THEN 
  IF (LHOOK) CALL DR_HOOK('FANION_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='FANION'
!
IF (ILPRFU.GE.1) THEN
  ILPREF=MIN (ILPRFU,INT (LEN (CLPREF), JPLIKB))
  CLPREF(1:ILPREF)=CDPREF(1:ILPREF)
ELSE
  ILPREF=8
  CLPREF(1:ILPREF)=FA%CHAINC(:ILPREF)
ENDIF
!
IF (ILSUFU.GE.1) THEN
  ILSUFF=MIN (ILSUFU,INT (LEN (CLSUFF), JPLIKB))
  CLSUFF(1:ILSUFF)=CDSUFF(1:ILSUFF)
ELSE
  ILSUFF=8
  CLSUFF(1:ILSUFF)=FA%CHAINC(:ILSUFF)
ENDIF
!
IF (.NOT.LLNOMU) THEN
  ILNOMU=MIN (ILPREF,FA%NCPCAD)
  CLNOMA(1:ILNOMU)=CLPREF(1:ILPREF)
ENDIF
!
WRITE (UNIT=CLMESS,                                             &
&       FMT='(''ARGUMENTS:'',I4,'','',I3,'','''''',A,            &
&       '''''','',I6,'','''''',A,'''''', LDEXIS= '',L1,          &
&       '', LDCOSP= '',L1,'', KNGRIB='',I2,'', KNBITS='',I3,     &
&       '',KSTRON='',I3,'',KPUILA='',I6)')                       &
&   KREP,KNUMER,CLPREF(1:ILPREF),KNIVAU,CLSUFF(1:ILSUFF),LDEXIS, &
&   LDCOSP,KNGRIB,KNBITS,KSTRON,KPUILA
CALL FAIPAR_MT_BB                                    &
&               (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS, &
&                CLNSPR,CLNOMA(1:ILNOMU),LLRLFI)
!
IF (LHOOK) CALL DR_HOOK('FANION_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"
#include "facom2.ixnvms.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FANION_BB                                     &
&           (KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, LDEXIS, &
&           LDCOSP, KNGRIB, KNBITS, KSTRON, KPUILA)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKB)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
LOGICAL                LDEXIS                                 !   OUT
LOGICAL                LDCOSP                                 !   OUT
INTEGER (KIND=JPLIKB)  KNGRIB                                 !   OUT
INTEGER (KIND=JPLIKB)  KNBITS                                 !   OUT
INTEGER (KIND=JPLIKB)  KSTRON                                 !   OUT
INTEGER (KIND=JPLIKB)  KPUILA                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FANION_MT_BB                                            &
&           (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, LDEXIS, &
&           LDCOSP, KNGRIB, KNBITS, KSTRON, KPUILA)

END SUBROUTINE

SUBROUTINE FANION_MM                                     &
&           (KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, LDEXIS, &
&           LDCOSP, KNGRIB, KNBITS, KSTRON, KPUILA)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
LOGICAL                LDEXIS                                 !   OUT
LOGICAL                LDCOSP                                 !   OUT
INTEGER (KIND=JPLIKM)  KNGRIB                                 !   OUT
INTEGER (KIND=JPLIKM)  KNBITS                                 !   OUT
INTEGER (KIND=JPLIKM)  KSTRON                                 !   OUT
INTEGER (KIND=JPLIKM)  KPUILA                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FANION_MT_MM                                            &
&           (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, LDEXIS, &
&           LDCOSP, KNGRIB, KNBITS, KSTRON, KPUILA)

END SUBROUTINE

SUBROUTINE FANION_MT_MM                                      &
&           (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, LDEXIS, &
&           LDCOSP, KNGRIB, KNBITS, KSTRON, KPUILA)
USE FA_MOD, ONLY : FA_COM,              &
&                   FA_COM_DEFAULT_INIT, &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
LOGICAL                LDEXIS                                 !   OUT
LOGICAL                LDCOSP                                 !   OUT
INTEGER (KIND=JPLIKM)  KNGRIB                                 !   OUT
INTEGER (KIND=JPLIKM)  KNBITS                                 !   OUT
INTEGER (KIND=JPLIKM)  KSTRON                                 !   OUT
INTEGER (KIND=JPLIKM)  KPUILA                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  INIVAU                                 ! IN   
INTEGER (KIND=JPLIKB)  INGRIB                                 !   OUT
INTEGER (KIND=JPLIKB)  INBITS                                 !   OUT
INTEGER (KIND=JPLIKB)  ISTRON                                 !   OUT
INTEGER (KIND=JPLIKB)  IPUILA                                 !   OUT
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)
INIVAU     = INT (    KNIVAU, JPLIKB)

CALL FANION_MT_BB                                            &
&           (FA, IREP, INUMER, CDPREF, INIVAU, CDSUFF, LDEXIS, &
&           LDCOSP, INGRIB, INBITS, ISTRON, IPUILA)

KREP       = INT (      IREP, JPLIKM)
KNGRIB     = INT (    INGRIB, JPLIKM)
KNBITS     = INT (    INBITS, JPLIKM)
KSTRON     = INT (    ISTRON, JPLIKM)
KPUILA     = INT (    IPUILA, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT 
!INTF KNUMER        IN    
!INTF CDPREF        IN    
!INTF KNIVAU        IN    
!INTF CDSUFF        IN    
!INTF LDEXIS          OUT 
!INTF LDCOSP          OUT 
!INTF KNGRIB          OUT 
!INTF KNBITS          OUT 
!INTF KSTRON          OUT 
!INTF KPUILA          OUT 
