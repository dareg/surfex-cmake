! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FAIPAG_MT_BB                                           &
&                     (FA,  KREP,   KRANG,  CDPREF, KNIVAU, CDSUFF, &
&                      KNIPAR )
USE FA_MOD, ONLY : FA_COM, JPNIIL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
!      Initialisation de quelques descripteurs de l'entete Gribex
!      (section 1) relatifs au parametre, a partir du nom FA du champ.
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KRANG  (Entree) ==> Rang de l'unite logique;
!                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
!                KNIVAU (Entree) ==> Niveau vertical eventuel;
!                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
!    ( Tableau ) KNIPAR (Sortie) ==> quelques descripteurs de la section 1 de
!                                    GRIBEX (KNIPAR(1)  =KSEC1(1),
!                                            KNIPAR(2:5)=KSEC1(6:9) ) et un
!                                    indicateur de type de champ (KNIPAR(6)=
!                                    KSEC1(18)):0->RAS; 2->min/max; 4->cumul
!
!     Original  : 06 juillet 2004, Denis Paradis DSI/DEV
!     --------
!
!     Modifications
!     -------------
!       R. El Ouaraini: 03-Oct-2006, enlever le commentaire de l'initialisation de KNIPAR(3) pour
!                  les types de niveaux :  hauteur, iso-tourb potent, isentrope et modele.
!
!*
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KRANG, KNIVAU, KNIPAR(6)
!
CHARACTER (LEN=*) CDPREF, CDSUFF
!
INTEGER (KIND=JPLIKB) INUMER, INIMES, J, JMEM
!
INTRINSIC LEN_TRIM
!
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     0.  -  INITIALISATIONS PREALABLES
!-----------------------------------------------------------------------
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FAIPAG_MT',0,ZHOOK_HANDLE)
KREP=0
IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
  KREP=-66
  GOTO 1001
ENDIF
!
!  DEFAUTS:
!
! Numero de version de la table de code parametres
KNIPAR(1)=1
! Indicateur de parametre (255=> valeur manquante)
KNIPAR(2)=255
! Indicateur de type de niveau (1=> surface)
KNIPAR(3)=1
! Niveau 1, Niveau 2 et type de champs
KNIPAR(4:6)=0
!**
!     1.  -  UTILISATION DE LA TABLE DE CORRESPONDANCE
!-----------------------------------------------------------------------
!
!
JMEM = 0
DO J = 1,FA%NBPARC
  IF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.               &
&      FA%CIPREF(J)(1:LEN_TRIM(FA%CIPREF(J))) .AND. &
&      CDSUFF(1:LEN_TRIM(CDSUFF)).EQ.               &
&      FA%CISUFF(J)(1:LEN_TRIM(FA%CISUFF(J)))) THEN
    JMEM = J
    EXIT
  ENDIF
ENDDO
IF (FA%LFAMOP.AND.JMEM.EQ.0) THEN
  WRITE (UNIT=FA%NULOUT,FMT=*)                                  &
&         'FAIPAG: WARNING, pas de reference GRIB pour ',        &
&         CDPREF(1:LEN_TRIM(CDPREF))//CDSUFF(1:LEN_TRIM(CDSUFF))
  WRITE (UNIT=FA%NULOUT,FMT=*)'       Les defauts seront utilises'
  GOTO 1001
ENDIF
IF (JMEM.NE.0) KNIPAR(1:6) = FA%NCODPA(JMEM,1:6)
!**
!     2.  -  INITIALISATION DU NIVEAU (AUTRE QUE 0)
!--------------------------------------------------
!
!     2.1 -  Champs sur un niveau isobare
!
IF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.'P') THEN
!  La pression est sur 5 chiffres: on la ramene a l'hPa
!     et on recree le niveau 1000 hPa
!
!  Si KNIVAU < 100, la pression fait moins d'un hPa et
!  on utilise une extension du GRIB introduite par le CEP:
!       KSEC1(7) = 210 (au lieu de 100)
!    et KSEC1(8) = pression en Pa
!
  IF (KNIVAU .GE. 100) THEN
    KNIPAR(3)=100
    KNIPAR(4)=KNIVAU/100
  ELSEIF (KNIVAU==0) THEN
    KNIPAR(3)=100
    KNIPAR(4)=1000
  ELSE
    KNIPAR(3)=210
    KNIPAR(4)=KNIVAU
  ENDIF
!
!     2.2 -  Champs sur un niveau hauteur
!
ELSEIF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.'H') THEN
 KNIPAR(3)=105
  KNIPAR(4)=KNIVAU
!
!     2.3 -  Champs sur un niveau iso-tourbillon-potentiel
!
!            ( unite SI = K m2 s-1 kg-1 = 10+6 PVU
!              mais l'unite retenu est le milliPVU: 10-9 SI)
ELSEIF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.'V') THEN
 KNIPAR(3)=117
! KNIVAU est exprime en 1/10 PVU
  KNIPAR(4)=KNIVAU*100
  IF (KNIVAU==0) KNIPAR(4)=1000
!
!     2.4 -  Champs sur un niveau isentrope
!
ELSEIF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.'T') THEN
 KNIPAR(3)=113
  KNIPAR(4)=KNIVAU
!
!     2.5 -  Champs sur un niveau modele
!
ELSEIF (CDPREF(1:LEN_TRIM(CDPREF)).EQ.'S') THEN
 KNIPAR(3)=109
  KNIPAR(4)=KNIVAU
ENDIF
!
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
LLFATA=LLMOER (KREP,KRANG)
!
IF (FA%LFAMOP.OR.LLFATA) THEN
  INIMES=2
  CLNSPR='FAIPAG'
  INUMER=JPNIIL
!
  WRITE (UNIT=FA%NULOUT,FMT=*)'FAIPAG: KNIPAR(1:6) = ',KNIPAR(1:6)
  WRITE (UNIT=FA%NULOUT,FMT=*)
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4, &
&       '', CDPREF='''''',A,'''''', KNIVAU='',I6,         &
&       '', CDSUFF='''''',A,'''''''')')                   &
&     KREP, KRANG, CDPREF(1:LEN_TRIM(CDPREF)), KNIVAU,    &
&     CDSUFF(1:LEN_TRIM(CDSUFF))
  CALL FAIPAR_MT_BB                                     &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CDPREF,.FALSE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FAIPAG_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FAIPAG_BB                                    &
&           (KREP, KRANG, CDPREF, KNIVAU, CDSUFF, KNIPAR)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKB)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
INTEGER (KIND=JPLIKB)  KNIPAR     (6)                         !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAIPAG_MT_BB                                           &
&           (FA, KREP, KRANG, CDPREF, KNIVAU, CDSUFF, KNIPAR)

END SUBROUTINE

SUBROUTINE FAIPAG_MM                                    &
&           (KREP, KRANG, CDPREF, KNIVAU, CDSUFF, KNIPAR)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIPAR     (6)                         !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAIPAG_MT_MM                                           &
&           (FA, KREP, KRANG, CDPREF, KNIVAU, CDSUFF, KNIPAR)

END SUBROUTINE

SUBROUTINE FAIPAG_MT_MM                                     &
&           (FA, KREP, KRANG, CDPREF, KNIVAU, CDSUFF, KNIPAR)
USE FA_MOD, ONLY : FA_COM,              &
&                   FA_COM_DEFAULT_INIT, &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIPAR     (6)                         !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  IRANG                                  ! IN   
INTEGER (KIND=JPLIKB)  INIVAU                                 ! IN   
INTEGER (KIND=JPLIKB)  INIPAR     (6)                         !   OUT
! Convert arguments

IRANG      = INT (     KRANG, JPLIKB)
INIVAU     = INT (    KNIVAU, JPLIKB)

CALL FAIPAG_MT_BB                                           &
&           (FA, IREP, IRANG, CDPREF, INIVAU, CDSUFF, INIPAR)

KREP       = INT (      IREP, JPLIKM)
KNIPAR     = INT (    INIPAR, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                               
!INTF KRANG         IN                                  
!INTF CDPREF        IN                                  
!INTF KNIVAU        IN                                  
!INTF CDSUFF        IN                                  
!INTF KNIPAR          OUT DIMS=6                        
