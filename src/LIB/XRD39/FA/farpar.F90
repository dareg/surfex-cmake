! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FARPAR_MT_BB                                        &
&                     (FA,  KREP, CDPREF, CDSUFF, KCODPA, KNUM )
USE FA_MOD, ONLY : FA_COM, JPNIIL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!      Sous-programme de reglage de la correspondance "nom d'article FA"
!              <-> "descripteurs GRIB du parametre+niveau"
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!    ( Tableau ) CDPREF (Entree) ==> Prefixe pour les KNUM noms d'article;
!    ( Tableau ) CDSUFF (Entree) ==> Suffixe pour les KNUM noms d'article;
!    ( Tableau ) KCODPA (Entree) ==> 6 descripteurs GRIB pour chacun
!                                    des KNUM parametres:
!      KCODPA(J,1) = KSEC1(1)  version de la table parametres
!      KCODPA(J,2) = KSEC1(6)  indicateur du parametre
!      KCODPA(J,3) = KSEC1(7)  indicateur du type de niveau
!      KCODPA(J,4) = KSEC1(8)  niveau
!      KCODPA(J,5) = KSEC1(9)  2ieme nv si couche, sinon 0
!      KCODPA(J,6) = KSEC1(18) indicateur du type de champ
!                              (0 sf si min/max:2 ou si cumul:4)
!
!                KNUM   (Entree  ==> Nombre de parametres a regler
!                          et        (dimension de CDPREF, CDSUFF et KCODPA)
!                        Sortie) ==> Nb de nouveaux parametres pouvant encore
!                                    etre definis lors d'un appel ulterieur.
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KNUM
INTEGER (KIND=JPLIKB) KCODPA(KNUM,6)
!
CHARACTER (LEN=*) CDPREF(KNUM), CDSUFF(KNUM)
!
INTEGER (KIND=JPLIKB) J, JJ, INUMER, INIMES, JMEM
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR

!
INTRINSIC LEN_TRIM
!
!**
!     0.  -  CONTROLES ET INITIALISATIONS PREALABLES
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FARPAR_MT',0,ZHOOK_HANDLE)
IF (KNUM.LT.1) THEN
  KREP=-129
  IF (FA%LFAMOP) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                         &
&           'FARPAR: Nb de parametres ',KNUM,' incorrect'
  ENDIF
  GOTO 1001
ENDIF
DO J = 1,KNUM
  IF ( INT (LEN_TRIM(CDPREF(J)), JPLIKB).LE.0 .OR.           &
&       INT (LEN_TRIM(CDPREF(J)), JPLIKB).GT.FA%JPXPRF ) THEN
    KREP=-129
    IF (FA%LFAMOP) THEN
      WRITE (UNIT=FA%NULOUT,FMT=*)                               &
&              'FARPAR: Longueur du prefixe ',CDPREF(J),          &
&              ' incorrecte : ',INT (LEN_TRIM(CDPREF(J)), JPLIKB)
    ENDIF
    GOTO 1001
  ENDIF
  IF ( INT (LEN_TRIM(CDSUFF(J)), JPLIKB).LE.0 .OR.     &
&       INT (LEN_TRIM(CDSUFF(J)), JPLIKB).GT.FA%JPXNOM- &
&       INT (LEN_TRIM(CDPREF(J)), JPLIKB) ) THEN
    KREP=-129
    IF (FA%LFAMOP) THEN
      WRITE (UNIT=FA%NULOUT,FMT=*)                               &
&              'FARPAR: Longueur du suffixe ',CDSUFF(J),          &
&              ' incorrecte : ',INT (LEN_TRIM(CDSUFF(J)), JPLIKB)
    ENDIF
    GOTO 1001
  ENDIF
  DO JJ = 1,3
    IF (KCODPA(J,JJ).LT.1 .OR. KCODPA(J,JJ).GT.255) THEN
      KREP=-129
      IF (FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)                    &
&                'FARPAR: descripteur GRIB num ',JJ,     &
&                ' pour le parametre num ',J,' ( ',      &
&                CDPREF(J)//CDSUFF(J),' ) incorrect : ', &
&                KCODPA(J,JJ)
      ENDIF
      GOTO 1001
    ENDIF
  ENDDO
  IF (KCODPA(J,6).LT.0 .OR. KCODPA(J,6).GT.124) THEN
    KREP=-129
    IF (FA%LFAMOP) THEN
      WRITE (UNIT=FA%NULOUT,FMT=*)                    &
&              'FARPAR: descripteur GRIB, KSEC1(18),', &
&              ' pour le parametre num ',J,' ( ',      &
&              CDPREF(J)//CDSUFF(J),' ) incorrect : ', &
&              KCODPA(J,6)
    ENDIF
    GOTO 1001
  ENDIF
  IF (KCODPA(J,4).LT.0) THEN
    KREP=-129
    IF (FA%LFAMOP) THEN
      WRITE (UNIT=FA%NULOUT,FMT=*)                    &
&              'FARPAR: descripteur GRIB, KSEC1(8),',  &
&              ' pour le parametre num ',J,' ( ',      &
&              CDPREF(J)//CDSUFF(J),' ) incorrect : ', &
&              KCODPA(J,4)
    ENDIF
    GOTO 1001
  ENDIF
ENDDO
!
!**
!     2.  -  Prise en compte des nouvelles correspondances
!---------------------------------------------------------
!
!
DO J = 1,KNUM
! Recherche prealable de l'eventuelle existence de la definition
! de ce parametre (il faudra alors l'ecraser).
  JMEM = 0
  DO JJ = 1,FA%NBPARC
    IF (CDPREF(J )(1:INT (LEN_TRIM(CDPREF(J )), JPLIKB)).EQ.     &
&    FA%CIPREF(JJ)(1:INT (LEN_TRIM(FA%CIPREF(JJ)), JPLIKB)) .AND. &
&    CDSUFF(J )(1:INT (LEN_TRIM(CDSUFF(J )), JPLIKB)).EQ.         &
&    FA%CISUFF(JJ)(1:INT (LEN_TRIM(FA%CISUFF(JJ)), JPLIKB))) THEN
      JMEM = JJ
      EXIT
    ENDIF
  ENDDO
  IF (JMEM==0) THEN
    IF (FA%NBPARC+1.GT.FA%JPXPAR) THEN
      KREP=-129
      IF (FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)                              &
&                'FARPAR: le nb lim de parametres definissables ', &
&                FA%JPXPAR
        WRITE (UNIT=FA%NULOUT,FMT=*)                          &
&                '        est depasse avec la definition de ', &
&                KNUM,' nouveaux parametres'
      ENDIF
      GOTO 1001
    ENDIF
    JMEM = FA%NBPARC + 1
    FA%NBPARC = FA%NBPARC + 1
  ENDIF
  FA%CIPREF(JMEM)     =                                 &
&         CDPREF(J)(1:INT (LEN_TRIM(CDPREF(J)), JPLIKB))
  FA%CISUFF(JMEM)     =                                 &
&         CDSUFF(J)(1:INT (LEN_TRIM(CDSUFF(J)), JPLIKB))
  FA%NCODPA(JMEM,1:6) = KCODPA(J,1:6)
  IF (FA%LFAMOP) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                              &
&           'FARPAR: Prise en compte de ',CDPREF(J)//CDSUFF(J)
    WRITE (UNIT=FA%NULOUT,FMT=*)                      &
&           '        associe a KSEC1(1,6:9 et 18) = ', &
&           FA%NCODPA(JMEM,1:6)
  ENDIF
ENDDO
KNUM = FA%JPXPAR - FA%NBPARC
!
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
!
IF (FA%LFAMOP) THEN
  INIMES=2
  CLNSPR='FARPAR'
  INUMER=JPNIIL
!
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I4)') KREP
  CALL FAIPAR_MT_BB                                     &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CLNSPR,.FALSE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FARPAR_MT',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FARPAR_BB                           &
&           (KREP, CDPREF, CDSUFF, KCODPA, KNUM)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KNUM                                   ! INOUT
CHARACTER (LEN=*)      CDPREF     (KNUM)                      ! IN   
CHARACTER (LEN=*)      CDSUFF     (KNUM)                      ! IN   
INTEGER (KIND=JPLIKB)  KCODPA     (KNUM,6)                    ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FARPAR_MT_BB                                  &
&           (FA, KREP, CDPREF, CDSUFF, KCODPA, KNUM)

END SUBROUTINE

SUBROUTINE FARPAR_MM                           &
&           (KREP, CDPREF, CDSUFF, KCODPA, KNUM)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUM                                   ! INOUT
CHARACTER (LEN=*)      CDPREF     (KNUM)                      ! IN   
CHARACTER (LEN=*)      CDSUFF     (KNUM)                      ! IN   
INTEGER (KIND=JPLIKM)  KCODPA     (KNUM,6)                    ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FARPAR_MT_MM                                  &
&           (FA, KREP, CDPREF, CDSUFF, KCODPA, KNUM)

END SUBROUTINE

SUBROUTINE FARPAR_MT_MM                            &
&           (FA, KREP, CDPREF, CDSUFF, KCODPA, KNUM)
USE FA_MOD, ONLY : FA_COM,              &
&                   FA_COM_DEFAULT_INIT, &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUM                                   ! INOUT
CHARACTER (LEN=*)      CDPREF     (KNUM)                      ! IN   
CHARACTER (LEN=*)      CDSUFF     (KNUM)                      ! IN   
INTEGER (KIND=JPLIKM)  KCODPA     (KNUM,6)                    ! IN   
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  ICODPA     (KNUM,6)                    ! IN   
INTEGER (KIND=JPLIKB)  INUM                                   ! INOUT
! Convert arguments

ICODPA     = INT (    KCODPA, JPLIKB)
INUM       = INT (      KNUM, JPLIKB)

CALL FARPAR_MT_BB                                  &
&           (FA, IREP, CDPREF, CDSUFF, ICODPA, INUM)

KREP       = INT (      IREP, JPLIKM)
KNUM       = INT (      INUM, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                               
!INTF CDPREF        IN    DIMS=KNUM                     
!INTF CDSUFF        IN    DIMS=KNUM                     
!INTF KCODPA        IN    DIMS=KNUM,6                   
!INTF KNUM          INOUT                               
