! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
! Sep-2012 P. Marguinaud Fix unitialized variables
SUBROUTINE FADECI_MT64                                            &
&                     (FA,  KREP,   KRANG,  CDNOMA, KVALCO, KLONGA, &
&                    KCHAMP, LDCOSP )
USE FA_MOD, ONLY : FA_COM, JPNIIL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
!      Controle de coherence et decodage d'un CHAMP HORIZONTAL
!      venant d'etre lu sur un fichier ARPEGE/ALADIN.
!       ( DECodage Interne d'un champ lu )
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KRANG  (Entree) ==> Rang de l'unite logique;
!                CDNOMA (Entree) ==> Nom d'article (prefabrique);
!    ( Tableau ) KVALCO (Entree) ==> Donnees issues de la lecture;
!                KLONGA (Entree) ==> Nombre de mots lus;
!    ( Tableau ) KCHAMP (Sortie) ==> Valeurs REELLES du champ lu;
!                LDCOSP (Entree) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux;
!*
!       En mode multi-taches, il doit y avoir verrouillage du fichier
!     concerne avant l'appel au sous-programme.
!
!     Modifications
!     -------------
!
!  Avril 2004, D. Paradis, DSI/DEV:
!
!    -Declaration ICHAUX en ALLOCATABLE, KCHAMP en profil implicite (gain mem.)
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KRANG, KLONGA
!
INTEGER (KIND=JPDBLE) KVALCO(*), KCHAMP(*)
!
LOGICAL LDCOSP
!
CHARACTER CDNOMA*(*)
!
REAL (KIND=JPDBLR) ZFOHYB (2)
!
INTEGER (KIND=JPLIKB) ILCHAM, ISTRIA, J, IDECAL, ICPACK
INTEGER (KIND=JPLIKB) IPUILA, IPOFIN
INTEGER (KIND=JPLIKB) ITRONC, IIND, ILOW, IHIGH, JTRON
INTEGER (KIND=JPLIKB) IDIMNC, INBITS
INTEGER (KIND=JPLIKB) IL, IADD, IRANGC, IILCHAM, INDECO
INTEGER (KIND=JPLIKB) IERR, INIMES
INTEGER (KIND=JPLIKB) IVALC3, IVALC4, IVALC5, IJLENV
INTEGER (KIND=JPLIKB) IJLENF, IDIZAI, IUNITE
INTEGER (KIND=JPLIKB) INUMER
!
INTEGER (KIND=JPDBLE), ALLOCATABLE :: ICHAUX(:)
INTEGER (KIND=JPLIKB) IB1PAR (FA%JPLB1P), IB2PAR (FA%JPLB2P)
!
LOGICAL LLARPE, LLMLAM, LLCOSP
!
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FADECI_MT',0,ZHOOK_HANDLE)
KREP=0
IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
  KREP=-66
  GOTO 1001
ENDIF
!
ISTRIA=0
INBITS=0
ICPACK=0
IPUILA=0
IVALC3=0
IVALC4=0
IVALC5=0
!**
!     2.  -  CONTROLE DES DONNEES DE L'ARTICLE
!-----------------------------------------------------------------------
!
IF (KVALCO(1).LT.-1.OR.KVALCO(1).GT.2.OR.                        &
&    KVALCO(2).LT.0.OR.KVALCO(2).GT.1.OR.                         &
&    (KVALCO(1).GT.0.AND.KVALCO(2).EQ.1.AND.KVALCO(4).LT.0)) THEN
  KREP=-91
  GOTO 1001
ELSE
  LLARPE=KVALCO(1).EQ.2
  LLCOSP=KVALCO(2).EQ.1
ENDIF
!
IF ((LLCOSP.AND..NOT.LDCOSP).OR.(.NOT.LLCOSP.AND.LDCOSP)) THEN
  KREP=-92
  GOTO 1001
ENDIF
!
IRANGC=FA%FICHIER(KRANG)%NUCADR
LLMLAM=FA%CADRE(IRANGC)%LIMLAM
!
IF (LDCOSP) THEN
  IF (LLMLAM) THEN
     ILCHAM=FA%CADRE(IRANGC)%NSFLAM
  ELSE    
    IF (KVALCO(1).EQ.-1) THEN
      ILCHAM=(1+FA%CADRE(IRANGC)%MTRONC)*(2+FA%CADRE(IRANGC)%MTRONC)
    ELSE
      ILCHAM=(1+FA%CADRE(IRANGC)%MTRONC)**2
    ENDIF
  ENDIF   
ELSE
  ILCHAM=FA%CADRE(IRANGC)%NVAPDG
ENDIF
!
!**
!     3.  -  DECODAGE DES DONNEES DE L'ARTICLE
!-----------------------------------------------------------------------
!
IF (KVALCO(1).EQ.-1 .OR. KVALCO(1).EQ.0) THEN
!
!          Cas ou il n'y a aucun codage... controle longueur d'article
!
  IF (KLONGA.LT.(ILCHAM+2)) THEN
    KREP=-93
    GOTO 1001
  ELSEIF (KLONGA.GT.(ILCHAM+2)) THEN
    KREP=-94
    IF (LLMOER(KREP,KRANG)) GOTO 1001
  ENDIF
!
!           Transfert du tableau d'entree a la suite des 2 valeurs
!     documentaires stockees en debut d'article.
!
  DO J=1,ILCHAM
  KCHAMP(J)=KVALCO(2+J)
  ENDDO
!
ELSE
!*
!     3.1 -  DECODAGE GRIB PROPREMENT DIT (STANDARD OU NON).
!-----------------------------------------------------------------------
!
  IDECAL=1+2*KVALCO(1)
  IF (LDCOSP) IDECAL=IDECAL+2
  IVALC3=KVALCO(3)
  IVALC4=KVALCO(4)
  IVALC5=KVALCO(5)
  IF (LDCOSP.AND.LLMLAM) THEN
!
    ALLOCATE (ICHAUX (ILCHAM))
!
    ITRONC=FA%CADRE(IRANGC)%MTRONC
    ISTRIA=FA%CADRE(IRANGC)%NOZPAR(4)-FA%CADRE(IRANGC)%NOZPAR(3)+1
    DO JTRON=1,ITRONC 
      IADD=4*(IVALC4+1-JTRON)
      IF (IADD.LE.0) IADD=4
      ISTRIA=ISTRIA+IADD
    ENDDO
    IILCHAM=ILCHAM-ISTRIA
    CALL FADECOGA(ICHAUX,IILCHAM,INBITS,FA%NBIMAC,IB1PAR,IB2PAR, &
&                 ZFOHYB(1),2_JPLIKB ,KVALCO(IDECAL+1),           &
&                 KLONGA-IDECAL,INDECO,IJLENV,IJLENF,ICPACK,      &
&                 IPUILA,IERR,KVALCO(IDECAL-1),KVALCO(IDECAL),    &
&                 LLARPE)
!
!  Controle de l'adequation entre nb de valeurs attendues/lues
!
    IF (IJLENF.LT.IILCHAM) THEN
      KREP=-93
      IF (FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)                               &
&                 'FADECI: erreur !!! Nbre de valeurs decodees = ', &
&                 IJLENF,' et nbre de valeurs attendues = ',IILCHAM
      ENDIF
      GOTO 1001
    ELSEIF (IJLENF.GT.IILCHAM) THEN
      KREP=-94
      IF (FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)                               &
&                 'FADECI: erreur !!! Nbre de valeurs decodees = ', &
&                 IJLENF,' et nbre de valeurs attendues = ',IILCHAM
      ENDIF
      IF (LLMOER(KREP,KRANG)) GOTO 1001
    ENDIF
    IIND=0
    DO JTRON=1,ITRONC      
    ILOW=2+2*JTRON+1
    IADD=4*(IVALC4+1-JTRON)
    IF (IADD.LE.0) IADD=4  
    DO J=FA%CADRE(IRANGC)%NOZPAR(ILOW)+IADD,FA%CADRE(IRANGC)%NOZPAR(ILOW+1)
    IIND=IIND+1        
    KCHAMP(J)=ICHAUX(IIND)
    ENDDO
    ENDDO
!
    IF (ALLOCATED( ICHAUX )) DEALLOCATE ( ICHAUX )
!
  ELSE                         
    CALL FADECOGA (KCHAMP,ILCHAM,INBITS,FA%NBIMAC,IB1PAR,IB2PAR, &
&                  ZFOHYB(1),2_JPLIKB ,KVALCO(IDECAL+1),          &
&                  KLONGA-IDECAL,INDECO,IJLENV,IJLENF,ICPACK,     &
&                  IPUILA,IERR,KVALCO(IDECAL-1),KVALCO(IDECAL),   &
&                  LLARPE)
!
!  Controle de l'adequation entre nb de valeurs attendues/lues
!
    IF (IJLENF.LT.ILCHAM) THEN
      KREP=-93
      IF (FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)                               &
&                 'FADECI: erreur !!! Nbre de valeurs decodees = ', &
&                 IJLENF,' et nbre de valeurs attendues = ',ILCHAM
      ENDIF
      GOTO 1001
    ELSEIF (IJLENF.GT.ILCHAM) THEN
      KREP=-94
      IF (FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)                               &
&                 'FADECI: erreur !!! Nbre de valeurs decodees = ', &
&                 IJLENF,' et nbre de valeurs attendues = ',ILCHAM
      ENDIF
      IF (LLMOER(KREP,KRANG)) GOTO 1001
    ENDIF
  ENDIF
!
  IF (IERR.EQ.-2) THEN
    KREP=-93
    GOTO 1001
  ELSEIF (IERR.NE.0) THEN
    KREP=-200+IERR
    GOTO 1001
  ELSEIF (IVALC3.NE.INBITS.OR.(LDCOSP.AND.        &
&          ((ICPACK.NE.IVALC4.AND..NOT.LLMLAM)     &
&   .OR.(.NOT.LLMLAM.AND.IPUILA.NE.IVALC5)))) THEN                 
    KREP=-95
    GOTO 1001
  ELSEIF (IB1PAR(4).GT.64) THEN
!
!     Controle effectue s'il y a un bloc 2 en retour du decodage.
!
    IDIZAI=IB2PAR(1)/10
    IUNITE=IB2PAR(1)-IDIZAI*10
!
    IF ((LDCOSP.AND..NOT.LLMLAM.AND.                        &
&     (IUNITE.NE.0.OR.IDIZAI.LT.5.OR.IDIZAI.GT.8))           &
&    .OR.(.NOT.LDCOSP.AND.                                   &
&         (IUNITE.NE.4.OR.IDIZAI.LT.0.OR.IDIZAI.GT.3))) THEN
      KREP=-95
      GOTO 1001
    ENDIF
!
  ENDIF
  IF (LDCOSP.AND.LLMLAM) THEN
    ICPACK=IVALC4
    IPUILA=IVALC5
  ENDIF          
!
  IF (LDCOSP) THEN
!
    IF (LLARPE) THEN
      IF (LLMLAM) THEN
        IPOFIN=IDECAL+INDECO+ISTRIA
      ELSE   
        IDIMNC=(1+ICPACK)**2
        IPOFIN=IDECAL+INDECO+IDIMNC
      ENDIF            
!
      IF (KLONGA.LT.IPOFIN) THEN
        KREP=-93
        GOTO 1001
      ELSEIF (KLONGA.GT.IPOFIN) THEN
        KREP=-94
        IF (LLMOER(KREP,KRANG)) GOTO 1001
      ENDIF
!*
!     3.2 -  TRANSFERT DES COEFFICIENTS SPECTRAUX NON COMPACTES.
!-----------------------------------------------------------------------
!        (et non fournis par DECOGA) stockes en fin d'article.
!
      IF (LLMLAM) THEN
        IIND=0
         DO JTRON=0,ITRONC
         IL=2+2*JTRON+1
         ILOW=FA%CADRE(IRANGC)%NOZPAR(IL)
         IF (JTRON.EQ.0) THEN
            IHIGH=FA%CADRE(IRANGC)%NOZPAR(IL+1)
         ELSE
            IHIGH=ILOW+4*(ICPACK+1-JTRON)-1
            IF (IHIGH.LE.ILOW) IHIGH=ILOW+3
         ENDIF     
         DO J=ILOW,IHIGH
         IIND=IIND+1
         KCHAMP(J)=KVALCO(IDECAL+INDECO+IIND)
         ENDDO
         ENDDO
      ELSE
        DO J=1,IDIMNC
        KCHAMP(J)=KVALCO(IDECAL+INDECO+J)
        ENDDO
      ENDIF
!
    ENDIF
!*
!     3.3 -  SI NECESSAIRE, RECONSTITUTION DU SPECTRE.
!-----------------------------------------------------------------------
!
    IF (IPUILA.NE.0) THEN
      CALL FARCIS_MT64                                    &
&                     (FA, KREP,KRANG,KCHAMP,ICPACK,IPUILA)
      IF (KREP.NE.0) GOTO 1001
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
LLFATA=LLMOER (KREP,KRANG)
!
IF (FA%LFAMOP.OR.LLFATA) THEN
  INIMES=2
  CLNSPR='FADECI'
  INUMER=JPNIIL
!
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4, &
&         '', CDNOMA='''''',A,'''''', LDCOSP= '',L1,      &
&         '', KLONGA='',I8)')                             &
&     KREP, KRANG, CDNOMA, LDCOSP, KLONGA
  CALL FAIPAR_MT64                                      &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CDNOMA,.FALSE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FADECI_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"

END SUBROUTINE




! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FADECI64                                     &
&           (KREP, KRANG, CDNOMA, KVALCO, KLONGA, KCHAMP, &
&           LDCOSP)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDNOMA                                 ! IN   
INTEGER (KIND=JPDBLE)  KVALCO     (*)                         ! IN   
INTEGER (KIND=JPLIKB)  KLONGA                                 ! IN   
INTEGER (KIND=JPDBLE)  KCHAMP     (*)                         !   OUT
LOGICAL                LDCOSP                                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FADECI_MT64                                            &
&           (FA, KREP, KRANG, CDNOMA, KVALCO, KLONGA, KCHAMP, &
&           LDCOSP)

END SUBROUTINE

SUBROUTINE FADECI                                       &
&           (KREP, KRANG, CDNOMA, KVALCO, KLONGA, KCHAMP, &
&           LDCOSP)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDNOMA                                 ! IN   
INTEGER (KIND=JPDBLE)  KVALCO     (*)                         ! IN   
INTEGER (KIND=JPLIKM)  KLONGA                                 ! IN   
INTEGER (KIND=JPDBLE)  KCHAMP     (*)                         !   OUT
LOGICAL                LDCOSP                                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FADECI_MT                                              &
&           (FA, KREP, KRANG, CDNOMA, KVALCO, KLONGA, KCHAMP, &
&           LDCOSP)

END SUBROUTINE

SUBROUTINE FADECI_MT                                        &
&           (FA, KREP, KRANG, CDNOMA, KVALCO, KLONGA, KCHAMP, &
&           LDCOSP)
USE FA_MOD, ONLY : FA_COM
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDNOMA                                 ! IN   
INTEGER (KIND=JPDBLE)  KVALCO     (*)                         ! IN   
INTEGER (KIND=JPLIKM)  KLONGA                                 ! IN   
INTEGER (KIND=JPDBLE)  KCHAMP     (*)                         !   OUT
LOGICAL                LDCOSP                                 ! IN   
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  IRANG                                  ! IN   
INTEGER (KIND=JPLIKB)  ILONGA                                 ! IN   
! Convert arguments

IRANG      = INT (     KRANG, JPLIKB)
ILONGA     = INT (    KLONGA, JPLIKB)

CALL FADECI_MT64                                            &
&           (FA, IREP, IRANG, CDNOMA, KVALCO, ILONGA, KCHAMP, &
&           LDCOSP)

KREP       = INT (      IREP, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                                                              
!INTF KRANG         IN                                                                 
!INTF CDNOMA        IN                                                                 
!INTF KVALCO        IN    DIMS=*                         KIND=JPDBLE                   
!INTF KLONGA        IN                                                                 
!INTF KCHAMP          OUT DIMS=*                         KIND=JPDBLE                   
!INTF LDCOSP        IN                                                                 
