! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe LFI
SUBROUTINE LFINAF_MT_BB                                  &
&                     (LFI, KREP, KNUMER, KNALDO, KNTROU,  &
&                      KNARES, KNAMAX )
USE LFIMOD, ONLY : LFICOM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!        SOUS-PROGRAMME DONNANT DES NOMBRES D'ARTICLES LOGIQUES DIVERS
!     ( DE DONNEES, "TROUS", POSSIBLE, MAXIMUM ) POUR UNE UNITE LOGIQUE
!     OUVERTE POUR LE LOGICIEL DE FICHIERS INDEXES *LFI* .
!**
!    ARGUMENTS : KREP   (SORTIE) ==> CODE-REPONSE DU SOUS-PROGRAMME;
!                KNUMER (ENTREE) ==> LFI%NUMERO DE L'UNITE LOGIQUE;
!                KNALDO (SORTIE) ==> NOMBRE D'ARTICLES LOGIQUES
!                                  *DE DONNEES* (TROUS EXCLUS) PRESENTS;
!                KNTROU (SORTIE) ==> NOMBRE D'ARTICLES LOGIQUES QUI SONT
!                                    DEVENUS DES "TROUS", PAR SUITE DE
!                                    REECRITURES D'ARTICLES PLUS LONGUES
!                                    QU'INITIALEMENT, ET QUI N'ONT PAS
!                                    PU ETRE (ENCORE) RECYCLES;
!                KNARES (SORTIE) ==> NOMBRE D'ARTICLES LOGIQUES POUVANT
!                                    ETRE ECRITS DANS LA PARTIE "PRERER-
!                                    VEE" DE L'INDEX (TROUS COMPRIS);
!                KNAMAX (SORTIE) ==> NOMBRE D'ARTICLES LOGIQUES MAXIMUM
!                                    POUVANT ETRE ECRITS SUR L'UNITE
!                                    LOGIQUE, EN "DEBORDANT" AU MAXIMUM
!                                   DE LA PARTIE PRERESERVEE DE L'INDEX.
!                                    ( TROUS COMPRIS )
!
!
TYPE(LFICOM) :: LFI
INTEGER (KIND=JPLIKB) KREP, KNUMER, KNALDO, KNTROU, KNARES, KNAMAX
INTEGER (KIND=JPLIKB) IRANG, IREP, IFACTM, ILARPH, INALPP, INIMES
!
CHARACTER(LEN=LFI%JPLSPX) CLNSPR
CHARACTER(LEN=LFI%JPLMES) CLMESS
CHARACTER(LEN=LFI%JPLFTX) CLACTI
LOGICAL LLFATA

!**
!     1.  -  CONTROLES DU PARAMETRE D'APPEL, PUIS INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('LFINAF_MT',0,ZHOOK_HANDLE)
CLACTI=''
CALL LFINUM_MT_BB                 &
&               (LFI, KNUMER,IRANG)
!
IF (IRANG.EQ.0) THEN
  KNTROU=0
  KNALDO=0
  KNARES=0
  KNAMAX=0
  IREP=-1
  GOTO 1001
ENDIF
!
 IF (LFI%LMULTI) CALL LFIVER_MT_BB                           &
&                                (LFI, LFI%VERRUE(IRANG),'ON')
!**
!     2.  -  CALCUL DIRECT DES ARGUMENTS DE SORTIE DU SOUS-PROGRAMME.
!-----------------------------------------------------------------------
!
IFACTM=LFI%MFACTM(IRANG)
ILARPH=LFI%JPLARD*IFACTM
INALPP=LFI%JPNAPP*IFACTM
KNTROU=LFI%MDES1D(IXM(LFI%JPNTRU,IRANG))+LFI%NBTROU(IRANG)
KNALDO=LFI%MDES1D(IXM(LFI%JPNALO,IRANG))-KNTROU
KNARES=INALPP*LFI%MDES1D(IXM(LFI%JPNPIR,IRANG))
KNAMAX=KNARES+INALPP*(ILARPH-LFI%JPLDOC)
IREP=0
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "LFIEMS" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
KREP=IREP
LLFATA=LLMOER (IREP,IRANG)
!
IF (IRANG.NE.0) THEN
  LFI%NDEROP(IRANG)=12
  LFI%NDERCO(IRANG)=IREP
   IF (LFI%LMULTI) CALL LFIVER_MT_BB                            &
&                                  (LFI, LFI%VERRUE(IRANG),'OFF')
ENDIF
!
IF (LLFATA.OR.IXNIMS (IRANG).EQ.2) THEN
  INIMES=2
ELSE
  IF (LHOOK) CALL DR_HOOK('LFINAF_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='LFINAF'
WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,            &
&       '', KNALDO='',I6,'', KNTROU='',I6,'', KNARES='',I6,         &
&       '', KNAMAX='',I6)') KREP,KNUMER,KNALDO,KNTROU,KNARES,KNAMAX
CALL LFIEMS_MT_BB                              &
&               (LFI, KNUMER,INIMES,IREP,LLFATA, &
&                CLMESS,CLNSPR,CLACTI)
!
IF (LHOOK) CALL DR_HOOK('LFINAF_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "lficom2.ixm.h"
#include "lficom2.ixnims.h"
#include "lficom2.llmoer.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE LFINAF_BB                                     &
&           (KREP, KNUMER, KNALDO, KNTROU, KNARES, KNAMAX)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  KNALDO                                 !   OUT
INTEGER (KIND=JPLIKB)  KNTROU                                 !   OUT
INTEGER (KIND=JPLIKB)  KNARES                                 !   OUT
INTEGER (KIND=JPLIKB)  KNAMAX                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFINAF_MT_BB                                             &
&           (LFI, KREP, KNUMER, KNALDO, KNTROU, KNARES, KNAMAX)

END SUBROUTINE

SUBROUTINE LFINAF_MM                                     &
&           (KREP, KNUMER, KNALDO, KNTROU, KNARES, KNAMAX)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KNALDO                                 !   OUT
INTEGER (KIND=JPLIKM)  KNTROU                                 !   OUT
INTEGER (KIND=JPLIKM)  KNARES                                 !   OUT
INTEGER (KIND=JPLIKM)  KNAMAX                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFINAF_MT_MM                                             &
&           (LFI, KREP, KNUMER, KNALDO, KNTROU, KNARES, KNAMAX)

END SUBROUTINE

SUBROUTINE LFINAF_MT_MM                                       &
&           (LFI, KREP, KNUMER, KNALDO, KNTROU, KNARES, KNAMAX)
USE LFIMOD, ONLY : LFICOM,              &
&                   LFICOM_DEFAULT_INIT, &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (LFICOM)          LFI                                    ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KNALDO                                 !   OUT
INTEGER (KIND=JPLIKM)  KNTROU                                 !   OUT
INTEGER (KIND=JPLIKM)  KNARES                                 !   OUT
INTEGER (KIND=JPLIKM)  KNAMAX                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  INALDO                                 !   OUT
INTEGER (KIND=JPLIKB)  INTROU                                 !   OUT
INTEGER (KIND=JPLIKB)  INARES                                 !   OUT
INTEGER (KIND=JPLIKB)  INAMAX                                 !   OUT
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)

CALL LFINAF_MT_BB                                             &
&           (LFI, IREP, INUMER, INALDO, INTROU, INARES, INAMAX)

KREP       = INT (      IREP, JPLIKM)
KNALDO     = INT (    INALDO, JPLIKM)
KNTROU     = INT (    INTROU, JPLIKM)
KNARES     = INT (    INARES, JPLIKM)
KNAMAX     = INT (    INAMAX, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT 
!INTF KNUMER        IN    
!INTF KNALDO          OUT 
!INTF KNTROU          OUT 
!INTF KNARES          OUT 
!INTF KNAMAX          OUT 
