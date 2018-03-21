! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FADECX_MT_BB                                       &
&                     (FA, KREP, KRANG, CDNOMA, KVALCO, KLONGA, &
&                      KCHAMP, LDCOSP )
USE FA_MOD, ONLY : FA_COM, JPNIIL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
!      Controle de coherence et decodage (GRIBEX) d'un CHAMP
!      HORIZONTAL venant d'etre lu sur un fichier ARPEGE/ALADIN.
!       ( DECodage d'un champ gribeX )
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
!  MODIFICATION :
!         JM AUDOIN  15/05/2007  Partie 3.1  Blindage controle changement unite
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
!#include "fagribexi.h"
!
REAL (KIND=JPDBLR) ZSEC2(10+2*(FA%JPXNIV+1)), ZSEC3(2)
REAL (KIND=JPDBLR), ALLOCATABLE ::  ZSEC4(:), ZCHAMP(:)
REAL (KIND=JPDBLR) ZPULAP
!
INTEGER (KIND=JPDBLE), ALLOCATABLE :: ICHAMP(:)
!
INTEGER (KIND=JPLIKB) ISEC0(2), ISEC1(FA%JPSEC1)
INTEGER (KIND=JPLIKB) ISEC2(FA%JPSEC2), ISEC3(2)
INTEGER (KIND=JPLIKB) ISEC4(FA%JPSEC4)
INTEGER (KIND=JPLIKB) ILCHAM, ISTRIA, IDECAL
INTEGER (KIND=JPLIKB) IPOFIN, ILONSEC2
INTEGER (KIND=JPLIKB) ITRONC, IIND, ILOW, IHIGH
INTEGER (KIND=JPLIKB) IL, IADD, IRANGC, IILCHAM
INTEGER (KIND=JPLIKB) INIMES
INTEGER (KIND=JPLIKB) IVALC3, IVALC4, IVALC5, IWORD
INTEGER (KIND=JPLIKB) INUMER, ILENG, IRET, IDX
INTEGER (KIND=JPLIKB) JN, JM, JLAT, JLON
INTEGER (KIND=JPLIKB) IFAORI, IFAMOD, INBIMO
INTEGER (KIND=JPLIKB) I7,I10,I16
!
LOGICAL LLMLAM, LLCOSP
!
CHARACTER(LEN=1) CLOPER
!
INTEGER (KIND=JPLIKB) DECF10
EXTERNAL DECF10
!
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FADECX_MT',0,ZHOOK_HANDLE)
KREP=0
IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
  KREP=-66
  GOTO 1001
ENDIF
!
CLOPER='D'
ISTRIA=0
!**
!     2.  -  CONTROLE DES DONNEES DE L'ARTICLE
!-----------------------------------------------------------------------
!
IF (KVALCO(1).NE.3.OR.                        &
&    KVALCO(2).LT.0.OR.KVALCO(2).GT.1.OR.      &
&    (KVALCO(2).EQ.1.AND.KVALCO(4).LT.0)) THEN
  KREP=-91
  GOTO 1001
ELSE
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
ITRONC=FA%CADRE(IRANGC)%MTRONC
!
IF (LDCOSP) THEN
  IF (LLMLAM) THEN
    ILCHAM=FA%CADRE(IRANGC)%NSFLAM
    ILONSEC2=21+ITRONC
  ELSE    
    ILCHAM=(1+ITRONC)*(2+ITRONC)
    ILONSEC2=22
  ENDIF   
ELSE
  ILCHAM=FA%CADRE(IRANGC)%NVAPDG
  IF (LLMLAM) THEN
    ILONSEC2=22
  ELSE
    ILONSEC2=22+FA%CADRE(IRANGC)%NLATIT
  ENDIF
ENDIF
!
ALLOCATE (ICHAMP(ILCHAM), ZCHAMP(ILCHAM), ZSEC4(ILCHAM))
!
!**
!     3.  -  DECODAGE GRIBEX DES DONNEES DE L'ARTICLE
!-----------------------------------------------------------------------
!
IDECAL=3
IVALC3=KVALCO(3)
IF (LDCOSP) THEN
  IDECAL=IDECAL+2
! IVALC4=ss-tronc non compactee
! IVALC5=puissance de laplacien
  IVALC4=KVALCO(4)
  IVALC5=KVALCO(5)
ENDIF
IILCHAM=ILCHAM
!
! Pour Aladin, le calcul du nb de coeff spectraux qui ont
! ete compactes est plus complexe (certains ont ete retires
! pour ne pas etre compactes: ss-tronc triangulaire).
!
IF (LDCOSP.AND.LLMLAM) THEN
  ISTRIA=4*(1+FA%CADRE(IRANGC)%NOZPAR(1)+FA%CADRE(IRANGC)%NOZPAR(2)+ &
&            IVALC4*(IVALC4-1)/2)
  IILCHAM=ILCHAM-ISTRIA
ENDIF
! ILENG=longueur disponible en entiers declares INTEGER dans KVALCO.
ILENG=(KIND(KVALCO)/4)*(KLONGA-IDECAL)
!
!     3.1 -  CHANGEMENT D'UNITE DE CERTAINS CHAMPS.
!            Il s'agit de champs dont les valeurs sont comprises
!            entre 0 et 1 dans le modele mais dont l'unite
!            conventionnelle dans le GRIB est le %.
!            Avant l'appel a GRIBEX, il faut leur redonner leurs
!            valeurs d'origine (comprises entre 0 et 1) en ajoutant
!            2 au facteur d'echelle decimal KSEC1(23) via DECF10.
!
I7 =MIN(7_JPLIKB ,INT (LEN(TRIM(CDNOMA)), JPLIKB))
I10=MIN(10_JPLIKB ,INT (LEN(TRIM(CDNOMA)), JPLIKB))
I16=MIN(16_JPLIKB ,INT (LEN(TRIM(CDNOMA)), JPLIKB))
IF (                                                             &
& CDNOMA(1:I10)=='SURFNEBUL.' .OR.                                &
& CDNOMA(1:I10)=='SURFALBEDO' .OR.                                &
& CDNOMA=='SURFPROP.VEGETAT' .OR.                                 &
& CDNOMA=='CLSHUMI.RELATIVE' .OR. CDNOMA=='CLSMAXI.HUMI.REL' .OR. &
& CDNOMA=='CLSMINI.HUMI.REL' .OR.                                 &
& (CDNOMA(1:1)=='P'.AND.CDNOMA(I7:I16)=='HUMI_RELAT').OR.         &
& (CDNOMA(1:1)=='H'.AND.CDNOMA(I7:I16)=='HUMI_RELAT')      ) THEN
  IADD = 2
  INBIMO = 32  ! Nombre de BIts par mot (un mot=INTEGER)
  IDX = DECF10 ( KVALCO(IDECAL+1), ILENG, IADD, &
&                   IFAORI, IFAMOD, INBIMO )
  IF (IDX==-1) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                   &
&           'FADECX: pas d''entete GRIB au debut !'
    KREP=-128
    GOTO 1001
  ELSEIF (IDX==-2) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                           &
&           'FADECX: edition du GRIB invalid pour DECF10 !'
    KREP=-128
    GOTO 1001
  ELSEIF (IDX > 0) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                             &
&           'FADECX: ERREUR dans appel a INXBIT par DECF10 !'
    WRITE (UNIT=FA%NULOUT,FMT=*)                       &
&           'FADECX: avec code retour de INXBIT = ',IDX
    KREP=-128
    GOTO 1001
  ELSEIF (IDX < -2) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                          &
&           'FADECX: code retour inconnu de DECF10 : ',IDX
    KREP=-126
    GOTO 1001
  ENDIF
ENDIF
!
!     3.2 -  APPEL A GRIBEX
!
IWORD=0
IRET=-1
!CALL FAGRIBEX(ISEC0,ISEC1,ISEC2,ZSEC2,ISEC3,ZSEC3,ISEC4,   &
!&              KCHAMP,IILCHAM,KVALCO(IDECAL+1),ILENG,IWORD, &
!&              CLOPER,IRET)
IF (FA%LFAMOP) THEN
  WRITE (UNIT=FA%NULOUT,FMT=*)                         &
&         ' FADECX: KLONGA, IDECAL, ILENG, IILCHAM = ', &
&         KLONGA, IDECAL, ILENG, IILCHAM
  WRITE (UNIT=FA%NULOUT,FMT=*) '       * ISEC0 = ',ISEC0
  WRITE (UNIT=FA%NULOUT,FMT=*) '       * ISEC1 = ',ISEC1
  WRITE (UNIT=FA%NULOUT,FMT=*)                                   &
&                     '       * ILONSEC2 ! ISEC2(1:ILONSEC2) = ', &
&                     ILONSEC2, ' ! ', ISEC2(1:ILONSEC2)
  WRITE (UNIT=FA%NULOUT,FMT=*) '       * ZSEC2(1:2) = ',ZSEC2(1:2)
  IF (ISEC2(12).GT.0) WRITE (UNIT=FA%NULOUT,FMT=*)           &
&          '       * ISEC2(12) ! ZSEC2(11:10+ISEC2(12)) = ',  &
&                    ISEC2(12), ' ! ', ZSEC2(11:10+ISEC2(12))
  WRITE (UNIT=FA%NULOUT,FMT=*) '       * FA%JPSEC4 ! ISEC4 = ', &
&                               FA%JPSEC4,' ! ',ISEC4
ENDIF
!*
!     3.1 -  CONTROLES DE COHERENCE
!-----------------------------------------------------------------------
!
IF (IRET.GT.0) THEN
! Erreur rapportee par GRIBEX
  KREP=-1000-IRET
  WRITE (UNIT=FA%NULOUT,FMT=*) ' FADECX: IRET, KREP = ',IRET, KREP 
  GOTO 1001
ELSEIF (IRET.LT.0) THEN
! Warning rapporte par GRIBEX
  WRITE (UNIT=FA%NULOUT,FMT=*)
  WRITE (UNIT=FA%NULOUT,FMT=*)                         &
&         '!------------------------------------------'
  WRITE (UNIT=FA%NULOUT,FMT=*)                         &
&         '!           FADECX:   WARNING !!!         !'
  WRITE (UNIT=FA%NULOUT,FMT=*)                         &
&         '!------------------------------------------'
  WRITE (UNIT=FA%NULOUT,FMT=*) ' Code retour de GRIBEX = ', &
&        IRET,' pour le champ: ',CDNOMA
  WRITE (UNIT=FA%NULOUT,FMT=*)
ENDIF
IF (ISEC4(1).LT.IILCHAM) THEN
  KREP=-93
  IF (FA%LFAMOP) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                                &
&         'FADECX: ERREUR !!! Nbre de valeurs decodees = ',      &
&            ISEC4(1),' et nbre de valeurs attendues = ',IILCHAM
  ENDIF
  GOTO 1001
ELSEIF (ISEC4(1).GT.IILCHAM) THEN
  KREP=-94
  IF (FA%LFAMOP) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                             &
&         'FADECX: ERREUR !!! Nbre de valeurs decodees = ',   &
&         ISEC4(1),' et nbre de valeurs attendues = ',IILCHAM
  ENDIF
  IF (LLMOER(KREP,KRANG)) GOTO 1001
ENDIF
!
IF (IVALC3.NE.ISEC4(2).AND.FA%LFAMOP) THEN
  WRITE (UNIT=FA%NULOUT,FMT=*)                                   &
&     ' FADECX: WARNING, le nb de bits de codage qui avait',      &
&     ' ete demande ( ',IVALC3,' ) est different de celui qui a', &
&          ' ete finalement retenu ( ',ISEC4(2),' ) par GRIBEX.'
  WRITE (UNIT=FA%NULOUT,FMT=*)                       &
&         ' => Gain de place sans perte de precision'
ENDIF
!
!  Dans le cas d'un champ spectral ARPEGE
!
IF (LDCOSP.AND..NOT.LLMLAM.AND.(ISEC4(18).NE.IVALC4 &
&    .OR.ISEC4(17).NE.IVALC5)) THEN                  
  IF (FA%LFAMOP) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                             &
&            'Ss-tronc non compactee dans GRIB = ',ISEC4(18), &
&            ' et on attend: ',IVALC4
    WRITE (UNIT=FA%NULOUT,FMT=*)                             &
&            'Puissance de laplacien dans GRIB = ',ISEC4(17), &
&            ' et on attend: ',IVALC5
  ENDIF
  KREP=-95
  GOTO 1001
ENDIF
!
! Controle de l'adequation entre le nb de mots lus par LFI et le detail:
! ( enrobage FA + message GRIBEX + eventuelles valeurs non-compactees ).
!
IWORD=1+(ISEC0(1)-1)/JPDBLE
IF (FA%LFAMOP) THEN
  WRITE (UNIT=FA%NULOUT,FMT=*) ' FADECX: IWORD = ',IWORD
ENDIF
IPOFIN=IDECAL+IWORD
IF (LDCOSP) THEN
  IF (LLMLAM) THEN
    IPOFIN=IPOFIN+ISTRIA
  ELSE
    IPOFIN=IPOFIN+(1+IVALC4)*(2+IVALC4)
  ENDIF
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
!     3.2 -  DEMODULATION DES COEFF. SPEC. ALADIN QUI ONT ETE COMPACTES
!-----------------------------------------------------------------------
!
IF (LDCOSP.AND.LLMLAM) THEN
!  Transfert des donnees decodees et modulees entieres en nombres reels
!  pour les demoduler. Comme KCHAMP est a profil implicite, on ne peut
!  s'en servir pour la fonction TRANSFER => il faut passer par ICHAMP!
  ICHAMP(1:IILCHAM) = KCHAMP(1:IILCHAM)
  ZSEC4 (1:IILCHAM) = TRANSFER(ICHAMP,ZSEC4,IILCHAM)
  ZCHAMP=0._JPDBLR
  ZPULAP=REAL(IVALC5,JPDBLR) * (-0.001_JPDBLR)
  IIND=0
  DO JM=1,FA%CADRE(IRANGC)%NOMPAR(2)
    ILOW=2+2*JM+1
    IADD=4*MAX(IVALC4+1-JM,1_JPLIKB )
!
    DO IDX=FA%CADRE(IRANGC)%NOMPAR(ILOW)+IADD,FA%CADRE(IRANGC)%NOMPAR(ILOW+1)
      IIND=IIND+1        
      JN=(IDX-FA%CADRE(IRANGC)%NOMPAR(ILOW))/4
      ZCHAMP(IDX)=ZSEC4(IIND) *                 &
&           ((REAL(JN**2+JM**2,JPDBLR))**ZPULAP)
    ENDDO
  ENDDO
!  Transfert des donnees decodees et demodulees reelles en nombres entiers
!  disposes aux bons endroits du tableau definitif.
  ICHAMP(1:ILCHAM) = TRANSFER(ZCHAMP,ICHAMP,ILCHAM)
  KCHAMP(1:ILCHAM) = ICHAMP(1:ILCHAM)
ENDIF
!*
!     3.3 -  TRANSFERT DES COEFFICIENTS SPECTRAUX NON COMPACTES.
!-----------------------------------------------------------------------
!        (et non fournis par GRIBEX) stockes en fin d'article.
!
IF (LDCOSP) THEN
  IF (LLMLAM) THEN
    IIND=0
    DO JM=0,FA%CADRE(IRANGC)%NOMPAR(2)
      IL=2+2*JM+1
      ILOW=FA%CADRE(IRANGC)%NOMPAR(IL)
!
      IF (JM.EQ.0) THEN
        IHIGH=FA%CADRE(IRANGC)%NOMPAR(IL+1)
      ELSE
        IHIGH=ILOW+4*(IVALC4+1-JM)-1
        IF (IHIGH.LE.ILOW) IHIGH=ILOW+3
      ENDIF
!
      DO IDX=ILOW,IHIGH
        IIND=IIND+1
        KCHAMP(IDX)=KVALCO(IDECAL+IWORD+IIND)
      ENDDO
    ENDDO
  ELSE
!
! Cas ARPEGE
!
    KCHAMP(1:2*(IVALC4+1))=                              &
&        KVALCO(IDECAL+IWORD+1:IDECAL+IWORD+2*(IVALC4+1))
    IIND=2*(IVALC4+1)-1
    IDX=2*(ITRONC+1)-1
    DO JM=1,IVALC4
    DO JN=JM,ITRONC
      IDX=IDX+2
      IF (JN.LE.IVALC4) THEN
        IIND=IIND+2
        KCHAMP(IDX) = KVALCO(IDECAL+IWORD+IIND)
        KCHAMP(IDX+1) = KVALCO(IDECAL+IWORD+IIND+1)
      ENDIF
    ENDDO
    ENDDO
!
  ENDIF
ENDIF
!*
!     3.4 - Renversement des valeurs en pts de grille des champs
!            lat-lon afin de les ranger Sud-Nord plutot que Nord-Sud
!            (on conserve le rangt W-E consecutif) a l'image du rangt
!            initial effectue par FullPos.
!-----------------------------------------------------------------------
!
IF ((ISEC2(1)==0.OR.ISEC2(1)==10.OR.ISEC2(1)==20.OR. &
&    ISEC2(1)==30) .AND. .NOT.LDCOSP) THEN
  IF (FA%LFAMOP) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                             &
&            ' FADECX: Grille LAT-LON issue BDAP -> ',        &
&            ' renversement des valeurs pour etre rangees SN'
  ENDIF
  DO JLAT=1,FA%CADRE(IRANGC)%NLATIT
  DO JLON=1,FA%CADRE(IRANGC)%NXLOPA
    JN=JLON+FA%CADRE(IRANGC)%NXLOPA*(JLAT-1)
    IDX=JLON+FA%CADRE(IRANGC)%NXLOPA*(FA%CADRE(IRANGC)%NLATIT-JLAT)
    ICHAMP(IDX) = KCHAMP(JN)
  ENDDO
  ENDDO
  KCHAMP(1:ILCHAM) = ICHAMP(1:ILCHAM)
ENDIF
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
IF (ALLOCATED(ICHAMP)) DEALLOCATE ( ICHAMP, ZCHAMP, ZSEC4 )
LLFATA=LLMOER (KREP,KRANG)
!
IF (FA%LFAMOP.OR.LLFATA) THEN
  INIMES=2
  CLNSPR='FADECX'
  INUMER=JPNIIL
!
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I5,'', KRANG='',I4, &
&         '', CDNOMA='''''',A,'''''', KLONGA= '',I8,      &
&         '', LDCOSP='',L1)')                             &
&     KREP, KRANG, CDNOMA, KLONGA, LDCOSP
  CALL FAIPAR_MT_BB                                     &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CDNOMA,.FALSE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FADECX_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FADECX_BB                                    &
&           (KREP, KRANG, CDNOMA, KVALCO, KLONGA, KCHAMP, &
&           LDCOSP)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
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

CALL FADECX_MT_BB                                           &
&           (FA, KREP, KRANG, CDNOMA, KVALCO, KLONGA, KCHAMP, &
&           LDCOSP)

END SUBROUTINE

SUBROUTINE FADECX_MM                                    &
&           (KREP, KRANG, CDNOMA, KVALCO, KLONGA, KCHAMP, &
&           LDCOSP)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
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

CALL FADECX_MT_MM                                           &
&           (FA, KREP, KRANG, CDNOMA, KVALCO, KLONGA, KCHAMP, &
&           LDCOSP)

END SUBROUTINE

SUBROUTINE FADECX_MT_MM                                     &
&           (FA, KREP, KRANG, CDNOMA, KVALCO, KLONGA, KCHAMP, &
&           LDCOSP)
USE FA_MOD, ONLY : FA_COM,              &
&                   FA_COM_DEFAULT_INIT, &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
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

CALL FADECX_MT_BB                                           &
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
