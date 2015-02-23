! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FACODX_MT64                                            &
&                     (FA,  KREP,   KRANG,  CDPREF, KNIVAU, CDSUFF, &
&                    PSEC4, LDCOSP, KVALCO, KLONGD )
USE FA_MOD, ONLY : FA_COM, JPNIIL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
!      PREPARATION (codage GRIBEX) d'un CHAMP HORIZONTAL
!      destine a etre ecrit sur un fichier ARPEGE/ALADIN.
!       ( CODage d'un champ a l'aide de gribeX )
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KRANG  (Entree) ==> Rang de l'unite logique;
!                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
!                KNIVAU (Entree) ==> Niveau vertical eventuel;
!                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
!    ( Tableau ) PSEC4  (Entree) ==> Valeurs REELLES du champ a ecrire;
!                LDCOSP (Entree) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux;
!    ( Tableau ) KVALCO (Sortie) ==> Donnees destinees a l'ecriture;
!                KLONGD (Sortie) ==> Nombre de mots a ecrire;
!*
!       En mode multi-taches, il doit y avoir verrouillage du fichier
!     concerne avant l'appel au sous-programme.
!
!      Modifications
!      -------------
!         R. El Ouaraini : 03-Oct-06, introduire la nouvelle geometrie pour tester ERPK
!
!         JM AUDOIN  :  15 mai 2007 partie 5 changement unite 
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KRANG, KNIVAU, KLONGD
!
INTEGER (KIND=JPDBLE) KVALCO(*)
REAL (KIND=JPDBLR) PSEC4(*)
!
LOGICAL LDCOSP
!
CHARACTER CDPREF*(*), CDSUFF*(*)
!
#include "fagribexr.h"
!
REAL (KIND=JPDBLR), ALLOCATABLE :: ZSEC4(:)
INTEGER (KIND=JPDBLE), ALLOCATABLE :: IVALCO(:)
REAL (KIND=JPDBLR) :: ZMIN, ZA
REAL (KIND=JPDBLR) :: ZSEC2(10+2*(FA%JPXNIV+1)), ZSEC3(2), ZPULAP
!
INTEGER (KIND=JPLIKB) ISEC0(2), ISEC1(FA%JPSEC1)
INTEGER (KIND=JPLIKB) ISEC2(FA%JPSEC2), ISEC3(2)
INTEGER (KIND=JPLIKB) ISEC4(FA%JPSEC4), ILONSEC2
INTEGER (KIND=JPLIKB) ILENG, IWORD, IRET, JM, IPULAP
INTEGER (KIND=JPLIKB) ILCHAM, JN, IDECAL, ICPACK
INTEGER (KIND=JPLIKB) ITRONC, ILOW, IHIGH, IDIMNC, INBITS
INTEGER (KIND=JPLIKB) IL, IADD, IRANGC, IILCHAM, INIMES
INTEGER (KIND=JPLIKB) INUMER,  IDX, JLAT, JLON, IDECOPT
INTEGER (KIND=JPLIKB) IFAORI, IFAMOD, INBIMO
!
LOGICAL LLMLAM
!
CHARACTER(LEN=1) CLOPER
!
INTEGER (KIND=JPLIKB) DECF10
EXTERNAL DECF10
INTRINSIC LEN_TRIM
!
CHARACTER(LEN=FA%JPXNOM) CLACTI 
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FACODX_MT',0,ZHOOK_HANDLE)
CLACTI=''
IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
  KREP=-66
  GOTO 1001
ENDIF
ICPACK=FA%FICHIER(KRANG)%NSTROF
IRANGC=FA%FICHIER(KRANG)%NUCADR
LLMLAM=FA%CADRE(IRANGC)%LIMLAM
ITRONC=FA%CADRE(IRANGC)%MTRONC
!
IF (LLMLAM) THEN
  IF (LDCOSP) THEN
    ILONSEC2=21+FA%CADRE(IRANGC)%NOMPAR(2)
  ELSE
    ILONSEC2=22
  ENDIF
ELSE
  IF (LDCOSP) THEN
    ILONSEC2=22
  ELSE
    ILONSEC2=22+FA%CADRE(IRANGC)%NLATIT
  ENDIF
ENDIF
!
KVALCO(1)=FA%FICHIER(KRANG)%NFGRIB
IDECAL=3
IF (LDCOSP) THEN
  IF (LLMLAM) THEN
    ILCHAM=FA%CADRE(IRANGC)%NSFLAM
  ELSE
    ILCHAM=(1+ITRONC)*(2+ITRONC)
  ENDIF
  KVALCO(2)=1
  INBITS=FA%FICHIER(KRANG)%NBFCSP 
  IDECAL=IDECAL+2
ELSE
  ILCHAM=FA%CADRE(IRANGC)%NVAPDG
  KVALCO(2)=0
  INBITS=FA%FICHIER(KRANG)%NBFPDG
ENDIF
KVALCO(3)=INBITS
IILCHAM = ILCHAM
IDECOPT = 0
!**
!     2.  -  PREPARATION DU TABLEAU DE DONNEES A ECRIRE SUR LE FICHIER.
!-----------------------------------------------------------------------
!
ALLOCATE (ZSEC4 (ILCHAM))
!
IF (LDCOSP .AND. LLMLAM) THEN
!
!       Champ ALADIN en coefficients spectraux... traitement particulier,
!     car non prevu dans GRIBEX (il y sera considere comme un champ lat-lon)
!     mais on a la possibilite de compacter une (pseudo-)puissance de
!     laplacien du champ a la place du champ, de maniere a augmenter
!     la precision du champ en "aplanissant" le spectre.
!
!     Determination de la puissance de Laplacien (en 1/1000 ieme)
!
  CALL FAPULA_MT64                                 &
&                 (FA,  KREP, KRANG, PSEC4, IPULAP )
  ZPULAP=REAL(IPULAP,JPDBLR)/1000._JPDBLR
!       ZPULAP=0.
!       IPULAP=0
  IF (FA%LFAMOP) THEN
    print *,'FACODX: puissance de laplacien selectionee ',ZPULAP, &
&          ' pour une sous-tronc de ',ICPACK
  ENDIF
  IF (KREP.NE.0) GOTO 1001
!
! Transfert des coeff spectraux devant etre compactes de PSEC4 a ZSEC4
! avec prise en compte du coefficient (n**2+m**2)**zpulap. Les coefficients
! concernes sont ceux inclus dans le quart de l'ellipse, hors axes (coeff
! nuls), et hors du triangle non-compacte (sous-troncature).
  IILCHAM=0
!
  DO JM=1,FA%CADRE(IRANGC)%NOMPAR(2)
    ILOW=2+2*JM+1
    IADD=4* MAX(ICPACK+1-JM,1_JPLIKB )
!
    DO IDX=FA%CADRE(IRANGC)%NOMPAR(ILOW)+IADD,FA%CADRE(IRANGC)%NOMPAR(ILOW+1)
      IILCHAM=IILCHAM+1
      JN=(IDX-FA%CADRE(IRANGC)%NOMPAR(ILOW))/4
      ZSEC4(IILCHAM)=PSEC4(IDX) *              &
&         ((REAL(JN**2+JM**2, JPDBLR))**ZPULAP)
    ENDDO
  ENDDO
! Number of elements in sub-triangle+axes:IDIMNC
  IDIMNC=ILCHAM-IILCHAM
! Recherche de l'amplitude et du min du champ
  ZMIN=ZSEC4(1)
  ZA=0._JPDBLR
  ZMIN = MINVAL(ZSEC4(1:IILCHAM))
  ZA   = MAXVAL(ZSEC4(1:IILCHAM)) - ZMIN
! Recherche du facteur decimal optimal pour utiliser
! au mieux les INBITS dans le codage de ce champ
  IF (FA%LFAMOP) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)'FACODX: traitement du champ: ', &
&          CDPREF,KNIVAU,CDSUFF
  ENDIF
  CALL FACDEC_MT64                                    &
&                 (FA, KREP, ZA, ZMIN, INBITS, IDECOPT)
  IF (KREP.NE.0) THEN
    IDECOPT = 0
    KREP = 0
  ENDIF
ELSEIF(LDCOSP .AND. .NOT.LLMLAM) THEN
!
!          Transfert du tableau d'entree dans un tableau local
!     de maniere a eviter l'ecrasement du tableau d'entree par "GRIBEX".
!
  ZSEC4(1:IILCHAM) = PSEC4(1:IILCHAM)
  IDIMNC=(1+ICPACK)*(2+ICPACK)
ELSE
!
!    CHAMPS NON SPECTRAUX: transfert du tableau d'entree dans un
!    tableau local de maniere a eviter son ecrasement par "GRIBEX".
!
!
  IDIMNC=0
! Tester si Nouvelle ou ancienne geometrie Aladin
IF (FA%CADRE(IRANGC)%SINLAT(1) .GE. 0) THEN
  IF (LLMLAM .AND. FA%CADRE(IRANGC)%SINLAT(10).LT.0) THEN
!  Parametre de projection negatif, donc pas de projection:
!  Il s'agit d'une grille lat-lon reguliere du type Full-Pos
!  (pour champ ARPEGE ou Aladin). Il faut donc renverser
!  le champ afin de ranger Nord-Sud les valeurs plutot que Sud-Nord
!  (on conserve le rangt W-E consecutif).
!  Le but est de satisfaire la BDAP qui attend un rangt NW-->SE.
!
    IF (FA%LFAMOP) THEN
      WRITE (UNIT=FA%NULOUT,FMT=*)                             &
&              ' FACODX: Grille LAT-LON pour BDAP -> ',         &
&              ' renversement des valeurs pour etre rangees NS'
    ENDIF
    DO JLAT=1,FA%CADRE(IRANGC)%NLATIT
    DO JLON=1,FA%CADRE(IRANGC)%NXLOPA
      JN=JLON+FA%CADRE(IRANGC)%NXLOPA*(JLAT-1)
      IDX=JLON+FA%CADRE(IRANGC)%NXLOPA*(FA%CADRE(IRANGC)%NLATIT-JLAT)
      ZSEC4(IDX) = PSEC4(JN)
    ENDDO
    ENDDO
  ELSE
    ZSEC4(1:IILCHAM) = PSEC4(1:IILCHAM)
  ENDIF
ELSE
  IF (LLMLAM .AND. FA%CADRE(IRANGC)%SINLAT(2).LT.0) THEN
    IF (FA%LFAMOP) THEN
      WRITE (UNIT=FA%NULOUT,FMT=*)                             &
&              ' FACODX: Grille LAT-LON pour BDAP -> ',         &
&              ' renversement des valeurs pour etre rangees NS'
    ENDIF
    DO JLAT=1,FA%CADRE(IRANGC)%NLATIT
    DO JLON=1,FA%CADRE(IRANGC)%NXLOPA
      JN=JLON+FA%CADRE(IRANGC)%NXLOPA*(JLAT-1)
      IDX=JLON+FA%CADRE(IRANGC)%NXLOPA*(FA%CADRE(IRANGC)%NLATIT-JLAT)
      ZSEC4(IDX) = PSEC4(JN)
    ENDDO
    ENDDO
  ELSE
    ZSEC4(1:IILCHAM) = PSEC4(1:IILCHAM)
  ENDIF
ENDIF
! Recherche de l'amplitude et du min du champ
  ZMIN=ZSEC4(1)
  ZA=0._JPDBLR
  ZMIN = MINVAL(ZSEC4(1:IILCHAM))
  ZA   = MAXVAL(ZSEC4(1:IILCHAM)) - ZMIN
! Recherche du facteur decimal optimal pour utiliser
! au mieux les INBITS dans le codage de ce champ
  IF (FA%LFAMOP) THEN 
  WRITE (UNIT=FA%NULOUT,FMT=*)'FACODX: traitement du champ: ', &
&          CDPREF,KNIVAU,CDSUFF
  ENDIF
  CALL FACDEC_MT64                                    &
&                 (FA, KREP, ZA, ZMIN, INBITS, IDECOPT)
  IF (KREP.NE.0) THEN
    IDECOPT = 0
    KREP = 0
  ENDIF
ENDIF
!*
!     3.  -  INITIALISATION DE L'ENROBAGE GRIB
!-----------------------------------------------------------------------
!
!     3.1 -  Sections 1, 2, 3 et 4 (sf la partie reelle pour 4)
!
CALL FAINIG_MT64                                                 &
&               (FA,  KREP, KRANG, CDPREF, KNIVAU, CDSUFF, LDCOSP, &
&              IILCHAM, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4)
IF (KREP.NE.0) THEN
  GOTO 1001
ENDIF
! Prise en compte du facteur decimal
ISEC1(23) = IDECOPT
!
!     3.2 -  Definition du type de codage
!
CLOPER='C'
IF (FA%FICHIER(KRANG)%NCOGRIF(1)==1) CLOPER='K'
!*
!     4.  -  CODAGE GRIB PROPREMENT DIT
!-----------------------------------------------------------------------
!
IRET=-1
! ILENG=longueur disponible en nb d'"entiers declares INTEGER" dans KVALCO.
! On part de l'hypothese ou le dimensionnement de KVALCO se fait
! dans la routine appelante a ILCHAM+2 (cas de l'absence de compactage).
ILENG=(KIND(KVALCO)/4)*(ILCHAM+2-IDECAL)
IWORD=0
!DP
!DP  TEST AVEC UNE PUISSANCE DE LAPLACIEN IMPOSEE
!DP
!DP   CALL GRSMKP(0)
!DP ISEC4(17) = 2000
!DP
IF (FA%LFAMOP) THEN
  WRITE (UNIT=FA%NULOUT,FMT=*)' FACODX: CLOPER = ',CLOPER
  WRITE (UNIT=FA%NULOUT,FMT=*)                                    &
&                    ' FACODX: IILCHAM, ILCHAM, IDECAL, ILENG = ', &
&                    IILCHAM, ILCHAM, IDECAL, ILENG
  WRITE (UNIT=FA%NULOUT,FMT=*)'       * ISEC1 = ',ISEC1
  WRITE (UNIT=FA%NULOUT,FMT=*)                                  &
&                    '       * ILONSEC2 ! ISEC2(1:ILONSEC2) = ', &
&                    ILONSEC2,' ! ', ISEC2(1:ILONSEC2)
  WRITE (UNIT=FA%NULOUT,FMT=*) '       * ZSEC2(1:2) = ',ZSEC2(1:2)
  IF (ISEC2(12).GT.0) WRITE (UNIT=FA%NULOUT,FMT=*)           &
&          '       * ISEC2(12) ! ZSEC2(11:10+ISEC2(12)) = ',  &
&                    ISEC2(12), ' ! ', ZSEC2(11:10+ISEC2(12))
  WRITE (UNIT=FA%NULOUT,FMT=*)'       * FA%JPSEC4 ! ISEC4 = ', &
&                              FA%JPSEC4,' ! ',ISEC4
  WRITE (UNIT=FA%NULOUT,FMT=*)'       * ZSEC4(1:20) = ', &
&                              ZSEC4(1:20)
ENDIF
!     WARNING GRIBEX ENLEVE 
CALL GRSDBG(0)
CALL GRSVCK(0)

CALL FAGRIBEX(ISEC0,ISEC1,ISEC2,ZSEC2,ISEC3,ZSEC3,ISEC4,  &
&              ZSEC4,IILCHAM,KVALCO(IDECAL+1),ILENG,IWORD, &
&              CLOPER,IRET)
!
IF (IRET.GT.0) THEN
! Erreur rapportee par GRIBEX
  KREP=-1000-IRET
  GOTO 1001
ELSEIF (IRET.LT.0) THEN
! Warning rapporte par GRIBEX
  WRITE (UNIT=FA%NULOUT,FMT=*)
  WRITE (UNIT=FA%NULOUT,FMT=*)                               &
&               '!------------------------------------------'
  WRITE (UNIT=FA%NULOUT,FMT=*)                               &
&               '!           FACODX:   WARNING !!!         !'
  WRITE (UNIT=FA%NULOUT,FMT=*)                               &
&               '!------------------------------------------'
  WRITE (UNIT=FA%NULOUT,FMT=*) ' Code retour de GRIBEX = ', &
&        IRET,' pour le champ: ',CDPREF,KNIVAU,CDSUFF
  WRITE (UNIT=FA%NULOUT,FMT=*)
ENDIF
!
! ISEC0(1) = nb d'octets dans le message GRIB
! IWORD    = nb de mots de JBDBLE octets (64 bits) du message GRIB
IWORD=1+(ISEC0(1)-1)/JPDBLE
KLONGD=IDECAL+IWORD+IDIMNC
IF (FA%LFAMOP) THEN
  WRITE (UNIT=FA%NULOUT,FMT=*)                                   &
&         ' FACODX: longueur du GRIB en nb octets et en mots = ', &
&         ISEC0(1), IWORD
  WRITE (UNIT=FA%NULOUT,FMT=*)                            &
&         ' FACODX: longueur de l''article FA en mots = ', &
&         KLONGD
  IF (ISEC4(4).EQ.64 .AND. ISEC4(3).EQ.128) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                         &
&           ' FACODX: complex packing with P=',ISEC4(17), &
&           ' and sub trunc = ',ISEC4(18)
  ENDIF
ENDIF
!
!  CAS D'UN DEPASSEMENT DE LA TAILLE MAX DE L'ARTICLE FINAL
!  On ramene ce cas a celui d'un tableau trop petit dans GRIBEX.
!
IF (KLONGD.GT.ILCHAM+2) THEN
  IF (FA%LFAMOP) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                          &
&            ' FACODX: article FA + long avec compactage', &
&            ' que sans => on le supprime'
  ENDIF
  IRET=710
  KREP=-1000-IRET
  GOTO 1001
ENDIF
!
!
!*
!     5.  -  CHANGEMENT D'UNITE DE CERTAINS CHAMPS.
!            Il s'agit de champs dont les valeurs sont comprises
!            entre 0 et 1 dans le modele mais dont l'unite
!            conventionnelle dans le GRIB est le %.
!---------------------------------------------------------------
!
IDX = 0
IF (CDPREF=='SURF') THEN
  IF (CDSUFF(1:6)=='NEBUL.'  .OR. CDSUFF(1:6)=='ALBEDO'  .OR. &
&      CDSUFF=='PROP.VEGETAT' ) THEN
    IDX = 1
  ENDIF
ELSEIF (CDPREF(1:LEN(TRIM(CDPREF)))=='CLS') THEN
  IF (CDSUFF=='HUMI.RELATIVE' .OR. CDSUFF=='MAXI.HUMI.REL' .OR. &
&      CDSUFF=='MINI.HUMI.REL' ) THEN
    IDX = 1
  ENDIF
ELSEIF (CDPREF(1:1)=='P'.OR.CDPREF(1:1)=='H') THEN
!         blindage 
  IF (CDSUFF(1:10)=='HUMI_RELAT') THEN
    IDX = 1
  ENDIF
ENDIF
IF (IDX==1) THEN
  IADD = -2
  INBIMO = 32  ! Nombre de BIts par mot (un mot=INTEGER)
  IDX = DECF10 ( KVALCO(IDECAL+1), ILENG, IADD, &
&                   IFAORI, IFAMOD, INBIMO )
  IF (IDX==-1) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                   &
&           'FACODX: pas d''entete GRIB au debut !'
    KREP=-128
    GOTO 1001
  ELSEIF (IDX==-2) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                           &
&           'FACODX: edition du GRIB invalid pour DECF10 !'
    KREP=-128
    GOTO 1001
  ELSEIF (IDX > 0) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                             &
&           'FACODX: ERREUR dans appel a INXBIT par DECF10 !'
    WRITE (UNIT=FA%NULOUT,FMT=*)                       &
&           '        avec code retour de INXBIT = ',IDX
    KREP=-128
    GOTO 1001
  ELSEIF (IDX < -2) THEN
    WRITE (UNIT=FA%NULOUT,FMT=*)                          &
&           'FACODX: code retour inconnu de DECF10 : ',IDX
    KREP=-126
    GOTO 1001
  ENDIF
ENDIF
!
!*
!     6.  -  TRANSFERT DES COEFFICIENTS SPECTRAUX NON COMPACTES.
!-----------------------------------------------------------------------
!        (et non traites par GRIBEX) en fin d'article.
!
IF (LDCOSP) THEN
  KVALCO(4)=ICPACK
  IF (LLMLAM) THEN
    KVALCO(5)=IPULAP
! Copy nonpacked part of PSEC4 (sub-triangle+axes) into KVALCO
    IILCHAM=0
    DO JM=0,FA%CADRE(IRANGC)%NOMPAR(2)
      IL=2+2*JM+1
      ILOW=FA%CADRE(IRANGC)%NOMPAR(IL)
!
      IF (JM.EQ.0) THEN
        IHIGH=FA%CADRE(IRANGC)%NOMPAR(IL+1)
      ELSE
        IHIGH=ILOW+4*(ICPACK+1-JM)-1
        IF (IHIGH.LE.ILOW) IHIGH=ILOW+3
      ENDIF     
!
      DO IDX=ILOW,IHIGH
        IILCHAM=IILCHAM+1
        ZSEC4(IILCHAM)=PSEC4(IDX)
      ENDDO
    ENDDO
    IF (IILCHAM.NE.IDIMNC) THEN
      WRITE (UNIT=FA%NULOUT,FMT='(A35,I10,A11,I10)')       &
&            'FACODX: incoherence entre IILCHAM= ',IILCHAM, &
&            'et IDIMNC= ',IDIMNC
      KREP=-126
      GOTO 1001
    ENDIF
  ELSE
    KVALCO(5)=ISEC4(17)
! Recuperation des coeff spectraux non compactes sachant que le
! rangement est fait par colonnes de JM=cst juxtaposees
    ZSEC4(1:2*(ICPACK+1))=PSEC4(1:2*(ICPACK+1))
    IILCHAM=2*(ICPACK+1)-1
    IDX=2*(ITRONC+1)-1
    DO JM=1,ICPACK
    DO JN=JM,ITRONC
      IDX=IDX+2
      IF (JN.LE.ICPACK) THEN
        IILCHAM=IILCHAM+2
        ZSEC4(IILCHAM) = PSEC4(IDX)
        ZSEC4(IILCHAM+1) = PSEC4(IDX+1)
      ENDIF
    ENDDO
    ENDDO
    IF (IILCHAM+1.NE.IDIMNC) THEN
      WRITE (UNIT=FA%NULOUT,FMT='(A35,I10,A11,I10)')           &
&            'FACODX: incoherence entre IILCHAM+1= ',IILCHAM+1, &
&            'et IDIMNC= ',IDIMNC
      KREP=-126
      GOTO 1001
    ENDIF
  ENDIF
! Les IDIMNC coeff spectraux non compactes doivent etre transferes
! sur le tableau d'entiers KVALCO apres le IDECAL+IWORD ieme elt.
!
!       KVALCO(IDECAL+IWORD+1:KLONGD)=TRANSFER(ZSEC4,KVALCO,IDIMNC)
  ALLOCATE (IVALCO(IDIMNC))
  IVALCO(1:IDIMNC)=TRANSFER(ZSEC4,IVALCO,IDIMNC)
  KVALCO(IDECAL+IWORD+1:KLONGD)=IVALCO(1:IDIMNC)
  DEALLOCATE (IVALCO)
ENDIF
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
IF (ALLOCATED(ZSEC4)) DEALLOCATE ( ZSEC4 )
!
! Cas particulier de l'erreur GRIBEX num 710: OUTPUT ARRAY TOO SMALL
! On s'en sert pour detecter un probleme de compactage lie a ce que
! le champ compacte+les descripteurs prennent plus de place que le
! champ non compacte...
! On sort donc du compactage (FACODX) pour demander un codage sans
! compactage (FACINE) avec rangement des valeurs selon le modele:
! FA%NFGRIB=-1.
!
IF (IRET==710) THEN
  CLNSPR='FACODX'
  INIMES=2
  INUMER=JPNIIL
!
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I5,'', KRANG='',I4, &
&         '', CDPREF='''''',A,'''''', KNIVAU='',I6,       &
&         '', CDSUFF='''''',A,'''''', LDCOSP= '',L1,      &
&         '', KLONGD='',I6)')                             &
&     KREP, KRANG, CDPREF(1:LEN_TRIM(CDPREF)), KNIVAU,    &
&     CDSUFF(1:LEN_TRIM(CDSUFF)), LDCOSP, KLONGD
  CALL FAIPAR_MT64                                      &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CLACTI,.FALSE.)
  CLMESS=                                                          &
& ' CAUTION: this field is not packed or it will occupy more space'
  CALL FAIPAR_MT64                                      &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CLACTI,.FALSE.)
  IF (LHOOK) CALL DR_HOOK('FACODX_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
!
!
LLFATA=LLMOER (KREP,KRANG)
!
IF (FA%LFAMOP.OR.LLFATA) THEN
  INIMES=2
  CLNSPR='FACODX'
  INUMER=JPNIIL
!
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I5,'', KRANG='',I4, &
&         '', CDPREF='''''',A,'''''', KNIVAU='',I6,       &
&         '', CDSUFF='''''',A,'''''', LDCOSP= '',L1,      &
&         '', KLONGD='',I6)')                             &
&     KREP, KRANG, CDPREF(1:LEN_TRIM(CDPREF)), KNIVAU,    &
&     CDSUFF(1:LEN_TRIM(CDSUFF)), LDCOSP, KLONGD
  CALL FAIPAR_MT64                                      &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CLACTI,.FALSE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FACODX_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FACODX64                                            &
&           (KREP, KRANG, CDPREF, KNIVAU, CDSUFF, PSEC4, LDCOSP, &
&           KVALCO, KLONGD)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKB)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
REAL (KIND=JPDBLR)     PSEC4      (*)                         ! IN   
LOGICAL                LDCOSP                                 ! IN   
INTEGER (KIND=JPDBLE)  KVALCO     (*)                         !   OUT
INTEGER (KIND=JPLIKB)  KLONGD                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FACODX_MT64                                           &
&           (FA, KREP, KRANG, CDPREF, KNIVAU, CDSUFF, PSEC4, &
&           LDCOSP, KVALCO, KLONGD)

END SUBROUTINE

SUBROUTINE FACODX                                              &
&           (KREP, KRANG, CDPREF, KNIVAU, CDSUFF, PSEC4, LDCOSP, &
&           KVALCO, KLONGD)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
REAL (KIND=JPDBLR)     PSEC4      (*)                         ! IN   
LOGICAL                LDCOSP                                 ! IN   
INTEGER (KIND=JPDBLE)  KVALCO     (*)                         !   OUT
INTEGER (KIND=JPLIKM)  KLONGD                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FACODX_MT                                             &
&           (FA, KREP, KRANG, CDPREF, KNIVAU, CDSUFF, PSEC4, &
&           LDCOSP, KVALCO, KLONGD)

END SUBROUTINE

SUBROUTINE FACODX_MT                                       &
&           (FA, KREP, KRANG, CDPREF, KNIVAU, CDSUFF, PSEC4, &
&           LDCOSP, KVALCO, KLONGD)
USE FA_MOD, ONLY : FA_COM
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
REAL (KIND=JPDBLR)     PSEC4      (*)                         ! IN   
LOGICAL                LDCOSP                                 ! IN   
INTEGER (KIND=JPDBLE)  KVALCO     (*)                         !   OUT
INTEGER (KIND=JPLIKM)  KLONGD                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  IRANG                                  ! IN   
INTEGER (KIND=JPLIKB)  INIVAU                                 ! IN   
INTEGER (KIND=JPLIKB)  ILONGD                                 !   OUT
! Convert arguments

IRANG      = INT (     KRANG, JPLIKB)
INIVAU     = INT (    KNIVAU, JPLIKB)

CALL FACODX_MT64                                           &
&           (FA, IREP, IRANG, CDPREF, INIVAU, CDSUFF, PSEC4, &
&           LDCOSP, KVALCO, ILONGD)

KREP       = INT (      IREP, JPLIKM)
KLONGD     = INT (    ILONGD, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                                                              
!INTF KRANG         IN                                                                 
!INTF CDPREF        IN                                                                 
!INTF KNIVAU        IN                                                                 
!INTF CDSUFF        IN                                                                 
!INTF PSEC4         IN    DIMS=*                                                       
!INTF LDCOSP        IN                                                                 
!INTF KVALCO          OUT DIMS=*                         KIND=JPDBLE                   
!INTF KLONGD          OUT                                                              
