! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
! Sep-2012 P. Marguinaud Fix uninitialized variables
SUBROUTINE FACINE_MT64                                            &
&                     (FA,  KREP,   KRANG,  CDNOMA, KCHAMP, LDCOSP, &
&                      KVALCO, KLONGD, KB1PAR )
USE FA_MOD, ONLY : FA_COM, JPNIIL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
!      PREPARATION (codage GRIB ou non) d'un CHAMP HORIZONTAL
!      destine a etre ecrit sur un fichier ARPEGE/ALADIN.
!       ( Codage Interne d'un (Nouveau ?) champ a Ecrire )
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KRANG  (Entree) ==> Rang de l'unite logique;
!                CDNOMA (Entree) ==> Nom d'article (prefabrique);
!    ( Tableau ) KCHAMP (Entree) ==> Valeurs REELLES du champ a ecrire;
!                LDCOSP (Entree) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux;
!    ( Tableau ) KVALCO (Sortie) ==> Donnees destinees a l'ecriture;
!                KLONGD (Sortie) ==> Nombre de mots a ecrire;
!    ( Tableau ) KB1PAR (Entree+ ==> Image des parametres de la section
!                        Sortie)     1 de GRIB.
!*
!       En mode multi-taches, il doit y avoir verrouillage du fichier
!     concerne avant l'appel au sous-programme.
!
!     Modifications
!     -------------
!  Juin 2001, D. Paradis, DT/DSI/DEV:
!
!    -retrait du compactage lorsqu'il conduit a un article de longueur
!     superieure a celle obtenue sans le compactage (permet aussi de
!     dimensionner KVALCO a ILCHAM+2 sans risquer un debordement)
!
!  Avril 2004, D. Paradis, DT/DSI/DEV:
!
!    -declaration de FA%ICHAMP et FA%ICHAUX en ALLOCATABLE (gain memoire)
!
! January 2010 Trygve Aspelien & Ryad El Khatib : 
!    - workaround against memory leak on IBM
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KRANG, KLONGD
!
INTEGER (KIND=JPDBLE) KCHAMP(*), KVALCO(*)
INTEGER (KIND=JPLIKB) KB1PAR (FA%JPLB1P)
!
LOGICAL LDCOSP
!
CHARACTER CDNOMA*(*)
!
INTEGER (KIND=JPLIKB) ILCHAM, ISTRIA, IVALC1, IVALC2
INTEGER (KIND=JPLIKB) J, IDECAL, ICPACK, IPUILA
INTEGER (KIND=JPLIKB) ITRONC, IIND, ILOW, IHIGH, JTRON
INTEGER (KIND=JPLIKB) IDIMNC, ILDISP, INBITS
INTEGER (KIND=JPLIKB) IL, IADD, IRANGC, IARR, IILCHAM
INTEGER (KIND=JPLIKB) INMOCC, IERR, INIMES
INTEGER (KIND=JPLIKB) INUMER, ITRONC2, ILONGFA
INTEGER (KIND=JPLIKB) ILONGSEC, ILONGDATA
INTEGER (KIND=JPLIKB) ILONGD
!
INTEGER (KIND=JPLIKB) IB2PAR (FA%JPLB2P)
!
LOGICAL LLARPE, LLMLAM
!
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FACINE_MT',0,ZHOOK_HANDLE)
ICPACK=0
IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
  KREP=-66
  GOTO 1001
ENDIF
!**
!     2.  -  FABRICATION DE L'ARTICLE A ECRIRE SUR LE FICHIER.
!-----------------------------------------------------------------------
!
IVALC1=FA%FICHIER(KRANG)%NFGRIB
LLARPE=IVALC1.EQ.2
IRANGC=FA%FICHIER(KRANG)%NUCADR
LLMLAM=FA%CADRE(IRANGC)%LIMLAM
!
!  Initialisation du nombre de valeurs du champ (ILCHAM)
!                 du type de champ (IVALC2): spectral/pdg
!                 du type de representation de donnees (IB2PAR(1))
!                 du nombre de bits utilises pour le compactage
!                 de la longueur (en mots de 64 bits) de l'enrobage de l'article
!                          FA (ILONGFA) + les donnees non compactees (en spectral)
!                 de la longueur (en octets) des sections 0, 1 et 5 du GRIB
!                          devant apparaitre dans l'article FA (ILONGSEC)
!                 de la longueur (en bits) de la section 4 du GRIB (ILONGDATA),
!                          devant etre un multiple de 16 bits.
!
! ILONGSEC= 4 (S0) + 24 (S1) + 4 (S5) pour le GRIB version 0
ILONGSEC=32
IF (LDCOSP) THEN
!
  IVALC2=1
  INBITS=FA%FICHIER(KRANG)%NBFCSP
!
  IF (LLMLAM) THEN
    IB2PAR(1)=34
    ILCHAM=FA%CADRE(IRANGC)%NSFLAM
    IF (IVALC1.GT.0) THEN
! calcul du nombre de coefficients non compactes ISTRIA
      ICPACK=FA%FICHIER(KRANG)%NSTROF
      ITRONC=FA%CADRE(IRANGC)%MTRONC
      ITRONC2=-FA%CADRE(IRANGC)%NTYPTR
      ISTRIA=4*(1+ITRONC+ITRONC2+(ICPACK*(ICPACK-1)/2))
!
      ILONGFA=3+2*IVALC1+ISTRIA
!  Les 88 bits correspondent aux 11 octets d'enrobage GRIB V0 de la section 4
      ILONGDATA=(ILCHAM-ISTRIA)*INBITS + 88
    ENDIF
  ELSE
    IF (IVALC1.EQ.-1) THEN
      ILCHAM=(1+FA%CADRE(IRANGC)%MTRONC)*(2+FA%CADRE(IRANGC)%MTRONC)
    ELSE
      ILCHAM=(1+FA%CADRE(IRANGC)%MTRONC)**2
    ENDIF
    IB2PAR(1)=80
    IF (IVALC1.GT.0) THEN
      ICPACK=FA%FICHIER(KRANG)%NSTROF
! calcul du nombre de coefficients non compactes IDIMNC
      IDIMNC=(1+ICPACK)**2
      ILONGFA=3+2*IVALC1+IDIMNC
!  Les 144 bits correspondent aux 18 octets d'enrobage GRIB V0 de la section 4
      ILONGDATA=ILCHAM*INBITS + IDIMNC*(32-INBITS) + 144
    ENDIF
  ENDIF
!
ELSE
!
  ILCHAM=FA%CADRE(IRANGC)%NVAPDG
  IVALC2=0
  IB2PAR(1)=34
  INBITS=FA%FICHIER(KRANG)%NBFPDG
  IF (IVALC1.GT.0) THEN
    ILONGFA=1+2*IVALC1
!  Les 88 bits correspondent aux 11 octets d'enrobage GRIB V0 de la section 4
    ILONGDATA=ILCHAM*INBITS + 88
  ENDIF
!
ENDIF
!
!  Retrait du compactage si celui-ci s'avere conduire
!  a un article plus long qu'en l'absence de compactage:
!
IF (IVALC1.GT.0) THEN
! Arrondi de ILONGDATA au premier multiple de 16 superieur ou egal
  ILONGDATA=16*(1+(ILONGDATA-1)/16)
! Calcul du nombre de mots (64 bits) de la partie GRIB
  ILONGD=1+(ILONGDATA+8*ILONGSEC-1)/64
! On ajoute l'enrobage FA et les eventuelles donnees non compactees
  ILONGD=ILONGD+ILONGFA
  IF (ILONGD.GT.ILCHAM+2) THEN
!          WRITE (FA%NULOUT,*)
!     S    '///// FACINE:  the packing of article ',CDNOMA,
!     S    ' is not done because it will generate'
!          WRITE (FA%NULOUT,*)
!     S    '         a size ( ',ILONGD,' words ) bigger than',
!     S    ' the one ( ',ILCHAM+2,' words ) obtained without packing.'
    IVALC1=0
!          WRITE (FA%NULOUT,*)
  ENDIF
ENDIF
!
ISTRIA = 0
!
IF (IVALC1.EQ.-1.OR.IVALC1.EQ.0) THEN
!
!          Cas ou il n'y a aucun codage...
!     transfert du tableau d'entree a la suite des 2 valeurs
!     documentaires stockees ci-dessus dans KVALCO.
!
  DO J=1,ILCHAM
  KVALCO(2+J)=KCHAMP(J)
  ENDDO
!
  KLONGD=2+ILCHAM
!
ELSE
!
!          Cas avec codage GRIB (standard ou non).
!
#ifndef RS6K
  IF ( ASSOCIATED(FA%ICHAMP) ) DEALLOCATE(FA%ICHAMP)
  ALLOCATE (FA%ICHAMP (ILCHAM))
#else
  IF (.NOT. ASSOCIATED(FA%ICHAMP)) THEN
    ALLOCATE (FA%ICHAMP (ILCHAM))
  ELSEIF ( ILCHAM .GT. SIZE(FA%ICHAMP) ) THEN
    DEALLOCATE(FA%ICHAMP)
    ALLOCATE (FA%ICHAMP (ILCHAM))
  ENDIF
#endif
!
  IDECAL=1+2*IVALC1
  KB1PAR(1)=98
  KB1PAR(2)=1
  KB1PAR(3)=254
  KB1PAR(4)=0
  KB1PAR(5)=255
  KB1PAR(9)=MOD (FA%FICHIER(KRANG)%MADATE(1),100_JPLIKB )
!
  DO J=2,FA%JPLDAT
  KB1PAR(8+J)=FA%FICHIER(KRANG)%MADATE(J)
  ENDDO
!
  IB2PAR(6)=2
  IPUILA=FA%FICHIER(KRANG)%NPUFLA
  ITRONC=FA%CADRE(IRANGC)%MTRONC
!
  IF (LDCOSP) THEN
!
!           Champ en coefficients spectraux... traitement particulier,
!     lie a la possibilite de compacter une (pseudo-)puissance de
!     laplacien du champ a la place du champ, de maniere a augmenter
!     la precision du champ en "aplanissant" le spectre.
!
    CALL FACSIM_MT64                                              &
&                   (FA, KREP,KRANG,KCHAMP,FA%ICHAMP,IPUILA,ICPACK)
    IF (FA%LFAMOP) THEN
      print *,'FACINE: puissance Dolby selectionnee ',IPUILA
    ENDIF
    IF (KREP.NE.0) GOTO 1001
    IF (LLMLAM) THEN
! Copy the elements to be compacted from FA%ICHAMP to a work array
! This is that part of the quarter-ellipse which is out of the triangle of no compacting.
! In addition, the axes of ellipse are also excluded because of zero-coefficients
#ifndef RS6K
      ALLOCATE (FA%ICHAUX (ILCHAM))
#else
      IF (.NOT. ASSOCIATED(FA%ICHAUX)) THEN
        ALLOCATE (FA%ICHAUX (ILCHAM))
      ELSEIF ( ILCHAM .GT. SIZE(FA%ICHAUX) ) THEN
        DEALLOCATE(FA%ICHAUX)
        ALLOCATE (FA%ICHAUX (ILCHAM))
      ENDIF
#endif
      IIND=0
!
      DO JTRON=1,ITRONC
      ILOW=2+2*JTRON+1
      IADD=4* MAX(ICPACK+1-JTRON,1_JPLIKB )
!
      DO J=FA%CADRE(IRANGC)%NOZPAR(ILOW)+IADD, &
&                FA%CADRE(IRANGC)%NOZPAR(ILOW+1)
      IIND=IIND+1
      FA%ICHAUX(IIND)=FA%ICHAMP(J)
      ENDDO
      ENDDO
! Number of elements in sub-triangle+axes:ISTRIA
      ISTRIA=ILCHAM-IIND
      IDIMNC=0
    ELSE
      ISTRIA=IDIMNC
    ENDIF
!
    IDECAL=IDECAL+2
    ILDISP=ILCHAM+2-IDECAL-(IVALC1-1)*ISTRIA
!
    IF (.NOT.LLARPE) THEN
!
!            Recopie des coefficients spectraux "non compactes",
!     qui sont codes en fait sur 32 bits dans le cas standard de GRIB.
!
      DO J=1,IDIMNC
      FA%ICHAMP(J)=KCHAMP(J)
      ENDDO
!
    ENDIF
!
  ELSE
!
!          Transfert du tableau d'entree dans un tableau local
!     de maniere a eviter l'ecrasement du tableau d'entree par "CODEGA".
!
    DO J=1,ILCHAM
    FA%ICHAMP(J)=KCHAMP(J)
    ENDDO
!
    IDIMNC=0
    ILDISP=ILCHAM+2-IDECAL
  ENDIF
!*
!     3.1 -  CODAGE GRIB PROPREMENT DIT.
!-----------------------------------------------------------------------
!
  IARR=0
!
  IF (LDCOSP.AND.LLMLAM) THEN
    IILCHAM=ILCHAM-ISTRIA
    CALL FACODEGA(FA%ICHAUX,IILCHAM,INBITS,FA%NBIMAC,KB1PAR, &
&                 IB2PAR,FA%CADRE(IRANGC)%SFOHYB(1,0),2_JPLIKB ,     &
&                 KVALCO(IDECAL+1),ILDISP,INMOCC,IARR,        &
&                 0_JPLIKB ,IPUILA,IERR,KVALCO(IDECAL-1),     &
&                 KVALCO(IDECAL),LLARPE)
  ELSE
    CALL FACODEGA(FA%ICHAMP,ILCHAM,INBITS,FA%NBIMAC,KB1PAR,IB2PAR, &
&                 FA%CADRE(IRANGC)%SFOHYB(1,0),2_JPLIKB ,                  &
&                 KVALCO(IDECAL+1),ILDISP,INMOCC,IARR,ICPACK,       &
&                 IPUILA,IERR,KVALCO(IDECAL-1),KVALCO(IDECAL),      &
&                 LLARPE)
  ENDIF
#ifndef RS6K
  IF (ASSOCIATED( FA%ICHAMP )) DEALLOCATE ( FA%ICHAMP )
#endif
!
  IF (IERR.NE.0) THEN
    KREP=-200+IERR
    GOTO 1001
  ELSEIF (LDCOSP) THEN
    KVALCO(4)=ICPACK
    KVALCO(5)=IPUILA
!
    IF (LLARPE) THEN
!*
!     3.2 -  TRANSFERT DES COEFFICIENTS SPECTRAUX NON COMPACTES.
!-----------------------------------------------------------------------
!        (et non traites par CODEGA) en fin d'article.
!
     IF (LLMLAM) THEN
! Copy nonpacked part of kchamp (sub-triangle+axes) into ichaux
       IIND=0
!
       DO JTRON=0,ITRONC
       IL=2+2*JTRON+1
       ILOW=FA%CADRE(IRANGC)%NOZPAR(IL)
!
       IF (JTRON.EQ.0) THEN
         IHIGH=FA%CADRE(IRANGC)%NOZPAR(IL+1)
       ELSE
         IHIGH=ILOW+4*(ICPACK+1-JTRON)-1
         IF (IHIGH.LE.ILOW) IHIGH=ILOW+3
       ENDIF
!
       DO J=ILOW,IHIGH
       IIND=IIND+1
       FA%ICHAUX(IIND)=KCHAMP(J)
       ENDDO
       ENDDO
!
       DO J=1,ISTRIA
       KVALCO(IDECAL+INMOCC+J)=FA%ICHAUX(J)
       ENDDO
!
#ifndef RS6K
       IF (ASSOCIATED( FA%ICHAUX )) DEALLOCATE ( FA%ICHAUX )
#endif
!
      ELSE
!
        DO J=1,IDIMNC
        KVALCO(IDECAL+INMOCC+J)=KCHAMP(J)
        ENDDO
!
      ENDIF
!
    ENDIF
!
  ENDIF
!
  KVALCO(3)=INBITS
!
  IF (LLMLAM) THEN
     KLONGD=IDECAL+INMOCC+ISTRIA
  ELSE
     KLONGD=IDECAL+INMOCC+IDIMNC
  ENDIF
!
ENDIF
!
KVALCO(1)=IVALC1
KVALCO(2)=IVALC2
IF (ASSOCIATED( FA%ICHAUX )) DEALLOCATE ( FA%ICHAUX )
IF (ASSOCIATED( FA%ICHAMP )) DEALLOCATE ( FA%ICHAMP )
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
  CLNSPR='FACINE'
  INUMER=JPNIIL
!
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4, &
&         '', CDNOMA='''''',A,'''''', LDCOSP= '',L1,      &
&         '', KLONGD='',I8)')                             &
&     KREP, KRANG, CDNOMA, LDCOSP, KLONGD
  CALL FAIPAR_MT64                                      &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR, CDNOMA,.FALSE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FACINE_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FACINE64                                     &
&           (KREP, KRANG, CDNOMA, KCHAMP, LDCOSP, KVALCO, &
&           KLONGD, KB1PAR)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDNOMA                                 ! IN   
INTEGER (KIND=JPDBLE)  KCHAMP     (*)                         ! IN   
LOGICAL                LDCOSP                                 ! IN   
INTEGER (KIND=JPDBLE)  KVALCO     (*)                         !   OUT
INTEGER (KIND=JPLIKB)  KLONGD                                 !   OUT
INTEGER (KIND=JPLIKB)  KB1PAR     (FA%JPLB1P)                 ! INOUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FACINE_MT64                                            &
&           (FA, KREP, KRANG, CDNOMA, KCHAMP, LDCOSP, KVALCO, &
&           KLONGD, KB1PAR)

END SUBROUTINE

SUBROUTINE FACINE                                       &
&           (KREP, KRANG, CDNOMA, KCHAMP, LDCOSP, KVALCO, &
&           KLONGD, KB1PAR)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDNOMA                                 ! IN   
INTEGER (KIND=JPDBLE)  KCHAMP     (*)                         ! IN   
LOGICAL                LDCOSP                                 ! IN   
INTEGER (KIND=JPDBLE)  KVALCO     (*)                         !   OUT
INTEGER (KIND=JPLIKM)  KLONGD                                 !   OUT
INTEGER (KIND=JPLIKM)  KB1PAR     (FA%JPLB1P)                 ! INOUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FACINE_MT                                              &
&           (FA, KREP, KRANG, CDNOMA, KCHAMP, LDCOSP, KVALCO, &
&           KLONGD, KB1PAR)

END SUBROUTINE

SUBROUTINE FACINE_MT                                        &
&           (FA, KREP, KRANG, CDNOMA, KCHAMP, LDCOSP, KVALCO, &
&           KLONGD, KB1PAR)
USE FA_MOD, ONLY : FA_COM
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
CHARACTER (LEN=*)      CDNOMA                                 ! IN   
INTEGER (KIND=JPDBLE)  KCHAMP     (*)                         ! IN   
LOGICAL                LDCOSP                                 ! IN   
INTEGER (KIND=JPDBLE)  KVALCO     (*)                         !   OUT
INTEGER (KIND=JPLIKM)  KLONGD                                 !   OUT
INTEGER (KIND=JPLIKM)  KB1PAR     (FA%JPLB1P)                 ! INOUT
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  IRANG                                  ! IN   
INTEGER (KIND=JPLIKB)  ILONGD                                 !   OUT
INTEGER (KIND=JPLIKB)  IB1PAR     (FA%JPLB1P)                 ! INOUT
! Convert arguments

IRANG      = INT (     KRANG, JPLIKB)
IB1PAR     = INT (    KB1PAR, JPLIKB)

CALL FACINE_MT64                                            &
&           (FA, IREP, IRANG, CDNOMA, KCHAMP, LDCOSP, KVALCO, &
&           ILONGD, IB1PAR)

KREP       = INT (      IREP, JPLIKM)
KLONGD     = INT (    ILONGD, JPLIKM)
KB1PAR     = INT (    IB1PAR, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                                                              
!INTF KRANG         IN                                                                 
!INTF CDNOMA        IN                                                                 
!INTF KCHAMP        IN    DIMS=*                         KIND=JPDBLE                   
!INTF LDCOSP        IN                                                                 
!INTF KVALCO          OUT DIMS=*                         KIND=JPDBLE                   
!INTF KLONGD          OUT                                                              
!INTF KB1PAR        INOUT DIMS=FA%JPLB1P                                               
