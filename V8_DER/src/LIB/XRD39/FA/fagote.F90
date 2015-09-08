! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FAGOTE_MT_BB                                          &
&                     (FA,  KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP,  &
&                      KSTRON, KPUILA, KDMOPL )
USE FA_MOD, ONLY : FA_COM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!        Ce sous-programme permet d'ajuster, pour un fichier ARPEGE
!     ouvert, les options liees au codage GRIB des champs.
!     CES OPTIONS NE SONT UTILISEES QUE POUR (RE)ECRIRE DES CHAMPS
!     codes en GRIB.
!       ( Grib, Options Techniques Effectives )
!**
!     Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                 KNUMER (Entree) ==> Numero d'Unite Logique concernee;
!                 KNGRIB (Entree) ==> Niveau de codage GRIB (-1,0,1,2,3);
!                 KNBPDG (Entree) ==> Nombre de bits par valeur point-
!                                     de-grille;
!                 KNBCSP (Entree) ==> Nombre de bits par partie reelle/
!                                     imaginaire de coeff. spectral;
!                 KSTRON (Entree) ==> Sous-troncature non compactee;
!                 KPUILA (Entree) ==> Puissance de laplacien;
!                 KDMOPL (Entree) ==> Degre de modulation de KPUILA.
!
!     Remarque:  KSTRON egal a -1 est accepte; dans ce cas on indexera
!                (pour chaque champ spectral ecrit) la sous-troncature
!                effective sur la troncature et sur le nombre de bits
!                par valeur compactee.
!                 
!     MODIF : 30/03/2007 JM AUDOIN  FA%LFAMOP Pour limiter Impression
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KNUMER, KNGRIB
INTEGER (KIND=JPLIKB) KNBPDG, KNBCSP, KSTRON, KPUILA
INTEGER (KIND=JPLIKB) KDMOPL
!
INTEGER (KIND=JPLIKB) IMINIM, IREP, IRANGC
INTEGER (KIND=JPLIKB) ITRONC, INIMES, IRANG, ITYPTR
!
LOGICAL LLVERF, LLMLAM
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
IF (LHOOK) CALL DR_HOOK('FAGOTE_MT',0,ZHOOK_HANDLE)
CLACTI=''
LLVERF=.FALSE.
IMINIM=MIN (2+KNGRIB,2+KNBPDG,2+KNBCSP,2+KSTRON,1+KDMOPL)
CALL FANUMU_MT_BB                &
&               (FA, KNUMER,IRANG)
!
IF (IRANG.EQ.0) THEN
  IREP=-51
  GOTO 1001
ELSEIF (IMINIM.LE.0) THEN
  IREP=-64
  GOTO 1001
ELSEIF (KNBPDG*KNBCSP.EQ.0 .AND. KNGRIB.GT.0) THEN
  IREP=-124
  GOTO 1001
ELSEIF (KNGRIB.GT.3) THEN
  IREP=-96
  GOTO 1001
ELSEIF (MAX (KNBPDG,KNBCSP).GT.FA%NBIMAX) THEN
  IREP=-97
  GOTO 1001
ELSEIF (ABS (KPUILA).GT.2**15-1) THEN
  IREP=-98
  GOTO 1001
ENDIF
!
!         Verrouillage eventuel du fichier.
!
IF (FA%LFAMUL) CALL LFIVER_MT_BB                             &
&                              (FA%LFI, FA%FICHIER(IRANG)%VRFICH,'ON')
LLVERF=FA%LFAMUL
!
IRANGC=FA%FICHIER(IRANG)%NUCADR
ITRONC=FA%CADRE(IRANGC)%MTRONC
ITYPTR=FA%CADRE(IRANGC)%NTYPTR
LLMLAM=FA%CADRE(IRANGC)%LIMLAM
!
IF (KSTRON.GE.ITRONC) THEN
  IREP=-99
  GOTO 1001
ELSEIF (ITYPTR.LT.0.AND.KSTRON.GE.(-ITYPTR)) THEN
  IREP=-99
  GOTO 1001
ENDIF
!**
!     2.  -  STOCKAGE DES NOUVEAUX PARAMETRES.
!-----------------------------------------------------------------------
!
IF (KPUILA.NE.FA%FICHIER(IRANG)%NPUFLA) THEN
  FA%FICHIER(IRANG)%NPUFLA=KPUILA
  FA%FICHIER(IRANG)%LIFLAP=.TRUE.
ENDIF
IF (KSTRON.NE.FA%FICHIER(IRANG)%NSTROF) THEN
  FA%FICHIER(IRANG)%NSTROF=KSTRON
  FA%FICHIER(IRANG)%LISC2F=.TRUE.
ENDIF
!
IF (FA%LFAMOP.AND.(FA%FICHIER(IRANG)%NFGRIB.EQ.3              &
&   .OR.FA%FICHIER(IRANG)%NFGRIB.EQ.-1)                        &
&  .AND.(KNGRIB.LT.3.AND.KNGRIB.GT.-1))           THEN
  WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
  WRITE (UNIT=FA%NULOUT,FMT=*)                                  &
&        'FAGOTE: WARNING!! Les champs spectraux NE devront',    &
&        ' PAS etre ranges comme dans le modele (rangt horiz.)', &
&        ' pour l''unite logique ',KNUMER
  WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
ENDIF
IF (FA%LFAMOP.AND.(FA%FICHIER(IRANG)%NFGRIB.LT.3             &
&   .AND.FA%FICHIER(IRANG)%NFGRIB.GT.-1)                      &
&  .AND.(KNGRIB.EQ.3.OR.KNGRIB.EQ.-1))           THEN
  WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
  WRITE (UNIT=FA%NULOUT,FMT=*)                                   &
&        'FAGOTE: WARNING!! Les champs spectraux devront',        &
&        ' etre ranges comme dans le modele (rangt verti.) pour', &
&        ' l''unite logique ',KNUMER
  WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
ENDIF
FA%FICHIER(IRANG)%NFGRIB=KNGRIB
FA%FICHIER(IRANG)%NBFPDG=KNBPDG
FA%FICHIER(IRANG)%NBFCSP=KNBCSP
FA%FICHIER(IRANG)%NMFDPL=KDMOPL
IREP=0
!
! Appel a FAINOC pour interpreter les eventuels defauts
! de -1 pris par FA%NBFPDG, FA%NBFCSP, FA%NSTROF et FA%NPUFLA en
! IRANG-ieme position.
!
CALL FAINOC_MT_BB           &
&               (FA,  IRANG )
!
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
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
  IF (LHOOK) CALL DR_HOOK('FAGOTE_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='FAGOTE'
!
!***** FAZZZZ - KREP=iiii, KNUMER=iii, KNGRIB=ii, KNBPDG=iii,   *****
!*****          KNBCSP=iii, KSTRON=ii, KPUILA=iii, KDMOPL=iii  *****
WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,      &
&       '', KNGRIB='',I2,'', KNBPDG='',I3,'', KNBCSP='',I3,   &
&       '', KSTRON='',I2,'', KPUILA='',I3,'', KDMOPL='',I3)') &
&   KREP,KNUMER,KNGRIB,KNBPDG,KNBCSP,KSTRON,KPUILA,KDMOPL
CALL FAIPAR_MT_BB                                    &
&               (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS, &
&                CLNSPR, CLACTI,.FALSE.)
!
IF (LHOOK) CALL DR_HOOK('FAGOTE_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"
#include "facom2.ixnvms.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FAGOTE_BB                                     &
&           (KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP, KSTRON, &
&           KPUILA, KDMOPL)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  KNGRIB                                 ! IN   
INTEGER (KIND=JPLIKB)  KNBPDG                                 ! IN   
INTEGER (KIND=JPLIKB)  KNBCSP                                 ! IN   
INTEGER (KIND=JPLIKB)  KSTRON                                 ! IN   
INTEGER (KIND=JPLIKB)  KPUILA                                 ! IN   
INTEGER (KIND=JPLIKB)  KDMOPL                                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAGOTE_MT_BB                                            &
&           (FA, KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP, KSTRON, &
&           KPUILA, KDMOPL)

END SUBROUTINE

SUBROUTINE FAGOTE_MM                                     &
&           (KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP, KSTRON, &
&           KPUILA, KDMOPL)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KNGRIB                                 ! IN   
INTEGER (KIND=JPLIKM)  KNBPDG                                 ! IN   
INTEGER (KIND=JPLIKM)  KNBCSP                                 ! IN   
INTEGER (KIND=JPLIKM)  KSTRON                                 ! IN   
INTEGER (KIND=JPLIKM)  KPUILA                                 ! IN   
INTEGER (KIND=JPLIKM)  KDMOPL                                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAGOTE_MT_MM                                            &
&           (FA, KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP, KSTRON, &
&           KPUILA, KDMOPL)

END SUBROUTINE

SUBROUTINE FAGOTE_MT_MM                                      &
&           (FA, KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP, KSTRON, &
&           KPUILA, KDMOPL)
USE FA_MOD, ONLY : FA_COM,              &
&                   FA_COM_DEFAULT_INIT, &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KNGRIB                                 ! IN   
INTEGER (KIND=JPLIKM)  KNBPDG                                 ! IN   
INTEGER (KIND=JPLIKM)  KNBCSP                                 ! IN   
INTEGER (KIND=JPLIKM)  KSTRON                                 ! IN   
INTEGER (KIND=JPLIKM)  KPUILA                                 ! IN   
INTEGER (KIND=JPLIKM)  KDMOPL                                 ! IN   
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  INGRIB                                 ! IN   
INTEGER (KIND=JPLIKB)  INBPDG                                 ! IN   
INTEGER (KIND=JPLIKB)  INBCSP                                 ! IN   
INTEGER (KIND=JPLIKB)  ISTRON                                 ! IN   
INTEGER (KIND=JPLIKB)  IPUILA                                 ! IN   
INTEGER (KIND=JPLIKB)  IDMOPL                                 ! IN   
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)
INGRIB     = INT (    KNGRIB, JPLIKB)
INBPDG     = INT (    KNBPDG, JPLIKB)
INBCSP     = INT (    KNBCSP, JPLIKB)
ISTRON     = INT (    KSTRON, JPLIKB)
IPUILA     = INT (    KPUILA, JPLIKB)
IDMOPL     = INT (    KDMOPL, JPLIKB)

CALL FAGOTE_MT_BB                                            &
&           (FA, IREP, INUMER, INGRIB, INBPDG, INBCSP, ISTRON, &
&           IPUILA, IDMOPL)

KREP       = INT (      IREP, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT 
!INTF KNUMER        IN    
!INTF KNGRIB        IN    
!INTF KNBPDG        IN    
!INTF KNBCSP        IN    
!INTF KSTRON        IN    
!INTF KPUILA        IN    
!INTF KDMOPL        IN    
