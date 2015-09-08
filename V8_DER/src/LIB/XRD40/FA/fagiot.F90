! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FAGIOT_MT64                                            &
&                     (FA,  KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA, &
&                    KDMOPL )
USE FA_MOD, ONLY : FA_COM, JPNIIL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!        Ce sous-programme permet de modifier les options implicites
!     liees au codage GRIB des champs.
!     CES OPTIONS NE SONT UTILISEES QUE POUR (RE)ECRIRE DES CHAMPS
!     codes en GRIB, et les nouvelles valeurs implicites ne serviront
!     que LORS d'une OUVERTURE de FICHIER ULTERIEURE.
!       ( Grib, Implicites Options Techniques )
!**
!     Arguments : KNGRIB (Entree) ==> Niveau de codage GRIB (-1,0,1,2,3);
!                 KNBPDG (Entree) ==> Nombre de bits par valeur point-
!                                     de-grille;
!                 KNBCSP (Entree) ==> Nombre de bits par partie reelle/
!                                     imaginaire de coeff. spectral;
!                 KSTRON (Entree) ==> Sous-troncature non compactee;
!                 KPUILA (Entree) ==> Puissance de laplacien;
!                 KDMOPL (Entree) ==> Degre de modulation de KPUILA.
!
!     N.B.: Il doit y avoir coherence vis-a-vis des cadres deja definis
!           et vis-a-vis des limites usagers.
!           ( ce qui en pratique, ne concerne que KSTRON )
!
!     Remarque:  KSTRON egal a -1 est accepte, et dans ce cas
!                on indexera (a chaque ouverture de fichier) la sous-
!                troncature effective sur la troncature.
!
!      MODIF 30/03/2007 JM AUDOIN FA%LFAMOP pour limiter IMPRESSION
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KNGRIB, KNBPDG, KNBCSP
INTEGER (KIND=JPLIKB) KSTRON, KPUILA, KDMOPL
!
INTEGER (KIND=JPLIKB) IMINIM, IREP, INIMES
INTEGER (KIND=JPLIKB) INUMER, J, IRANGC
!
LOGICAL LLVERG
!
!
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
IF (LHOOK) CALL DR_HOOK('FAGIOT_MT',0,ZHOOK_HANDLE)
CLACTI=''
IF (FA%FAGIOT_LLPREA) THEN
!
!          A la premiere utilisation, appel au sous-programme "FARINE".
!
  CALL FARINE_MT64             &
&                (FA, 2_JPLIKB )
  FA%FAGIOT_LLPREA=.FALSE.
ENDIF
!
LLVERG=.FALSE.
IMINIM=MIN (2+KNGRIB,2+KNBPDG,2+KNBCSP,2+KSTRON,1+KDMOPL)
!
IF (IMINIM.LE.0) THEN
  IREP=-64
  GOTO 1001
ELSEIF (KNBPDG*KNBCSP.EQ.0) THEN
  IREP=-64
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
!         Verrouillage global eventuel.
!
IF (FA%LFAMUL) CALL LFIVER_MT64                        &
&                              (FA%LFI, FA%VRGLAS,'ON')
LLVERG=FA%LFAMUL
!
IF (KSTRON.GE.FA%NXTRON) THEN
  IREP=-113
  GOTO 1001
ENDIF
!
!        Coherence de "KSTRON" vis-a-vis de la troncature des cadres
!     deja definis.
!
DO J=1,FA%NCADEF
IRANGC=FA%NCAIND(J)
!
IF (KSTRON.GE.FA%CADRE(IRANGC)%MTRONC) THEN
  IREP=-99
  GOTO 1001
ENDIF
!
ENDDO
!**
!     2.  -  STOCKAGE DES NOUVEAUX PARAMETRES.
!-----------------------------------------------------------------------
!
IF (FA%LFAMOP.AND.(FA%NIGRIB.EQ.-1.OR.FA%NIGRIB.EQ.3).AND. &
&    (KNGRIB.GT.-1.AND.KNGRIB.LT.3))     THEN
  WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
  WRITE (UNIT=FA%NULOUT,FMT=*)                                   &
&          'FAGIOT: CAUTION!! Les champs spectraux ARPEGE ne',    &
&          ' devront pas etre ranges comme dans le modele ARPEGE'
  WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
ENDIF
IF (FA%LFAMOP.AND.(KNGRIB.EQ.-1.OR.KNGRIB.EQ.3).AND. &
&    (FA%NIGRIB.GT.-1.AND.FA%NIGRIB.LT.3))     THEN
  WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
  WRITE (UNIT=FA%NULOUT,FMT=*)                                     &
&          'FAGIOT: CAUTION!! Les champs spectraux ARPEGE devront', &
&          ' etre ranges comme dans le modele ARPEGE'
  WRITE (UNIT=FA%NULOUT,FMT=*)'-----------------'
ENDIF
FA%NIGRIB=KNGRIB
FA%NBIPDG=KNBPDG
FA%NBICSP=KNBCSP
FA%NSTROI=KSTRON
FA%NPUILA=KPUILA
FA%NMIDPL=KDMOPL
IREP=0
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
LLFATA=LLMOER (IREP,0_JPLIKB )
!
!        Deverrouillage global eventuel.
!
IF (LLVERG) CALL LFIVER_MT64                         &
&                           (FA%LFI, FA%VRGLAS,'OFF')
!
IF (LLFATA) THEN
  INIMES=2
ELSE
  INIMES=FA%NIMSGA
ENDIF
!
IF (.NOT.LLFATA.AND.INIMES.NE.2)  THEN 
  IF (LHOOK) CALL DR_HOOK('FAGIOT_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='FAGIOT'
!
WRITE (UNIT=CLMESS,FMT='(''KNGRIB='',I2,'', KNBPDG='',I3,  &
&       '', KNBCSP='',I3,'', KSTRON='',I2,'', KPUILA='',I3, &
&       '', KDMOPL='',I3)')                                 &
&   KNGRIB,KNBPDG,KNBCSP,KSTRON,KPUILA,KDMOPL
INUMER=JPNIIL
CALL FAIPAR_MT64                                     &
&               (FA, INUMER,INIMES,IREP,LLFATA,CLMESS, &
&                CLNSPR,CLACTI,.FALSE.)
!
IF (LHOOK) CALL DR_HOOK('FAGIOT_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FAGIOT64                                        &
&           (KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA, KDMOPL)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KNGRIB                                 ! IN   
INTEGER (KIND=JPLIKB)  KNBPDG                                 ! IN   
INTEGER (KIND=JPLIKB)  KNBCSP                                 ! IN   
INTEGER (KIND=JPLIKB)  KSTRON                                 ! IN   
INTEGER (KIND=JPLIKB)  KPUILA                                 ! IN   
INTEGER (KIND=JPLIKB)  KDMOPL                                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAGIOT_MT64                                               &
&           (FA, KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA, KDMOPL)

END SUBROUTINE

SUBROUTINE FAGIOT                                          &
&           (KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA, KDMOPL)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KNGRIB                                 ! IN   
INTEGER (KIND=JPLIKM)  KNBPDG                                 ! IN   
INTEGER (KIND=JPLIKM)  KNBCSP                                 ! IN   
INTEGER (KIND=JPLIKM)  KSTRON                                 ! IN   
INTEGER (KIND=JPLIKM)  KPUILA                                 ! IN   
INTEGER (KIND=JPLIKM)  KDMOPL                                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAGIOT_MT                                                 &
&           (FA, KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA, KDMOPL)

END SUBROUTINE

SUBROUTINE FAGIOT_MT                                           &
&           (FA, KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA, KDMOPL)
USE FA_MOD, ONLY : FA_COM
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KNGRIB                                 ! IN   
INTEGER (KIND=JPLIKM)  KNBPDG                                 ! IN   
INTEGER (KIND=JPLIKM)  KNBCSP                                 ! IN   
INTEGER (KIND=JPLIKM)  KSTRON                                 ! IN   
INTEGER (KIND=JPLIKM)  KPUILA                                 ! IN   
INTEGER (KIND=JPLIKM)  KDMOPL                                 ! IN   
! Local integers
INTEGER (KIND=JPLIKB)  INGRIB                                 ! IN   
INTEGER (KIND=JPLIKB)  INBPDG                                 ! IN   
INTEGER (KIND=JPLIKB)  INBCSP                                 ! IN   
INTEGER (KIND=JPLIKB)  ISTRON                                 ! IN   
INTEGER (KIND=JPLIKB)  IPUILA                                 ! IN   
INTEGER (KIND=JPLIKB)  IDMOPL                                 ! IN   
! Convert arguments

INGRIB     = INT (    KNGRIB, JPLIKB)
INBPDG     = INT (    KNBPDG, JPLIKB)
INBCSP     = INT (    KNBCSP, JPLIKB)
ISTRON     = INT (    KSTRON, JPLIKB)
IPUILA     = INT (    KPUILA, JPLIKB)
IDMOPL     = INT (    KDMOPL, JPLIKB)

CALL FAGIOT_MT64                                               &
&           (FA, INGRIB, INBPDG, INBCSP, ISTRON, IPUILA, IDMOPL)


END SUBROUTINE

!INTF KNGRIB        IN    
!INTF KNBPDG        IN    
!INTF KNBCSP        IN    
!INTF KSTRON        IN    
!INTF KPUILA        IN    
!INTF KDMOPL        IN    
