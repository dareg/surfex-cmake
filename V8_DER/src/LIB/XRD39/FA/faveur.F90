! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FAVEUR_MT_BB                                          &
&                     (FA,  KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP,  &
&                    KSTRON, KPUILA, KDMOPL )
USE FA_MOD, ONLY : FA_COM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!        Ce sous-programme permet d'obtenir, pour un fichier ARPEGE
!     ouvert, les options courantes liees au codage GRIB des champs.
!     CES OPTIONS NE SONT UTILISEES QUE POUR (RE)ECRIRE DES CHAMPS
!     codes en GRIB.
!       ( Visualisation (?) options Effectives pour l'UtilisateuR )
!**
!     Arguments : KREP   ==> Code-reponse du sous-programme;
!     (tous de    KNUMER ==> Numero d'Unite Logique concernee;
!      SORTIE)    KNGRIB ==> Niveau de codage GRIB (-1,0,1,2,3);
!                 KNBPDG ==> Nombre de bits par valeur point-de-grille;
!                 KNBCSP ==> Nombre de bits par partie reelle/imaginaire
!                            de coefficient spectral;
!                 KSTRON ==> Sous-troncature non compactee;
!                 KPUILA ==> Puissance de laplacien;
!                 KDMOPL ==> Degre de modulation de KPUILA.
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KNUMER, KNGRIB
INTEGER (KIND=JPLIKB) KNBPDG, KNBCSP, KSTRON, KPUILA
INTEGER (KIND=JPLIKB) KDMOPL
!
INTEGER (KIND=JPLIKB) IREP, IRANG, INIMES
!
LOGICAL LLVERF
!
CHARACTER(LEN=FA%JPXNOM) CLACTI 
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     1.  -  INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FAVEUR_MT',0,ZHOOK_HANDLE)
CLACTI=''
LLVERF=.FALSE.
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
!**
!     2.  -  RECOPIE DES VALEURS EN COMMON DANS LES ARGUMENTS.
!-----------------------------------------------------------------------
!
KNGRIB=FA%FICHIER(IRANG)%NFGRIB
KNBPDG=FA%FICHIER(IRANG)%NBFPDG
KNBCSP=FA%FICHIER(IRANG)%NBFCSP
KSTRON=FA%FICHIER(IRANG)%NSTROF
KPUILA=FA%FICHIER(IRANG)%NPUFLA
KDMOPL=FA%FICHIER(IRANG)%NMFDPL
IREP=0
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
  IF (LHOOK) CALL DR_HOOK('FAVEUR_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='FAVEUR'
!
WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,      &
&       '', KNGRIB='',I2,'', KNBPDG='',I3,'', KNBCSP='',I3,   &
&       '', KSTRON='',I2,'', KPUILA='',I3,'', KDMOPL='',I3)') &
&   KREP,KNUMER,KNGRIB,KNBPDG,KNBCSP,KSTRON,KPUILA,KDMOPL
CALL FAIPAR_MT_BB                                    &
&               (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS, &
&                CLNSPR,CLACTI,.FALSE.)
!
IF (LHOOK) CALL DR_HOOK('FAVEUR_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"
#include "facom2.ixnvms.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FAVEUR_BB                                     &
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
INTEGER (KIND=JPLIKB)  KNGRIB                                 !   OUT
INTEGER (KIND=JPLIKB)  KNBPDG                                 !   OUT
INTEGER (KIND=JPLIKB)  KNBCSP                                 !   OUT
INTEGER (KIND=JPLIKB)  KSTRON                                 !   OUT
INTEGER (KIND=JPLIKB)  KPUILA                                 !   OUT
INTEGER (KIND=JPLIKB)  KDMOPL                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAVEUR_MT_BB                                            &
&           (FA, KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP, KSTRON, &
&           KPUILA, KDMOPL)

END SUBROUTINE

SUBROUTINE FAVEUR_MM                                     &
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
INTEGER (KIND=JPLIKM)  KNGRIB                                 !   OUT
INTEGER (KIND=JPLIKM)  KNBPDG                                 !   OUT
INTEGER (KIND=JPLIKM)  KNBCSP                                 !   OUT
INTEGER (KIND=JPLIKM)  KSTRON                                 !   OUT
INTEGER (KIND=JPLIKM)  KPUILA                                 !   OUT
INTEGER (KIND=JPLIKM)  KDMOPL                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAVEUR_MT_MM                                            &
&           (FA, KREP, KNUMER, KNGRIB, KNBPDG, KNBCSP, KSTRON, &
&           KPUILA, KDMOPL)

END SUBROUTINE

SUBROUTINE FAVEUR_MT_MM                                      &
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
INTEGER (KIND=JPLIKM)  KNGRIB                                 !   OUT
INTEGER (KIND=JPLIKM)  KNBPDG                                 !   OUT
INTEGER (KIND=JPLIKM)  KNBCSP                                 !   OUT
INTEGER (KIND=JPLIKM)  KSTRON                                 !   OUT
INTEGER (KIND=JPLIKM)  KPUILA                                 !   OUT
INTEGER (KIND=JPLIKM)  KDMOPL                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  INGRIB                                 !   OUT
INTEGER (KIND=JPLIKB)  INBPDG                                 !   OUT
INTEGER (KIND=JPLIKB)  INBCSP                                 !   OUT
INTEGER (KIND=JPLIKB)  ISTRON                                 !   OUT
INTEGER (KIND=JPLIKB)  IPUILA                                 !   OUT
INTEGER (KIND=JPLIKB)  IDMOPL                                 !   OUT
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)

CALL FAVEUR_MT_BB                                            &
&           (FA, IREP, INUMER, INGRIB, INBPDG, INBCSP, ISTRON, &
&           IPUILA, IDMOPL)

KREP       = INT (      IREP, JPLIKM)
KNGRIB     = INT (    INGRIB, JPLIKM)
KNBPDG     = INT (    INBPDG, JPLIKM)
KNBCSP     = INT (    INBCSP, JPLIKM)
KSTRON     = INT (    ISTRON, JPLIKM)
KPUILA     = INT (    IPUILA, JPLIKM)
KDMOPL     = INT (    IDMOPL, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT 
!INTF KNUMER        IN    
!INTF KNGRIB          OUT 
!INTF KNBPDG          OUT 
!INTF KNBCSP          OUT 
!INTF KSTRON          OUT 
!INTF KPUILA          OUT 
!INTF KDMOPL          OUT 
