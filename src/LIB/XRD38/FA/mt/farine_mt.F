C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FARINE_MT (FA,  KOPTIO )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme est charge des INITIALISATIONS du logiciel
C     de Fichiers ARPEGE FA ( Routine d'INitialisation )
C**
C        Argument : KOPTIO  ==> OPTION concernant le mode d'utilisation.
C                  (Entree)     (MULTI-TACHE ou NON)
C     VALEURS POSSIBLES : 0 ==> Mode MONO-Tache prescrit;
C                         1 ==> Mode MULTI-Taches prescrit;
C                         2 ==> Utilisation du mode par defaut si c'est
C                               le premier appel; sinon on garde le mode
C                               prescrit anterieurement .
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KOPTIO
C
      INTEGER IREP, J, JN, IDEBUT, IFIN, INDM, INDN, INIMES, INUMER
      INTEGER JJPUIS, INCHIF
C
      LOGICAL LLNMUL, LLASGN, LLREL
C
C
#include "facom2.h"
#include "facom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FARINE_MT',0,ZHOOK_HANDLE)
      IF (KOPTIO.LT.0.OR.KOPTIO.GT.2) THEN
        IREP=-52
        GOTO 1001
      ENDIF
C
      IF (FA%FARINE_LLPREA) THEN
C
C         C'EST LE PREMIER APPEL AU sous-programme - INITIALISATIONS .
C
        FA%NFIOUV=0
        FA%NCADEF=0
        FA%NRFAGA=1
        FA%NIMSGA=1
        FA%NBIPDG=24
        FA%NBICSP=24
        FA%NPUILA=1
        FA%NSTROI=10
        FA%NMIDPL=5
        FA%NIGRIB=2
        FA%LFAMOP=.FALSE.
        FA%LERRFA(0)=.TRUE.
        FA%NIVOMS(0)=0
        FA%NCPCAD=LEN (FA%CNOMCA(1))
        FA%SPSMIN=45000.
        FA%SPSMAX=110000.
        FA%MPRESX=100000
        FA%NBIMAC=64
        FA%NBIMAX=31
        FA%LIGARD=.FALSE.
        FA%NTYPTX=2
        FA%NXNIVV=FA%JPXNIV
        FA%NXTRON=FA%JPXTRO
        FA%NXLATI=FA%JPXLAT
        FA%NXLONG=FA%JPXLON
C
        DO 101 J=1,FA%JPNXFA
          FA%NULOGI(J)=FA%JPNIIL
          FA%NRASHO(J)=0
          FA%NRASVE(J)=0
  101   CONTINUE
C
        DO 102 J=1,FA%JPNXCA
          FA%CNOMCA(J)=' '
  102   CONTINUE
C
        DO 1020 J=1,FA%JPXNOM
          FA%CHAINC(J:J)='?'
 1020   CONTINUE
C
C        Descripteurs lies aux types de niveaux verticaux:
C        prefixes reconnus, extrema possibles du niveau dans le cas
C        d'une coordonnee verticale, premiers elements de xB1PAR.
C
        DO 111 J=0,FA%JPTNIV
C
C       Initialisation par defaut pour la serie d'affectations qui suit.
C
          FA%NIVDSC(0,J)=0
          FA%NIVDSC(1,J)=0
          FA%NIVDSC(2,J)=0
C
C        Tous les cas sont reputes "niveau vrai" (pas des couches).
C
          FA%NIVDSC(4,J)=0
  111   CONTINUE
C
C           Cas du type non reconnu.
C
        FA%NIVDSC(0,0)=0
        FA%NIVDSC(3,0)=200
C
C           Niveau hybride.
C
        FA%CTNPRF(1)='S'
        FA%NIVDSC(0,1)=3
        FA%NIVDSC(2,1)=FA%JPXNIV
        FA%NIVDSC(3,1)=109
C
C           Niveau isobare.
C
        FA%CTNPRF(2)='P'
        FA%NIVDSC(0,2)=5
        FA%NIVDSC(1,2)=1
        FA%NIVDSC(2,2)=10**5
        FA%NIVDSC(3,2)=100
C
C           Niveau iso-hauteur (au-dessus d'un relief de reference).
C
        FA%CTNPRF(3)='H'
        FA%NIVDSC(0,3)=5
        FA%NIVDSC(2,3)=10**5-1
        FA%NIVDSC(3,3)=105
C
C           Niveau iso-tourbillon_potentiel.
C
        FA%CTNPRF(4)='V'
        FA%NIVDSC(0,4)=3
        FA%NIVDSC(2,4)=10**3-1
        FA%NIVDSC(3,4)=117
C
C           Niveau iso-temperature_potentielle.
C
        FA%CTNPRF(5)='T'
        FA%NIVDSC(0,5)=3
        FA%NIVDSC(2,5)=10**3-1
        FA%NIVDSC(3,5)=113
C
C           Niveau surface.
C
        FA%CTNPRF(6)='SURF'
        FA%NIVDSC(3,6)=1
C
C           Niveau de vent max (jet).
C
        FA%CTNPRF(7)='JET'
        FA%NIVDSC(3,7)=6
C
C           Niveau surface.
C
        FA%CTNPRF(8)='TROPO'
        FA%NIVDSC(3,8)=7
C
C           Niveau moyen de la mer.
C
        FA%CTNPRF(9)='MER'
        FA%NIVDSC(3,9)=102
C
C           Niveau hybride bis pour MOCAGE
C
        FA%CTNPRF(10)='L'
        FA%NIVDSC(0,10)=3
        FA%NIVDSC(2,10)=FA%JPXNIV
C  
C           Niveau iso-temperature.
C
        FA%CTNPRF(11)='KB'
        FA%NIVDSC(0,11)=3
        FA%NIVDSC(2,11)=10**3-1
        FA%NIVDSC(3,11)=113
C
        FA%CTNPRF(12)='KT'
        FA%NIVDSC(0,12)=3
        FA%NIVDSC(2,12)=10**3-1
        FA%NIVDSC(3,12)=113

C
C  Initialisations pour la mise en oeuvre de GRIBEX
C
C  1/ On force GRIBEX a calculer la puissance de laplacien
C        CALL GRSMKP(1)
C  2/ On retire l'arrondi du message GRIB a un multiple de 120 octets
C        CALL GRSRND(0)
C
C  3/ Creation de la correspondance "nom article FA" et
C                                   "descripteurs GRIBEX"
        CALL FAICOR_MT (FA)
C
C
C  4/ Definition du codage GRIBEX par defaut:
C
C  Il s'agit d'une compression "APAC1" (meilleure solution entre
C  l'absence de compression et la compression ligne a ligne)
C  associee a la compression "general extended", a la differentiation
C  spatiale (-1: calcul dynamique de l'ordre) et au rearrangement
C  boustrophedonique.
C
C  "APAC1"
        FA%NCODGRI(1) = 1
        FA%NCODGRI(2) = 16
        FA%NCODGRI(3) = 0
        FA%NCODGRI(4) = 0
        FA%NCODGRI(5) = 16
        FA%NCODGRI(6) = 0
C  compression "general extended"
        FA%NCODGRI(7) = 8
C  arrangement boustrophedonique
        FA%NCODGRI(8) = 4
C  differentiation spatiale
        FA%NCODGRI( 9)= 0
        FA%NCODGRI(10)= -1
C
C  5/ Initialisation de logiques pilotant des initialisations ulterieures
C
C  Il faudra initialiser XLAPxDx
        FA%LIXLAP=.TRUE.
C  Il ne faut pas initialiser FLAP1Dx(): on attend l'ouverture du fichier
        FA%LIFLAP(:)=.FALSE.
C  Il faudra initialiser le tableau FA%NSEC1 (section 1 GRIBEX) via FAISC1
        FA%LISEC1(:)=.TRUE.
C  Il faudra initialiser les tableaux NSEC2xxx et FA%XSEC2
C  (section 2 GRIBEX) via FAISC2
        FA%LISEC2(:)=.TRUE.
C  Il faudra initialiser le tableau FA%NSC2ALF
C  (section 2 GRIBEX) via FAIS2F
        FA%LISC2F(:)=.TRUE.
C
        FA%FARINE_LLPREA=.FALSE.
        LLNMUL=(KOPTIO.EQ.1).OR.(KOPTIO.EQ.2.AND.FA%FARINE_LLDEFM)
        LLASGN=LLNMUL
        LLREL=.FALSE.
        CALL LFIINI_MT(FA%LFI, KOPTIO)
C
      ELSEIF (KOPTIO.EQ.2) THEN
C
C         CE N'EST PAS LE PREMIER APPEL, MAIS COMME L'ARGUMENT VAUT 2,
C         ON LAISSE LES CHOSES EN PLACE .
C
        LLNMUL=FA%LFAMUL
        LLASGN=.FALSE.
        LLREL =.FALSE.
C
      ELSE
        LLNMUL=KOPTIO.EQ.1
        LLASGN=LLNMUL.AND.(.NOT.FA%LFAMUL)
        LLREL =(.NOT.LLNMUL).AND.FA%LFAMUL
C
C     CE N'EST PAS LE PREMIER APPEL ET LE MODE EST PASSE 'EXPLICITEMENT'
C
        IF ((LLASGN.OR.LLREL).AND.FA%NFIOUV.NE.0) THEN
          IREP=-54
          GOTO 1001
        ENDIF
C
        CALL LFIINI_MT(FA%LFI, KOPTIO)
C
      ENDIF
C
      FA%LFAMUL=LLNMUL
      IREP=0
C
      IF (LLASGN) THEN
        CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ASGN')
      ELSEIF (LLREL) THEN
        CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'REL')
      ENDIF
C
 1001 CONTINUE
C
C        MESSAGERIE EVENTUELLE, AVEC ABORT SI NECESSAIRE .
C
      LLFATA=IREP.NE.0.AND.FA%NRFAGA.NE.2
C
      IF (LLFATA.OR.FA%LFAMOP) THEN
        INIMES=2
      ELSEIF (IREP.NE.0) THEN
        INIMES=0
      ELSEIF (FA%NIMSGA.EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('FARINE_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FARINE'
      INUMER=FA%JPNIIL
C
      IF (MAX (INIMES,FA%NIMSGA).EQ.2) THEN
        WRITE (UNIT=CLMESS,
     S         FMT='(''KOPTIO='',I5,'', CODE INTERNE='',I4)'
     S         ) KOPTIO,IREP
        IF (INIMES.NE.2) CALL FAIPAR_MT (FA, INUMER,FA%NIMSGA,IREP,
     S                                .FALSE.,CLMESS,
     S                                CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S             CLNSPR,CLACTI,.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FARINE_MT',1,ZHOOK_HANDLE)
      END

