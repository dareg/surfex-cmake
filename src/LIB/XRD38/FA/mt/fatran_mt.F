C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FATRAN_MT (FA,  KREP,  KNUMER,  PCHAME, PCHAMS, LDOPT )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme du logiciel de Fichiers ARPEGE permettant la
C      TRANsposition d'un champ spectral ARPEGE ou ALADIN,
C      d'un rangement des coeff selon le MODELE vers un rangement
C      des coeff selon FA+GRIB_version0 et inversement.
C
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Rang de l'unite logique;
C    ( Tableau ) PCHAME (Entree) ==> Valeurs du champ a transposer;
C    ( Tableau ) PCHAMS (Sortie) ==> Valeurs du champ transpose;
C                LDOPT   (Entree) ==> Option de transposition;
C                                    si .TRUE.  alors PCHAME range comme MODELE
C                                                     (soit "verticalement")
C                                                     PCHAMS range comme FA-GRIB0
C                                                     (soit "horizontalement")
C                                    si .FALSE. alors PCHAME range comme FA-GRIB0
C                                                     PCHAMS range comme MODELE
C*
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER
C
      LOGICAL LDOPT
C
      REAL (KIND=JPDBLR) PCHAME(*), PCHAMS(*)
C
      INTEGER JN, JM, J, INDEX, ILOW, IHIGH, IRANGC, IRANG
      INTEGER INIMES, ITRONC, IMSMAX
      INTEGER, ALLOCATABLE :: IND(:,:)
C
      LOGICAL LLMLAM
C
#include "facom2.h"
#include "facom_mt.h"
C
C
C
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FATRAN_MT',0,ZHOOK_HANDLE)
      KREP=0
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        KREP=-51
        GOTO 1001
      ENDIF
      IRANGC=FA%NUCADR(IRANG)
      LLMLAM=FA%LIMLAM(IRANGC)
      ITRONC=FA%MTRONC(IRANGC)
      IF (LLMLAM) THEN
        IMSMAX = FA%NOMPAR(2,IRANGC)
      ENDIF
C
C Initialisation de l'indirection pour IRANGC, si ce n'est pas deja fait.
C Si ARPEGE, FA%FATRAN_INDIRECT(J,IRANGC)=INDEX signifie que les indices J dans
C   le tableau "FA+GRIB0" et INDEX dans le tableau "modele ARPEGE"
C    designent un meme coeff spectral.
C Si Aladin, FA%FATRAN_INDIRECT(JM*(ITRONC+1)+JN+1,IRANGC)=J signifie que l'indice
C   J dans le tableau "FA+GRIB0" est le premier coeff associe au couple (JM,JN)
C   ou JM est le nombre d'onde zonal et JN le nombre d'onde meridien.
C   4 coeff spectraux sont associes a chaque couple (JM,JN) car JM varie
C   de 0 a IMSMAX et JN varie de 0 a ITRONC (soit 1/4 de l'ellipse).
C
C  CAS ARPEGE
C
      IF (FA%FATRAN_LLINDIR(IRANGC).AND..NOT.LLMLAM) THEN
        ALLOCATE (IND(0:FA%JPXTRO,-FA%JPXTRO:FA%JPXTRO))
        DO JN=0,ITRONC
          ILOW=JN**2+1
          IHIGH=(JN+1)**2
          JM=-JN-1
          DO J=ILOW,IHIGH
            JM=JM+1
            IND(JN,JM)=J
          ENDDO
        ENDDO
C
        INDEX=-1
        DO JM=0,ITRONC
        DO JN=JM,ITRONC
          INDEX=INDEX+2
          FA%FATRAN_INDIRECT  (IND(JN, JM),IRANGC)=INDEX
          IF (JM.NE.0) THEN
            FA%FATRAN_INDIRECT(IND(JN,-JM),IRANGC)=INDEX+1
          ENDIF
        ENDDO
        ENDDO
        FA%FATRAN_LLINDIR(IRANGC)=.FALSE.
        DEALLOCATE (IND)
      ENDIF
C
C  CAS ALADIN
C
      IF (FA%FATRAN_LLINDIR(IRANGC).AND.LLMLAM) THEN
        DO JN=0,ITRONC
        DO J=FA%NOZPAR(2*JN+3,IRANGC), FA%NOZPAR(2*JN+4,IRANGC), 4
          JM=(J-FA%NOZPAR(2*JN+3,IRANGC)) / 4
          FA%FATRAN_INDIRECT(JM*(ITRONC+1)+JN+1,IRANGC) = J
        ENDDO
        ENDDO
        FA%FATRAN_LLINDIR(IRANGC)=.FALSE.
      ENDIF
C**
C     2.  -  TRANSPOSITION DES DONNEES
C-----------------------------------------------------------------------
C
C  CAS ALADIN
C
      IF (LLMLAM) THEN
          IF (LDOPT) THEN
C PCHAME range comme MODELE (soit "verticalement")
C
          DO JM=0,IMSMAX
          DO INDEX=FA%NOMPAR(2*JM+3,IRANGC), FA%NOMPAR(2*JM+4,IRANGC), 4
            JN = (INDEX-FA%NOMPAR(2*JM+3,IRANGC)) / 4
            J  = FA%FATRAN_INDIRECT(JM*(ITRONC+1)+JN+1,IRANGC)
            PCHAMS(J  )=PCHAME(INDEX  )  
            PCHAMS(J+1)=PCHAME(INDEX+1)
            PCHAMS(J+2)=PCHAME(INDEX+2)
            PCHAMS(J+3)=PCHAME(INDEX+3)
          ENDDO
          ENDDO
        ELSE
C PCHAME range comme FA+GRIB0 (soit "horizontalement")
C
          DO JM=0,IMSMAX
          DO INDEX=FA%NOMPAR(2*JM+3,IRANGC), FA%NOMPAR(2*JM+4,IRANGC), 4
            JN = (INDEX-FA%NOMPAR(2*JM+3,IRANGC)) / 4
            J  = FA%FATRAN_INDIRECT(JM*(ITRONC+1)+JN+1,IRANGC)
            PCHAMS(INDEX  )=PCHAME(J  )  
            PCHAMS(INDEX+1)=PCHAME(J+1)
            PCHAMS(INDEX+2)=PCHAME(J+2)
            PCHAMS(INDEX+3)=PCHAME(J+3)
          ENDDO
          ENDDO
        ENDIF
      ELSE
C
C  CAS ARPEGE
C
C  1/ Passage du rangement des coeff. spectraux du type modele ARPEGE
C  a celui de FA associe a GRIB version0 (et pas associe a GRIBEX qui
C  reprend la structure de tableau de ARPEGE).
C
        IF (LDOPT) THEN
          DO JN=0,ITRONC
            ILOW=JN**2+1
            IHIGH=(JN+1)**2
            DO J=ILOW,IHIGH
              PCHAMS(J)=PCHAME(FA%FATRAN_INDIRECT(J,IRANGC))
            ENDDO
          ENDDO
C
C  2/ Passage du rangement des coeff. spectraux du type FA associe
C  a GRIB version0 (et pas associe a GRIBEX qui reprend la structure de
C  tableau de ARPEGE) a celui du type modele ARPEGE.
C
        ELSE
C  Initialisation de la partie "JM=0" a zero, pour y introduire
C  ensuite uniquement les coeff reels correspondant dans PCHAME
C  (les coeff imaginaires etant donc crees et mis a zero).
          PCHAMS(1:2*(ITRONC+1))=0.
C
          DO JN=0,ITRONC
            ILOW=JN**2+1
            IHIGH=(JN+1)**2
            DO J=ILOW,IHIGH
              PCHAMS(FA%FATRAN_INDIRECT(J,IRANGC))=PCHAME(J)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      LLFATA=LLMOER (KREP,IRANG)
C
      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=IXNVMS(IRANG)
      ENDIF
C
      IF (INIMES.EQ.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('FATRAN_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FATRAN'
C
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', IRANG='',I4,
     S         '', LDOPT='',L)')  KREP, IRANG, LDOPT
      CALL FAIPAR_MT (FA, KNUMER,INIMES,KREP,LLFATA,CLMESS,
     S               CLNSPR,CLNSPR,.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FATRAN_MT',1,ZHOOK_HANDLE)
      END

