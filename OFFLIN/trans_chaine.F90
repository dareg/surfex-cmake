!----------------------------------------------------------
SUBROUTINE TRANS_CHAINE(CHAINE,ENTIER,OPTION)
!----------------------------------------------------------
!
! transforme l'entier 23416 en la chaine de caracteres '23416' 
! option est utilise si superieur a 0. Dans ce cas si l'entier
! a un nombre de digits inferieur a option, la difference est
! remplie avec des '0' dans chaine.
!
!
! ex: entier=256 option=0  ===>   chaine='256'
!     entier=256 option=2  ===>   chaine='256'
!     entier=256 option=7  ===>   chaine='0000256'

!
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE PARKIND1  ,ONLY : JPRB
!
      IMPLICIT NONE
      DOUBLE PRECISION REEL1
      CHARACTER(LEN=*) CHAINE
      INTEGER OPTION
      CHARACTER              :: C
      INTEGER                :: RESTE , NBDIGIT , NUM , ENTIER, I, J
      DOUBLE PRECISION       :: DIVI
      REAL(KIND=JPRB) :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('TRANS_CHAINE',0,ZHOOK_HANDLE)
      REEL1  = ENTIER/1.
      CHAINE = ''
      RESTE  = ENTIER
      NBDIGIT= INT(LOG(REEL1)/LOG(10.)+0.001)+1
      DIVI   = 1

      DIG_LOOP: DO I = 1 , NBDIGIT - 1
      !-------------------------------

          DIVI = DIVI * 10

      END DO DIG_LOOP
      !---------------

      IF (OPTION.GT.0) THEN
         IF (OPTION.GT.NBDIGIT) THEN

            OPT_LOOP: DO I = 1 , OPTION - NBDIGIT
            !------------------------------------

               CHAINE(I:I)='0'

            END DO OPT_LOOP
            !---------------
          
         ENDIF
      ENDIF

      DIGI_LOOP: DO I = 1 , NBDIGIT
      !----------------------------

         NUM=INT(RESTE/DIVI)
         C=CHAR(NUM+48)
         IF (OPTION.GT.0) J=OPTION-NBDIGIT+I
         IF (OPTION.EQ.0) J=I
         CHAINE(J:J)=C
         RESTE=NINT(RESTE-NUM*DIVI)
         DIVI=DIVI/10

      END DO DIGI_LOOP
IF (LHOOK) CALL DR_HOOK('TRANS_CHAINE',1,ZHOOK_HANDLE)
      !---------------

END SUBROUTINE TRANS_CHAINE
