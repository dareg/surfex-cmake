!
SUBROUTINE TRANS_CHAINE(HCHAINE,KENTIER,KOPTION)
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
!
INTEGER,          INTENT(IN ) :: KENTIER, KOPTION
CHARACTER(LEN=*), INTENT(OUT) :: HCHAINE
!
CHARACTER :: YC
INTEGER :: IDIVI
INTEGER :: IRESTE, INBDIGIT, INUM
INTEGER:: I, J
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('TRANS_CHAINE',0,ZHOOK_HANDLE)
!
HCHAINE = ''
!
IRESTE   = KENTIER
!
IF (KENTIER==0) THEN
  INBDIGIT = 1
ELSE
  INBDIGIT = INT( LOG(KENTIER*1.)/LOG(10.)+0.001 ) + 1
ENDIF
!
IDIVI = 10**(INBDIGIT-1)
!
IF ( KOPTION.GT.MAX(0,INBDIGIT) ) THEN
  !
  DO I = 1, KOPTION-INBDIGIT
    HCHAINE(I:I) = '0'
  END DO
  !     
ENDIF
!
DO I = 1,INBDIGIT
!
  INUM = INT(IRESTE*1./IDIVI)
  YC = CHAR(INUM+48)
  IF (KOPTION.GT.0) THEN
    J = KOPTION-INBDIGIT+I
  ELSEIF (KOPTION.EQ.0) THEN
    J = I
  ENDIF
  HCHAINE(J:J) = YC
  IRESTE = IRESTE - INUM*IDIVI
  IDIVI = IDIVI/10
!
END DO
!
IF (LHOOK) CALL DR_HOOK('TRANS_CHAINE',1,ZHOOK_HANDLE)
!---------------
END SUBROUTINE TRANS_CHAINE
