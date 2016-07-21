!     #########
SUBROUTINE OL_READ_ATM_ASCII (KFORC_STEP,                                 &
                              PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,    &
                              PRAIN,PPS,PCO2,PIMPWET,PIMPDRY,PDIR                         )  
!**************************************************************************
!
!!    PURPOSE
!!    -------
!         Read in the ascii file the atmospheric forcing for the actual time
!         step KFORC_STEP, and for the next one.
!         The two time step are needed for the time interpolation of the
!         forcing.
!         If the end of the file  is reached, set the two step to the last
!         values.
!         Return undef value if the variable is not present
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      A. Lemonsu  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     03/2008       
!
USE MODD_IO_SURF_OL, ONLY : XCOUNT
USE MODI_READ_SURF_ATM
USE MODN_IO_OFFLINE,  ONLY : NIMPUROF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
! global variables
REAL, DIMENSION(:,:),INTENT(OUT) :: PTA
REAL, DIMENSION(:,:),INTENT(OUT) :: PQA
REAL, DIMENSION(:,:),INTENT(OUT) :: PWIND
REAL, DIMENSION(:,:),INTENT(OUT) :: PDIR_SW
REAL, DIMENSION(:,:),INTENT(OUT) :: PSCA_SW
REAL, DIMENSION(:,:),INTENT(OUT) :: PLW
REAL, DIMENSION(:,:),INTENT(OUT) :: PSNOW
REAL, DIMENSION(:,:),INTENT(OUT) :: PRAIN
REAL, DIMENSION(:,:),INTENT(OUT) :: PPS
REAL, DIMENSION(:,:),INTENT(OUT) :: PCO2
REAL, DIMENSION(:,:,:),INTENT(OUT) :: PIMPWET
REAL, DIMENSION(:,:,:),INTENT(OUT) :: PIMPDRY
REAL, DIMENSION(:,:),INTENT(OUT) :: PDIR
INTEGER,INTENT(IN)               :: KFORC_STEP
! local variables
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER  :: JIMP
!
! read data
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_ASCII',0,ZHOOK_HANDLE)
 CALL READ_SURF_ATM('ASCII ',PTA    (:,1:XCOUNT),KFORC_STEP,XCOUNT,22)
 CALL READ_SURF_ATM('ASCII ',PQA    (:,1:XCOUNT),KFORC_STEP,XCOUNT,23)
 CALL READ_SURF_ATM('ASCII ',PWIND  (:,1:XCOUNT),KFORC_STEP,XCOUNT,24)
 CALL READ_SURF_ATM('ASCII ',PLW    (:,1:XCOUNT),KFORC_STEP,XCOUNT,25)
 CALL READ_SURF_ATM('ASCII ',PDIR_SW(:,1:XCOUNT),KFORC_STEP,XCOUNT,26)
 CALL READ_SURF_ATM('ASCII ',PSCA_SW(:,1:XCOUNT),KFORC_STEP,XCOUNT,27)
 CALL READ_SURF_ATM('ASCII ',PRAIN  (:,1:XCOUNT),KFORC_STEP,XCOUNT,28)
 CALL READ_SURF_ATM('ASCII ',PSNOW  (:,1:XCOUNT),KFORC_STEP,XCOUNT,29)
 CALL READ_SURF_ATM('ASCII ',PPS    (:,1:XCOUNT),KFORC_STEP,XCOUNT,30)
 CALL READ_SURF_ATM('ASCII ',PDIR   (:,1:XCOUNT),KFORC_STEP,XCOUNT,31)
 CALL READ_SURF_ATM('ASCII ',PCO2   (:,1:XCOUNT),KFORC_STEP,XCOUNT,32)
 DO JIMP=1,NIMPUROF
  CALL READ_SURF_ATM('ASCII ',PIMPWET   (:,JIMP,1:XCOUNT),KFORC_STEP,XCOUNT,33)
  CALL READ_SURF_ATM('ASCII ',PIMPDRY   (:,JIMP,1:XCOUNT),KFORC_STEP,XCOUNT,34)
 ENDDO
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_ASCII',1,ZHOOK_HANDLE)
!
END SUBROUTINE OL_READ_ATM_ASCII
