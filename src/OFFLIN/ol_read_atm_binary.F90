!     #########
SUBROUTINE OL_READ_ATM_BINARY(KFORC_STEP,                                 &
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
REAL, DIMENSION(:,:),INTENT(INOUT) :: PTA
REAL, DIMENSION(:,:),INTENT(INOUT) :: PQA
REAL, DIMENSION(:,:),INTENT(INOUT) :: PWIND
REAL, DIMENSION(:,:),INTENT(INOUT) :: PDIR_SW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PSCA_SW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PLW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PSNOW
REAL, DIMENSION(:,:),INTENT(INOUT) :: PRAIN
REAL, DIMENSION(:,:),INTENT(INOUT) :: PPS
REAL, DIMENSION(:,:),INTENT(INOUT) :: PCO2
REAL, DIMENSION(:,:,:),INTENT(OUT) :: PIMPWET
REAL, DIMENSION(:,:,:),INTENT(OUT) :: PIMPDRY
REAL, DIMENSION(:,:),INTENT(INOUT) :: PDIR
INTEGER,INTENT(IN)               :: KFORC_STEP
! local variables
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER  :: JIMP
! 
! read data
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_BINARY',0,ZHOOK_HANDLE)
 CALL READ_SURF_ATM('BINARY',PTA    (:,1:XCOUNT),KFORC_STEP,XCOUNT,22)
 CALL READ_SURF_ATM('BINARY',PQA    (:,1:XCOUNT),KFORC_STEP,XCOUNT,23)
 CALL READ_SURF_ATM('BINARY',PWIND  (:,1:XCOUNT),KFORC_STEP,XCOUNT,24)
 CALL READ_SURF_ATM('BINARY',PLW    (:,1:XCOUNT),KFORC_STEP,XCOUNT,25)
 CALL READ_SURF_ATM('BINARY',PDIR_SW(:,1:XCOUNT),KFORC_STEP,XCOUNT,26)
 CALL READ_SURF_ATM('BINARY',PSCA_SW(:,1:XCOUNT),KFORC_STEP,XCOUNT,27)
 CALL READ_SURF_ATM('BINARY',PRAIN  (:,1:XCOUNT),KFORC_STEP,XCOUNT,28)
 CALL READ_SURF_ATM('BINARY',PSNOW  (:,1:XCOUNT),KFORC_STEP,XCOUNT,29)
 CALL READ_SURF_ATM('BINARY',PPS    (:,1:XCOUNT),KFORC_STEP,XCOUNT,30)
 CALL READ_SURF_ATM('BINARY',PDIR   (:,1:XCOUNT),KFORC_STEP,XCOUNT,31)
 CALL READ_SURF_ATM('BINARY',PCO2   (:,1:XCOUNT),KFORC_STEP,XCOUNT,32)
 DO JIMP=1,NIMPUROF
  CALL READ_SURF_ATM('BINARY ',PIMPWET   (:,JIMP,1:XCOUNT),KFORC_STEP,XCOUNT,33)
  CALL READ_SURF_ATM('BINARY ',PIMPDRY   (:,JIMP,1:XCOUNT),KFORC_STEP,XCOUNT,34)
 ENDDO
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_BINARY',1,ZHOOK_HANDLE)
!
END SUBROUTINE OL_READ_ATM_BINARY
