!     #########
SUBROUTINE OL_READ_ATM_ASCII (HSURF_FILETYPE, KFORC_STEP,                 &
                                PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,    &
                                PRAIN,PPS,PCO2,PDIR                         )  
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
!!	A. Lemonsu  *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     03/2008

!          
!
USE MODD_IO_SURF_OL, ONLY : XCOUNT
USE MODI_READ_SURF_ATM
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
REAL, DIMENSION(:,:),INTENT(OUT) :: PDIR
INTEGER,INTENT(IN)               :: KFORC_STEP
CHARACTER(LEN=6)    ,INTENT(IN)  :: HSURF_FILETYPE

! local variables
INTEGER                          :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!

   
! read data
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_ASCII',0,ZHOOK_HANDLE)
CALL READ_SURF_ATM('ASCII ','Forc_TA.txt    ',PTA    (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,22)
CALL READ_SURF_ATM('ASCII ','Forc_QA.txt    ',PQA    (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,23)
CALL READ_SURF_ATM('ASCII ','Forc_WIND.txt  ',PWIND  (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,24)
CALL READ_SURF_ATM('ASCII ','Forc_LW.txt    ',PLW    (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,25)
CALL READ_SURF_ATM('ASCII ','Forc_DIR_SW.txt',PDIR_SW(1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,26)
CALL READ_SURF_ATM('ASCII ','Forc_SCA_SW.txt',PSCA_SW(1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,27)
CALL READ_SURF_ATM('ASCII ','Forc_RAIN.txt  ',PRAIN  (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,28)
CALL READ_SURF_ATM('ASCII ','Forc_SNOW.txt  ',PSNOW  (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,29)
CALL READ_SURF_ATM('ASCII ','Forc_PS.txt    ',PPS    (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,30)
CALL READ_SURF_ATM('ASCII ','Forc_DIR.txt   ',PDIR   (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,31)
CALL READ_SURF_ATM('ASCII ','Forc_CO2.txt   ',PCO2   (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,32)
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_ASCII',1,ZHOOK_HANDLE)


END SUBROUTINE OL_READ_ATM_ASCII
