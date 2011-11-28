!     #########
SUBROUTINE OL_READ_ATM_BINARY(HSURF_FILETYPE, KFORC_STEP,                 &
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
REAL, DIMENSION(:,:),INTENT(INOUT) :: PDIR
INTEGER,INTENT(IN)               :: KFORC_STEP
CHARACTER(LEN=6)    ,INTENT(IN)  :: HSURF_FILETYPE

! local variables
INTEGER                          :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!

   
! read data

IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_BINARY',0,ZHOOK_HANDLE)
CALL READ_SURF_ATM('BINARY','Forc_TA.bin    ',PTA    (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,22)
CALL READ_SURF_ATM('BINARY','Forc_QA.bin    ',PQA    (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,23)
CALL READ_SURF_ATM('BINARY','Forc_WIND.bin  ',PWIND  (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,24)
CALL READ_SURF_ATM('BINARY','Forc_LW.bin    ',PLW    (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,25)
CALL READ_SURF_ATM('BINARY','Forc_DIR_SW.bin',PDIR_SW(1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,26)
CALL READ_SURF_ATM('BINARY','Forc_SCA_SW.bin',PSCA_SW(1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,27)
CALL READ_SURF_ATM('BINARY','Forc_RAIN.bin  ',PRAIN  (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,28)
CALL READ_SURF_ATM('BINARY','Forc_SNOW.bin  ',PSNOW  (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,29)
CALL READ_SURF_ATM('BINARY','Forc_PS.bin    ',PPS    (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,30)
CALL READ_SURF_ATM('BINARY','Forc_DIR.bin   ',PDIR   (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,31)
CALL READ_SURF_ATM('BINARY','Forc_CO2.bin   ',PCO2   (1:XCOUNT,:),KFORC_STEP,XCOUNT,IRET,32)
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_BINARY',1,ZHOOK_HANDLE)


END SUBROUTINE OL_READ_ATM_BINARY
