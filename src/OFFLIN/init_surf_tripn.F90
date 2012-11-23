!#########################################################
SUBROUTINE INIT_SURF_TRIP_n(HPROGRAM, KI, KSW, ORESTART, KYEAR,     &
                              KMONTH, PDURATION, KTRIP_MONTH,        &
                              KTRIP_COUNT, PZENITH, PSW_BANDS,       &
                              PEMIS, PTSRAD, PDIR_ALB, PSCA_ALB      )  
!#########################################################
!
!!****  *INIT_SURF_TRIP_n* - routine to initialize the SURFACE-TRIP coupling
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
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
!!	B. Decharme   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2008
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS,       ONLY : XDAY
!
USE MODD_SURFEX_MPI, ONLY : NPROC
USE MODD_SURFEX_OMP, ONLY : NBLOCKTOT
!
USE MODI_GET_LUOUT
USE MODI_GET_CONF_ISBA_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_INI_CSTS
USE MODI_INIT_COUPLING_SURF_TRIP_n
USE MODI_INIT_TRIP_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM    ! program calling surf. schemes
LOGICAL,                          INTENT(IN)  :: ORESTART    
INTEGER,                          INTENT(IN)  :: KI         ! Surfex grid dimension
INTEGER,                          INTENT(IN)  :: KSW        ! Number of spectral bands
INTEGER,                          INTENT(IN)  :: KYEAR       ! current year (UTC)
INTEGER,                          INTENT(IN)  :: KMONTH      ! current month (UTC)
REAL,                             INTENT(IN)  :: PDURATION
INTEGER,                          INTENT(OUT) :: KTRIP_MONTH ! current output month (UTC)
INTEGER,                          INTENT(OUT) :: KTRIP_COUNT ! current TRIP counter
!
REAL,             DIMENSION(KI),  INTENT(IN)  :: PZENITH   ! solar zenithal angle
REAL,             DIMENSION(KI), INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),  INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSRAD    ! radiative temperature
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
CHARACTER(LEN=10) :: YGRID
!
LOGICAL :: LTRIP
LOGICAL :: LFLOOD
!
INTEGER :: IDIMTAB
INTEGER :: ILUOUT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_SURF_TRIP_N',0,ZHOOK_HANDLE)
CALL INI_CSTS
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
! * 1. Get ISBA configuration
!      
CALL GET_CONF_ISBA_n(LTRIP,LFLOOD,YGRID,IDIMTAB)
!
! * 2. Initialyse TRIP
!
IF(.NOT.LTRIP)THEN
!        
  IF(LFLOOD)THEN
    CALL ABOR1_SFX('INIT_SURF_TRIPN: LFLOOD=T BUT LTRIP=F')
  ENDIF
!
ELSE
!
 IF (NPROC>1) CALL ABOR1_SFX('INIT_SURF_TRIPN: TRIP CANNOT RUN WITH MORE THAN 1 MPI TASK') 
 IF (NBLOCKTOT>1) CALL ABOR1_SFX("INIT_SURF_TRIPN: TRIP CANNOT RUN WITH NUMEROUS OPENMP BLOCKS")
!
 KTRIP_MONTH=0
 IF(PDURATION/XDAY<=31.)THEN
    KTRIP_MONTH=KMONTH
 ELSEIF(PDURATION/XDAY>366.)THEN
   CALL ABOR1_SFX('Trip output time can not be superior to one year per run')
 ENDIF
!
  KTRIP_COUNT = 0
!  
  CALL INIT_TRIP_n(HPROGRAM,KYEAR,KTRIP_MONTH,ORESTART)
!
! * 3. Test and initialyse Surface-TRIP coupling
!
  CALL INIT_COUPLING_SURF_TRIP_n(HPROGRAM,KI,KSW,LFLOOD,YGRID,IDIMTAB,PZENITH,&
                                    PSW_BANDS,PEMIS,PTSRAD,PDIR_ALB,PSCA_ALB    )  
!          
ENDIF
IF (LHOOK) CALL DR_HOOK('INIT_SURF_TRIP_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_SURF_TRIP_n
