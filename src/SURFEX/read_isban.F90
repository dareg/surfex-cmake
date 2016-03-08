!     #########
      SUBROUTINE READ_ISBA_n (DTCO, IO, IR, PCLAY, PLAI, PBSLAI, U, HPROGRAM)
!     ##################################
!
!!****  *READ_ISBA_n* - routine to initialise ISBA variables
!!                         
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!!
!!      READ_SURF for general reading : 08/2003 (S.Malardel)
!!      B. Decharme  2008    : Floodplains
!!      B. Decharme  01/2009 : Optional Arpege deep soil temperature read
!!      A.L. Gibelin   03/09 : modifications for CENTURY model 
!!      A.L. Gibelin    04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin    06/2009 : Soil carbon variables for CNT option
!!      B. Decharme  09/2012 : suppress NWG_LAYER (parallelization problems)
!!      T. Aspelien  08/2013 : Read diagnostics for assimilation
!!      P. Samuelsson   10/2014 : MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_CO2V_PAR,       ONLY : XANFMINIT, XCONDCTMIN
!                          
USE MODD_ASSIM,          ONLY : LASSIM,CASSIM_ISBA,XAT2M_ISBA,XAHU2M_ISBA,&
                                XAZON10M_ISBA,XAMER10M_ISBA,NIFIC,NVAR, &
                                COBS,NOBSTYPE,CVAR,LPRT,XTPRT,NIVAR,CBIO, &
                                XADDINFL,NENS,XSIGMA,NIE
!                                
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_SNOW_PAR,       ONLY : XZ0SN
!
USE MODI_READ_SURF
!
USE MODI_READ_GR_SNOW
USE MODI_ABOR1_SFX
USE MODI_IO_BUFF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_TYPE_DIM_n
USE MODE_RANDOM
USE MODE_EKF
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
REAL, DIMENSION(:,:), INTENT(IN) :: PCLAY
REAL, DIMENSION(:,:), INTENT(INOUT) :: PLAI
REAL, DIMENSION(:,:), INTENT(INOUT) :: PBSLAI
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! calling program
!
!*       0.2   Declarations of local variables
!              -------------------------------
INTEGER           :: ILU          ! 1D physical dimension
!
INTEGER           :: IRESP          ! Error code after redding
!
CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
!
CHARACTER(LEN=4)  :: YLVL
!
REAL, DIMENSION(:,:,:),ALLOCATABLE :: ZLAI
REAL, DIMENSION(:,:),ALLOCATABLE  :: ZWORK      ! 2D array to write data in file
REAL, DIMENSION(:), ALLOCATABLE :: ZCOFSWI
!
REAL,DIMENSION(IO%NPATCH) :: ZVLAIMIN
REAL :: ZCOEF
!
INTEGER :: IWORK   ! Work integer
!
INTEGER :: JP, JL, JNBIOMASS, JNLITTER, JNSOILCARB, JNLITTLEVS  ! loop counter on layers
INTEGER :: JVAR, JI
!
INTEGER           :: IVERSION       ! surface version
INTEGER           :: IBUGFIX
INTEGER           :: IIVAR
INTEGER           :: IOBS
INTEGER           :: IBSUP
INTEGER           :: ISIZE_LMEB_PATCH
!
LOGICAL :: GKNOWN
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_NATURE'
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     'NATURE',ILU)
!
!
!*       2.     Prognostic fields:
!               -----------------
!
ALLOCATE(ZWORK(ILU,IO%NPATCH))
!* soil temperatures
!
IF(IO%LTEMP_ARP)THEN
  IWORK=IO%NTEMPLAYER_ARP
ELSEIF(IO%CISBA=='DIF')THEN
  IWORK=IO%NGROUND_LAYER
ELSE
  IWORK=2 !Only 2 temperature layer in ISBA-FR
ENDIF
!
IF ( TRIM(CASSIM_ISBA)=="ENKF") THEN
  ALLOCATE(IR%XRED_NOISE(ILU,IO%NPATCH,NVAR))
  IR%XRED_NOISE(:,:,:) = 0.
  ALLOCATE(ZCOFSWI(ILU))
  CALL COFSWI(PCLAY(:,1),ZCOFSWI)
ELSE
  ALLOCATE(IR%XRED_NOISE(0,0,0))
  ALLOCATE(ZCOFSWI(0))
ENDIF
!
ALLOCATE(IR%XTG(ILU,IWORK,IO%NPATCH))
IR%XTG(:,:,:)=XUNDEF
!
DO JL=1,IWORK
  WRITE(YLVL,'(I4)') JL
  YRECFM='TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  IR%XTG(:,JL,:)=ZWORK(:,:)
END DO
!
! Perturb value if requested
IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
  !
  DO JL=1,IWORK
  ! read in control variable
    IF ( (TRIM(CVAR(NIVAR))=="TG1" .AND. JL==1) .OR. &
         (TRIM(CVAR(NIVAR))=="TG2" .AND. JL==2) ) THEN
      WHERE ( IR%XTG(:,JL,:)/=XUNDEF )
        IR%XTG(:,JL,:) = IR%XTG(:,JL,:) + XTPRT(NIVAR)*IR%XTG(:,JL,:)
      ENDWHERE
    ENDIF
  END DO
  !
ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 ) THEN
  !
  CALL MAKE_ENS_ENKF(IWORK,ILU,"TG ",ZCOFSWI,IR%XTG,IR%XRED_NOISE)
  !
ENDIF
!
!
!* soil liquid and ice water contents
!
ALLOCATE(IR%XWG (ILU,IO%NGROUND_LAYER,IO%NPATCH))
ALLOCATE(IR%XWGI(ILU,IO%NGROUND_LAYER,IO%NPATCH))
!
IR%XWG (:,:,:)=XUNDEF
IR%XWGI(:,:,:)=XUNDEF
!
DO JL=1,IO%NGROUND_LAYER
  WRITE(YLVL,'(I4)') JL
  YRECFM='WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
   CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
   IR%XWG(:,JL,:)=ZWORK(:,:)
END DO
!
! Perturb value if requested
IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
   !
   DO JL=1,IO%NGROUND_LAYER
    ! read in control variable
    IF ( (TRIM(CVAR(NIVAR))=="WG1" .AND. JL==1) .OR. & 
         (TRIM(CVAR(NIVAR))=="WG2" .AND. JL==2) ) THEN
      WHERE ( IR%XWG(:,JL,:)/=XUNDEF )
        IR%XWG(:,JL,:) = IR%XWG(:,JL,:) + XTPRT(NIVAR)*IR%XWG(:,JL,:)
      ENDWHERE
    ENDIF
   END DO
   !
ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 ) THEN
  !
  CALL MAKE_ENS_ENKF(IWORK,ILU,"WG ",ZCOFSWI,IR%XWG,IR%XRED_NOISE)
  !
ENDIF
!
IF(IO%CISBA=='DIF')THEN
  IWORK=IO%NGROUND_LAYER
ELSE
  IWORK=2 !Only 2 soil ice layer in ISBA-FR
ENDIF
!
DO JL=1,IWORK
  WRITE(YLVL,'(I4)') JL
  YRECFM='WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  IR%XWGI(:,JL,:)=ZWORK(:,:)
END DO
!
!* water intercepted on leaves
!
ALLOCATE(IR%XWR(ILU,IO%NPATCH))
!
YRECFM = 'WR'
 CALL READ_SURF(HPROGRAM,YRECFM,IR%XWR(:,:),IRESP)
!
!* Leaf Area Index
!
IF (IO%CPHOTO=='LAI' .OR. IO%CPHOTO=='LST' .OR. IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
  YRECFM = 'LAI'
  CALL READ_SURF(HPROGRAM,YRECFM,PLAI(:,:),IRESP)
  IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
    !
    ! read in control variable
    IF ( TRIM(CVAR(NIVAR))=="LAI" ) THEN
      WHERE ( PLAI(:,:)/=XUNDEF ) 
        PLAI(:,:) = PLAI(:,:) + XTPRT(NIVAR)*PLAI(:,:)
      ENDWHERE
    ENDIF
    !
  ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 ) THEN
    !
    IF (IO%NPATCH==12) THEN
      ZVLAIMIN = (/0.3,0.3,0.3,0.3,1.0,1.0,0.3,0.3,0.3,0.3,0.3,0.3/)
    ELSE
      ZVLAIMIN = (/0.3/)
    ENDIF
    !
    ALLOCATE(ZLAI(ILU,1,IO%NPATCH))
    ZLAI(:,1,:) = PLAI(:,:)
    CALL MAKE_ENS_ENKF(1,ILU,"LAI",ZCOFSWI,ZLAI,IR%XRED_NOISE)
    DO JP = 1,IO%NPATCH
      PLAI(:,JP) = MAX(ZVLAIMIN(JP),ZLAI(:,1,JP))
    ENDDO
    DEALLOCATE(ZLAI)
    !    
  ENDIF  
END IF
!
!* snow mantel
!
 CALL READ_GR_SNOW(HPROGRAM,'VEG','     ',ILU,IO%NPATCH,IR%TSNOW  )
!
YRECFM='VERSION'
 CALL READ_SURF(HPROGRAM,YRECFM,IVERSION,IRESP)
!
YRECFM='BUG'
 CALL READ_SURF(HPROGRAM,YRECFM,IBUGFIX,IRESP)
!
IF(IO%LGLACIER)THEN
  ALLOCATE(IR%XICE_STO(ILU,IO%NPATCH))
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=2) THEN
    YRECFM = 'ICE_STO'
    CALL READ_SURF(HPROGRAM,YRECFM,IR%XICE_STO(:,:),IRESP)
  ELSE
    IR%XICE_STO(:,:) = 0.0
  ENDIF
ELSE
  ALLOCATE(IR%XICE_STO(0,0))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.  MEB Prognostic or Semi-prognostic variables
!            -------------------------------------------
!
ISIZE_LMEB_PATCH=COUNT(IO%LMEB_PATCH(:))
!
IF (ISIZE_LMEB_PATCH>0) THEN
!
!* water intercepted on litter

 ALLOCATE(IR%XWRL(ILU,IO%NPATCH))
 YRECFM = 'WRL'
 CALL READ_SURF(HPROGRAM,YRECFM,IR%XWRL(:,:),IRESP)

 ALLOCATE(IR%XWRLI(ILU,IO%NPATCH))
 YRECFM = 'WRLI'
 CALL READ_SURF(HPROGRAM,YRECFM,IR%XWRLI(:,:),IRESP)
!
!* snow intercepted on vegetation canopy leaves
!
  ALLOCATE(IR%XWRVN(ILU,IO%NPATCH))
  YRECFM = 'WRVN'
  CALL READ_SURF(HPROGRAM,YRECFM,IR%XWRVN(:,:),IRESP)
!
!* vegetation canopy temperature
!
  ALLOCATE(IR%XTV(ILU,IO%NPATCH))
  YRECFM = 'TV'
  CALL READ_SURF(HPROGRAM,YRECFM,IR%XTV(:,:),IRESP)
!
!* litter temperature
!
  ALLOCATE(IR%XTL(ILU,IO%NPATCH))
  YRECFM = 'TL'
  CALL READ_SURF(HPROGRAM,YRECFM,IR%XTL(:,:),IRESP)
!
!* vegetation canopy air temperature
!
  ALLOCATE(IR%XTC(ILU,IO%NPATCH))
  YRECFM = 'TC'
  CALL READ_SURF(HPROGRAM,YRECFM,IR%XTC(:,:),IRESP)
!
!* vegetation canopy air specific humidity
!
  ALLOCATE(IR%XQC(ILU,IO%NPATCH))
  YRECFM = 'QC'
  CALL READ_SURF(HPROGRAM,YRECFM,IR%XQC(:,:),IRESP)
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.  Semi-prognostic variables
!            -------------------------
!
ALLOCATE(IR%XRESA(ILU,IO%NPATCH))
ALLOCATE(IR%XLE  (ILU,IO%NPATCH))
IF (IO%CPHOTO/='NON') THEN
  ALLOCATE(IR%XANFM  (ILU,IO%NPATCH))
  ALLOCATE(IR%XAN    (ILU,IO%NPATCH))
  ALLOCATE(IR%XANDAY (ILU,IO%NPATCH))
END IF
!
IF(IO%CPHOTO/='NON') THEN
  ALLOCATE(IR%XBIOMASS         (ILU,IO%NNBIOMASS,IO%NPATCH))
  ALLOCATE(IR%XRESP_BIOMASS    (ILU,IO%NNBIOMASS,IO%NPATCH))
  IR%XBIOMASS(:,:,:) = 0.
  IR%XRESP_BIOMASS(:,:,:) = 0.    
END IF
!
!
!* aerodynamical resistance
!
YRECFM = 'RESA'
IR%XRESA(:,:) = 100.
 CALL READ_SURF(HPROGRAM,YRECFM,IR%XRESA(:,:),IRESP)
!
!* patch averaged radiative temperature (K)
!
ALLOCATE(IR%XTSRAD_NAT(ILU))
IF (IVERSION<6) THEN
  IR%XTSRAD_NAT(:)=0.
  DO JP=1,IO%NPATCH
    IR%XTSRAD_NAT(:)=IR%XTSRAD_NAT(:)+IR%XTG(:,1,JP)
  ENDDO
  IR%XTSRAD_NAT(:)=IR%XTSRAD_NAT(:)/IO%NPATCH
ELSE
  YRECFM='TSRAD_NAT'
  CALL READ_SURF(HPROGRAM,YRECFM,IR%XTSRAD_NAT(:),IRESP)
ENDIF
!
IR%XLE(:,:) = XUNDEF
!
!*       5. ISBA-AGS variables
!
IF (IO%CPHOTO/='NON') THEN
  YRECFM = 'AN'
  IR%XAN(:,:) = 0.
  CALL READ_SURF(HPROGRAM,YRECFM,IR%XAN(:,:),IRESP)
  !
  YRECFM = 'ANDAY'
  IR%XANDAY(:,:) = 0.
  CALL READ_SURF(HPROGRAM,YRECFM,IR%XANDAY(:,:),IRESP)
  !
  YRECFM = 'ANFM'
  IR%XANFM(:,:) = XANFMINIT
  CALL READ_SURF(HPROGRAM,YRECFM,IR%XANFM(:,:),IRESP)
  !
  YRECFM = 'LE_AGS'
  IR%XLE(:,:) = 0.
  CALL READ_SURF(HPROGRAM,YRECFM,IR%XLE(:,:),IRESP)
END IF
!
IF (IO%CPHOTO=='LAI' .OR. IO%CPHOTO=='LST') THEN
  !
  IR%XBIOMASS(:,1,:) = PBSLAI(:,:) * PLAI(:,:)

ELSEIF (IO%CPHOTO=='NIT'.OR.IO%CPHOTO=='NCB') THEN
  !
  DO JNBIOMASS=1,IO%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
      YRECFM='BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='BIOMASS'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF
    CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
      ! read in control variable
      IF ( TRIM(CVAR(NIVAR)) == "LAI" .AND. TRIM(CBIO)==TRIM(YRECFM) ) THEN
        WHERE ( ZWORK(:,:)/=XUNDEF ) 
          ZWORK(:,:) = ZWORK(:,:) + XTPRT(NIVAR)*ZWORK(:,:)
        ENDWHERE
      ENDIF
    ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 .AND. .NOT.LASSIM ) THEN
      !
      IF ( TRIM(CBIO)==TRIM(YRECFM) ) THEN
        DO JVAR = 1,NVAR
          IF (TRIM(CVAR(JVAR)) == "LAI") THEN
            DO JI = 1,ILU
              DO JP = 1,IO%NPATCH
                ZWORK(JI,JP) = ZWORK(JI,JP) + XADDINFL(JVAR)*RANDOM_NORMAL()
              ENDDO
            ENDDO
            EXIT
          ENDIF
        ENDDO
      ENDIF
      !      
    ENDIF     
    IR%XBIOMASS(:,JNBIOMASS,:)=ZWORK
  END DO
!
  IWORK=0
  IF(IO%CPHOTO=='NCB'.OR.IVERSION<8)IWORK=2
!
  DO JNBIOMASS=2,IO%NNBIOMASS-IWORK
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. (IVERSION==7 .AND. IBUGFIX>=3)) THEN
      YRECFM='RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='RESP_BIOM'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF    
    CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
      ! read in control variable
      IF ( TRIM(CVAR(NIVAR)) == "LAI" .AND. TRIM(CBIO)==TRIM(YRECFM) ) THEN
        WHERE ( ZWORK(:,:)/=XUNDEF ) 
          ZWORK(:,:) = ZWORK(:,:) + XTPRT(NIVAR)*ZWORK(:,:)
        ENDWHERE
    ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 .AND. .NOT.LASSIM ) THEN
      !
      IF ( TRIM(CBIO)==TRIM(YRECFM) ) THEN
        DO JVAR = 1,NVAR
          IF (TRIM(CVAR(JVAR)) == "LAI") THEN
            DO JI = 1,ILU
              DO JP = 1,IO%NPATCH
                ZWORK(JI,JP) = ZWORK(JI,JP) + XADDINFL(JVAR)*RANDOM_NORMAL()
              ENDDO
            ENDDO
            EXIT
          ENDIF
        ENDDO
      ENDIF
      !  
      ENDIF
    ENDIF      
    IR%XRESP_BIOMASS(:,JNBIOMASS,:)=ZWORK
  END DO
  !
ENDIF
!
DEALLOCATE(ZCOFSWI)
!
!*       6. Soil carbon
!
!
IF (IO%CRESPSL=='CNT') THEN
  !
  ALLOCATE(IR%XLITTER          (ILU,IO%NNLITTER,IO%NNLITTLEVS,IO%NPATCH))
  ALLOCATE(IR%XSOILCARB        (ILU,IO%NNSOILCARB,IO%NPATCH))
  ALLOCATE(IR%XLIGNIN_STRUC    (ILU,IO%NNLITTLEVS,IO%NPATCH))
  !
  IR%XLITTER(:,:,:,:) = 0.
  DO JNLITTER=1,IO%NNLITTER
    DO JNLITTLEVS=1,IO%NNLITTLEVS
      WRITE(YLVL,'(I1,A1,I1)') JNLITTER,'_',JNLITTLEVS
      YRECFM='LITTER'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
      IR%XLITTER(:,JNLITTER,JNLITTLEVS,:)=ZWORK
    END DO
  END DO

  IR%XSOILCARB(:,:,:) = 0.
  DO JNSOILCARB=1,IO%NNSOILCARB
    WRITE(YLVL,'(I4)') JNSOILCARB
    YRECFM='SOILCARB'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    IR%XSOILCARB(:,JNSOILCARB,:)=ZWORK
  END DO
!
  IR%XLIGNIN_STRUC(:,:,:) = 0.
  DO JNLITTLEVS=1,IO%NNLITTLEVS
    WRITE(YLVL,'(I4)') JNLITTLEVS
    YRECFM='LIGNIN_STR'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    IR%XLIGNIN_STRUC(:,JNLITTLEVS,:)=ZWORK
  END DO
!
ENDIF

IF ( LASSIM ) THEN
  IF ( TRIM(CASSIM_ISBA) == "OI" ) THEN
    IF ( IO%NPATCH /= 1 ) CALL ABOR1_SFX ('Reading of diagnostical values for'&
                       & //'assimilation at the moment only works for one patch for OI')          
    ! Diagnostic fields for assimilation
    IF ( .NOT. ALLOCATED(XAT2M_ISBA)) ALLOCATE(XAT2M_ISBA(ILU,1))
    XAT2M_ISBA=XUNDEF
    YRECFM='T2M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(HPROGRAM,YRECFM,XAT2M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAHU2M_ISBA)) ALLOCATE(XAHU2M_ISBA(ILU,1))
    XAHU2M_ISBA=XUNDEF
    YRECFM='HU2M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(HPROGRAM,YRECFM,XAHU2M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAZON10M_ISBA)) ALLOCATE(XAZON10M_ISBA(ILU,1))
    XAZON10M_ISBA=XUNDEF
    YRECFM='ZON10M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(HPROGRAM,YRECFM,XAZON10M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAMER10M_ISBA)) ALLOCATE(XAMER10M_ISBA(ILU,1))
    XAMER10M_ISBA=XUNDEF
    YRECFM='MER10M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(HPROGRAM,YRECFM,XAMER10M_ISBA(:,1),IRESP)
  ELSEIF ( NIFIC/=NVAR+2 ) THEN
    ! Diagnostic fields for EKF assimilation ("observations")
    DO IOBS = 1,NOBSTYPE
     SELECT CASE (TRIM(COBS(IOBS)))
       CASE("T2M")
         IF ( .NOT. ALLOCATED(XAT2M_ISBA)) ALLOCATE(XAT2M_ISBA(ILU,1))
         XAT2M_ISBA=XUNDEF
         YRECFM='T2M'
         CALL IO_BUFF(YRECFM,'R',GKNOWN)
         CALL READ_SURF(HPROGRAM,YRECFM,XAT2M_ISBA(:,1),IRESP)
       CASE("HU2M")
         IF ( .NOT. ALLOCATED(XAHU2M_ISBA)) ALLOCATE(XAHU2M_ISBA(ILU,1))
         XAHU2M_ISBA=XUNDEF
         YRECFM='HU2M'
         CALL IO_BUFF(YRECFM,'R',GKNOWN)
         CALL READ_SURF(HPROGRAM,YRECFM,XAHU2M_ISBA(:,1),IRESP)
       CASE("WG1")
         ! This is already read above
       CASE("LAI")
         ! This is already read above   
       CASE("SWE")
         ! This is handled independently 
       CASE DEFAULT
         CALL ABOR1_SFX("Mapping of "//TRIM(COBS(IOBS))//" is not defined in READ_ISBA_n!")
     END SELECT
    ENDDO
  ENDIF
ENDIF
!
DEALLOCATE(ZWORK)
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_N',1,ZHOOK_HANDLE)
!
CONTAINS
!
SUBROUTINE MAKE_ENS_ENKF(KWORK,KLU,HREC,PCOFSWI,PVAR,PRED_NOISE)
!
USE MODD_ASSIM, ONLY : LENS_GEN, XADDTIMECORR, XADDINFL, XASSIM_WINH
!
USE MODI_ADD_NOISE
USE MODE_RANDOM
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KWORK
INTEGER, INTENT(IN) :: KLU
 CHARACTER(LEN=3), INTENT(IN) :: HREC
REAL, DIMENSION(:), INTENT(IN) :: PCOFSWI
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PVAR
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRED_NOISE
!
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=4) :: YLVL
 CHARACTER(LEN=3) :: YVAR
REAL :: ZWHITE_NOISE, ZVAR0
INTEGER :: JL, JI, JP, IVAR
LOGICAL :: GPASS
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_N:MAKE_ENS_ENKF',0,ZHOOK_HANDLE)
!
!
DO JL=1,KWORK
  !
  IF (KWORK>1) THEN
    WRITE(YLVL,'(I4)') JL
    YRECFM = TRIM(HREC)//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ELSE
    YRECFM = TRIM(HREC)
  ENDIF
  !
  IVAR = 0
  DO JVAR = 1,NVAR
    GPASS = ( TRIM(CVAR(JVAR))==TRIM(YRECFM) )
    IF (GPASS) THEN
      IVAR = JVAR
      EXIT
    ENDIF
  ENDDO
  !
  IF ( GPASS ) THEN
    !
    IF (XADDINFL(IVAR)>0.) THEN
      !
      IF (LASSIM) THEN
        !
        WRITE(YVAR,'(I3)') IVAR
        YRECFM='RED_NOISE'//ADJUSTL(YVAR(:LEN_TRIM(YVAR)))
        CALL READ_SURF(HPROGRAM,YRECFM,PRED_NOISE(:,:,IVAR),IRESP)
        !
      ELSEIF (.NOT.LENS_GEN .AND. XADDTIMECORR(IVAR)>0. ) THEN
        !
        WRITE(YVAR,'(I3)') IVAR
        YRECFM='RED_NOISE'//ADJUSTL(YVAR(:LEN_TRIM(YVAR)))
        CALL READ_SURF(HPROGRAM,YRECFM,PRED_NOISE(:,:,IVAR),IRESP)
        !
        DO JI = 1,KLU
          DO JP = 1,IO%NPATCH
            ZWHITE_NOISE = XADDINFL(IVAR)*PCOFSWI(JI)*RANDOM_NORMAL()
            CALL ADD_NOISE(XADDTIMECORR(IVAR),XASSIM_WINH,ZWHITE_NOISE,PRED_NOISE(JI,JP,IVAR))
         ENDDO
         ENDDO
        !
        ZCOEF = XASSIM_WINH/24.
        !
      ELSE
        !
        DO JI = 1,ILU
          DO JP = 1,IO%NPATCH 
            PRED_NOISE(JI,JP,IVAR) = XADDINFL(IVAR)*PCOFSWI(JI)*RANDOM_NORMAL()
          ENDDO
        ENDDO
        !
        ZCOEF = 1. 
        !
      ENDIF
      !
      IF (.NOT.LASSIM) THEN
        !
        DO JI = 1,ILU
          DO JP = 1,IO%NPATCH
            IF ( PVAR(JI,JL,JP)/=XUNDEF ) THEN
              !
              ZVAR0 = PVAR(JI,JL,JP)
              !
              PVAR(JI,JL,JP) = PVAR(JI,JL,JP) + ZCOEF * PRED_NOISE(JI,JP,IVAR)
              !
              IF (PVAR(JI,JL,JP) < 0.) THEN
                IF (LENS_GEN) THEN
                  PVAR(JI,JL,JP) = ABS(PVAR(JI,JL,JP))
                ELSE
                  PVAR(JI,JL,JP) = ZVAR0
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        !
      ENDIF
      !
    ENDIF
    !
  ENDIF
  !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_N:MAKE_ENS_ENKF',1,ZHOOK_HANDLE)
!
END SUBROUTINE MAKE_ENS_ENKF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_ISBA_n
