!     #########
      SUBROUTINE READ_ISBA_n (DTCO, I, U, &
                              HPROGRAM)
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
!
!
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_ISBA_n, ONLY : ISBA_t
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
TYPE(ISBA_t), INTENT(INOUT) :: I
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
REAL,DIMENSION(I%O%NPATCH) :: ZVLAIMIN
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
ALLOCATE(ZWORK(ILU,I%O%NPATCH))
!* soil temperatures
!
IF(I%O%LTEMP_ARP)THEN
  IWORK=I%O%NTEMPLAYER_ARP
ELSEIF(I%O%CISBA=='DIF')THEN
  IWORK=I%O%NGROUND_LAYER
ELSE
  IWORK=2 !Only 2 temperature layer in ISBA-FR
ENDIF
!
IF ( TRIM(CASSIM_ISBA)=="ENKF") THEN
  ALLOCATE(I%R%XRED_NOISE(ILU,I%O%NPATCH,NVAR))
  I%R%XRED_NOISE(:,:,:) = 0.
  ALLOCATE(ZCOFSWI(ILU))
  CALL COFSWI(I%P%XCLAY(:,1),ZCOFSWI)
ELSE
  ALLOCATE(I%R%XRED_NOISE(0,0,0))
  ALLOCATE(ZCOFSWI(0))
ENDIF
!
ALLOCATE(I%R%XTG(ILU,IWORK,I%O%NPATCH))
I%R%XTG(:,:,:)=XUNDEF
!
DO JL=1,IWORK
  WRITE(YLVL,'(I4)') JL
  YRECFM='TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  I%R%XTG(:,JL,:)=ZWORK(:,:)
END DO
!
! Perturb value if requested
IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
  !
  DO JL=1,IWORK
  ! read in control variable
    IF ( (TRIM(CVAR(NIVAR))=="TG1" .AND. JL==1) .OR. &
         (TRIM(CVAR(NIVAR))=="TG2" .AND. JL==2) ) THEN
      WHERE ( I%R%XTG(:,JL,:)/=XUNDEF )
        I%R%XTG(:,JL,:) = I%R%XTG(:,JL,:) + XTPRT(NIVAR)*I%R%XTG(:,JL,:)
      ENDWHERE
    ENDIF
  END DO
  !
ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 ) THEN
  !
  CALL MAKE_ENS_ENKF(IWORK,ILU,"TG ",ZCOFSWI,I%R%XTG,I%R%XRED_NOISE)
  !
ENDIF
!
!
!* soil liquid and ice water contents
!
ALLOCATE(I%R%XWG (ILU,I%O%NGROUND_LAYER,I%O%NPATCH))
ALLOCATE(I%R%XWGI(ILU,I%O%NGROUND_LAYER,I%O%NPATCH))
!
I%R%XWG (:,:,:)=XUNDEF
I%R%XWGI(:,:,:)=XUNDEF
!
DO JL=1,I%O%NGROUND_LAYER
  WRITE(YLVL,'(I4)') JL
  YRECFM='WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
   CALL READ_SURF(&
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
   I%R%XWG(:,JL,:)=ZWORK(:,:)
END DO
!
! Perturb value if requested
IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
   !
   DO JL=1,I%O%NGROUND_LAYER
    ! read in control variable
    IF ( (TRIM(CVAR(NIVAR))=="WG1" .AND. JL==1) .OR. & 
         (TRIM(CVAR(NIVAR))=="WG2" .AND. JL==2) ) THEN
      WHERE ( I%R%XWG(:,JL,:)/=XUNDEF )
        I%R%XWG(:,JL,:) = I%R%XWG(:,JL,:) + XTPRT(NIVAR)*I%R%XWG(:,JL,:)
      ENDWHERE
    ENDIF
   END DO
   !
ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 ) THEN
  !
  CALL MAKE_ENS_ENKF(IWORK,ILU,"WG ",ZCOFSWI,I%R%XWG,I%R%XRED_NOISE)
  !
ENDIF
!
IF(I%O%CISBA=='DIF')THEN
  IWORK=I%O%NGROUND_LAYER
ELSE
  IWORK=2 !Only 2 soil ice layer in ISBA-FR
ENDIF
!
DO JL=1,IWORK
  WRITE(YLVL,'(I4)') JL
  YRECFM='WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  I%R%XWGI(:,JL,:)=ZWORK(:,:)
END DO
!
!* water intercepted on leaves
!
ALLOCATE(I%R%XWR(ILU,I%O%NPATCH))
!
YRECFM = 'WR'
 CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%R%XWR(:,:),IRESP)
!
!* Leaf Area Index
!
IF (I%O%CPHOTO=='LAI' .OR. I%O%CPHOTO=='LST' .OR. I%O%CPHOTO=='NIT' .OR. I%O%CPHOTO=='NCB') THEN
  YRECFM = 'LAI'
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%M%T%XLAI(:,:),IRESP)
  IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
    !
    ! read in control variable
    IF ( TRIM(CVAR(NIVAR))=="LAI" ) THEN
      WHERE ( I%M%T%XLAI(:,:)/=XUNDEF ) 
        I%M%T%XLAI(:,:) = I%M%T%XLAI(:,:) + XTPRT(NIVAR)*I%M%T%XLAI(:,:)
      ENDWHERE
    ENDIF
    !
  ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 ) THEN
    !
    IF (I%O%NPATCH==12) THEN
      ZVLAIMIN = (/0.3,0.3,0.3,0.3,1.0,1.0,0.3,0.3,0.3,0.3,0.3,0.3/)
    ELSE
      ZVLAIMIN = (/0.3/)
    ENDIF
    !
    ALLOCATE(ZLAI(ILU,1,I%O%NPATCH))
    ZLAI(:,1,:) = I%M%T%XLAI(:,:)
    CALL MAKE_ENS_ENKF(1,ILU,"LAI",ZCOFSWI,ZLAI,I%R%XRED_NOISE)
    DO JP = 1,I%O%NPATCH
      I%M%T%XLAI(:,JP) = MAX(ZVLAIMIN(JP),ZLAI(:,1,JP))
    ENDDO
    DEALLOCATE(ZLAI)
    !    
  ENDIF  
END IF
!
!* snow mantel
!
 CALL READ_GR_SNOW(&
                   HPROGRAM,'VEG','     ',ILU,I%O%NPATCH,I%R%TSNOW  )
!
YRECFM='VERSION'
 CALL READ_SURF(&
                 HPROGRAM,YRECFM,IVERSION,IRESP)
!
YRECFM='BUG'
 CALL READ_SURF(&
                 HPROGRAM,YRECFM,IBUGFIX,IRESP)
!
IF(I%O%LGLACIER)THEN
  ALLOCATE(I%R%XICE_STO(ILU,I%O%NPATCH))
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=2) THEN
    YRECFM = 'ICE_STO'
    CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%R%XICE_STO(:,:),IRESP)
  ELSE
    I%R%XICE_STO(:,:) = 0.0
  ENDIF
ELSE
  ALLOCATE(I%R%XICE_STO(0,0))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.  MEB Prognostic or Semi-prognostic variables
!            -------------------------------------------
!
ISIZE_LMEB_PATCH=COUNT(I%O%LMEB_PATCH(:))
!
IF (ISIZE_LMEB_PATCH>0) THEN
!
!* water intercepted on litter

 ALLOCATE(I%R%XWRL(ILU,I%O%NPATCH))
 YRECFM = 'WRL'
 CALL READ_SURF(HPROGRAM,YRECFM,I%R%XWRL(:,:),IRESP)

 ALLOCATE(I%R%XWRLI(ILU,I%O%NPATCH))
 YRECFM = 'WRLI'
 CALL READ_SURF(HPROGRAM,YRECFM,I%R%XWRLI(:,:),IRESP)
!
!* snow intercepted on vegetation canopy leaves
!
  ALLOCATE(I%R%XWRVN(ILU,I%O%NPATCH))
  YRECFM = 'WRVN'
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%R%XWRVN(:,:),IRESP)
!
!* vegetation canopy temperature
!
  ALLOCATE(I%R%XTV(ILU,I%O%NPATCH))
  YRECFM = 'TV'
  CALL READ_SURF(HPROGRAM,YRECFM,I%R%XTV(:,:),IRESP)
!
!* litter temperature
!
  ALLOCATE(I%R%XTL(ILU,I%O%NPATCH))
  YRECFM = 'TL'
  CALL READ_SURF(HPROGRAM,YRECFM,I%R%XTL(:,:),IRESP)
!
!* vegetation canopy air temperature
!
  ALLOCATE(I%R%XTC(ILU,I%O%NPATCH))
  YRECFM = 'TC'
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%R%XTC(:,:),IRESP)
!
!* vegetation canopy air specific humidity
!
  ALLOCATE(I%R%XQC(ILU,I%O%NPATCH))
  YRECFM = 'QC'
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%R%XQC(:,:),IRESP)
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.  Semi-prognostic variables
!            -------------------------
!
ALLOCATE(I%R%XRESA(ILU,I%O%NPATCH))
ALLOCATE(I%R%XLE  (ILU,I%O%NPATCH))
IF (I%O%CPHOTO/='NON') THEN
  ALLOCATE(I%R%XANFM  (ILU,I%O%NPATCH))
  ALLOCATE(I%R%XAN    (ILU,I%O%NPATCH))
  ALLOCATE(I%R%XANDAY (ILU,I%O%NPATCH))
END IF
!
IF(I%O%CPHOTO/='NON') THEN
  ALLOCATE(I%R%XBIOMASS         (ILU,I%O%NNBIOMASS,I%O%NPATCH))
  ALLOCATE(I%R%XRESP_BIOMASS    (ILU,I%O%NNBIOMASS,I%O%NPATCH))
  I%R%XBIOMASS(:,:,:) = 0.
  I%R%XRESP_BIOMASS(:,:,:) = 0.    
END IF
!
!
!* aerodynamical resistance
!
YRECFM = 'RESA'
I%R%XRESA(:,:) = 100.
 CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%R%XRESA(:,:),IRESP)
!
!* patch averaged radiative temperature (K)
!
ALLOCATE(I%R%XTSRAD_NAT(ILU))
IF (IVERSION<6) THEN
  I%R%XTSRAD_NAT(:)=0.
  DO JP=1,I%O%NPATCH
    I%R%XTSRAD_NAT(:)=I%R%XTSRAD_NAT(:)+I%R%XTG(:,1,JP)
  ENDDO
  I%R%XTSRAD_NAT(:)=I%R%XTSRAD_NAT(:)/I%O%NPATCH
ELSE
  YRECFM='TSRAD_NAT'
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%R%XTSRAD_NAT(:),IRESP)
ENDIF
!
I%R%XLE(:,:) = XUNDEF
!
!*       5. ISBA-AGS variables
!
IF (I%O%CPHOTO/='NON') THEN
  YRECFM = 'AN'
  I%R%XAN(:,:) = 0.
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%R%XAN(:,:),IRESP)
  !
  YRECFM = 'ANDAY'
  I%R%XANDAY(:,:) = 0.
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%R%XANDAY(:,:),IRESP)
  !
  YRECFM = 'ANFM'
  I%R%XANFM(:,:) = XANFMINIT
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%R%XANFM(:,:),IRESP)
  !
  YRECFM = 'LE_AGS'
  I%R%XLE(:,:) = 0.
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,I%R%XLE(:,:),IRESP)
END IF
!
IF (I%O%CPHOTO=='LAI' .OR. I%O%CPHOTO=='LST') THEN
  !
  I%R%XBIOMASS(:,1,:) = I%M%T%XBSLAI(:,:) * I%M%T%XLAI(:,:)

ELSEIF (I%O%CPHOTO=='NIT'.OR.I%O%CPHOTO=='NCB') THEN
  !
  DO JNBIOMASS=1,I%O%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
      YRECFM='BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='BIOMASS'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF
    CALL READ_SURF(&
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
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
              DO JP = 1,I%O%NPATCH
                ZWORK(JI,JP) = ZWORK(JI,JP) + XADDINFL(JVAR)*RANDOM_NORMAL()
              ENDDO
            ENDDO
            EXIT
          ENDIF
        ENDDO
      ENDIF
      !      
    ENDIF     
    I%R%XBIOMASS(:,JNBIOMASS,:)=ZWORK
  END DO
!
  IWORK=0
  IF(I%O%CPHOTO=='NCB'.OR.IVERSION<8)IWORK=2
!
  DO JNBIOMASS=2,I%O%NNBIOMASS-IWORK
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. (IVERSION==7 .AND. IBUGFIX>=3)) THEN
      YRECFM='RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='RESP_BIOM'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF    
    CALL READ_SURF(&
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
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
              DO JP = 1,I%O%NPATCH
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
    I%R%XRESP_BIOMASS(:,JNBIOMASS,:)=ZWORK
  END DO
  !
ENDIF
!
DEALLOCATE(ZCOFSWI)
!
!*       6. Soil carbon
!
!
IF (I%O%CRESPSL=='CNT') THEN
  !
  ALLOCATE(I%R%XLITTER          (ILU,I%O%NNLITTER,I%O%NNLITTLEVS,I%O%NPATCH))
  ALLOCATE(I%R%XSOILCARB        (ILU,I%O%NNSOILCARB,I%O%NPATCH))
  ALLOCATE(I%R%XLIGNIN_STRUC    (ILU,I%O%NNLITTLEVS,I%O%NPATCH))
  !
  I%R%XLITTER(:,:,:,:) = 0.
  DO JNLITTER=1,I%O%NNLITTER
    DO JNLITTLEVS=1,I%O%NNLITTLEVS
      WRITE(YLVL,'(I1,A1,I1)') JNLITTER,'_',JNLITTLEVS
      YRECFM='LITTER'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      CALL READ_SURF(&
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
      I%R%XLITTER(:,JNLITTER,JNLITTLEVS,:)=ZWORK
    END DO
  END DO

  I%R%XSOILCARB(:,:,:) = 0.
  DO JNSOILCARB=1,I%O%NNSOILCARB
    WRITE(YLVL,'(I4)') JNSOILCARB
    YRECFM='SOILCARB'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL READ_SURF(&
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    I%R%XSOILCARB(:,JNSOILCARB,:)=ZWORK
  END DO
!
  I%R%XLIGNIN_STRUC(:,:,:) = 0.
  DO JNLITTLEVS=1,I%O%NNLITTLEVS
    WRITE(YLVL,'(I4)') JNLITTLEVS
    YRECFM='LIGNIN_STR'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL READ_SURF(&
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    I%R%XLIGNIN_STRUC(:,JNLITTLEVS,:)=ZWORK
  END DO
!
ENDIF

IF ( LASSIM ) THEN
  IF ( TRIM(CASSIM_ISBA) == "OI" ) THEN
    IF ( I%O%NPATCH /= 1 ) CALL ABOR1_SFX ('Reading of diagnostical values for'&
                       & //'assimilation at the moment only works for one patch for OI')          
    ! Diagnostic fields for assimilation
    IF ( .NOT. ALLOCATED(XAT2M_ISBA)) ALLOCATE(XAT2M_ISBA(ILU,1))
    XAT2M_ISBA=XUNDEF
    YRECFM='T2M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(&
                 HPROGRAM,YRECFM,XAT2M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAHU2M_ISBA)) ALLOCATE(XAHU2M_ISBA(ILU,1))
    XAHU2M_ISBA=XUNDEF
    YRECFM='HU2M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(&
                 HPROGRAM,YRECFM,XAHU2M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAZON10M_ISBA)) ALLOCATE(XAZON10M_ISBA(ILU,1))
    XAZON10M_ISBA=XUNDEF
    YRECFM='ZON10M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(&
                 HPROGRAM,YRECFM,XAZON10M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAMER10M_ISBA)) ALLOCATE(XAMER10M_ISBA(ILU,1))
    XAMER10M_ISBA=XUNDEF
    YRECFM='MER10M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(&
                 HPROGRAM,YRECFM,XAMER10M_ISBA(:,1),IRESP)
  ELSEIF ( NIFIC/=NVAR+2 ) THEN
    ! Diagnostic fields for EKF assimilation ("observations")
    DO IOBS = 1,NOBSTYPE
     SELECT CASE (TRIM(COBS(IOBS)))
       CASE("T2M")
         IF ( .NOT. ALLOCATED(XAT2M_ISBA)) ALLOCATE(XAT2M_ISBA(ILU,1))
         XAT2M_ISBA=XUNDEF
         YRECFM='T2M'
         CALL IO_BUFF(YRECFM,'R',GKNOWN)
         CALL READ_SURF(&
                 HPROGRAM,YRECFM,XAT2M_ISBA(:,1),IRESP)
       CASE("HU2M")
         IF ( .NOT. ALLOCATED(XAHU2M_ISBA)) ALLOCATE(XAHU2M_ISBA(ILU,1))
         XAHU2M_ISBA=XUNDEF
         YRECFM='HU2M'
         CALL IO_BUFF(YRECFM,'R',GKNOWN)
         CALL READ_SURF(&
                 HPROGRAM,YRECFM,XAHU2M_ISBA(:,1),IRESP)
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
          DO JP = 1,I%O%NPATCH
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
          DO JP = 1,I%O%NPATCH 
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
          DO JP = 1,I%O%NPATCH
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
