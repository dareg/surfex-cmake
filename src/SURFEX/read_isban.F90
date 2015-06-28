!     #########
      SUBROUTINE READ_ISBA_n(HPROGRAM)
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
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_CO2V_PAR,       ONLY : XANFMINIT, XCONDCTMIN
USE MODD_ISBA_n, ONLY : I => ISBA
!                          
USE MODD_ASSIM,          ONLY : LASSIM,CASSIM_ISBA,XAT2M_ISBA,XAHU2M_ISBA,&
                              & XAZON10M_ISBA,XAMER10M_ISBA,NIPERT,NVAR, &
                              & COBS,NOBSTYPE,CVAR,LPRT,XTPRT,NIVAR,CBIO
!                                
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_SNOW_PAR,       ONLY : XZ0SN
!
USE MODI_READ_SURF
!
USE MODI_READ_GR_SNOW
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_TYPE_DIM_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
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
REAL, DIMENSION(:,:),ALLOCATABLE  :: ZWORK      ! 2D array to write data in file
!
INTEGER :: IWORK   ! Work integer
!
INTEGER :: JP, JL, JNBIOMASS, JNLITTER, JNSOILCARB, JNLITTLEVS  ! loop counter on layers
!
INTEGER           :: IVERSION       ! surface version
INTEGER           :: IBUGFIX
INTEGER           :: IIVAR
INTEGER           :: IOBS
INTEGER           :: IBSUP
INTEGER           :: ISIZE_LMEB_PATCH
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
ALLOCATE(ZWORK(ILU,I%NPATCH))
!* soil temperatures
!
IF(I%LTEMP_ARP)THEN
  IWORK=I%NTEMPLAYER_ARP
ELSEIF(I%CISBA=='DIF')THEN
  IWORK=I%NGROUND_LAYER
ELSE
  IWORK=2 !Only 2 temperature layer in ISBA-FR
ENDIF
!
ALLOCATE(I%XTG(ILU,IWORK,I%NPATCH))
I%XTG(:,:,:)=XUNDEF
!
DO JL=1,IWORK
  WRITE(YLVL,'(I4)') JL
  YRECFM='TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  I%XTG(:,JL,:)=ZWORK(:,:)
END DO
!
! Perturb value if requested
IF ( LPRT ) THEN
  DO JL=1,IWORK
  ! read in control variable
    IF ( (TRIM(CVAR(NIVAR))=="TG1" .AND. JL==1) .OR. &
         (TRIM(CVAR(NIVAR))=="TG2" .AND. JL==2) ) THEN
      WHERE ( I%XTG(:,JL,:)/=XUNDEF )
        I%XTG(:,JL,:) = I%XTG(:,JL,:) + XTPRT(NIVAR)*I%XTG(:,JL,:)
      ENDWHERE
    ENDIF
  END DO
ENDIF
!
!
!* soil liquid and ice water contents
!
ALLOCATE(I%XWG (ILU,I%NGROUND_LAYER,I%NPATCH))
ALLOCATE(I%XWGI(ILU,I%NGROUND_LAYER,I%NPATCH))
!
I%XWG (:,:,:)=XUNDEF
I%XWGI(:,:,:)=XUNDEF
!
DO JL=1,I%NGROUND_LAYER
  WRITE(YLVL,'(I4)') JL
  YRECFM='WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
   CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
   I%XWG(:,JL,:)=ZWORK(:,:)
END DO
!
! Perturb value if requested
IF ( LPRT ) THEN
   DO JL=1,I%NGROUND_LAYER
    ! read in control variable
    IF ( (TRIM(CVAR(NIVAR))=="WG1" .AND. JL==1) .OR. & 
         (TRIM(CVAR(NIVAR))=="WG2" .AND. JL==2) ) THEN
      WHERE ( I%XWG(:,JL,:)/=XUNDEF )
        I%XWG(:,JL,:) = I%XWG(:,JL,:) + XTPRT(NIVAR)*I%XWG(:,JL,:)
      ENDWHERE
    ENDIF
   END DO
ENDIF
!
IF(I%CISBA=='DIF')THEN
  IWORK=I%NGROUND_LAYER
ELSE
  IWORK=2 !Only 2 soil ice layer in ISBA-FR
ENDIF
!
DO JL=1,IWORK
  WRITE(YLVL,'(I4)') JL
  YRECFM='WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  I%XWGI(:,JL,:)=ZWORK(:,:)
END DO
!
!* water intercepted on leaves
!
ALLOCATE(I%XWR(ILU,I%NPATCH))
!
YRECFM = 'WR'
 CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XWR(:,:),IRESP)
!
!* Leaf Area Index
!
IF (I%CPHOTO=='LAI' .OR. I%CPHOTO=='LST' .OR. I%CPHOTO=='NIT' .OR. I%CPHOTO=='NCB') THEN
  YRECFM = 'LAI'
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XLAI(:,:),IRESP)
  IF ( LPRT ) THEN
    ! read in control variable
    IF ( TRIM(CVAR(NIVAR))=="LAI" ) THEN
      WHERE ( I%XLAI(:,:)/=XUNDEF ) 
        I%XLAI(:,:) = I%XLAI(:,:) + XTPRT(NIVAR)*I%XLAI(:,:)
      ENDWHERE
    ENDIF
  ENDIF  
END IF
!
!* snow mantel
!
 CALL READ_GR_SNOW(IOB, &
                   HPROGRAM,'VEG','     ',ILU,I%NPATCH,I%TSNOW  )
!
YRECFM='VERSION'
 CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,IVERSION,IRESP)
!
YRECFM='BUG'
 CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,IBUGFIX,IRESP)
!
IF(I%LGLACIER)THEN
  ALLOCATE(I%XICE_STO(ILU,I%NPATCH))
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=2) THEN
    YRECFM = 'ICE_STO'
    CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XICE_STO(:,:),IRESP)
  ELSE
    I%XICE_STO(:,:) = 0.0
  ENDIF
ELSE
  ALLOCATE(I%XICE_STO(0,0))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.  MEB Prognostic or Semi-prognostic variables
!            -------------------------------------------
!
ISIZE_LMEB_PATCH=COUNT(I%LMEB_PATCH(:))
!
IF (ISIZE_LMEB_PATCH>0) THEN
!
!* water intercepted on vegetation canopy leaves
!
  !nogv ALLOCATE(XWRV(ILU,NPATCH))
  !nogv YRECFM = 'WRV'
  !nogv CALL READ_SURF(HPROGRAM,YRECFM,XWRV(:,:),IRESP)
!
!* snow intercepted on vegetation canopy leaves
!
  ALLOCATE(I%XWRVN(ILU,I%NPATCH))
  YRECFM = 'WRVN'
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XWRVN(:,:),IRESP)
!
!* vegetation canopy temperature
!
  ALLOCATE(I%XTV(ILU,I%NPATCH))
  YRECFM = 'TV'
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XTV(:,:),IRESP)
!
!* vegetation canopy air temperature
!
  ALLOCATE(I%XTC(ILU,I%NPATCH))
  YRECFM = 'TC'
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XTC(:,:),IRESP)
!
!* vegetation canopy air specific humidity
!
  ALLOCATE(I%XQC(ILU,I%NPATCH))
  YRECFM = 'QC'
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XQC(:,:),IRESP)
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.  Semi-prognostic variables
!            -------------------------
!
ALLOCATE(I%XRESA(ILU,I%NPATCH))
ALLOCATE(I%XLE  (ILU,I%NPATCH))
IF (I%CPHOTO/='NON') THEN
  ALLOCATE(I%XANFM  (ILU,I%NPATCH))
  ALLOCATE(I%XAN    (ILU,I%NPATCH))
  ALLOCATE(I%XANDAY (ILU,I%NPATCH))
END IF
!
IF(I%CPHOTO/='NON') THEN
  ALLOCATE(I%XBIOMASS         (ILU,I%NNBIOMASS,I%NPATCH))
  ALLOCATE(I%XRESP_BIOMASS    (ILU,I%NNBIOMASS,I%NPATCH))
END IF
!
!
!* aerodynamical resistance
!
YRECFM = 'RESA'
I%XRESA(:,:) = 100.
 CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XRESA(:,:),IRESP)
!
!* patch averaged radiative temperature (K)
!
ALLOCATE(I%XTSRAD_NAT(ILU))
IF (IVERSION<6) THEN
  I%XTSRAD_NAT(:)=0.
  DO JP=1,I%NPATCH
    I%XTSRAD_NAT(:)=I%XTSRAD_NAT(:)+I%XTG(:,1,JP)
  ENDDO
  I%XTSRAD_NAT(:)=I%XTSRAD_NAT(:)/I%NPATCH
ELSE
  YRECFM='TSRAD_NAT'
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XTSRAD_NAT(:),IRESP)
ENDIF
!
I%XLE(:,:) = XUNDEF
!
!*       5. ISBA-AGS variables
!
IF (I%CPHOTO/='NON') THEN
  YRECFM = 'AN'
  I%XAN(:,:) = 0.
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XAN(:,:),IRESP)
  !
  YRECFM = 'ANDAY'
  I%XANDAY(:,:) = 0.
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XANDAY(:,:),IRESP)
  !
  YRECFM = 'ANFM'
  I%XANFM(:,:) = XANFMINIT
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XANFM(:,:),IRESP)
  !
  YRECFM = 'LE_AGS'
  I%XLE(:,:) = 0.
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,I%XLE(:,:),IRESP)
END IF
!
IF (I%CPHOTO=='AGS' .OR. I%CPHOTO=='AST') THEN
  !
  I%XBIOMASS(:,:,:) = 0.
  I%XRESP_BIOMASS(:,:,:) = 0.

ELSEIF (I%CPHOTO=='LAI' .OR. I%CPHOTO=='LST') THEN
  !
  I%XBIOMASS(:,1,:) = I%XBSLAI(:,:) * I%XLAI(:,:)
  I%XRESP_BIOMASS(:,:,:) = 0.

ELSEIF (I%CPHOTO=='NIT'.OR.I%CPHOTO=='NCB') THEN
  !
  I%XBIOMASS(:,:,:) = 0.
  DO JNBIOMASS=1,I%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
      YRECFM='BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='BIOMASS'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF
    CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    IF ( LPRT ) THEN
      ! read in control variable
      IF ( TRIM(CVAR(NIVAR)) == "LAI" .AND. TRIM(CBIO)==TRIM(YRECFM) ) THEN
        WHERE ( ZWORK(:,:)/=XUNDEF ) 
          ZWORK(:,:) = ZWORK(:,:) + XTPRT(NIVAR)*ZWORK(:,:)
        ENDWHERE
      ENDIF
    ENDIF     
    I%XBIOMASS(:,JNBIOMASS,:)=ZWORK
  END DO
!
  IWORK=0
  IF(I%CPHOTO=='NCB'.OR.IVERSION<8)IWORK=2
!
  I%XRESP_BIOMASS(:,:,:) = 0.
  DO JNBIOMASS=2,I%NNBIOMASS-IWORK
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. (IVERSION==7 .AND. IBUGFIX>=3)) THEN
      YRECFM='RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='RESP_BIOM'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF    
    CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    IF ( LPRT ) THEN
      ! read in control variable
      IF ( TRIM(CVAR(NIVAR)) == "LAI" .AND. TRIM(CBIO)==TRIM(YRECFM) ) THEN
        WHERE ( ZWORK(:,:)/=XUNDEF ) 
          ZWORK(:,:) = ZWORK(:,:) + XTPRT(NIVAR)*ZWORK(:,:)
        ENDWHERE
      ENDIF
    ENDIF      
    I%XRESP_BIOMASS(:,JNBIOMASS,:)=ZWORK
  END DO
  !
ENDIF
!
!*       6. Soil carbon
!
!
IF (I%CRESPSL=='CNT') THEN
  !
  ALLOCATE(I%XLITTER          (ILU,I%NNLITTER,I%NNLITTLEVS,I%NPATCH))
  ALLOCATE(I%XSOILCARB        (ILU,I%NNSOILCARB,I%NPATCH))
  ALLOCATE(I%XLIGNIN_STRUC    (ILU,I%NNLITTLEVS,I%NPATCH))
  !
  I%XLITTER(:,:,:,:) = 0.
  DO JNLITTER=1,I%NNLITTER
    DO JNLITTLEVS=1,I%NNLITTLEVS
      WRITE(YLVL,'(I1,A1,I1)') JNLITTER,'_',JNLITTLEVS
      YRECFM='LITTER'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
      I%XLITTER(:,JNLITTER,JNLITTLEVS,:)=ZWORK
    END DO
  END DO

  I%XSOILCARB(:,:,:) = 0.
  DO JNSOILCARB=1,I%NNSOILCARB
    WRITE(YLVL,'(I4)') JNSOILCARB
    YRECFM='SOILCARB'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    I%XSOILCARB(:,JNSOILCARB,:)=ZWORK
  END DO
!
  I%XLIGNIN_STRUC(:,:,:) = 0.
  DO JNLITTLEVS=1,I%NNLITTLEVS
    WRITE(YLVL,'(I4)') JNLITTLEVS
    YRECFM='LIGNIN_STR'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    I%XLIGNIN_STRUC(:,JNLITTLEVS,:)=ZWORK
  END DO
!
ENDIF

IF ( LASSIM ) THEN
  IF ( TRIM(CASSIM_ISBA) == "OI" ) THEN
    IF ( I%NPATCH /= 1 ) CALL ABOR1_SFX ('Reading of diagnostical values for'&
                       & //'assimilation at the moment only works for one patch')          
    ! Diagnostic fields for assimilation
    IF ( .NOT. ALLOCATED(XAT2M_ISBA)) ALLOCATE(XAT2M_ISBA(ILU,1))
    XAT2M_ISBA=XUNDEF
    YRECFM='T2M'
    CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,XAT2M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAHU2M_ISBA)) ALLOCATE(XAHU2M_ISBA(ILU,1))
    XAHU2M_ISBA=XUNDEF
    YRECFM='HU2M'
    CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,XAHU2M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAZON10M_ISBA)) ALLOCATE(XAZON10M_ISBA(ILU,1))
    XAZON10M_ISBA=XUNDEF
    YRECFM='ZON10M'
    CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,XAZON10M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAMER10M_ISBA)) ALLOCATE(XAMER10M_ISBA(ILU,1))
    XAMER10M_ISBA=XUNDEF
    YRECFM='MER10M'
    CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,XAMER10M_ISBA(:,1),IRESP)
  ELSEIF ( NIPERT/=NVAR+2 ) THEN
    ! Diagnostic fields for EKF assimilation ("observations")
    DO IOBS = 1,NOBSTYPE
     SELECT CASE (TRIM(COBS(IOBS)))
       CASE("T2M")
         IF ( .NOT. ALLOCATED(XAT2M_ISBA)) ALLOCATE(XAT2M_ISBA(ILU,I%NPATCH))
         XAT2M_ISBA=XUNDEF
         IF (I%NPATCH>1) THEN
           YRECFM='T2M_P'
         ELSE
           YRECFM='T2M_ISBA'
         ENDIF
         CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,XAT2M_ISBA(:,:),IRESP)
       CASE("HU2M")
         IF ( .NOT. ALLOCATED(XAHU2M_ISBA)) ALLOCATE(XAHU2M_ISBA(ILU,I%NPATCH))
         XAHU2M_ISBA=XUNDEF
         IF (I%NPATCH>1) THEN
           YRECFM='HU2M_P'
         ELSE
           YRECFM='HU2M_ISBA'
         ENDIF
         CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,XAHU2M_ISBA(:,:),IRESP)
       CASE("WG1")
         ! This is already read above
       CASE("LAI")
         ! This is already read above         
       CASE DEFAULT
         CALL ABOR1_SFX("Mapping of "//COBS(IOBS)//" is not defined in READ_ISBA_n!")
     END SELECT
    ENDDO
  ENDIF
ENDIF

!
!
DEALLOCATE(ZWORK)
IF (LHOOK) CALL DR_HOOK('READ_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_ISBA_n
