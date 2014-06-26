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
!!	V. Masson   *Meteo France*	
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
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_CO2V_PAR,       ONLY : XANFMINIT, XCONDCTMIN
USE MODD_ISBA_n,         ONLY : NGROUND_LAYER, NPATCH, NNBIOMASS,   &
                                  NNLITTER, NNLITTLEVS, NNSOILCARB,   &
                                  CPHOTO, CRESPSL, XTSRAD_NAT,        &
                                  XTG, XWG, XWGI, XWR, XLAI, TSNOW,   &
                                  XRESA, XANFM, XAN, XLE, XANDAY,     &
                                  XBSLAI, XBIOMASS, XRESP_BIOMASS,    &
                                  XLITTER, XSOILCARB, XLIGNIN_STRUC,  &
                                  LFLOOD, XZ0_FLOOD, LTEMP_ARP,       &
                                  NTEMPLAYER_ARP, LGLACIER, XICE_STO  
USE MODD_ASSIM,          ONLY : LASSIM,CASSIM_ISBA,XAT2M_ISBA,XAHU2M_ISBA,&
                              & XAZON10M_ISBA,XAMER10M_ISBA,              &
                              & COBS,NOBSTYPE,CVAR,LPRT,XTPRT,NIVAR
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
 CALL GET_TYPE_DIM_n('NATURE',ILU)
!
!
!*       2.     Prognostic fields:
!               -----------------
!
ALLOCATE(ZWORK(ILU,NPATCH))
!* soil temperatures
!
IF(LTEMP_ARP)THEN
  IWORK=NTEMPLAYER_ARP
ELSE
  IWORK=NGROUND_LAYER
ENDIF
!
ALLOCATE(XTG(ILU,IWORK,NPATCH))
!
DO JL=1,IWORK
  WRITE(YLVL,'(I4)') JL
  YRECFM='TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  XTG(:,JL,:)=ZWORK

  ! Perturb value if requested
  IF ( LPRT ) THEN
    ! read in control variable
    IF (( TRIM(CVAR(NIVAR)) == "TG1" ) .AND. ( JL == 1 ) ) THEN
      XTG(:,JL,:) = XTG(:,JL,:) + XTPRT(NIVAR)*XTG(:,JL,:)
    ENDIF
    IF (( TRIM(CVAR(NIVAR)) == "TG2" ) .AND. ( JL == 2 ) ) THEN
      XTG(:,JL,:) = XTG(:,JL,:) + XTPRT(NIVAR)*XTG(:,JL,:)
    ENDIF
  ENDIF
END DO
!
!
!* soil liquid and ice water contents
!
ALLOCATE(XWG (ILU,NGROUND_LAYER,NPATCH))
ALLOCATE(XWGI(ILU,NGROUND_LAYER,NPATCH))
!
XWG (:,:,:)=XUNDEF
XWGI(:,:,:)=XUNDEF
!
DO JL=1,NGROUND_LAYER
  WRITE(YLVL,'(I4)') JL
  YRECFM='WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
   CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
   XWG(:,JL,:)=ZWORK
   
   ! Perturb value if requested
  IF ( LPRT ) THEN
    ! read in control variable
    IF (( TRIM(CVAR(NIVAR)) == "WG1" ) .AND. ( JL == 1 ) ) THEN
      XWG(:,JL,:) = XWG(:,JL,:) + XTPRT(NIVAR)*XWG(:,JL,:)
    ENDIF
    IF (( TRIM(CVAR(NIVAR)) == "WG2" ) .AND. ( JL == 2 ) ) THEN
      XWG(:,JL,:) = XWG(:,JL,:) + XTPRT(NIVAR)*XWG(:,JL,:)
    ENDIF
  ENDIF
END DO
!
DO JL=1,NGROUND_LAYER
  WRITE(YLVL,'(I4)') JL
  YRECFM='WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  XWGI(:,JL,:)=ZWORK
END DO
!
!* water intercepted on leaves
!
ALLOCATE(XWR(ILU,NPATCH))
!
YRECFM = 'WR'
 CALL READ_SURF(HPROGRAM,YRECFM,XWR(:,:),IRESP)
!
!* roughness length of Flood water
!
IF(LFLOOD)THEN
  ALLOCATE(XZ0_FLOOD(ILU,NPATCH))
  YRECFM = 'Z0_FLOOD'
  CALL READ_SURF(HPROGRAM,YRECFM,XZ0_FLOOD(:,:),IRESP)
ENDIF
!
!* Leaf Area Index
!
IF (CPHOTO=='LAI' .OR. CPHOTO=='LST' .OR. CPHOTO=='NIT' .OR. CPHOTO=='NCB') THEN
  YRECFM = 'LAI'
  CALL READ_SURF(HPROGRAM,YRECFM,XLAI(:,:),IRESP)
  IF ( LPRT ) THEN
    ! read in control variable
    IF ( TRIM(CVAR(NIVAR)) == "LAI" ) THEN
      XLAI(:,:) = XLAI(:,:) + XTPRT(NIVAR)*XLAI(:,:)
    ENDIF
  ENDIF  
END IF
!
!* snow mantel
!
 CALL READ_GR_SNOW(HPROGRAM,'VEG','     ',ILU,NPATCH,TSNOW  )
!
YRECFM='VERSION'
 CALL READ_SURF(HPROGRAM,YRECFM,IVERSION,IRESP)
!
YRECFM='BUG'
 CALL READ_SURF(HPROGRAM,YRECFM,IBUGFIX,IRESP)
!
IF(LGLACIER)THEN
  ALLOCATE(XICE_STO(ILU,NPATCH))
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=2) THEN
    YRECFM = 'ICE_STO'
    CALL READ_SURF(HPROGRAM,YRECFM,XICE_STO(:,:),IRESP)
  ELSE
    XICE_STO(:,:) = 0.0
  ENDIF
ELSE
  ALLOCATE(XICE_STO(0,0))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.  Semi-prognostic variables
!            -------------------------
!
ALLOCATE(XRESA(ILU,NPATCH))
ALLOCATE(XLE  (ILU,NPATCH))
IF (CPHOTO/='NON') THEN
  ALLOCATE(XANFM  (ILU,NPATCH))
  ALLOCATE(XAN    (ILU,NPATCH))
  ALLOCATE(XANDAY (ILU,NPATCH))
END IF
!
IF(CPHOTO/='NON') THEN
  ALLOCATE(XBIOMASS         (ILU,NNBIOMASS,NPATCH))
  ALLOCATE(XRESP_BIOMASS    (ILU,NNBIOMASS,NPATCH))
END IF
!
!
!* aerodynamical resistance
!
YRECFM = 'RESA'
XRESA(:,:) = 100.
 CALL READ_SURF(HPROGRAM,YRECFM,XRESA(:,:),IRESP)
!
!* patch averaged radiative temperature (K)
!
ALLOCATE(XTSRAD_NAT(ILU))
IF (IVERSION<6) THEN
  XTSRAD_NAT(:)=0.
  DO JP=1,NPATCH
    XTSRAD_NAT(:)=XTSRAD_NAT(:)+XTG(:,1,JP)
  ENDDO
  XTSRAD_NAT(:)=XTSRAD_NAT(:)/NPATCH
ELSE
  YRECFM='TSRAD_NAT'
  CALL READ_SURF(HPROGRAM,YRECFM,XTSRAD_NAT(:),IRESP)
ENDIF
!
XLE(:,:) = XUNDEF
!
!*       5. ISBA-AGS variables
!
IF (CPHOTO/='NON') THEN
  YRECFM = 'AN'
  XAN(:,:) = 0.
  CALL READ_SURF(HPROGRAM,YRECFM,XAN(:,:),IRESP)
  !
  YRECFM = 'ANDAY'
  XANDAY(:,:) = 0.
  CALL READ_SURF(HPROGRAM,YRECFM,XANDAY(:,:),IRESP)
  !
  YRECFM = 'ANFM'
  XANFM(:,:) = XANFMINIT
  CALL READ_SURF(HPROGRAM,YRECFM,XANFM(:,:),IRESP)
  !
  YRECFM = 'LE_AGS'
  XLE(:,:) = 0.
  CALL READ_SURF(HPROGRAM,YRECFM,XLE(:,:),IRESP)
END IF
!
IF (CPHOTO=='AGS' .OR. CPHOTO=='AST') THEN
  !
  XBIOMASS(:,:,:) = 0.
  XRESP_BIOMASS(:,:,:) = 0.

ELSEIF (CPHOTO=='LAI' .OR. CPHOTO=='LST') THEN
  !
  XBIOMASS(:,1,:) = XBSLAI(:,:) * XLAI(:,:)
  XRESP_BIOMASS(:,:,:) = 0.

ELSEIF (CPHOTO=='NIT') THEN
  !
  XBIOMASS(:,:,:) = 0.
  DO JNBIOMASS=1,NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
      YRECFM='BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='BIOMASS'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF
    CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    XBIOMASS(:,JNBIOMASS,:)=ZWORK
  END DO

  XRESP_BIOMASS(:,:,:) = 0.
  DO JNBIOMASS=2,NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
      YRECFM='RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='RESP_BIOM'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF    
    CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    XRESP_BIOMASS(:,JNBIOMASS,:)=ZWORK
  END DO

ELSEIF (CPHOTO=='NCB') THEN
  !
  XBIOMASS(:,:,:) = 0.
  DO JNBIOMASS=1,NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
      YRECFM='BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='BIOMASS'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF    
    CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    XBIOMASS(:,JNBIOMASS,:)=ZWORK
  END DO

  XRESP_BIOMASS(:,:,:) = 0.
  DO JNBIOMASS=2,NNBIOMASS-2
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
      YRECFM='RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='RESP_BIOM'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF    
    CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    XRESP_BIOMASS(:,JNBIOMASS,:)=ZWORK
  END DO
  !
ENDIF
!
!*       6. Soil carbon
!
!
IF (CRESPSL=='CNT') THEN
  ALLOCATE(XLITTER          (ILU,NNLITTER,NNLITTLEVS,NPATCH))
  ALLOCATE(XSOILCARB        (ILU,NNSOILCARB,NPATCH))
  ALLOCATE(XLIGNIN_STRUC    (ILU,NNLITTLEVS,NPATCH))
END IF
!
IF (CRESPSL=='CNT') THEN
  !
  XLITTER(:,:,:,:) = 0.
  DO JNLITTER=1,NNLITTER
    DO JNLITTLEVS=1,NNLITTLEVS
      WRITE(YLVL,'(I1,A1,I1)') JNLITTER,'_',JNLITTLEVS
      YRECFM='LITTER'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
      XLITTER(:,JNLITTER,JNLITTLEVS,:)=ZWORK
    END DO
  END DO

  XSOILCARB(:,:,:) = 0.
  DO JNSOILCARB=1,NNSOILCARB
    WRITE(YLVL,'(I4)') JNSOILCARB
    YRECFM='SOILCARB'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    XSOILCARB(:,JNSOILCARB,:)=ZWORK
  END DO
!
  XLIGNIN_STRUC(:,:,:) = 0.
  DO JNLITTLEVS=1,NNLITTLEVS
    WRITE(YLVL,'(I4)') JNLITTLEVS
    YRECFM='LIGNIN_STR'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
    XLIGNIN_STRUC(:,JNLITTLEVS,:)=ZWORK
  END DO
!
ENDIF

IF ( LASSIM ) THEN
  IF ( NPATCH /= 1 ) CALL ABOR1_SFX ('Reading of diagnostical values for'&
                     & //'assimilation at the moment only works for one patch')
  IF ( TRIM(CASSIM_ISBA) == "OI" ) THEN
    ! Diagnostic fields for assimilation
    IF ( .NOT. ALLOCATED(XAT2M_ISBA)) ALLOCATE(XAT2M_ISBA(ILU,1))
    XAT2M_ISBA=XUNDEF
    YRECFM='T2M'
    CALL READ_SURF(HPROGRAM,YRECFM,XAT2M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAHU2M_ISBA)) ALLOCATE(XAHU2M_ISBA(ILU,1))
    XAHU2M_ISBA=XUNDEF
    YRECFM='HU2M'
    CALL READ_SURF(HPROGRAM,YRECFM,XAHU2M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAZON10M_ISBA)) ALLOCATE(XAZON10M_ISBA(ILU,1))
    XAZON10M_ISBA=XUNDEF
    YRECFM='ZON10M'
    CALL READ_SURF(HPROGRAM,YRECFM,XAZON10M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAMER10M_ISBA)) ALLOCATE(XAMER10M_ISBA(ILU,1))
    XAMER10M_ISBA=XUNDEF
    YRECFM='MER10M'
    CALL READ_SURF(HPROGRAM,YRECFM,XAMER10M_ISBA(:,1),IRESP)
  ELSE
    ! Diagnostic fields for EKF assimilation ("observations")
    DO IOBS = 1,NOBSTYPE
     SELECT CASE (TRIM(COBS(IOBS)))
       CASE("T2M")
         IF ( .NOT. ALLOCATED(XAT2M_ISBA)) ALLOCATE(XAT2M_ISBA(ILU,1))
         XAT2M_ISBA=XUNDEF
         YRECFM='T2M_ISBA'
         CALL READ_SURF(HPROGRAM,YRECFM,XAT2M_ISBA(:,1),IRESP)
       CASE("HU2M")
         IF ( .NOT. ALLOCATED(XAHU2M_ISBA)) ALLOCATE(XAHU2M_ISBA(ILU,1))
         XAHU2M_ISBA=XUNDEF
         YRECFM='HU2M_ISBA'
         CALL READ_SURF(HPROGRAM,YRECFM,XAHU2M_ISBA(:,1),IRESP)
       CASE("WG1")
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
