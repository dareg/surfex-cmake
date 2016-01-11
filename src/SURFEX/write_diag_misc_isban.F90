!     #########
      SUBROUTINE WRITE_DIAG_MISC_ISBA_n (DTCO, DGU, U, DGI, DGIP, DGMI, DGMIP, I, &
                                         HPROGRAM)
!     #################################
!
!!****  *WRITE_DIAG_MISC_ISBA* - writes the ISBA diagnostic fields
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2004
!!      B. Decharme    2008  Total Albedo, Total SWI and Floodplains
!!      B. Decharme 06/2009  key to write (or not) patch result
!!      A.L. Gibelin 04/09 : Add respiration diagnostics
!!      A.L. Gibelin 05/09 : Add carbon spinup
!!      A.L. Gibelin 07/09 : Suppress RDK and transform GPP as a diagnostic
!!      D. Carrer    04/11 : Add FAPAR and effective LAI
!!      B. Decharme  09/2012 : suppress NWG_LAYER (parallelization problems)
!!      B. Decharme  09/12 : Carbon fluxes in diag_evap
!!      B. Decharme  09/12   New diag for DIF:
!!                           F2 stress
!!                           Root zone swi, wg and wgi
!!                           swi, wg and wgi comparable to ISBA-FR-DG2 and DG3 layers
!!                           active layer thickness over permafrost
!!                           frozen layer thickness over non-permafrost
!!      B. Decharme  06/13   All snow outputs noted SN
!!                           XTSRAD_NAT instead of XAVG_TSRAD
!!                           delete NWG_SIZE
!!                           water table depth
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_DIAG_n, ONLY : DIAG_t, DIAG_PATCH_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t, DIAG_MISC_ISBA_PATCH_t
USE MODD_ISBA_n, ONLY : ISBA_t
!
USE MODD_SURF_PAR,        ONLY :   NUNDEF, XUNDEF
!
USE MODD_ASSIM, ONLY : LASSIM, CASSIM_ISBA, NVAR
!                                 
USE MODD_AGRI,            ONLY :   LAGRIP
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_t), INTENT(INOUT) :: DGU
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(DIAG_t), INTENT(INOUT) :: DGI
TYPE(DIAG_PATCH_t), INTENT(INOUT) :: DGIP
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
TYPE(DIAG_MISC_ISBA_PATCH_t), INTENT(INOUT) :: DGMIP
TYPE(ISBA_t), INTENT(INOUT) :: I
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=1) :: YVAR
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
 CHARACTER(LEN=2)  :: YLVL
 CHARACTER(LEN=20) :: YFORM
!
INTEGER           :: JLAYER, JJ, IDEPTH, JVAR, JPATCH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_MISC_ISBA_N',0,ZHOOK_HANDLE)
CALL INIT_IO_SURF_n(DTCO, DGU, U, &
                     HPROGRAM,'NATURE','ISBA  ','WRITE')
!
!-------------------------------------------------------------------------------
!
IF (DGMI%LSURF_MISC_BUDGET) THEN
  !
  !*       2.     Miscellaneous fields :
  !
  !-------------------------------------------------------------------------------
  !
  !        2.1    Halstead coefficient
  !               --------------------
  !
  YRECFM='HV_ISBA'
  YCOMMENT='Halstead coefficient averaged over tile nature (-)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XHV(:),IRESP,HCOMMENT=YCOMMENT)
  !
  !        2.2    Snow fractions
  !               --------------
  !
  YRECFM='PSNG_ISBA'
  YCOMMENT='snow fraction over ground averaged over tile nature (-)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XPSNG(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='PSNV_ISBA'
  YCOMMENT='snow fraction over vegetation averaged over tile nature (-)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XPSNV(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='PSN_ISBA'
  YCOMMENT='total snow fraction averaged over tile nature (-)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XPSN(:),IRESP,HCOMMENT=YCOMMENT)
  !
  !        2.3    Total Albedo and surface temperature
  !               ------------------------------------
  !
  IF (I%R%TSNOW%SCHEME=='3-L' .OR. I%R%TSNOW%SCHEME=='CRO') THEN
    !        
    YRECFM='TS_ISBA'
    YCOMMENT='total surface temperature (isba+snow) over tile nature'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGI%XTS(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='TSRAD_ISBA'
    YCOMMENT='total radiative surface temperature (isba+snow) over tile nature'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,I%R%XTSRAD_NAT(:),IRESP,HCOMMENT=YCOMMENT)
    !
  END IF
  !
  !        2.4    Soil Wetness Index, Water content and active layer depth
  !               --------------------------------------------------------
  !  
  IF(I%O%CISBA=='DIF')THEN
    DO JLAYER = 1,I%O%NGROUND_LAYER
     DO JJ=1,SIZE(I%M%X%NWG_LAYER,1)
        IDEPTH=MAXVAL(I%M%X%NWG_LAYER(JJ,:),I%M%X%NWG_LAYER(JJ,:)/=NUNDEF)
        IF(JLAYER>IDEPTH)THEN  
          DGMI%XSWI (JJ,JLAYER) = XUNDEF
          DGMI%XTSWI(JJ,JLAYER) = XUNDEF
        ENDIF
      ENDDO 
    ENDDO
  ENDIF         
  !
  DO JLAYER=1,I%O%NGROUND_LAYER
    !
    WRITE(YLVL,'(I2)') JLAYER
    !
    YRECFM='SWI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YRECFM=YRECFM(:LEN_TRIM(YRECFM))//'_ISBA'
    YFORM='(A29,I1.1,A4)'
    IF (JLAYER >= 10)  YFORM='(A29,I2.2,A4)'
    WRITE(YCOMMENT,FMT=YFORM) 'soil wetness index for layer ',JLAYER,' (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XSWI(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='TSWI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YRECFM=YRECFM(:LEN_TRIM(YRECFM))//'_ISBA'
    YFORM='(A29,I1.1,A4)'
    IF (JLAYER >= 10)  YFORM='(A29,I2.2,A4)'
    WRITE(YCOMMENT,FMT=YFORM) 'total swi (liquid+solid) for layer ',JLAYER,' (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XTSWI(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
    !
  END DO
  !
  YRECFM='SWI_T_ISBA'
  YCOMMENT='soil wetness index over the soil column (-)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XSOIL_SWI(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='TSWI_T_ISBA'
  YCOMMENT='total soil wetness index over the soil column (-)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XSOIL_TSWI(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='WGTOT_T_ISBA'
  YCOMMENT='total water content (liquid+solid) over the soil column (kg/m2)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XSOIL_TWG(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='WGI_T_ISBA'
  YCOMMENT='total ice content (solid) over the soil column (kg/m2)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XSOIL_TWGI(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='WGTOT_ISBA'
  YCOMMENT='total volumetric water content (liquid+solid) over the soil column (m3/m3)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XSOIL_WG(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='WGI_ISBA'
  YCOMMENT='total volumetric ice content (solid) over the soil column (m3/m3)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XSOIL_WGI(:),IRESP,HCOMMENT=YCOMMENT)
  !
  IF(I%O%CISBA=='DIF') THEN
    !
    IF (DGMI%LSURF_MISC_DIF)THEN
      !
      YRECFM='TSWI_D2_ISBA'
      YCOMMENT='total soil wetness index over comparable FR-DG2 reservoir (-)'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XFRD2_TSWI(:),IRESP,HCOMMENT=YCOMMENT)
      !
      YRECFM='WG_D2_ISBA'
      YCOMMENT='liquid water content over comparable FR-DG2 reservoir (m3/m3)'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XFRD2_TWG(:),IRESP,HCOMMENT=YCOMMENT)
      !
      YRECFM='WGI_D2_ISBA'
      YCOMMENT='ice content over comparable FR-DG2 reservoir (m3/m3)'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XFRD2_TWGI(:),IRESP,HCOMMENT=YCOMMENT)  
      !
      YRECFM='TSWI_D3_ISBA'
      YCOMMENT='total soil wetness index over comparable FR-DG3 reservoir (-)'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XFRD3_TSWI(:),IRESP,HCOMMENT=YCOMMENT)
      !
      YRECFM='WG_D3_ISBA'
      YCOMMENT='liquid water content over comparable FR-DG3 reservoir (m3/m3)'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XFRD3_TWG(:),IRESP,HCOMMENT=YCOMMENT)
      !
      YRECFM='WGI_D3_ISBA'
      YCOMMENT='ice content over comparable FR-DG3 reservoir (m3/m3)'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XFRD3_TWGI(:),IRESP,HCOMMENT=YCOMMENT)  
      !
    ENDIF
    !
    YRECFM='ALT_ISBA'
    YCOMMENT='active layer thickness over permafrost (m)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XALT(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='FLT_ISBA'
    YCOMMENT='frozen layer thickness over non-permafrost (m)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XFLT(:),IRESP,HCOMMENT=YCOMMENT)
    !
  ENDIF
  !
  !        2.5    Snow outputs
  !               -------------
  !
  YRECFM='WSN_T_ISBA'
  YCOMMENT='Total_snow_reservoir (kg/m2)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XTWSNOW(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='DSN_T_ISBA'
  YCOMMENT='Total_snow_depth (m)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XTDSNOW(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='TSN_T_ISBA'
  YCOMMENT='Total_snow_temperature (K)'
  CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XTTSNOW(:),IRESP,HCOMMENT=YCOMMENT)
  !
  !        2.6    SGH scheme
  !               ----------
  !
  IF(I%O%CRUNOFF=='SGH '.OR.I%O%CRUNOFF=='DT92')THEN     
    YRECFM='FSAT_ISBA'
    YCOMMENT='Soil saturated fraction (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XFSAT(:),IRESP,HCOMMENT=YCOMMENT)
  ENDIF
  !
  IF(I%O%CRAIN=='SGH ')THEN
    YRECFM='MUF_ISBA'
    YCOMMENT='fraction of the grid cell reached by the rainfall (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,I%I%XMUF(:),IRESP,HCOMMENT=YCOMMENT)
  ENDIF
  !
  !        2.7    Flooding scheme
  !               ---------------
  !
  IF(I%O%LFLOOD)THEN
    !
    YRECFM='FFG_ISBA'
    YCOMMENT='flood fraction over ground averaged over tile nature (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XFFG(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='FFV_ISBA'
    YCOMMENT='flood fraction over vegetation averaged over tile nature (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XFFV(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='FF_ISBA'
    YCOMMENT='total flood fraction averaged over tile nature (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XFF(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='FFLOOD_ISBA'
    YCOMMENT='Grdi-cell potential flood fraction (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,I%I%XFFLOOD(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='PIFLOOD_ISBA'
    YCOMMENT='Grdi-cell Potential_floodplain_infiltration (kg/m2)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,I%I%XPIFLOOD(:),IRESP,HCOMMENT=YCOMMENT)
    !
  ENDIF
  !
  !        2.8    Total LAI
  !               ---------
  !
  IF(I%O%CPHOTO/='NON'.OR.I%O%NPATCH>1)THEN        
    YRECFM='LAI_ISBA'
    YCOMMENT='leaf area index (m2/m2)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,DGMI%XLAI(:),IRESP,HCOMMENT=YCOMMENT)
  ENDIF
  !
  !        2.9    Water table depth
  !               -----------------
  !
  IF(I%O%LWTD)THEN
    !
    YRECFM='FWTD_ISBA'
    YCOMMENT='grid-cell fraction of water table to rise'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,I%IP%XFWTD(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='WTD_ISBA'
    YCOMMENT='water table depth from RRM model or observation (m)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,I%IP%XWTD(:),IRESP,HCOMMENT=YCOMMENT)
    !
  ENDIF
  !*       3.     Miscellaneous fields for each patch :
  !               -------------------------------------
  !
  !----------------------------------------------------------------------------
  !User wants (or not) patch output
  IF(DGI%LPATCH_BUDGET)THEN
    !----------------------------------------------------------------------------
    !
    !        3.1    Soil Wetness Index and active layer depth
    !               -----------------------------------------   
    !
    ALLOCATE(ZWORK(SIZE(DGMIP%AL(1)%XSWI(:,1)),I%O%NPATCH))
    !
    DO JLAYER=1,I%O%NGROUND_LAYER
      !
      WRITE(YLVL,'(I2)') JLAYER
      !
      YRECFM='SWI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      YFORM='(A39,I1.1,A4)'
      IF (JLAYER >= 10)  YFORM='(A39,I2.2,A4)'
      WRITE(YCOMMENT,FMT=YFORM) 'soil wetness index per patch for layer ',JLAYER,' (-)'
      DO JPATCH=1,I%O%NPATCH
        ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XSWI(:,JLAYER)
      ENDDO
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
      !
      YRECFM='TSWI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      YFORM='(A39,I1.1,A4)'
      IF (JLAYER >= 10)  YFORM='(A39,I2.2,A4)'
      WRITE(YCOMMENT,FMT=YFORM) 'total swi (liquid+solid) per patch for layer ',JLAYER,' (-)'
      DO JPATCH=1,I%O%NPATCH
        ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XTSWI(:,JLAYER)
      ENDDO      
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
      !
    END DO
    !
    IF(I%O%CISBA=='DIF')THEN
      !
      DO JPATCH=1,I%O%NPATCH
        ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XALT(:)
      ENDDO      
      YRECFM='ALT_P'
      YCOMMENT='active layer thickness over permafrost per patch (m)'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
      !
      DO JPATCH=1,I%O%NPATCH
        ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XFLT(:)
      ENDDO
      YRECFM='FLT_P'
      YCOMMENT='frozen layer thickness over non-permafrost per patch (m)'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT) 
      !
    ENDIF
    !    
    !        3.2    Snow fractions
    !               --------------
    !
    DO JPATCH=1,I%O%NPATCH
      ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XPSNG(:)
    ENDDO    
    YRECFM='PSNG_P'
    YCOMMENT='snow fraction per patch over ground '
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    DO JPATCH=1,I%O%NPATCH
      ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XPSNV(:)
    ENDDO
    YRECFM='PSNV_P'
    YCOMMENT='snow fraction per patch over vegetation'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    DO JPATCH=1,I%O%NPATCH
      ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XPSN(:)
    ENDDO
    YRECFM='PSN_P'
    YCOMMENT='total snow fraction per patch'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    !        3.3    SGH scheme
    !               ----------
    !
    IF(I%O%CRUNOFF=='DT92')THEN 
      DO JPATCH=1,I%O%NPATCH
        ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XFSAT(:)
      ENDDO      
      YRECFM='FSAT_P'
      YCOMMENT='Soil saturated fraction per patch (-)'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    ENDIF
    !
    !        3.3    Flood fractions
    !               --------------
    !
    IF(I%O%LFLOOD)THEN
      !        
      DO JPATCH=1,I%O%NPATCH
        ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XFFG(:)
      ENDDO
      YRECFM='FFG_P'
      YCOMMENT='flood fraction per patch over ground '
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
      !
      DO JPATCH=1,I%O%NPATCH
        ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XFFV(:)
      ENDDO
      YRECFM='FFV_P'
      YCOMMENT='flood fraction per patch over vegetation'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
      !
      DO JPATCH=1,I%O%NPATCH
        ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XFF(:)
      ENDDO
      YRECFM='FF_P'
      YCOMMENT='total flood fraction per patch'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
      !
    ENDIF
    !
    !        3.4    Total Albedo
    !               ------------
    !
    !
    IF (I%R%TSNOW%SCHEME=='3-L' .OR. I%R%TSNOW%SCHEME=='CRO') THEN
      DO JPATCH=1,I%O%NPATCH
        ZWORK(:,JPATCH) = DGIP%AL(JPATCH)%XTS(:)
      ENDDO
      YRECFM='TS_P'
      YCOMMENT='total surface temperature (isba+snow) per patch'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK(:,:),IRESP,HCOMMENT=YCOMMENT)
      DO JPATCH=1,I%O%NPATCH
        ZWORK(:,JPATCH) = DGIP%AL(JPATCH)%XTSRAD(:)
      ENDDO
      YRECFM='TSRAD_P'
      YCOMMENT='total radiative surface temperature (isba+snow) per patch'
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK(:,:),IRESP,HCOMMENT=YCOMMENT)
    ENDIF
    !
    !        3.5    Halstead coefficient
    !               --------------------
    !
    DO JPATCH=1,I%O%NPATCH
       ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XHV(:)
    ENDDO
    YRECFM='HV'
    YCOMMENT='Halstead coefficient per patch'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    !        3.6  Snow outputs 
    !        -----------------
    !
    DO JPATCH=1,I%O%NPATCH
       ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XTWSNOW(:)
    ENDDO
    YRECFM='WSN_T_P'
    YCOMMENT='X_Y_WSNOW_TOT (kg/m2) per patch'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    DO JPATCH=1,I%O%NPATCH
       ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XTDSNOW(:)
    ENDDO
    YRECFM='DSN_T_P'
    YCOMMENT='X_Y_DSNOW_TOT (m) per patch'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    DO JPATCH=1,I%O%NPATCH
       ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XTTSNOW(:)
    ENDDO
    YRECFM='TSN_T_P'
    YCOMMENT='X_Y_TSNOW_TOT (k) per patch'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    IF (I%R%TSNOW%SCHEME=='3-L' .OR. I%R%TSNOW%SCHEME=='CRO') THEN
      !
      DO JLAYER=1,I%R%TSNOW%NLAYER
        !
        WRITE(YLVL,'(I2)') JLAYER
        !
        YRECFM='SNOWLIQ'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
        YFORM='(A17,I1.1,A4)'
        IF (JLAYER >= 10)  YFORM='(A17,I2.2,A4)'
        WRITE(YCOMMENT,FMT=YFORM) 'snow liquid water',JLAYER,' (m)'
        DO JPATCH=1,I%O%NPATCH
          ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XSNOWLIQ(:,JLAYER)
        ENDDO
        CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
        !
        YRECFM='SNOWTEMP'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
        YFORM='(A16,I1.1,A4)'
        IF (JLAYER >= 10)  YFORM='(A16,I2.2,A4)'
        WRITE(YCOMMENT,FMT=YFORM) 'snow temperature',JLAYER,' (K)'
        DO JPATCH=1,I%O%NPATCH
          ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XSNOWTEMP(:,JLAYER)
        ENDDO
        CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
        !
      END DO
      !        
    ENDIF
    !
    DEALLOCATE(ZWORK)
    !
  END IF
  !
  IF (LAGRIP) THEN
    !
    !        2.8    Irrigation threshold
    !               --------------------
    !
    ALLOCATE(ZWORK(SIZE(DGMIP%AL(1)%XSEUIL(:)),I%O%NPATCH))
    DO JPATCH=1,I%O%NPATCH
      ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XSEUIL(:)
    ENDDO
    YRECFM='IRRISEUIL'
    YCOMMENT='irrigation threshold per patch'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    DEALLOCATE(ZWORK)
    !
  ENDIF
  !
  IF (I%O%LTR_ML) THEN
    !
    ALLOCATE(ZWORK(SIZE(DGMIP%AL(1)%XFAPAR(:)),I%O%NPATCH))
    !
    DO JPATCH=1,I%O%NPATCH
      ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XFAPAR(:)
    ENDDO
    YRECFM='FAPAR'
    YCOMMENT='FAPAR (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    DO JPATCH=1,I%O%NPATCH
      ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XFAPIR(:)
    ENDDO
    YRECFM='FAPIR'
    YCOMMENT='FAPIR (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    DO JPATCH=1,I%O%NPATCH
      ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XFAPAR_BS(:)
    ENDDO
    YRECFM='FAPAR_BS'
    YCOMMENT='FAPAR_BS (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    DO JPATCH=1,I%O%NPATCH
      ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XFAPIR_BS(:)
    ENDDO
    YRECFM='FAPIR_BS'
    YCOMMENT='FAPIR_BS (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    DO JPATCH=1,I%O%NPATCH
      ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XDFAPARC(:)
    ENDDO
    YRECFM='DFAPARC'
    YCOMMENT='DFAPARC (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    DO JPATCH=1,I%O%NPATCH
      ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XDFAPIRC(:)
    ENDDO
    YRECFM='DFAPIRC'
    YCOMMENT='DFAPIRC (-)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    DO JPATCH=1,I%O%NPATCH
      ZWORK(:,JPATCH) = DGMIP%AL(JPATCH)%XDLAI_EFFC(:)
    ENDDO
    YRECFM='DLAI_EFFC'
    YCOMMENT='DLAI_EFFC (m2/m2)'
    CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
    !
    DEALLOCATE(ZWORK)
    !
  ENDIF
  ! 
  IF (LASSIM .AND. CASSIM_ISBA=="EKF  ") THEN
    !
    DO JVAR = 1,NVAR
      WRITE(YVAR,FMT='(I1.1)') JVAR
      YRECFM="ANAL_INCR"//YVAR
      YCOMMENT="by patch"
      CALL WRITE_SURF(DGU, U, HPROGRAM,YRECFM,I%R%XINCR(:,I%O%NPATCH*(JVAR-1)+1:I%O%NPATCH*JVAR),IRESP,HCOMMENT=YCOMMENT)
    ENDDO
    !
  ENDIF
  !
ENDIF
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_MISC_ISBA_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_DIAG_MISC_ISBA_n
