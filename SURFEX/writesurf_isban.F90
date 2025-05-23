!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE WRITESURF_ISBA_n (HSELECT, OSNOWDIMNC, CHI, NDST, &
                                   IO, S, NP, NPE, KI, HPROGRAM, OLAND_USE)
!     #####################################
!
!!****  *WRITESURF_ISBA_n* - writes ISBA prognostic fields
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
!
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
!!      P. LeMoigne 12/2004 : correct dimensionning if more than 10 layers in
!!                            the soil (diffusion version)
!!      B. Decharme  2008    : Floodplains
!!      B. Decharme  01/2009 : Optional Arpege deep soil temperature write
!!      A.L. Gibelin   03/09 : modifications for CENTURY model 
!!      A.L. Gibelin 04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin 06/2009 : Soil carbon variables for CNT option
!!      B. Decharme  07/2011 : land_use semi-prognostic variables
!!      B. Decharme  09/2012 : suppress NWG_LAYER (parallelization problems)
!!      B. Decharme  09/2012 : write some key for prep_read_external
!!      B. Decharme  04/2013 : Only 2 temperature layer in ISBA-FR
!!      P. Samuelsson 10/2014: MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURFEX_MPI, ONLY : NRANK
!
USE MODN_PREP_SURF_ATM, ONLY : LWRITE_EXTERN
USE MODD_WRITE_SURF_ATM,  ONLY : LSPLIT_PATCH
!
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_DST_n, ONLY : DST_NP_t
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_NP_t, ISBA_NPE_t, ISBA_S_t
!
USE MODD_SURF_PAR, ONLY : NUNDEF
!
USE MODD_ASSIM, ONLY : LASSIM, CASSIM, CASSIM_ISBA, NIE, NENS, &
                       XADDTIMECORR, LENS_GEN, NVAR
!
USE MODD_DST_SURF
!
USE MODI_WRITE_FIELD_1D_PATCH
USE MODI_WRITE_SURF
USE MODI_WRITESURF_GR_SNOW
USE MODI_ALLOCATE_GR_SNOW
USE MODI_DEALLOC_GR_SNOW
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT 
LOGICAL, INTENT(IN) :: OSNOWDIMNC  
!
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(DST_NP_t), INTENT(INOUT) :: NDST
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t), INTENT(INOUT) :: NPE
INTEGER, INTENT(IN) :: KI
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
LOGICAL,           INTENT(IN)  :: OLAND_USE !
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=4 ) :: YLVL
 CHARACTER(LEN=3 ) :: YVAR
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
 CHARACTER(LEN=25) :: YFORM          ! Writing format
  CHARACTER(LEN=2) :: YPAT
!
INTEGER :: JJ, JL, JP, JNB, JNL, JNS, JNLV  ! loop counter on levels
INTEGER :: IWORK   ! Work integer
INTEGER :: JSV
INTEGER :: ISIZE_LMEB_PATCH
INTEGER :: JVAR
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
!*       2.     Prognostic fields:
!               -----------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_ISBA_N',0,ZHOOK_HANDLE)
!
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
DO JL=1,IWORK
  WRITE(YLVL,'(I4)') JL
  YRECFM='TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A6,I1.1,A4)'
  IF (JL >= 10)  YFORM='(A6,I2.2,A4)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_TG',JL,' (K)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                NP%AL(JP)%NR_P,NPE%AL(JP)%XTG(:,JL),KI,S%XWORK_WR)
  ENDDO
END DO
!
!* soil liquid water contents
!
DO JL=1,IO%NGROUND_LAYER
  WRITE(YLVL,'(I4)') JL     
  YRECFM='WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A6,I1.1,A8)'
  IF (JL >= 10)  YFORM='(A6,I2.2,A8)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_WG',JL,' (m3/m3)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                NP%AL(JP)%NR_P,NPE%AL(JP)%XWG(:,JL),KI,S%XWORK_WR)  
  ENDDO
END DO
!
!* soil ice water contents
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
  YFORM='(A7,I1.1,A8)'
  IF (JL >= 10)  YFORM='(A7,I2.2,A8)'
  WRITE(YCOMMENT,YFORM) 'X_Y_WGI',JL,' (m3/m3)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                NP%AL(JP)%NR_P,NPE%AL(JP)%XWGI(:,JL),KI,S%XWORK_WR)  
  ENDDO
END DO
!
!* water intercepted on leaves
!
YRECFM='WR'
YCOMMENT='X_Y_WR (kg/m2)'
DO JP = 1,IO%NPATCH
  CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                NP%AL(JP)%NR_P,NPE%AL(JP)%XWR(:),KI,S%XWORK_WR)    
ENDDO
!
!* Glacier ice storage
!
YRECFM = 'GLACIER'
YCOMMENT='LGLACIER key for external prep'   
CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%LGLACIER,IRESP,HCOMMENT=YCOMMENT)
!
IF(IO%LGLACIER)THEN
  YRECFM='ICE_STO'
  YCOMMENT='X_Y_ICE_STO (kg/m2)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                NP%AL(JP)%NR_P,NPE%AL(JP)%XICE_STO(:),KI,S%XWORK_WR)   
  ENDDO
ENDIF
!
!* Leaf Area Index
!
IF (IO%CPHOTO/='NON' .AND. IO%CPHOTO/='AST') THEN
  !
  YRECFM='LAI'
  !
  YCOMMENT='X_Y_LAI (m2/m2)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                NP%AL(JP)%NR_P,NPE%AL(JP)%XLAI(:),KI,S%XWORK_WR)    
  ENDDO
  !
END IF
!
IF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. (LASSIM .OR. NIE/=0) ) THEN
  DO JVAR = 1,NVAR
    IF ( XADDTIMECORR(JVAR)>0. ) THEN
      WRITE(YVAR,'(I3)') JVAR
      YCOMMENT = 'Red_Noise_Enkf'
      YRECFM='RD_NS'//ADJUSTL(YVAR(:LEN_TRIM(YVAR)))
      DO JP = 1,IO%NPATCH
        CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                     NP%AL(JP)%NR_P,NP%AL(JP)%XRED_NOISE(:,JVAR),KI,S%XWORK_WR)       
      ENDDO
    ENDIF
  ENDDO
ENDIF
!
!* snow mantel
!
DO JP = 1,IO%NPATCH
  CALL WRITESURF_GR_SNOW(OSNOWDIMNC, HSELECT, HPROGRAM, 'VEG', '     ', KI, &
           NP%AL(JP)%NR_P, JP, NPE%AL(JP)%TSNOW, S%XWSN_WR, S%XRHO_WR, &
           S%XHEA_WR, S%XAGE_WR, S%XSG1_WR, S%XSG2_WR, S%XHIS_WR, S%XALB_WR, S%XIMP_WR)
ENDDO
!
!* key and/or field usefull to make an external prep
!
IF(IO%CISBA=='DIF')THEN
!
  YRECFM = 'SOC'
  YCOMMENT='SOC key for external prep'
  CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%LSOC,IRESP,HCOMMENT=YCOMMENT)
!
ELSE
!
  YRECFM = 'TEMPARP'
  YCOMMENT='LTEMP_ARP key for external prep'
  CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%LTEMP_ARP,IRESP,HCOMMENT=YCOMMENT)
!
  IF(IO%LTEMP_ARP)THEN
    YRECFM = 'NTEMPLARP'
    YCOMMENT='NTEMPLAYER_ARP for external prep'
    CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%NTEMPLAYER_ARP,IRESP,HCOMMENT=YCOMMENT)
  ENDIF
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.  MEB Prognostic or Semi-prognostic variables
!            -------------------------------------------
!
!
ISIZE_LMEB_PATCH=COUNT(IO%LMEB_PATCH(:))
!
IF (ISIZE_LMEB_PATCH>0) THEN
!
!* water intercepted on canopy vegetation leaves
!
  YRECFM='WRL'
  YCOMMENT='X_Y_WRL (kg/m2)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                  NP%AL(JP)%NR_P,NPE%AL(JP)%XWRL(:),KI,S%XWORK_WR)    
  ENDDO
!
!* ice on litter
!
  YRECFM='WRLI'
  YCOMMENT='X_Y_WRLI (kg/m2)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                  NP%AL(JP)%NR_P,NPE%AL(JP)%XWRLI(:),KI,S%XWORK_WR)    
  ENDDO
!
!* snow intercepted on canopy vegetation leaves
!
  YRECFM='WRVN'
  YCOMMENT='X_Y_WRVN (kg/m2)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                  NP%AL(JP)%NR_P,NPE%AL(JP)%XWRVN(:),KI,S%XWORK_WR)    
  ENDDO

!
!* canopy vegetation temperature
!
  YRECFM='TV'
  YCOMMENT='X_Y_TV (K)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                  NP%AL(JP)%NR_P,NPE%AL(JP)%XTV(:),KI,S%XWORK_WR)    
  ENDDO
!
!* litter temperature
!
  YRECFM='TL'
  YCOMMENT='X_Y_TL (K)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                  NP%AL(JP)%NR_P,NPE%AL(JP)%XTL(:),KI,S%XWORK_WR)    
  ENDDO
!
!* vegetation canopy air temperature
!
  YRECFM='TC'
  YCOMMENT='X_Y_TC (K)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                  NP%AL(JP)%NR_P,NPE%AL(JP)%XTC(:),KI,S%XWORK_WR)    
  ENDDO
!
!* vegetation canopy air specific humidity
!
  YRECFM='QC'
  YCOMMENT='X_Y_QC (kg/kg)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                  NP%AL(JP)%NR_P,NPE%AL(JP)%XQC(:),KI,S%XWORK_WR)    
  ENDDO
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.  Semi-prognostic variables
!            -------------------------
!
!
!* Fraction for each patch
!
YRECFM='PATCH'
YCOMMENT='fraction for each patch (-)'
DO JP = 1,IO%NPATCH
  CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                  NP%AL(JP)%NR_P,NP%AL(JP)%XPATCH(:),KI,S%XWORK_WR)    
ENDDO
!
!* patch averaged radiative temperature (K)
!
YRECFM='TSRAD_NAT'
YCOMMENT='X_TSRAD_NAT (K)'
 CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,S%XTSRAD_NAT(:),IRESP,HCOMMENT=YCOMMENT)
!
!* aerodynamical resistance
!
YRECFM='RESA'
YCOMMENT='X_Y_RESA (s/m)'
DO JP = 1,IO%NPATCH
  CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                  NP%AL(JP)%NR_P,NPE%AL(JP)%XRESA(:),KI,S%XWORK_WR)    
ENDDO
!
!* Land use variables
!
IF(OLAND_USE .OR. LWRITE_EXTERN)THEN
!
  DO JL=1,IO%NGROUND_LAYER
    WRITE(YLVL,'(I4)') JL
    YFORM='(A6,I1.1,A8)'
    IF (JL >= 10)  YFORM='(A6,I2.2,A8)'
    IF (OLAND_USE) THEN
      YRECFM='OLD_DG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      WRITE(YCOMMENT,FMT=YFORM) 'X_Y_OLD_DG',JL,' (m)'
    ELSE
      YRECFM='DG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      WRITE(YCOMMENT,FMT=YFORM) 'X_Y_DG',JL,' (m)'
    ENDIF
    DO JP = 1,IO%NPATCH
      CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NP%AL(JP)%XDG(:,JL),KI,S%XWORK_WR)    
    ENDDO
  END DO
!
ENDIF
!
!* ISBA-AGS variables
!
IF (IO%CPHOTO/='NON') THEN
  YRECFM='AN'
  YCOMMENT='X_Y_AN (kgCO2/kgair m/s)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                    NP%AL(JP)%NR_P,NPE%AL(JP)%XAN(:),KI,S%XWORK_WR)    
  ENDDO
!
  YRECFM='ANDAY'
  YCOMMENT='X_Y_ANDAY (kgCO2/m2/day)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                    NP%AL(JP)%NR_P,NPE%AL(JP)%XANDAY(:),KI,S%XWORK_WR)    
  ENDDO  
!
  YRECFM='ANFM'
  YCOMMENT='X_Y_ANFM (kgCO2/kgair m/s)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                    NP%AL(JP)%NR_P,NPE%AL(JP)%XANFM(:),KI,S%XWORK_WR)    
  ENDDO    
!
  YRECFM='LE_AGS'
  YCOMMENT='X_Y_LE_AGS (W/m2)'
  DO JP = 1,IO%NPATCH
    CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                    NP%AL(JP)%NR_P,NPE%AL(JP)%XLE(:),KI,S%XWORK_WR)    
  ENDDO    
END IF
!
!
IF (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
  !
  DO JNB=1,IO%NNBIOMASS
    WRITE(YLVL,'(I1)') JNB
    YRECFM='BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YFORM='(A11,I1.1,A10)'
    WRITE(YCOMMENT,FMT=YFORM) 'X_Y_BIOMASS',JNB,' (kgDM/m2)'
    DO JP = 1,IO%NPATCH
      CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XBIOMASS(:,JNB),KI,S%XWORK_WR)    
    ENDDO      
  END DO
  !
  !
  DO JNB=2,IO%NNBIOMASS
    WRITE(YLVL,'(I1)') JNB
    YRECFM='RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YFORM='(A16,I1.1,A10)'
    WRITE(YCOMMENT,FMT=YFORM) 'X_Y_RESP_BIOMASS',JNB,' (kg/m2/s)'
    DO JP = 1,IO%NPATCH
      CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XRESP_BIOMASS(:,JNB),KI,S%XWORK_WR)    
    ENDDO      
  END DO
  !
END IF
!
!* Soil carbon
!
YRECFM = 'RESPSL'
YCOMMENT=YRECFM
 CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%CRESPSL,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='NLITTER'
YCOMMENT=YRECFM
 CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%NNLITTER,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='NLITTLEVS'
YCOMMENT=YRECFM
 CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%NNLITTLEVS,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='NSOILCARB'
YCOMMENT=YRECFM
 CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%NNSOILCARB,IRESP,HCOMMENT=YCOMMENT)
!
IF(IO%LSPINUPCARBS.OR.IO%LSPINUPCARBW)THEN
  YRECFM='NBYEARSOLD'
  YCOMMENT='yrs'
  CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%NNBYEARSOLD,IRESP,HCOMMENT=YCOMMENT)
ENDIF
!
IF (IO%CRESPSL=='CNT') THEN
  !
  DO JNL=1,IO%NNLITTER
    DO JNLV=1,IO%NNLITTLEVS
      WRITE(YLVL,'(I1,A1,I1)') JNL,'_',JNLV
      YRECFM='LITTER'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      YFORM='(A10,I1.1,A1,I1.1,A8)'
      WRITE(YCOMMENT,FMT=YFORM) 'X_Y_LITTER',JNL,' ',JNLV,' (gC/m2)'
      DO JP = 1,IO%NPATCH
        CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                        NP%AL(JP)%NR_P,NPE%AL(JP)%XLITTER(:,JNL,JNLV),KI,S%XWORK_WR)    
      ENDDO        
    END DO
  END DO

  DO JNS=1,IO%NNSOILCARB
    WRITE(YLVL,'(I4)') JNS
    YRECFM='SOILCARB'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YFORM='(A8,I1.1,A8)'
    WRITE(YCOMMENT,FMT=YFORM) 'X_Y_SOILCARB',JNS,' (gC/m2)'
    DO JP = 1,IO%NPATCH
      CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XSOILCARB(:,JNS),KI,S%XWORK_WR)    
    ENDDO     
  END DO
!
  DO JNLV=1,IO%NNLITTLEVS
    WRITE(YLVL,'(I4)') JNLV
    YRECFM='LIGN_STR'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YFORM='(A12,I1.1,A8)'
    WRITE(YCOMMENT,FMT=YFORM) 'X_Y_LIGNIN_STRUC',JNLV,' (-)'
    DO JP = 1,IO%NPATCH
      CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XLIGNIN_STRUC(:,JNLV),KI,S%XWORK_WR)    
    ENDDO       
  END DO
!
ENDIF
!
!
IF (CHI%SVI%NDSTEQ > 0)THEN
  DO JSV = 1,NDSTMDE ! for all dust modes
    WRITE(YRECFM,'(A6,I3.3)')'F_DSTM',JSV
    YCOMMENT='X_Y_'//YRECFM//' (kg/m2)'
    DO JP = 1,IO%NPATCH
      CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NDST%AL(JP)%XSFDSTM(:,JSV),KI,S%XWORK_WR)    
    ENDDO     
  END DO
ENDIF
!
!-------------------------------------------------------------------------------

!*       5.  Time
!            ----
!
YRECFM='DTCUR'
YCOMMENT='s'
 CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,S%TTIME,IRESP,HCOMMENT=YCOMMENT)
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_ISBA_n
