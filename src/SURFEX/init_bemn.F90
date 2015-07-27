!     #############################################################
      SUBROUTINE INIT_BEM_n (IOB, &
                              B, BOP, BDD, DTB, DTCO, DTT, UG, U, TG, T, TOP, &
                             KLUOUT)
!     #############################################################
!
!!****  *INIT_BEM_n* - routine to initialize Building Energy Model
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
!!      G. Pigeon   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2012
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
!
!
!
!
!
USE MODD_CH_EMIS_FIELD_n, ONLY : CHE => CH_EMIS_FIELD
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
USE MODD_CH_SEAFLUX_n, ONLY : CHS => CH_SEAFLUX
USE MODD_CH_SNAP_n, ONLY : CHN => CH_EMIS_SNAP
USE MODD_CH_SURF_n, ONLY : CHU => CH_SURF
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
USE MODD_CH_WATFLUX_n, ONLY : CHW => CH_WATFLUX
USE MODD_DATA_SEAFLUX_n, ONLY : DTS => DATA_SEAFLUX
USE MODD_DATA_TSZ0_n, ONLY : DTZ => DATA_TSZ0
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DGEI => DIAG_EVAP_ISBA
USE MODD_DIAG_FLAKE_n, ONLY : DGF => DIAG_FLAKE
USE MODD_DIAG_ISBA_n, ONLY : DGI => DIAG_ISBA
USE MODD_DIAG_MISC_ISBA_n, ONLY : DGMI => DIAG_MISC_ISBA
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
USE MODD_DIAG_OCEAN_n, ONLY : DGO => DIAG_OCEAN
USE MODD_DIAG_SEAFLUX_n, ONLY : DGS => DIAG_SEAFLUX
USE MODD_DIAG_SEAICE_n, ONLY : DGSI => DIAG_SEAICE
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
USE MODD_DIAG_UTCI_TEB_n, ONLY : DGUT => DIAG_UTCI_TEB
USE MODD_DIAG_WATFLUX_n, ONLY : DGW => DIAG_WATFLUX
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_FLAKE_SBL_n, ONLY : FSB => FLAKE_SBL
USE MODD_GR_BIOG_n, ONLY : GB => GR_BIOG
USE MODD_ISBA_CANOPY_n, ONLY : ICP => ISBA_CANOPY
USE MODD_ISBA_n, ONLY : I => ISBA
USE MODD_OCEAN_n, ONLY : O => OCEAN
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
USE MODD_SEAFLUX_SBL_n, ONLY : SSB => SEAFLUX_SBL
USE MODD_SV_n, ONLY : SV => SV
USE MODD_TEB_CANOPY_n, ONLY : TCP => TEB_CANOPY
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GREENROOF_n, ONLY : TGR => TEB_GREENROOF
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
USE MODD_WATFLUX_SBL_n, ONLY : WSB => WATFLUX_SBL
!
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
!
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_t
USE MODD_DATA_BEM_n, ONLY : DATA_BEM_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
!
USE MODI_CONVERT_PATCH_TEB
USE MODI_WINDOW_DATA
USE MODI_HVAC_AUTOSIZE
USE MODI_BEM_MORPHO
USE MODI_STORES_HVAC_AUTOSIZE
!
USE MODI_ABOR1_SFX
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
!
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
!
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BLD_DESC_t), INTENT(INOUT) :: BDD
TYPE(DATA_BEM_t), INTENT(INOUT) :: DTB
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_TEB_t), INTENT(INOUT) :: DTT
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
!
INTEGER, INTENT(IN) :: KLUOUT ! logical unit of output listing
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER                         :: JJ               ! counter
INTEGER                         :: ILU              ! sizes of TEB arrays
LOGICAL                         :: GPRINT           ! flag for warning prints in output file
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!

IF (LHOOK) CALL DR_HOOK('INIT_BEM_N',0,ZHOOK_HANDLE)
!
!
!*       3.     Physiographic data fields from land cover:
!               -----------------------------------------
!
ILU = SIZE(TOP%XCOVER,1)
IF (TOP%CBEM=='DEF') ILU=0
!
ALLOCATE(B%CUR%XHC_FLOOR    (ILU,BOP%NFLOOR_LAYER))
ALLOCATE(B%CUR%XTC_FLOOR    (ILU,BOP%NFLOOR_LAYER))
ALLOCATE(B%CUR%XD_FLOOR     (ILU,BOP%NFLOOR_LAYER))
!
ALLOCATE(B%CUR%XTCOOL_TARGET(ILU))
ALLOCATE(B%CUR%XTHEAT_TARGET(ILU))
ALLOCATE(B%CUR%XEFF_HEAT    (ILU))
ALLOCATE(B%CUR%XSHGC        (ILU))
ALLOCATE(B%CUR%XQIN         (ILU))
ALLOCATE(B%CUR%XQIN_FRAD    (ILU))
ALLOCATE(B%CUR%XSHGC_SH     (ILU))
ALLOCATE(B%CUR%XU_WIN       (ILU))
ALLOCATE(B%CUR%XTRAN_WIN    (ILU))
ALLOCATE(B%CUR%XFLOOR_HEIGHT(ILU))
ALLOCATE(B%CUR%XINF         (ILU))
!
ALLOCATE(B%CUR%XQIN_FLAT    (ILU))
ALLOCATE(B%CUR%XHR_TARGET   (ILU))
ALLOCATE(B%CUR%XV_VENT      (ILU))
ALLOCATE(B%CUR%XCAP_SYS_HEAT(ILU))
ALLOCATE(B%CUR%XCAP_SYS_RAT (ILU))
ALLOCATE(B%CUR%XT_ADP       (ILU))
ALLOCATE(B%CUR%XM_SYS_RAT   (ILU))
ALLOCATE(B%CUR%XCOP_RAT     (ILU))
ALLOCATE(B%CUR%XT_SIZE_MAX  (ILU))
ALLOCATE(B%CUR%XT_SIZE_MIN  (ILU))
ALLOCATE(B%CUR%XF_WATER_COND(ILU))
ALLOCATE(B%CUR%CNATVENT     (ILU))
ALLOCATE(B%CUR%XNATVENT     (ILU))
!
ALLOCATE(B%CUR%XABS_WIN (ILU))
ALLOCATE(B%CUR%XUGG_WIN (ILU))
ALLOCATE(B%CUR%LSHADE   (ILU))
ALLOCATE(B%CUR%XSHADE   (ILU))
ALLOCATE(B%CUR%LSHAD_DAY(ILU))
ALLOCATE(B%CUR%LNATVENT_NIGHT(ILU))
ALLOCATE(B%CUR%XAUX_MAX    (ILU))
ALLOCATE(B%CUR%XN_FLOOR(ILU))
ALLOCATE(B%CUR%XGLAZ_O_BLD(ILU))
ALLOCATE(B%CUR%XMASS_O_BLD(ILU))
ALLOCATE(B%CUR%XFLOOR_HW_RATIO(ILU))
ALLOCATE(B%CUR%XF_FLOOR_MASS(ILU))
ALLOCATE(B%CUR%XF_FLOOR_WALL(ILU))
ALLOCATE(B%CUR%XF_FLOOR_WIN(ILU))
ALLOCATE(B%CUR%XF_FLOOR_ROOF(ILU))
ALLOCATE(B%CUR%XF_WALL_FLOOR(ILU))
ALLOCATE(B%CUR%XF_WALL_MASS(ILU))
ALLOCATE(B%CUR%XF_WALL_WIN(ILU))
ALLOCATE(B%CUR%XF_WIN_FLOOR(ILU))
ALLOCATE(B%CUR%XF_WIN_MASS(ILU))
ALLOCATE(B%CUR%XF_WIN_WALL(ILU))
ALLOCATE(B%CUR%XF_WIN_WIN(ILU))
ALLOCATE(B%CUR%XF_MASS_FLOOR(ILU))
ALLOCATE(B%CUR%XF_MASS_WALL(ILU))
ALLOCATE(B%CUR%XF_MASS_WIN(ILU))

SELECT CASE(TOP%CBEM)
!----------
CASE("DEF")
!-----------
   !parameters that needs to be 0 for calculation
   B%CUR%XGR  (:)         = 0.
   B%CUR%XF_WASTE_CAN(:)  = 0.
!----------
CASE("BEM")
!----------

  B%CUR%XAUX_MAX(:) = 5.
  CALL CONVERT_PATCH_TEB(BDD, DTB, DTCO, DTT, TOP, &
                         TOP%XCOVER,TOP%LCOVER,0.,                                            &
                      PHC_FLOOR=B%CUR%XHC_FLOOR, PTC_FLOOR=B%CUR%XTC_FLOOR, PD_FLOOR=B%CUR%XD_FLOOR,    &
                      PTCOOL_TARGET=B%CUR%XTCOOL_TARGET, PTHEAT_TARGET=B%CUR%XTHEAT_TARGET,       &
                      PF_WASTE_CAN=B%CUR%XF_WASTE_CAN, PEFF_HEAT=B%CUR%XEFF_HEAT, PQIN=B%CUR%XQIN,      &
                      PQIN_FRAD=B%CUR%XQIN_FRAD, PSHGC=B%CUR%XSHGC, PU_WIN=B%CUR%XU_WIN, PGR=B%CUR%XGR,       &
                      PSHGC_SH=B%CUR%XSHGC_SH, PFLOOR_HEIGHT=B%CUR%XFLOOR_HEIGHT, PINF=B%CUR%XINF,      &
                      PF_WATER_COND=B%CUR%XF_WATER_COND, PQIN_FLAT=B%CUR%XQIN_FLAT,               &
                      PHR_TARGET=B%CUR%XHR_TARGET, PV_VENT=B%CUR%XV_VENT,                         &
                      PCAP_SYS_HEAT=B%CUR%XCAP_SYS_HEAT, PCAP_SYS_RAT=B%CUR%XCAP_SYS_RAT,         &
                      PT_ADP=B%CUR%XT_ADP, PM_SYS_RAT=B%CUR%XM_SYS_RAT, PCOP_RAT=B%CUR%XCOP_RAT,        &
                      PT_SIZE_MAX=B%CUR%XT_SIZE_MAX, PT_SIZE_MIN=B%CUR%XT_SIZE_MIN,               &
                      PSHADE=B%CUR%XSHADE, PNATVENT=B%CUR%XNATVENT)
   !
   !
   ! *.     indoor relative surf. and view factors
   !        --------------------------------------
   !
   CALL BEM_MORPHO(T%CUR%XBLD, T%CUR%XWALL_O_HOR, T%CUR%XBLD_HEIGHT, B%CUR%XFLOOR_HEIGHT, B%CUR%XGR,         &
                   B%CUR%XN_FLOOR, T%CUR%XWALL_O_BLD, B%CUR%XGLAZ_O_BLD, B%CUR%XMASS_O_BLD,            &
                   B%CUR%XFLOOR_HW_RATIO, B%CUR%XF_FLOOR_MASS, B%CUR%XF_FLOOR_WALL, B%CUR%XF_FLOOR_WIN,&
                   B%CUR%XF_FLOOR_ROOF, B%CUR%XF_WALL_FLOOR, B%CUR%XF_WALL_MASS,                 &
                   B%CUR%XF_WALL_WIN, B%CUR%XF_WIN_FLOOR, B%CUR%XF_WIN_MASS, B%CUR%XF_WIN_WALL,        &
                   B%CUR%XF_MASS_FLOOR, B%CUR%XF_MASS_WALL, B%CUR%XF_MASS_WIN, B%CUR%XF_WASTE_CAN,     &
                   B%CUR%XF_WIN_WIN      )
   !
   ! *.     Window optical and thermal data
   !        -------------------------------
   !
   CALL WINDOW_DATA(ILU, B%CUR%XSHGC, B%CUR%XU_WIN, B%CUR%XALB_WIN, B%CUR%XABS_WIN, B%CUR%XUGG_WIN, B%CUR%XTRAN_WIN)
   GPRINT = .FALSE.
   DO JJ=1,SIZE(B%CUR%XSHADE)
      IF (B%CUR%XSHADE(JJ) >= 0.0 .AND. B%CUR%XSHADE(JJ) < 0.5) THEN
         B%CUR%LSHADE(JJ) = .FALSE.
      ELSEIF (B%CUR%XSHADE(JJ) >= 0.5 .AND. B%CUR%XSHADE(JJ) <= 1.0) THEN
         B%CUR%LSHADE(JJ) = .TRUE.
      ELSE
       GPRINT = .TRUE.
       B%CUR%LSHADE(JJ) = .FALSE.
      ENDIF
   ENDDO
   IF (GPRINT) WRITE(KLUOUT,*) &
   'TEB-BEM : Error in specifying shading devices for at least one point, no shading device for these points'
   B%CUR%LSHAD_DAY(:) = .FALSE.
   !
   ! *.     Nocturnal surventilation
   !        ------------------------
   GPRINT = .FALSE.
   DO JJ=1,SIZE(B%CUR%XNATVENT)
      IF (B%CUR%XNATVENT(JJ) >= 0.0 .AND. B%CUR%XNATVENT(JJ) < 0.5) THEN
        B%CUR%CNATVENT(JJ) = 'NONE'
      ELSEIF (B%CUR%XNATVENT(JJ) >= 0.5 .AND. B%CUR%XNATVENT(JJ) < 1.5) THEN
        B%CUR%CNATVENT(JJ) = 'MANU'
      ELSEIF (B%CUR%XNATVENT(JJ) >= 1.5 .AND. B%CUR%XNATVENT(JJ) <= 2.5) THEN
        B%CUR%CNATVENT(JJ) = 'AUTO'        
      ELSEIF (B%CUR%XNATVENT(JJ) >= 2.5 .AND. B%CUR%XNATVENT(JJ) <= 3.5) THEN
        B%CUR%CNATVENT(JJ) = 'MECH'        
      ELSE
        GPRINT = .TRUE.
        B%CUR%CNATVENT(JJ) = 'NONE'        
      ENDIF
    ENDDO
    IF (GPRINT) WRITE(KLUOUT,*) 'TEB-BEM : Chosen option for surventilation is not yet implemented; None venting is kept instead'

   B%CUR%LNATVENT_NIGHT(:) = .FALSE.
   !
END SELECT
!
!-------------------------------------------------------------------------------
!
!*       8.     Building HVAC automatic sizing:
!               -------------------------------  
IF (TOP%CBEM=='BEM' .AND. BOP%LAUTOSIZE) THEN
  CALL HVAC_AUTOSIZE(BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, DTT, &
                           DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, &
                           DGUT, DGW, F, FSB, GB, ICP, I, O, S, SSB, SV, &
                           TCP, TGD, TGDO, TGR, TGRO, TVG, W, WSB, &
                     IOB, &
                     B, BOP, UG, U, TG, T, TOP, &
                     ILU,KLUOUT)
  !* stores the real systems characteristics in physiographic data 
  !  for further use
  CALL STORES_HVAC_AUTOSIZE(B, BOP, DTB)
ENDIF
!
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('INIT_BEM_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE INIT_BEM_n
