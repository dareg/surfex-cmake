!     #############################################################
      SUBROUTINE INIT_BEM_n(KLUOUT)
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
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS 
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_BEM_n, ONLY : B => BEM
!
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
ALLOCATE(B%XHC_FLOOR    (ILU,BOP%NFLOOR_LAYER))
ALLOCATE(B%XTC_FLOOR    (ILU,BOP%NFLOOR_LAYER))
ALLOCATE(B%XD_FLOOR     (ILU,BOP%NFLOOR_LAYER))
!
ALLOCATE(B%XTCOOL_TARGET(ILU))
ALLOCATE(B%XTHEAT_TARGET(ILU))
ALLOCATE(B%XEFF_HEAT    (ILU))
ALLOCATE(B%XSHGC        (ILU))
ALLOCATE(B%XQIN         (ILU))
ALLOCATE(B%XQIN_FRAD    (ILU))
ALLOCATE(B%XSHGC_SH     (ILU))
ALLOCATE(B%XU_WIN       (ILU))
ALLOCATE(B%XTRAN_WIN    (ILU))
ALLOCATE(B%XFLOOR_HEIGHT(ILU))
ALLOCATE(B%XINF         (ILU))
!
ALLOCATE(B%XQIN_FLAT    (ILU))
ALLOCATE(B%XHR_TARGET   (ILU))
ALLOCATE(B%XV_VENT      (ILU))
ALLOCATE(B%XCAP_SYS_HEAT(ILU))
ALLOCATE(B%XCAP_SYS_RAT (ILU))
ALLOCATE(B%XT_ADP       (ILU))
ALLOCATE(B%XM_SYS_RAT   (ILU))
ALLOCATE(B%XCOP_RAT     (ILU))
ALLOCATE(B%XT_SIZE_MAX  (ILU))
ALLOCATE(B%XT_SIZE_MIN  (ILU))
ALLOCATE(B%XF_WATER_COND(ILU))
ALLOCATE(B%CNATVENT     (ILU))
ALLOCATE(B%XNATVENT     (ILU))
!
ALLOCATE(B%XABS_WIN (ILU))
ALLOCATE(B%XUGG_WIN (ILU))
ALLOCATE(B%LSHADE   (ILU))
ALLOCATE(B%XSHADE   (ILU))
ALLOCATE(B%LSHAD_DAY(ILU))
ALLOCATE(B%LNATVENT_NIGHT(ILU))
ALLOCATE(B%XAUX_MAX    (ILU))
ALLOCATE(B%XN_FLOOR(ILU))
ALLOCATE(B%XGLAZ_O_BLD(ILU))
ALLOCATE(B%XMASS_O_BLD(ILU))
ALLOCATE(B%XFLOOR_HW_RATIO(ILU))
ALLOCATE(B%XF_FLOOR_MASS(ILU))
ALLOCATE(B%XF_FLOOR_WALL(ILU))
ALLOCATE(B%XF_FLOOR_WIN(ILU))
ALLOCATE(B%XF_FLOOR_ROOF(ILU))
ALLOCATE(B%XF_WALL_FLOOR(ILU))
ALLOCATE(B%XF_WALL_MASS(ILU))
ALLOCATE(B%XF_WALL_WIN(ILU))
ALLOCATE(B%XF_WIN_FLOOR(ILU))
ALLOCATE(B%XF_WIN_MASS(ILU))
ALLOCATE(B%XF_WIN_WALL(ILU))
ALLOCATE(B%XF_WIN_WIN(ILU))
ALLOCATE(B%XF_MASS_FLOOR(ILU))
ALLOCATE(B%XF_MASS_WALL(ILU))
ALLOCATE(B%XF_MASS_WIN(ILU))

SELECT CASE(TOP%CBEM)
!----------
CASE("DEF")
!-----------
   !parameters that needs to be 0 for calculation
   B%XGR  (:)         = 0.
   B%XF_WASTE_CAN(:)  = 0.
!----------
CASE("BEM")
!----------

  B%XAUX_MAX(:) = 5.
  CALL CONVERT_PATCH_TEB(TOP%XCOVER,TOP%LCOVER,0.,                                            &
                      PHC_FLOOR=B%XHC_FLOOR, PTC_FLOOR=B%XTC_FLOOR, PD_FLOOR=B%XD_FLOOR,    &
                      PTCOOL_TARGET=B%XTCOOL_TARGET, PTHEAT_TARGET=B%XTHEAT_TARGET,       &
                      PF_WASTE_CAN=B%XF_WASTE_CAN, PEFF_HEAT=B%XEFF_HEAT, PQIN=B%XQIN,      &
                      PQIN_FRAD=B%XQIN_FRAD, PSHGC=B%XSHGC, PU_WIN=B%XU_WIN, PGR=B%XGR,       &
                      PSHGC_SH=B%XSHGC_SH, PFLOOR_HEIGHT=B%XFLOOR_HEIGHT, PINF=B%XINF,      &
                      PF_WATER_COND=B%XF_WATER_COND, PQIN_FLAT=B%XQIN_FLAT,               &
                      PHR_TARGET=B%XHR_TARGET, PV_VENT=B%XV_VENT,                         &
                      PCAP_SYS_HEAT=B%XCAP_SYS_HEAT, PCAP_SYS_RAT=B%XCAP_SYS_RAT,         &
                      PT_ADP=B%XT_ADP, PM_SYS_RAT=B%XM_SYS_RAT, PCOP_RAT=B%XCOP_RAT,        &
                      PT_SIZE_MAX=B%XT_SIZE_MAX, PT_SIZE_MIN=B%XT_SIZE_MIN,               &
                      PSHADE=B%XSHADE, PNATVENT=B%XNATVENT)
   !
   !
   ! *.     indoor relative surf. and view factors
   !        --------------------------------------
   !
   CALL BEM_MORPHO(T%XBLD, T%XWALL_O_HOR, T%XBLD_HEIGHT, B%XFLOOR_HEIGHT, B%XGR,         &
                   B%XN_FLOOR, T%XWALL_O_BLD, B%XGLAZ_O_BLD, B%XMASS_O_BLD,            &
                   B%XFLOOR_HW_RATIO, B%XF_FLOOR_MASS, B%XF_FLOOR_WALL, B%XF_FLOOR_WIN,&
                   B%XF_FLOOR_ROOF, B%XF_WALL_FLOOR, B%XF_WALL_MASS,                 &
                   B%XF_WALL_WIN, B%XF_WIN_FLOOR, B%XF_WIN_MASS, B%XF_WIN_WALL,        &
                   B%XF_MASS_FLOOR, B%XF_MASS_WALL, B%XF_MASS_WIN, B%XF_WASTE_CAN,     &
                   B%XF_WIN_WIN      )
   !
   ! *.     Window optical and thermal data
   !        -------------------------------
   !
   CALL WINDOW_DATA(ILU, B%XSHGC, B%XU_WIN, B%XALB_WIN, B%XABS_WIN, B%XUGG_WIN, B%XTRAN_WIN)
   GPRINT = .FALSE.
   DO JJ=1,SIZE(B%XSHADE)
      IF (B%XSHADE(JJ) >= 0.0 .AND. B%XSHADE(JJ) < 0.5) THEN
         B%LSHADE(JJ) = .FALSE.
      ELSEIF (B%XSHADE(JJ) >= 0.5 .AND. B%XSHADE(JJ) <= 1.0) THEN
         B%LSHADE(JJ) = .TRUE.
      ELSE
       GPRINT = .TRUE.
       B%LSHADE(JJ) = .FALSE.
      ENDIF
   ENDDO
   IF (GPRINT) WRITE(KLUOUT,*) &
   'TEB-BEM : Error in specifying shading devices for at least one point, no shading device for these points'
   B%LSHAD_DAY(:) = .FALSE.
   !
   ! *.     Nocturnal surventilation
   !        ------------------------
   GPRINT = .FALSE.
   DO JJ=1,SIZE(B%XNATVENT)
      IF (B%XNATVENT(JJ) >= 0.0 .AND. B%XNATVENT(JJ) < 0.5) THEN
        B%CNATVENT(JJ) = 'NONE'
      ELSEIF (B%XNATVENT(JJ) >= 0.5 .AND. B%XNATVENT(JJ) < 1.5) THEN
        B%CNATVENT(JJ) = 'MANU'
      ELSEIF (B%XNATVENT(JJ) >= 1.5 .AND. B%XNATVENT(JJ) <= 2.5) THEN
        B%CNATVENT(JJ) = 'AUTO'        
      ELSEIF (B%XNATVENT(JJ) >= 2.5 .AND. B%XNATVENT(JJ) <= 3.5) THEN
        B%CNATVENT(JJ) = 'MECH'        
      ELSE
        GPRINT = .TRUE.
        B%CNATVENT(JJ) = 'NONE'        
      ENDIF
    ENDDO
    IF (GPRINT) WRITE(KLUOUT,*) 'TEB-BEM : Chosen option for surventilation is not yet implemented; None venting is kept instead'

   B%LNATVENT_NIGHT(:) = .FALSE.
   !
END SELECT
!
!-------------------------------------------------------------------------------
!
!*       8.     Building HVAC automatic sizing:
!               -------------------------------  
IF (TOP%CBEM=='BEM' .AND. BOP%LAUTOSIZE) THEN
  CALL HVAC_AUTOSIZE(ILU,KLUOUT)
  !* stores the real systems characteristics in physiographic data 
  !  for further use
  CALL STORES_HVAC_AUTOSIZE
ENDIF
!
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('INIT_BEM_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE INIT_BEM_n
