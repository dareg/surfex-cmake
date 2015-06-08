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
USE MODD_DATA_BEM_n, ONLY : DTB => DATA_BEM
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
  CALL CONVERT_PATCH_TEB(TOP%XCOVER,TOP%LCOVER,0.,                                            &
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
  CALL HVAC_AUTOSIZE(ILU,KLUOUT)
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
