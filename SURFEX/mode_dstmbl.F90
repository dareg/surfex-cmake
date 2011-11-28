!     #########
MODULE MODE_DSTMBL

  !PURPOSE: 
  !Take input from ISBA land surface model and
  !calculate a dust flux which is consistent with the input.

  !THEORY:
  !Based on Marticorena/Bergametti, 1995 and Zender et al 2003 (JGR)

  !CODE HISTORY
  !Code is a modified version of dstmbl.F90 in the DEAD model
  !Original version was downloaded from the DEAD homepage
  !http://dust.ess.uci.edu/dead/ on January 10th 2005

  !AUTHOR (or rather "code modifyer")
  !Alf Grini <alf.grini@cnrm.meteo.fr>

    USE MODD_DST_SURF, ONLY :  XFLX_MSS_FDG_FCT
!
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
  USE PARKIND1  ,ONLY : JPRB
!
  implicit none
  public

contains

  subroutine dustflux_get(          &
         PUSTAR,                     &!I [m/s] Wind friction speed 
         PRHOA,                     &!I [kg/m3] air density at 2m height 
         PWG,                       &!I [m3/m3] volumetric water content 
         PZ0,                       &!I [m] roughness length of surface
         PWSAT,                     &!I [m3 m-3] saturation liquid water content
         PCLAY,                     &!I [frc] mass fraction clay
         PSAND,                     &!I [frc] mass fraction sand
         PWIND10M,                  &!I [m/s] wind at 10m altitude
         PSFDST,                    &!O [kg/m2/sec] Vertical dust flux
         KSIZE                     &!I [nbr] number of points for calculation
         )  

    USE MODE_DSTMBLUTL                     !Dust mobilization subroutines
    USE MODD_DST_SURF, ONLY :  XFLX_MSS_FDG_FCT

    implicit none
    
    !INPUT, set their dimensions to their passed lengths or to KSIZE ?
    integer, intent(in)                  :: KSIZE    ![nbr] length of passed arrays
    real, intent(in), dimension(KSIZE)   :: PUSTAR   ![m/s] wind friction speed
    real, intent(in), dimension(KSIZE)   :: PRHOA    ![kg/m3] air density
    real, intent(in), dimension(KSIZE)   :: PCLAY    ![frc] mass fraction clay
    real, intent(in), dimension(KSIZE)   :: PSAND    ![frc] mass fraction sand
    real, intent(in), dimension(KSIZE)   :: PWG      ![m3 m-3] volumetric water fraction
    real, intent(in), dimension(KSIZE)   :: PWSAT    ![m3 m-3] saturation water content
    real, intent(in), dimension(KSIZE)   :: PZ0      ![m] surface roughness length
    real, intent(in), dimension(KSIZE)   :: PWIND10M ![m/s] wind at 10m altitude

    !OUTPUT the flux of dust
    real, intent(out), dimension(KSIZE)  :: PSFDST   ! [kg m-2 s-1] Output flux of atmospheric dust

!!!!!!!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!!!!!!

!#ifdef AlG01
!    real,parameter::XFLX_MSS_FDG_FCT=28. ! [frc] Global mass flux tuning factor (a posteriori)
!#else
!    real,parameter::XFLX_MSS_FDG_FCT=7.0e-4 ! [frc] Global mass flux tuning factor (a posteriori)
!    real,parameter::XFLX_MSS_FDG_FCT=21.0e-4 ! [frc] Global mass flux tuning factor (a posteriori)
    !real,parameter::XFLX_MSS_FDG_FCT=12.0e-4 ! [frc] values used in Masdev47
!    real,parameter::flx_mss_fdg_fctm=4.0e-4 ! [frc] Global mass flux tuning factor (a posteriori) (proposez by Pierre)

!#endif
    real,parameter::hgt_rfr=10.0 ! [m] Reference height for mobilization processes
    real,parameter::hgt_zpd_mbl=0.0 ! [m] Zero plane displacement for erodible surfaces
    real,parameter::rgh_mmn_mbl=100.0e-6 ! [m] Roughness length momentum for erodible surfaces MaB95 p. 16420, GMB98 p. 6205
    ! fxm: rgh_mmn_smt set to 33.3e-6 um, MaB95 p. 16426 recommend 10.0e-6
    real,parameter::rgh_mmn_smt=33.3e-6 ! [m] Smooth roughness length MaB95 p. 16426, MaB97 p. 4392, GMB98 p. 6207
    real,parameter::wnd_min_mbl=1.0 ! [m s-1] Minimum windspeed used for mobilization 
    real,parameter::wnd_frc_rsl=0.95d0 ! [frc] Fraction of wind PDF to resolve 

    !Define local variables:
    logical :: flg_mbl(KSIZE)          ![frc] Mobilization candidate flag
    real  :: frc_thr_ncr_drg(KSIZE)  ![frc] fraction by which drag partitioning increases threshold wind
    real  :: frc_thr_ncr_wtr(KSIZE)  ![frc] Fraction by which soil wetness increases threshold wind
    real  :: gwc_sfc(KSIZE)          ![kg/kg] Gravimetric water content
    real  :: wnd_frc_thr_slt(KSIZE)  ![m/s] Threshold wind friction speed when all effects taken into account
    real  :: wnd_frc_slt(KSIZE)         ![m/s] wind friction speed after modified for saltation feedbacks
    real  :: flx_mss_hrz_slt_ttl_wbn(KSIZE) ![kg m-1 s-1] Vertically integrated horizontal saltation soil flux for a wind bin 
    real  :: flx_mss_vrt_dst_ttl_wbn(KSIZE)     ![kg m-2 s-1]
    real  :: wnd_rfr_thr_slt(KSIZE)             ![m s-1] Threshold wind speed at reference level
    real  :: mbl_bsn_fct(KSIZE)                 ![frc] enhancement factor for grid cells with higher erodibility
    real  :: dst_slt_flx_rat_ttl(KSIZE)         ![m-1] ratio of vertical to horizontal flux (alpha in several papers)
    real  :: ZCLAY(KSIZE)                       ![frc] dummy for fraction of clay
    real  :: wnd_rfr(KSIZE)                     ![m s-1] wind speed at reference level    

    integer             :: i                   !Counter for number of points (used in loops)
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    
    !Allocate the local variables
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBL:DUSTFLUX_GET',0,ZHOOK_HANDLE)

    !Set fraction of clay to 20, this does not work anyway
    ZCLAY(:)=0.2

    !Initialize mobilization candidate flag
    flg_mbl(:)=.TRUE.
    
    !Initialize vertical dust flux
    flx_mss_vrt_dst_ttl_wbn(:)=0.d0

    !fxm: Get erodibility limitation factor, use something connected to amount of sand
    !Discuss with Valery Masson
    mbl_bsn_fct(:)=PSAND(:)
!
!!! utilisé dans le calcul de l'effet Owen 
    wnd_rfr(:) = PWIND10M(:)                    ! proposed
!
    ! Factor by which surface roughness increases threshold friction velocity 
    !++grini: fxm: USE WHOLE ARRAY OF Z0 INSTEAD OF ONLY RGH_MMN_MBL AS IN OLD CODE
    call frc_thr_ncr_drg_get( &
           frc_thr_ncr_drg,  &! O [frc] Factor by which surface roughness increases threshold friction velocity
           PZ0,             &! I [m] Roughness length momentum for this (erodible) surface
           rgh_mmn_smt,   &!I [m] Smooth roughness length
           KSIZE,          &!I [nbr] loop size
           KSIZE )         !I [nbr] dimension size  

    ! Convert volumetric water content to gravimetric water content
    call vwc2gwc( &
           flg_mbl,  &! I [flg] Mobilization candidate flag
           gwc_sfc,  &! O [kg kg-1] Gravimetric water content
           PWSAT,    &! I [m3 m-3] Saturated volumetric water content (sand-dependent)
           PWG,      &! I [m3 m-3] Volumetric water content
           KSIZE,    &! I [nbr] loop size
           KSIZE)     ! I [nbr] array dimension  
!
    ! Factor by which soil moisture increases threshold friction velocity 
    call frc_thr_ncr_wtr_get( &
           flg_mbl,  &! I [flg] Mobilization candidate flag
           frc_thr_ncr_wtr,  &! O [frc] Factor by which moisture increases threshold friction velocity
           PCLAY,  &! I [kg kg-1] Mass fraction of clay
           gwc_sfc, &! I [kg kg-1] Gravimetric water content
           KSIZE ,  &!I [nbr] loop size
           KSIZE )   !I [nbr] dimension size  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! fxm: Use surface density not midlayer density
    call wnd_frc_thr_slt_get(                      &!++grini: Removed variables which were not used
           PRHOA,    &! I [kg m-3] Midlayer density
           wnd_frc_thr_slt,      &! O [m s-1] Threshold friction velocity for saltation
           KSIZE,                &! I [nbr] loop dimension
           KSIZE )                ! I [nbr] Array dimension  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i=1,KSIZE
       wnd_frc_thr_slt(i)=  &! [m s-1] Threshold friction velocity for saltation
              wnd_frc_thr_slt(i)*  &! [m s-1] Threshold for dry, flat ground
              frc_thr_ncr_wtr(i)  !*  ! [frc] Adjustment for moisture 
       !frc_thr_ncr_drg(i) ! [frc] Adjustment for roughness
     end do                    ! end loop over elements 

    ! Threshold saltation wind speed
    !Needed for the saltation feedback roughening effect
    !do i=1,KSIZE
    !   if (flg_mbl(i)) then
    !      wnd_rfr_thr_slt(i)= & ! [m s-1] Threshold 10 m wind speed for saltation
    !           wnd_rfr(i)*wnd_frc_thr_slt(i)/PUSTAR(i) !++alfgr 
    !   endif                  ! endif flg_mbl
    !end do                    ! end loop over lon


    !CHECK IF THIS CAN BE USED EASILY
    !NEEDS 10M WIND SPEED WHICH IS MAYBE KNOWN MAYBE NOT !
    ! Saltation increases friction speed by roughening surface
    !call wnd_frc_slt_get( &
    !    flg_mbl, & ! I [flg] Mobilization candidate flag
    !     PUSTAR, & ! I [m s-1] Surface friction velocity
    !     wnd_frc_slt, & ! O [m s-1] Saltating friction velocity
    !     wnd_rfr, & ! I [m s-1] Wind speed at reference height
    !     wnd_rfr_thr_slt) ! I [m s-1] Threshold 10 m wind speed for saltation

    !Skip the roughening of surface effect for now, and 
    !just use the wind friction speed as it is modified
    !by drag partitioning
    wnd_frc_slt(:)=PUSTAR(:)/frc_thr_ncr_drg(:) 

    ! Horizontal streamwise mass flux for old "bulk" formulation
    call flx_mss_hrz_slt_ttl_Whi79_get( &
           PRHOA,  &! I [kg m-3] Midlayer density
           flg_mbl,  &! I [flg] Mobilization candidate flag
           flx_mss_hrz_slt_ttl_wbn,  &! O [kg m-1 s-1] Vertically integrated streamwise mass flux in wind bin
           wnd_frc_slt,  &! I [m s-1] Saltating friction velocity
           wnd_frc_thr_slt,  &! I [m s-1] Threshold friction speed for saltation
           KSIZE,            &! I [nbr] loop size
           KSIZE )            ! I [nbr] array dimension size  


    ! Apply land surface and vegetation limitations and global tuning factor
    do i=1,KSIZE
       ! Vertically integrated streamwise mass flux in wind bin
       flx_mss_hrz_slt_ttl_wbn(i)=flx_mss_hrz_slt_ttl_wbn(i)  &! [kg m-2 s-1]
            !*lnd_frc_mbl(i) & ! [frc] Bare ground fraction
              *mbl_bsn_fct(i)  &! [frc] Erodibility factor
              *XFLX_MSS_FDG_FCT  ! [frc] Global mass flux tuning factor (empirical)  
    end do ! end loop over lon
    
    ! Vertical dust mass flux
    call flx_mss_vrt_dst_ttl_MaB95_get( &
           dst_slt_flx_rat_ttl,  &! O [m-1] Ratio of vertical dust flux to streamwise mass flux
           flg_mbl,  &! I [flg] Mobilization candidate flag
           flx_mss_hrz_slt_ttl_wbn,  &! I [kg m-1 s-1] Vertically integrated streamwise mass flux
           flx_mss_vrt_dst_ttl_wbn,  &! O [kg m-2 s-1] Total vertical mass flux of dust
           ZCLAY, &! I [kg kg-1] Mass fraction clay
           KSIZE,   &!I [nbr] loop size
           KSIZE )   !I [nbr] array dimension size  
!
    !Assign the output vertical dust flux to the value calculated
    !PSFDST(:) = flx_mss_vrt_dst_ttl_wbn(:)

    PSFDST(:) = dst_slt_flx_rat_ttl(:)  * flx_mss_hrz_slt_ttl_wbn(:)

  IF (LHOOK) CALL DR_HOOK('MODE_DSTMBL:DUSTFLUX_GET',1,ZHOOK_HANDLE)

  end subroutine dustflux_get

END MODULE MODE_DSTMBL
