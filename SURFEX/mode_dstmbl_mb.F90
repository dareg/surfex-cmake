!     #########
MODULE MODE_DSTMBL_MB

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
  !Mohamed MOKHTARI <m.mokhtari@meteo.dz  , m_morad06@yahoo.fr>

    USE MODD_DST_SURF, ONLY :  CVERMOD, XFLX_MSS_FDG_FCT

  implicit none

contains

  subroutine dustflux_get_MB(          &
       PUSTAR                   &  !I [m/s] Wind friction speed 
       ,PRHOA                   &  !I [kg/m3] air density at 2m height 
       ,PWG                     &  !I [m3/m3] volumetric water content 
       ,PZ0                     &  !I [m] roughness length of surface
       ,PWSAT                   &  !I [m3 m-3] saturation liquid water content
       ,PCLAY                   &  !I [frc] mass fraction clay
       ,PSAND                   &  !I [frc] mass fraction sand
       ,PWIND10M                &  !I [m/s] wind at 10m altitude
       ,PSFDST                  &  !O [kg/m2/sec] Vertical dust flux
       ,KSIZE                   &  !I [nbr] number of points for calculation
       )

    USE MODE_DSTMBLUTL                     !Dust mobilization subroutines
    USE MODD_DST_SURF, ONLY :  CVERMOD, XFLX_MSS_FDG_FCT
    USE MODD_CSTS, only: XPI

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
    real,parameter::flx_mss_fdg_fctm=0.04d0  ! [frc] Global mass flux tuning factor (a posteriori)

!#endif
    real,parameter::hgt_rfr=10.0 ! [m] Reference height for mobilization processes
    real,parameter::hgt_zpd_mbl=0.0 ! [m] Zero plane displacement for erodible surfaces
    real,parameter::rgh_mmn_mbl=100.0e-6 ! [m] Roughness length momentum for erodible surfaces MaB95 p. 16420, GMB98 p. 6205
    ! fxm: rgh_mmn_smt set to 33.3e-6 um, MaB95 p. 16426 recommend 10.0e-6
    real,parameter::rgh_mmn_smt=33.3e-6 ! [m] Smooth roughness length MaB95 p. 16426, MaB97 p. 4392, GMB98 p. 6207
    real,parameter::wnd_min_mbl=1.0 ! [m s-1] Minimum windspeed used for mobilization 
    real,parameter::wnd_frc_rsl=0.95d0 ! [frc] Fraction of wind PDF to resolve 
!
    real,parameter::DLNDP=0.1d0        ! [µm]  Dln(DP) 
    integer,parameter::NTEX=12         ! number of texture
    integer,parameter::NMODE=3         ! number of mode
    integer,parameter::NDP=100         ! number of paticle
!
    !Define local variables:
    logical  :: flg_mbl(KSIZE)          ![frc] Mobilization candidate flag
    real  :: frc_thr_ncr_drg(KSIZE)  ![frc] fraction by which drag partitioning increases threshold wind
    real  :: frc_thr_ncr_wtr(KSIZE)  ![frc] Fraction by which soil wetness increases threshold wind
    real  :: gwc_sfc(KSIZE)          ![kg/kg] Gravimetric water content
    real  :: wnd_frc_thr_slt(KSIZE)  ![m/s] Threshold wind friction speed when all effects taken into account
    real  :: wnd_frc_slt(KSIZE)      ![m/s] Threshold wind friction speed when all effects taken into account    
    real  :: flx_mss_hrz_slt_ttl_wbn(KSIZE)     ![kg m-1 s-1] Vertically integrated horizontal saltation soil flux for a wind bin 
    real  :: flx_mss_vrt_dst_ttl_wbn(KSIZE)     ![kg m-2 s-1]
    real  :: mbl_bsn_fct(KSIZE)                 ![frc] enhancement factor for grid cells with higher erodibility
    real  :: dst_slt_flx_rat_ttl(KSIZE)         ![m-1] ratio of vertical to horizontal flux (alpha in several papers)
    real  :: wnd_rfr(KSIZE)                     ![m s-1] wind speed at reference level    
    real  :: wnd_rfr_thr_slt(KSIZE)             ![m s-1] Threshold wind speed at reference level

    real  :: ZSILT(KSIZE)                       ![frc] dummy for fraction of silt 
    real  :: ZWPRM(KSIZE)                       !threshold soil wetness
!    real :: flx_mss_hrz_slt_ttl_tot(KSIZE) ![kg m-1 s-1] Vertically integrated horizontal saltation soil flux for a wind bin 
    real  :: flx_mss_vrt_dst_tot(KSIZE)         ![kg m-2 s-1] Vertically integrated flux for a wind bin
!    
    real  :: ZDSRLV(NTEX,NDP)                 ! Surface relative des grains du sol
    real  :: ZZS0(NTEX)                       ! rugosité de la surface lisse
    real  :: ZDP(NDP)                         ! [µm] diamètre du particule  
    real  :: ZDPLN(NDP)                       ! ln(Dp)
!    
    integer  :: ZTEXT(KSIZE)                    !soil texture
!
    real ZP_DP ! [5m] 

!
    integer             :: i                             !Counter for number of points (used in loops)
    integer             :: JMOD                          !Counter for number of mode
    integer             :: IDP                           !Counter for number of particle
    integer             :: ITEX                          !Counter for number of texture
    
    !Initialize mobilization candidate flag
    flg_mbl(:)=.TRUE.
    
    !Initialize vertical dust flux
    flx_mss_vrt_dst_ttl_wbn(:)=0.d0
!    flx_mss_hrz_slt_ttl_tot(:)= 0.d0
    flx_mss_vrt_dst_tot(:)=0.d0
    !fxm: Get erodibility limitation factor, use something connected to amount of sand
    !Discuss with Valery Masson
    mbl_bsn_fct(:)=PSAND(:)
!
!!! utilisï¿½ dans le calcul de l'effet Owen 
    wnd_rfr(:) = PWIND10M(:)                    ! proposed
  
  do i = 1, KSIZE
    ZSILT(i) = 1.0 - PCLAY(i) - PSAND(i)    
     if (ZSILT(i)<=0.) ZSILT(i) = 0.0
!      
      if (PCLAY(i) < 7E-2) then
          if ((PSAND(i) <=1.).and.(PSAND(i)>=(0.3*PCLAY(i)+0.87))) then
          ZTEXT(i)=1       !Sand
           ZWPRM(i)= 1.5E-2
          elseif ((PSAND(i) <(0.3*PCLAY(i)+0.87)).and.(PSAND(i)>=(PCLAY(i)+0.7))) then
          ZTEXT(i)=2       !Loamy Sand
          ZWPRM(i)= 2.5E-2
          elseif ((PSAND(i) <(PCLAY(i)+0.7)).and.(PSAND(i)>=(0.5 -PCLAY(i)))) then
          ZTEXT(i)=3       !Sandy Loam
          ZWPRM(i)= 3.0E-2
          elseif ((ZSILT(i) <0.8).and.(ZSILT(i)>=0.5)) then
          ZTEXT(i)=4       !Silt Loam
          ZWPRM(i)= 5.E-2
          elseif ((ZSILT(i) <1.).and.(ZSILT(i)>=0.8)) then
          ZTEXT(i)=12      !Silt
          ZWPRM(i)= 2.5E-2
          endif       
      endif
!  
      if ((PCLAY(i) < 10E-2).and.(PCLAY(i)>=7E-2)) then
          if ((PSAND(i) <0.93).and.(PSAND(i)>=(0.3*PCLAY(i)+0.87))) then
          ZTEXT(i)=1       !Sand
          ZWPRM(i)= 1.5E-2
          elseif ((PSAND(i) <(0.3*PCLAY(i)+0.87)).and.(PSAND(i)>=(PCLAY(i)+0.7))) then
          ZTEXT(i)=2       !Loamy Sand
          ZWPRM(i)= 2.5E-2
          elseif ((PSAND(i) <(PCLAY(i)+0.7)).and.(PSAND(i)>=0.52)) then
          ZTEXT(i)=3       !Sandy Loam
          ZWPRM(i)= 3.0E-2
          elseif ((PSAND(i) <0.52).and.(PSAND(i)>=(0.5 -PCLAY(i)))) then
          ZTEXT(i)=5       !Loam
          ZWPRM(i)= 4E-2
          elseif ((ZSILT(i) <0.8).and.(ZSILT(i)>=0.5)) then
          ZTEXT(i)=4       !Silt Loam
          ZWPRM(i)= 5E-2
          elseif ((ZSILT(i) <1.).and.(ZSILT(i)>=0.8)) then
          ZTEXT(i)=12       !Silt
          ZWPRM(i)= 2.5E-2
          endif       
      endif
!
      if ((PCLAY(i) < 12E-2).and.(PCLAY(i)>=10E-2)) then
          if ((PSAND(i) <0.9).and.(PSAND(i)>=(PCLAY(i)+0.7))) then
          ZTEXT(i)=2       !Loamy Sand
          ZWPRM(i)= 2.5E-2
          elseif ((PSAND(i) <(PCLAY(i)+0.7)).and.(PSAND(i)>=0.52)) then
          ZTEXT(i)=3       !Sandy Loam
          ZWPRM(i)= 3.0E-2
          elseif ((PSAND(i) <0.52).and.(PSAND(i)>=0.5 - PCLAY(i))) then
          ZTEXT(i)=5       !Loam
          ZWPRM(i)= 4E-2
          elseif ((ZSILT(i) <0.8).and.(ZSILT(i)>=0.5)) then
          ZTEXT(i)=4       !Silt Loam
          ZWPRM(i)= 5E-2
          elseif ((ZSILT(i) <1.).and.(ZSILT(i)>=0.8)) then
          ZTEXT(i)=12       !Silt
          ZWPRM(i)= 2.5E-2
          endif       
      endif
!      
      if ((PCLAY(i) < 15E-2).and.(PCLAY(i)>=12E-2)) then
          if ((PSAND(i) <0.88).and.(PSAND(i)>=(PCLAY(i)+0.7))) then
          ZTEXT(i)=2       !Loamy Sand
          ZWPRM(i)= 2.5E-2
          elseif ((PSAND(i) <(PCLAY(i)+0.7)).and.(PSAND(i)>=0.52)) then
          ZTEXT(i)=3       !Sandy Loam
          ZWPRM(i)= 3.0E-2
          elseif ((PSAND(i) <0.52).and.(PSAND(i)>=(0.5 - PCLAY(i)))) then
          ZTEXT(i)=5       !Loam
          ZWPRM(i)= 4E-2
          elseif ((ZSILT(i) <0.88).and.(ZSILT(i)>=0.5)) then
          ZTEXT(i)=4       !Silt Loam
          ZWPRM(i)= 5E-2
          endif       
      endif
!
      if ((PCLAY(i) < 20E-2).and.(PCLAY(i)>=15E-2)) then
          if ((PSAND(i) <0.85).and.(PSAND(i)>=0.52)) then
          ZTEXT(i)=3       !Sandy Loam
          ZWPRM(i)= 3.0E-2
          elseif ((PSAND(i) <0.52).and.(PSAND(i)>=(0.5 - PCLAY(i)))) then
          ZTEXT(i)=5       !Loam
          ZWPRM(i)= 4E-2
          elseif ((ZSILT(i) <0.85).and.(ZSILT(i)>=0.5)) then
          ZTEXT(i)=4       !Silt Loam
          ZWPRM(i)= 5E-2
          endif       
      endif
!
      if ((PCLAY(i) < 28E-2).and.(PCLAY(i)>=20E-2)) then
          if ((ZSILT(i) <0.28).and.(ZSILT(i)>=0.)) then
          ZTEXT(i)=6       !Sandy Clay Loam
          ZWPRM(i)= 6.E-2
          elseif ((ZSILT(i) <0.5).and.(ZSILT(i)>=0.28)) then
          ZTEXT(i)=5       !Loam
          ZWPRM(i)= 4E-2
          elseif ((ZSILT(i) <0.8).and.(ZSILT(i)>=0.5)) then
          ZTEXT(i)=4       !Silt Loam
          ZWPRM(i)= 5E-2
          endif       
      endif
!
      if ((PCLAY(i) < 36E-2).and.(PCLAY(i)>=28E-2)) then
          if ((PSAND(i) <0.72).and.(PSAND(i)>=0.45)) then
          ZTEXT(i)=6       !Sandy Clay Loam
          ZWPRM(i)= 6.E-2
          elseif ((PSAND(i) <0.45).and.(PSAND(i)>=0.2)) then
          ZTEXT(i)=8       !Clay Loam
          ZWPRM(i)= 6.8E-2
          elseif ((PSAND(i) <0.2).and.(PSAND(i)>=0.)) then
          ZTEXT(i)=7       !Silty Clay Loam
          ZWPRM(i)= 6.8E-2
          endif       
      endif
!      
      if ((PCLAY(i) < 40E-2).and.(PCLAY(i)>=36E-2)) then
          if ((PSAND(i) <0.64).and.(PSAND(i)>=0.45)) then
          ZTEXT(i)=9       !Sandy Clay
          ZWPRM(i)= 10E-2
          elseif ((PSAND(i) <0.45).and.(PSAND(i)>=0.2)) then
          ZTEXT(i)=8       !Clay Loam
          ZWPRM(i)= 6.8E-2
          elseif ((PSAND(i) <0.2).and.(PSAND(i)>=0.)) then
          ZTEXT(i)=7       !Silty Clay Loam
          ZWPRM(i)= 6.8E-2
          endif       
      endif
!
      if ((PCLAY(i) < 55E-2).and.(PCLAY(i)>=40E-2)) then
          if ((PSAND(i) <0.60).and.(PSAND(i)>=0.45)) then
          ZTEXT(i)=9       !Sandy Clay
          ZWPRM(i)= 10.E-2
          elseif ((PSAND(i) <0.45).and.(PSAND(i)>=0.05).and.(ZSILT(i) <0.40).and.(ZSILT(i)>=0.)) then
          ZTEXT(i)=11       !Clay
          ZWPRM(i)= 11.5E-2
          elseif ((ZSILT(i) <0.60).and.(ZSILT(i)>=0.40)) then
          ZTEXT(i)=10       !Silty Clay
          ZWPRM(i)= 10.5E-2
          endif       
      endif
!
      if ((PCLAY(i) < 60E-2).and.(PCLAY(i)>=55E-2)) then
          if ((ZSILT(i) <0.40).and.(ZSILT(i)>=0.)) then
          ZTEXT(i)=11       !Clay
          ZWPRM(i)= 11.5E-2
          elseif ((ZSILT(i) <0.45).and.(ZSILT(i)>=0.40)) then
          ZTEXT(i)=10       !Silty Clay
          ZWPRM(i)= 10.5E-2
          endif       
      endif
!
      if ((PCLAY(i) < 100E-2).and.(PCLAY(i)>=60E-2)) then
          if ((ZSILT(i) <0.40).and.(ZSILT(i)>=0.)) then
          ZTEXT(i)=11       !Clay
          ZWPRM(i)= 11.5E-2
          endif       
      endif
!
  enddo  
!&&&&&&&&&&&&&&&&&&&
   call distribution(         &
        KSIZE                 &  !I loop over lon 
       ,NDP                   &  !I loop over diameter
       ,NTEX                  &  !I loop over texture 
       ,NMODE                 &  !I loop over mode 
       ,DLNDP                 &  !I [µm] dln(Dp)
       ,ZTEXT                 &  !I texture
       ,ZDSRLV                &  !O surface relative occupée par chaque particule
       ,ZZS0)                    !O rugosité de la surface lisse  
!
    ! Factor by which surface roughness increases threshold friction velocity 
    !++grini: fxm: USE WHOLE ARRAY OF Z0 INSTEAD OF ONLY RGH_MMN_MBL AS IN OLD CODE
    call frc_thr_ncr_drg_get_mm( &
         frc_thr_ncr_drg,        & ! O [frc] Factor by which surface roughness increases threshold friction velocity
         PZ0,                    & ! I [m] Roughness length momentum for this (erodible) surface
         ZZS0,                   & !I [m] Smooth roughness length
         ZTEXT,                  & !I Soil texture
         KSIZE,                  & !I [nbr] loop size
         KSIZE,                  & !I [nbr] dimension size
	 NTEX)                  !  I [nbr] dimension texture size
!	 
    ! Convert volumetric water content to gravimetric water content
    call vwc2gwc( &
         flg_mbl, & ! I [flg] Mobilization candidate flag
         gwc_sfc, & ! O [kg kg-1] Gravimetric water content
         PWSAT, &   ! I [m3 m-3] Saturated volumetric water content (sand-dependent)
         PWG,   &   ! I [m3 m-3] Volumetric water content
         KSIZE, &   ! I [nbr] loop size
         KSIZE)     ! I [nbr] array dimension
         
  
    ! Factor by which soil moisture increases threshold friction velocity 
frc_thr_ncr_wtr(:) = 1.0d0
    do i=1,KSIZE
       if (flg_mbl(i)) then
          ! Adjust threshold velocity for inhibition by moisture
          if (gwc_sfc(i) > ZWPRM(i)) &
               frc_thr_ncr_wtr(i)=sqrt(1.0d0+1.21d0*(100.0d0*(gwc_sfc(i)-ZWPRM(i)))**0.68d0) ! [frc] FMB99 p. 155 (15)
          frc_thr_ncr_wtr(i)=max(frc_thr_ncr_wtr(i),1.0)
       endif ! endif flg_mbl
    end do   ! end loop over lon

!
!&&&&&&&&& 
ZDP (1)  = 0.1d0                  ! µm 
ZDPLN(1) = log(ZDP(1))           ! µm
DO  IDP = 2 , NDP

  ZDPLN(IDP) = ZDPLN(IDP-1) + DLNDP  
  ZDP(IDP)   = exp(ZDPLN(IDP))
  ZP_DP   = ZDP(IDP)
    !Initialize vertical dust flux
    flx_mss_vrt_dst_ttl_wbn(:)=0.d0
    flx_mss_hrz_slt_ttl_wbn(:)=0.d0
  if (ZDP(IDP) < 1.0d0)  then
           
	   wnd_frc_thr_slt(:)= 0.7d0 &                ! [m s-1] Threshold friction velocity for saltation dry ground
	                       * frc_thr_ncr_wtr(:) & ! [m s-1] Threshold friction velocity for saltation
			       * frc_thr_ncr_drg(:)   ! I [frc] Adjustment for roughness 
	   
  elseif ((ZDP(IDP) >= 1.0d0) .and. (ZDP(IDP) < 10.0d0)) then 

	   wnd_frc_thr_slt(:)=0.5d0 &                 ! [m s-1] Threshold friction velocity for saltation dry ground
	                       * frc_thr_ncr_wtr(:) & ! [m s-1] Threshold friction velocity for saltation
			       * frc_thr_ncr_drg(:)   ! I [frc] Adjustment for roughness 
  
  elseif ((ZDP(IDP) >= 10.0d0) .and. (ZDP(IDP) < 40.0d0)) then 

	   wnd_frc_thr_slt(:)=0.3d0 &                 ! [m s-1] Threshold friction velocity for saltation dry ground
	                       * frc_thr_ncr_wtr(:) & ! [m s-1] Threshold friction velocity for saltation
			       * frc_thr_ncr_drg(:)   ! I [frc] Adjustment for roughness 
       	   
  elseif ((ZDP(IDP) >= 40.0d0) .and. (ZDP(IDP) < 60.0d0)) then 

	   wnd_frc_thr_slt(:)=0.2d0 &                ! [m s-1] Threshold friction velocity for saltation dry ground
	                       * frc_thr_ncr_wtr(:) & ! [m s-1] Threshold friction velocity for saltation
			       * frc_thr_ncr_drg(:)   ! I [frc] Adjustment for roughness 

  else
    
    call wnd_thr_slt_get_mm( &                     !++grini: Removed variables which were not used
         PRHOA,           &     ! I [kg m-3] Midlayer density
         wnd_frc_thr_slt, &     ! O [m s-1] Threshold friction velocity for saltation
	 frc_thr_ncr_wtr, &     ! I [frc] Adjustment for moisture
	 frc_thr_ncr_drg, &     ! I [frc] Adjustment for roughness
	 ZP_DP,           &     ! I [5m]  DP  
         KSIZE,           &     ! I [nbr] loop dimension
         KSIZE)                 ! I [nbr] Array dimension 
         
  endif
!
!  print*,'ap_wnd_thr_get=',MINVAL(wnd_frc_thr_slt(:)),&
!                   MAXVAL(wnd_frc_thr_slt(:))
    do i=1,KSIZE
        if (flg_mbl(i)) then
!    
!           wnd_frc_thr_slt(i)=     &    ! [m s-1] Threshold friction velocity for saltation
!               wnd_frc_thr_slt(i)* &    ! [m s-1] Threshold for dry, flat ground
!               frc_thr_ncr_wtr(i)* &    ! [frc] Adjustment for moisture
!               frc_thr_ncr_drg(i)       ! [frc] Adjustment for roughness            ! proposed
!
    ! Threshold saltation wind speed
    !Needed for the saltation feedback roughening effect
    
          wnd_rfr_thr_slt(i)= & ! [m s-1] Threshold 10 m wind speed for saltation
               wnd_rfr(i)*wnd_frc_thr_slt(i)/PUSTAR(i) !++alfgr 
        endif                  ! endif flg_mbl
    end do                              ! end loop over elements

!

    !CHECK IF THIS CAN BE USED EASILY
    !NEEDS 10M WIND SPEED WHICH IS MAYBE KNOWN MAYBE NOT !
    ! Saltation increases friction speed by roughening surface
    call wnd_frc_slt_get( &
        flg_mbl,          & ! I [flg] Mobilization candidate flag
         PUSTAR,          & ! I [m s-1] Surface friction velocity
         wnd_frc_slt,     & ! O [m s-1] Saltating friction velocity
         wnd_rfr,         & ! I [m s-1] Wind speed at reference height
         wnd_rfr_thr_slt, & ! I [m s-1] Threshold 10 m wind speed for saltation
         KSIZE,           & !I [nbr] loop size
         KSIZE           & !I [nbr] array dimension size
         )              !I [nbr] array dimension size
	 
!
   call flx_mss_hrz_slt_get_mm( &
         PRHOA,                   & ! I [kg m-3] Midlayer density
         flg_mbl,                 & ! I [flg] Mobilization candidate flag
         flx_mss_hrz_slt_ttl_wbn, & ! O [kg m-1 s-1] Vertically integrated streamwise mass flux in wind bin
         wnd_frc_slt,             & ! I [m s-1] Saltating friction velocity
         wnd_frc_thr_slt,         & ! I [m s-1] Threshold friction speed for saltation
         KSIZE,                   & ! I [nbr] loop size
         KSIZE                   & ! I [nbr] array dimension size
         )                   ! I [nbr] array dimension size
!    
    do i=1,KSIZE
       ! Vertically integrated streamwise mass flux in wind bin
       ITEX = ZTEXT(i)
       flx_mss_hrz_slt_ttl_wbn(i)=flx_mss_hrz_slt_ttl_wbn(i) & ! [kg m-2 s-1]
             *ZDSRLV(ITEX,IDP)*mbl_bsn_fct(i)   & ! surface basale relative
             *flx_mss_fdg_fctm  ! [frc] Global mass flux tuning factor (empirical)
    end do ! end loop over lon
	 
!&&&&&&&&
    ! Vertical dust mass flux
    call flx_mss_vrt_dst_Aust_get_mm( &
         dst_slt_flx_rat_ttl, & ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
         flg_mbl, & ! I [flg] Mobilization candidate flag
         flx_mss_hrz_slt_ttl_wbn, &  ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
         flx_mss_vrt_dst_ttl_wbn, &  ! O [kg m-2 s-1] Total vertical mass flux of dust
         PRHOA,               &  ! I                      ! proposed
         wnd_frc_thr_slt,     &  ! I [m s-1] Threshold friction speed for saltation
         KSIZE,               &  ! I [nbr] loop size
         KSIZE               &  ! I [nbr] array dimension size
         )                   ! I [nbr] array dimension size

    do i=1,KSIZE
       ! Vertically integrated streamwise mass flux in wind bin
       flx_mss_vrt_dst_tot(i)=flx_mss_vrt_dst_tot(i)+flx_mss_vrt_dst_ttl_wbn(i)
    end do ! end loop over lon
! 
ENDDO      ! end loop over NDP

!    
    PSFDST(:) = flx_mss_vrt_dst_tot(:)
! 	 
  end subroutine dustflux_get_MB

END MODULE MODE_DSTMBL_MB
