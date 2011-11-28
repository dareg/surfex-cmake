!     ######spl
MODULE MODE_DSTMBLUTL ! [mdl] Mobilization utilities
!
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
  USE PARKIND1  ,ONLY : JPRB
!
  implicit none
  public ! [stt] Symbols are public unless individually qualified as private

contains
  
  subroutine wnd_frc_thr_slt_get( &
       dns_mdp, & ! I [kg m-3] Midlayer density
       wnd_frc_thr_slt,& ! O [m s-1] Threshold friction velocity for saltation
       PLON,           & ! I [nbr] number of points to calculate 
       PLOND )           ! I [nbr] dimension of arrays

    ! Purpose: Compute dry threshold friction velocity for saltation
    !++grini use dstgrd ! [mdl] Dust grid sizes
    !++grini use pmgrid ! [mdl] Spatial resolution parameters
    !++grini use dstcst ! [mdl] Physical constants for dust routines
    USE MODD_CSTS, only: XG
    implicit none
    ! Parameters
    real,parameter::dmt_slt_opt=75.0e-6 ! [m] Optimal diameter for saltation, IvW82 p. 117 Fgr. 8, Pye87 p. 31, MBA97 p. 4388, SRL96 (2)
    real,parameter::dns_slt=2650.0 ! [kg m-3] Density of optimal saltation particles, MBA97 p. 4388 
    ! Input
    !++grini real,intent(in)::dns_aer(dst_nbr) ! [kg m-3] Particle density
    !++grini real,intent(in)::dmt_aer(dst_nbr) ! [m] Particle diameter
    integer, intent(in) :: PLON         ! [nbr] loop size
    integer, intent(in) :: PLOND        ! [nbr] dimension size
    real,intent(in)::dns_mdp(plond) ! [kg m-3] Midlayer density
    ! Output
    real,intent(out)::wnd_frc_thr_slt(plond) ! [m s-1] Threshold friction velocity for saltation
    ! Local
    integer lon_idx ! [idx] Counting index
    real ryn_nbr_frc_thr_prx_opt ! [frc] Threshold friction Reynolds number approximation for optimal size
    real ryn_nbr_frc_thr_opt_fnc ! [frc] Threshold friction Reynolds factor for saltation calculation
    real dns_fct ! Density ratio factor for saltation calculation
    real icf_fct ! Interparticle cohesive forces factor for saltation calculation
    real tmp1 ! Factor in saltation computation
    real grv_sfc   !Gravitation (m s-2) 
    REAL(KIND=JPRB) :: ZHOOK_HANDLE   

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:WND_FRC_THR_SLT_GET',0,ZHOOK_HANDLE)

    ! Main Code
    !++alfgr set gravitation to a real number which is close to 9.8 on earth in SI units
    grv_sfc=XG

    ! Initialize some variables
    ! MaB95 pzn. for Re*t(D_opt) circumvents iterative solution
    ryn_nbr_frc_thr_prx_opt=0.38d0+1331.0d0*(100.0d0*dmt_slt_opt)**1.56d0 ! [frc] "B" MaB95 p. 16417 (5)
    ! Given Re*t(D_opt), compute time independent factors contributing to u*t
    icf_fct=1.0d0+6.0d-07/(dns_slt*grv_sfc*(dmt_slt_opt**2.5d0)) ! [frc] IvW82 p. 115 (6) MaB95 p. 16417 (4) Interparticle cohesive forces
    dns_fct=dns_slt*grv_sfc*dmt_slt_opt ! IvW82 p. 115 (6) MaB95 p. 16417 (4)
    if (ryn_nbr_frc_thr_prx_opt < 0.03d0) then
       stop 'dst: wnd_frc_thr_slt_get() reports ryn_nbr_frc_thr_prx_opt < 0.03'
    else if (ryn_nbr_frc_thr_prx_opt < 10.0d0) then
       ryn_nbr_frc_thr_opt_fnc=-1.0d0+1.928d0*(ryn_nbr_frc_thr_prx_opt**0.0922d0) ! [frc] IvW82 p. 114 (3), MaB95 p. 16417 (6)
       ryn_nbr_frc_thr_opt_fnc=0.1291d0*0.1291d0/ryn_nbr_frc_thr_opt_fnc ! [frc] 
    else 
       ryn_nbr_frc_thr_opt_fnc=1.0d0-0.0858d0*exp(-0.0617d0*(ryn_nbr_frc_thr_prx_opt-10.0d0)) ! [frc] IvW82 p. 114 (3), MaB95 p. 16417 (7)
       ryn_nbr_frc_thr_opt_fnc=0.120d0*0.120d0*ryn_nbr_frc_thr_opt_fnc*ryn_nbr_frc_thr_opt_fnc! [frc]
    endif ! endif
    ! This method minimizes the number of square root computations performed
    tmp1=sqrt(icf_fct*dns_fct*ryn_nbr_frc_thr_opt_fnc)
    do lon_idx=1,plon
       wnd_frc_thr_slt(lon_idx)=tmp1/sqrt(dns_mdp(lon_idx)) ! [m s-1] Threshold friction velocity for saltation dry ground 
    end do ! end loop over lon
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:WND_FRC_THR_SLT_GET',1,ZHOOK_HANDLE)
  end subroutine wnd_frc_thr_slt_get ! end wnd_frc_thr_slt_get()
!&&&&&
  subroutine wnd_thr_slt_get_mm( &
       dns_mdp, & ! I [kg m-3] Midlayer density
       wnd_frc_thr_slt,& ! O [m s-1] Threshold friction velocity for saltation
       frc_thr_ncr_wtr,& ! I [frc] Adjustment for moisture
       frc_thr_ncr_drg,& ! I[frc] Adjustment for roughness
       PDP,            & ! I [5m] 
       PLON,           & ! I [nbr] number of points to calculate 
       PLOND)            ! I [nbr] dimension of arrays
    ! Purpose: Compute dry threshold friction velocity for saltation
    !++grini use dstgrd ! [mdl] Dust grid sizes
    !++grini use pmgrid ! [mdl] Spatial resolution parameters
    !++grini use dstcst ! [mdl] Physical constants for dust routines
    USE MODD_CSTS, only: XG
    implicit none
    ! Parameters
    real,parameter:: dns_slt=2650.0 ! [kg m-3] Density of optimal saltation particles, MBA97 p. 4388 
    real,parameter:: x=1.56d0  
    real,parameter:: a=1331.d0  
    real,parameter:: b=0.38d0  
    
    ! Input
    !++grini real,intent(in)::dns_aer(dst_nbr) ! [kg m-3] Particle density
    !++grini real,intent(in)::dmt_aer(dst_nbr) ! [m] Particle diameter
    integer, intent(in) :: PLON         ! [nbr] loop size
    integer, intent(in) :: PLOND        ! [nbr] dimension size
    
    real,intent(in)::dns_mdp(plond) ! [kg m-3] Midlayer density
    real,intent(in)::PDP            ![5m] 
    real,intent(in)::frc_thr_ncr_wtr(plond) ! [frc] Factor by which moisture increases threshold friction velocity
    real,intent(in)::frc_thr_ncr_drg(plond) ! [frc] Factor by which roughness increases threshold friction velocity
    
    ! Output
    real,intent(out)::wnd_frc_thr_slt(plond) ! [m s-1] Threshold friction velocity for saltation
    ! Local
    integer lon_idx ! [idx] Counting index
    real ryn, ryn1     ! [frc] Threshold friction Reynolds number approximation for optimal size
    real dp, k, k1, k2 ! [frc] Threshold friction Reynolds factor for saltation calculation
    real dns_fct ! Density ratio factor for saltation calculation
    real icf_fct ! Interparticle cohesive forces factor for saltation calculation
    real tmp     ! Factor in saltation computation
    real grv_sfc   !Gravitation (m s-2) 
    real wnd_thr   ! [m s-1] Threshold friction velocity for saltation dry ground 

    REAL(KIND=JPRB) :: ZHOOK_HANDLE     
    
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:WND_THR_SLT_GET_MM',0,ZHOOK_HANDLE)

    ! Main Code
    !++alfgr set gravitation to a real number which is close to 9.8 on earth in SI units
    grv_sfc=XG
       
    do lon_idx=1,plon
      
        dp=1E-6*PDP
        ryn=b + a*(100.0*dp)**x
        k1=dns_slt*grv_sfc*dp
        k2 =1.0d0 + 6.0E-7/(dns_slt*grv_sfc*(dp**2.5))

        IF (ryn.lt.0.03) THEN 
        STOP
        ELSEIF (ryn.lt.10.0) THEN 
        ryn1=-1.0+1.928*(ryn**0.0922)
        ryn1 = 0.1291*0.1291/ryn1
        ELSE
        ryn1=1.0-0.0858*exp(-0.0617*(ryn-10.0))
        ryn1=0.129*0.129*ryn1*ryn1
        END IF
         tmp=sqrt(k2*k1*ryn1)
         wnd_thr=tmp/sqrt(dns_mdp(lon_idx)) ! [m s-1] Threshold friction velocity for saltation dry ground
	 wnd_frc_thr_slt(lon_idx)=wnd_thr*frc_thr_ncr_wtr(lon_idx)* &    ! [frc] Adjustment for moisture
	                                  frc_thr_ncr_drg(lon_idx)       ! [frc] Adjustment for roughness 
    end do
!
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:WND_THR_SLT_GET_MM',1,ZHOOK_HANDLE)
  end subroutine wnd_thr_slt_get_mm 
!
!&&&&&  
  subroutine wnd_rfr_thr_slt_get( &
       wnd_frc, & ! I [m s-1] Surface friction velocity
       wnd_frc_thr_slt, & ! I [m s-1] Threshold friction velocity for saltation
       wnd_mdp, & ! I [m s-1] Surface layer mean wind speed
       wnd_rfr, & ! I [m s-1] Wind speed at reference height
       wnd_rfr_thr_slt, & ! O [m s-1] Threshold 10 m wind speed for saltation
       PLON,            & ! I [nbr] loop size
       PLOND            & ! I [nbr] dimension of arrays
       )
    ! Purpose: Compute threshold 10 m wind speed
    !++ alfgr use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    ! Input
    integer, intent(in) :: PLON     ![nbr] loop size
    integer, intent(in) :: PLOND    ![nbr] dimension size
    real,intent(in)::wnd_frc(plond) ! [m s-1] Friction velocity
    real,intent(in)::wnd_frc_thr_slt(plond) ! [m s-1] Threshold friction velocity for saltation
    real,intent(in)::wnd_mdp(plond) ! [m s-1] Surface layer mean wind speed
    real,intent(in)::wnd_rfr(plond) ! [m s-1] Wind speed at reference height
    ! Output
    real,intent(out)::wnd_rfr_thr_slt(plond) ! [m s-1] Threshold 10 m wind speed for saltation
    ! Local
    integer lon_idx ! [idx] Counting index
    REAL(KIND=JPRB) :: ZHOOK_HANDLE       

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:WND_RFR_THR_SLT_GET',0,ZHOOK_HANDLE)
    
    ! Main Code
    ! Compute threshold horizontal wind speed at reference height
    do lon_idx=1,plon
       ! A more complicated procedure would recompute mno_lng for wnd_frc_thr,
       ! and then integrate vertically from rgh_mmn+hgt_zpd to hgt_rfr
       ! wnd_crc_fct is (1/k)*[ln(z-D)/z0 - psi(zeta2) + psi(zeta1)]
       wnd_rfr_thr_slt(lon_idx)=wnd_frc_thr_slt(lon_idx)*wnd_rfr(lon_idx)/wnd_frc(lon_idx) ! [m s-1]
    end do ! end loop over lon
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:WND_RFR_THR_SLT_GET',1,ZHOOK_HANDLE)    
  end subroutine wnd_rfr_thr_slt_get ! end wnd_rfr_thr_slt_get()
  
  subroutine vwc2gwc( &
       flg_mbl, & ! I [flg] Mobilization candidate flag
       gwc_sfc, & ! O [kg kg-1] Gravimetric water content
       vwc_sat, & ! I [m3 m-3] Saturated volumetric water content (sand-dependent)
       vwc_sfc, & ! I [m3 m-3] Volumetric water content
       PLON,    &  !I [nbr] loop size
       PLOND    )  !I [nbr] dimension size
    ! Purpose: Convert volumetric water content to gravimetric water content
    !++alfgr use pmgrid ! [mdl] Spatial resolution parameters
    use MODD_CSTS, only : XRHOLI    ! [kg/m3] density of liquid water
    !++alfgr use dstblm,only:dns_H2O_lqd_std ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    implicit none
    ! Parameters
    real,parameter::dns_prt_sfc=2650.0 ! [kg m-3] Dry density of soil particles (excluding pores)
    ! Input
    integer, intent(in) :: PLON        ! [nbr] loop size
    integer, intent(in) :: PLOND       ! [nbr] dimension size
    logical,intent(in)::flg_mbl(plond) ! [flg] Mobilization candidate flag
    real,intent(in)::vwc_sat(plond) ! [m3 m-3] Saturated volumetric water content (sand-dependent)
    real,intent(in)::vwc_sfc(plond) ! [m3 m-3] Volumetric water content
    ! Output
    real,intent(out)::gwc_sfc(plond) ! [kg kg-1] Gravimetric water content
    ! Local
    integer lon_idx             ! [idx] Counting index
    real dns_blk_dry(plond) ! [kg m-3] Bulk density of dry surface soil
    real dns_H2O_lqd_std    ! [kg m-3] Density of water
    REAL(KIND=JPRB) :: ZHOOK_HANDLE     
    ! Main Code

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:VWC2GWC',0,ZHOOK_HANDLE)
    !Set liquid water density
    dns_H2O_lqd_std=XRHOLI

    ! Initialize output
    do lon_idx=1,plon
       if (flg_mbl(lon_idx)) then
          ! Assume volume of air pores when dry equals saturated VWC
          ! This implies air pores are completely filled by water in saturated soil
          dns_blk_dry(lon_idx)=dns_prt_sfc*(1.0-vwc_sat(lon_idx)) ! [kg m-3] Bulk density of dry surface soil
          gwc_sfc(lon_idx)=vwc_sfc(lon_idx)*dns_H2O_lqd_std/dns_blk_dry(lon_idx) ! O [kg kg-1] Gravimetric water content
       endif ! endif flg_mbl
    enddo ! end loop over lon
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:VWC2GWC',1,ZHOOK_HANDLE)
  end subroutine vwc2gwc
  
  subroutine frc_thr_ncr_wtr_get( &
       flg_mbl, & ! I [flg] Mobilization candidate flag
       frc_thr_ncr_wtr, & ! O [frc] Factor by which moisture increases threshold friction velocity
       mss_frc_cly, & ! I [frc] Mass fraction of clay
       gwc_sfc, & ! I [kg kg-1] Gravimetric water content
       PLON,    & ! I [nbr] loop size 
       PLOND)     ! I [nbr] dimension size
   
    ! Purpose: Compute factor by which soil moisture increases threshold friction velocity 
    ! This parameterization is based on FMB99
    !++alfgr use pmgrid ! [mdl] Spatial resolution parameters


    implicit none
    ! Parameters
    real,parameter::awet=3. ! facteur qui tient compte du comportement de mod?le
                            ! dans le calcul de l'humidit? du sol(a posteriori)
    ! Input
    integer, intent(in) :: PLON        !I [nbr] loop size
    integer, intent(in) :: PLOND       !I [nbr] dimension size
    logical,intent(in)::flg_mbl(plond) ! I [flg] Mobilization candidate flag
    real,intent(in)::mss_frc_cly(plond) ! [frc] Mass fraction of clay
    real,intent(in)::gwc_sfc(plond) ! [kg kg-1] Gravimetric water content
    ! Output
    real,intent(out)::frc_thr_ncr_wtr(plond) ! [frc] Factor by which moisture increases threshold friction velocity
    ! Local
    integer lon_idx ! [idx] Counting index
    real gwc_thr(plond) ! [kg kg-1] Threshold gravimetric water content
    real wet_thr1(plond) ! [kg kg-1] Threshold gravimetric water content 
    real wet_thr2(plond) ! [kg kg-1] Threshold gravimetric water content 
    REAL(KIND=JPRB) :: ZHOOK_HANDLE 

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FRC_THR_NCR_WTR_GET',0,ZHOOK_HANDLE)
    ! Main Code
    ! Initialize output
    frc_thr_ncr_wtr(:)=1.0d0 ! [frc] Factor by which moisture increases threshold friction velocity
    do lon_idx=1,plon
      if (flg_mbl(lon_idx)) then
          ! Adjust threshold velocity for inhibition by moisture
          ! frc_thr_ncr_wtr(lon_idx)=exp(22.7d0*vwc_sfc(lon_idx)) ! [frc] SRL96
          ! Compute threshold soil moisture based on clay content

          ! Old Alf Grini code: bug on DEAD detected ? 
          ! Modification proposed by M. Mokhtari.. Accepted for all cases
          !  if (CVERMOD=='CMDVER') then
              wet_thr1(lon_idx)=awet*mss_frc_cly(lon_idx)*(0.17d0+0.14d0*mss_frc_cly(lon_idx)) ! proposed.
             wet_thr2(lon_idx)=max(0.02,wet_thr1(lon_idx))                  ! proposed
             gwc_thr(lon_idx)=min(0.14,wet_thr2(lon_idx))                   ! proposed  
          !  else 
             ! gwc_thr=mss_frc_cly*(0.17d0+0.14d0*mss_frc_cly) ! [m3 m-3] FMB99 p. 155 (14)
             ! fxm: 19991105 remove factor of mss_frc_cly from gwc_thr to improve large scale behavior
          ! Begin Old Alf code 
          !   gwc_thr(lon_idx)=0.17d0+0.14d0*mss_frc_cly(lon_idx) ! [m3 m-3] 
          !  endif
          ! End Old Alf code 
          if (gwc_sfc(lon_idx) > gwc_thr(lon_idx)) &
               frc_thr_ncr_wtr(lon_idx)=sqrt(1.0d0+1.21d0*(100.0d0*(gwc_sfc(lon_idx)-gwc_thr(lon_idx)))**0.68d0) ! [frc] FMB99 p. 155 (15)
       endif ! endif flg_mbl
    enddo ! end loop over lon
    ! Uncomment following line to remove all dependence on gwc_sfc
    ! frc_thr_ncr_wtr(lon_idx)=1.0d0 ! [frc]
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FRC_THR_NCR_WTR_GET',1,ZHOOK_HANDLE)
  end subroutine frc_thr_ncr_wtr_get
  
  !++alfgr fxm: Fix this so that we can actually use rgh_mmn_mbl different for grid cells
  subroutine frc_thr_ncr_drg_get( &
       frc_thr_ncr_drg, & ! O [frc] Factor by which surface roughness increases threshold friction velocity
       rgh_mmn_mbl, & ! I [m] Roughness length momentum for erodible surfaces
       rgh_mmn_smt, & ! I [m] Smooth roughness length
       PLON,        & ! I [nbr] loop size
       PLOND )        ! I [nbr] dimension size 
    ! Purpose: Compute factor by which surface roughness increases threshold friction velocity
    ! This parameterization is based on MaB95 and GMB98
    !++alfgr use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    ! Input
    integer, intent(in) :: PLON  ! [nbr] loop size
    integer, intent(in) :: PLOND ! [nbr] dimension size
    real,intent(in)::rgh_mmn_mbl(PLOND) ! [m] Roughness length momentum for erodible surfaces
    real,intent(in)::rgh_mmn_smt ! [m] Smooth roughness length
    ! Output
    real,intent(out)::frc_thr_ncr_drg(plond) ! [frc] Factor by which roughness increases threshold friction velocity
    ! Local
    integer lon_idx ! [idx] Counting index
    real :: wnd_frc_fsh_frc(PLOND) ! [frc] Efficient fraction of wind friction
    real :: wnd_frc_fsh_frc_rcp(PLOND) ! [frc] Reciprocal of wnd_frc_fsh_frc
    REAL(KIND=JPRB) :: ZHOOK_HANDLE 

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FRC_THR_NCR_DRG_GET',0,ZHOOK_HANDLE)    
    ! Main Code
    
    ! Adjust threshold velocity for inhibition by roughness elements
    do lon_idx=1,plon
       wnd_frc_fsh_frc(lon_idx)= & ! [frc] MaB95 p. 16420, GMB98 p. 6207
            +1.0d0-log(rgh_mmn_mbl(lon_idx)/rgh_mmn_smt)/log(0.35d0*((0.1d0/rgh_mmn_smt)**0.8d0))
       wnd_frc_fsh_frc(lon_idx)=min(1.,wnd_frc_fsh_frc(lon_idx)) !exclude larger values than 1
       wnd_frc_fsh_frc(lon_idx)=max(1.e-6,wnd_frc_fsh_frc(lon_idx)) !exclude smaller values than 1.d-6
       if (wnd_frc_fsh_frc(lon_idx) <= 0.0d0.or.wnd_frc_fsh_frc(lon_idx) > 1.0d0) stop &
            'dst: ERROR frc_thr_ncr_drg_get() reports wnd_frc_fsh_frc out of range'
       wnd_frc_fsh_frc_rcp(lon_idx)=1.0d0/wnd_frc_fsh_frc(lon_idx) ! [frc] 
       frc_thr_ncr_drg(lon_idx)=wnd_frc_fsh_frc_rcp(lon_idx) ! [frc]
       ! fxm: 19991012 Set frc_thr_ncr_drg=1.0, equivalent to assuming mobilization takes place at smooth roughness length
    enddo
    !++alfgr removed this frc_thr_ncr_drg=1.0d0
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FRC_THR_NCR_DRG_GET',1,ZHOOK_HANDLE) 
  end subroutine frc_thr_ncr_drg_get ! end frc_thr_ncr_drg_get()
!
  subroutine frc_thr_ncr_drg_get_mm( &
       frc_thr_ncr_drg, & ! O [frc] Factor by which surface roughness increases threshold friction velocity
       rgh_mmn_mbl, & ! I [m] Roughness length momentum for erodible surfaces
       rgh_mmn_smt, & ! I [m] Smooth roughness length
       PTEXT      , & ! I Soil texture
       PLON,        & ! I [nbr] loop size
       PLOND,       & ! I [nbr] dimension size 
       NTEX)          ! I [nbr] Texture size
    ! Purpose: Compute factor by which surface roughness increases threshold friction velocity
    ! This parameterization is based on MaB95 and GMB98
    !++alfgr use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    ! Input
    integer, intent(in) :: PLON  ! [nbr] loop size
    integer, intent(in) :: PLOND ! [nbr] dimension size
    integer, intent(in) :: NTEX  ! [nbr] dimension texture size
    real,intent(in)::rgh_mmn_mbl(PLOND) ! [m] Roughness length momentum for erodible surfaces
    real,intent(in)::rgh_mmn_smt(NTEX)  ! [m] Smooth roughness length
    integer,intent(in)::PTEXT(PLOND)       ! Soil Texture
    
    ! Output
    real,intent(out)::frc_thr_ncr_drg(plond) ! [frc] Factor by which roughness increases threshold friction velocity
    ! Local
    integer lon_idx ! [idx] Counting index
    integer ITEX ! [idx] Counting index texture

    real rgh_z0  !  [m] Roughness length momentum for erodible surfaces
    
    real :: wnd_frc_fsh_frc(PLOND) ! [frc] Efficient fraction of wind friction
    real :: wnd_frc_fsh_frc_rcp(PLOND) ! [frc] Reciprocal of wnd_frc_fsh_frc
    REAL(KIND=JPRB) :: ZHOOK_HANDLE 

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FRC_THR_NCR_DRG_GET_MM',0,ZHOOK_HANDLE)     
    ! Main Code
    
    ! Adjust threshold velocity for inhibition by roughness elements
    do lon_idx=1,plon
       ITEX = PTEXT(lon_idx)
       rgh_z0=rgh_mmn_mbl(lon_idx)
       if (rgh_z0<=rgh_mmn_smt(ITEX)) rgh_z0=rgh_mmn_smt(ITEX)
       wnd_frc_fsh_frc(lon_idx)= & ! [frc] MaB95 p. 16420, GMB98 p. 6207
            +1.0d0-log(rgh_z0/rgh_mmn_smt(ITEX))/log(0.35d0*((0.1d0/rgh_mmn_smt(ITEX))**0.8d0))
       wnd_frc_fsh_frc(lon_idx)=min(1.,wnd_frc_fsh_frc(lon_idx)) !exclude larger values than 1
       wnd_frc_fsh_frc(lon_idx)=max(1.e-6,wnd_frc_fsh_frc(lon_idx)) !exclude smaller values than 1.d-6
       if (wnd_frc_fsh_frc(lon_idx) <= 0.0d0.or.wnd_frc_fsh_frc(lon_idx) > 1.0d0) stop &
            'dst: ERROR frc_thr_ncr_drg_get() reports wnd_frc_fsh_frc out of range'
       wnd_frc_fsh_frc_rcp(lon_idx)=1.0d0/wnd_frc_fsh_frc(lon_idx) ! [frc] 
       frc_thr_ncr_drg(lon_idx)=wnd_frc_fsh_frc_rcp(lon_idx) ! [frc]
       ! fxm: 19991012 Set frc_thr_ncr_drg=1.0, equivalent to assuming mobilization takes place at smooth roughness length
    enddo
    !++alfgr removed this frc_thr_ncr_drg=1.0d0
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FRC_THR_NCR_DRG_GET_MM',1,ZHOOK_HANDLE)     
  end subroutine frc_thr_ncr_drg_get_mm ! end frc_thr_ncr_drg_get()
!  
  subroutine wnd_frc_slt_get( &
       flg_mbl, & ! I [flg] Mobilization candidate flag
       wnd_frc, & ! I [m s-1] Surface friction velocity
       wnd_frc_slt, & ! O [m s-1] Saltating friction velocity
       wnd_rfr, & ! I [m s-1] Wind speed at reference height
       wnd_rfr_thr_slt, & ! I [m s-1] Threshold 10 m wind speed for saltation
       PLON,           &  ! I [nbr] loop size
       PLOND)             ! I [nbr] dimension size
    ! Purpose: Compute the saltating friction velocity
    ! Saltation increases friction speed by roughening surface, AKA "Owen's effect"
    ! This acts as a positive feedback to the friction speed
    ! GMB98 parameterized this feedback in terms of 10 m windspeeds
    !++alfgr use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    ! Input
    integer, intent(in) :: PLON        ! I [nbr] loop size
    integer, intent(in) :: PLOND       ! I [nbr] dimension size
    logical,intent(in)::flg_mbl(plond) ! I [flg] Mobilization candidate flag
    real,intent(in)::wnd_frc(plond) ! I [m s-1] Surface friction velocity
    real,intent(in)::wnd_rfr(plond) ! I [m s-1] Wind speed at reference height
    real,intent(in)::wnd_rfr_thr_slt(plond) ! I [m s-1] Threshold 10 m wind speed for saltation
    ! Output
    real,intent(out)::wnd_frc_slt(plond) ! O [m s-1] Saltating friction velocity
    ! Local
    integer lon_idx ! [idx] Counting index
    real wnd_rfr_dlt ! [m s-1] Reference windspeed excess over threshold
    real wnd_frc_slt_dlt ! [m s-1] Friction velocity increase from saltation
    REAL(KIND=JPRB) :: ZHOOK_HANDLE 

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:WND_FRC_SLT_GET',0,ZHOOK_HANDLE)   
    ! Main Code
    ! Compute saltating friction velocity, accounting for "Owen's effect"
    wnd_frc_slt=wnd_frc ! [m s-1] Saltating friction velocity
    do lon_idx=1,plon
       if (flg_mbl(lon_idx).and.wnd_rfr(lon_idx) >= wnd_rfr_thr_slt(lon_idx)) then
          ! Saltation roughens the boundary layer, AKA "Owen's effect"
          ! GMB98 p. 6206 Fig. 1 shows observed/computed u* dependence on observed U(1 m)
          ! GMB98 p. 6209 (12) has u* in cm s-1 and U, Ut in m s-1, personal communication, D. Gillette, 19990529
          ! With everything in MKS, the 0.3 coefficient in GMB98 (12) becomes 0.003 
          ! Increase in friction velocity due to saltation varies as square of 
          ! difference between reference wind speed and reference threshold speed 
          wnd_rfr_dlt=wnd_rfr(lon_idx)-wnd_rfr_thr_slt(lon_idx)
          wnd_frc_slt_dlt=0.003d0*wnd_rfr_dlt*wnd_rfr_dlt ! [m s-1] Friction velocity increase from saltation GMB98 p. 6209
          wnd_frc_slt(lon_idx)=wnd_frc(lon_idx)+wnd_frc_slt_dlt ! [m s-1] Saltating friction velocity
       endif ! endif wnd_frc_mbl > wnd_frc_thr_slt
    end do ! end loop over lon
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:WND_FRC_SLT_GET',1,ZHOOK_HANDLE)
  end subroutine wnd_frc_slt_get ! end wnd_frc_slt_get()
!&&&&&&&&&&&&&&&&&&&

  subroutine wnd_frc_slt_get_mm( &
       flg_mbl,         & ! I [flg] Mobilization candidate flag
       wnd_frc,         & ! I [m s-1] Surface friction velocity
       wnd_frc_thr_slt, & ! I [m s-1] Threshold friction velocity for saltation
       wnd_frc_slt,     & ! O [m s-1] Saltating friction velocity
       wnd_rfr,         & ! I [m s-1] Wind speed at reference height
       PLON,            & ! I [nbr] loop size
       PLOND)             ! I [nbr] dimension size
       
    ! Purpose: Compute the saltating friction velocity
    ! Saltation increases friction speed by roughening surface, AKA "Owen's effect"
    ! This acts as a positive feedback to the friction speed
    ! GMB98 parameterized this feedback in terms of 10 m windspeeds
    !++alfgr use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    ! Input
    integer, intent(in) :: PLON        ! I [nbr] loop size
    integer, intent(in) :: PLOND       ! I [nbr] dimension size
    logical,intent(in)::flg_mbl(plond) ! I [flg] Mobilization candidate flag
    real,intent(in)::wnd_frc(plond)    ! I [m s-1] Surface friction velocity
    real,intent(in)::wnd_frc_thr_slt(plond)    ! I [m s-1] Threshold surface friction velocity
    real,intent(in)::wnd_rfr(plond)            ! I [m s-1] Wind speed at reference height
    ! Output
    real,intent(out)::wnd_frc_slt(plond)    ! O [m s-1] Saltating friction velocity
    ! Local
    integer lon_idx ! [idx] Counting index
    real wnd_rfr_dlt ! [m s-1] Reference windspeed excess over threshold
    real wnd_frc_slt_dlt ! [m s-1] Friction velocity increase from saltation
    real :: wnd_rfr_thr_slt(PLOND) 
    REAL(KIND=JPRB) :: ZHOOK_HANDLE 

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:WND_FRC_SLT_GET_MM',0,ZHOOK_HANDLE)       
    ! Main Code
    ! Compute saltating friction velocity, accounting for "Owen's effect"
    
    do lon_idx=1,plon
    
          wnd_rfr_thr_slt(lon_idx)=     &                 ! [m s-1] Threshold 10 m wind speed for saltation
               wnd_rfr(lon_idx)         &                 ! [m s-1] 10 m wind speed 
	      *wnd_frc_thr_slt(lon_idx) &                 ! [m s-1] Threshold wind speed for saltation
	      /wnd_frc(lon_idx)                           ! [m s-1] wind speed for saltation++alfgr 
	       
          wnd_frc_slt(lon_idx)= wnd_frc(lon_idx) ! [m s-1] Saltating friction velocity
	          
       if (flg_mbl(lon_idx).and.(wnd_rfr(lon_idx) >= wnd_rfr_thr_slt(lon_idx))) then
          ! Saltation roughens the boundary layer, AKA "Owen's effect"
          ! GMB98 p. 6206 Fig. 1 shows observed/computed u* dependence on observed U(1 m)
          ! GMB98 p. 6209 (12) has u* in cm s-1 and U, Ut in m s-1, personal communication, D. Gillette, 19990529
          ! With everything in MKS, the 0.3 coefficient in GMB98 (12) becomes 0.003 
          ! Increase in friction velocity due to saltation varies as square of 
          ! difference between reference wind speed and reference threshold speed 
          wnd_rfr_dlt=wnd_rfr(lon_idx)-wnd_rfr_thr_slt(lon_idx)
          wnd_frc_slt_dlt=0.003d0*wnd_rfr_dlt*wnd_rfr_dlt ! [m s-1] Friction velocity increase from saltation GMB98 p. 6209
          wnd_frc_slt(lon_idx)=wnd_frc(lon_idx)+wnd_frc_slt_dlt ! [m s-1] Saltating friction velocity
       endif ! endif wnd_frc_mbl > wnd_frc_thr_slt
    end do   ! end loop over lon
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:WND_FRC_SLT_GET_MM',1,ZHOOK_HANDLE) 
  end subroutine wnd_frc_slt_get_mm 
!
!&&&&&&&&&&&&&&&&&&&  
  subroutine flx_mss_hrz_slt_ttl_Whi79_get( &
       dns_mdp, & ! I [kg m-3] Midlayer density
       flg_mbl, & ! I [flg] Mobilization candidate flag
       flx_mss_hrz_slt_ttl, & ! O [kg m-1 s-1] Vertically integrated streamwise mass flux
       wnd_frc, & ! I [m s-1] Surface friction velocity
       wnd_frc_thr_slt, & ! I [m s-1] Threshold friction speed for saltation
       PLON,            & ! I [nbr] loop size
       PLOND)             ! I [nbr] dimension size 
    ! Purpose: Compute vertically integrated streamwise mass flux of particles
    ! Theory: Uses method proposed by White (1979)
    ! fxm: use surface air density not midlayer density
    !++alfgr use pmgrid ! [mdl] Spatial resolution parameters
    !++alfgr use dstcst ! [mdl] Physical constants for dust routines
    use MODD_CSTS, only : XG  ! Gravitation constant
    implicit none
    ! Parameters
    real,parameter::cst_slt=2.61 ! [frc] Saltation constant Whi79 p. 4648, MaB97 p. 16422 
    ! Input
    integer, intent(in) :: PLON        ! I [nbr] loop size
    integer, intent(in) :: PLOND       ! I [nbr] dimension size
    logical,intent(in)::flg_mbl(plond) ! I [flg] Mobilization candidate flag
    real,intent(in)::dns_mdp(plond) ! I [kg m-3] Midlayer density
    real,intent(in)::wnd_frc(plond) ! I [m s-1] Surface friction velocity
    real,intent(in)::wnd_frc_thr_slt(plond) ! I [m s-1] Threshold friction speed for saltation
    ! Output
    real,intent(out)::flx_mss_hrz_slt_ttl(plond) ! O [kg m-1 s-1] Vertically integrated streamwise mass flux
    integer rcd ! [rcd] Return success code
    ! Local
    real wnd_frc_rat ! [frc] Ratio of wind friction threshold to wind friction
    real grv_sfc     ! [m s-2] Gravitation constant
    integer lon_idx ! [idx] Counting index for lon
    REAL(KIND=JPRB) :: ZHOOK_HANDLE 

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FLX_MSS_HRZ_SLT_TTL_WHI79_GET',0,ZHOOK_HANDLE)       
    ! Main Code
    
    !++alfgr Set gravitation.
    grv_sfc = XG

    ! Initialize output
    flx_mss_hrz_slt_ttl(:)=0.0d0 ! [kg m-1 s-1]
    
    do lon_idx=1,plon
       if (flg_mbl(lon_idx).and.wnd_frc(lon_idx) > wnd_frc_thr_slt(lon_idx)) then
          wnd_frc_rat=wnd_frc_thr_slt(lon_idx)/wnd_frc(lon_idx) ! [frc]
          flx_mss_hrz_slt_ttl(lon_idx)= & ! [kg m-1 s-1] 
               cst_slt*dns_mdp(lon_idx)*(wnd_frc(lon_idx)**3.0d0)* &
               (1.0d0-wnd_frc_rat)*(1.0d0+wnd_frc_rat)*(1.0d0+wnd_frc_rat)/grv_sfc ! Whi79 p. 4648 (19), MaB97 p. 16422 (28)
       endif ! endif 
    end do ! end loop over lon
    
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FLX_MSS_HRZ_SLT_TTL_WHI79_GET',1,ZHOOK_HANDLE) 
  end subroutine flx_mss_hrz_slt_ttl_Whi79_get ! end flx_mss_hrz_slt_ttl_Whi79_get()
!
  subroutine flx_mss_hrz_slt_get_mm( &
       dns_mdp,             & ! I [kg m-3] Midlayer density
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       flx_mss_hrz_slt_ttl, & ! O [kg m-1 s-1] Vertically integrated streamwise mass flux
       wnd_frc,             & ! I [m s-1] Surface friction velocity
       wnd_frc_thr_slt,     & ! I [m s-1] Threshold friction speed for saltation
       PLON,                & ! I [nbr] loop size
       PLOND                & ! I [nbr] dimension size 
       )               ! I [nbr] dimension size 
       
    ! Purpose: Compute vertically integrated streamwise mass flux of particles
    ! Theory: Uses method proposed by White (1979)
    ! fxm: use surface air density not midlayer density
    !++alfgr use pmgrid ! [mdl] Spatial resolution parameters
    !++alfgr use dstcst ! [mdl] Physical constants for dust routines
    use MODD_CSTS, only : XG  ! Gravitation constant
    implicit none
    ! Parameters
!    real,parameter::cst_slt=2.61 ! [frc] Saltation constant Whi79 p. 4648, MaB97 p. 16422 
    ! Input
    integer, intent(in) :: PLON                 ! I [nbr] loop size
    integer, intent(in) :: PLOND                ! I [nbr] dimension size
    logical,intent(in)::flg_mbl(plond)          ! I [flg] Mobilization candidate flag
    real,intent(in)::dns_mdp(plond)             ! I [kg m-3] Midlayer density
    real,intent(in)::wnd_frc(plond)         ! I [m s-1] Surface friction velocity
    real,intent(in)::wnd_frc_thr_slt(plond) ! I [m s-1] Threshold friction speed for saltation
    ! Output
    real,intent(out)::flx_mss_hrz_slt_ttl(plond) ! O [kg m-1 s-1] Vertically integrated streamwise mass flux
    integer rcd ! [rcd] Return success code
    ! Local
    real wnd_frc_rat ! [frc] Ratio of wind friction threshold to wind friction
    real grv_sfc     ! [m s-2] Gravitation constant
    integer lon_idx ! [idx] Counting index for lon
    integer pdp_idp ! [idx] Counting index for lon
    REAL(KIND=JPRB) :: ZHOOK_HANDLE 

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FLX_MSS_HRZ_SLT_GET_MM',0,ZHOOK_HANDLE)  
    ! Main Code
    
    !++alfgr Set gravitation.
    grv_sfc = XG

    ! Initialize output
    flx_mss_hrz_slt_ttl(:)=0.0d0 ! [kg m-1 s-1]
    
    do lon_idx=1,plon
       if (flg_mbl(lon_idx).and.wnd_frc(lon_idx) > wnd_frc_thr_slt(lon_idx)) then
          wnd_frc_rat=wnd_frc_thr_slt(lon_idx)/wnd_frc(lon_idx) ! [frc]
          flx_mss_hrz_slt_ttl(lon_idx)= & ! [kg m-1 s-1] 
               dns_mdp(lon_idx)*(wnd_frc(lon_idx)**3.0d0)* &
               (1.0d0-wnd_frc_rat)*(1.0d0+wnd_frc_rat)*(1.0d0+wnd_frc_rat)/grv_sfc ! Whi79 p. 4648 (19), MaB97 p. 16422 (28)
       endif ! endif 
    end do ! end loop over lon
    
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FLX_MSS_HRZ_SLT_GET_MM',1,ZHOOK_HANDLE)      
  end subroutine flx_mss_hrz_slt_get_mm ! end flx_mss_hrz_slt_ttl_Whi79_get()
!  
  subroutine flx_mss_vrt_dst_ttl_MaB95_get( &
       dst_slt_flx_rat_ttl, & ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
       flg_mbl, & ! I [flg] Mobilization candidate flag
       flx_mss_hrz_slt_ttl, & ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
       flx_mss_vrt_dst_ttl, & ! O [kg m-2 s-1] Total vertical mass flux of dust
       mss_frc_cly,    & ! I [frc] Mass fraction clay 
       PLON,           & ! I [nbr] loop size
       PLOND)            ! I [nbr] dimension size
    
    ! Purpose: Diagnose total vertical mass flux of dust from vertically integrated streamwise mass flux
    ! Theory: Uses clay-based method proposed by Marticorena & Bergametti (1995)
    ! Their parameterization is based only on data for mss_frc_cly < 0.20
    ! For clayier soils, dst_slt_flx_rat_ttl may behave dramatically differently
    ! Whether this behavior changes when mss_frc_cly > 0.20 is unknown
    ! Anecdotal evidence suggests vertical flux decreases for mss_frc_cly > 0.20
    ! Thus we use min[mss_frc_cly,0.20] in MaB95 parameterization
    implicit none
    ! Parameters
    ! Input
    integer, intent(in) :: plon        ! I [nbr] loop size
    integer, intent(in) :: plond       ! I [nbr] dimension size
    logical,intent(in)::flg_mbl(plond) ! I [flg] Mobilization candidate flag
    real,intent(in)::flx_mss_hrz_slt_ttl(plond) ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
    real,intent(in)::mss_frc_cly(plond) ! I [frc] Mass fraction clay 
    ! Output
    real,intent(out)::dst_slt_flx_rat_ttl(plond) ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
    real,intent(out)::flx_mss_vrt_dst_ttl(plond) ! O [kg m-2 s-1] Total vertical mass flux of dust
    integer rcd ! [rcd] Return success code
    ! Local
    integer lon_idx ! [idx] Counting index for lon
    real mss_frc_cly_vld ! [frc] Mass fraction clay limited to 0.20
    real ln10 ! [frc] Natural log of 10
    REAL(KIND=JPRB) :: ZHOOK_HANDLE 

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FLX_MSS_VRT_DST_TTL_MAB95_GET',0,ZHOOK_HANDLE)      
    ! Initialize some variables
    ln10=log(10.0d0) ! [frc] Natural log of 10
    do lon_idx=1,plon
       if (flg_mbl(lon_idx)) then
          ! 19990603: fxm: Dust production is EXTREMELY sensitive to this parameter, which changes flux by 3 orders of magnitude in 0.0 < mss_frc_cly < 0.20
          mss_frc_cly_vld=min(mss_frc_cly(lon_idx),0.2) ! [frc]
          dst_slt_flx_rat_ttl(lon_idx)= & ! [m-1]
               100.0d0*exp(ln10*(13.4d0*mss_frc_cly_vld-6.0d0)) ! MaB95 p. 16423 (47)
          flx_mss_vrt_dst_ttl(lon_idx)=flx_mss_hrz_slt_ttl(lon_idx)*dst_slt_flx_rat_ttl(lon_idx) ! [kg m-1 s-1] 
       endif ! endif flg_mbl
    end do ! end loop over lon
    
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FLX_MSS_VRT_DST_TTL_MAB95_GET',1,ZHOOK_HANDLE)     
  end subroutine flx_mss_vrt_dst_ttl_MaB95_get ! end flx_mss_vrt_dst_ttl_MaB95_get()
!
  subroutine flx_mss_vrt_dst_Aust_get_mm( &  ! Shao and al. 1993
       dst_slt_flx_rat_ttl, & ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
       flg_mbl, & ! I [flg] Mobilization candidate flag
       flx_mss_hrz_slt_ttl, & ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
       flx_mss_vrt_dst_ttl, & ! O [kg m-2 s-1] Total vertical mass flux of dust
       dns_mdp, & ! I [kg m-3] Midlayer density 
       wnd_frc_thr_slt, & ! I [m s-1] Threshold friction speed for saltation       
       PLON,            & ! I [nbr] loop size
       PLOND            & ! I [nbr] dimension size
       )               ! I [nbr] dimension size
    
    ! Purpose: Diagnose total vertical mass flux of dust from vertically integrated streamwise mass flux
 !    
    use MODD_CSTS, only : XG  ! Gravitation constant
!    
    implicit none
    ! Parameters
    real,parameter::dmt_slt_opt=75.0e-3 ! [mm] 
    real,parameter::dmt_ero_opt=6.7e-3  ! [mm] 
    real,parameter::dns_prt_sfc=2650.0 ! [kg m-3] Dry density of soil particles (excluding pores)
    real,parameter::gama =2.5          ! dimensionless speed factor
    
    ! Input
    integer, intent(in) :: plon        ! I [nbr] loop size
    integer, intent(in) :: plond       ! I [nbr] dimension size
    logical,intent(in)::flg_mbl(plond) ! I [flg] Mobilization candidate flag
    real,intent(in)::flx_mss_hrz_slt_ttl(plond) ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
    real,intent(in)::dns_mdp(plond) ! I [kg m-3] Midlayer density 
    real,intent(in)::wnd_frc_thr_slt(plond) ! II [m s-1] Threshold friction speed for saltation
    !
    ! Output
    real,intent(out)::dst_slt_flx_rat_ttl(plond) ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
    real,intent(out)::flx_mss_vrt_dst_ttl(plond) ! O [kg m-2 s-1] Total vertical mass flux of dust
    integer rcd ! [rcd] Return success code
    ! Local
    integer lon_idx ! [idx] Counting index for lon
    integer pdp_idp ! [idx] Counting index for pdp
    
    real expdd ! 
    real lnds  ! 
    real beta  ! 
    real bgxg  ! 
    
    real wnd_frc_rat ! [frc] Ratio of wind friction threshold to wind friction
    real grv_sfc     ! [m s-2] Gravitation constant
    real :: rop_roa(PLOND) ! [frc] ratio of particle densite to Midlayer density 
    REAL(KIND=JPRB) :: ZHOOK_HANDLE 

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FLX_MSS_VRT_DST_AUST_GET_MM',0,ZHOOK_HANDLE)  
   
    ! Initialize some variables
    expdd=exp(-140.7d0*dmt_ero_opt + 0.37d0) 
    lnds = (0.328*1.0d-4) + (0.125d0*1.0d-4*log(dmt_slt_opt))    
    beta =  lnds * expdd   ! [mm]  
    grv_sfc = XG
    bgxg = beta * gama* grv_sfc
!   
    do lon_idx=1,plon
       if (flg_mbl(lon_idx)) then
          rop_roa(lon_idx)=dns_prt_sfc/dns_mdp(lon_idx)
         dst_slt_flx_rat_ttl(lon_idx)= & ! [mm/m]
          (2.0d0/3.0d0)*bgxg*rop_roa(lon_idx)/(wnd_frc_thr_slt(lon_idx)**2.0d0)
          flx_mss_vrt_dst_ttl(lon_idx)=1.0d-3 * &
	  flx_mss_hrz_slt_ttl(lon_idx)*dst_slt_flx_rat_ttl(lon_idx) ! [kg m-1 s-1] 
       endif ! endif flg_mbl
    end do ! end loop over lon
!    
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:FLX_MSS_VRT_DST_AUST_GET_MM',1,ZHOOK_HANDLE) 
  end subroutine flx_mss_vrt_dst_Aust_get_mm ! end flx_mss_vrt_dst_ttl_Aust_get()


!
!&&&&&&&&&&&& 
  subroutine distribution(    &
        KSIZE                 &  !I loop over lon 
       ,NDP                   &  !I loop over diameter
       ,NTEX                  &  !I loop over texture 
       ,NMODE                 &  !I loop over mode 
       ,PDLNDP                &  !I [?m] delta Dp
       ,PTEXT                 &  !I texture
       ,PDSRLV                &  !O surface relative occup?e par chaque particule
       ,PZS0)                     !O rugosit? de lka surface lisse  
    USE MODD_CSTS, only: XPI

    implicit none
    
    !INPUT, set their dimensions to their passed lengths or to KSIZE ?
    integer, intent(in)                  :: KSIZE    ![nbr] length of passed arrays
    integer, intent(in)                  :: NDP      ![nbr] loop over diameter
    integer, intent(in)                  :: NTEX     ![nbr] loop over texture
    integer, intent(in)                  :: NMODE    ![nbr] loop over mode
    real,    intent(in)                  :: PDLNDP   ! delta Dp [?m]
    integer, intent(in), dimension(KSIZE)   :: PTEXT    !texture
    !OUTPUT surface relative
    real, intent(out), dimension(NTEX,NDP)  :: PDSRLV   ! [--] Output surface relative
    real, intent(out), dimension(NTEX)  :: PZS0         ! [m] Output rugosit? de la surface lisse (Dp/30)
    
    !OUTPUT wind friction speed after modified for saltation feedbacks
!============
!

    real,parameter::dns_slt=2650.d0    ! [kg m-3] Density of optimal saltation particles, MBA97 p. 4388 
    !Define local variables:
    real   :: ZSIGMA(NMODE,NTEX)                      !d?viation g?om?trique standard du mode
    real   :: ZMFRAC(NMODE,NTEX)                    !fraction massique du mode
    real   :: ZDMED(NMODE,NTEX)                       ! diametre median du mode
    real   :: ZDELDLN(NDP)
    real   :: ZDPLN(NDP)
    real   :: ZDS(NDP,NTEX)
    real   :: ZSP(NDP)
    real   :: ZDM(NDP,NTEX)
    real   :: ZDMLN(NDP)
    real   :: ZDP(NDP)
    
!
    integer             :: i                             !Counter for number of points (used in loops)
    integer             :: JMOD                          !Counter for number of mode
    integer             :: IDP                           !Counter for number of particle
    integer             :: ITEX                          !Counter for number of texture
    integer             :: IMASK                         !Counter NDP*KSIZE
    
    real                :: DM1, DM2
    REAL(KIND=JPRB) :: ZHOOK_HANDLE 

    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:DISTRIBUTION',0,ZHOOK_HANDLE)  
        
!================
do ITEX = 1 , NTEX

  if (ITEX==1) then      !Sand

        ZMFRAC(1,ITEX)=0.0
 	ZMFRAC(2,ITEX)=0.1
	ZMFRAC(3,ITEX)=0.9
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.7
	ZSIGMA(3,ITEX)=1.6
!	
	ZDMED(1,ITEX)=10.0
	ZDMED(2,ITEX)=100.0
	ZDMED(3,ITEX)=1000.0
!
	PZS0(ITEX) =1E-6*ZDMED(3,ITEX)/30.0
!	
  elseif (ITEX== 2) then      !Loamy Sand

        ZMFRAC(1,ITEX)=0.1
 	ZMFRAC(2,ITEX)=0.3
	ZMFRAC(3,ITEX)=0.6
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.7
	ZSIGMA(3,ITEX)=1.6
!	
	ZDMED(1,ITEX)=10.0
	ZDMED(2,ITEX)=100.0
	ZDMED(3,ITEX)=690.0
!	
	PZS0(ITEX) = 1E-6*ZDMED(3,ITEX)/30.0
!	
  elseif (ITEX== 3) then      !Sandy Loam

        ZMFRAC(1,ITEX)=0.1
 	ZMFRAC(2,ITEX)=0.3
	ZMFRAC(3,ITEX)=0.6
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.7
	ZSIGMA(3,ITEX)=1.6
!	
	ZDMED(1,ITEX)=5.0
	ZDMED(2,ITEX)=100.0
	ZDMED(3,ITEX)=520.0
!	
	PZS0(ITEX) = 1E-6*ZDMED(3,ITEX)/30.0
!	
  elseif (ITEX == 4) then     !Silt Loam

        ZMFRAC(1,ITEX)=0.15
 	ZMFRAC(2,ITEX)=0.35
	ZMFRAC(3,ITEX)=0.50
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.7
	ZSIGMA(3,ITEX)=1.6
!	
	ZDMED(1,ITEX)=5.0
	ZDMED(2,ITEX)=100.0
	ZDMED(3,ITEX)=520.0
!	
	PZS0(ITEX) = 1E-6*ZDMED(3,ITEX)/30.0
!	
  elseif (ITEX== 5) then     !Loam

        ZMFRAC(1,ITEX)=0.15
 	ZMFRAC(2,ITEX)=0.50
	ZMFRAC(3,ITEX)=0.35
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.7
	ZSIGMA(3,ITEX)=1.6
!	
	ZDMED(1,ITEX)=2.5
	ZDMED(2,ITEX)=75.0
	ZDMED(3,ITEX)=520.0
!	
	PZS0(ITEX) = 1E-6*ZDMED(3,ITEX)/30.0
!	
  elseif (ITEX== 6) then     !Sandy Clay Loam

        ZMFRAC(1,ITEX)=0.20
 	ZMFRAC(2,ITEX)=0.50
	ZMFRAC(3,ITEX)=0.30
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.7
	ZSIGMA(3,ITEX)=1.7
!	
	ZDMED(1,ITEX)=2.5
	ZDMED(2,ITEX)=75.0
	ZDMED(3,ITEX)=210.0
!	
	PZS0(ITEX) = 1E-6*ZDMED(3,ITEX)/30.0
!	
  elseif (ITEX== 7) then     !Silty Clay Loam

        ZMFRAC(1,ITEX)=0.20
 	ZMFRAC(2,ITEX)=0.50
	ZMFRAC(3,ITEX)=0.30
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.7
	ZSIGMA(3,ITEX)=1.7
!	
	ZDMED(1,ITEX)=2.5
	ZDMED(2,ITEX)=50.0
	ZDMED(3,ITEX)=210.0
!	
	PZS0(ITEX) = 1E-6*ZDMED(3,ITEX)/30.0
!	
  elseif (ITEX== 8) then     !Clay Loam

        ZMFRAC(1,ITEX)=0.30
 	ZMFRAC(2,ITEX)=0.50
	ZMFRAC(3,ITEX)=0.20
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.7
	ZSIGMA(3,ITEX)=1.7
!	
	ZDMED(1,ITEX)=1.0
	ZDMED(2,ITEX)=50.0
	ZDMED(3,ITEX)=125.0
!	
	PZS0(ITEX) = 1E-6*ZDMED(3,ITEX)/30.0
!	
  elseif (ITEX== 9) then     !Sandy Clay

        ZMFRAC(1,ITEX)=0.35
 	ZMFRAC(2,ITEX)=0.00
	ZMFRAC(3,ITEX)=0.65
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.8
	ZSIGMA(3,ITEX)=1.8
!	
	ZDMED(1,ITEX)=1.0
	ZDMED(2,ITEX)=10.0
	ZDMED(3,ITEX)=100.0
!	
	PZS0(ITEX) = 1E-6*ZDMED(3,ITEX)/30.0
!	
  elseif (ITEX== 10) then    !Silty Clay

        ZMFRAC(1,ITEX)=0.40
 	ZMFRAC(2,ITEX)=0.00
	ZMFRAC(3,ITEX)=0.60
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.8
	ZSIGMA(3,ITEX)=1.8
!	
	ZDMED(1,ITEX)=0.5
	ZDMED(2,ITEX)=10.0
	ZDMED(3,ITEX)=100.0
!	
	PZS0(ITEX) = 1E-6*ZDMED(3,ITEX)/30.0
!	
  elseif (ITEX== 11) then     !Clay

        ZMFRAC(1,ITEX)=0.50
 	ZMFRAC(2,ITEX)=0.00
	ZMFRAC(3,ITEX)=0.50
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.8
	ZSIGMA(3,ITEX)=1.8
!	
	ZDMED(1,ITEX)=0.5
	ZDMED(2,ITEX)=10.0
	ZDMED(3,ITEX)=100.0
!	
	PZS0(ITEX) = 1E-6*ZDMED(3,ITEX)/30.0
!	
  else                        !Silt

        ZMFRAC(1,ITEX)=0.15
 	ZMFRAC(2,ITEX)=0.4
	ZMFRAC(3,ITEX)=0.45
!	
	ZSIGMA(1,ITEX)=1.8
	ZSIGMA(2,ITEX)=1.7
	ZSIGMA(3,ITEX)=1.6
!	
	ZDMED(1,ITEX)=2.5
	ZDMED(2,ITEX)=75.0
	ZDMED(3,ITEX)=520.0
!	
	PZS0(ITEX) = 1E-6*ZDMED(3,ITEX)/30.0
!	
  endif  
end do
!
! calcul pour chaque particule 
! Boucle NDP 
! Initialisation
! surface basale totale
 ZSP(:) = 0.d0
!
! Initialisation pour IDP = 1
 
 ZDP(1) = 0.1                     ! ?m
 ZDPLN(1) = log(ZDP(1))           ! ?m
! 
DO IDP=2,NDP
    ZDPLN(IDP) = ZDPLN(IDP-1) + PDLNDP
    ZDP(IDP)=exp(ZDPLN(IDP))
END DO
!
! Calcul de la distribution massique des particules (dm(Dp)/dln(Dp)) pour chaque particule de diametre Dp

do ITEX = 1, NTEX
!
  DO IDP=1,NDP
    ZDMLN(IDP)=0.0d0
    DO JMOD=1, NMODE        
      DM1=ZMFRAC(JMOD,ITEX)/(sqrt(2.d0*XPI)*log(ZSIGMA(JMOD,ITEX)))
      DM2=exp((log(ZDP(IDP))-log(ZDMED(JMOD,ITEX)))**2./(-2.d0*(log(ZSIGMA(JMOD,ITEX)))**2.d0))
      ZDMLN(IDP)= ZDMLN(IDP) + DM1*DM2
    END DO
  ZDM (IDP,ITEX) = ZDMLN(IDP) * PDLNDP  
! Calcul de la surface basale occup?e par chaque particule de diam?tre Dp  
  ZDS(IDP,ITEX)= 3.d0*ZDM(IDP,ITEX)/(2.*dns_slt*ZDP(IDP)) 
              
! Calcul de la surface basale totale Sp= (somme (ds(Dp)):
  ZSP(ITEX) = ZDS(IDP,ITEX) + ZSP(ITEX) 
       
  ENDDO
!  
end do       ! end loop over ntex
!
! Calcul de la distribution de surface relative des particules
DO IDP = 1, NDP
  do ITEX = 1, NTEX
    PDSRLV(ITEX,IDP) = ZDS(IDP,ITEX)/ZSP(ITEX)
  end do  
END DO 
! 
    IF (LHOOK) CALL DR_HOOK('MODE_DSTMBLUTL:DISTRIBUTION',1,ZHOOK_HANDLE)  
  end subroutine distribution
!    
END MODULE MODE_DSTMBLUTL
