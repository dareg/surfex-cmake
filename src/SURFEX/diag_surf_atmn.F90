!     #########
SUBROUTINE DIAG_SURF_ATM_n (DGEI, DGEIC, DGF, DGFC, DGL, DGLC, DGI, DGIC, DGS, DGSC, &
                            DGU, DGUP, DGUC, DGUPC, DGT, DGW, DGWC, U, USS, &
                            HPROGRAM)
!     #################################################################################
!
!!****  *DIAG_SURF_ATM_n * - Chooses the surface schemes for diagnostics
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified    01/2006 : sea flux parameterization.
!!      Modified    08/2008 : cumulated fluxes
!       B. decharme 04/2013 : Add EVAP and SUBL diag
!!------------------------------------------------------------------
!

!
!
!
!
!
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_n, ONLY : DIAG_t, DIAG_PATCH_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
!
USE MODD_SURF_CONF,      ONLY : CPROGNAME
USE MODD_DATA_COVER_PAR, ONLY : NTILESFC
!
USE MODI_DIAG_NATURE_n 
USE MODI_DIAG_SEA_n 
USE MODI_DIAG_INLAND_WATER_n 
USE MODI_DIAG_TOWN_n 
USE MODI_AVERAGE_DIAG
!
USE MODI_MINZS_VERT_SHIFT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIC
TYPE(DIAG_t), INTENT(INOUT) :: DGF
TYPE(DIAG_t), INTENT(INOUT) :: DGFC
TYPE(DIAG_t), INTENT(INOUT) :: DGL
TYPE(DIAG_t), INTENT(INOUT) :: DGLC
TYPE(DIAG_t), INTENT(INOUT) :: DGI
TYPE(DIAG_t), INTENT(INOUT) :: DGIC
TYPE(DIAG_t), INTENT(INOUT) :: DGS
TYPE(DIAG_t), INTENT(INOUT) :: DGSC
TYPE(DIAG_t), INTENT(INOUT) :: DGU
TYPE(DIAG_PATCH_t), INTENT(INOUT) :: DGUP
TYPE(DIAG_t), INTENT(INOUT) :: DGUC
TYPE(DIAG_PATCH_t), INTENT(INOUT) :: DGUPC
TYPE(DIAG_t), INTENT(INOUT) :: DGT
TYPE(DIAG_t), INTENT(INOUT) :: DGW
TYPE(DIAG_t), INTENT(INOUT) :: DGWC
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
!
!
!*      0.2    declarations of local variables
!
INTEGER :: JTILE                        ! loop on type of surface
LOGICAL :: GNATURE, GTOWN, GWATER, GSEA ! .T. if the corresponding surface is represented
INTEGER :: JSW                          ! number of spectral whort wave bands
!
REAL, DIMENSION(SIZE(U%XSEA),NTILESFC) :: ZFRAC_TILE! fraction of each tile
INTEGER, DIMENSION(5) :: IFACT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
! Preliminaries: Tile related operations
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_N',0,ZHOOK_HANDLE)
CPROGNAME = HPROGRAM
!
! FLAGS for the various surfaces:
!
GSEA      = U%NDIM_SEA    >0
GWATER    = U%NDIM_WATER  >0
GTOWN     = U%NDIM_TOWN   >0
GNATURE   = U%NDIM_NATURE >0
!
! Tile counter:
!
JTILE     = 0 
!
! Fractions for each tile:
!
ZFRAC_TILE(:,:)    = 0.0
!
! Number of spectral short wave bands for detailed radiation budget
JSW = SIZE(DGUP%AL(1)%XSWBD,2)
!
!
CALL GET_DIMS(IFACT)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! SEA Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
! first, pack vector...then call ALMA routine
!
JTILE               = JTILE + 1
!
IF(GSEA)THEN
! 
  ZFRAC_TILE(:,JTILE) = U%XSEA(:)
!
  CALL TREAT_SURF(JTILE,U%NSIZE_SEA,U%NR_SEA,IFACT)
!
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! INLAND WATER Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE               = JTILE + 1
!
IF(GWATER)THEN
!
  ZFRAC_TILE(:,JTILE) = U%XWATER(:)
!
  CALL TREAT_SURF(JTILE,U%NSIZE_WATER,U%NR_WATER,IFACT)
!
ENDIF 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! NATURAL SURFACE Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE               = JTILE + 1
!
IF(GNATURE)THEN
!
    ZFRAC_TILE(:,JTILE) = U%XNATURE(:)
!
  CALL TREAT_SURF(JTILE,U%NSIZE_NATURE,U%NR_NATURE,IFACT)  
!
ENDIF 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! URBAN Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE               = JTILE + 1
!
IF(GTOWN)THEN
!
    ZFRAC_TILE(:,JTILE) = U%XTOWN(:)
!
  CALL TREAT_SURF(JTILE,U%NSIZE_TOWN,U%NR_TOWN,IFACT)  
!
ENDIF 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Grid box average fluxes/properties:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
CALL AVERAGE_DIAG(ZFRAC_TILE, DGU, DGUP, DGUC, DGUPC)              
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Quantities at 2 meters above the minimum orography of the grid mesh
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF (DGU%L2M_MIN_ZS) CALL GET_2M
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_N',1,ZHOOK_HANDLE)
CONTAINS
!=======================================================================================
SUBROUTINE GET_2M
!
REAL, DIMENSION(SIZE(U%XSEA)) :: ZPS         ! surface air pressure
REAL, DIMENSION(SIZE(U%XSEA)) :: ZRHOA       ! surface air density
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_n:GET_2M',0,ZHOOK_HANDLE)
!
CALL MINZS_VERT_SHIFT(U%XZS,USS%XMIN_ZS,DGU%XT2M,DGU%XQ2M,DGU%XPS,DGU%XRHOA, &
                      DGU%XT2M_MIN_ZS,DGU%XQ2M_MIN_ZS,ZPS,ZRHOA)  
DGU%XHU2M_MIN_ZS = DGU%XHU2M
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_n:GET_2M',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_2M
!
!=======================================================================================
!
SUBROUTINE GET_DIMS(KFACT)
!
IMPLICIT NONE
!
INTEGER, DIMENSION(5), INTENT(OUT) :: KFACT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_n:GET_DIMS',0,ZHOOK_HANDLE)
!
KFACT(:)=0
!
IF (DGU%LSURF_BUDGET) KFACT(1)=1
!
IF (DGU%LSURF_BUDGETC) KFACT(2)=1
!
IF (DGU%N2M>=1) KFACT(3)=1
!
IF (DGU%LCOEF) KFACT(4)=1
!
IF (DGU%LSURF_VARS) KFACT(5)=1
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_n:GET_DIMS',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_DIMS
!
!=======================================================================================
!
SUBROUTINE TREAT_SURF(KTILE,KSIZE,KMASK,KFACT)
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)               :: KTILE
INTEGER, INTENT(IN)               :: KSIZE
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
INTEGER, DIMENSION(5), INTENT(IN) :: KFACT
!
REAL, DIMENSION(KSIZE) :: ZP_TS       ! surface temperature (K)
!
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_RN       ! Net radiation (W/m2)
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_H        ! sensible heat flux (W/m2)
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_LE       ! total latent heat flux (W/m2)
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_LEI      ! sublimation latent heat flux (W/m2)
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_GFLUX    ! storage flux (W/m2)
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_EVAP     ! total evapotranspiration (kg/m2/s)
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_SUBL     ! sublimation (kg/m2/s)
!
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_SWD      ! short wave incoming radiation (W/m2)
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_SWU      ! short wave outgoing radiation (W/m2)
REAL, DIMENSION(KSIZE*KFACT(1),JSW*KFACT(1)) :: ZP_SWBD   ! short wave incoming radiation by spectral band (W/m2)
REAL, DIMENSION(KSIZE*KFACT(1),JSW*KFACT(1)) :: ZP_SWBU   ! short wave outgoing radiation by spectral band (W/m2)
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_LWD      ! long wave incoming radiation (W/m2)
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_LWU      ! long wave outgoing radiation (W/m2)
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_FMU      ! zonal friction
REAL, DIMENSION(KSIZE*KFACT(1)) :: ZP_FMV      ! meridian friction 
!
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_RNC      ! Cumulated Net radiation (W/m2)
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_HC       ! Cumulated sensible heat flux (W/m2)
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_LEC      ! Cumulated total latent heat flux (W/m2)
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_LEIC     ! Cumulated sublimation latent heat flux (W/m2)
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_GFLUXC   ! Cumulated storage flux (W/m2)
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_EVAPC    ! Cumulated total evapotranspiration (kg/m2)
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_SUBLC    ! Cumulated sublimation (kg/m2)
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_SWDC     ! Cumulated short wave incoming radiation (W/m2)
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_SWUC     ! Cumulated short wave outgoing radiation (W/m2)
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_LWDC     ! Cumulated long wave incoming radiation (W/m2)
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_LWUC     ! Cumulated long wave outgoing radiation (W/m2)
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_FMUC     ! Cumulated zonal friction
REAL, DIMENSION(KSIZE*KFACT(2)) :: ZP_FMVC     ! Cumulated meridian friction 
!
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_RI       ! Richardson number
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_T2M      ! air temperature at 2 meters (K)
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_T2M_MIN  ! Minimum air temperature at 2 meters (K)
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_T2M_MAX  ! Maximum air temperature at 2 meters (K)
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_Q2M      ! air humidity at 2 meters (kg/kg)
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_HU2M     ! air relative humidity at 2 meters (-)
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_HU2M_MIN ! Minimum air relative humidity at 2 meters (-)
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_HU2M_MAX ! Maximum air relative humidity at 2 meters (-)
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_ZON10M   ! zonal wind at 10 meters (m/s)
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_MER10M   ! meridian wind at 10 meters (m/s)
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_WIND10M  ! wind at 10 meters (m/s)
REAL, DIMENSION(KSIZE*KFACT(3)) :: ZP_WIND10M_MAX ! Maximum wind at 10 meters (m/s)
!
REAL, DIMENSION(KSIZE*KFACT(4)) :: ZP_CD       ! drag coefficient for wind
REAL, DIMENSION(KSIZE*KFACT(4)) :: ZP_CH       ! drag coefficient for heat
REAL, DIMENSION(KSIZE*KFACT(4)) :: ZP_CE       ! drag coefficient for evaporation
REAL, DIMENSION(KSIZE*KFACT(4)) :: ZP_Z0       ! roughness length for momentum
REAL, DIMENSION(KSIZE*KFACT(4)) :: ZP_Z0H      ! roughness length for heat
!
REAL, DIMENSION(KSIZE*KFACT(5)) :: ZP_QS       ! specific humidity
!
INTEGER :: JJ, JJSW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_n:TREAT_SURF',0,ZHOOK_HANDLE)
!
IF (KTILE==1) THEN
  !
  CALL DIAG_SEA_n(DGL, DGLC, DGS, DGSC, U, &
                  HPROGRAM,                             &
                  ZP_RN, ZP_H, ZP_LE, ZP_LEI, ZP_GFLUX, &
                  ZP_RI, ZP_CD, ZP_CH, ZP_CE,           &
                  ZP_QS, ZP_Z0, ZP_Z0H,                 &
                  ZP_T2M, ZP_TS, ZP_Q2M, ZP_HU2M,       &
                  ZP_ZON10M, ZP_MER10M,                 &
                  ZP_SWD, ZP_SWU, ZP_SWBD, ZP_SWBU,     &
                  ZP_LWD, ZP_LWU, ZP_FMU, ZP_FMV,       &
                  ZP_RNC, ZP_HC, ZP_LEC, ZP_GFLUXC,     &
                  ZP_SWDC, ZP_SWUC, ZP_LWDC, ZP_LWUC,   &
                  ZP_FMUC, ZP_FMVC, ZP_T2M_MIN,         &
                  ZP_T2M_MAX, ZP_LEIC, ZP_HU2M_MIN,     &
                  ZP_HU2M_MAX, ZP_WIND10M,              &
                  ZP_WIND10M_MAX,                       &
                  ZP_EVAP, ZP_EVAPC, ZP_SUBL, ZP_SUBLC  )   
  !
ELSEIF (KTILE==2) THEN
  !
  CALL DIAG_INLAND_WATER_n(DGF, DGFC, DGL, DGLC, DGW, DGWC, U, &
                           HPROGRAM,                            &
                           ZP_RN, ZP_H, ZP_LE, ZP_LEI, ZP_GFLUX,&
                           ZP_RI, ZP_CD, ZP_CH, ZP_CE,          &
                           ZP_QS, ZP_Z0, ZP_Z0H,                &
                           ZP_T2M, ZP_TS, ZP_Q2M, ZP_HU2M,      &
                           ZP_ZON10M, ZP_MER10M,                &
                           ZP_SWD, ZP_SWU, ZP_SWBD, ZP_SWBU,    &
                           ZP_LWD, ZP_LWU, ZP_FMU, ZP_FMV,      &
                           ZP_RNC, ZP_HC, ZP_LEC, ZP_GFLUXC,    &
                           ZP_SWDC, ZP_SWUC, ZP_LWDC, ZP_LWUC,  &
                           ZP_FMUC, ZP_FMVC, ZP_T2M_MIN,        &
                           ZP_T2M_MAX, ZP_LEIC, ZP_HU2M_MIN,    &
                           ZP_HU2M_MAX, ZP_WIND10M,             &
                           ZP_WIND10M_MAX,                      &
                           ZP_EVAP, ZP_EVAPC, ZP_SUBL, ZP_SUBLC  )   
  !
ELSEIF (KTILE==3) THEN
  !
  CALL DIAG_NATURE_n(DGEI, DGEIC, DGL, DGLC, DGI, DGIC, U, &
                     HPROGRAM,                            &
                     ZP_RN, ZP_H, ZP_LE, ZP_LEI, ZP_GFLUX,&
                     ZP_RI, ZP_CD, ZP_CH, ZP_CE,          &
                     ZP_QS, ZP_Z0, ZP_Z0H,                &
                     ZP_T2M, ZP_TS, ZP_Q2M, ZP_HU2M,      &
                     ZP_ZON10M, ZP_MER10M,                &
                     ZP_SWD, ZP_SWU, ZP_SWBD, ZP_SWBU,    &
                     ZP_LWD, ZP_LWU, ZP_FMU, ZP_FMV,      &
                     ZP_RNC, ZP_HC, ZP_LEC, ZP_GFLUXC,    &
                     ZP_SWDC, ZP_SWUC, ZP_LWDC, ZP_LWUC,  &
                     ZP_FMUC, ZP_FMVC, ZP_T2M_MIN,        &
                     ZP_T2M_MAX, ZP_LEIC, ZP_HU2M_MIN,    &
                     ZP_HU2M_MAX, ZP_WIND10M,             &
                     ZP_WIND10M_MAX,                      &
                     ZP_EVAP, ZP_EVAPC, ZP_SUBL, ZP_SUBLC )   
  !
ELSEIF (KTILE==4) THEN
  !
  CALL DIAG_TOWN_n(DGL, DGLC, DGT, U, &
                   HPROGRAM,                            &
                   ZP_RN, ZP_H, ZP_LE, ZP_LEI, ZP_GFLUX,&
                   ZP_RI, ZP_CD, ZP_CH, ZP_CE,          &
                   ZP_QS, ZP_Z0, ZP_Z0H,                &
                   ZP_T2M, ZP_TS, ZP_Q2M, ZP_HU2M,      &
                   ZP_ZON10M, ZP_MER10M,                &
                   ZP_SWD, ZP_SWU, ZP_SWBD, ZP_SWBU,    &
                   ZP_LWD, ZP_LWU, ZP_FMU, ZP_FMV,      &
                   ZP_RNC, ZP_HC, ZP_LEC, ZP_GFLUXC,    &
                   ZP_SWDC, ZP_SWUC, ZP_LWDC, ZP_LWUC,  &
                   ZP_FMUC, ZP_FMVC, ZP_T2M_MIN,        &
                   ZP_T2M_MAX, ZP_LEIC, ZP_HU2M_MIN,    &
                   ZP_HU2M_MAX, ZP_WIND10M,             &
                   ZP_WIND10M_MAX,                      &
                   ZP_EVAP, ZP_EVAPC, ZP_SUBL, ZP_SUBLC )  
  !
ENDIF
!
!----------------------------------------------------------------------
IF (DGU%LSURF_BUDGET) THEN
  DO JJ=1,KSIZE
   DGUP%AL(KTILE)%XRN      (KMASK(JJ))  = ZP_RN       (JJ)
   DGUP%AL(KTILE)%XH       (KMASK(JJ))  = ZP_H        (JJ)
   DGUP%AL(KTILE)%XLE      (KMASK(JJ))  = ZP_LE       (JJ)
   DGUP%AL(KTILE)%XLEI     (KMASK(JJ))  = ZP_LEI      (JJ)
   DGUP%AL(KTILE)%XGFLUX   (KMASK(JJ))  = ZP_GFLUX    (JJ)
   DGUP%AL(KTILE)%XEVAP    (KMASK(JJ))  = ZP_EVAP     (JJ)
   DGUP%AL(KTILE)%XSUBL    (KMASK(JJ))  = ZP_SUBL     (JJ)
   DGUP%AL(KTILE)%XSWD     (KMASK(JJ))  = ZP_SWD      (JJ)
   DGUP%AL(KTILE)%XSWU     (KMASK(JJ))  = ZP_SWU      (JJ)
   DGUP%AL(KTILE)%XLWD     (KMASK(JJ))  = ZP_LWD      (JJ)
   DGUP%AL(KTILE)%XLWU     (KMASK(JJ))  = ZP_LWU      (JJ)
   DGUP%AL(KTILE)%XFMU     (KMASK(JJ))  = ZP_FMU      (JJ)
   DGUP%AL(KTILE)%XFMV     (KMASK(JJ))  = ZP_FMV      (JJ)
   DO JJSW=1, SIZE(DGUP%AL(KTILE)%XSWBD,2)
      DGUP%AL(KTILE)%XSWBD    (KMASK(JJ),JJSW) = ZP_SWBD     (JJ,JJSW)
      DGUP%AL(KTILE)%XSWBU    (KMASK(JJ),JJSW) = ZP_SWBU     (JJ,JJSW)
   ENDDO
  ENDDO
END IF
!
IF (DGU%LSURF_BUDGETC) THEN
  DO JJ=1,KSIZE
   DGUPC%AL(KTILE)%XRN      (KMASK(JJ))  = ZP_RNC       (JJ)
   DGUPC%AL(KTILE)%XH       (KMASK(JJ))  = ZP_HC        (JJ)
   DGUPC%AL(KTILE)%XLE      (KMASK(JJ))  = ZP_LEC       (JJ)
   DGUPC%AL(KTILE)%XLEI     (KMASK(JJ))  = ZP_LEIC      (JJ)
   DGUPC%AL(KTILE)%XGFLUX   (KMASK(JJ))  = ZP_GFLUXC    (JJ)
   DGUPC%AL(KTILE)%XEVAP    (KMASK(JJ))  = ZP_EVAPC     (JJ)
   DGUPC%AL(KTILE)%XSUBL    (KMASK(JJ))  = ZP_SUBLC     (JJ)
   DGUPC%AL(KTILE)%XSWD     (KMASK(JJ))  = ZP_SWDC      (JJ)
   DGUPC%AL(KTILE)%XSWU     (KMASK(JJ))  = ZP_SWUC      (JJ)
   DGUPC%AL(KTILE)%XLWD     (KMASK(JJ))  = ZP_LWDC      (JJ)
   DGUPC%AL(KTILE)%XLWU     (KMASK(JJ))  = ZP_LWUC      (JJ)
   DGUPC%AL(KTILE)%XFMU     (KMASK(JJ))  = ZP_FMUC      (JJ)
   DGUPC%AL(KTILE)%XFMV     (KMASK(JJ))  = ZP_FMVC      (JJ)
  ENDDO
END IF
!
DO JJ=1,KSIZE
   DGUP%AL(KTILE)%XTS       (KMASK(JJ))  = ZP_TS      (JJ)
ENDDO
!
IF (DGU%N2M>=1) THEN
  DO JJ=1,KSIZE
   DGUP%AL(KTILE)%XRI      (KMASK(JJ))  = ZP_RI       (JJ)
   DGUP%AL(KTILE)%XT2M     (KMASK(JJ))  = ZP_T2M      (JJ)
   DGUP%AL(KTILE)%XT2M_MIN (KMASK(JJ))  = ZP_T2M_MIN  (JJ)
   DGUP%AL(KTILE)%XT2M_MAX (KMASK(JJ))  = ZP_T2M_MAX  (JJ)
   DGUP%AL(KTILE)%XQ2M     (KMASK(JJ))  = ZP_Q2M      (JJ)
   DGUP%AL(KTILE)%XHU2M    (KMASK(JJ))  = ZP_HU2M     (JJ)
   DGUP%AL(KTILE)%XHU2M_MIN(KMASK(JJ))  = ZP_HU2M_MIN (JJ)
   DGUP%AL(KTILE)%XHU2M_MAX(KMASK(JJ))  = ZP_HU2M_MAX (JJ)
   DGUP%AL(KTILE)%XZON10M  (KMASK(JJ))  = ZP_ZON10M   (JJ)
   DGUP%AL(KTILE)%XMER10M  (KMASK(JJ))  = ZP_MER10M   (JJ)
   DGUP%AL(KTILE)%XWIND10M (KMASK(JJ))  = ZP_WIND10M   (JJ)
   DGUP%AL(KTILE)%XWIND10M_MAX (KMASK(JJ))  = ZP_WIND10M_MAX   (JJ)
  ENDDO
END IF
!
IF (DGU%LCOEF) THEN
  DO JJ=1,KSIZE
   DGUP%AL(KTILE)%XCD      (KMASK(JJ))  = ZP_CD       (JJ)
   DGUP%AL(KTILE)%XCH      (KMASK(JJ))  = ZP_CH       (JJ)
   DGUP%AL(KTILE)%XCE      (KMASK(JJ))  = ZP_CE       (JJ)
   DGUP%AL(KTILE)%XQS      (KMASK(JJ))  = ZP_QS       (JJ)
   DGUP%AL(KTILE)%XZ0      (KMASK(JJ))  = ZP_Z0       (JJ)
   DGUP%AL(KTILE)%XZ0H     (KMASK(JJ))  = ZP_Z0H      (JJ)
  ENDDO
END IF
!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_n:TREAT_SURF',1,ZHOOK_HANDLE)
!
END SUBROUTINE TREAT_SURF
!=======================================================================================
END SUBROUTINE DIAG_SURF_ATM_n
