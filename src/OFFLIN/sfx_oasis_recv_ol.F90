!#########
SUBROUTINE SFX_OASIS_RECV_OL (F, I, S, U, W, &
                              HPROGRAM,KI,KSW,PTIMEC,PTSTEP_SURF,   &
                              PZENITH,PSW_BANDS,          &
                              PTSRAD,PDIR_ALB,PSCA_ALB,PEMIS,PTSURF )
!#############################################
!
!!****  *SFX_OASIS_RECV_OL* - Offline driver that receive coupling fields from oasis
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_OFF_SURFEX_n, ONLY : GOTO_MODEL
!
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODI_SFX_OASIS_RECV
USE MODI_PUT_SFXCPL_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef AIX64
!$ USE OMP_LIB
#endif
!
IMPLICIT NONE
!
#ifndef AIX64
!$ INCLUDE 'omp_lib.h'
#endif
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(WATFLUX_t), INTENT(INOUT) :: W
!
CHARACTER(LEN=6),       INTENT(IN)  :: HPROGRAM    ! program calling surf. schemes
!
INTEGER,                INTENT(IN)  :: KI          ! number of points on this proc
INTEGER,                INTENT(IN)  :: KSW         ! number of short-wave spectral bands
REAL,                   INTENT(IN)  :: PTIMEC      ! Cumulated run time step (s)
REAL,                   INTENT(IN)  :: PTSTEP_SURF ! Surfex time step
!
REAL, DIMENSION(KI),    INTENT(IN)  :: PZENITH   ! zenithal angle       (radian from the vertical)
REAL, DIMENSION(KSW),   INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
!
REAL, DIMENSION(KI),    INTENT(OUT) :: PTSRAD    ! radiative temperature                 (K)
REAL, DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each spectral band  (-)
REAL, DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each spectral band (-)
REAL, DIMENSION(KI),    INTENT(OUT) :: PEMIS     ! emissivity                            (-)
REAL, DIMENSION(KI),    INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(KI) :: ZLAND_WTD     ! Land water table depth (m)
REAL, DIMENSION(KI) :: ZLAND_FWTD    ! Land grid-cell fraction of water table rise (-)
REAL, DIMENSION(KI) :: ZLAND_FFLOOD  ! Land Floodplains fraction (-)
REAL, DIMENSION(KI) :: ZLAND_PIFLOOD ! Land Potential flood infiltration (kg/m2)
REAL, DIMENSION(KI) :: ZSEA_SST      ! Sea surface temperature (K)
REAL, DIMENSION(KI) :: ZSEA_UCU      ! Sea u-current stress (Pa)
REAL, DIMENSION(KI) :: ZSEA_VCU      ! Sea v-current stress (Pa)
REAL, DIMENSION(KI) :: ZSEAICE_SIT   ! Sea-ice Temperature (K)
REAL, DIMENSION(KI) :: ZSEAICE_CVR   ! Sea-ice cover (-)
REAL, DIMENSION(KI) :: ZSEAICE_ALB   ! Sea-ice albedo (-)
!
LOGICAL  :: GOASIS_PUT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_RECV_OL',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*       1.     init coupling fields:
!               ----------------------------------
!
ZLAND_WTD    (:) = XUNDEF
ZLAND_FWTD   (:) = XUNDEF
ZLAND_FFLOOD (:) = XUNDEF
ZLAND_PIFLOOD(:) = XUNDEF
!
ZSEA_SST(:) = XUNDEF
ZSEA_UCU(:) = XUNDEF
ZSEA_VCU(:) = XUNDEF
!
ZSEAICE_SIT(:) = XUNDEF
ZSEAICE_CVR(:) = XUNDEF
ZSEAICE_ALB(:) = XUNDEF
!
!
!*       2.     Receive fields to other models proc by proc:
!               --------------------------------------------
!
GOASIS_PUT = .FALSE.
!
CALL SFX_OASIS_RECV(HPROGRAM,KI,KSW,PTIMEC-PTSTEP_SURF,&
                    GOASIS_PUT,                        &
                    ZLAND_WTD    (:),ZLAND_FWTD   (:), &
                    ZLAND_FFLOOD (:),ZLAND_PIFLOOD(:), &
                    ZSEA_SST     (:),ZSEA_UCU     (:), &
                    ZSEA_VCU     (:),ZSEAICE_SIT  (:), &
                    ZSEAICE_CVR  (:),ZSEAICE_ALB  (:)  )
!
!
!*       3.     Put definitions for exchange of coupling fields :
!               -------------------------------------------------
!
IF(GOASIS_PUT)THEN
!
  CALL PUT_SFXCPL_n(F, I, S, U, W, &
                     HPROGRAM,KI,KSW,PSW_BANDS,       &
                     PZENITH      (:),ZLAND_WTD   (:),&
                     ZLAND_FWTD   (:),ZLAND_FFLOOD(:),&
                     ZLAND_PIFLOOD(:),ZSEA_SST    (:),&
                     ZSEA_UCU     (:),ZSEA_VCU    (:),&                     
                     ZSEAICE_SIT  (:),ZSEAICE_CVR (:),&
                     ZSEAICE_ALB  (:),PTSRAD      (:),&
                     PDIR_ALB     (:,:),PSCA_ALB    (:,:),&
                     PEMIS        (:),PTSURF      (:) )
!
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_RECV_OL',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SFX_OASIS_RECV_OL
