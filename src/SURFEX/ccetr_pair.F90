!     ######
!
   SUBROUTINE CCETR_PAIR(KNIV, PABC, PABC_SUP, PSSA_SUP, PSSA_INF,          &
                        PB_SUP, PB_INF, PIA, PXMUS, PLAI, PALB_VEG, PALB_SOIL, &
                        PFD_SKY, PFD_VEG, PTR, PXIA, PLAI_EFF                  )
   
!
!!***	*CCETR_PAIR* ***
!!
!!    PURPOSE
!!    -------
!!    Calculates radiative transfer within the canopy
!!
!!**  METHOD
!!    ------
!!    Carrer et al. 
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------      
!!    USE MODD_CO2V_PAR
!!
!!    REFERENCE
!!    ---------
!!    Carrer et al. ??
!!      
!!    AUTHOR
!!    ------
!!	  D. Carrer           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/04/11 
!!
!-------------------------------------------------------------------------------
!
USE MODD_CO2V_PAR,   ONLY : XK_SUP, XK_INF, XXSI_SUP, XXSI_INF 
USE MODD_CSTS,       ONLY : XI0
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!*       0.     DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER, INTENT(IN)               :: KNIV
REAL,    INTENT(IN)               :: PABC, PABC_SUP
!                              PABC    = abscissa needed for integration
!                                       of net assimilation and stomatal 
!                                       conductance over canopy depth
REAL,    INTENT(IN)               :: PSSA_SUP, PSSA_INF
REAL, DIMENSION(:), INTENT(IN)    :: PB_SUP, PB_INF
REAL, DIMENSION(:), INTENT(IN)    :: PIA, PXMUS, PLAI
!	                       PIA   = absorbed PAR / PIR
!                              PXMUS = cosine of solar zenith angle
!                              PLAI  = leaf area index
REAL, DIMENSION(:), INTENT(IN)    :: PALB_VEG, PALB_SOIL
REAL, DIMENSION(:), INTENT(IN)    :: PFD_SKY
REAL, DIMENSION(:), INTENT(INOUT) :: PFD_VEG, PTR
REAL, DIMENSION(:), INTENT(OUT)   :: PXIA
!                              PXIA  = abs. radiation of veg
REAL, DIMENSION(:), INTENT(OUT)   :: PLAI_EFF
!
!*      0.2    declarations of local variables

!
REAL, DIMENSION(SIZE(PLAI)) :: ZSLAI_TRU, ZB_DR, ZOMEGA_DR, ZOMEGA_DF, &
                               ZFD_VEG, ZTDF, ZIDR, ZIDF
!                              ZIDF  = interception of diffusion
!                              ZIDR  = direct interception
!                              XB_DR = DH albedo of upper/lower layers
REAL    :: ZB_DF
REAL    :: ZSUP, ZINF, ZGT_SUP, ZGT_INF, ZGT
INTEGER :: I
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CCETR_PAIR',0,ZHOOK_HANDLE)
!
ZSUP = - 0.461 * XXSI_SUP + 3.8
ZINF = - 0.461 * XXSI_INF + 3.8
!
PLAI_EFF(:) = 0.
!
!Angular projection of the leaves 
!  0.5                    : spherical distribution 
!  (2./!PI)*sin(zs*!Dtor) : vertical distribution
!  cos(zs*!Dtor)          : horizontal distribution
ZGT_SUP = 0.5 
ZGT_INF = 0.5
!
IF (PABC.GT.0.8) THEN
  ZGT = ZGT_SUP
  CALL SET_PARAM(PSSA_SUP, XK_SUP, ZSUP, PB_SUP, ZOMEGA_DR, ZB_DR, ZOMEGA_DF, ZB_DF)
ELSE
  ZGT = ZGT_INF
  CALL SET_PARAM(PSSA_INF, XK_INF, ZINF, PB_INF, ZOMEGA_DR, ZB_DR, ZOMEGA_DF, ZB_DF)
ENDIF
!
!
DO I=1,SIZE(PIA)
  IF (PIA(I)>0.) THEN
    ZSLAI_TRU(I) = (PABC_SUP-PABC)*PLAI(I)
    PLAI_EFF(I) = ZOMEGA_DR(I)*ZSLAI_TRU(I)
    ! transmittance of direct beam
    ZIDR(I) = EXP(-ZGT*ZB_DR(I)*ZOMEGA_DR(I)*ZSLAI_TRU(I)/PXMUS(I))
    ! transmittance of diffuse beam
    ZIDF(I) = EXP(-ZB_DF*ZOMEGA_DF(I)*ZSLAI_TRU(I))
    !
    PTR(I) = ((1.-PFD_VEG(I))*ZIDR(I) + PFD_VEG(I)*ZIDF(I))*PTR(I)
  ENDIF
ENDDO
!
!
IF (PABC.GT.0.8) THEN
  DO I=1,SIZE(PIA)
    IF (PIA(I)>0.) THEN
      ! diffuse fraction due to vegetation
      ZFD_VEG(I) = EXP(-(1.-PABC)*ZOMEGA_DR(I)*PLAI(I))
      ZFD_VEG(I) = (1. - ZFD_VEG(I)) / (1. - (1.-PXMUS(I))*ZFD_VEG(I))
      PFD_VEG(I) = MIN(ZFD_VEG(I) + PFD_SKY(I),1.)
    ENDIF
  ENDDO
ENDIF
!
! transmissivity of upper layers
!
PXIA(:) = 0.
WHERE (PIA(:)>0.) PXIA(:) = (1-PALB_VEG(:))*(1.-PTR(:))*PIA(:)
!
IF (KNIV .EQ. 1) THEN
  DO I=1,SIZE(PIA)
    IF (PIA(I)>0.) THEN  
      ! -- reflection of surface ---   
      ! transmittance diffuse up - all layer
      ZTDF(I) = EXP(-ZB_DF*ZOMEGA_DF(I)*(1.-PABC)*PLAI(I))
      PXIA(I)= PXIA(I) + (1.-PALB_VEG(I))**2*PALB_SOIL(I)*(1.-ZTDF(I))*PTR(I)*PIA(I)
    ENDIF
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('CCETR_PAIR',1,ZHOOK_HANDLE)
!
CONTAINS
!
SUBROUTINE SET_PARAM(PSSA, PK, PP, PB, POMEGA_DR, PB_DR, POMEGA_DF, PB_DF)
!
REAL, INTENT(IN) :: PSSA, PK, PP
REAL, DIMENSION(:), INTENT(IN) :: PB
REAL, DIMENSION(:), INTENT(OUT) :: POMEGA_DR, PB_DR
REAL, DIMENSION(:), INTENT(OUT) :: POMEGA_DF
REAL, INTENT(OUT) :: PB_DF
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CCETR_PAIR:SET_PARAM',0,ZHOOK_HANDLE)
!
DO I=1,SIZE(PB_DR)
  IF (PIA(I).NE.0) THEN
    ! direct case
    ! Directional albedo of upper/lower layer
    PB_DR(I) = 1.-(1.-SQRT(1.-PSSA))/(1.+2.*PXMUS(I)*SQRT(1.-PSSA)) 
    ! CLUMPING INDEX 
    POMEGA_DR(I) = 1. / (1.+ PB(I)*EXP(-PK*(ACOS(PXMUS(I)))**PP))
    ! diffus case
    ! CLUMPING INDEX
    POMEGA_DF(I) = (1.+PB(I)/2.)/(1.+PB(I))
  ENDIF
ENDDO
!Bihemispherical albedo of upper/lower layer - diffuse case
PB_DF = 1.-(1.-SQRT(1.-PSSA))/(1.+ SQRT(1.-PSSA))
!
IF (LHOOK) CALL DR_HOOK('CCETR_PAIR:SET_PARAM',1,ZHOOK_HANDLE)
END SUBROUTINE SET_PARAM
!
END SUBROUTINE CCETR_PAIR
 
