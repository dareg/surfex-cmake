!#############################################################
SUBROUTINE INIT_ISBA_LANDUSE (DTCO, IG, I, UG, U, &
                               HPROGRAM)  
!#############################################################
!
!!****  *INIT_ISBA_LANDUSE* - routine to initialize land use for ISBA field
!!
!!    PURPOSE
!!    -------
!     Extrapolation from existing surounding cells with same patch properties:
!!      (1) IPTS=n  interpol field with n pts
!!      (2) IPTS=0  conserve cells mass  
!!   Case 2 : simple extrapolation based on the inside cell informations.
!!             this is donne before conserving cell or global mass
!!
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
!!      Original    07/2011
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_TYPE_SNOW
USE MODD_SURF_PAR,ONLY : XUNDEF                 
!
USE MODI_GET_LUOUT
USE MODI_INI_VAR_FROM_PATCH
USE MODI_CONSERV_GLOBAL_MASS
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
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(SIZE(I%M%X%XDG,1),SIZE(I%M%X%XDG,2),SIZE(I%M%X%XDG,3)) :: ZZDG     ! Actual layer thicknesses
REAL, DIMENSION(SIZE(I%M%X%XDG,1),SIZE(I%M%X%XDG,2),SIZE(I%M%X%XDG,3)) :: ZZDG_OLD ! Old layer thicknesses
REAL, DIMENSION(SIZE(I%M%X%XDG,1),SIZE(I%M%X%XDG,2),SIZE(I%M%X%XDG,3)) :: ZWG_OLD  ! Old XWG
REAL, DIMENSION(SIZE(I%M%X%XDG,1),SIZE(I%M%X%XDG,2),SIZE(I%M%X%XDG,3)) :: ZWGI_OLD ! Old XWGI
REAL, DIMENSION(SIZE(I%M%X%XDG,1),1,SIZE(I%M%X%XDG,3)) :: ZTEST
!
INTEGER :: ILUOUT
INTEGER :: JLAYER, JNBIOMASS, JNLITTER, JNLITTLEVS, JNSOILCARB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_ISBA_LANDUSE',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------
!
IF(ALL(I%M%X%XDG(:,I%O%NGROUND_LAYER,:)==I%M%X%XDG_OLD(:,I%O%NGROUND_LAYER,:)))THEN
  IF (LHOOK) CALL DR_HOOK('INIT_ISBA_LANDUSE',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
!-------------------------------------------------------------------------------
! Conserve mass in the cell
!-------------------------------------------------------------------------------
!
 CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'WR      ', I%R%XWR     (:,:),0)

IF (I%O%LGLACIER) CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'ICE_STO ', I%R%XICE_STO(:,:),0)
!
DO JLAYER=1,SIZE(I%R%XTG,2)
   CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'TEMP GRO', I%R%XTG(:,JLAYER,:),0)
END DO
!
!
 CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'ALBSNOW ', I%R%TSNOW%ALB(:,:),0)
!
IF (I%R%TSNOW%SCHEME=='1-L'  .OR. I%R%TSNOW%SCHEME=='3-L' .OR. I%R%TSNOW%SCHEME=='CRO') THEN
   CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'EMISSNOW', I%R%TSNOW%EMIS(:,:),0)    
   CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'TSSNOW  ', I%R%TSNOW%TS  (:,:),0)
ENDIF
!
DO JLAYER=1,I%R%TSNOW%NLAYER
   !
   CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'WSNOW   ', I%R%TSNOW%WSNOW(:,JLAYER,:),0)
   !
   IF (I%R%TSNOW%SCHEME=='3-L' .OR. I%R%TSNOW%SCHEME=='CRO') THEN            
      CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'TEMPSNOW', I%R%TSNOW%TEMP(:,JLAYER,:),0)
      CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'HEATSNOW', I%R%TSNOW%HEAT(:,JLAYER,:),0)     
      CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'AGESNOW ', I%R%TSNOW%AGE (:,JLAYER,:),0)
   ENDIF
   !
   IF (I%R%TSNOW%SCHEME=='1-L') THEN
      CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'TSNOW   ', I%R%TSNOW%T(:,JLAYER,:),0)
   ENDIF
   !
   IF(I%R%TSNOW%SCHEME=='CRO') THEN
      CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'GRANSNOW', I%R%TSNOW%GRAN1(:,JLAYER,:),0)
      CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'GRANSNOW', I%R%TSNOW%GRAN2(:,JLAYER,:),0)
      CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'HISTSNOW', I%R%TSNOW%HIST (:,JLAYER,:),0)
   ENDIF
   !
ENDDO
!
!-------------------------------------------------------------------------------
! Conserve mass globaly because soil depth change
!-------------------------------------------------------------------------------
!
ZWG_OLD(:,:,:) =I%R%XWG (:,:,:)
ZWGI_OLD(:,:,:)=I%R%XWGI(:,:,:)
!
DO JLAYER=1,I%O%NGROUND_LAYER
   CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'WG      ', I%R%XWG (:,JLAYER,:),0)
   CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'WGI     ', I%R%XWGI(:,JLAYER,:),0)
ENDDO
!
ZZDG    (:,1,:)=I%M%X%XDG    (:,1,:)
ZZDG_OLD(:,1,:)=I%M%X%XDG_OLD(:,1,:)
IF(I%O%CISBA=='DIF')THEN
  DO JLAYER=2,I%O%NGROUND_LAYER
     ZZDG    (:,JLAYER,:)=I%M%X%XDG    (:,JLAYER,:)-I%M%X%XDG    (:,JLAYER-1,:)
     ZZDG_OLD(:,JLAYER,:)=I%M%X%XDG_OLD(:,JLAYER,:)-I%M%X%XDG_OLD(:,JLAYER-1,:)
  ENDDO
ELSE     
  ZZDG    (:,2,:)=I%M%X%XDG    (:,2,:)
  ZZDG_OLD(:,2,:)=I%M%X%XDG_OLD(:,2,:)
  IF(I%O%CISBA=='3-L' )THEN
    ZZDG    (:,3,:)=I%M%X%XDG    (:,3,:)-I%M%X%XDG    (:,2,:)
    ZZDG_OLD(:,3,:)=I%M%X%XDG_OLD(:,3,:)-I%M%X%XDG_OLD(:,2,:)
  ENDIF 
ENDIF
!
WHERE(ZZDG(:,:,:)    >1.E+10)ZZDG    (:,:,:)=0.
WHERE(ZZDG_OLD(:,:,:)>1.E+10)ZZDG_OLD(:,:,:)=0.
!
 CALL CONSERV_GLOBAL_MASS(DTCO, IG, I, U, &
                          ILUOUT,ZZDG,ZZDG_OLD,I%R%XWG, ZWG_OLD )
 CALL CONSERV_GLOBAL_MASS(DTCO, IG, I, U, &
                          ILUOUT,ZZDG,ZZDG_OLD,I%R%XWGI,ZWGI_OLD)
!
!-------------------------------------------------------------------------------
! Extrapolation with 3 pts 
!-------------------------------------------------------------------------------
!
 CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'RESA    ', I%R%XRESA(:,:),3)
!
DO JLAYER=1,I%R%TSNOW%NLAYER
   CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'RHOSNOW ', I%R%TSNOW%RHO  (:,JLAYER,:),3)
ENDDO
!
IF (I%O%CPHOTO/='NON') THEN
   !
   CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'AN      ', I%R%XAN   (:,:),3)
   CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'ANDAY   ', I%R%XANDAY(:,:),3)   
   CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'ANFM    ', I%R%XANFM (:,:),3)
   CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'LE      ', I%R%XLE   (:,:),3)
   !
   DO JNBIOMASS=1,I%O%NNBIOMASS
      CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'RESPBIOM', I%R%XRESP_BIOMASS(:,JNBIOMASS,:),3)
      CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'BIOMASS ', I%R%XBIOMASS     (:,JNBIOMASS,:),3)
   ENDDO
   !
   IF (I%O%CRESPSL=='CNT') THEN
      !
      DO JNLITTLEVS=1,I%O%NNLITTLEVS
         CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'LIGNINST',I%R%XLIGNIN_STRUC(:,JNLITTLEVS,:),3)
         DO JNLITTER=1,I%O%NNLITTER
            CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'LITTER  ',I%R%XLITTER(:,JNLITTER,JNLITTLEVS,:),3)
         ENDDO
      ENDDO
      !
      DO JNSOILCARB=1,I%O%NNSOILCARB
         CALL INI_VAR_FROM_PATCH(DTCO, I, UG, U, &
                         HPROGRAM,ILUOUT,'SOILCARB',I%R%XSOILCARB(:,JNSOILCARB,:),3)
      ENDDO
      !
   ENDIF
   !
ENDIF
!
!-------------------------------------------------------------------------------
!  
IF (LHOOK) CALL DR_HOOK('INIT_ISBA_LANDUSE',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_ISBA_LANDUSE
