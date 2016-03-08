!#############################################################
SUBROUTINE INIT_ISBA_LANDUSE (DTCO, UG, U, IO, IR, PMESH_SIZE, PDG, PDG_OLD, PLAI, &
                              PPATCH, PPATCH_OLD, HPROGRAM)  
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
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
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
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
REAL, DIMENSION(:), INTENT(IN) :: PMESH_SIZE
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDG
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDG_OLD
REAL, DIMENSION(:,:), INTENT(IN) :: PLAI
REAL, DIMENSION(:,:), INTENT(IN) :: PPATCH
REAL, DIMENSION(:,:), INTENT(IN) :: PPATCH_OLD
!
 CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(SIZE(PDG,1),SIZE(PDG,2),SIZE(PDG,3)) :: ZZDG     ! Actual layer thicknesses
REAL, DIMENSION(SIZE(PDG,1),SIZE(PDG,2),SIZE(PDG,3)) :: ZZDG_OLD ! Old layer thicknesses
REAL, DIMENSION(SIZE(PDG,1),SIZE(PDG,2),SIZE(PDG,3)) :: ZWG_OLD  ! Old XWG
REAL, DIMENSION(SIZE(PDG,1),SIZE(PDG,2),SIZE(PDG,3)) :: ZWGI_OLD ! Old XWGI
REAL, DIMENSION(SIZE(PDG,1),1,SIZE(PDG,3)) :: ZTEST
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
IF(ALL(PDG(:,IO%NGROUND_LAYER,:)==PDG_OLD(:,IO%NGROUND_LAYER,:)))THEN
  IF (LHOOK) CALL DR_HOOK('INIT_ISBA_LANDUSE',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
!-------------------------------------------------------------------------------
! Conserve mass in the cell
!-------------------------------------------------------------------------------
!
 CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'WR      ', IR%XWR     (:,:),0)

IF (IO%LGLACIER) CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'ICE_STO ', IR%XICE_STO(:,:),0)
!
DO JLAYER=1,SIZE(IR%XTG,2)
   CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'TEMP GRO', IR%XTG(:,JLAYER,:),0)
END DO
!
!
 CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'ALBSNOW ', IR%TSNOW%ALB(:,:),0)
!
IF (IR%TSNOW%SCHEME=='1-L'  .OR. IR%TSNOW%SCHEME=='3-L' .OR. IR%TSNOW%SCHEME=='CRO') THEN
   CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'EMISSNOW', IR%TSNOW%EMIS(:,:),0)    
   CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'TSSNOW  ', IR%TSNOW%TS  (:,:),0)
ENDIF
!
DO JLAYER=1,IR%TSNOW%NLAYER
   !
   CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'WSNOW   ', IR%TSNOW%WSNOW(:,JLAYER,:),0)
   !
   IF (IR%TSNOW%SCHEME=='3-L' .OR. IR%TSNOW%SCHEME=='CRO') THEN            
      CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'TEMPSNOW', IR%TSNOW%TEMP(:,JLAYER,:),0)
      CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'HEATSNOW', IR%TSNOW%HEAT(:,JLAYER,:),0)     
      CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'AGESNOW ', IR%TSNOW%AGE (:,JLAYER,:),0)
   ENDIF
   !
   IF (IR%TSNOW%SCHEME=='1-L') THEN
      CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'TSNOW   ', IR%TSNOW%T(:,JLAYER,:),0)
   ENDIF
   !
   IF(IR%TSNOW%SCHEME=='CRO') THEN
      CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'GRANSNOW', IR%TSNOW%GRAN1(:,JLAYER,:),0)
      CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'GRANSNOW', IR%TSNOW%GRAN2(:,JLAYER,:),0)
      CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'HISTSNOW', IR%TSNOW%HIST (:,JLAYER,:),0)
   ENDIF
   !
ENDDO
!
!-------------------------------------------------------------------------------
! Conserve mass globaly because soil depth change
!-------------------------------------------------------------------------------
!
ZWG_OLD(:,:,:) =IR%XWG (:,:,:)
ZWGI_OLD(:,:,:)=IR%XWGI(:,:,:)
!
DO JLAYER=1,IO%NGROUND_LAYER
   CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'WG      ', IR%XWG (:,JLAYER,:),0)
   CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'WGI     ', IR%XWGI(:,JLAYER,:),0)
ENDDO
!
ZZDG    (:,1,:)=PDG    (:,1,:)
ZZDG_OLD(:,1,:)=PDG_OLD(:,1,:)
IF(IO%CISBA=='DIF')THEN
  DO JLAYER=2,IO%NGROUND_LAYER
     ZZDG    (:,JLAYER,:)=PDG    (:,JLAYER,:)-PDG    (:,JLAYER-1,:)
     ZZDG_OLD(:,JLAYER,:)=PDG_OLD(:,JLAYER,:)-PDG_OLD(:,JLAYER-1,:)
  ENDDO
ELSE     
  ZZDG    (:,2,:)=PDG    (:,2,:)
  ZZDG_OLD(:,2,:)=PDG_OLD(:,2,:)
  IF(IO%CISBA=='3-L' )THEN
    ZZDG    (:,3,:)=PDG    (:,3,:)-PDG    (:,2,:)
    ZZDG_OLD(:,3,:)=PDG_OLD(:,3,:)-PDG_OLD(:,2,:)
  ENDIF 
ENDIF
!
WHERE(ZZDG(:,:,:)    >1.E+10)ZZDG    (:,:,:)=0.
WHERE(ZZDG_OLD(:,:,:)>1.E+10)ZZDG_OLD(:,:,:)=0.
!
 CALL CONSERV_GLOBAL_MASS(DTCO, U, PMESH_SIZE, PPATCH, PPATCH_OLD, &
                          ILUOUT,ZZDG,ZZDG_OLD,IR%XWG, ZWG_OLD )
 CALL CONSERV_GLOBAL_MASS(DTCO, U, PMESH_SIZE, PPATCH, PPATCH_OLD, &
                          ILUOUT,ZZDG,ZZDG_OLD,IR%XWGI,ZWGI_OLD)
!
!-------------------------------------------------------------------------------
! Extrapolation with 3 pts 
!-------------------------------------------------------------------------------
!
 CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'RESA    ', IR%XRESA(:,:),3)
!
DO JLAYER=1,IR%TSNOW%NLAYER
   CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'RHOSNOW ', IR%TSNOW%RHO  (:,JLAYER,:),3)
ENDDO
!
IF (IO%CPHOTO/='NON') THEN
   !
   CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'AN      ', IR%XAN   (:,:),3)
   CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'ANDAY   ', IR%XANDAY(:,:),3)   
   CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'ANFM    ', IR%XANFM (:,:),3)
   CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'LE      ', IR%XLE   (:,:),3)
   !
   DO JNBIOMASS=1,IO%NNBIOMASS
      CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'RESPBIOM', IR%XRESP_BIOMASS(:,JNBIOMASS,:),3)
      CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'BIOMASS ', IR%XBIOMASS     (:,JNBIOMASS,:),3)
   ENDDO
   !
   IF (IO%CRESPSL=='CNT') THEN
      !
      DO JNLITTLEVS=1,IO%NNLITTLEVS
         CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'LIGNINST',IR%XLIGNIN_STRUC(:,JNLITTLEVS,:),3)
         DO JNLITTER=1,IO%NNLITTER
            CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'LITTER  ',IR%XLITTER(:,JNLITTER,JNLITTLEVS,:),3)
         ENDDO
      ENDDO
      !
      DO JNSOILCARB=1,IO%NNSOILCARB
         CALL INI_VAR_FROM_PATCH(DTCO, UG, U, PPATCH, PPATCH_OLD, PLAI, &
                         HPROGRAM,ILUOUT,'SOILCARB',IR%XSOILCARB(:,JNSOILCARB,:),3)
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
