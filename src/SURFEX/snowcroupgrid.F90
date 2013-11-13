!     #########
SUBROUTINE SNOWCROUPGRID(PSNOW,PSNOWDZ,PSNOWDZN,PSNOWRHO,          &
                         PSNOWHEAT,PSNOWGRAN1, PSNOWGRAN2,         &
                         PSNOWHIST,PSNOWHMASS,OSNOWFALL, PSNOWBIS, &
                         PSNOWDZBIS,                               &
                         PSNOWHEATBIS,PSNOWRHOBIS, PSNOWGRAN1BIS,  &
                         PSNOWGRAN2BIS,PSNOWHISTBIS,PTSTEP,PSR,    &
                         PTA,PVMOD                                 )
!
USE MODD_SNOW_PAR, ONLY : XSNOWCRITD
USE MODE_SNOW3L
!
! modifs_EB: transformation de grille uniquement si chute de neige, une couche
! trop petite ou HTN < 3 cm (==>omodifgrid=.T.)
! modifs pour traiter les grains comme les variables d'origine
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!       0.1 declarations of arguments        
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ,PSNOWRHO,PSNOWDZN,PSNOWHEAT   
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZBIS,PSNOWHEATBIS,PSNOWRHOBIS, &
                                       PSNOWGRAN1BIS,PSNOWGRAN2BIS,PSNOWHISTBIS                     
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST    
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOW,PSNOWBIS,PSNOWHMASS 
!
REAL, DIMENSION(:),      INTENT(IN) :: PSR, PTA, PVMOD
!
LOGICAL,DIMENSION(:),    INTENT(IN) :: OSNOWFALL 
!
REAL, INTENT(IN)  :: PTSTEP 
!
!       0.2 declaration of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWHEATN,ZSNOWRHON  
! ajout EB
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWGRAN1N, &
                                                      ZSNOWGRAN2N,ZSNOWHISTN
REAL, DIMENSION(SIZE(PSNOW))   :: ZSUMHEAT, ZSUMSWE,  ZSNOWMIX_DELTA, ZSNOWHMASS         
REAL, DIMENSION(SIZE(PSNOW))   :: ZNDENT, ZNVIEU
REAL, DIMENSION(SIZE(PSNOW))   :: ZSNOWMIN
!
REAL :: ZTSTEP
!
INTEGER :: INLVLS
INTEGER :: JJ,JST,JIT
!
LOGICAL :: GMODIFGRID
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!       0.3 initialization
!
IF (LHOOK) CALL DR_HOOK('SNOWCROUPGRID',0,ZHOOK_HANDLE)
!
INLVLS = SIZE(PSNOWDZ(:,:),2)
!      
!       1. update grid for snowpack>3cm
!
!
! ajout EB
IF( PSNOWGRAN1(1,3)>80.0 .AND. PSNOWGRAN2(1,3)>2.E-2 ) THEN
  WRITE(*,*) 'pb1 psnowgran',PSNOWGRAN1(1,3),PSNOWGRAN2(1,3),&
                             PSNOWRHO(1,3),PSNOWDZ(1,3)  
ENDIF
!
IF( PSNOWGRAN1BIS(1,3)>80.0 .AND. PSNOWGRAN2BIS(1,3)>2.E-2 ) THEN
  WRITE(*,*) 'pb1 psnowgranBIS',PSNOWGRAN1BIS(1,3),PSNOWGRAN2BIS(1,3)
ENDIF
!
GMODIFGRID = .FALSE.
!         
DO JJ = 1,SIZE(PSNOW(:))
  !
  IF( OSNOWFALL(JJ) .AND. GMODIFGRID ) THEN
! IF(OSNOWFALL(JJ)) THEN
    !
    WRITE(*,*) 'ALERTE'
    ! ajout EB
    GMODIFGRID = .TRUE.
    !
    ZTSTEP = PTSTEP/10.0
    PSNOWHMASS(JJ) = 0.
    !
    DO JIT = 1,10
      !
      PSNOWHMASS(JJ) = 0.
      !
      CALL SNOWCROADDSNOW(ZTSTEP,PSR(JJ),PTA(JJ),PVMOD(JJ),                       &
                          PSNOW(JJ),PSNOWRHO(JJ,1),PSNOWDZ(JJ,1),PSNOWHEAT(JJ,1), &
                          ZSNOWHMASS(JJ),PSNOWGRAN1(JJ,1),PSNOWGRAN2(JJ,1),       &
                          PSNOWHIST(JJ,1),INLVLS                                  )
      !
      PSNOWHMASS (JJ)   = PSNOWHMASS(JJ) + ZSNOWHMASS(JJ)
      !   
      ZSNOWRHON  (JJ,:) = PSNOWRHO  (JJ,:)
      ZSNOWHEATN (JJ,:) = PSNOWHEAT (JJ,:)      
      !ajout EB           
      ZSNOWGRAN1N(JJ,:) = PSNOWGRAN1(JJ,:)      
      ZSNOWGRAN2N(JJ,:) = PSNOWGRAN2(JJ,:)      
      ZSNOWHISTN (JJ,:) = PSNOWHIST (JJ,:)      
      !
      CALL SNOW3LGRID(PSNOWDZN(JJ,:),PSNOW(JJ))
      !
      CALL SNOWCROTRANSF_1D(PSNOW(JJ),PSNOWDZ(JJ,:),PSNOWDZN(JJ,:),             &
                            ZSNOWRHON(JJ,:),ZSNOWHEATN(JJ,:),ZSNOWGRAN1N(JJ,:), &
                            ZSNOWGRAN2N(JJ,:),ZSNOWHISTN(JJ,:)                  )
      !
    END DO
    !
  ELSE 
    !
    ZSNOWRHON  (JJ,:) = PSNOWRHOBIS  (JJ,:)
    ZSNOWHEATN (JJ,:) = PSNOWHEATBIS (JJ,:)
    !ajout EB           
    ZSNOWGRAN1N(JJ,:) = PSNOWGRAN1BIS(JJ,:)
    ZSNOWGRAN2N(JJ,:) = PSNOWGRAN2BIS(JJ,:)
    ZSNOWHISTN (JJ,:) = PSNOWHISTBIS (JJ,:)
    !
    ! ajout EB         
    ! on change de grille seulement si il y a une trop petite couche
    ZSNOWMIN(JJ) = PSNOW(JJ)
    !
    DO JST = 1,INLVLS
      IF ( PSNOWDZ(JJ,JST)<ZSNOWMIN(JJ) ) ZSNOWMIN(JJ) = PSNOWDZ(JJ,JST)
    ENDDO
    !
    IF ( ZSNOWMIN(JJ)<MIN( 0.001, PSNOW(JJ)/(2*INLVLS) ) .OR. PSR(JJ)>0.0 ) THEN
      !  
      GMODIFGRID=.TRUE.
      !
      ! WRITE (*,*) 'avant recalcul: snowmin=',psnowmin(jj1),'psr=', psr(jj1) 
      ! WRITE (*,*) PSNOW(JJ), PSNOWBIS(JJ)
      ! WRITE (*,*) PSNOWRHO(JJ,1), PSNOWRHOBIS(JJ,1)
      ! WRITE (*,*) PSNOWHEAT(JJ,1), PSNOWHEATBIS(JJ,1)
      ! WRITE (*,*) PSNOWDZ(JJ,1), PSNOWDZBIS(JJ,1)
      CALL SNOW3LGRID(PSNOWDZN(JJ,:),PSNOWBIS(JJ))
      !
      WRITE(*,*) 'psr, psnowmin,psnow', PSR(JJ), ZSNOWMIN(JJ),PSNOW(JJ)
      !
      !  CALL SNOW3LTRANSF_1D(PSNOWBIS(JJ),PSNOWDZBIS(JJ,:),PSNOWDZN(JJ,:),  &
      !                 ZSNOWRHON(JJ,:),ZSNOWHEATN(JJ,:),ZSNOWGRAN1N(JJ,:),  &
      !                 ZSNOWGRAN2N(JJ,:),ZSNOWHISTN(JJ,:),OSNOW_METAMO)
      CALL SNOWNLTRANSFGRID_1D(PSNOWBIS(JJ),PSNOWDZBIS(JJ,:),PSNOWDZN(JJ,:),       &
                               ZSNOWRHON(JJ,:),ZSNOWHEATN(JJ,:),ZSNOWGRAN1N(JJ,:), &
                               ZSNOWGRAN2N(JJ,:),ZSNOWHISTN(JJ,:)                  )  
      !  WRITE(*,*) 'après'
      !  WRITE (*,*) PSNOW(JJ), PSNOWBIS(JJ)
      !  WRITE (*,*) PSNOWRHO(JJ,1), PSNOWRHOBIS(JJ,1)
      !  WRITE (*,*) PSNOWHEAT(JJ,1), PSNOWHEATBIS(JJ,1)
      !  WRITE (*,*) PSNOWDZ(JJ,1), PSNOWDZBIS(JJ,1)
      !
    ENDIF
    !
    PSNOW     (JJ)   = PSNOWBIS     (JJ)  
    PSNOWRHO  (JJ,:) = PSNOWRHOBIS  (JJ,:)
    PSNOWHEAT (JJ,:) = PSNOWHEATBIS (JJ,:) 
    PSNOWDZ   (JJ,:) = PSNOWDZBIS   (JJ,:)                       
    !          
    PSNOWGRAN1(JJ,:) = PSNOWGRAN1BIS(JJ,:)
    PSNOWGRAN2(JJ,:) = PSNOWGRAN2BIS(JJ,:) 
    PSNOWHIST (JJ,:) = PSNOWHISTBIS (JJ,:)
    !
  ENDIF  
  !
  IF ( PSNOWGRAN1(1,3)>80.0 .AND. PSNOWGRAN2(1,3)>2.E-2 ) THEN
    WRITE(*,*) 'pb2 psnowgran',PSNOWGRAN1(1,3),PSNOWGRAN2(1,3)
  ENDIF
  !
  IF( PSNOWGRAN1BIS(1,3)>80.0 .AND. PSNOWGRAN2BIS(1,3)>2.E-2 ) THEN
    WRITE(*,*) 'pb2 psnowgranBIS',PSNOWGRAN1BIS(1,3),PSNOWGRAN2BIS(1,3)
  ENDIF
  !
  IF( ZSNOWGRAN1N(1,3)>80.0 .AND.  ZSNOWGRAN2N(1,3)>2.E-2 ) THEN
    WRITE(*,*) 'pb2 psnowgranN',ZSNOWGRAN1N(1,3),ZSNOWGRAN2N(1,3), &
                                ZSNOWHISTN(1,3),PSNOW(1),INLVLS
  ENDIF
  !
END DO       
!
!      2. update grid for thin snowpack<3 cm 
!
ZSUMHEAT      (:) = 0.0
ZSUMSWE       (:) = 0.0
ZSNOWMIX_DELTA(:) = 0.0
ZNDENT        (:) = 0.0
ZNVIEU        (:) = 0.0      
!        
DO JJ = 1,SIZE(PSNOWHEAT,1)
  !
  IF( PSNOW(JJ)<XSNOWCRITD ) THEN
    !
    ZSNOWMIX_DELTA(JJ) = 1.0  
    !
    ! ajout EB                
    GMODIFGRID = .TRUE.
    !
    DO JST = 1,INLVLS
      !
      ZSUMHEAT(JJ) = ZSUMHEAT(JJ) + ZSNOWHEATN(JJ,JST)
      ZSUMSWE (JJ) = ZSUMSWE (JJ) + ZSNOWRHON (JJ,JST) * PSNOWDZ(JJ,JST)
      !           
      IF ( PSNOWGRAN1(JJ,JST)<0. ) THEN       ! Dendritic snow
        ZNDENT(JJ) = ZNDENT(JJ) + 1.0
      ELSE                             ! Non dendritic snow
        ZNVIEU(JJ) = ZNVIEU(JJ) + 1.0
      ENDIF
      !
    END DO
    !
  ENDIF
  !
END DO
!
! Average properties for grains : determine which grain type is the most
! important in the snowpack.
!IF(OSNOW_METAMO)THEN
! modifs EB pour changer variables de la subroutine suivante        
! IF((PSNOWGRAN1(1,3)>80.0).and.(PSNOWGRAN2(1,3)>2.E-2))then
!       WRITE(*,*) 'pb2b psnowgran',PSNOWGRAN1(1,3),PSNOWGRAN2(1,3) &
!              , PSNOWHIST(1,3),PSNOW(1),zndent(1),znvieu(1),inlvls
!endif
! ajout EB suppression temporaire de cet appel pour vérif TRANSFGRID 
 CALL SNOW3LAVGRAIN(PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,                 &
                    ZSNOWGRAN1N,ZSNOWGRAN2N,ZSNOWHISTN,ZNDENT,ZNVIEU )   
!   IF((ZSNOWGRAN1N(1,3)>80.0).and.(ZSNOWGRAN2N(1,3)>2.E-2))then
!         WRITE(*,*) 'pb2b psnowgranN',ZSNOWGRAN1N(1,3),ZSNOWGRAN2N(1,3) &
!                , ZSNOWHISTN(1,3),PSNOW(1),zndent(1),znvieu(1),inlvls
! endif
!ENDIF
!
!IF(ZSNOWMIX_DELTA(1)>0.)THEN
 !       WRITE(*,*) 'PG1_1',PSNOWGRAN1(1,1),'PG2_1',PSNOWGRAN2(1,1)
!ENDIF
!
DO JST = 1,INLVLS
  !
  ZSNOWHEATN(:,JST) = ZSNOWMIX_DELTA(:)        * ZSUMHEAT(:)/INLVLS + &
                     ( 1.0-ZSNOWMIX_DELTA(:) ) * ZSNOWHEATN(:,JST)  
  !
  PSNOWDZN  (:,JST) = ZSNOWMIX_DELTA(:)         * PSNOW(:)/INLVLS + &
                      ( 1.0-ZSNOWMIX_DELTA(:) ) * PSNOWDZN(:,JST)  
  !
  ZSNOWRHON (:,JST) = ZSNOWMIX_DELTA(:)         * ZSUMSWE(:)/PSNOW(:) + &
                      ( 1.0-ZSNOWMIX_DELTA(:) ) * ZSNOWRHON(:,JST)  
  !
ENDDO
!                      
!     3. Update mass (density and thickness) and heat:
!
! ajout EB pour ne faire cet update que pour les couches fines 
IF( ZSNOWGRAN1N(1,3)>80.0 .AND. ZSNOWGRAN2N(1,3)>2.E-2 ) THEN
  WRITE(*,*) 'pb3 psnowgran',ZSNOWGRAN1N(1,3),ZSNOWGRAN2N(1,3)
ENDIF
!
IF( PSNOWGRAN1BIS(1,3)>80.0 .AND. PSNOWGRAN2BIS(1,3)>2.E-2 ) THEN
  WRITE(*,*) 'pb3 psnowgranBIS',PSNOWGRAN1BIS(1,3),PSNOWGRAN2BIS(1,3)
ENDIF
!
IF ( GMODIFGRID ) THEN
  ! WRITE(*,*) omodifgrid, psnowmin(1), psr(1), zsnowmix_delta(1)        
  PSNOWRHO  (:,:) = ZSNOWRHON  (:,:)
  PSNOWDZ   (:,:) = PSNOWDZN   (:,:)
  PSNOWHEAT (:,:) = ZSNOWHEATN (:,:) 
  !
  PSNOWGRAN1(:,:) = ZSNOWGRAN1N(:,:)
  PSNOWGRAN2(:,:) = ZSNOWGRAN2N(:,:) 
  PSNOWHIST (:,:) = ZSNOWHISTN (:,:)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SNOWCROUPGRID',1,ZHOOK_HANDLE)
! 
!#############################################################
!#############################################################
!#############################################################
!#############################################################
!
CONTAINS
!
SUBROUTINE SNOWCROADDSNOW(PTSTEP,PSR,PTA,PVMOD,                        &
                          PSNOW,PSNOWRHO,PSNOWDZ,PSNOWHEAT,PSNOWHMASS, &
                          PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,KNLVLS       )  
!!    PURPOSE
!!    -------
!     Add snow to snowpack 
!     Update mass and heat content of uppermost layer.
!
USE MODD_CSTS,     ONLY : XLMTT, XTT, XCI
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES, XSNOWDMIN, XANSMAX, XSNOWCRITD, &
                          XSNOWFALL_A_SN, XSNOWFALL_B_SN, XSNOWFALL_C_SN
!
USE MODD_SNOW_METAMO, ONLY : XNDEN1, XNDEN2, XNDEN3, XGRAN,  &
                             XNSPH1, XNSPH2, XNSPH3, XNSPH4, &
                             XDIAET, XDIAFP, XDIAGF
!
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)    :: PTSTEP
!
REAL, INTENT(IN)    :: PSR, PTA, PVMOD
!
REAL, INTENT(INOUT) :: PSNOW
!
REAL, INTENT(INOUT) :: PSNOWRHO, PSNOWDZ, PSNOWHEAT
!
REAL, INTENT(OUT)   :: PSNOWHMASS
!
REAL, INTENT(INOUT) :: PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST
!   
INTEGER, INTENT(IN) :: KNLVLS  
!
!*      0.2    declarations of local variables
!
REAL :: ZSNOWFALL, ZRHOSNEW, ZSNOW, ZSNOWTEMP, ZSNOWFALL_DELTA, &
        ZSCAP, ZSDEN, ZSPHE, ZDIAMD, ZDIAMV, ZDIAMN, ZSPHERD,   &
        ZSPHERV, ZSPHERN, ZDENT, ZSPHERN2, ZCOEF1, ZCOEF2   
!
INTEGER :: JJ
!
! ISBA-ES CROCUS (Pahaut 1976): snowfall density coefficients:
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCROADDSNOW',0,ZHOOK_HANDLE)
!
ZRHOSNEW   = XRHOSMIN_ES
ZSNOWFALL  = 0.0
ZSCAP      = 0.0
ZSNOW      = PSNOW
ZSDEN      = 0.0
ZSPHE      = 0.0
!
PSNOWHMASS = 0.0
!
! 1. Add snow into snowpack:
! --------------------------
!
! Heat content of newly fallen snow (J/m2):
! NOTE for now we assume the snowfall has
! the temperature of the snow surface upon reaching the snow.
! This is done as opposed to using the air temperature since
! this flux is quite small and has little to no impact
! on the time scales of interest. If we use the above assumption
! then, then the snowfall advective heat flux is zero.
!
ZSCAP = PSNOWRHO * XCI
!
ZSNOWTEMP = XTT + ( PSNOWHEAT +  XLMTT*PSNOWRHO*PSNOWDZ ) / &
                  ( ZSCAP * MAX( XSNOWDMIN/KNLVLS, PSNOWDZ ) )  
ZSNOWTEMP = MIN( XTT, ZSNOWTEMP )
!
PSNOWHMASS = PSR * ( XCI * (ZSNOWTEMP-XTT) - XLMTT ) * PTSTEP
!
! Snowfall density: Following CROCUS (Pahaut 1976)
!
ZRHOSNEW = MAX( XRHOSMIN_ES, XSNOWFALL_A_SN +               &
                             XSNOWFALL_B_SN * ( PTA-XTT ) + &
                             XSNOWFALL_C_SN * SQRT(PVMOD)   )   
!
! Augment total pack depth:
!
ZSNOWFALL  = PSR * PTSTEP / ZRHOSNEW    ! snowfall thickness (m)
!
! Fresh snowfall changes the snowpack
! density, increases the total liquid water
! equivalent: in uppermost snow layer:
!
!
ZSDEN = MAX( MIN( XNDEN1*PVMOD-XNDEN2, XNDEN3 ), -XGRAN )
ZSPHE = MIN( MAX( XNSPH1*PVMOD+XNSPH2, XNSPH3 ), XNSPH4 )
!
ZCOEF1 = PSNOWDZ * PSNOWRHO
ZCOEF2 = ZSNOWFALL * ZRHOSNEW
!
IF ( PSNOWGRAN1<0.0 ) THEN
  !
  PSNOWGRAN1 = ( ZCOEF1 * PSNOWGRAN1 + ZCOEF2 * ZSDEN ) / ( ZCOEF1 + ZCOEF2 )  
  PSNOWGRAN2 = ( ZCOEF1 * PSNOWGRAN2 + ZCOEF2 * ZSPHE ) / ( ZCOEF1 + ZCOEF2 )
  !
ELSE
  !
  ZSPHERD =  ZSPHE / XGRAN
  ZDIAMD  = -ZSDEN / XGRAN        * XDIAET + &
            ( 1 + ZSDEN / XGRAN ) * ( ZSPHERD * XDIAGF + ( 1 - ZSPHERD ) * XDIAFP )    
  ! 
  ZSPHERV = PSNOWGRAN1 / XGRAN  
  ZDIAMV  = PSNOWGRAN2
  !
  ZDIAMN  = ( ZCOEF1 * ZDIAMV   + ZCOEF2 * ZDIAMD  ) / ( ZCOEF1 + ZCOEF2 )  
  ZSPHERN = ( ZCOEF1 * ZSPHERV  + ZCOEF2 * ZSPHERD ) / ( ZCOEF1 + ZCOEF2 )  
  !
  ZSPHERN2 = ZSPHERN * XDIAGF + (1-ZSPHERN) * XDIAFP
  IF( ZDIAMN < ZSPHERN2 ) THEN
    !
    ZDENT = ( ZDIAMN - ZSPHERN2 ) / ( XDIAET - ZSPHERN2 )
    PSNOWGRAN1 = -XGRAN * ZDENT
    PSNOWGRAN2 =  XGRAN * ZSPHERN  
    !
  ELSE
    !
    PSNOWGRAN1 = XGRAN * ZSPHERN
    PSNOWGRAN2 = ZDIAMN
    !
  ENDIF 
  !
ENDIF
!
PSNOW = PSNOW + ZSNOWFALL
!
PSNOWRHO = ( PSNOWDZ * PSNOWRHO + ZSNOWFALL * ZRHOSNEW ) / &
           ( PSNOWDZ            + ZSNOWFALL )  
!
PSNOWDZ  = PSNOWDZ + ZSNOWFALL
!
! Add energy of snowfall to snowpack:
! Update heat content (J/m2) (therefore the snow temperature
! and liquid content):
!
PSNOWHEAT  = PSNOWHEAT + PSNOWHMASS
!
PSNOWHIST = 0.0    
!
IF (LHOOK) CALL DR_HOOK('SNOWCROADDSNOW',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROADDSNOW
!
!################################################################################
!################################################################################
!################################################################################
!
SUBROUTINE SNOWCROTRANSF_1D(PSNOW,PSNOWDZ,PSNOWDZN,        &
                            PSNOWRHO,PSNOWHEAT,PSNOWGRAN1, &
                            PSNOWGRAN2, PSNOWHIST          )  
!
!!    PURPOSE
!!    -------
!     Snow mass and heat redistibution due to grid thickness
!     configuration resetting. Total mass and heat content
!     of the overall snowpack unchanged/conserved within this routine.
!
USE MODD_SNOW_PAR, ONLY : XSNOWCRITD
USE MODD_SNOW_METAMO
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                  :: PSNOW
!
REAL, DIMENSION(:), INTENT(INOUT) :: PSNOWHEAT, PSNOWRHO, PSNOWDZ,     &
                                     PSNOWDZN, PSNOWGRAN1, PSNOWGRAN2, &
                                     PSNOWHIST    
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1)) :: ZSNOWRHON,ZSNOWGRAN1,ZSNOWGRAN2,  &
                                     ZSNOWHEATN,ZSNOWHIST                                         
!
REAL, DIMENSION(SIZE(PSNOWRHO,1)-1) :: ZSNOWZO, ZSNOWZN, ZSNOWDDZ, ZDELTA, &
                                       ZRHO, ZHEAT
!
REAL :: ZHTN_NEW, ZHTN_OLD
!
INTEGER :: JJ, JLAYER1, JLAYER2
INTEGER :: INLVLS
!
! ajout EB
!integer :: jjj
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCROTRANSF_1D',0,ZHOOK_HANDLE)
!
INLVLS = SIZE(PSNOWDZ(:),1)

ZSNOWGRAN1(:) = PSNOWGRAN1(:)
ZSNOWGRAN2(:) = PSNOWGRAN2(:)
ZSNOWHIST (:) = PSNOWHIST (:)
!ajout EB pour vérifier égalité des HTN

ZHTN_NEW = 0.
ZHTN_OLD = 0.
!
DO JJ = 1,INLVLS
  ZHTN_NEW = ZHTN_NEW + PSNOWDZ (JJ)
  ZHTN_OLD = ZHTN_OLD + PSNOWDZN(JJ)
ENDDO
!
IF ( ABS( ZHTN_NEW-ZHTN_OLD ) > 1.E-7 ) WRITE(*,*) 'DIFHTN:', ZHTN_NEW, ZHTN_OLD
!
!
! 1. Calculate vertical grid depths (m):
! --------------------------------------
!
ZSNOWZO(1) = PSNOWDZ (1)
ZSNOWZN(1) = PSNOWDZN(1)
!
DO JJ = 2,INLVLS-1
  ZSNOWZO(JJ) = ZSNOWZO(JJ-1) + PSNOWDZ (JJ)
  ZSNOWZN(JJ) = ZSNOWZN(JJ-1) + PSNOWDZN(JJ)
ENDDO
!
! 2. Calculate thickness changes (m):
! -----------------------------------
!
ZSNOWDDZ(:) = ZSNOWZN(:) - ZSNOWZO(:)
!
! Calculate the delta function:
!
ZDELTA(:) = 0.0
WHERE( ZSNOWDDZ(:)>0.0 ) ZDELTA(:) = 1.0
!
! 3. Calculate mass and heat transfers due to grid adjustment/changes:
! --------------------------------------------------------------------
! Do transfers at boundaries first:
! Upper layer:
!
DO JJ = 2,INLVLS
  ZRHO (JJ-1) = ZDELTA(JJ-1) * PSNOWRHO(JJ) + (1.0-ZDELTA(JJ-1)) * PSNOWRHO(JJ-1)
  ZHEAT(JJ-1) = ZDELTA(JJ-1)       * PSNOWHEAT(JJ)  /PSNOWDZ(JJ) + &
                (1.0-ZDELTA(JJ-1)) * PSNOWHEAT(JJ-1)/PSNOWDZ(JJ-1)
ENDDO
!
!WRITE(*,*) 'AV_AGREG1', PSNOWGRAN1(1)
!WRITE(*,*) 'AP_AGREG1', PSNOWGRAN1(1)                       
! 
!IF(ZSNOWGRAN1(1)>99.OR.ZSNOWGRAN1(1)<-99.)THEN
!     WRITE(*,*) 'ZG1',ZSNOWGRAN1(1)
!     read(*,*)
!ENDIF           
! Lowest layer:
!
! Update interior layer mass and heat :
!
!WRITE(*,*) 'AV_AGREG8', PSNOWGRAN1(8),'ZD7',ZDELTA(7),'ZD9',ZDELTA(9)
DO JJ = 1,INLVLS
  !
  IF ( JJ>1 ) THEN
    IF ( ZDELTA(JJ-1)==0.0 ) THEN
      JLAYER1 = JJ
      JLAYER2 = JJ - 1
      CALL SNOW3LAGREG(PSNOWDZN,PSNOWDZ,PSNOWRHO,PSNOWGRAN1,PSNOWGRAN2, &
                       PSNOWHIST,ZSNOWGRAN1,ZSNOWGRAN2,ZSNOWHIST,       &
                       JLAYER1,JLAYER2,ZSNOWDDZ                         )
!           if (zsnowgran1(2) > 900.or. psnowgran1(2)>900.) then
!                   WRITE (*,*) 'agr',PSNOWGRAN1(2),PSNOWGRAN1(2), &
!                   zSNOWGRAN1(2),zSNOWGRAN1(2), JLAYER1,JLAYER2
!                   WRITE(*,*) 'dzold:',(psnowdz(jjj),jjj=1, INLVLS)
!                   WRITE(*,*) 'dznew:',(psnowdzn(jjj),jjj=1, INLVLS)
!                   stop
!           endif
    ENDIF
  ENDIF
!
!code initial vincent IF(OSNOW_METAMO.AND.ZDELTA(JJ+1)==1.0)THEN
!code initial vincent         JLAYER1=JJ
!code initial vincent         JLAYER2=JJ+1
!code initial vincent         CALL SNOW3LAGREG(PSNOWDZN,PSNOWDZ,PSNOWRHO,PSNOWGRAN1,PSNOWGRAN2,    &
!code initial vincent                        PSNOWHIST,ZSNOWGRAN1,ZSNOWGRAN2,ZSNOWHIST,            &
!code initial vincent                        JLAYER1,JLAYER2,ZSNOWDDZ)
!code initial vincent ENDIF

!plm
  IF ( JJ<INLVLS .AND. ZDELTA(JJ)==1.0 ) THEN
    JLAYER1 = JJ - 1
    IF ( JJ==1 ) JLAYER1 = JLAYER1 + 1
    JLAYER2 = JJ
    IF ( JJ==1 ) JLAYER2 = JLAYER2 + 1
    CALL SNOW3LAGREG(PSNOWDZN,PSNOWDZ,PSNOWRHO,PSNOWGRAN1,PSNOWGRAN2, &
                     PSNOWHIST,ZSNOWGRAN1,ZSNOWGRAN2,ZSNOWHIST,       &
                     JLAYER1,JLAYER2,ZSNOWDDZ                         )  
  ENDIF
  !plm
  !
  ZSNOWRHON (JJ) = ( PSNOWDZ(JJ) * PSNOWRHO(JJ) ) / PSNOWDZN(JJ)
  ZSNOWHEATN(JJ) = PSNOWHEAT(JJ)
  !
  IF ( JJ>1 ) THEN
    ZSNOWRHON (JJ) = ZSNOWRHON (JJ) - ( ZSNOWDDZ(JJ-1) * ZRHO (JJ-1) ) / PSNOWDZN(JJ) 
    ZSNOWHEATN(JJ) = ZSNOWHEATN(JJ) -   ZSNOWDDZ(JJ-1) * ZHEAT(JJ-1)
  ENDIF
  !
  IF ( JJ<INLVLS ) THEN       
    ZSNOWRHON (JJ) = ZSNOWRHON (JJ) + ( ZSNOWDDZ(JJ) * ZRHO (JJ) ) / PSNOWDZN(JJ)            
    ZSNOWHEATN(JJ) = ZSNOWHEATN(JJ) +   ZSNOWDDZ(JJ) * ZHEAT(JJ) 
  ENDIF
  ! 
ENDDO
!WRITE(*,*) 'AP_AGREG8', PSNOWGRAN1(8)
!
!
! 5. Update mass (density and thickness) and heat:
! ------------------------------------------------
!
PSNOWRHO  (:) = ZSNOWRHON (:)
!PSNOWDZ   (:) = PSNOWDZN  (:)
PSNOWHEAT (:) = ZSNOWHEATN(:)
PSNOWGRAN1(:) = ZSNOWGRAN1(:)
PSNOWGRAN2(:) = ZSNOWGRAN2(:)
PSNOWHIST (:) =  ZSNOWHIST(:)
!
IF (LHOOK) CALL DR_HOOK('SNOWCROTRANSF_1D',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROTRANSF_1D
!################################################################################
!################################################################################
!################################################################################
!
!
SUBROUTINE SNOWNLTRANSFGRID_1D(PSNOW,PSNOWDZ,PSNOWDZN,        &
                               PSNOWRHO,PSNOWHEAT,PSNOWGRAN1, &
                               PSNOWGRAN2, PSNOWHIST          )  
!
!!    PURPOSE
!!    -------
!     Snow mass and heat redistibution due to grid thickness
!     configuration resetting. Total mass and heat content
!     of the overall snowpack unchanged/conserved within this routine.
!
USE MODD_SNOW_PAR, ONLY : XSNOWCRITD, XD1, XD2, XD3, XX, XVALB6
USE MODD_SNOW_METAMO
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                  :: PSNOW
!
REAL, DIMENSION(:), INTENT(INOUT) :: PSNOWHEAT, PSNOWRHO, PSNOWDZ,     &
                                     PSNOWDZN, PSNOWGRAN1, PSNOWGRAN2, &
                                     PSNOWHIST  
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1)) :: ZSNOWRHON,ZSNOWGRAN1N,ZSNOWGRAN2N,   &
                                     ZSNOWHEATN,ZSNOWHISTN,ZSNOWAGEO,ZSNOWAGEN
!
REAL, DIMENSION(SIZE(PSNOWRHO,1)) :: ZSNOWZTOP_OLD, ZSNOWZTOP_NEW, &
                                     ZSNOWZBOT_OLD, ZSNOWZBOT_NEW
!
REAL :: ZPSNOW_OLD, ZPSNOW_NEW
REAL :: ZMASTOTN, ZSNOWHEAN
!
INTEGER :: INLVLS,INLVLS_OLD,INLVLS_NEW
INTEGER :: JST
!                                                          
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWNLTRANSFGRID_1D',0,ZHOOK_HANDLE)
!
INLVLS = SIZE(PSNOWDZ(:),1)
! a ce stade, les couches restent identiques
!
INLVLS_OLD = INLVLS
INLVLS_NEW = INLVLS
ZPSNOW_OLD = PSNOW
ZPSNOW_NEW = PSNOW
!
! 1. Calculate vertical grid limits (m):
! --------------------------------------
!
ZSNOWZTOP_OLD(1) = ZPSNOW_OLD
ZSNOWZTOP_NEW(1) = ZPSNOW_NEW
!
DO JST = 1,INLVLS_OLD
  IF ( JST>1 ) ZSNOWZTOP_OLD(JST) = ZSNOWZBOT_OLD(JST-1)
  ZSNOWZBOT_OLD(JST) = ZSNOWZTOP_OLD(JST) - PSNOWDZ(JST)
ENDDO
!
DO JST = 1,INLVLS_NEW
  IF ( JST>1 ) ZSNOWZTOP_NEW(JST) = ZSNOWZBOT_NEW(JST-1)
  ZSNOWZBOT_NEW(JST) = ZSNOWZTOP_NEW(JST) - PSNOWDZN(JST)
ENDDO
!
! 3. Calculate mass and heat transfers due to grid adjustment/changes:
! --------------------------------------------------------------------
!
! on boucle sur les couches de la nouvelle grille et pour chaque couche
! on somme ou on moyenne les quantités des couches totales ou partielles
! des couches anciennes qui la constituent
!
ZSNOWAGEO(:) = 0.
 CALL GET_MASS_HEAT(1,INLVLS_NEW,INLVLS_OLD,                                   &
                    ZSNOWZTOP_OLD,ZSNOWZTOP_NEW,ZSNOWZBOT_OLD,ZSNOWZBOT_NEW, &
                    PSNOWRHO,PSNOWDZ,PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,        &
                    ZSNOWAGEO,PSNOWHEAT,                                     &
                    ZSNOWRHON,PSNOWDZN,ZSNOWGRAN1N,ZSNOWGRAN2N,ZSNOWHISTN,   &
                    ZSNOWAGEN,ZSNOWHEATN                                     )
!     
!verifs
WRITE(*,*) 'verifs chgt grille INLVLS=',INLVLS
!
ZSNOWHEAN = 0.
ZMASTOTN  = 0.

DO JST=1,INLVLS
  ZSNOWHEAN = ZSNOWHEAN + PSNOWHEAT(JST) - ZSNOWHEATN(JST)
  ZMASTOTN  = ZMASTOTN  + PSNOWRHO(JST)*PSNOWDZ(JST) - ZSNOWRHON(JST)*PSNOWDZN(JST)
  WRITE(*,*) JST,'DZ', PSNOWDZ   (JST), PSNOWDZN   (JST)
  WRITE(*,*)    'RHO', PSNOWRHO  (JST), ZSNOWRHON  (JST)
  WRITE(*,*)   'HEAT', PSNOWHEAT (JST), ZSNOWHEATN (JST)
  WRITE(*,*)    'GR1', PSNOWGRAN1(JST), ZSNOWGRAN1N(JST)
  WRITE(*,*)    'GR2', PSNOWGRAN2(JST), ZSNOWGRAN2N(JST)
  WRITE(*,*)   'HIST', PSNOWHIST (JST), ZSNOWHISTN (JST)
ENDDO
!
WRITE(*,*) 'diff', ZSNOWHEAN,ZMASTOTN
!
! 5. Update mass (density and thickness) and heat:
! ------------------------------------------------
!
PSNOWRHO  (:) = ZSNOWRHON  (:)
PSNOWHEAT (:) = ZSNOWHEATN (:)
PSNOWGRAN1(:) = ZSNOWGRAN1N(:)
PSNOWGRAN2(:) = ZSNOWGRAN2N(:)
PSNOWHIST (:) = ZSNOWHISTN (:)
!
IF (LHOOK) CALL DR_HOOK('SNOWNLTRANSFGRID_1D',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWNLTRANSFGRID_1D
!
!############################################################################
END SUBROUTINE SNOWCROUPGRID
