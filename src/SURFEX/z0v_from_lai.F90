!#######################
MODULE MODI_Z0V_FROM_LAI
!#######################
!
INTERFACE Z0V_FROM_LAI
!
    FUNCTION Z0V_FROM_LAI_0D(PLAI,PH_TREE,PVEGTYPE,OAGRI_TO_GRASS) RESULT(PZ0)
!
REAL,                 INTENT(IN) :: PLAI         ! Leaf area Index
REAL,                 INTENT(IN) :: PH_TREE      ! height of trees
REAL,   DIMENSION(:), INTENT(IN) :: PVEGTYPE     ! type of vegetation
LOGICAL,              INTENT(IN) :: OAGRI_TO_GRASS
!
REAL                             :: PZ0          ! vegetation roughness
!
END FUNCTION Z0V_FROM_LAI_0D
!
!
    FUNCTION Z0V_FROM_LAI_1D(PLAI,PH_TREE,PVEGTYPE,OAGRI_TO_GRASS) RESULT(PZ0)
!
REAL,   DIMENSION(:),   INTENT(IN) :: PLAI         ! Leaf area Index
REAL,   DIMENSION(:),   INTENT(IN) :: PH_TREE      ! height of trees
REAL,   DIMENSION(:,:), INTENT(IN) :: PVEGTYPE     ! type of vegetation
LOGICAL,                INTENT(IN) :: OAGRI_TO_GRASS
!
REAL,   DIMENSION(SIZE(PLAI))      :: PZ0          ! vegetation roughness
!
END FUNCTION Z0V_FROM_LAI_1D
!
!
    FUNCTION Z0V_FROM_LAI_2D(PLAI,PH_TREE,PVEGTYPE,OAGRI_TO_GRASS) RESULT(PZ0)
!
REAL,   DIMENSION(:,:),   INTENT(IN) :: PLAI         ! Leaf area Index
REAL,   DIMENSION(:,:),   INTENT(IN) :: PH_TREE      ! height of trees
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PVEGTYPE     ! type of vegetation
LOGICAL,                  INTENT(IN) :: OAGRI_TO_GRASS
!
REAL,   DIMENSION(SIZE(PLAI,1),SIZE(PLAI,2)) :: PZ0  ! vegetation roughness
!
END FUNCTION Z0V_FROM_LAI_2D
!
    FUNCTION Z0V_FROM_LAI_PATCH(PLAI,PH_TREE,PVEGTYPE,OAGRI_TO_GRASS) RESULT(PZ0)
!
REAL,   DIMENSION(:),   INTENT(IN) :: PLAI         ! Leaf area Index
REAL,   DIMENSION(:),   INTENT(IN) :: PH_TREE      ! height of trees
REAL,   DIMENSION(:),   INTENT(IN) :: PVEGTYPE     ! type of vegetation
LOGICAL,                INTENT(IN) :: OAGRI_TO_GRASS
!
REAL,   DIMENSION(SIZE(PLAI)) :: PZ0  ! vegetation roughness
!
END FUNCTION Z0V_FROM_LAI_PATCH
!
END INTERFACE
!
END MODULE MODI_Z0V_FROM_LAI
!

!   ###########################################################
    FUNCTION Z0V_FROM_LAI_0D(PLAI,PH_TREE,PVEGTYPE,OAGRI_TO_GRASS) RESULT(PZ0)
!   ###########################################################
!!
!!    PURPOSE
!!    -------
!
!     Calculates vegetation roughness from leaf
!    area index and type of vegetation
!    (most of types; forest and vineyards; grassland)
!              
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    none
!!
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Masson and A. Boone          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/03/99
!!      P. Aumond   10/10/2009     Because drag force applied in atmospheric 
!!                                 model, Z0tree -> z0grass
!!      R. Alkama    05/2012   : Extantion from 12 to 19 vegtypes 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_DATA_COVER_PAR, ONLY : NVT_NO, NVT_ROCK, NVT_SNOW, NVT_TEBD,     & 
                                  NVT_BONE, NVT_TRBE, NVT_C3, NVT_C4,     &
                                  NVT_IRR, NVT_GRAS, NVT_TROG,NVT_PARK,   &
                                  NVT_TRBD, NVT_TEBE, NVT_TENE, NVT_BOBD, &
                                  NVT_BOND, NVT_BOGR, NVT_SHRB  
USE MODD_TREEDRAG,       ONLY : LTREEDRAG
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,                 INTENT(IN) :: PLAI         ! Leaf area Index
REAL,                 INTENT(IN) :: PH_TREE      ! height of trees
REAL,   DIMENSION(:), INTENT(IN) :: PVEGTYPE     ! type of vegetation
LOGICAL,              INTENT(IN) :: OAGRI_TO_GRASS
!
REAL                             :: PZ0          ! vegetation roughness
!
!*      0.2    declarations of local variables
!
REAL                            :: ZALLEN_H    ! Allen formula for height
REAL                            :: ZLAI        ! LAI for vegetated areas
!
REAL, DIMENSION(SIZE(PVEGTYPE)) :: ZH          ! height for each type
REAL                            :: ZAVG_H      ! averaged height
REAL                            :: ZZREF       ! reference height        
!
INTEGER                         :: JTYPE       ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_Z0V_FROM_LAI:Z0V_FROM_LAI_0D',0,ZHOOK_HANDLE)
ZLAI = PLAI
IF ( PVEGTYPE(NVT_NO  ) + PVEGTYPE(NVT_ROCK) + PVEGTYPE(NVT_SNOW) < 1.) THEN
  ZLAI = PLAI / (1.-PVEGTYPE(NVT_NO)-PVEGTYPE(NVT_ROCK)-PVEGTYPE(NVT_SNOW))
END IF
!
ZALLEN_H = 0.
IF ( PLAI /= XUNDEF) THEN
  ZALLEN_H = EXP((ZLAI-3.5)/(1.3))
END IF
!
ZH(NVT_PARK) = ZLAI / 6.                    ! irr. grassland
IF (LTREEDRAG) THEN
  ZH(NVT_TEBD) = ZLAI / 6.                  ! forest
  ZH(NVT_BONE) = ZLAI / 6.                  ! forest
  ZH(NVT_TRBE) = ZLAI / 6.                  ! forest
  ZH(NVT_TRBD) = ZLAI / 6.                  ! forest
  ZH(NVT_TEBE) = ZLAI / 6.                  ! forest
  ZH(NVT_TENE) = ZLAI / 6.                  ! forest
  ZH(NVT_BOBD) = ZLAI / 6.                  ! forest
  ZH(NVT_BOND) = ZLAI / 6.                  ! forest
  ZH(NVT_SHRB) = ZLAI / 6.                  ! forest  
ELSE
  ZH(NVT_TEBD) = PH_TREE                  ! forest
  ZH(NVT_BONE) = PH_TREE                  ! forest
  ZH(NVT_TRBE) = PH_TREE                  ! forest
  ZH(NVT_TRBD) = PH_TREE                  ! forest
  ZH(NVT_TEBE) = PH_TREE                  ! forest
  ZH(NVT_TENE) = PH_TREE                  ! forest
  ZH(NVT_BOBD) = PH_TREE                  ! forest
  ZH(NVT_BOND) = PH_TREE                  ! forest
  ZH(NVT_SHRB) = PH_TREE                  ! forest  
END IF
ZH(NVT_GRAS) = ZLAI / 6.                    ! grassland
ZH(NVT_BOGR) = ZLAI / 6.                    ! boreal grassland
ZH(NVT_TROG) = ZLAI / 6.                    ! tropical grassland
IF(OAGRI_TO_GRASS)THEN
  ZH(NVT_C3  ) = ZLAI / 6.
  ZH(NVT_C4  ) = ZLAI / 6.
  ZH(NVT_IRR ) = ZLAI / 6.
ELSE
  ZH(NVT_C3  ) = MIN(1. , ZALLEN_H )          ! cultures
  ZH(NVT_C4  ) = MIN(2.5, ZALLEN_H )          ! C4 types
  ZH(NVT_IRR ) = MIN(2.5, ZALLEN_H )          ! irrigated crops (as C4)
ENDIF
ZH(NVT_NO  ) = 0.1                          ! no vegetation (smooth)
ZH(NVT_ROCK) = 1.                           ! no vegetation (rocks)
ZH(NVT_SNOW) = 0.01                         ! no vegetation (snow)
!
ZH(:) = MAX(ZH(:),0.001)
!
ZZREF = 10.
ZAVG_H = 0.
DO JTYPE=1,SIZE(PVEGTYPE)
  ZAVG_H = ZAVG_H + PVEGTYPE(JTYPE) / (LOG(0.13*ZH(JTYPE)/ZZREF))**2
END DO
ZAVG_H = MAX(ZAVG_H,0.00001)

ZAVG_H = ZZREF / 0.13 * EXP (-1./SQRT(ZAVG_H))
!
PZ0  = MAX(0.001, 0.13*ZAVG_H)
IF (LHOOK) CALL DR_HOOK('MODI_Z0V_FROM_LAI:Z0V_FROM_LAI_0D',1,ZHOOK_HANDLE)
!
!-----------------------------------------------------------------
!
END FUNCTION Z0V_FROM_LAI_0D
!
!   ###########################################################
    FUNCTION Z0V_FROM_LAI_1D(PLAI,PH_TREE,PVEGTYPE,OAGRI_TO_GRASS) RESULT(PZ0)
!   ###########################################################
!!
!!    PURPOSE
!!    -------
!
!     Calculates vegetation roughness from leaf
!    area index and type of vegetation
!    (most of types; forest and vineyards; grassland)
!              
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    none
!!
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Masson and A. Boone          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/03/99
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_DATA_COVER_PAR, ONLY : NVT_NO, NVT_ROCK, NVT_SNOW, NVT_TEBD,     & 
                                  NVT_BONE, NVT_TRBE, NVT_C3, NVT_C4,     &
                                  NVT_IRR, NVT_GRAS, NVT_TROG,NVT_PARK,   &
                                  NVT_TRBD, NVT_TEBE, NVT_TENE, NVT_BOBD, &
                                  NVT_BOND, NVT_BOGR, NVT_SHRB  
USE MODD_TREEDRAG,       ONLY : LTREEDRAG
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:),   INTENT(IN) :: PLAI         ! Leaf area Index
REAL,   DIMENSION(:),   INTENT(IN) :: PH_TREE      ! height of trees
REAL,   DIMENSION(:,:), INTENT(IN) :: PVEGTYPE     ! type of vegetation
LOGICAL,                INTENT(IN) :: OAGRI_TO_GRASS
!
REAL,   DIMENSION(SIZE(PLAI))      :: PZ0          ! vegetation roughness
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PLAI))                  :: ZALLEN_H    ! Allen formula for height
REAL, DIMENSION(SIZE(PLAI))                  :: ZLAI        ! LAI for vegetated areas
!
REAL, DIMENSION(SIZE(PLAI),SIZE(PVEGTYPE,2)) :: ZH          ! height for each type
REAL, DIMENSION(SIZE(PLAI))                  :: ZAVG_H      ! averaged height
REAL                                         :: ZZREF       ! reference height        
!
INTEGER                                      :: JTYPE       ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_Z0V_FROM_LAI:Z0V_FROM_LAI_1D',0,ZHOOK_HANDLE)
ZH(:,:) = XUNDEF
!
ZLAI(:) = PLAI(:)
WHERE ( PVEGTYPE(:,NVT_NO  ) + PVEGTYPE(:,NVT_ROCK) + PVEGTYPE(:,NVT_SNOW) < 1.) 
  ZLAI(:) = PLAI(:) / (1.-PVEGTYPE(:,NVT_NO)-PVEGTYPE(:,NVT_ROCK)-PVEGTYPE(:,NVT_SNOW))
END WHERE
!
ZALLEN_H(:) = 0.
WHERE (PLAI(:) /= XUNDEF)
  ZALLEN_H(:) = EXP((ZLAI(:)-3.5)/(1.3))
END WHERE
!
!
ZH(:,NVT_PARK) = ZLAI(:) / 6.                 ! irr. grassland
IF (LTREEDRAG) THEN
  ZH(:,NVT_TEBD) = ZLAI(:) / 6.         ! forest
  ZH(:,NVT_BONE) = ZLAI(:) / 6.         ! forest
  ZH(:,NVT_TRBE) = ZLAI(:) / 6.         ! forest
  ZH(:,NVT_TRBD) = ZLAI(:) / 6.         ! forest
  ZH(:,NVT_TEBE) = ZLAI(:) / 6.         ! forest
  ZH(:,NVT_TENE) = ZLAI(:) / 6.         ! forest
  ZH(:,NVT_BOBD) = ZLAI(:) / 6.         ! forest
  ZH(:,NVT_BOND) = ZLAI(:) / 6.         ! forest
  ZH(:,NVT_SHRB) = ZLAI(:) / 6.         ! forest  
ELSE
  ZH(:,NVT_TEBD) = PH_TREE(:)           ! forest
  ZH(:,NVT_BONE) = PH_TREE(:)           ! forest
  ZH(:,NVT_TRBE) = PH_TREE(:)           ! forest
  ZH(:,NVT_TRBD) = PH_TREE(:)           ! forest
  ZH(:,NVT_TEBE) = PH_TREE(:)           ! forest
  ZH(:,NVT_TENE) = PH_TREE(:)           ! forest
  ZH(:,NVT_BOBD) = PH_TREE(:)           ! forest
  ZH(:,NVT_BOND) = PH_TREE(:)           ! forest
  ZH(:,NVT_SHRB) = PH_TREE(:)           ! forest  
END IF
ZH(:,NVT_GRAS) = ZLAI(:) / 6.           ! grassland
ZH(:,NVT_BOGR) = ZLAI(:) / 6.           ! boreal grassland
ZH(:,NVT_TROG) = ZLAI(:) / 6.           ! tropical grassland
IF(OAGRI_TO_GRASS)THEN
  ZH(:,NVT_C3  ) = ZLAI(:) / 6.
  ZH(:,NVT_C4  ) = ZLAI(:) / 6.
  ZH(:,NVT_IRR ) = ZLAI(:) / 6.
ELSE
  ZH(:,NVT_C3  ) = MIN(1. , ZALLEN_H(:) )          ! cultures
  ZH(:,NVT_C4  ) = MIN(2.5, ZALLEN_H(:) )          ! C4 types
  ZH(:,NVT_IRR ) = MIN(2.5, ZALLEN_H(:) )          ! irrigated crops (as C4)
ENDIF
ZH(:,NVT_NO  ) = 0.1                    ! no vegetation (smooth)
ZH(:,NVT_ROCK) = 1.                     ! no vegetation (rocks)
ZH(:,NVT_SNOW) = 0.01                   ! no vegetation (snow)
!
ZH(:,:) = MAX(ZH(:,:),0.001)
!
ZZREF = 10.
ZAVG_H(:) = 0.
DO JTYPE=1,SIZE(PVEGTYPE,2)
  ZAVG_H(:) = ZAVG_H(:) + PVEGTYPE(:,JTYPE) / (LOG(0.13*ZH(:,JTYPE)/ZZREF))**2
END DO

ZAVG_H = MAX(ZAVG_H,0.00001)

ZAVG_H(:) = ZZREF / 0.13 * EXP (-1./SQRT(ZAVG_H(:)))
!
PZ0 (:) = MAX(0.001, 0.13*ZAVG_H(:))
IF (LHOOK) CALL DR_HOOK('MODI_Z0V_FROM_LAI:Z0V_FROM_LAI_1D',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
END FUNCTION Z0V_FROM_LAI_1D
!
!   ###########################################################
    FUNCTION Z0V_FROM_LAI_2D(PLAI,PH_TREE,PVEGTYPE,OAGRI_TO_GRASS) RESULT(PZ0)
!   ###########################################################
!!
!!    PURPOSE
!!    -------
!
!     Calculates vegetation roughness from leaf
!    area index and type of vegetation
!    (most of types; forest and vineyards; grassland)
!              
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    none
!!
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Masson and A. Boone          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/03/99
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_DATA_COVER_PAR, ONLY : NVT_NO, NVT_ROCK, NVT_SNOW, NVT_TEBD,     & 
                                  NVT_BONE, NVT_TRBE, NVT_C3, NVT_C4,     &
                                  NVT_IRR, NVT_GRAS, NVT_TROG,NVT_PARK,   &
                                  NVT_TRBD, NVT_TEBE, NVT_TENE, NVT_BOBD, &
                                  NVT_BOND, NVT_BOGR, NVT_SHRB  
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_TREEDRAG,       ONLY : LTREEDRAG
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:,:),   INTENT(IN) :: PLAI         ! Leaf area Index
REAL,   DIMENSION(:,:),   INTENT(IN) :: PH_TREE      ! height of trees
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PVEGTYPE     ! type of vegetation
LOGICAL,                  INTENT(IN) :: OAGRI_TO_GRASS
!
REAL,   DIMENSION(SIZE(PLAI,1),SIZE(PLAI,2)) :: PZ0  ! vegetation roughness
!
!*      0.2    declarations of local variables
!

REAL, DIMENSION(SIZE(PLAI,1),SIZE(PLAI,2))                  :: ZALLEN_H ! Allen formula for height
REAL, DIMENSION(SIZE(PLAI,1),SIZE(PLAI,2))                  :: ZLAI     ! LAI for vegetated areas
!
REAL, DIMENSION(SIZE(PLAI,1),SIZE(PLAI,2),SIZE(PVEGTYPE,3)) :: ZH       ! height for each type
REAL, DIMENSION(SIZE(PLAI,1),SIZE(PLAI,2))                  :: ZAVG_H   ! averaged height
REAL                                                        :: ZZREF    ! reference height        
!
INTEGER                                                     :: JTYPE    ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_Z0V_FROM_LAI:Z0V_FROM_LAI_2D',0,ZHOOK_HANDLE)
ZH(:,:,:)=XUNDEF
!
ZLAI(:,:) = PLAI(:,:)
WHERE ( PVEGTYPE(:,:,NVT_NO  ) + PVEGTYPE(:,:,NVT_ROCK) + PVEGTYPE(:,:,NVT_SNOW) < 1.) 
  ZLAI(:,:) = PLAI(:,:) / (1.-PVEGTYPE(:,:,NVT_NO)-PVEGTYPE(:,:,NVT_ROCK)-PVEGTYPE(:,:,NVT_SNOW))
END WHERE
!
ZALLEN_H(:,:) = 0.
WHERE(PLAI(:,:)/=XUNDEF)
  ZALLEN_H(:,:) = EXP((ZLAI(:,:)-3.5)/(1.3))
END WHERE
!
!
ZH(:,:,NVT_PARK) = ZLAI(:,:) / 6.               ! irr. grassland
IF (LTREEDRAG) THEN
  ZH(:,:,NVT_TEBD) = ZLAI(:,:) / 6.         ! forest
  ZH(:,:,NVT_BONE) = ZLAI(:,:) / 6.         ! forest
  ZH(:,:,NVT_TRBE) = ZLAI(:,:) / 6.         ! forest
  ZH(:,:,NVT_TRBD) = ZLAI(:,:) / 6.         ! forest
  ZH(:,:,NVT_TEBE) = ZLAI(:,:) / 6.         ! forest
  ZH(:,:,NVT_TENE) = ZLAI(:,:) / 6.         ! forest
  ZH(:,:,NVT_BOBD) = ZLAI(:,:) / 6.         ! forest
  ZH(:,:,NVT_BOND) = ZLAI(:,:) / 6.         ! forest
  ZH(:,:,NVT_SHRB) = ZLAI(:,:) / 6.         ! forest  
ELSE
  ZH(:,:,NVT_TEBD) = PH_TREE(:,:)           ! forest
  ZH(:,:,NVT_BONE) = PH_TREE(:,:)           ! forest
  ZH(:,:,NVT_TRBE) = PH_TREE(:,:)           ! forest
  ZH(:,:,NVT_TRBD) = PH_TREE(:,:)           ! forest
  ZH(:,:,NVT_TEBE) = PH_TREE(:,:)           ! forest
  ZH(:,:,NVT_TENE) = PH_TREE(:,:)           ! forest
  ZH(:,:,NVT_BOBD) = PH_TREE(:,:)           ! forest
  ZH(:,:,NVT_BOND) = PH_TREE(:,:)           ! forest
  ZH(:,:,NVT_SHRB) = PH_TREE(:,:)           ! forest   
END IF
ZH(:,:,NVT_GRAS) = ZLAI(:,:) / 6.               ! grassland
ZH(:,:,NVT_BOGR) = ZLAI(:,:) / 6.               ! boreal grassland
ZH(:,:,NVT_TROG) = ZLAI(:,:) / 6.               ! tropical grassland
IF(OAGRI_TO_GRASS)THEN
  ZH(:,:,NVT_C3  ) = ZLAI(:,:) / 6.
  ZH(:,:,NVT_C4  ) = ZLAI(:,:) / 6.
  ZH(:,:,NVT_IRR ) = ZLAI(:,:) / 6.
ELSE
  ZH(:,:,NVT_C3  ) = MIN(1. , ZALLEN_H(:,:) )          ! cultures
  ZH(:,:,NVT_C4  ) = MIN(2.5, ZALLEN_H(:,:) )          ! C4 types
  ZH(:,:,NVT_IRR ) = MIN(2.5, ZALLEN_H(:,:) )          ! irrigated crops (as C4)
ENDIF
ZH(:,:,NVT_NO  ) = 0.1                          ! no vegetation (smooth)
ZH(:,:,NVT_ROCK) = 1.                           ! no vegetation (rocks)
ZH(:,:,NVT_SNOW) = 0.01                         ! no vegetation (snow)
!
ZH(:,:,:) = MAX(ZH(:,:,:),0.001)
!
ZZREF = 10.
ZAVG_H(:,:) = 0.
DO JTYPE=1,SIZE(PVEGTYPE,3)
  ZAVG_H(:,:) = ZAVG_H(:,:) + PVEGTYPE(:,:,JTYPE) / (LOG(0.13*ZH(:,:,JTYPE)/ZZREF))**2
END DO
ZAVG_H(:,:) = MAX(ZAVG_H(:,:),0.00001)
ZAVG_H(:,:) = ZZREF / 0.13 * EXP (-1./SQRT(ZAVG_H(:,:)))
!
PZ0 (:,:) = MAX(0.001, 0.13*ZAVG_H(:,:))
!
WHERE (PLAI(:,:) == XUNDEF)
  PZ0(:,:) = XUNDEF
END WHERE
IF (LHOOK) CALL DR_HOOK('MODI_Z0V_FROM_LAI:Z0V_FROM_LAI_2D',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
END FUNCTION Z0V_FROM_LAI_2D
!
!
!
!   ###########################################################
    FUNCTION Z0V_FROM_LAI_PATCH(PLAI,PH_TREE,PVEGTYPE,OAGRI_TO_GRASS) RESULT(PZ0)
!   ###########################################################
!!
!!    PURPOSE
!!    -------
!
!     Calculates vegetation roughness from leaf
!    area index and type of vegetation for each patch
!    (most of types; forest and vineyards; grassland)
!              
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    none
!!
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!        F.Solmon
!!	V. Masson and A. Boone          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/03/99
!!      
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_DATA_COVER_PAR, ONLY : NVT_NO, NVT_ROCK, NVT_SNOW, NVT_TEBD,     & 
                                  NVT_BONE, NVT_TRBE, NVT_C3, NVT_C4,     &
                                  NVT_IRR, NVT_GRAS, NVT_TROG,NVT_PARK,   &
                                  NVT_TRBD, NVT_TEBE, NVT_TENE, NVT_BOBD, &
                                  NVT_BOND, NVT_BOGR, NVT_SHRB
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_TREEDRAG,       ONLY : LTREEDRAG
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:),   INTENT(IN) :: PLAI         ! Leaf area Index
REAL,   DIMENSION(:),   INTENT(IN) :: PH_TREE      ! height of trees
REAL,   DIMENSION(:),   INTENT(IN) :: PVEGTYPE     ! type of vegetation
LOGICAL,                INTENT(IN) :: OAGRI_TO_GRASS
!
REAL,   DIMENSION(SIZE(PLAI))      :: PZ0          ! vegetation roughness
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PLAI)) :: ZALLEN_H    ! Allen formula for height
!
REAL, DIMENSION(SIZE(PLAI)) :: ZH          ! height for each type
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_Z0V_FROM_LAI:Z0V_FROM_LAI_PATCH',0,ZHOOK_HANDLE)
ZH(:) = XUNDEF
!
WHERE (PLAI(:)/= XUNDEF)
  ZALLEN_H(:) = EXP((PLAI(:)-3.5)/(1.3))
END WHERE
!
!
IF (PVEGTYPE(NVT_PARK)>0.) ZH(NVT_PARK) = PLAI(NVT_PARK) / 6.          ! irr. grasslands
IF (LTREEDRAG) THEN
  IF (PVEGTYPE(NVT_TEBD)>0.) ZH(NVT_TEBD) = PLAI(NVT_TEBD) / 6.        ! broadleaf forest
  IF (PVEGTYPE(NVT_BONE)>0.) ZH(NVT_BONE) = PLAI(NVT_BONE) / 6.        ! coniferous forest
  IF (PVEGTYPE(NVT_TRBE)>0.) ZH(NVT_TRBE) = PLAI(NVT_TRBE) / 6.        ! euqatorial forest
  IF (PVEGTYPE(NVT_TRBD)>0.) ZH(NVT_TRBD) = PLAI(NVT_TRBD) / 6.        ! broadleaf forest
  IF (PVEGTYPE(NVT_TEBE)>0.) ZH(NVT_TEBE) = PLAI(NVT_TEBE) / 6.        ! coniferous forest
  IF (PVEGTYPE(NVT_TENE)>0.) ZH(NVT_TENE) = PLAI(NVT_TENE) / 6.        ! euqatorial forest
  IF (PVEGTYPE(NVT_BOBD)>0.) ZH(NVT_BOBD) = PLAI(NVT_BOBD) / 6.        ! broadleaf forest
  IF (PVEGTYPE(NVT_BOND)>0.) ZH(NVT_BOND) = PLAI(NVT_BOND) / 6.        ! coniferous forest
  IF (PVEGTYPE(NVT_SHRB)>0.) ZH(NVT_SHRB) = PLAI(NVT_SHRB) / 6.        ! euqatorial forest  
ELSE
  IF (PVEGTYPE(NVT_TEBD)>0.) ZH(NVT_TEBD) = PH_TREE(NVT_TEBD)          ! broadleaf forest
  IF (PVEGTYPE(NVT_BONE)>0.) ZH(NVT_BONE) = PH_TREE(NVT_BONE)          ! coniferous forest
  IF (PVEGTYPE(NVT_TRBE)>0.) ZH(NVT_TRBE) = PH_TREE(NVT_TRBE)          ! euqatorial forest
  IF (PVEGTYPE(NVT_TRBD)>0.) ZH(NVT_TRBD) = PH_TREE(NVT_TRBD)          ! broadleaf forest
  IF (PVEGTYPE(NVT_TEBE)>0.) ZH(NVT_TEBE) = PH_TREE(NVT_TEBE)          ! coniferous forest
  IF (PVEGTYPE(NVT_TENE)>0.) ZH(NVT_TENE) = PH_TREE(NVT_TENE)          ! euqatorial forest
  IF (PVEGTYPE(NVT_BOBD)>0.) ZH(NVT_BOBD) = PH_TREE(NVT_BOBD)          ! broadleaf forest
  IF (PVEGTYPE(NVT_BOND)>0.) ZH(NVT_BOND) = PH_TREE(NVT_BOND)          ! coniferous forest
  IF (PVEGTYPE(NVT_SHRB)>0.) ZH(NVT_SHRB) = PH_TREE(NVT_SHRB)          ! euqatorial forest  
END IF
IF (PVEGTYPE(NVT_GRAS)>0.) ZH(NVT_GRAS) = PLAI(NVT_GRAS) / 6.          ! grassland
IF (PVEGTYPE(NVT_BOGR)>0.) ZH(NVT_BOGR) = PLAI(NVT_BOGR) / 6.          ! boreal grassland
IF (PVEGTYPE(NVT_TROG)>0.) ZH(NVT_TROG) = PLAI(NVT_TROG) / 6.          ! tropical grassland
IF(OAGRI_TO_GRASS)THEN
  IF (PVEGTYPE(NVT_C3  )>0.) ZH(NVT_C3  ) = PLAI(NVT_C3)  / 6.  ! cultures
  IF (PVEGTYPE(NVT_C4  )>0.) ZH(NVT_C4  ) = PLAI(NVT_C4)  / 6.  ! C4 types
  IF (PVEGTYPE(NVT_IRR )>0.) ZH(NVT_IRR ) = PLAI(NVT_IRR) / 6.  ! irrigated crops (as C4)
ELSE
  IF (PVEGTYPE(NVT_C3  )>0.) ZH(NVT_C3  ) = MIN(1. , ZALLEN_H(NVT_C3) )  ! cultures
  IF (PVEGTYPE(NVT_C4  )>0.) ZH(NVT_C4  ) = MIN(2.5, ZALLEN_H(NVT_C4) )  ! C4 types
  IF (PVEGTYPE(NVT_IRR )>0.) ZH(NVT_IRR ) = MIN(2.5, ZALLEN_H(NVT_IRR) ) ! irrigated crops (as C4)
ENDIF
IF (PVEGTYPE(NVT_NO  )>0.) ZH(NVT_NO  ) = 0.1                          ! no vegetation (smooth)
IF (PVEGTYPE(NVT_ROCK)>0.) ZH(NVT_ROCK) = 1.                           ! no vegetation (rocks)
IF (PVEGTYPE(NVT_SNOW)>0.) ZH(NVT_SNOW) = 0.01                         ! no vegetation (snow)
!
ZH(:) = MAX(ZH(:),0.001)
!
PZ0 (:) = MAX(0.001, 0.13*ZH(:)) ! rugosite pour chaque vegtype
!-----------------------------------------------------------------
!
WHERE (PLAI(:) == XUNDEF)
  PZ0(:) = XUNDEF
END WHERE
IF (LHOOK) CALL DR_HOOK('MODI_Z0V_FROM_LAI:Z0V_FROM_LAI_PATCH',1,ZHOOK_HANDLE)
!
!
END FUNCTION Z0V_FROM_LAI_PATCH
!
