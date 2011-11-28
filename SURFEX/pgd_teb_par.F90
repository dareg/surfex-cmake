!     #########
      SUBROUTINE PGD_TEB_PAR(HPROGRAM,KROOF_LAYER,KROAD_LAYER,KWALL_LAYER,OGARDEN)
!     ##############################################################
!
!!**** *PGD_TEB_PAR* monitor for averaging and interpolations of cover fractions
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    10/12/97
!!
!!       Modified 08/12/05, P. Le Moigne: user defined fields
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_TEB_GRID_n, ONLY : NDIM
USE MODD_DATA_TEB_n, ONLY : XPAR_Z0_TOWN, XPAR_BLD, XPAR_ALB_ROOF,       &
                              XPAR_EMIS_ROOF, XPAR_HC_ROOF, XPAR_TC_ROOF,  &
                              XPAR_D_ROOF, XPAR_ALB_ROAD, XPAR_EMIS_ROAD,  &
                              XPAR_HC_ROAD, XPAR_TC_ROAD, XPAR_D_ROAD,     &
                              XPAR_ALB_WALL, XPAR_EMIS_WALL, XPAR_HC_WALL, &
                              XPAR_TC_WALL, XPAR_D_WALL, XPAR_BLD_HEIGHT,  &
                              XPAR_WALL_O_HOR,                             &
                              XPAR_H_TRAFFIC, XPAR_LE_TRAFFIC,           &
                              XPAR_H_INDUSTRY, XPAR_LE_INDUSTRY ,        &
                              XPAR_VEG_ROOF, XPAR_GARDEN, XPAR_URBTYPE,  &
                              LDATA_Z0_TOWN, LDATA_BLD, LDATA_ALB_ROOF,       &
                              LDATA_EMIS_ROOF, LDATA_HC_ROOF, LDATA_TC_ROOF,  &
                              LDATA_D_ROOF, LDATA_ALB_ROAD, LDATA_EMIS_ROAD,  &
                              LDATA_HC_ROAD, LDATA_TC_ROAD, LDATA_D_ROAD,     &
                              LDATA_ALB_WALL, LDATA_EMIS_WALL, LDATA_HC_WALL, &
                              LDATA_TC_WALL, LDATA_D_WALL, LDATA_BLD_HEIGHT,  &
                              LDATA_WALL_O_HOR,                             &
                              LDATA_H_TRAFFIC, LDATA_LE_TRAFFIC,           &
                              LDATA_H_INDUSTRY, LDATA_LE_INDUSTRY ,        &
                              LDATA_VEG_ROOF, LDATA_GARDEN, LDATA_URBTYPE
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_INI_VAR_FROM_DATA_0D
USE MODI_INI_VAR_FROM_DATA
!
USE MODE_POS_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=6),    INTENT(IN)    :: HPROGRAM     ! Type of program
INTEGER,             INTENT(INOUT) :: KROOF_LAYER  ! number of roof layers
INTEGER,             INTENT(INOUT) :: KROAD_LAYER  ! number of road layers
INTEGER,             INTENT(INOUT) :: KWALL_LAYER  ! number of wall layers
LOGICAL,             INTENT(IN)    :: OGARDEN      ! T if urban green areas
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER               :: ILUOUT    ! output listing logical unit
INTEGER               :: ILUNAM    ! namelist file  logical unit
LOGICAL               :: GFOUND    ! true if namelist is found
!
INTEGER               :: JLAYER    ! loop counter on layers
!
!*    0.3    Declaration of namelists
!            ------------------------
!
INTEGER            :: NROOF_LAYER
INTEGER            :: NROAD_LAYER
INTEGER            :: NWALL_LAYER
INTEGER, PARAMETER :: NROOF_MAX  = 9
INTEGER, PARAMETER :: NROAD_MAX  = 9
INTEGER, PARAMETER :: NWALL_MAX  = 9
!
! Geometric Parameters:
!
REAL                                    :: XUNIF_URBTYPE
CHARACTER(LEN=28)                       :: CFNAM_URBTYPE
CHARACTER(LEN=6)                        :: CFTYP_URBTYPE
!
REAL                                    :: XUNIF_BLD          ! fraction of buildings            (-)
REAL                                    :: XUNIF_BLD_HEIGHT   ! buildings height 'h'             (m)
REAL                                    :: XUNIF_WALL_O_HOR   ! wall surf. / hor. surf.          (-)
REAL                                    :: XUNIF_Z0_TOWN      ! roughness length for momentum    (m)
REAL                                    :: XUNIF_GARDEN       ! fraction of veg in the streets   (-)
REAL                                    :: XUNIF_VEG_ROOF     ! fraction of veg on roofs         (-)
CHARACTER(LEN=28)                       :: CFNAM_BLD          ! file name for BLD 
CHARACTER(LEN=28)                       :: CFNAM_BLD_HEIGHT   ! file name for BLD_HEIGHT
CHARACTER(LEN=28)                       :: CFNAM_WALL_O_HOR   ! file name for WALL_O_HOR
CHARACTER(LEN=28)                       :: CFNAM_Z0_TOWN      ! file name for Z0_TOWN
CHARACTER(LEN=28)                       :: CFNAM_GARDEN       ! file name for GARDEN  
CHARACTER(LEN=28)                       :: CFNAM_VEG_ROOF     ! file name for VEG_ROOF
CHARACTER(LEN=6)                        :: CFTYP_BLD          ! file type for BLD 
CHARACTER(LEN=6)                        :: CFTYP_BLD_HEIGHT   ! file type for BLD_HEIGHT
CHARACTER(LEN=6)                        :: CFTYP_WALL_O_HOR   ! file type for WALL_O_HOR
CHARACTER(LEN=6)                        :: CFTYP_Z0_TOWN      ! file type for Z0_TOWN
CHARACTER(LEN=6)                        :: CFTYP_GARDEN       ! file type for GARDEN  
CHARACTER(LEN=6)                        :: CFTYP_VEG_ROOF     ! file type for VEG_ROOF
!
! Roof parameters
!
REAL                                    :: XUNIF_ALB_ROOF     ! roof albedo                      (-)
REAL                                    :: XUNIF_EMIS_ROOF    ! roof emissivity                  (-)
CHARACTER(LEN=28)                       :: CFNAM_ALB_ROOF     ! file name for ALB_ROOF
CHARACTER(LEN=28)                       :: CFNAM_EMIS_ROOF    ! file name for EMIS_ROOF
CHARACTER(LEN=6)                        :: CFTYP_ALB_ROOF     ! file name for ALB_ROOF   
CHARACTER(LEN=6)                        :: CFTYP_EMIS_ROOF    ! file name for EMIS_ROOF  
REAL, DIMENSION(NROOF_MAX)              :: XUNIF_HC_ROOF      ! roof layers heat capacity        (J/K/m3)
REAL, DIMENSION(NROOF_MAX)              :: XUNIF_TC_ROOF      ! roof layers thermal conductivity (W/K/m)
REAL, DIMENSION(NROOF_MAX)              :: XUNIF_D_ROOF       ! depth of roof layers             (m)
CHARACTER(LEN=28), DIMENSION(NROOF_MAX) :: CFNAM_HC_ROOF      ! file name for HC_ROOF   
CHARACTER(LEN=28), DIMENSION(NROOF_MAX) :: CFNAM_TC_ROOF      ! file name for TC_ROOF
CHARACTER(LEN=28), DIMENSION(NROOF_MAX) :: CFNAM_D_ROOF       ! file name for D_ROOF
CHARACTER(LEN=6), DIMENSION(NROOF_MAX) :: CFTYP_HC_ROOF      ! file type for HC_ROOF   
CHARACTER(LEN=6), DIMENSION(NROOF_MAX) :: CFTYP_TC_ROOF      ! file type for TC_ROOF
CHARACTER(LEN=6), DIMENSION(NROOF_MAX) :: CFTYP_D_ROOF       ! file type for D_ROOF
!
!
! Road parameters
!
REAL                                    :: XUNIF_ALB_ROAD     ! road albedo                      (-)
REAL                                    :: XUNIF_EMIS_ROAD    ! road emissivity                  (-)
CHARACTER(LEN=28)                       :: CFNAM_ALB_ROAD     ! file name for ALB_ROAD
CHARACTER(LEN=28)                       :: CFNAM_EMIS_ROAD    ! file name for EMIS_ROAD
CHARACTER(LEN=6)                        :: CFTYP_ALB_ROAD     ! file type for ALB_ROAD
CHARACTER(LEN=6)                        :: CFTYP_EMIS_ROAD    ! file type for EMIS_ROAD
REAL, DIMENSION(NROAD_MAX)              :: XUNIF_HC_ROAD      ! road layers heat capacity        (J/K/m3)
REAL, DIMENSION(NROAD_MAX)              :: XUNIF_TC_ROAD      ! road layers thermal conductivity (W/K/m)
REAL, DIMENSION(NROAD_MAX)              :: XUNIF_D_ROAD       ! depth of road layers             (m)
CHARACTER(LEN=28), DIMENSION(NROAD_MAX) :: CFNAM_HC_ROAD      ! file name for HC_ROAD   
CHARACTER(LEN=28), DIMENSION(NROAD_MAX) :: CFNAM_TC_ROAD      ! file name for TC_ROAD
CHARACTER(LEN=28), DIMENSION(NROAD_MAX) :: CFNAM_D_ROAD       ! file name for D_ROAD
CHARACTER(LEN=6), DIMENSION(NROAD_MAX) :: CFTYP_HC_ROAD      ! file type for HC_ROAD   
CHARACTER(LEN=6), DIMENSION(NROAD_MAX) :: CFTYP_TC_ROAD      ! file type for TC_ROAD
CHARACTER(LEN=6), DIMENSION(NROAD_MAX) :: CFTYP_D_ROAD       ! file type for D_ROAD
!
! Wall parameters
!
REAL                                    :: XUNIF_ALB_WALL     ! wall albedo                      (-)
REAL                                    :: XUNIF_EMIS_WALL    ! wall emissivity                  (-)
CHARACTER(LEN=28)                       :: CFNAM_ALB_WALL     ! file name for ALB_WALL
CHARACTER(LEN=28)                       :: CFNAM_EMIS_WALL    ! file name for EMIS_WALL
CHARACTER(LEN=6)                        :: CFTYP_ALB_WALL     ! file type for ALB_WALL
CHARACTER(LEN=6)                        :: CFTYP_EMIS_WALL    ! file type for EMIS_WALL
REAL, DIMENSION(NWALL_MAX)              :: XUNIF_HC_WALL      ! wall layers heat capacity        (J/K/m3)
REAL, DIMENSION(NWALL_MAX)              :: XUNIF_TC_WALL      ! wall layers thermal conductivity (W/K/m)
REAL, DIMENSION(NWALL_MAX)              :: XUNIF_D_WALL       ! depth of wall layers             (m)
CHARACTER(LEN=28), DIMENSION(NWALL_MAX) :: CFNAM_HC_WALL      ! file name for HC_WALL   
CHARACTER(LEN=28), DIMENSION(NWALL_MAX) :: CFNAM_TC_WALL      ! file name for TC_WALL
CHARACTER(LEN=28), DIMENSION(NWALL_MAX) :: CFNAM_D_WALL       ! file name for D_WALL
CHARACTER(LEN=6), DIMENSION(NWALL_MAX) :: CFTYP_HC_WALL      ! file type for HC_WALL   
CHARACTER(LEN=6), DIMENSION(NWALL_MAX) :: CFTYP_TC_WALL      ! file type for TC_WALL
CHARACTER(LEN=6), DIMENSION(NWALL_MAX) :: CFTYP_D_WALL       ! file type for D_WALL
!
! anthropogenic fluxes
!
REAL                                    :: XUNIF_H_TRAFFIC    ! anthropogenic sensible
!                                                             ! heat fluxes due to traffic       (W/m2)
REAL                                    :: XUNIF_LE_TRAFFIC   ! anthropogenic latent
!                                                             ! heat fluxes due to traffic       (W/m2)
REAL                                    :: XUNIF_H_INDUSTRY   ! anthropogenic sensible                   
!                                                             ! heat fluxes due to factories     (W/m2)
REAL                                    :: XUNIF_LE_INDUSTRY  ! anthropogenic latent
!                                                             ! heat fluxes due to factories     (W/m2)
CHARACTER(LEN=28)                       :: CFNAM_H_TRAFFIC    ! file name for H_TRAFFIC
CHARACTER(LEN=28)                       :: CFNAM_LE_TRAFFIC   ! file name for LE_TRAFFIC
CHARACTER(LEN=28)                       :: CFNAM_H_INDUSTRY   ! file name for H_INDUSTRY
CHARACTER(LEN=28)                       :: CFNAM_LE_INDUSTRY  ! file name for LE_INDUSTRY
CHARACTER(LEN=6)                        :: CFTYP_H_TRAFFIC    ! file type for H_TRAFFIC
CHARACTER(LEN=6)                        :: CFTYP_LE_TRAFFIC   ! file type for LE_TRAFFIC
CHARACTER(LEN=6)                        :: CFTYP_H_INDUSTRY   ! file type for H_INDUSTRY
CHARACTER(LEN=6)                        :: CFTYP_LE_INDUSTRY  ! file type for LE_INDUSTRY
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!

NAMELIST/NAM_DATA_TEB/        NROOF_LAYER, NROAD_LAYER, NWALL_LAYER,   &
                              XUNIF_URBTYPE, CFNAM_URBTYPE, CFTYP_URBTYPE,    &
                              XUNIF_ALB_ROOF,                                 &
                              XUNIF_EMIS_ROOF, XUNIF_HC_ROOF, XUNIF_TC_ROOF,  &
                              XUNIF_D_ROOF, XUNIF_ALB_ROAD, XUNIF_EMIS_ROAD,  &
                              XUNIF_HC_ROAD, XUNIF_TC_ROAD, XUNIF_D_ROAD,     &
                              XUNIF_ALB_WALL, XUNIF_EMIS_WALL, XUNIF_HC_WALL, &
                              XUNIF_TC_WALL, XUNIF_D_WALL,                    &
                              XUNIF_Z0_TOWN, XUNIF_BLD, XUNIF_BLD_HEIGHT,     &
                              XUNIF_WALL_O_HOR,                               &
                              XUNIF_H_TRAFFIC, XUNIF_LE_TRAFFIC,              &
                              XUNIF_H_INDUSTRY, XUNIF_LE_INDUSTRY,            &
                              XUNIF_GARDEN, XUNIF_VEG_ROOF,                   &
                              CFNAM_ALB_ROOF,                                 &
                              CFNAM_EMIS_ROOF, CFNAM_HC_ROOF, CFNAM_TC_ROOF,  &
                              CFNAM_D_ROOF, CFNAM_ALB_ROAD, CFNAM_EMIS_ROAD,  &
                              CFNAM_HC_ROAD, CFNAM_TC_ROAD, CFNAM_D_ROAD,     &
                              CFNAM_ALB_WALL, CFNAM_EMIS_WALL, CFNAM_HC_WALL, &
                              CFNAM_TC_WALL, CFNAM_D_WALL,                    &
                              CFNAM_Z0_TOWN, CFNAM_BLD, CFNAM_BLD_HEIGHT,     &
                              CFNAM_WALL_O_HOR,                               &
                              CFNAM_H_TRAFFIC, CFNAM_LE_TRAFFIC,              &
                              CFNAM_H_INDUSTRY, CFNAM_LE_INDUSTRY,            &
                              CFNAM_GARDEN, CFNAM_VEG_ROOF,                   &
                              CFTYP_ALB_ROOF,                                 &
                              CFTYP_EMIS_ROOF, CFTYP_HC_ROOF, CFTYP_TC_ROOF,  &
                              CFTYP_D_ROOF, CFTYP_ALB_ROAD, CFTYP_EMIS_ROAD,  &
                              CFTYP_HC_ROAD, CFTYP_TC_ROAD, CFTYP_D_ROAD,     &
                              CFTYP_ALB_WALL, CFTYP_EMIS_WALL, CFTYP_HC_WALL, &
                              CFTYP_TC_WALL, CFTYP_D_WALL,                    &
                              CFTYP_Z0_TOWN, CFTYP_BLD, CFTYP_BLD_HEIGHT,     &
                              CFTYP_WALL_O_HOR,                               &
                              CFTYP_H_TRAFFIC, CFTYP_LE_TRAFFIC,              &
                              CFTYP_H_INDUSTRY, CFTYP_LE_INDUSTRY,            &
                              CFTYP_GARDEN, CFTYP_VEG_ROOF  
!
!-------------------------------------------------------------------------------
!
!*    1.      Initializations
!             ---------------
!
IF (LHOOK)   CALL DR_HOOK('PGD_TEB_PAR',0,ZHOOK_HANDLE)
NROOF_LAYER=3
NROAD_LAYER=3
NWALL_LAYER=3
XUNIF_URBTYPE      = XUNDEF
XUNIF_BLD          = XUNDEF
XUNIF_BLD_HEIGHT   = XUNDEF
XUNIF_WALL_O_HOR   = XUNDEF
XUNIF_Z0_TOWN      = XUNDEF
XUNIF_ALB_ROOF     = XUNDEF
XUNIF_EMIS_ROOF    = XUNDEF
XUNIF_HC_ROOF      = XUNDEF
XUNIF_TC_ROOF      = XUNDEF
XUNIF_D_ROOF       = XUNDEF
XUNIF_ALB_ROAD     = XUNDEF
XUNIF_EMIS_ROAD    = XUNDEF
XUNIF_HC_ROAD      = XUNDEF
XUNIF_TC_ROAD      = XUNDEF
XUNIF_D_ROAD       = XUNDEF
XUNIF_ALB_WALL     = XUNDEF
XUNIF_EMIS_WALL    = XUNDEF
XUNIF_HC_WALL      = XUNDEF
XUNIF_TC_WALL      = XUNDEF
XUNIF_D_WALL       = XUNDEF
XUNIF_H_TRAFFIC    = XUNDEF
XUNIF_LE_TRAFFIC   = XUNDEF
XUNIF_H_INDUSTRY   = XUNDEF
XUNIF_LE_INDUSTRY  = XUNDEF
XUNIF_VEG_ROOF     = XUNDEF
XUNIF_GARDEN       = XUNDEF

CFNAM_URBTYPE      = '                            '
CFNAM_BLD          = '                            '
CFNAM_BLD_HEIGHT   = '                            '
CFNAM_WALL_O_HOR   = '                            '
CFNAM_Z0_TOWN      = '                            '

CFNAM_ALB_ROOF (:) = '                            '
CFNAM_EMIS_ROOF(:) = '                            '
CFNAM_HC_ROOF  (:) = '                            '
CFNAM_TC_ROOF  (:) = '                            '
CFNAM_D_ROOF   (:) = '                            '
CFNAM_ALB_ROAD (:) = '                            '
CFNAM_EMIS_ROAD(:) = '                            '
CFNAM_HC_ROAD  (:) = '                            '
CFNAM_TC_ROAD  (:) = '                            '
CFNAM_D_ROAD   (:) = '                            '
CFNAM_ALB_WALL (:) = '                            '
CFNAM_EMIS_WALL(:) = '                            '
CFNAM_HC_WALL  (:) = '                            '
CFNAM_TC_WALL  (:) = '                            '
CFNAM_D_WALL   (:) = '                            '

CFNAM_H_TRAFFIC    = '                            '
CFNAM_LE_TRAFFIC   = '                            '
CFNAM_H_INDUSTRY   = '                            '
CFNAM_LE_INDUSTRY  = '                            '

CFNAM_GARDEN       = '                            '
CFNAM_VEG_ROOF     = '                            '

CFTYP_URBTYPE      = '      '
CFTYP_BLD          = '      '
CFTYP_BLD_HEIGHT   = '      '
CFTYP_WALL_O_HOR   = '      '
CFTYP_Z0_TOWN      = '      '
CFTYP_ALB_ROOF(:)  = '      '
CFTYP_EMIS_ROOF(:) = '      '
CFTYP_HC_ROOF(:)   = '      '
CFTYP_TC_ROOF(:)   = '      '
CFTYP_D_ROOF(:)    = '      '
CFTYP_ALB_ROAD(:)  = '      '
CFTYP_EMIS_ROAD(:) = '      '
CFTYP_HC_ROAD(:)   = '      '
CFTYP_TC_ROAD(:)   = '      '
CFTYP_D_ROAD(:)    = '      '
CFTYP_ALB_WALL(:)  = '      '
CFTYP_EMIS_WALL(:) = '      '
CFTYP_HC_WALL(:)   = '      '
CFTYP_TC_WALL(:)   = '      '
CFTYP_D_WALL(:)    = '      '
CFTYP_H_TRAFFIC    = '      '
CFTYP_LE_TRAFFIC   = '      '
CFTYP_H_INDUSTRY   = '      '
CFTYP_LE_INDUSTRY  = '      '
CFTYP_GARDEN       = '      '
CFTYP_VEG_ROOF     = '      '
!
!-------------------------------------------------------------------------------
!
!*    2.      Input file for cover types
!             --------------------------
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
CALL OPEN_NAMELIST(HPROGRAM,ILUNAM)
!
CALL POSNAM(ILUNAM,'NAM_DATA_TEB',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_DATA_TEB)
!
CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
!
!-------------------------------------------------------------------------------
KROOF_LAYER = NROOF_LAYER
KROAD_LAYER = NROAD_LAYER
KWALL_LAYER = NWALL_LAYER
!-------------------------------------------------------------------------------
ALLOCATE(XPAR_URBTYPE     (NDIM,1))
ALLOCATE(XPAR_Z0_TOWN     (NDIM))
ALLOCATE(XPAR_ALB_ROOF    (NDIM))
ALLOCATE(XPAR_EMIS_ROOF   (NDIM))
ALLOCATE(XPAR_ALB_ROAD    (NDIM))
ALLOCATE(XPAR_EMIS_ROAD   (NDIM))
ALLOCATE(XPAR_ALB_WALL    (NDIM))
ALLOCATE(XPAR_EMIS_WALL   (NDIM))
ALLOCATE(XPAR_BLD         (NDIM))
ALLOCATE(XPAR_BLD_HEIGHT  (NDIM))
ALLOCATE(XPAR_WALL_O_HOR  (NDIM))
ALLOCATE(XPAR_H_TRAFFIC   (NDIM))
ALLOCATE(XPAR_LE_TRAFFIC  (NDIM))
ALLOCATE(XPAR_H_INDUSTRY  (NDIM))
ALLOCATE(XPAR_LE_INDUSTRY (NDIM))
ALLOCATE(XPAR_GARDEN      (NDIM))
ALLOCATE(XPAR_VEG_ROOF    (NDIM))
!
ALLOCATE(XPAR_HC_ROOF     (NDIM,NROOF_LAYER))
ALLOCATE(XPAR_TC_ROOF     (NDIM,NROOF_LAYER))
ALLOCATE(XPAR_D_ROOF      (NDIM,NROOF_LAYER))
ALLOCATE(XPAR_HC_ROAD     (NDIM,NROAD_LAYER))
ALLOCATE(XPAR_TC_ROAD     (NDIM,NROAD_LAYER))
ALLOCATE(XPAR_D_ROAD      (NDIM,NROAD_LAYER))
ALLOCATE(XPAR_HC_WALL     (NDIM,NWALL_LAYER))
ALLOCATE(XPAR_TC_WALL     (NDIM,NWALL_LAYER))
ALLOCATE(XPAR_D_WALL      (NDIM,NWALL_LAYER))
!
!-------------------------------------------------------------------------------
IF (NROOF_MAX < NROOF_LAYER) THEN
  WRITE(ILUOUT,*) '---------------------------------------------'
  WRITE(ILUOUT,*) 'Please update pgd_teb_par.f90 routine :      '
  WRITE(ILUOUT,*) 'The maximum number of ROOF LAYER             '
  WRITE(ILUOUT,*) 'in the declaration of the namelist variables '
  WRITE(ILUOUT,*) 'must be increased to : ', NROOF_LAYER
  WRITE(ILUOUT,*) '---------------------------------------------'
  CALL ABOR1_SFX('PGD_TEB_PAR: MAXIMUM NUMBER OF NROOF_LAYER MUST BE INCREASED')
ENDIF
!-------------------------------------------------------------------------------
IF (NROAD_MAX < NROAD_LAYER) THEN
  WRITE(ILUOUT,*) '---------------------------------------------'
  WRITE(ILUOUT,*) 'Please update pgd_teb_par.f90 routine :      '
  WRITE(ILUOUT,*) 'The maximum number of ROAD LAYER             '
  WRITE(ILUOUT,*) 'in the declaration of the namelist variables '
  WRITE(ILUOUT,*) 'must be increased to : ', NROAD_LAYER
  WRITE(ILUOUT,*) '---------------------------------------------'
  CALL ABOR1_SFX('PGD_TEB_PAR: MAXIMUM NUMBER OF NROAD_LAYER MUST BE INCREASED')
ENDIF
!-------------------------------------------------------------------------------
IF (NWALL_MAX < NWALL_LAYER) THEN
  WRITE(ILUOUT,*) '---------------------------------------------'
  WRITE(ILUOUT,*) 'Please update pgd_teb_par.f90 routine :      '
  WRITE(ILUOUT,*) 'The maximum number of WALL LAYER             '
  WRITE(ILUOUT,*) 'in the declaration of the namelist variables '
  WRITE(ILUOUT,*) 'must be increased to : ', NWALL_LAYER
  WRITE(ILUOUT,*) '---------------------------------------------'
  CALL ABOR1_SFX('PGD_TEB_PAR: MAXIMUM NUMBER OF NWALL_LAYER MUST BE INCREASED')
ENDIF

!-------------------------------------------------------------------------------
!
!*    3.      user defined fields are prescribed
!             ----------------------------------
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','URBTYPE    ','TWN', CFNAM_URBTYPE,CFTYP_URBTYPE,XUNIF_URBTYPE,&
        XPAR_URBTYPE(:,1),LDATA_URBTYPE )
IF (.NOT.LDATA_URBTYPE) DEALLOCATE(XPAR_URBTYPE)
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','BLD        ','TWN', CFNAM_BLD,CFTYP_BLD,XUNIF_BLD,XPAR_BLD,LDATA_BLD )
IF (.NOT.LDATA_BLD) DEALLOCATE(XPAR_BLD)
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','BLD_HEIGHT ','TWN',CFNAM_BLD_HEIGHT,CFTYP_BLD_HEIGHT,XUNIF_BLD_HEIGHT,&
        XPAR_BLD_HEIGHT,LDATA_BLD_HEIGHT)
IF (.NOT.LDATA_BLD_HEIGHT) DEALLOCATE(XPAR_BLD_HEIGHT)
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','WALL_O_HOR ','TWN',CFNAM_WALL_O_HOR,CFTYP_WALL_O_HOR,XUNIF_WALL_O_HOR,&
        XPAR_WALL_O_HOR,LDATA_WALL_O_HOR)
IF (.NOT.LDATA_WALL_O_HOR) DEALLOCATE(XPAR_WALL_O_HOR)
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','Z0_TOWN    ','TWN',CFNAM_Z0_TOWN,CFTYP_Z0_TOWN,XUNIF_Z0_TOWN,&
        XPAR_Z0_TOWN,LDATA_Z0_TOWN)
IF (.NOT.LDATA_Z0_TOWN) DEALLOCATE(XPAR_Z0_TOWN)
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','ALB_ROOF   ','TWN',CFNAM_ALB_ROOF,CFTYP_ALB_ROOF,XUNIF_ALB_ROOF  ,&
        XPAR_ALB_ROOF,LDATA_ALB_ROOF)
IF (.NOT.LDATA_ALB_ROOF) DEALLOCATE(XPAR_ALB_ROOF)
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','EMIS_ROOF  ','TWN',CFNAM_EMIS_ROOF,CFTYP_EMIS_ROOF,XUNIF_EMIS_ROOF ,&
        XPAR_EMIS_ROOF,LDATA_EMIS_ROOF)
IF (.NOT.LDATA_EMIS_ROOF) DEALLOCATE(XPAR_EMIS_ROOF)
!
CALL INI_VAR_FROM_DATA(HPROGRAM,'INV','HC_ROOF  ','TWN',CFNAM_HC_ROOF,CFTYP_HC_ROOF, &
        XUNIF_HC_ROOF,XPAR_HC_ROOF,LDATA_HC_ROOF ) 
IF (.NOT.LDATA_HC_ROOF) DEALLOCATE(XPAR_HC_ROOF)
! 
CALL INI_VAR_FROM_DATA(HPROGRAM,'ARI','TC_ROOF  ','TWN',CFNAM_TC_ROOF,CFTYP_TC_ROOF, &
                 XUNIF_TC_ROOF ,XPAR_TC_ROOF, LDATA_TC_ROOF ) 
IF (.NOT.LDATA_TC_ROOF) DEALLOCATE(XPAR_TC_ROOF)
! 
CALL INI_VAR_FROM_DATA(HPROGRAM,'ARI','D_ROOF   ','TWN',CFNAM_D_ROOF,CFTYP_D_ROOF, &
                 XUNIF_D_ROOF  ,XPAR_D_ROOF , LDATA_D_ROOF ) 
IF (.NOT.LDATA_D_ROOF) DEALLOCATE(XPAR_D_ROOF)
! 
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','ALB_ROAD   ','TWN',CFNAM_ALB_ROAD  ,CFTYP_ALB_ROAD  ,XUNIF_ALB_ROAD  ,&
        XPAR_ALB_ROAD, LDATA_ALB_ROAD  )
IF (.NOT.LDATA_ALB_ROAD) DEALLOCATE(XPAR_ALB_ROAD)
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','EMIS_ROAD  ','TWN',CFNAM_EMIS_ROAD ,CFTYP_EMIS_ROAD ,XUNIF_EMIS_ROAD ,&
        XPAR_EMIS_ROAD, LDATA_EMIS_ROAD )
IF (.NOT.LDATA_EMIS_ROAD) DEALLOCATE(XPAR_EMIS_ROAD)
!
CALL INI_VAR_FROM_DATA(HPROGRAM,'INV','HC_ROAD  ','TWN',CFNAM_HC_ROAD ,CFTYP_HC_ROAD , &
                   XUNIF_HC_ROAD ,XPAR_HC_ROAD, LDATA_HC_ROAD  )  
IF (.NOT.LDATA_HC_ROAD) DEALLOCATE(XPAR_HC_ROAD)
!
CALL INI_VAR_FROM_DATA(HPROGRAM,'ARI','TC_ROAD  ','TWN',CFNAM_TC_ROAD ,CFTYP_TC_ROAD , &
                   XUNIF_TC_ROAD ,XPAR_TC_ROAD, LDATA_TC_ROAD  )  
IF (.NOT.LDATA_TC_ROAD) DEALLOCATE(XPAR_TC_ROAD)
!
CALL INI_VAR_FROM_DATA(HPROGRAM,'ARI','D_ROAD   ','TWN',CFNAM_D_ROAD  ,CFTYP_D_ROAD  , &
                   XUNIF_D_ROAD  ,XPAR_D_ROAD , LDATA_D_ROAD  )
IF (.NOT.LDATA_D_ROAD) DEALLOCATE(XPAR_D_ROAD)
!  
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','ALB_WALL   ','TWN',CFNAM_ALB_WALL  ,CFTYP_ALB_WALL  ,XUNIF_ALB_WALL  ,&
        XPAR_ALB_WALL, LDATA_ALB_WALL   )
IF (.NOT.LDATA_ALB_WALL) DEALLOCATE(XPAR_ALB_WALL)
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','EMIS_WALL  ','TWN',CFNAM_EMIS_WALL ,CFTYP_EMIS_WALL ,XUNIF_EMIS_WALL ,&
        XPAR_EMIS_WALL, LDATA_EMIS_WALL  )
IF (.NOT.LDATA_EMIS_WALL) DEALLOCATE(XPAR_EMIS_WALL)
!
CALL INI_VAR_FROM_DATA(HPROGRAM,'INV','HC_WALL  ','TWN',CFNAM_HC_WALL ,CFTYP_HC_WALL , &
                   XUNIF_HC_WALL ,XPAR_HC_WALL, LDATA_HC_WALL  ) 
IF (.NOT.LDATA_HC_WALL) DEALLOCATE(XPAR_HC_WALL)
! 
CALL INI_VAR_FROM_DATA(HPROGRAM,'ARI','TC_WALL  ','TWN',CFNAM_TC_WALL ,CFTYP_TC_WALL , &
                   XUNIF_TC_WALL ,XPAR_TC_WALL, LDATA_TC_WALL  ) 
IF (.NOT.LDATA_TC_WALL) DEALLOCATE(XPAR_TC_WALL)
! 
CALL INI_VAR_FROM_DATA(HPROGRAM,'ARI','D_WALL   ','TWN',CFNAM_D_WALL  ,CFTYP_D_WALL  , &
                   XUNIF_D_WALL  ,XPAR_D_WALL , LDATA_D_WALL  ) 
IF (.NOT.LDATA_D_WALL) DEALLOCATE(XPAR_D_WALL)
! 
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','H_TRAFFIC  ','TWN',CFNAM_H_TRAFFIC  ,CFTYP_H_TRAFFIC  ,XUNIF_H_TRAFFIC  ,&
        XPAR_H_TRAFFIC, LDATA_H_TRAFFIC   )
IF (.NOT.LDATA_H_TRAFFIC) DEALLOCATE(XPAR_H_TRAFFIC)
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','LE_TRAFFIC ','TWN',CFNAM_LE_TRAFFIC ,CFTYP_LE_TRAFFIC ,XUNIF_LE_TRAFFIC ,&
        XPAR_LE_TRAFFIC, LDATA_LE_TRAFFIC  )
IF (.NOT.LDATA_LE_TRAFFIC) DEALLOCATE(XPAR_LE_TRAFFIC)
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','H_INDUSTRY ','TWN',CFNAM_H_INDUSTRY ,CFTYP_H_INDUSTRY ,XUNIF_H_INDUSTRY ,&
        XPAR_H_INDUSTRY, LDATA_H_INDUSTRY  )
IF (.NOT.LDATA_H_INDUSTRY) DEALLOCATE(XPAR_H_INDUSTRY)
!
CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','LE_INDUSTRY','TWN',CFNAM_LE_INDUSTRY,CFTYP_LE_INDUSTRY,XUNIF_LE_INDUSTRY,&
        XPAR_LE_INDUSTRY, LDATA_LE_INDUSTRY )
IF (.NOT.LDATA_LE_INDUSTRY) DEALLOCATE(XPAR_LE_INDUSTRY)
!
!-------------------------------------------------------------------------------
!
!* gardens
!
IF (OGARDEN) THEN
  CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','GARDEN     ','TWN',CFNAM_GARDEN    ,CFTYP_GARDEN    ,XUNIF_GARDEN    ,&
        XPAR_GARDEN, LDATA_GARDEN    )
  IF (.NOT.LDATA_GARDEN) DEALLOCATE(XPAR_GARDEN)
ELSE IF ( (XUNIF_GARDEN/=0. .AND. XUNIF_GARDEN/=XUNDEF) .OR. LEN_TRIM(CFNAM_GARDEN)/=0) THEN
  WRITE(ILUOUT,*) '---------------------------------------------'
  WRITE(ILUOUT,*) ' You chose not to include gardens in urban areas : LGARDEN=.FALSE.     '
  WRITE(ILUOUT,*) ' But            '
  IF (XUNIF_GARDEN/=0. .AND. XUNIF_GARDEN/=XUNDEF) THEN
    WRITE(ILUOUT,*) ' You also chose a garden fraction that is not zero : XUNIF_GARDEN=',XUNIF_GARDEN
  ELSE
    WRITE(ILUOUT,*) ' You also chose a garden fraction that is not zero : CFNAM_GARDEN=',CFNAM_GARDEN
  END IF
  WRITE(ILUOUT,*) '- - - - - - - - - - - - - - - - - - - - - - -'
  WRITE(ILUOUT,*) ' Please choose either:'
  WRITE(ILUOUT,*) ' LGARDEN=.TRUE. or set GARDEN fraction to zero (XUNIF_GARDEN=0.) in namelist PGD_TEB_PAR'
  WRITE(ILUOUT,*) '- - - - - - - - - - - - - - - - - - - - - - -'
  WRITE(ILUOUT,*) ' Beware that in this case, it may be necessary to change the'
  WRITE(ILUOUT,*) ' road fraction if you want to keep the same canyon aspect ratio'
  WRITE(ILUOUT,*) '---------------------------------------------'
  CALL ABOR1_SFX('PGD_TEB_PAR: GARDEN flag and GARDEN fraction not coherent')
END IF
!
!* vegetalized roofs (not yet implemented)
!
IF (.FALSE.) THEN
  CALL INI_VAR_FROM_DATA_0D(HPROGRAM,'ARI','VEG_ROOF   ','TWN',CFNAM_VEG_ROOF  ,CFTYP_VEG_ROOF  ,XUNIF_VEG_ROOF  ,XPAR_VEG_ROOF, &
        LDATA_VEG_ROOF  )
  IF (.NOT.LDATA_VEG_ROOF) DEALLOCATE(XPAR_VEG_ROOF)
ENDIF
!
!-------------------------------------------------------------------------------
!
!* coherence checks
!
!*    0.999 > building fraction > 0.001
!
IF (LDATA_BLD) THEN
  XPAR_BLD = MAX(XPAR_BLD,0.001)
  XPAR_BLD = MIN(XPAR_BLD,0.999)
END IF
!
!* gardens + building fraction < 0.999
!
IF (LDATA_BLD .AND. LDATA_GARDEN) THEN
  WHERE (XPAR_GARDEN+XPAR_BLD>=0.999)
    XPAR_GARDEN = 0.999 - XPAR_BLD
  END WHERE
END IF
!
!-------------------------------------------------------------------------------
IF (LHOOK)   CALL DR_HOOK('PGD_TEB_PAR',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PGD_TEB_PAR
