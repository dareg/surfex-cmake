!     #########
SUBROUTINE PREP_SNOW_EXTERN(HPROGRAM,HSURF,HFILE,HFILETYPE,KLUOUT,PFIELD,OSNOW_IDEAL,KLAYER)
!     #################################################################################
!
!
USE MODD_TYPE_SNOW
!
USE MODI_PREP_GRID_EXTERN
USE MODI_OPEN_AUX_IO_SURF
USE MODI_CLOSE_AUX_IO_SURF
USE MODE_SNOW3L
USE MODI_INTERP_GRID
USE MODI_READ_GR_SNOW
USE MODI_READ_SURF
USE MODI_SNOW_T_WLIQ_TO_HEAT
!
USE MODD_PREP,           ONLY : CINGRID_TYPE, CINTERP_TYPE
USE MODD_PREP_SNOW,      ONLY : XGRID_SNOW, NGRID_LEVEL
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_CSTS,           ONLY : XTT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_PUT_ON_ALL_VEGTYPES
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=10),  INTENT(IN)  :: HSURF     ! type of field
CHARACTER(LEN=28),  INTENT(IN)  :: HFILE     ! name of file
CHARACTER(LEN=6),   INTENT(IN)  :: HFILETYPE ! type of file
INTEGER,            INTENT(IN)  :: KLUOUT    ! logical unit of output listing
REAL,DIMENSION(:,:,:), POINTER  :: PFIELD    ! field to interpolate horizontally
LOGICAL,            INTENT(IN)  :: OSNOW_IDEAL
INTEGER,            INTENT(IN)  :: KLAYER
!
!*      0.2    declarations of local variables
!
TYPE(SURF_SNOW)                    :: TZSNOW ! snow characteristics

REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZFIELD       ! work field on input snow grid
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZFIELD_FINE  ! work field on fine snow grid
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTEMP        ! snow temperature
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWLIQ        ! liquid water snow pack content
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZD           ! total snow depth
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDEPTH       ! thickness of each layer (m)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZGRID        ! normalized input grid

CHARACTER(LEN=16)                 :: YRECFM         ! record name
INTEGER                           :: IRESP          ! error return code
INTEGER                           :: IVEGTYPE       ! actual number of vegtypes
INTEGER                           :: JLAYER         ! loop on snow vertical grids
INTEGER                           :: INI
CHARACTER(LEN=4)                  :: YAREA          ! area treated ('ROOF','ROAD','VEG ')
INTEGER                           :: IPATCH         ! number of input patch
INTEGER                           :: JPATCH         ! loop on patch
CHARACTER(LEN=6)                  :: YMASK          ! type of tile mask
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!*      3.     Area being treated
!              ------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_SNOW_EXTERN',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
YAREA = HSURF(7:10)
!
IF (YAREA=='VEG ') THEN
  IVEGTYPE = NVEGTYPE
  YMASK = 'NATURE'
ELSE
  IVEGTYPE = 1
  YMASK    = 'TOWN  '
  IPATCH   = 1
END IF
!
!*      1.     Preparation of IO for reading in the file
!              -----------------------------------------
!
!* Note that all points are read, even those without physical meaning.
!  These points will not be used during the horizontal interpolation step.
!  Their value must be defined as XUNDEF.
!
CALL OPEN_AUX_IO_SURF(HFILE,HFILETYPE,YMASK)
!
IF (YAREA=='VEG ') THEN
  YRECFM = 'PATCH_NUMBER'
  CALL READ_SURF(HFILETYPE,YRECFM,IPATCH,IRESP)
END IF
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading of grid
!              ---------------
!
CALL PREP_GRID_EXTERN(HFILETYPE,KLUOUT,CINGRID_TYPE,CINTERP_TYPE,INI)
!
!-------------------------------------------------------------------------------------
!
!*      4.     Reading of snow data
!              ---------------------
!
CALL READ_GR_SNOW(HFILETYPE,TRIM(YAREA),INI,IPATCH,TZSNOW,HDIR='A')
!
IF (TZSNOW%NLAYER.GT.KLAYER) THEN
  TZSNOW%NLAYER=KLAYER
ELSEIF (TZSNOW%NLAYER.LT.KLAYER) THEN
  CALL ABOR1_SFX("PREP_SNOW_EXTERN: SNOW NLAYER IN EXTERN FILE MUST BE THE SAME AS CURRENT NLAYER")
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      5.     Total snow content
!              ------------------
!
SELECT CASE (HSURF(1:3))
  CASE ('WWW')
    IF (OSNOW_IDEAL) THEN
      ALLOCATE(ZFIELD(INI,TZSNOW%NLAYER,IPATCH))
      ZFIELD(:,:,:) = TZSNOW%WSNOW(:,1:TZSNOW%NLAYER,:)
      ALLOCATE(PFIELD(INI,TZSNOW%NLAYER,IVEGTYPE))
      CALL PUT_ON_ALL_VEGTYPES(INI,TZSNOW%NLAYER,IPATCH,IVEGTYPE,ZFIELD,PFIELD)
    ELSE
      ALLOCATE(ZFIELD(INI,1,IPATCH))
      ZFIELD = 0.
      DO JLAYER=1,TZSNOW%NLAYER
        ZFIELD(:,1,:) = ZFIELD(:,1,:) + TZSNOW%WSNOW(:,JLAYER,:)
      END DO
      ALLOCATE(PFIELD(INI,1,IVEGTYPE))
      CALL PUT_ON_ALL_VEGTYPES(INI,1,IPATCH,IVEGTYPE,ZFIELD,PFIELD)
    ENDIF
    DEALLOCATE(ZFIELD)
!
!-------------------------------------------------------------------------------------
!
!*      6.     Albedo
!              ------
!
  CASE ('ALB')
    ALLOCATE(ZFIELD(INI,1,IPATCH))
    ZFIELD(:,1,:) = TZSNOW%ALB(:,:)
    !
    ALLOCATE(PFIELD(INI,1,IVEGTYPE))
    CALL PUT_ON_ALL_VEGTYPES(INI,1,IPATCH,IVEGTYPE,ZFIELD,PFIELD)
    DEALLOCATE(ZFIELD)
!
!-------------------------------------------------------------------------------------
!
!*      7.     Total depth
!              -----------
!
  CASE ('DEP')
    IF (OSNOW_IDEAL) THEN
      ALLOCATE(ZDEPTH(INI,TZSNOW%NLAYER,IPATCH))
      ZDEPTH(:,:,:) = TZSNOW%WSNOW(:,1:TZSNOW%NLAYER,:)/TZSNOW%RHO(:,1:TZSNOW%NLAYER,:)
      WHERE(TZSNOW%WSNOW(:,1:TZSNOW%NLAYER,:)==XUNDEF) ZDEPTH(:,:,:)=XUNDEF
      ALLOCATE(PFIELD(INI,TZSNOW%NLAYER,IVEGTYPE))
      CALL PUT_ON_ALL_VEGTYPES(INI,TZSNOW%NLAYER,IPATCH,IVEGTYPE,ZDEPTH,PFIELD)
    ELSE
      ALLOCATE(ZDEPTH(INI,1,IPATCH))
      ZDEPTH = 0.
      DO JPATCH=1,IPATCH
        DO JLAYER=1,TZSNOW%NLAYER
          ZDEPTH(:,1,JPATCH) = ZDEPTH(:,1,JPATCH) + &
                               TZSNOW%WSNOW(:,JLAYER,JPATCH)/TZSNOW%RHO(:,JLAYER,JPATCH)
        END DO
        WHERE(TZSNOW%WSNOW(:,1,JPATCH)==XUNDEF) ZDEPTH(:,1,JPATCH)=XUNDEF
      END DO
      ALLOCATE(PFIELD(INI,1,IVEGTYPE))
      CALL PUT_ON_ALL_VEGTYPES(INI,1,IPATCH,IVEGTYPE,ZDEPTH,PFIELD)
    ENDIF
    DEALLOCATE(ZDEPTH)
!
!-------------------------------------------------------------------------------------
!
!*      8.     Density or heat profile
!              -----------------------
!
  CASE ('RHO','HEA')
    ALLOCATE(ZFIELD     (INI,TZSNOW%NLAYER,IPATCH))

    SELECT CASE (TZSNOW%SCHEME)
      CASE ('D95','1-L','EBA')
        ALLOCATE(ZFIELD_FINE(INI,NGRID_LEVEL,IPATCH))      
        !* computes output physical variable
        IF (HSURF(1:1)=='R') ZFIELD(:,1,:) = TZSNOW%RHO(:,1,:)
        IF (HSURF(1:1)=='H') THEN
          ALLOCATE(ZTEMP(INI,TZSNOW%NLAYER,IPATCH))
          ALLOCATE(ZWLIQ(INI,TZSNOW%NLAYER,IPATCH))
          IF (TZSNOW%SCHEME=='D95'.OR.TZSNOW%SCHEME=='EBA') ZTEMP(:,1,:) = XTT
          IF (TZSNOW%SCHEME=='1-L') ZTEMP(:,1,:) = TZSNOW%T(:,1,:)
          ZWLIQ(:,:,:) = 0.
          CALL SNOW_T_WLIQ_TO_HEAT(ZFIELD,TZSNOW%RHO,ZTEMP,ZWLIQ)
          DEALLOCATE(ZTEMP)
          DEALLOCATE(ZWLIQ)
        END IF
        !* put profile on fine snow grid
        DO JLAYER=1,NGRID_LEVEL
          ZFIELD_FINE(:,JLAYER,:) = ZFIELD(:,1,:)
        END DO
        ALLOCATE(PFIELD(INI,NGRID_LEVEL,IVEGTYPE))
        CALL PUT_ON_ALL_VEGTYPES(INI,NGRID_LEVEL,IPATCH,IVEGTYPE,ZFIELD_FINE,PFIELD)

      CASE ('3-L','CRO')
        !* input physical variable
        IF (HSURF(1:1)=='R') ZFIELD(:,:,:) = TZSNOW%RHO (:,1:TZSNOW%NLAYER,:)
        IF (HSURF(1:1)=='H') ZFIELD(:,:,:) = TZSNOW%HEAT(:,1:TZSNOW%NLAYER,:)
        !
        IF (OSNOW_IDEAL) THEN
          ALLOCATE(ZFIELD_FINE(INI,TZSNOW%NLAYER,IPATCH))
          IF (HSURF(1:1)=='R') ZFIELD_FINE(:,:,:) = ZFIELD(:,:,:)
          IF (HSURF(1:1)=='H') ZFIELD_FINE(:,:,:) = ZFIELD(:,:,:)
          ALLOCATE(PFIELD(INI,TZSNOW%NLAYER,IVEGTYPE))
          CALL PUT_ON_ALL_VEGTYPES(INI,TZSNOW%NLAYER,IPATCH,IVEGTYPE,ZFIELD_FINE,PFIELD)
        ELSE
          ALLOCATE(ZFIELD_FINE(INI,NGRID_LEVEL,IPATCH))      
          !* total depth
          ALLOCATE(ZD(INI,IPATCH))
          ZD = 0.
          DO JLAYER=1,TZSNOW%NLAYER
            ZD(:,:) = ZD(:,:) + TZSNOW%WSNOW(:,JLAYER,:)/TZSNOW%RHO(:,JLAYER,:)
          END DO
          !* input snow layer thickness
          ALLOCATE(ZDEPTH(INI,TZSNOW%NLAYER,IPATCH))
          DO JPATCH=1,IPATCH
            CALL SNOW3LGRID(ZDEPTH(:,:,JPATCH),ZD(:,JPATCH))
          END DO
          !* input normalized grid
          ALLOCATE(ZGRID(INI,TZSNOW%NLAYER,IPATCH))
          WHERE(ZD(:,:)>0.)
            ZGRID(:,1,:) = 0.5 * ZDEPTH(:,1,:) / ZD(:,:)
          ELSEWHERE
            ZGRID(:,1,:) = 0.5 / FLOAT(TZSNOW%NLAYER)
          END WHERE
          DO JLAYER = 2,TZSNOW%NLAYER
            WHERE (ZD(:,:)>0.)
              ZGRID(:,JLAYER,:) = ZGRID(:,JLAYER-1,:)                     &
                                 + 0.5 * ZDEPTH(:,JLAYER-1,:) / ZD(:,:)   &
                                 + 0.5 * ZDEPTH(:,JLAYER  ,:) / ZD(:,:)
            ELSEWHERE
              ZGRID(:,JLAYER,:) = (FLOAT(JLAYER)-0.5) / FLOAT(TZSNOW%NLAYER)
            END WHERE
          END DO
          DEALLOCATE(ZD)
          DEALLOCATE(ZDEPTH)
          !* interpolation of profile onto fine normalized snow grid
          DO JPATCH=1,IPATCH
            CALL INTERP_GRID(ZGRID(:,:,JPATCH),ZFIELD(:,:,JPATCH),    &
                             XGRID_SNOW(:),    ZFIELD_FINE(:,:,JPATCH))
          END DO
          DEALLOCATE(ZGRID)
          ALLOCATE(PFIELD(INI,NGRID_LEVEL,IVEGTYPE))
          CALL PUT_ON_ALL_VEGTYPES(INI,NGRID_LEVEL,IPATCH,IVEGTYPE,ZFIELD_FINE,PFIELD)
        ENDIF
      END SELECT
    !
    !* put field form patch to all vegtypes
    DEALLOCATE(ZFIELD)
    DEALLOCATE(ZFIELD_FINE)
!
!*      9.     Crocus specific fields
!              -----------------------
!
  CASE ('SG1','SG2','HIS','AGE')
    ALLOCATE(ZFIELD(INI,TZSNOW%NLAYER,IPATCH))
    IF (HSURF(1:3)=='SG1') ZFIELD(:,:,:) = TZSNOW%GRAN1(:,1:TZSNOW%NLAYER,:)
    IF (HSURF(1:3)=='SG2') ZFIELD(:,:,:) = TZSNOW%GRAN2(:,1:TZSNOW%NLAYER,:)
    IF (HSURF(1:3)=='HIS') ZFIELD(:,:,:) = TZSNOW%HIST(:,1:TZSNOW%NLAYER,:)
    IF (HSURF(1:3)=='AGE') ZFIELD(:,:,:) = TZSNOW%AGE(:,1:TZSNOW%NLAYER,:)
    ALLOCATE(PFIELD(INI,TZSNOW%NLAYER,IVEGTYPE))
    CALL PUT_ON_ALL_VEGTYPES(INI,TZSNOW%NLAYER,IPATCH,IVEGTYPE,ZFIELD,PFIELD)
    DEALLOCATE(ZFIELD)
  !    
END SELECT
!
!-------------------------------------------------------------------------------------
!
!*      9.     End of IO
!              ---------
!
CALL CLOSE_AUX_IO_SURF(HFILE,HFILETYPE)
IF (LHOOK) CALL DR_HOOK('PREP_SNOW_EXTERN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
END SUBROUTINE PREP_SNOW_EXTERN
