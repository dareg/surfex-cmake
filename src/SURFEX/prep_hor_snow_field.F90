!     #########
SUBROUTINE PREP_HOR_SNOW_FIELD (DTCO, &
                                 IG, U, &
                                HPROGRAM,                       &
                                HFILE,HFILETYPE,                &
                                HFILEPGD,HFILEPGDTYPE,          &
                                KLUOUT,OUNIF,HSNSURF,KPATCH,    &
                                KTEB_PATCH, &
                                KL,TPSNOW, TPTIME,              &
                                PUNIF_WSNOW, PUNIF_RSNOW,       &
                                PUNIF_TSNOW, PUNIF_LWCSNOW,     &
                                PUNIF_ASNOW, OSNOW_IDEAL,       &
                                PUNIF_SG1SNOW, PUNIF_SG2SNOW,   &
                                PUNIF_HISTSNOW,PUNIF_AGESNOW,   &                                
                                PF,PDEPTH,PVEGTYPE,             &
                                PVEGTYPE_PATCH,PPATCH           )
!     #######################################################
!
!!****  *PREP_HOR_SNOW_FIELD* - reads, interpolates and prepares a snow field
!!
!!    PURPOSE
!!    -------
!!
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
!!      P. Le Moigne 10/2005, Phasage Arome
!!      B. Decharme  10/2013, Phasage Arpège-Climat
!!      M. Lafaysse 11/2012, snow liquid water content
!!      B. Decharme  04/2014, external init with FA files
!!                            new init for ES
!!------------------------------------------------------------------
!
!
!
USE MODD_GRID_GRIB, ONLY : CINMODEL
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
!
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_TYPE_SNOW
USE MODD_TYPE_DATE_SURF, ONLY : DATE_TIME
!
USE MODD_CSTS,           ONLY : XTT
USE MODD_PREP_SNOW,      ONLY : XGRID_SNOW
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE, NVT_SNOW
USE MODD_PREP,           ONLY : LINTERP,CINTERP_TYPE,CINGRID_TYPE
!
USE MODD_SNOW_PAR, ONLY : XANSMAX
!
USE MODI_PREP_GRIB_GRID
USE MODI_PREP_SNOW_GRIB
USE MODI_PREP_SNOW_UNIF
USE MODI_PREP_SNOW_EXTERN
USE MODI_PREP_SNOW_BUFFER
USE MODI_HOR_INTERPOL
USE MODI_VEGTYPE_GRID_TO_PATCH_GRID
USE MODI_SNOW_HEAT_TO_T_WLIQ
USE MODI_VEGTYPE_TO_PATCH
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*      0.1    declarations of arguments
!
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
!
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=28),  INTENT(IN)  :: HFILE     ! file name
 CHARACTER(LEN=6),   INTENT(IN)  :: HFILETYPE ! file type
 CHARACTER(LEN=28),  INTENT(IN)  :: HFILEPGD     ! file name
 CHARACTER(LEN=6),   INTENT(IN)  :: HFILEPGDTYPE ! file type
INTEGER,            INTENT(IN)  :: KLUOUT    ! logical unit of output listing
LOGICAL,            INTENT(IN)  :: OUNIF     ! flag for prescribed uniform field
 CHARACTER(LEN=10)               :: HSNSURF   ! type of field
INTEGER,            INTENT(IN)  :: KPATCH    ! patch number for output scheme
INTEGER,            INTENT(IN) :: KTEB_PATCH
INTEGER,            INTENT(IN)  :: KL        ! number of points
TYPE(SURF_SNOW)                 :: TPSNOW    ! snow fields
TYPE(DATE_TIME),    INTENT(IN)  :: TPTIME    ! date and time
REAL, DIMENSION(:), INTENT(IN)  :: PUNIF_WSNOW ! prescribed snow content (kg/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PUNIF_RSNOW ! prescribed density (kg/m3)
REAL, DIMENSION(:), INTENT(IN)  :: PUNIF_TSNOW ! prescribed temperature (K)
REAL, DIMENSION(:), INTENT(IN)  :: PUNIF_LWCSNOW ! prescribed snow liquid water content (kg/m3)
REAL,               INTENT(IN)  :: PUNIF_ASNOW ! prescribed albedo (-)
LOGICAL,            INTENT(IN)  :: OSNOW_IDEAL
REAL, DIMENSION(:), INTENT(IN)  :: PUNIF_SG1SNOW ! 
REAL, DIMENSION(:), INTENT(IN)  :: PUNIF_SG2SNOW ! 
REAL, DIMENSION(:), INTENT(IN)  :: PUNIF_HISTSNOW ! 
REAL, DIMENSION(:), INTENT(IN)  :: PUNIF_AGESNOW ! 

REAL,DIMENSION(:,:,:),  INTENT(OUT),OPTIONAL :: PF     ! output field (x,kpatch)
REAL,DIMENSION(:,:,:),INTENT(IN), OPTIONAL :: PDEPTH ! thickness of each snow layer
REAL,DIMENSION(:,:),  INTENT(IN), OPTIONAL :: PVEGTYPE ! fraction of each vegtype
REAL,DIMENSION(:,:,:),  INTENT(IN), OPTIONAL :: PVEGTYPE_PATCH ! fraction of each vegtype per patch
REAL,DIMENSION(:,:),  INTENT(IN), OPTIONAL :: PPATCH ! fraction of each patch
!
!
!*      0.2    declarations of local variables
!
REAL, POINTER, DIMENSION(:,:,:)     :: ZFIELDIN=>NULL()  ! field to interpolate horizontally
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZFIELDOUTP ! field interpolated   horizontally
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZFIELDOUTV !
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZD        ! snow depth (x, kpatch)
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZW        ! work array (x, fine   snow grid, kpatch)
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZHEAT     ! work array (x, output snow grid, kpatch)
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZGRID     ! grid array (x, output snow grid, kpatch)
!
TYPE (DATE_TIME)              :: TZTIME_GRIB    ! current date and time
LOGICAL                       :: GSNOW_IDEAL
INTEGER                       :: JPATCH    ! loop on patches
INTEGER                       :: JVEGTYPE  ! loop on vegtypes
INTEGER                       :: JLAYER    ! loop on layers
INTEGER :: INFOMPI, INL, INP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
!*      1.     Does the field exist?
!
!
IF (LHOOK) CALL DR_HOOK('PREP_HOR_SNOW_FIELD',0,ZHOOK_HANDLE)
!
GSNOW_IDEAL = .FALSE.
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      2.     Reading of input  configuration (Grid and interpolation type)
!
IF (OUNIF) THEN
  GSNOW_IDEAL = OSNOW_IDEAL
  CALL PREP_SNOW_UNIF(KLUOUT,HSNSURF,ZFIELDIN, TPTIME, GSNOW_IDEAL,       &
                      PUNIF_WSNOW, PUNIF_RSNOW, PUNIF_TSNOW,              &
                      PUNIF_LWCSNOW, PUNIF_ASNOW, PUNIF_SG1SNOW,          &
                      PUNIF_SG2SNOW, PUNIF_HISTSNOW, PUNIF_AGESNOW,       &
                      TPSNOW%NLAYER                                       )
ELSE IF (HFILETYPE=='GRIB  ') THEN
  CALL PREP_GRIB_GRID(HFILE,KLUOUT,CINMODEL,CINGRID_TYPE,CINTERP_TYPE,TZTIME_GRIB)            
   IF (NRANK==NPIO) CALL PREP_SNOW_GRIB(HPROGRAM,HSNSURF,HFILE,KLUOUT,TPSNOW%NLAYER,ZFIELDIN)        
ELSE IF (HFILETYPE=='MESONH' .OR. HFILETYPE=='ASCII ' .OR. HFILETYPE=='LFI   '.OR. HFILETYPE=='FA    ') THEN
  GSNOW_IDEAL = OSNOW_IDEAL
  CALL PREP_SNOW_EXTERN(&
                        HPROGRAM,HSNSURF,HFILE,HFILETYPE,HFILEPGD,HFILEPGDTYPE,&
                        KLUOUT,ZFIELDIN,GSNOW_IDEAL,TPSNOW%NLAYER,KTEB_PATCH)
ELSE IF (HFILETYPE=='BUFFER') THEN
  CALL PREP_SNOW_BUFFER(IG, U, &
                        HPROGRAM,HSNSURF,KLUOUT,TPSNOW%NLAYER,ZFIELDIN)
ELSE
  CALL ABOR1_SFX('PREP_HOR_SNOW_FIELD: data file type not supported : '//HFILETYPE)
END IF
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      3.     Horizontal interpolation
!
IF (NRANK==NPIO) THEN
  INL = SIZE(ZFIELDIN,2)
  INP = SIZE(ZFIELDIN,3)
ELSEIF (.NOT.ASSOCIATED(ZFIELDIN)) THEN
 ALLOCATE(ZFIELDIN(0,0,0))
ENDIF
!
IF (NPROC>1) THEN
#ifdef SFX_MPI
  CALL MPI_BCAST(INL,KIND(INL)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
  CALL MPI_BCAST(INP,KIND(INP)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
#endif
ENDIF
!    
ALLOCATE(ZFIELDOUTP(KL,INL,INP))
!
DO JVEGTYPE = 1, INP
  IF (INP==NVEGTYPE) THEN
    JPATCH = 1
    IF (KPATCH>1) JPATCH = VEGTYPE_TO_PATCH(JVEGTYPE,KPATCH)
    !* does not interpolates snow caracteristics on points without snow
    IF (PRESENT(PDEPTH)) LINTERP(:) = ( PDEPTH(:,1,JPATCH) /= 0. .AND. PDEPTH(:,1,JPATCH) /= XUNDEF )
    IF (PRESENT(PPATCH)) LINTERP(:) = (LINTERP(:) .AND. PPATCH(:,JPATCH)>0.)
  ELSE
    IF (PRESENT(PDEPTH)) THEN
      LINTERP(:) = .FALSE.
      DO JPATCH = 1,KPATCH
        LINTERP(:) = LINTERP(:) .OR. (PDEPTH(:,1,JPATCH)/=0. .AND. PDEPTH(:,1,JPATCH)/=XUNDEF)
      ENDDO
    ENDIF
  ENDIF
  !* horizontal interpolation
  CALL HOR_INTERPOL(DTCO, U, &
                    KLUOUT,ZFIELDIN(:,:,JVEGTYPE),ZFIELDOUTP(:,:,JVEGTYPE))
  !
  LINTERP(:) = .TRUE.
END DO
!
DEALLOCATE(ZFIELDIN )
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      4.     Transformation from vegtype grid to patch grid, if any
!
ALLOCATE(ZW (KL,INL,KPATCH))
!
ZW = 0.
IF (KPATCH/=INP) THEN
  !
  ALLOCATE(ZFIELDOUTV(KL,INL,NVEGTYPE))
  CALL PUT_ON_ALL_VEGTYPES(KL,INL,INP,NVEGTYPE,ZFIELDOUTP,ZFIELDOUTV)
  !
  !*      6.     Transformation from vegtype grid to patch grid
  !
  CALL VEGTYPE_GRID_TO_PATCH_GRID(KPATCH,PVEGTYPE_PATCH,PPATCH,ZFIELDOUTV,ZW)
  DEALLOCATE(ZFIELDOUTV)
  !
ELSE
  !
  ZW(:,:,:) = ZFIELDOUTP(:,:,:)
  !
ENDIF
!
DEALLOCATE(ZFIELDOUTP)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      5.     Defines normalized output grid, if depths of snow layers are present
!
IF (PRESENT(PDEPTH) .AND. .NOT.GSNOW_IDEAL) THEN
!
!* total snow depth
!
  ALLOCATE(ZD(SIZE(TPSNOW%WSNOW,1),KPATCH))
  ZD(:,:)=0.
  DO JPATCH=1,KPATCH
    DO JLAYER=1,TPSNOW%NLAYER
      WHERE (PDEPTH(:,JLAYER,JPATCH)/=XUNDEF) ZD(:,JPATCH) = ZD(:,JPATCH) + PDEPTH(:,JLAYER,JPATCH)
    END DO
  END DO
!
!* grid at center of layers
!
  ALLOCATE(ZGRID(SIZE(ZW,1),TPSNOW%NLAYER,KPATCH))
  ZGRID(:,1,:) = PDEPTH(:,1,:)
  IF(TPSNOW%NLAYER>1)THEN
    DO JPATCH=1,KPATCH
      DO JLAYER=2,TPSNOW%NLAYER
        ZGRID(:,JLAYER,JPATCH) = ZGRID(:,JLAYER-1,JPATCH) + PDEPTH(:,JLAYER,JPATCH)
      ENDDO
    ENDDO
  ENDIF
!
!* normalized grid
!
  DO JPATCH=1,KPATCH
    DO JLAYER=1,TPSNOW%NLAYER
      WHERE (ZD(:,JPATCH)/=0.)
        ZGRID(:,JLAYER,JPATCH) = ZGRID(:,JLAYER,JPATCH) / ZD(:,JPATCH)
      ELSEWHERE
        ZGRID(:,JLAYER,JPATCH) = 1.0
      END WHERE
    END DO
  END DO
!
  DEALLOCATE(ZD)
!
ELSEIF (.NOT.GSNOW_IDEAL) THEN
  IF (HSNSURF(1:3)=='RHO' .OR. HSNSURF(1:3)=='HEA') THEN
    WRITE(KLUOUT,*) 'when interpolation profiles of snow pack quantities,'
    WRITE(KLUOUT,*) 'depth of snow layers must be given'
    CALL ABOR1_SFX('PREP_HOR_SNOW_FIELD: DEPTH OF SNOW LAYERS NEEDED')
  END IF
END IF
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      6.     Return to historical variable
!
SELECT CASE (HSNSURF(1:3))
  !
  CASE('WWW')  ! total snow content
    !
    IF (GSNOW_IDEAL) THEN
      PF(:,:,:) = ZW(:,:,:)
    ELSE
      DO JLAYER=1,SIZE(PF,2)
        PF(:,JLAYER,:) = ZW(:,1,:)
      ENDDO
    ENDIF
    !
    IF (PRESENT(PPATCH)) THEN
      DO JLAYER = 1,TPSNOW%NLAYER
        WHERE(PPATCH(:,:)==0.)
          PF(:,JLAYER,:) = XUNDEF
        END WHERE
      ENDDO
    ENDIF
  !
  CASE('DEP')  ! snow thickness
    !
    IF (GSNOW_IDEAL) THEN
      PF(:,:,:) = ZW(:,:,:)
    ELSE
      DO JLAYER=1,SIZE(PF,2)
        PF(:,JLAYER,:) = ZW(:,1,:)
      ENDDO
    ENDIF
    !
    IF (PRESENT(PPATCH)) THEN
      DO JLAYER = 1,TPSNOW%NLAYER
        WHERE(PPATCH(:,:)==0.)
          PF(:,JLAYER,:) = XUNDEF
        END WHERE
      ENDDO
    ENDIF
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  CASE('RHO') 
    !
    IF (GSNOW_IDEAL) THEN
      TPSNOW%RHO(:,:,:) = ZW(:,:,:)
    ELSEIF(SIZE(ZW,2)==1) THEN
      DO JLAYER = 1,TPSNOW%NLAYER
        TPSNOW%RHO(:,JLAYER,:) = ZW(:,1,:)
      ENDDO      
    ELSEIF(SIZE(ZW,2)==TPSNOW%NLAYER)THEN
      TPSNOW%RHO(:,:,:) = ZW(:,:,:)
    ELSE
      !* interpolation on snow levels
      CALL INIT_FROM_REF_GRID(XGRID_SNOW,ZW,ZGRID,TPSNOW%RHO)
    ENDIF
    !
    !* mask for areas where there is no snow
    DO JPATCH=1,KPATCH
      DO JLAYER=1,TPSNOW%NLAYER
        WHERE(PDEPTH(:,JLAYER,JPATCH)==0. .OR. PDEPTH(:,JLAYER,JPATCH)==XUNDEF) TPSNOW%RHO(:,JLAYER,JPATCH) = XUNDEF
      END DO
    END DO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  CASE('ALB')
    !
    DO JPATCH=1,KPATCH
      TPSNOW%ALB(:,JPATCH) = ZW(:,1,JPATCH)
    END DO
    !
    !* mask for areas where there is no snow
    DO JPATCH=1,KPATCH
      WHERE(PDEPTH(:,1,JPATCH)==0. .OR. PDEPTH(:,1,JPATCH)==XUNDEF)  TPSNOW%ALB(:,JPATCH) = XUNDEF
    END DO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  CASE('HEA') 
    !
    IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
      !
      IF (GSNOW_IDEAL) THEN
        TPSNOW%HEAT(:,:,:) = ZW(:,:,:)
      ELSEIF(SIZE(ZW,2)==1) THEN
        DO JLAYER = 1,TPSNOW%NLAYER
          TPSNOW%HEAT(:,JLAYER,:) = ZW(:,1,:)
        ENDDO        
      ELSEIF(SIZE(ZW,2)==TPSNOW%NLAYER)THEN
        TPSNOW%HEAT(:,:,:) = ZW(:,:,:)
      ELSE
        !* interpolation of heat on snow levels
        CALL INIT_FROM_REF_GRID(XGRID_SNOW,ZW,ZGRID,TPSNOW%HEAT)
      ENDIF
      !
      !* mask for areas where there is no snow
      DO JPATCH=1,KPATCH
        DO JLAYER=1,TPSNOW%NLAYER
          WHERE(PDEPTH(:,JLAYER,JPATCH)==0. .OR. PDEPTH(:,JLAYER,JPATCH)==XUNDEF) TPSNOW%HEAT(:,JLAYER,JPATCH) = XUNDEF
        END DO
      END DO
      !
    ELSE IF (TPSNOW%SCHEME=='1-L') THEN
      !* interpolation of heat on snow levels
      ALLOCATE(ZHEAT(KL,TPSNOW%NLAYER,KPATCH))
      !
      IF (GSNOW_IDEAL) THEN
        ZHEAT(:,:,:) = ZW(:,:,:)
      ELSEIF(SIZE(ZW,2)==1) THEN
        DO JLAYER = 1,TPSNOW%NLAYER
          ZHEAT(:,JLAYER,:) = ZW(:,1,:)
        ENDDO        
      ELSEIF(SIZE(ZW,2)==TPSNOW%NLAYER)THEN
        ZHEAT(:,:,:) = ZW(:,:,:)
      ELSE
        CALL INIT_FROM_REF_GRID(XGRID_SNOW,ZW,ZGRID,ZHEAT)
      ENDIF
      !
      !* transformation from heat to temperature
      CALL SNOW_HEAT_TO_T_WLIQ(ZHEAT,TPSNOW%RHO,TPSNOW%T)
      WHERE (TPSNOW%T>XTT) TPSNOW%T = XTT
      DEALLOCATE(ZHEAT)
      !
      !* mask for areas where there is no snow
      DO JPATCH=1,KPATCH
        DO JLAYER=1,TPSNOW%NLAYER
          WHERE(PDEPTH(:,JLAYER,JPATCH)==0. .OR. PDEPTH(:,JLAYER,JPATCH)==XUNDEF) TPSNOW%T(:,JLAYER,JPATCH) = XUNDEF
        END DO
      END DO
      !
    END IF
  !
  !
  CASE('SG1')
    !
    IF (GSNOW_IDEAL) THEN
      TPSNOW%GRAN1(:,:,:) = ZW(:,:,:)
    ELSEIF(SIZE(ZW,2)==1) THEN
      DO JLAYER = 1,TPSNOW%NLAYER
        TPSNOW%GRAN1(:,JLAYER,:) = ZW(:,1,:)
      ENDDO      
    ELSEIF(SIZE(ZW,2)==TPSNOW%NLAYER)THEN
      TPSNOW%GRAN1(:,:,:) = ZW(:,:,:)
    ELSE
      !* interpolation of heat on snow levels
      CALL INIT_FROM_REF_GRID(XGRID_SNOW,ZW,ZGRID,TPSNOW%GRAN1)
    ENDIF
    !
    !* mask for areas where there is no snow
    DO JPATCH=1,KPATCH
      DO JLAYER=1,TPSNOW%NLAYER
        WHERE(PDEPTH(:,JLAYER,JPATCH)==0. .OR. PDEPTH(:,JLAYER,JPATCH)==XUNDEF) TPSNOW%GRAN1(:,JLAYER,JPATCH) = XUNDEF
      END DO
    END DO
    !
  CASE('SG2')
    !
    IF (GSNOW_IDEAL) THEN
      TPSNOW%GRAN2(:,:,:) = ZW(:,:,:)
    ELSEIF(SIZE(ZW,2)==1) THEN
      DO JLAYER = 1,TPSNOW%NLAYER
        TPSNOW%GRAN2(:,JLAYER,:) = ZW(:,1,:)
      ENDDO      
    ELSEIF(SIZE(ZW,2)==TPSNOW%NLAYER)THEN
      TPSNOW%GRAN2(:,:,:) = ZW(:,:,:)
    ELSE
      !* interpolation of heat on snow levels
      CALL INIT_FROM_REF_GRID(XGRID_SNOW,ZW,ZGRID,TPSNOW%GRAN2)
    ENDIF
    !
    !* mask for areas where there is no snow
    DO JPATCH=1,KPATCH
      DO JLAYER=1,TPSNOW%NLAYER
        WHERE(PDEPTH(:,JLAYER,JPATCH)==0. .OR. PDEPTH(:,JLAYER,JPATCH)==XUNDEF) TPSNOW%GRAN2(:,JLAYER,JPATCH) = XUNDEF
      END DO
    END DO
    !
  CASE('HIS')
    !
    IF (GSNOW_IDEAL) THEN
      TPSNOW%HIST(:,:,:) = ZW(:,:,:)
    ELSEIF(SIZE(ZW,2)==1) THEN
      DO JLAYER = 1,TPSNOW%NLAYER
        TPSNOW%HIST(:,JLAYER,:) = ZW(:,1,:)
      ENDDO      
    ELSEIF(SIZE(ZW,2)==TPSNOW%NLAYER)THEN
      TPSNOW%HIST(:,:,:) = ZW(:,:,:)
    ELSE
      !* interpolation of heat on snow levels
      CALL INIT_FROM_REF_GRID(XGRID_SNOW,ZW,ZGRID,TPSNOW%HIST)
    ENDIF
    !
    !* mask for areas where there is no snow
    DO JPATCH=1,KPATCH
      DO JLAYER=1,TPSNOW%NLAYER
        WHERE(PDEPTH(:,JLAYER,JPATCH)==0. .OR. PDEPTH(:,JLAYER,JPATCH)==XUNDEF) TPSNOW%HIST(:,JLAYER,JPATCH) = XUNDEF
      END DO
    END DO
    !
  CASE('AGE')
    !
    IF (TPSNOW%SCHEME=='3-L'.AND.(.NOT.GSNOW_IDEAL).AND.(.NOT.OUNIF))THEN
       TPSNOW%AGE(:,:,:) = 0.0
    ELSE
      IF (GSNOW_IDEAL) THEN
        TPSNOW%AGE(:,:,:) = ZW(:,:,:)
      ELSEIF(SIZE(ZW,2)==1) THEN
        DO JLAYER = 1,TPSNOW%NLAYER
          TPSNOW%AGE(:,JLAYER,:) = ZW(:,1,:)
        ENDDO        
      ELSEIF(SIZE(ZW,2)==TPSNOW%NLAYER)THEN
        TPSNOW%AGE(:,:,:) = ZW(:,:,:)
      ELSE
        !* interpolation of heat on snow levels
        CALL INIT_FROM_REF_GRID(XGRID_SNOW,ZW,ZGRID,TPSNOW%AGE)
      ENDIF
    ENDIF
    !
    !* mask for areas where there is no snow
    DO JPATCH=1,KPATCH
      DO JLAYER=1,TPSNOW%NLAYER
        WHERE(PDEPTH(:,JLAYER,JPATCH)==0. .OR. PDEPTH(:,JLAYER,JPATCH)==XUNDEF) TPSNOW%AGE(:,JLAYER,JPATCH) = XUNDEF
      END DO
    END DO
    !
END SELECT
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      7.     Deallocations
!
IF (PRESENT(PDEPTH) .AND. .NOT.GSNOW_IDEAL) DEALLOCATE(ZGRID    )
DEALLOCATE(ZW       )
IF (LHOOK) CALL DR_HOOK('PREP_HOR_SNOW_FIELD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------------
!
SUBROUTINE INIT_FROM_REF_GRID(PGRID1,PT1,PD2,PT2)
!
USE MODI_INTERP_GRID_NAT
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PT1    ! variable profile
REAL, DIMENSION(:),     INTENT(IN)  :: PGRID1 ! normalized grid
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PD2    ! output layer thickness
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PT2    ! variable profile
!
INTEGER                                  :: JL, JL1  ! loop counter
REAL, DIMENSION(SIZE(PT1,1),SIZE(PGRID1),SIZE(PT1,3)) :: ZT1
REAL, DIMENSION(SIZE(PT1,1),SIZE(PGRID1)) :: ZD1 ! input grid
REAL, DIMENSION(SIZE(PD2,1),SIZE(PD2,2)) :: ZD2 ! output grid
INTEGER                       :: JPATCH    ! loop on patches
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (LHOOK) CALL DR_HOOK('INIT_FROM_REF_GRID',0,ZHOOK_HANDLE)
DO JPATCH=1,KPATCH
  ZD2(:,:) = 0.
  DO JL=1,SIZE(ZD2,2)
    ZD2(:,JL) = PD2(:,JL,JPATCH)
  END DO
  !
  DO JL=1,SIZE(PGRID1)
    JL1 = MIN(JL,SIZE(PT1,2))
    ZT1(:,JL,:) = PT1(:,JL1,:)
    ZD1(:,JL) = PGRID1(JL)
  END DO
  !
  CALL INTERP_GRID_NAT(ZD1,ZT1(:,:,JPATCH),ZD2,PT2(:,:,JPATCH))
END DO
IF (LHOOK) CALL DR_HOOK('INIT_FROM_REF_GRID',1,ZHOOK_HANDLE)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
END SUBROUTINE INIT_FROM_REF_GRID
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_HOR_SNOW_FIELD
