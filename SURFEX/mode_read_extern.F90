!     #####################
MODULE MODE_READ_EXTERN
!     #####################
!-------------------------------------------------------------------
!
USE MODI_READ_LECOCLIMAP
!
USE MODI_PUT_ON_ALL_VEGTYPES
USE MODI_OLD_NAME
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
CONTAINS
!
!---------------------------------------------------------------------------------------
!
!     #######################
      SUBROUTINE READ_EXTERN_DEPTH(HPROGRAM,KLUOUT,HISBA,HNAT,KNI,KLAYER,KPATCH,PDEPTH)
!     #######################
!
USE MODD_SURF_PAR,       ONLY : NUNDEF, XUNDEF
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER, NVEGTYPE
!
USE MODI_READ_SURF
USE MODI_CONVERT_COVER_ISBA
USE MODI_GARDEN_SOIL_DEPTH
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
!
CHARACTER(LEN=6),     INTENT(IN)    :: HPROGRAM  ! program calling surf. schemes
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
CHARACTER(LEN=3),     INTENT(IN)    :: HISBA     ! type of ISBA soil scheme
CHARACTER(LEN=3),     INTENT(IN)    :: HNAT      ! type of surface (nature, gardens)
INTEGER,              INTENT(IN)    :: KNI       ! number of points
INTEGER,              INTENT(IN)    :: KLAYER    ! number of layers
INTEGER,              INTENT(IN)    :: KPATCH    ! number of patch
REAL, DIMENSION(:,:,:), POINTER     :: PDEPTH    ! middle depth of each layer
!
!
!* local variables
!  ---------------
!
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
INTEGER           :: IRESP          ! reading return code
INTEGER           :: ILAYER         ! number of soil layers
INTEGER           :: JLAYER         ! loop counter
INTEGER           :: JPATCH         ! loop counter
!
CHARACTER(LEN=3)                     :: YISBA
LOGICAL, DIMENSION(JPCOVER)          :: GCOVER ! flag to read the covers
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZCOVER ! cover fractions
REAL,  DIMENSION(:,:,:), ALLOCATABLE :: ZD     ! depth of each inter-layer
REAL,  DIMENSION(:,:,:), ALLOCATABLE :: ZDG    ! depth of each inter-layer
REAL,  DIMENSION(:,:,:), ALLOCATABLE :: ZDEPTH ! middle of each layer for each patch
REAL,  DIMENSION(:,:),   ALLOCATABLE :: ZWORK  ! work array
REAL,  DIMENSION(KNI)                :: ZHVEG  ! high vegetation fraction
REAL,  DIMENSION(KNI)                :: ZLVEG  ! low  vegetation fraction
REAL,  DIMENSION(KNI)                :: ZNVEG  ! no   vegetation fraction
CHARACTER(LEN=4)                     :: YHVEG  ! type of high vegetation
CHARACTER(LEN=4)                     :: YLVEG  ! type of low  vegetation
CHARACTER(LEN=4)                     :: YNVEG  ! type of no   vegetation
LOGICAL                              :: GECOCLIMAP ! T if ecoclimap is used
LOGICAL                              :: GPAR_GARDEN! T if garden data are used
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_EXTERN:READ_EXTERN_DEPTH',0,ZHOOK_HANDLE)
!
IF (HNAT=='NAT') THEN
  CALL READ_LECOCLIMAP(HPROGRAM,GECOCLIMAP)
ELSE
  CALL READ_SURF(HPROGRAM,'PAR_GARDEN',GPAR_GARDEN,IRESP)
  GECOCLIMAP = .NOT. GPAR_GARDEN
END IF
!
!------------------------------------------------------------------------------
!
ALLOCATE(ZDG   (KNI,KLAYER,KPATCH))
!
IF (GECOCLIMAP) THEN
  !
  !* reading of the cover to obtain the depth of inter-layers
  !
  CALL OLD_NAME(HPROGRAM,'COVER_LIST      ',YRECFM)
  CALL READ_SURF(HPROGRAM,YRECFM,GCOVER(:),IRESP,HDIR='-')
  !
  ALLOCATE(ZCOVER(KNI,JPCOVER))
  YRECFM='COVER'
  CALL READ_SURF(HPROGRAM,YRECFM,ZCOVER(:,:),GCOVER(:),IRESP,HDIR='A')  
  !
  !* computes soil layers
  !
  CALL CONVERT_COVER_ISBA(HISBA,NUNDEF,ZCOVER,'   ',HNAT,PDG=ZDG)
  !
  DEALLOCATE(ZCOVER)
!-------------------------------------------------------------------
ELSE IF (HNAT=='NAT') THEN
  !
  !* directly read soil layers in the file for nature ISBA soil layers
  !
  ALLOCATE(ZWORK(KNI,KPATCH))
  DO JLAYER=1,KLAYER
    IF (JLAYER<10)  WRITE(YRECFM,FMT='(A4,I1.1)') 'D_DG',JLAYER
    IF (JLAYER>=10) WRITE(YRECFM,FMT='(A4,I2.2)') 'D_DG',JLAYER
    CALL READ_SURF(HPROGRAM,YRECFM,ZWORK,IRESP,HDIR='A')
    DO JPATCH=1,KPATCH
      ZDG(:,JLAYER,JPATCH) = ZWORK(:,JPATCH)
    END DO
  END DO
  DEALLOCATE(ZWORK)
!-------------------------------------------------------------------
ELSE IF (HNAT=='GRD') THEN
  !
  !* computes soil layers from vegetation fractions read in the file
  !
  CALL READ_SURF(HPROGRAM,'D_TYPE_HVEG',YHVEG,IRESP)
  CALL READ_SURF(HPROGRAM,'D_TYPE_LVEG',YLVEG,IRESP)
  CALL READ_SURF(HPROGRAM,'D_TYPE_NVEG',YNVEG,IRESP)
  CALL READ_SURF(HPROGRAM,'D_FRAC_HVEG',ZHVEG,IRESP,HDIR='A')
  CALL READ_SURF(HPROGRAM,'D_FRAC_LVEG',ZLVEG,IRESP,HDIR='A')
  CALL READ_SURF(HPROGRAM,'D_FRAC_NVEG',ZNVEG,IRESP,HDIR='A')
  ! Ground layers
  CALL GARDEN_SOIL_DEPTH(YNVEG,YLVEG,YHVEG,ZNVEG,ZLVEG,ZHVEG,ZDG)
  !
END IF
!
!-------------------------------------------------------------------
!
!* In force-restore ISBA, adds a layer at bottom of surface layer and a layer
!  between root and deep layers.
!
IF (HISBA=='2-L' .OR. HISBA=='3-L') THEN
  ILAYER = KLAYER + 1
  IF (HISBA=='3-L') ILAYER = ILAYER + 1
  ALLOCATE(ZD    (KNI,ILAYER,KPATCH))
  DO JPATCH=1,KPATCH
    ! for interpolations, middle of surface layer must be at least at 1cm
    ZD(:,1,JPATCH) = MIN(3.*ZDG(:,1,JPATCH),MAX(ZDG(:,1,JPATCH),0.02))
    ! new layer below surface layer. This layer will be at root depth layer humidity
    ZD(:,2,JPATCH) = MIN(4.*ZDG(:,1,JPATCH),0.5*(ZDG(:,1,JPATCH)+ZDG(:,2,JPATCH)))
    ! root layer
    ZD(:,3,JPATCH) = ZDG(:,2,JPATCH)
    IF (HISBA=='3-L') THEN
      ! between root and deep layers. This layer will have deep soil humidity.
      WHERE (ZDG(:,2,JPATCH)<ZDG(:,3,JPATCH))
        ZD(:,4,JPATCH) = 0.75 * ZDG(:,2,JPATCH) + 0.25 * ZDG(:,3,JPATCH)
      ELSEWHERE
        ZD(:,4,JPATCH) = ZDG(:,3,JPATCH)
      END WHERE
      ! deep layer
      ZD(:,5,JPATCH) = ZDG(:,3,JPATCH)
    END IF
  END DO
ELSE
  ILAYER = KLAYER
  ALLOCATE(ZD    (KNI,ILAYER,KPATCH))
  ZD = ZDG
END IF
!
DEALLOCATE(ZDG)
!
!-------------------------------------------------------------------
!* recovers middle layer depth (from the surface)
ALLOCATE(ZDEPTH    (KNI,ILAYER,KPATCH))
ZDEPTH = XUNDEF
DO JPATCH=1,KPATCH
  WHERE(ZD(:,1,JPATCH)/=XUNDEF) &
    ZDEPTH    (:,1,JPATCH)=ZD(:,1,JPATCH)/2.  
  DO JLAYER=2,ILAYER
    WHERE(ZD(:,1,JPATCH)/=XUNDEF) &
      ZDEPTH    (:,JLAYER,JPATCH) = (ZD(:,JLAYER-1,JPATCH) + ZD(:,JLAYER,JPATCH))/2.  
  END DO
END DO
DEALLOCATE(ZD)
!-------------------------------------------------------------------
!
ALLOCATE(PDEPTH    (KNI,ILAYER,NVEGTYPE))
CALL PUT_ON_ALL_VEGTYPES(KNI,ILAYER,KPATCH,NVEGTYPE,ZDEPTH,PDEPTH)
DEALLOCATE(ZDEPTH)

IF (LHOOK) CALL DR_HOOK('MODE_READ_EXTERN:READ_EXTERN_DEPTH',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------
!
END SUBROUTINE READ_EXTERN_DEPTH
!
!
!-------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!
!     #######################
      SUBROUTINE READ_EXTERN_ISBA(HPROGRAM,KLUOUT,KNI,HFIELD,PFIELD,PDEPTH)
!     #######################
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODI_READ_SURF
USE MODE_SOIL
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
!
CHARACTER(LEN=6),     INTENT(IN)    :: HPROGRAM  ! program calling surf. schemes
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
INTEGER,              INTENT(IN)    :: KNI       ! number of points
CHARACTER(LEN=7),     INTENT(IN)    :: HFIELD    ! field name
REAL, DIMENSION(:,:,:), POINTER       :: PFIELD    ! field to initialize
REAL, DIMENSION(:,:,:), POINTER       :: PDEPTH    ! middle depth of each layer
!
!
!* local variables
!  ---------------
!
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4)  :: YLVL
CHARACTER(LEN=3)  :: YISBA          ! type of ISBA soil scheme
CHARACTER(LEN=3)  :: YNAT           ! type of surface (nature, garden)
CHARACTER(LEN=2)  :: YDIF           ! type of Theta(psi) function
CHARACTER(LEN=4)  :: YPEDOTF        ! type of pedo-transfert function
INTEGER           :: IRESP          ! reading return code
INTEGER           :: ILAYER         ! number of layers
INTEGER           :: JLAYER         ! loop counter
INTEGER           :: IPATCH         ! number of patch
INTEGER           :: JPATCH         ! loop counter
INTEGER           :: JVEGTYPE       ! loop counter
!
REAL,  DIMENSION(:,:,:), ALLOCATABLE :: ZFIELD ! field read, one level, all patches
REAL,  DIMENSION(:,:),   ALLOCATABLE :: ZWORK  ! field read, one level, all patches
!
REAL,  DIMENSION(:,:,:), ALLOCATABLE :: ZVAR      ! profile of physical variable
REAL,  DIMENSION(:),   ALLOCATABLE   :: ZCLAY     ! clay fraction
REAL,  DIMENSION(:),   ALLOCATABLE   :: ZSAND     ! sand fraction
REAL,  DIMENSION(:),   ALLOCATABLE   :: ZWWILT    ! wilting point
REAL,  DIMENSION(:),   ALLOCATABLE   :: ZWFC      ! field capacity
REAL,  DIMENSION(:),   ALLOCATABLE   :: ZWSAT     ! saturation
!
INTEGER :: IVERSION   ! surface version
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_EXTERN:READ_EXTERN_ISBA',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,*) ' | Reading ',HFIELD,' in externalized file'
!
!------------------------------------------------------------------------------
!
!* Read number of soil layers
!
YRECFM='GROUND_LAYER'
IF (HFIELD(1:3)=='TWN') YRECFM='TWN_LAYER'
CALL READ_SURF(HPROGRAM,YRECFM,ILAYER,IRESP)
!
!* number of tiles
!
IPATCH=1
IF (HFIELD(1:3)/='TWN') THEN
  YRECFM='PATCH_NUMBER'
  CALL READ_SURF(HPROGRAM,YRECFM,IPATCH,IRESP)
END IF
!
!* soil scheme
!
YRECFM='ISBA'
IF (HFIELD(1:3)=='TWN') YRECFM='TWN_ISBA'
CALL READ_SURF(HPROGRAM,YRECFM,YISBA,IRESP)
!
YRECFM='VERSION'
CALL READ_SURF(HPROGRAM,YRECFM,IVERSION,IRESP)
!
IF (IVERSION>=7) THEN
  !
  !* Theta(psi) function
  !
  YRECFM='THETA_PSI'
  IF (HFIELD(1:3)=='TWN') YRECFM='TWN_THETAPSI'
  CALL READ_SURF(HPROGRAM,YRECFM,YDIF,IRESP)
  !
  !* Pedo-transfert function
  !
  YRECFM='PEDOTF'
  IF (HFIELD(1:3)=='TWN') YRECFM='TWN_PEDOTF'
  CALL READ_SURF(HPROGRAM,YRECFM,YPEDOTF,IRESP)
  !
ELSE
  YDIF = 'BC'
  YPEDOTF = 'CH78'
ENDIF
!
!Only Brook and Corey with Force-Restore scheme
IF(YISBA/='DIF')THEN
  YDIF='BC'
  YPEDOTF='CH78'
ENDIF
!
!* Allocate soil variable profile
!  ------------------------------
!
ALLOCATE(ZVAR(KNI,ILAYER,IPATCH))
ALLOCATE(ZWORK(KNI,IPATCH))
!
! *.  Read soil variable profile
!     --------------------------
!
DO JLAYER=1,ILAYER
  WRITE(YLVL,'(I4)') JLAYER
  YRECFM=TRIM(HFIELD)//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP,HDIR='A')
  DO JPATCH=1,IPATCH
    ZVAR(:,JLAYER,JPATCH)=ZWORK(:,JPATCH)
  END DO
END DO
!
DEALLOCATE(ZWORK)
!
! *.  Read soil grid
!     --------------
!
IF ((HFIELD=='TG    ' .OR. HFIELD=='TWN_TG') .AND. (YISBA=='2-L' .OR. YISBA=='3-L')) THEN
  ALLOCATE(PDEPTH    (KNI,ILAYER,NVEGTYPE))
  DO JVEGTYPE=1,NVEGTYPE
    PDEPTH(:,1,JVEGTYPE) = 0.
    PDEPTH(:,2,JVEGTYPE) = 0.2
    IF (ILAYER==3) PDEPTH(:,3,JVEGTYPE) = 3.
  END DO
ELSE
  YNAT='NAT'
  IF (HFIELD(1:3)=='TWN') YNAT='GRD'
  CALL READ_EXTERN_DEPTH(HPROGRAM,KLUOUT,YISBA,YNAT,KNI,ILAYER,IPATCH,PDEPTH)
END IF
!
!-------------------------------------------------------------------------------
!
! *.  Read clay fraction
!     ------------------
!
ALLOCATE(ZCLAY(KNI))
YRECFM='CLAY'
IF (HFIELD(1:3)=='TWN') YRECFM='TWN_CLAY'
CALL READ_SURF(HPROGRAM,YRECFM,ZCLAY(:),IRESP,HDIR='A')
!
!-------------------------------------------------------------------------------
!
! *.  Read sand fraction
!     ------------------
!
ALLOCATE(ZSAND(KNI))
YRECFM='SAND'
IF (HFIELD(1:3)=='TWN') YRECFM='TWN_SAND'
CALL READ_SURF(HPROGRAM,YRECFM,ZSAND(:),IRESP,HDIR='A')
!
!-------------------------------------------------------------------------------
!
! *.  Compute relative humidity from units kg/m^2 (SWI)
!     ------------------------------------------------
!
!* In case of force-restore ISBA, adds one layer at bottom of surface layer
IF ((HFIELD=='WG    ' .OR. HFIELD=='TWN_WG  ' .OR. HFIELD=='TWN_WGI ' .OR. HFIELD=='WGI   ') &
     .AND. (YISBA=='2-L' .OR. YISBA=='3-L')) THEN
  ALLOCATE(ZFIELD(KNI,ILAYER,IPATCH))
  ZFIELD(:,:,:) = ZVAR(:,:,:)
  DEALLOCATE(ZVAR)
  !
  ILAYER = ILAYER + 1
  IF ( YISBA=='3-L' ) ILAYER = ILAYER + 1
  ALLOCATE(ZVAR(KNI,ILAYER,IPATCH))
  DO JPATCH=1,IPATCH
    ZVAR(:,1,JPATCH)=ZFIELD(:,1,JPATCH)
    ZVAR(:,2,JPATCH)=ZFIELD(:,2,JPATCH)  ! new layer at root layer humidity but below surface layer
    ZVAR(:,3,JPATCH)=ZFIELD(:,2,JPATCH)
    IF ( YISBA=='3-L' ) THEN
      ZVAR(:,4,JPATCH)=ZFIELD(:,3,JPATCH)
      ZVAR(:,5,JPATCH)=ZFIELD(:,3,JPATCH)
    END IF
  END DO
  DEALLOCATE(ZFIELD)
END IF
!
ALLOCATE(ZFIELD(KNI,ILAYER,IPATCH))
ZFIELD = ZVAR
!
IF (HFIELD=='WG    ' .OR. HFIELD=='WGI   ' &
  .OR. HFIELD=='TWN_WG  ' .OR. HFIELD=='TWN_WGI ') THEN
  !
  ! Compute ISBA model constants
  !
  ALLOCATE (ZWFC  (KNI))
  ALLOCATE (ZWWILT(KNI))
  ALLOCATE (ZWSAT (KNI))
  !
  IF(YDIF=='BC')THEN
     ZWSAT (:) = WSAT_FUNC (ZCLAY(:),ZSAND(:),YPEDOTF)
     ZWWILT(:) = WWILT_FUNC(ZCLAY(:),ZSAND(:),YPEDOTF)
     ZWFC  (:) = WFC_FUNC  (ZCLAY(:),ZSAND(:),YPEDOTF)
  ELSE
     WRITE  (KLUOUT,*) 'No yet implemented with Van Genuchten function'
     STOP
  ENDIF
  !
  DEALLOCATE (ZSAND)
  DEALLOCATE (ZCLAY)

  ZFIELD(:,:,:) = XUNDEF
  !
  IF (HFIELD=='WG    '.OR. HFIELD=='TWN_WG  ') THEN
    DO JPATCH=1,IPATCH
      DO JLAYER=1,ILAYER
        WHERE(ZVAR(:,JLAYER,JPATCH)/=XUNDEF)
          ZVAR(:,JLAYER,JPATCH) = MAX(MIN(ZVAR(:,JLAYER,JPATCH),ZWSAT(:)),0.)
          !
          ZFIELD(:,JLAYER,JPATCH) = (ZVAR(:,JLAYER,JPATCH) - ZWWILT(:)) / (ZWFC(:) - ZWWILT(:))
        END WHERE
      END DO
    END DO
  ELSE IF (HFIELD=='WGI   '.OR. HFIELD=='TWN_WGI ') THEN
    DO JPATCH=1,IPATCH
      DO JLAYER=1,ILAYER
        WHERE(ZVAR(:,JLAYER,JPATCH)/=XUNDEF) &
          ZFIELD(:,JLAYER,JPATCH) = ZVAR(:,JLAYER,JPATCH) / ZWSAT(:)  
      END DO
    END DO
  END IF
!
  DEALLOCATE (ZWSAT)
  DEALLOCATE (ZWWILT)
  DEALLOCATE (ZWFC)
!
!
END IF
!
DEALLOCATE(ZVAR)
!-------------------------------------------------------------------------------
!
! *.  Set the field on all vegtypes
!     -----------------------------
!
ALLOCATE(PFIELD(KNI,ILAYER,NVEGTYPE))
CALL PUT_ON_ALL_VEGTYPES(KNI,ILAYER,IPATCH,NVEGTYPE,ZFIELD,PFIELD)
DEALLOCATE(ZFIELD)
IF (LHOOK) CALL DR_HOOK('MODE_READ_EXTERN:READ_EXTERN_ISBA',1,ZHOOK_HANDLE)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE READ_EXTERN_ISBA
!
!------------------------------------------------------------------------------
!
END MODULE MODE_READ_EXTERN                       
