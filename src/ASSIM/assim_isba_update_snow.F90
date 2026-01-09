!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE ASSIM_ISBA_UPDATE_SNOW (IO, S, NP, NPE, HPROGRAM, KI, PSWE, HTEST )

! ------------------------------------------------------------------------------------------
!  *****************************************************************************************
!
!  Routine to update snow field for ISBA
!  Trygve Aspelien, Separating IO  06/2013
!
!
! ******************************************************************************************
! ------------------------------------------------------------------------------------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_NP_t, ISBA_NPE_t, ISBA_S_t, ISBA_PE_t, ISBA_P_t
!
USE MODD_CSTS,        ONLY : XTT
USE MODD_ASSIM,       ONLY : LSWE, LSDPSINI, XSDPSINI, LSDPSMIN, XSDPSMIN
USE MODD_DATA_COVER_PAR, ONLY : NVT_SNOW
USE MODD_SURF_PAR,    ONLY : XUNDEF
USE MODD_SNOW_PAR,    ONLY : XANSMIN, XANSMAX, XRHOSMIN, XRHOSMAX, XSNOWDMIN
USE MODE_SNOW3L,      ONLY : SNOW3LGRID

USE MODI_SNOW_T_WLIQ_TO_HEAT
USE MODI_PACK_SAME_RANK
USE MODI_ABOR1_SFX
USE MODI_ASSIM_GATHER_WRITE_INCREMENTS
!
USE YOMHOOK,          ONLY : LHOOK,DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t),       INTENT(INOUT) :: S
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t),     INTENT(INOUT) :: NPE
!
CHARACTER(LEN=6),    INTENT(IN)    :: HPROGRAM  ! program calling surf. schemes
INTEGER,             INTENT(IN)    :: KI
REAL, DIMENSION(KI), INTENT(IN)    :: PSWE
CHARACTER(LEN=2),    INTENT(IN)    :: HTEST     ! must be equal to 'OK'
!
!    Declarations of local variables
!
TYPE(ISBA_P_t),      POINTER       :: PK
TYPE(ISBA_PE_t),     POINTER       :: PEK
!
CHARACTER(LEN=2)                :: CP
REAL, ALLOCATABLE, DIMENSION(:) :: ZSWE       ! Patch value of updated snow
REAL, ALLOCATABLE, DIMENSION(:) :: ZSWE_ORIG  ! Patch value of 1. guess
REAL, ALLOCATABLE, DIMENSION(:) :: ZSWEINC 
REAL, ALLOCATABLE, DIMENSION(:) :: ZD_IN
REAL, ALLOCATABLE, DIMENSION(:) :: ZSNA,ZSNA_ORIG
REAL, ALLOCATABLE, DIMENSION(:) :: ZSNR,ZSNR_ORIG
REAL, ALLOCATABLE, DIMENSION(:) :: ZD,ZHEAT,ZD_ORIG, ZSNOWD
REAL, ALLOCATABLE, DIMENSION(:,:) :: ZSNOWDZN
REAL, ALLOCATABLE, DIMENSION(:,:) :: ZSNOWRHO, ZSNOWAGE, ZSNOWHEAT
REAL, ALLOCATABLE, DIMENSION(:)   :: ZSNOWALB
LOGICAL,ALLOCATABLE,DIMENSION(:)  :: OMASK, OSNOWDEF
INTEGER, ALLOCATABLE,DIMENSION(:) :: JS_MASK
INTEGER                           :: JS_POINTS

INTEGER  :: JL,JP,JI,JJ
REAL :: ZSWEMIN,ZSWEMAX,ZSWEMEAN,ZINCMIN,ZINCMAX,ZINCMEAN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
! ----------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_UPDATE_SNOW',0,ZHOOK_HANDLE)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_ISBA_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF

! Loop patches
DO JP=1,IO%NPATCH

  PK  => NP%AL(JP)
  PEK => NPE%AL(JP)
  ALLOCATE(ZSWE(PK%NSIZE_P))
  ALLOCATE(ZSWE_ORIG(PK%NSIZE_P))
  ALLOCATE(ZD_ORIG(PK%NSIZE_P))
  ALLOCATE(ZSWEINC(PK%NSIZE_P))
  ALLOCATE(ZD_IN(PK%NSIZE_P))
  ALLOCATE(ZSNA(PK%NSIZE_P))
  ALLOCATE(ZSNA_ORIG(PK%NSIZE_P))
  ALLOCATE(ZSNR(PK%NSIZE_P))
  ALLOCATE(ZSNR_ORIG(PK%NSIZE_P))
  ALLOCATE(OMASK(PK%NSIZE_P))
  ALLOCATE(OSNOWDEF(PK%NSIZE_P))

  ZSWE(:)=0.
  ! Pack ISBA fields to this patch
  CALL PACK_SAME_RANK(PK%NR_P,PSWE,ZD_IN)

  ! Make sure we have positive values
  WHERE( ZD_IN(:) < 0.0 )
    ZD_IN(:) = 0.0
  ENDWHERE
  ! Set mask on defined input snow field
  WHERE( ZD_IN(:) /= -999. )
    OSNOWDEF(:) = .TRUE.
  ELSEWHERE
    OSNOWDEF(:) = .FALSE.
  ENDWHERE

  ! Default snow depth
  ZD_ORIG = 0.
  ZD_ORIG(:) = SUM(PEK%TSNOW%WSNOW(:,:)/PEK%TSNOW%RHO(:,:),MASK=(PEK%TSNOW%RHO /= XUNDEF),DIM=2)
  ! Default SWE
  ZSWE_ORIG(:) = SUM(PEK%TSNOW%WSNOW(:,:),DIM=2)
  ! Find average snow pack density if more than one layer
  IF ( PEK%TSNOW%NLAYER > 1 ) THEN
    IF ( PEK%TSNOW%SCHEME == '3-L' ) THEN
      ZSNR(:) = XUNDEF
      WHERE ( ZD_ORIG(:) /= 0. )
        ZSNR(:) = ZSWE_ORIG(:)/ZD_ORIG(:)
      ENDWHERE
    ELSE
       CALL ABOR1_SFX('ASSIM_ISBA_n: THE SNOW SCHEME '//TRIM(PEK%TSNOW%SCHEME)// &
 & ' IS A MULTILAYER SCHEME BUT NO CONVERSION FROM DEPTH TO METER IS IMPLEMENTED')
    ENDIF
  ELSE
    ZSNR(:)=PEK%TSNOW%RHO(:,1)
  ENDIF

  ! Convert analysed SWE to SD
  IF ( LSWE ) THEN
    WHERE( OSNOWDEF )
      WHERE( ZD_IN(:) /= 0. )
        WHERE( ZSNR(:) /= XUNDEF )
          ZD_IN(:) = ZD_IN(:) / ZSNR(:)
        ELSEWHERE
          ZD_IN(:) = ZD_IN(:) / ( 0.5 * ( XRHOSMIN + XRHOSMAX ))
        ENDWHERE
      ENDWHERE
    ENDWHERE
  ENDIF

  ! Set snow=0 where 1. guess = 0 and Ts>0, to avoid that the snow analysis introduce snow where it is no snow.
  WHERE ( ZSWE_ORIG(:)<1.0E-10 .AND. PEK%XTG(:,1) >XTT )
    ZD_IN(:)   = 0.0
  END WHERE

  ! Set initial SD = XSDPSMIN on vegetation type PERMANENT SNOW
  IF ( LSDPSINI ) THEN
    WHERE ( OSNOWDEF )
      WHERE (PK%XVEGTYPE_PATCH(:,NVT_SNOW) > 0.5)
         ZD_IN(:) = XSDPSINI*PK%XVEGTYPE_PATCH(:,NVT_SNOW)
      ENDWHERE
    ENDWHERE
  ENDIF

  ! Make sure that the SD >= XSDPSMIN on vegetation type PERMANENT SNOW
  IF ( LSDPSMIN ) THEN
    WHERE ( OSNOWDEF )
      WHERE ( PK%XVEGTYPE_PATCH(:,NVT_SNOW) > 0.75 .AND. ZD_IN(:) < (XSDPSMIN*PK%XVEGTYPE_PATCH(:,NVT_SNOW))/50.0)
         ZD_IN(:) = XSDPSMIN*PK%XVEGTYPE_PATCH(:,NVT_SNOW)
      ENDWHERE
    ENDWHERE
  ENDIF

  ! Snow albedo and density are given initial values in points  
  ! which get initial snow in the snow analysis
  OMASK=.FALSE.
  ZSNA(:)=PEK%TSNOW%ALB(:)
  ZSNA_ORIG(:) = ZSNA(:)
  ZSNR_ORIG(:) = SUM(PEK%TSNOW%RHO(:,:), MASK=(PEK%TSNOW%RHO/=XUNDEF), DIM=2)
  WHERE ( OSNOWDEF )
    WHERE ( ZSWE_ORIG(:) < 1.0E-10 .AND. ZD_IN(:)>= 1.0E-10 )
      ZSNA(:)    = 0.5 * ( XANSMIN + XANSMAX )
      ZSNR(:)    = 0.5 * ( XRHOSMIN + XRHOSMAX )
      ! Create a mask for new snow points
      OMASK(:)=.TRUE.
    END WHERE
  END WHERE

  IF ( JP < 10 ) THEN
   WRITE(CP(1:1),'(I1)') JP
   CP(2:2)=" "
  ELSE
   WRITE(CP(1:2),'(I2)') JP 
  ENDIF

  IF ( PEK%TSNOW%SCHEME == 'D95' ) THEN

    WHERE( OSNOWDEF )
      ! Update snow/albedo/density
      PEK%TSNOW%ALB(:)=ZSNA(:)

      ! Only modify density for new snow
      WHERE ( OMASK(:) )
        PEK%TSNOW%RHO(:,1)=ZSNR(:)
      ENDWHERE

      ! Set SWE
      PEK%TSNOW%WSNOW(:,1) = ZD_IN(:) * ZSNR(:)
    ENDWHERE

    ! Find modified snow depth in grid points
    ALLOCATE(ZD(PK%NSIZE_P))
    ZD=0.
    WHERE ( PEK%TSNOW%RHO(:,1) /= XUNDEF )
      ZD(:) = PEK%TSNOW%WSNOW(:,1) / PEK%TSNOW%RHO(:,1)
    ENDWHERE

    ! Print snow depth increments
    CALL ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,"DSN_VEG_P"//CP,ZD_ORIG,ZD,XUNDEF,LSTAT=.TRUE.)
    DEALLOCATE(ZD)

    ! Print increments
    ZSNA(:) = PEK%TSNOW%ALB(:)
    ZSWE(:) = SUM(PEK%TSNOW%WSNOW(:,:),DIM=2)
    ZSNR(:) = SUM(PEK%TSNOW%RHO(:,:),MASK=(PEK%TSNOW%RHO/=XUNDEF),DIM=2)
    CALL ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,"ASN_VEG_P"//CP,ZSNA_ORIG,ZSNA,XUNDEF,LSTAT=.TRUE.)
    CALL ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,"WSN_VEG_P"//CP,ZSWE_ORIG,ZSWE,XUNDEF,LSTAT=.TRUE.)
    CALL ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,"RSN_VEG_P"//CP,ZSNR_ORIG,ZSNR,XUNDEF,LSTAT=.TRUE.)

  ! ISBA-ES need heat and possible rescaling of layers
  ELSE IF ( PEK%TSNOW%SCHEME == '3-L' .OR. PEK%TSNOW%SCHEME =='1-L' ) THEN

    ZD_ORIG=0.
    ZD_ORIG(:) = SUM(PEK%TSNOW%WSNOW(:,:)/PEK%TSNOW%RHO(:,:),MASK=(PEK%TSNOW%RHO /= XUNDEF),DIM=2)
    ZSWE_ORIG(:) = SUM(PEK%TSNOW%WSNOW(:,:),DIM=2)

    ! Update snow/albedo/density for new snow
    ALLOCATE(ZHEAT(PK%NSIZE_P))

    WHERE( OSNOWDEF )
      WHERE ( OMASK(:))
        PEK%TSNOW%ALB(:)=ZSNA(:)
      ENDWHERE
    ENDWHERE

    ! Initialize density and heat for new snow
    DO JL=1,PEK%TSNOW%NLAYER
      WHERE( OSNOWDEF )
        WHERE ( OMASK(:))
          PEK%TSNOW%RHO(:,JL) = 0.5 * ( XRHOSMIN + XRHOSMAX )
        ENDWHERE
      ENDWHERE
      CALL SNOW_T_WLIQ_TO_HEAT(ZHEAT,PEK%TSNOW%RHO(:,JL),PEK%XTG(:,1))
      WHERE( OSNOWDEF )
        WHERE ( OMASK(:))
          PEK%TSNOW%HEAT(:,JL) = ZHEAT(:)
          PEK%TSNOW%AGE(:,JL) = 0
        ENDWHERE
      ENDWHERE
    ENDDO
    DEALLOCATE(ZHEAT)

    ! Need to work only for snow points for ISBA-ES routines
    ! Set snow mask
    JS_POINTS = 0
    DO JI=1,PK%NSIZE_P
      IF ( OSNOWDEF(JI) ) THEN
        IF ( ZD_IN(JI) > XSNOWDMIN ) THEN
          JS_POINTS = JS_POINTS + 1
        ENDIF
      ENDIF
    ENDDO

    ALLOCATE(JS_MASK(JS_POINTS))
    JS_MASK = -1
    JJ = 0
    DO JI=1,PK%NSIZE_P
      IF ( OSNOWDEF(JI) ) THEN
        IF ( ZD_IN(JI) > XSNOWDMIN ) THEN
          JJ = JJ + 1
          JS_MASK(JJ) = JI
        ENDIF
      ENDIF
    ENDDO

    IF ( JJ /= JS_POINTS ) THEN
       CALL ABOR1_SFX('Mismatch in snow points')
    ENDIF

    ! Find modified snow depth in grid points
    ALLOCATE(ZSNOWDZN(JS_POINTS,PEK%TSNOW%NLAYER))
    ALLOCATE(ZSNOWD(JS_POINTS))
    ALLOCATE(ZSNOWRHO(JS_POINTS,PEK%TSNOW%NLAYER))
    ALLOCATE(ZSNOWHEAT(JS_POINTS,PEK%TSNOW%NLAYER))
    ALLOCATE(ZSNOWAGE(JS_POINTS,PEK%TSNOW%NLAYER))
    ALLOCATE(ZSNOWALB(JS_POINTS))

    DO JI=1,JS_POINTS
      ZSNOWALB(JI)    = PEK%TSNOW%ALB(JS_MASK(JI))
      ZSNOWRHO(JI,:)  = PEK%TSNOW%RHO(JS_MASK(JI),:)
      ZSNOWHEAT(JI,:) = PEK%TSNOW%HEAT(JS_MASK(JI),:)
      ZSNOWAGE(JI,:)  = PEK%TSNOW%AGE(JS_MASK(JI),:)
      ZSNOWD(JI) = ZD_IN(JS_MASK(JI))
    ENDDO

    !* Rescale snow pack after update for multilayer scheme
    IF ( PEK%TSNOW%SCHEME == '3-L' ) THEN
      ! Re grid
       CALL SNOW3LGRID(ZSNOWDZN(:,:),ZSNOWD(:))

      ! First initialize SWE to zero unless undefined
       DO JI = 1,PK%NSIZE_P
          IF ( OSNOWDEF(JI) ) THEN
            PEK%TSNOW%WSNOW(JI,:) = 0.
            PEK%TSNOW%RHO(JI,:)   = XUNDEF
            PEK%TSNOW%HEAT(JI,:)  = XUNDEF
            PEK%TSNOW%AGE(JI,:)   = XUNDEF
            PEK%TSNOW%ALB(JI)     = XUNDEF
          ENDIF
      ENDDO
      ! Map back to full patch field with mask
      DO JI = 1,JS_POINTS
        PEK%TSNOW%WSNOW(JS_MASK(JI),:) = ZSNOWDZN(JI,:)*ZSNOWRHO(JI,:)
        PEK%TSNOW%RHO(JS_MASK(JI),:)   = ZSNOWRHO(JI,:)
        PEK%TSNOW%HEAT(JS_MASK(JI),:)  = ZSNOWHEAT(JI,:)
        PEK%TSNOW%AGE(JS_MASK(JI),:)   = ZSNOWAGE(JI,:)
        PEK%TSNOW%ALB(JS_MASK(JI))     = ZSNOWALB(JI)
      ENDDO
    ENDIF
    DEALLOCATE(ZSNOWDZN)
    DEALLOCATE(JS_MASK)
    DEALLOCATE(ZSNOWD)
    DEALLOCATE(ZSNOWRHO)
    DEALLOCATE(ZSNOWHEAT)
    DEALLOCATE(ZSNOWAGE)
    DEALLOCATE(ZSNOWALB)

    ! Finished snow pack
    ALLOCATE(ZD(PK%NSIZE_P))
    ZD = 0.
    ZD(:) = SUM(PEK%TSNOW%WSNOW(:,:)/PEK%TSNOW%RHO(:,:),MASK=(PEK%TSNOW%RHO/=XUNDEF),DIM=2)

    ! Print snow depth increments
    CALL ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,"DSN_VEG_P"//CP,ZD_ORIG,ZD,XUNDEF,LSTAT=.TRUE.)
    DEALLOCATE(ZD)

    ! Print increments
    ZSNA(:) = PEK%TSNOW%ALB(:)
    ZSWE(:) = SUM(PEK%TSNOW%WSNOW(:,:),DIM=2)
    ZSNR(:) = SUM(PEK%TSNOW%RHO(:,:),MASK=(PEK%TSNOW%RHO/=XUNDEF),DIM=2)
    CALL ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,"ASN_VEG_P"//CP,ZSNA_ORIG,ZSNA,XUNDEF,LSTAT=.TRUE.)
    CALL ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,"WSN_VEG_P"//CP,ZSWE_ORIG,ZSWE,XUNDEF,LSTAT=.TRUE.)
    CALL ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,"RSN_VEG_P"//CP,ZSNR_ORIG,ZSNR,XUNDEF,LSTAT=.TRUE.)

  ELSE
    CALL ABOR1_SFX('ASSIM_ISBA_n: THE SNOW SCHEME '//TRIM(PEK%TSNOW%SCHEME)//&
       & ' HAS NO SD UPDATE IMPLEMENTED')
  ENDIF

  DEALLOCATE(OMASK)
  DEALLOCATE(ZSWE)
  DEALLOCATE(ZSWE_ORIG)
  DEALLOCATE(ZD_ORIG)
  DEALLOCATE(ZSWEINC)
  DEALLOCATE(ZD_IN)
  DEALLOCATE(ZSNA)
  DEALLOCATE(ZSNA_ORIG)
  DEALLOCATE(ZSNR)
  DEALLOCATE(ZSNR_ORIG)
  DEALLOCATE(OSNOWDEF)
ENDDO


!
! -------------------------------------------------------------------------------------
 IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_UPDATE_SNOW',1,ZHOOK_HANDLE)
 END SUBROUTINE ASSIM_ISBA_UPDATE_SNOW

