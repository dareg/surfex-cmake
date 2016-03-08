!     #########
SUBROUTINE PREP_HOR_ISBA_CC_FIELD (DTCO, U, KLAT, IP, O, R, PLAI, PVEGTYPE, &
                                   HPROGRAM,HSURF,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!     #################################################################################
!
!!****  *PREP_HOR_ISBA_CC_FIELD* - reads, interpolates and prepares an ISBA-CC field
!                                   only external case implemeted
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     B. Decharme
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2014
!!------------------------------------------------------------------
!
!
!
!
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODD_CO2V_PAR,  ONLY : XCA_NIT, XCC_NIT
!
USE MODD_PREP,      ONLY : LINTERP, CMASK
!

USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF,NUNDEF
!
USE MODI_READ_PREP_ISBA_CONF
USE MODI_ABOR1_SFX
USE MODI_HOR_INTERPOL
USE MODI_VEGTYPE_GRID_TO_PATCH_GRID
USE MODI_GET_LUOUT
USE MODI_PREP_ISBA_CC_EXTERN
USE MODI_PUT_ON_ALL_VEGTYPES
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
INTEGER, INTENT(IN) :: KLAT
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: O
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
REAL, DIMENSION(:,:), INTENT(INOUT) :: PLAI
REAL, DIMENSION(:,:), INTENT(INOUT) :: PVEGTYPE
TYPE(ISBA_PROG_t), INTENT(INOUT) :: R
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=8),   INTENT(IN)  :: HSURF     ! type of field
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
!
!*      0.2    declarations of local variables
!
 CHARACTER(LEN=6)              :: YFILETYPE ! type of input file
 CHARACTER(LEN=28)             :: YFILE     ! name of file
 CHARACTER(LEN=6)              :: YFILEPGDTYPE ! type of input file
 CHARACTER(LEN=28)             :: YFILEPGD     ! name of file
REAL, POINTER, DIMENSION(:,:,:)     :: ZFIELDIN  ! field to interpolate horizontally
REAL, POINTER, DIMENSION(:,:)       :: ZFIELD    ! field to interpolate horizontally
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZFIELDOUTP ! field interpolated   horizontally
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZFIELDOUTV ! field interpolated   horizontally
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZW        ! work array (x, fine   soil grid, npatch)
!
INTEGER                       :: ILUOUT    ! output listing logical unit
!
LOGICAL                       :: GUNIF     ! flag for prescribed uniform field
LOGICAL                       :: GPREP_AGS ! flag to prepare ags field (only external case implemeted)
!
INTEGER                       :: JPATCH    ! loop on patches
INTEGER                       :: JVEGTYPE  ! loop on vegtypes
INTEGER                       :: INL, INP, JJ, JL ! Work integer
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_HOR_ISBA_CC_FIELD',0,ZHOOK_HANDLE)
!
!*      1.     Reading of input file name and type
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
CALL READ_PREP_ISBA_CONF(HPROGRAM,HSURF,YFILE,YFILETYPE,YFILEPGD,YFILEPGDTYPE,   &
                         HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,ILUOUT,GUNIF)
!
CMASK = 'NATURE'
!
GPREP_AGS = .TRUE.
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading of input  configuration (Grid and interpolation type)
!
IF (GUNIF) THEN
   GPREP_AGS = .FALSE.
ELSE IF (YFILETYPE=='ASCLLV') THEN
   GPREP_AGS = .FALSE.
ELSE IF (YFILETYPE=='GRIB  ') THEN
   GPREP_AGS = .FALSE.
ELSE IF (YFILETYPE=='MESONH' .OR. YFILETYPE=='ASCII ' .OR. YFILETYPE=='LFI   '.OR.YFILETYPE=='FA    ') THEN
   CALL PREP_ISBA_CC_EXTERN(HPROGRAM,HSURF,YFILE,YFILETYPE,YFILEPGD,YFILEPGDTYPE,ILUOUT,ZFIELDIN,GPREP_AGS)
ELSE IF (YFILETYPE=='BUFFER') THEN
   GPREP_AGS = .FALSE.
ELSE IF (YFILETYPE=='NETCDF') THEN
   GPREP_AGS = .FALSE.
ELSE
   CALL ABOR1_SFX('PREP_HOR_ISBA_CC_FIELD: data file type not supported : '//YFILETYPE)
END IF
!
!-------------------------------------------------------------------------------------
!
!*      3.     Horizontal interpolation
!
IF(GPREP_AGS)THEN
!
  INL = SIZE(ZFIELDIN,2)
  INP = SIZE(ZFIELDIN,3)
!
  ALLOCATE(ZFIELDOUTP(KLAT,INL,INP))
  ALLOCATE(ZFIELD(SIZE(ZFIELDIN,1),INL))
!
  DO JPATCH = 1, INP
    ZFIELD(:,:)=ZFIELDIN(:,:,JPATCH)
    IF (INP==NVEGTYPE) THEN
       LINTERP = (PVEGTYPE(:,JPATCH) > 0.)
    ELSEIF(INP==O%NPATCH)THEN
       LINTERP = (IP%XPATCH(:,JPATCH) > 0.)
    ENDIF
    CALL HOR_INTERPOL(DTCO, U, ILUOUT,ZFIELD,ZFIELDOUTP(:,:,JPATCH))
    LINTERP = .TRUE.
  END DO
!
  DEALLOCATE(ZFIELD)
  DEALLOCATE(ZFIELDIN)
!
  ALLOCATE(ZFIELDOUTV(KLAT,INL,NVEGTYPE))
!
  CALL PUT_ON_ALL_VEGTYPES(KLAT,INL,INP,NVEGTYPE,ZFIELDOUTP,ZFIELDOUTV)
!
  DEALLOCATE(ZFIELDOUTP)
!
!
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      6.     Transformation from vegtype grid to patch grid
!
IF(GPREP_AGS)THEN
!
  ALLOCATE(ZW (KLAT,SIZE(ZFIELDOUTV,2),O%NPATCH))
!
  ZW(:,:,:) = 0.
  CALL VEGTYPE_GRID_TO_PATCH_GRID(O%NPATCH,IP%XVEGTYPE_PATCH,IP%XPATCH,ZFIELDOUTV,ZW)
!
ELSE
!
  SELECT CASE (HSURF)
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  
    CASE('BIOMASS') 
     ALLOCATE(ZW(KLAT,O%NNBIOMASS,O%NPATCH))
     ZW(:,:,:) = 0.
     WHERE(PLAI(:,:)/=XUNDEF)
       ZW(:,1,:) = PLAI(:,:) * IP%XBSLAI_NITRO(:,:)
     ENDWHERE
     ZW(:,2,:) = MAX( 0., (ZW(:,1,:)/ (XCC_NIT/10.**XCA_NIT))  &
                          **(1.0/(1.0-XCA_NIT)) - ZW(:,1,:) )
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    CASE('LITTER') 
     ALLOCATE(ZW(KLAT,O%NNLITTER*O%NNLITTLEVS,O%NPATCH))
     ZW(:,:,:) = 0.0
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    CASE('SOILCARB') 
     ALLOCATE(ZW(KLAT,O%NNSOILCARB,O%NPATCH))
     ZW(:,:,:) = 0.0
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    CASE('LIGNIN') 
     ALLOCATE(ZW(KLAT,O%NNLITTLEVS,O%NPATCH))
     ZW(:,:,:) = 0.0
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
  END SELECT
!
ENDIF
!-------------------------------------------------------------------------------------
!
!*      7.     Return to historical variable
!
!
SELECT CASE (HSURF)
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('BIOMASS') 
  ALLOCATE(R%XBIOMASS(KLAT,O%NNBIOMASS,O%NPATCH))
  INL=MIN(O%NNBIOMASS,SIZE(ZW,2))
  DO JL=1,INL
     WHERE(ZW(:,JL,:)/=XUNDEF)
       R%XBIOMASS(:,JL,:) = ZW(:,JL,:)
     ELSEWHERE
       R%XBIOMASS(:,JL,:) = 0.0
     ENDWHERE
  ENDDO
  IF(O%NNBIOMASS>INL)THEN
    DO JL=INL+1,O%NNBIOMASS
       WHERE(ZW(:,JL,:)/=XUNDEF)
         R%XBIOMASS(:,JL,:) = ZW(:,INL,:)
       ELSEWHERE
         R%XBIOMASS(:,JL,:) = 0.0
       ENDWHERE
    ENDDO          
  ENDIF
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('LITTER') 
  ALLOCATE(R%XLITTER(KLAT,O%NNLITTER,O%NNLITTLEVS,O%NPATCH))
  DO JPATCH=1,O%NPATCH
    INL=0
    DO JJ=1,O%NNLITTER
       DO JL=1,O%NNLITTLEVS
          INL=INL+1
          WHERE(ZW(:,INL,JPATCH)/=XUNDEF)
             R%XLITTER(:,JJ,JL,JPATCH) = ZW(:,INL,JPATCH)
          ELSEWHERE
             R%XLITTER(:,JJ,JL,JPATCH) = 0.0
          ENDWHERE
       ENDDO
    ENDDO
  END DO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('SOILCARB') 
  ALLOCATE(R%XSOILCARB(KLAT,O%NNSOILCARB,O%NPATCH))
  WHERE(ZW(:,:,:)/=XUNDEF)
    R%XSOILCARB(:,:,:) = ZW(:,:,:)
  ELSEWHERE
    R%XSOILCARB(:,:,:) = 0.0
  ENDWHERE
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('LIGNIN') 
  ALLOCATE(R%XLIGNIN_STRUC(KLAT,O%NNLITTLEVS,O%NPATCH))
  WHERE(ZW(:,:,:)/=XUNDEF)
    R%XLIGNIN_STRUC(:,:,:) = ZW(:,:,:)
  ELSEWHERE
    R%XLIGNIN_STRUC(:,:,:) = 0.0
  ENDWHERE
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
END SELECT
!
DEALLOCATE(ZW)
!-------------------------------------------------------------------------------------
!
!*      8.     Deallocations
!
IF (ALLOCATED(ZFIELDOUTV)) DEALLOCATE(ZFIELDOUTV)
!
IF (LHOOK) CALL DR_HOOK('PREP_HOR_ISBA_CC_FIELD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!
!
END SUBROUTINE PREP_HOR_ISBA_CC_FIELD
