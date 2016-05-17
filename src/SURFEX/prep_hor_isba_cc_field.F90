!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE PREP_HOR_ISBA_CC_FIELD (DTCO, U, KLAT, IO, S, NK, NP, NPE,  &
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
USE MODD_ISBA_n, ONLY : ISBA_NK_t, ISBA_NP_t, ISBA_NPE_t, ISBA_K_t, ISBA_S_t, &
                        ISBA_PE_t, ISBA_P_t
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
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_NK_t), INTENT(INOUT) :: NK
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t), INTENT(INOUT) :: NPE
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
TYPE FOUT
  REAL, DIMENSION(:,:), ALLOCATABLE :: ZOUT
END TYPE FOUT
TYPE NFOUT
  TYPE(FOUT), DIMENSION(:), ALLOCATABLE :: AL
END TYPE NFOUT
TYPE (NFOUT) :: ZW
TYPE(NFOUT) :: ZF
!
TYPE(ISBA_P_t), POINTER :: PK
TYPE(ISBA_PE_t), POINTER :: PEK
!
 CHARACTER(LEN=6)              :: YFILETYPE ! type of input file
 CHARACTER(LEN=28)             :: YFILE     ! name of file
 CHARACTER(LEN=6)              :: YFILEPGDTYPE ! type of input file
 CHARACTER(LEN=28)             :: YFILEPGD     ! name of file
REAL, POINTER, DIMENSION(:,:,:)     :: ZFIELDIN  ! field to interpolate horizontally
REAL, POINTER, DIMENSION(:,:)       :: ZFIELD    ! field to interpolate horizontally
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZFIELDOUTP ! field interpolated   horizontally
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZFIELDOUTV ! field interpolated   horizontally
!
INTEGER                       :: ILUOUT    ! output listing logical unit
!
LOGICAL                       :: GUNIF     ! flag for prescribed uniform field
LOGICAL                       :: GPREP_AGS ! flag to prepare ags field (only external case implemeted)
!
INTEGER                       :: JPATCH    ! loop on patches
INTEGER                       :: JVEGTYPE  ! loop on vegtypes
INTEGER                       :: INL, INP, JJ, JL, JP ! Work integer
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
       LINTERP = (S%XVEGTYPE(:,JPATCH) > 0.)
    ELSEIF(INP==IO%NPATCH)THEN
       LINTERP = (S%XPATCH(:,JPATCH) > 0.)
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
ALLOCATE(ZW%AL(IO%NPATCH))
!
IF(GPREP_AGS)THEN
!
  DO JP = 1,IO%NPATCH
    PK => NP%AL(JP)
    !
    ALLOCATE(ZW%AL(JP)%ZOUT(PK%NSIZE_P,INL))
    !
    CALL VEGTYPE_GRID_TO_PATCH_GRID(JP,IO%NPATCH,PK%XVEGTYPE_PATCH,PK%XPATCH,&
                                    PK%NR_P,ZFIELDOUTV,ZW%AL(JP)%ZOUT)
  ENDDO
  !
  DEALLOCATE(ZFIELDOUTV)
  !
ELSE
!
  DO JP = 1,IO%NPATCH
    PK => NP%AL(JP)
    PEK => NPE%AL(JP)
    !
    SELECT CASE (HSURF)
      !
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  
      CASE('BIOMASS') 
        ALLOCATE(ZW%AL(JP)%ZOUT(PK%NSIZE_P,IO%NNBIOMASS))
        ZW%AL(JP)%ZOUT(:,:) = 0.
        WHERE(PEK%XLAI(:)/=XUNDEF)
          ZW%AL(JP)%ZOUT(:,1) = PEK%XLAI(:) * PK%XBSLAI_NITRO(:)
        ENDWHERE
        ZW%AL(JP)%ZOUT(:,2) = MAX( 0., (ZW%AL(JP)%ZOUT(:,1)/ (XCC_NIT/10.**XCA_NIT))  &
                                      **(1.0/(1.0-XCA_NIT)) - ZW%AL(JP)%ZOUT(:,1) )
      !
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !
      CASE('LITTER') 
        ALLOCATE(ZW%AL(JP)%ZOUT(PK%NSIZE_P,IO%NNLITTER*IO%NNLITTLEVS))
        ZW%AL(JP)%ZOUT(:,:) = 0.0
      !
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !
      CASE('SOILCARB') 
        ALLOCATE(ZW%AL(JP)%ZOUT(PK%NSIZE_P,IO%NNSOILCARB))
        ZW%AL(JP)%ZOUT(:,:) = 0.0
      !
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !
      CASE('LIGNIN') 
        ALLOCATE(ZW%AL(JP)%ZOUT(PK%NSIZE_P,IO%NNLITTLEVS))
        ZW%AL(JP)%ZOUT(:,:) = 0.0
      !
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !
    END SELECT
    !
  ENDDO
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
   DO JP = 1,IO%NPATCH
     PEK => NPE%AL(JP)
     PK => NP%AL(JP)

     ALLOCATE(PEK%XBIOMASS(PK%NSIZE_P,IO%NNBIOMASS))
     INL=MIN(IO%NNBIOMASS,SIZE(ZW%AL(JP)%ZOUT,2))
     DO JL=1,INL
       WHERE(ZW%AL(JP)%ZOUT(:,JL)/=XUNDEF)
         PEK%XBIOMASS(:,JL) = ZW%AL(JP)%ZOUT(:,JL)
       ELSEWHERE
         PEK%XBIOMASS(:,JL) = 0.0
       ENDWHERE
     ENDDO
     IF(IO%NNBIOMASS>INL)THEN
       DO JL=INL+1,IO%NNBIOMASS
         WHERE(ZW%AL(JP)%ZOUT(:,JL)/=XUNDEF)
           PEK%XBIOMASS(:,JL) = ZW%AL(JP)%ZOUT(:,INL)
         ELSEWHERE
           PEK%XBIOMASS(:,JL) = 0.0
         ENDWHERE
       ENDDO          
     ENDIF
   ENDDO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('LITTER')
   DO JP = 1,IO%NPATCH
     PEK => NPE%AL(JP)
     PK => NP%AL(JP)

     ALLOCATE(PEK%XLITTER(PK%NSIZE_P,IO%NNLITTER,IO%NNLITTLEVS))
     INL=0
     DO JJ=1,IO%NNLITTER
       DO JL=1,IO%NNLITTLEVS
         INL=INL+1
         WHERE(ZW%AL(JP)%ZOUT(:,INL)/=XUNDEF)
           PEK%XLITTER(:,JJ,JL) = ZW%AL(JP)%ZOUT(:,INL)
         ELSEWHERE
           PEK%XLITTER(:,JJ,JL) = 0.0
         ENDWHERE
       ENDDO
     ENDDO
   END DO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('SOILCARB') 
   DO JP = 1,IO%NPATCH
     PEK => NPE%AL(JP)
     PK => NP%AL(JP)

     ALLOCATE(PEK%XSOILCARB(PK%NSIZE_P,IO%NNSOILCARB))
     WHERE(ZW%AL(JP)%ZOUT(:,:)/=XUNDEF)
       PEK%XSOILCARB(:,:) = ZW%AL(JP)%ZOUT(:,:)
     ELSEWHERE
       PEK%XSOILCARB(:,:) = 0.0
     ENDWHERE
   ENDDO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('LIGNIN') 
   DO JP = 1,IO%NPATCH
     PEK => NPE%AL(JP)
     PK => NP%AL(JP)

     ALLOCATE(PEK%XLIGNIN_STRUC(PK%NSIZE_P,IO%NNLITTLEVS))
     WHERE(ZW%AL(JP)%ZOUT(:,:)/=XUNDEF)
       PEK%XLIGNIN_STRUC(:,:) = ZW%AL(JP)%ZOUT(:,:)
     ELSEWHERE
       PEK%XLIGNIN_STRUC(:,:) = 0.0
     ENDWHERE
   ENDDO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
END SELECT
!
DO JP = 1,IO%NPATCH
  DEALLOCATE(ZW%AL(JP)%ZOUT)
ENDDO
DEALLOCATE(ZW%AL)
!-------------------------------------------------------------------------------------
!
!*      8.     Deallocations
!
!
IF (LHOOK) CALL DR_HOOK('PREP_HOR_ISBA_CC_FIELD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!
!
END SUBROUTINE PREP_HOR_ISBA_CC_FIELD
