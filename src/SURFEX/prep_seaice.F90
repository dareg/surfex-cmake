!     #########
SUBROUTINE PREP_SEAICE(HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!     #################################################################################
!
!!****  *PREP_SEAICE* - prepares variables for SEAICE scheme (for now : Gelato only)
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
!!     S. Sénési 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2014
!!------------------------------------------------------------------
!
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODI_GET_LUOUT
USE MODI_GET_TYPE_DIM_N
!
USE MODD_TYPES_GLT,   ONLY : T_GLT
!
USE MODN_PREP_SEAFLUX,   ONLY : CPREP_SEAICE_SCHEME => CSEAICE_SCHEME
USE MODD_SEAFLUX_n,      ONLY : CSEAICE_SCHEME,LHANDLE_SIC, TGLT, & 
                                XSSS, XSIC,                       &
                                LINTERPOL_SSS, CINTERPOL_SSS,     &
                                LINTERPOL_SIC, CINTERPOL_SIC,     &
                                LINTERPOL_SIT, CINTERPOL_SIT,     &
                                XSSS_MTH, XSIC_MTH, XSIT_MTH
USE MODI_PREP_HOR_SEAFLUX_FIELD
!
! Will be use later for interpolating input fields :
! USE MODD_SEAFLUX_GRID_n, ONLY : CGRID, XGRID_PAR, XLAT, XLON
!
USE MODD_GLT_PARAM, ONLY : nl, nt, nx, ny, nxglo, nyglo 
USE MODI_GLTOOLS_ALLOC
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
!
!*      0.2    declarations of local variables
!
INTEGER :: IK,IL  ! loop counter on ice categories and layers 
INTEGER :: JMTH,INMTH
INTEGER :: ILUOUT
LOGICAL :: GFOUND         ! Return code when searching namelist
INTEGER :: ILUNAM         ! logical unit of namelist file
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------------------
!
!*      0.     Default of configuration
!
!
IF (LHOOK) CALL DR_HOOK('PREP_SEAICE',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------------
!
!*      1.     Interpret namelist
!
CSEAICE_SCHEME=CPREP_SEAICE_SCHEME
!
LHANDLE_SIC = .FALSE.
IF(TRIM(CSEAICE_SCHEME)/='NONE' .OR. TRIM(CINTERPOL_SIC)/='NONE' )THEN
  LHANDLE_SIC=.TRUE.
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading and horizontal interpolations
!
!*      2.1    Salinity
!
CALL PREP_HOR_SEAFLUX_FIELD(HPROGRAM,'SSS    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.2    Seaice cover
!
IF (LHANDLE_SIC) THEN 
   CALL PREP_HOR_SEAFLUX_FIELD(HPROGRAM,'SIC    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
ENDIF
!
!
!
!*      4'     Optional preparation of interpolation of monthly Sea Surface salinity
!
LINTERPOL_SSS=.FALSE.
IF (TRIM(CSEAICE_SCHEME) /= 'NONE') THEN
   IF(TRIM(CINTERPOL_SSS)/='NONE')THEN
      LINTERPOL_SSS=.TRUE.
      !
      ! Precedent, Current and Next Monthly SSS
      INMTH=3
      ! Precedent, Current and Next Annual Monthly SSS
      IF(TRIM(CINTERPOL_SSS)=='ANNUAL')INMTH=14
      !
      ALLOCATE(XSSS_MTH(SIZE(XSSS),INMTH))
      DO JMTH=1,INMTH
         XSSS_MTH(:,JMTH)=XSSS(:)
      ENDDO
      !
   ENDIF
ENDIF
!-------------------------------------------------------------------------------------
!
!*      4''     Optional preparation of interpolation of monthly sea ice cover and sea 
!               ice thickness 
!
LINTERPOL_SIC=.FALSE.
IF(TRIM(CINTERPOL_SIC)/='NONE')THEN
   LINTERPOL_SIC=.TRUE.
ENDIF
!
IF(TRIM(CINTERPOL_SIT)/='NONE')THEN
   LINTERPOL_SIT=.TRUE.
ENDIF
!
IF(LINTERPOL_SIC)THEN
   !
   ! Precedent, Current and Next Monthly SIC
   INMTH=3
   ! Precedent, Current and Next Annual Monthly SIC
   IF(TRIM(CINTERPOL_SIC)=='ANNUAL')INMTH=14
   !
   ALLOCATE(XSIC_MTH(SIZE(XSSS),INMTH))
   DO JMTH=1,INMTH
      XSIC_MTH(:,JMTH)=XSIC(:)
   ENDDO
!
ENDIF
!
IF(LINTERPOL_SIT)THEN
   !
   ! Precedent, Current and Next Monthly SIT
   INMTH=3
   ! Precedent, Current and Next Annual Monthly SIT
   IF(TRIM(CINTERPOL_SIT)=='ANNUAL')INMTH=14
   !
   ALLOCATE(XSIT_MTH(SIZE(XSSS),INMTH))
   DO JMTH=1,INMTH
      XSIT_MTH(:,JMTH)=XUNDEF
   ENDDO
!
ENDIF
!-------------------------------------------------------------------------------------
!
!*      Creating default initial state for Gelato 
!
!WRITE(ILUOUT,*) ' NO FILE PROVIDED FOR SEAICE MODEL PROGNOSTIC FIELDS  -> Creating a default initial state for Gelato'
!
CALL GET_TYPE_DIM_n('SEA   ',nx)
ny=1
nyglo=1
nxglo=nx
CALL GLTOOLS_ALLOC(TGLT)
!
!*       G1    Prognostic fields with only space dimension(s) :
!
TGLT%ust(:,1)=0.
!
!*       G2     Prognostic fields with space and ice-category dimension(s) :
!
! sea ice age 
TGLT%sit(:,:,1)%age=0.
! melt pond volume 
TGLT%sit(:,:,1)%vmp=0.
! sea ice surface albedo 
TGLT%sit(:,:,1)%asn=0.
! sea ice fraction 
TGLT%sit(:,:,1)%fsi=0.
! sea ice thickness 
TGLT%sit(:,:,1)%hsi=1.*TGLT%sit(:,:,1)%fsi
! sea ice salinity 
TGLT%sit(:,:,1)%ssi=0.
! sea ice surface temperature 
TGLT%sit(:,:,1)%tsf=260.
! snow thickness 
TGLT%sit(:,:,1)%hsn=0.
! snow density 
TGLT%sit(:,:,1)%rsn=100.
!
!*       G3     Prognostic fields with space, ice-category and layer dimensions :
!
! sea ice vertical gltools_enthalpy profile for all types and levels
TGLT%sil(:,:,:,1)%ent=-1000. 
!
IF (LHOOK) CALL DR_HOOK('PREP_SEAICE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_SEAICE
