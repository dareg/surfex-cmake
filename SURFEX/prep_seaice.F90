!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE PREP_SEAICE (UG, DTCO, DTS, O, OR, KLAT, S, U, GCP, &
                        HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
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
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_t
USE MODD_OCEAN_n, ONLY : OCEAN_t
USE MODD_OCEAN_REL_n, ONLY : OCEAN_REL_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_GRID_CONF_PROJ_n, ONLY : GRID_CONF_PROJ_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODI_GET_TYPE_DIM_N
!
USE MODD_PREP,       ONLY : CINGRID_TYPE, CINTERP_TYPE
!
USE MODN_PREP_SEAFLUX,   ONLY : CPREP_SEAICE_SCHEME => CSEAICE_SCHEME
USE MODI_PREP_HOR_SEAFLUX_FIELD
USE MODI_READ_PREP_SEAFLUX_CONF
USE MODD_PREP_SEAFLUX,   ONLY : XSIC_UNIF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: DTS
TYPE(OCEAN_t), INTENT(INOUT) :: O
TYPE(OCEAN_REL_t), INTENT(INOUT) :: OR
INTEGER, INTENT(IN) :: KLAT
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(GRID_CONF_PROJ_t),INTENT(INOUT) :: GCP
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
!
!*      0.2    declarations of local variables
!
INTEGER :: JMTH,INMTH
INTEGER :: INP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_SEAICE',0,ZHOOK_HANDLE)
!
!*      0.     Default of configuration
!
!
!S%CSEAICE_SCHEME=CPREP_SEAICE_SCHEME ! ???

S%LHANDLE_SIC = .FALSE.
IF(TRIM(S%CSEAICE_SCHEME)/='NONE' .OR. TRIM(S%CINTERPOL_SIC)/='NONE' )THEN
  S%LHANDLE_SIC=.TRUE.
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading and horizontal interpolations of Seaice cover
!
IF (S%LHANDLE_SIC) THEN 
   CALL PREP_HOR_SEAFLUX_FIELD(DTCO, UG, U, GCP, DTS, O, OR, KLAT, S, &
                               HPROGRAM,'SIC      ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      3.     Optional preparation of interpolation of monthly sea ice cover and sea 
!              ice thickness 
!
S%LINTERPOL_SIC=.FALSE.
IF(TRIM(S%CINTERPOL_SIC)/='NONE')THEN
   S%LINTERPOL_SIC=.TRUE.
ENDIF
!
IF(TRIM(S%CINTERPOL_SIT)/='NONE')THEN
   S%LINTERPOL_SIT=.TRUE.
ENDIF
!
IF(S%LINTERPOL_SIC)THEN
   !
   ! Precedent, Current, Next, and Second-next Monthly SIC
   INMTH=4
   !
   ALLOCATE(S%XSIC_MTH(SIZE(S%XSIC),INMTH))
   DO JMTH=1,INMTH
      S%XSIC_MTH(:,JMTH)=S%XSIC(:)
   ENDDO
!
ENDIF
!
IF(S%LINTERPOL_SIT)THEN
   !
   !Precedent, Current, Next, and Second-next Monthly SIT
   INMTH=4
   !
   ALLOCATE(S%XSIT_MTH(SIZE(S%XSIC),INMTH))
   DO JMTH=1,INMTH
      S%XSIT_MTH(:,JMTH)=XUNDEF
   ENDDO
!
ENDIF

CALL GET_TYPE_DIM_n(DTCO, U, 'SEA   ', INP)
CALL S%ICE%PREP(DTCO, U, GCP, INP, KLAT, &
                HPROGRAM, HATMFILE, HATMFILETYPE, HPGDFILE, HPGDFILETYPE)
!
IF (LHOOK) CALL DR_HOOK('PREP_SEAICE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_SEAICE
