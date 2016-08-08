!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE SFX_XIOS_DECLARE_FIELD(HREC, HDOMAIN, HAXIS, KLEV, HCOMMENT , KFREQOP) 
!!
!!
!!     PURPOSE
!!     --------
!!
!!     Declare field HREC and some attributes to XIOS if needed
!!     If 'units' or 'name' attribute is not defined using XIOS config 
!!        files , use HCOMMENT to declare it
!!     Same for domain and other axis, either using relevant args or
!!     with default values
!!  
!!     IMPLICIT ARGUMENTS :
!!     -------------------- 
!!
!!     EXTERNAL
!!     --------
!!
!!     XIOS LIBRARY
!!
!!     REFERENCE
!!     ---------
!!
!!     XIOS Reference guide - Yann Meurdesoif - 10/10/2014 - 
!!     svn co -r 515 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0 <dir> 
!!       cd <dir>/doc ; ....
!!
!!     AUTHOR
!!     ------
!!
!!     S.Sénési, CNRM
!!
!!     MODIFICATION
!!     --------------
!!
!!     Original    03/2016
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_XIOS,       ONLY     : LXIOS_DEF_CLOSED, YPATCH_DIM_NAME, NBASE_XIOS_FREQ
USE MODD_SURF_PAR,   ONLY     : XUNDEF

#ifdef WXIOS
USE XIOS
#endif
!
USE MODI_SET_AXIS
USE MODI_ABOR1_SFX
!
USE YOMHOOK, ONLY  : LHOOK,   DR_HOOK
USE PARKIND1, ONLY : JPRB, JPIM
!
IMPLICIT NONE

!
!   Arguments
!
 CHARACTER(LEN=*)   ,INTENT(IN)            :: HREC     ! field id
 CHARACTER(LEN=*)   ,INTENT(IN), OPTIONAL  :: HDOMAIN  ! name of the horiz domain
 CHARACTER(LEN=*)   ,INTENT(IN), OPTIONAL  :: HAXIS    ! name of the additional axis
INTEGER             ,INTENT(IN), OPTIONAL  :: KLEV     ! Axis size 
 CHARACTER(LEN=*)   ,INTENT(IN), OPTIONAL  :: HCOMMENT ! Comment string a la Surfex
INTEGER(KIND=JPIM)  ,INTENT(IN), OPTIONAL  :: KFREQOP  ! Sampling frequency, in minutes
!
 CHARACTER(1000)    :: YLDOMAIN
 CHARACTER(1000)    :: YLCOMMENT
 CHARACTER(1000)    :: YAXIS
 CHARACTER(3)       :: YIDIM
!
INTEGER(KIND=JPIM) :: IFREQOP  ! Sampling frequency, in minutes
INTEGER(KIND=JPIM) :: IIDIM, ILEV
LOGICAL            :: GISDEF, LVALID_AXIS, LLOOP
!
REAL(KIND=JPRB)    :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SFX_XIOS_DECLARE_FIELD',0,ZHOOK_HANDLE)
!
#ifdef WXIOS
!
! ----------------------------------------------------------------------
!   If XIOS init phase is over, just returns
! ----------------------------------------------------------------------
!
IF (LXIOS_DEF_CLOSED) THEN
   IF (LHOOK) CALL DR_HOOK('SFX_XIOS_DECLARE_FIELD_INTERNAL',1,ZHOOK_HANDLE)
   RETURN
ENDIF
!
YLDOMAIN='FULL'
IF (PRESENT(HDOMAIN)) YLDOMAIN=TRIM(HDOMAIN)
YLCOMMENT=''
IF (PRESENT(HCOMMENT)) YLCOMMENT=TRIM(HCOMMENT)
IFREQOP=0
IF (PRESENT(KFREQOP)) IFREQOP=KFREQOP
!
! HANDLE ADDITIONAL AXIS
!
!
IF (.NOT. PRESENT(HAXIS)) THEN
   !
   CALL SFX_XIOS_DECLARE_FIELD_INTERNAL(HREC, YLDOMAIN, YLCOMMENT, IFREQOP)
   !
ELSE 
   !
   LLOOP=.FALSE. ! By default, try to declare the field with its axis
   IF (TRIM(HAXIS)==TRIM(YPATCH_DIM_NAME)) THEN
      ! Enable the writing of the set of individual arrays 
      IF ( .NOT. PRESENT(KLEV)) THEN 
         CALL XIOS_GET_AXIS_ATTR(HAXIS, n_glo=ILEV)
      ELSE
         ILEV=KLEV
      ENDIF
      DO IIDIM=1,KLEV
         IF ( IIDIM < 10 ) THEN 
            WRITE(YIDIM,'(I1)') IIDIM
         ELSE
            WRITE(YIDIM,'(I2)') IIDIM
         ENDIF
         !write(0,*) '<field id="'//trim(HREC)//'_'//TRIM(YIDIM)//'", domain_ref="'//trim(YLDOMAIN)//'" />'
         CALL SFX_XIOS_DECLARE_FIELD_INTERNAL(TRIM(HREC)//'_'//TRIM(YIDIM), YLDOMAIN, YLCOMMENT, IFREQOP)
      END DO
      !
   ELSE 
      ! 
      CALL SFX_XIOS_DECLARE_FIELD_INTERNAL(HREC, YLDOMAIN, YLCOMMENT, IFREQOP)
      CALL XIOS_IS_DEFINED_FIELD_ATTR(HREC,axis_ref=GISDEF) 
      ! CALL XIOS_IS_DEFINED_FIELD_ATTR(HREC,grid_ref=LGRIDDEF) 
      !IF ( .NOT. GISDEF .AND. .NOT. LGRIDDEF ) THEN 
      IF ( .NOT. GISDEF ) THEN 
         IF ( TRIM(HAXIS) == '')  THEN
            LVALID_AXIS=.FALSE.
            YAXIS='dim_for_'//TRIM(HREC)
         ELSE
            LVALID_AXIS=XIOS_IS_VALID_AXIS(trim(HAXIS))
            YAXIS=TRIM(HAXIS)
         ENDIF
         IF (.NOT. LVALID_AXIS) THEN 
            IF ( PRESENT(KLEV)) THEN 
               CALL SET_AXIS(TRIM(YAXIS),KSIZE=KLEV)
            ELSE
               CALL ABOR1_SFX('SFX_XIOS_DECLARE_FIELD : MUST PROVIDE KLEV OR AN ALREADY DECLARED HAXIS for '//HREC)
            ENDIF
         ENDIF
         CALL XIOS_SET_FIELD_ATTR(HREC, axis_ref=TRIM(YAXIS))
      ENDIF
   ENDIF
   !
ENDIF
!
#endif

IF (LHOOK) CALL DR_HOOK('SFX_XIOS_DECLARE_FIELD',1,ZHOOK_HANDLE)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


CONTAINS 

SUBROUTINE SFX_XIOS_DECLARE_FIELD_INTERNAL(HREC, HDOMAIN, HCOMMENT , KFREQOP) 
USE MODD_XIOS,       ONLY     : COUTPUT_DEFAULT
USE MODD_SURF_PAR,   ONLY     : XUNDEF
USE YOMHOOK    , ONLY         : LHOOK,   DR_HOOK
USE PARKIND1   , ONLY         : JPRB, JPIM
!
#ifdef WXIOS
USE XIOS
#endif
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!   Arguments
!
CHARACTER(LEN=*)   ,INTENT(IN)  :: HREC     ! field id
CHARACTER(LEN=*)   ,INTENT(IN)  :: HDOMAIN  ! name of the horiz domain
CHARACTER(LEN=*)   ,INTENT(IN)  :: HCOMMENT ! Comment string a la Surfex
INTEGER(KIND=JPIM) ,INTENT(IN)  :: KFREQOP  ! Sampling frequency, in minutes
!
!  Local variables
!
LOGICAL            :: GISDEF
INTEGER            :: IPO,IPF
!
REAL(KIND=JPRB)    :: ZHOOK_HANDLE
!
#ifdef WXIOS
TYPE(xios_field)      :: field_hdl, other_field_hdl
TYPE(xios_fieldgroup) :: fieldgroup_hdl
TYPE(xios_file)       :: file_hdl
#endif
!
IF (LHOOK) CALL DR_HOOK('SFX_XIOS_DECLARE_FIELD',0,ZHOOK_HANDLE)
!
#ifdef WXIOS
!
!$OMP SINGLE
!
! ----------------------------------------------------------------------
!  We are still in the XIOS init phase =>  Define field if necessary 
! ----------------------------------------------------------------------
!
IF (.NOT. XIOS_IS_VALID_FIELD(HREC))  THEN

  CALL XIOS_GET_HANDLE("field_definition",fieldgroup_hdl)
  CALL XIOS_ADD_CHILD(fieldgroup_hdl,field_hdl,HREC)
  !IF (.NOT. XIOS_IS_VALID_FIELD("default_field")) &
  !    CALL ABOR1_SFX('sfx_xios_check_field:cannot output field '//HREC//' : no default_field is defined')
  CALL XIOS_SET_ATTR(field_hdl,name=HREC)
  !
  ! ----------------------------------------------------------------------
  ! If default_ouput file is defined and enabled, add this field to it
  ! ----------------------------------------------------------------------
  !
  IF ( XIOS_IS_VALID_FILE(COUTPUT_DEFAULT)) THEN 

    CALL XIOS_GET_HANDLE(COUTPUT_DEFAULT,file_hdl)
    !CALL XIOS_IS_DEFINED_FILE_ATTR(COUTPUT_DEFAULT,enabled=GISDEF) 
    !IF (GISDEF ) CALL XIOS_GET_FILE_ATTR(COUTPUT_DEFAULT,enabled=GISDEF)
    !IF (GISDEF) THEN 
    CALL XIOS_ADD_CHILD(file_hdl,field_hdl)
    CALL XIOS_SET_ATTR(field_hdl,field_ref=HREC)
    !ENDIF
    !ELSE
    !   CALL ABOR1_SFX('sfx_xios_check_field : cannot output field '//HREC//' : no default_output file is defined')

  ENDIF

ENDIF
!
! ----------------------------------------------------------------------
! If field enabling is not defined, set it to TRUE
! ----------------------------------------------------------------------
!
 CALL XIOS_IS_DEFINED_FIELD_ATTR(HREC,enabled=GISDEF) 
!
IF ( .NOT. GISDEF ) THEN
  CALL XIOS_SET_FIELD_ATTR(HREC, enabled=.TRUE.)
ENDIF
!
! ----------------------------------------------------------------------
! If field attribute 'domain' is not defined, set it
! ----------------------------------------------------------------------
!
!CALL XIOS_IS_DEFINED_FIELD_ATTR(HREC,grid_ref=LGRIDDEF) 
CALL XIOS_IS_DEFINED_FIELD_ATTR(HREC,domain_ref=GISDEF)
!
!IF ( .NOT. GISDEF .AND. .NOT. LGRIDDEF ) THEN 
IF ( .NOT. GISDEF ) THEN 
   IF (TRIM(HDOMAIN)=='') &
        CALL ABOR1_SFX('SFX_XIOS_DECLARE_FIELD_INTERNAL : MUST PROVIDE HDOMAIN '//HREC)
  CALL XIOS_SET_FIELD_ATTR(HREC, domain_ref=TRIM(HDOMAIN))
ENDIF
!
! ----------------------------------------------------------------------
! If Freq_op  is not defined , set it to the provided value (def : timestep)
! ----------------------------------------------------------------------
!
 CALL XIOS_IS_DEFINED_FIELD_ATTR(HREC,freq_op=GISDEF) 
IF ( .NOT. GISDEF ) THEN 
  IF (KFREQOP /= 0) THEN 
    ! Il faudrait tester la coherence avec NBASE_XIOS_FREQ*XIOS_TIMESTEP ...
    CALL XIOS_SET_FIELD_ATTR(HREC, freq_op=KFREQOP*XIOS_MINUTE)
  !ELSE
  !   CALL XIOS_SET_FIELD_ATTR(HREC, freq_op=XIOS_TIMESTEP)
  ENDIF
ENDIF
!
! ----------------------------------------------------------------------
! If Freq_offset  is not defined , set it to 0
! ----------------------------------------------------------------------
!
 CALL XIOS_IS_DEFINED_FIELD_ATTR(HREC,freq_offset=GISDEF) 
IF ( .NOT. GISDEF ) THEN 
  CALL XIOS_SET_FIELD_ATTR(HREC, freq_offset=0.*XIOS_SECOND)
ENDIF
!
! ----------------------------------------------------------------------
! If NetCDF variable name is not defined , set it
! ----------------------------------------------------------------------
!
 CALL XIOS_IS_DEFINED_FIELD_ATTR(HREC,name=GISDEF)
IF ( .NOT. GISDEF ) THEN
  CALL XIOS_SET_FIELD_ATTR(HREC, name=trim(HREC))
ENDIF
!
!
! ----------------------------------------------------------------------
! If operation is not defined , set it to "instant"
! ----------------------------------------------------------------------
!
 CALL XIOS_IS_DEFINED_FIELD_ATTR(HREC,operation=GISDEF) 
IF ( .NOT. GISDEF ) THEN 
  CALL XIOS_SET_FIELD_ATTR(HREC, operation="instant")
ENDIF
!
!
! ----------------------------------------------------------------------
! If field attribute 'long_name' is not defined or empty, set it 
! ----------------------------------------------------------------------
!
 CALL XIOS_IS_DEFINED_FIELD_ATTR(HREC,long_name=GISDEF) 
IF ( .NOT. GISDEF .AND. (TRIM(HCOMMENT) /= '') ) THEN 
  !CALL XIOS_SET_FIELD_ATTR(HREC,long_name=TRIM(HCOMMENT))
ENDIF
!
!
! ------------------------------------------------------------------------
! If field attribute 'unit' is not defined or empty, try to guess a value 
! from HCOMMENT (using rightmost string between parenthesis)
! ------------------------------------------------------------------------
!
 CALL XIOS_IS_DEFINED_FIELD_ATTR(HREC,unit=GISDEF)  
IF ( .NOT. GISDEF ) THEN 
  IPO = INDEX(HCOMMENT,"(",.TRUE.)
  IPF = INDEX(HCOMMENT,")",.TRUE.)
  IF ( (IPO > 0) .AND. (IPF>IPO+1) ) THEN
    CALL XIOS_SET_FIELD_ATTR(HREC,unit=HCOMMENT(IPO+1:IPF-1))
  ENDIF
ENDIF
!
! ----------------------------------------------------------------------
! Set default value to Surfex's one
! ----------------------------------------------------------------------
!
 CALL XIOS_SET_FIELD_ATTR(HREC,default_value=XUNDEF)
!
!$OMP END SINGLE
#endif
!
IF (LHOOK) CALL DR_HOOK('SFX_XIOS_DECLARE_FIELD_INTERNAL',1,ZHOOK_HANDLE)
! ----------------------------------------------------------------------
!
END SUBROUTINE SFX_XIOS_DECLARE_FIELD_INTERNAL





END SUBROUTINE SFX_XIOS_DECLARE_FIELD
