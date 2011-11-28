!     #########
      SUBROUTINE READ_PREP_WATFLUX_CONF(HPROGRAM,HVAR,HFILE,HFILETYPE,HATMFILE,HATMFILETYPE,KLUOUT,OUNIF)
!     #######################################################
!
!!****  *READ_PREP_WATFLUX_CONF* - routine to read the configuration for
!!                                 WATFLUX fields preparation
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      P. Le Moigne 10/2005, Phasage Arome
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_READ_PREP_SURF_ATM_CONF
!
USE MODN_PREP_WATFLUX
USE MODD_PREP_WATFLUX, ONLY : CFILE_WATFLX, CTYPE, XTS_WATER_UNIF
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling ISBA
CHARACTER(LEN=7),  INTENT(IN)  :: HVAR     ! variable treated
CHARACTER(LEN=28), INTENT(OUT) :: HFILE    ! file name
CHARACTER(LEN=6),  INTENT(OUT) :: HFILETYPE! file type
CHARACTER(LEN=28), INTENT(IN)  :: HATMFILE    ! atmospheric file name
CHARACTER(LEN=6),  INTENT(IN)  :: HATMFILETYPE! atmospheric file type
INTEGER,           INTENT(IN)  :: KLUOUT   ! logical unit of output listing
LOGICAL,           INTENT(OUT) :: OUNIF    ! flag for prescribed uniform field

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears 
                                    ! at the open of the file in LFI  routines 
INTEGER           :: ILUNAM         ! Logical unit of namelist file
!
CHARACTER(LEN=28) :: YNAMELIST      ! namelist file
!
LOGICAL           :: GFOUND         ! Return code when searching namelist
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('READ_PREP_WATFLUX_CONF',0,ZHOOK_HANDLE)
HFILE = '                         '
HFILETYPE = '      '
!
OUNIF     = .FALSE.
!
!-------------------------------------------------------------------------------
!
!* choice of input file
!  --------------------
!
IF (LEN_TRIM(HFILE)==0 .AND. LEN_TRIM(CFILE_WATFLX)>0 .AND. LEN_TRIM(CTYPE)>0) THEN
  HFILE     = CFILE_WATFLX
  HFILETYPE = CTYPE
END IF
!
!! If no file name in the scheme namelist,
!! try to find a name in NAM_SURF_ATM
!
IF (LEN_TRIM(HFILE)==0) THEN
!
CALL READ_PREP_SURF_ATM_CONF(HPROGRAM,HFILE,HFILETYPE,HATMFILE,HATMFILETYPE,KLUOUT)
!
END IF
!-------------------------------------------------------------------------------
!
!* Is an uniform field prescribed?
!  ------------------------------
!
    OUNIF = (XTS_WATER_UNIF/=XUNDEF) 
!
!-------------------------------------------------------------------------------
!
!* If no file and no uniform field is prescribed: error
!  ---------------------------------------------
!
IF (HVAR=='DATE   ' .OR. HVAR=='ZS     ' .AND. LHOOK) CALL DR_HOOK('READ_PREP_WATFLUX_CONF',1,ZHOOK_HANDLE)
IF (HVAR=='DATE   ' .OR. HVAR=='ZS     ') RETURN
!
IF (LEN_TRIM(HFILETYPE)==0 .AND. .NOT. OUNIF) THEN
   CALL ABOR1_SFX('READ_PREP_WATFLUX_CONF: AN INPUT VALUE IS REQUIRED FOR '//HVAR)
END IF
IF (LHOOK) CALL DR_HOOK('READ_PREP_WATFLUX_CONF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_PREP_WATFLUX_CONF
