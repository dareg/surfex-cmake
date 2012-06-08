!     #########
      SUBROUTINE BUILD_PRONOSLIST_n(TPEMISS,TPPRONOS,KCH,KLUOUT,KVERB)
!!    #######################################################################
!!
!!*** *BUILD_PRONOSLIST*
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
!!
!!
!!    AUTHOR
!!    ------
!!    D. Gazen
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 01/02/00
!!    C. Mari  30/10/00 call to MODD_TYPE_EFUTIL
!!    D. Gazen 01/12/03 change emissions handling for surf. externalization
!!    P. Tulet 01/05/05 aerosols primary emission
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_OPEN_INPUTB
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_TYPE_EFUTIL
USE MODD_SV_n,  ONLY: NBEQ,  CSV, NSV_CHSBEG, NSV_CHSEND,NSV_AERBEG, NSV_AEREND
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!
TYPE(EMISSVAR_T), DIMENSION(:),INTENT(IN)  :: TPEMISS
TYPE(PRONOSVAR_T),             POINTER     :: TPPRONOS
INTEGER,                       INTENT(IN)  :: KCH     ! logical unit of input chemistry file
INTEGER,                       INTENT(IN)  :: KLUOUT  ! output listing channel
INTEGER,                       INTENT(IN)  :: KVERB   ! verbose level
!
!*       0.2  declaration of local variables
!
CHARACTER(LEN=256) :: YINPLINE ! input agregation line read from Namelist
INTEGER :: INDX     ! 
INTEGER :: INBCOEFF ! Numer of agregations coeff for one species
INTEGER :: JI       ! loop index
INTEGER :: INDX_PRO ! index of the pronostic variable in CNAMES array
INTEGER :: IERR
CHARACTER(LEN=32) :: YPRO_NAME, YEMIS_NAME ! Name of the pronostic & emission species
LOGICAL :: GFOUND
CHARACTER(LEN=6), DIMENSION(:),POINTER :: CNAMES
TYPE(PRONOSVAR_T),             POINTER :: HEAD,CURRENT
INTEGER :: IEQ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
!*    EXECUTABLE STATEMENTS
!     ---------------------
!
IF (LHOOK) CALL DR_HOOK('BUILD_PRONOSLIST_N',0,ZHOOK_HANDLE)
IF (KVERB >= 6) THEN
  WRITE(KLUOUT,*) '********     SUBROUTINE (CHIMIE): BUILD_PRONOSLIST     ********'
END IF
!
! CNAMES points on chemical variables name
IF (NSV_AEREND > 0) THEN
CNAMES => CSV(NSV_CHSBEG:NSV_AEREND)
IEQ = NSV_AEREND - NSV_CHSBEG + 1
ELSE
CNAMES => CSV(NSV_CHSBEG:NSV_CHSEND)
IEQ = NSV_CHSEND  - NSV_CHSBEG + 1
END IF

!
! Namelist is opened and the agregation eq. are reached
!
CALL CH_OPEN_INPUTB("AGREGATION", KCH , KLUOUT)
!
! Parse each eq. line and build the TPPRONOS list
!
NULLIFY(HEAD)
NULLIFY(CURRENT)
DO 
!
! Read a line and convert 'tab' to 'space' characters
! until the keyword 'END_AGREGATION' is reached
  READ(KCH,'(A)',IOSTAT=IERR) YINPLINE
  IF (IERR /= 0) EXIT
  YINPLINE = TRIM(ADJUSTL(YINPLINE))
  IF (LEN_TRIM(YINPLINE) == 0) CYCLE ! skip blank line
  IF (YINPLINE == 'END_AGREGATION') EXIT
  CALL TAB2SPACE(YINPLINE)
!
!
!Extract pronostic variable name
  INDX = INDEX(YINPLINE,' ')
  YPRO_NAME = YINPLINE(1:INDX-1)
  WRITE(KLUOUT,*) 'YPRO_NAME = ',YPRO_NAME
!
! search the variable in CNAMES, STOP if not FOUND
  GFOUND = .FALSE.
  DO JI=1,IEQ
    IF (CNAMES(JI) == YPRO_NAME) THEN 
      INDX_PRO = JI
      GFOUND = .TRUE.
      EXIT
    END IF
  END DO
  IF (.NOT. GFOUND) THEN
    WRITE(KLUOUT,*) 'BUILD_PRONOSLIST ERROR : ',TRIM(YPRO_NAME),&
            ' not found in pronostic variables list !'  
    CALL ABOR1_SFX('CH_BUILDPRONOSN: VARIABLE NOT FOUND')
  END IF
!
! If YPRO_NAME variable already encountered : append the new equation (coeffs)
  GFOUND = .FALSE.
  INBCOEFF = 0
  CURRENT=>HEAD
  DO WHILE(ASSOCIATED(CURRENT))
    IF (CURRENT%NAMINDEX == INDX_PRO) THEN
      INBCOEFF = CURRENT%NBCOEFF
      GFOUND   = .TRUE.
      EXIT
    END IF
    CURRENT=>CURRENT%NEXT
  END DO
  IF (.NOT. GFOUND) THEN
!   New pronostic cell is created
    ALLOCATE(CURRENT)
    CURRENT%NAMINDEX = INDX_PRO
    CURRENT%NEXT     => HEAD
    HEAD => CURRENT
  END IF
!
!
! Extract the agregation coeffs
  DO
! get REAL coeff
    YINPLINE = ADJUSTL(YINPLINE(INDX:))
    INDX = INDEX(YINPLINE,' ')
    IF (INDX == 1) EXIT
    INBCOEFF = INBCOEFF+1
    IF (INBCOEFF > JPNBCOEFFMAX) THEN
      WRITE(KLUOUT,*) 'FATAL ERROR : Number of aggregation coefficients for ',&
             TRIM(YPRO_NAME),' exceeds constant JPNBCOEFFMAX = ',JPNBCOEFFMAX  
      WRITE(KLUOUT,*) '=> You should increase the JPNBCOEFFMAX value in modd_type_efutil.f90'
      CALL ABOR1_SFX('CH_BUILDPRONOSN: NUMBER OF AGGREGATION COEFFICIENTS TOO BIG')
    END IF
    READ(YINPLINE(1:INDX-1),*) CURRENT%XCOEFF(INBCOEFF)
!
! get EMIS species name
    YINPLINE = ADJUSTL(YINPLINE(INDX:))
    INDX = INDEX(YINPLINE,' ')
    YEMIS_NAME = YINPLINE(1:INDX-1)
!
! check EMIS species name
    GFOUND = .FALSE.
    DO JI=1,SIZE(TPEMISS)
      IF (TPEMISS(JI)%CNAME == YEMIS_NAME) THEN
        GFOUND = .TRUE.
        CURRENT%NEFINDEX(INBCOEFF) = JI
        EXIT
      END IF
    END DO
    IF (.NOT. GFOUND) THEN
      WRITE(KLUOUT,*) 'ERROR : ',TRIM(YEMIS_NAME),&
              ' not found in emission variables list !'  
      CALL ABOR1_SFX('CH_BUILDPRONOSN: UNKNOWN EMISSION VARIABLE')
    END IF
  END DO
  CURRENT%NBCOEFF = INBCOEFF
END DO

!
! Update TPPRONOS pointer with head of list
TPPRONOS => HEAD
!
IF (KVERB >= 6) THEN
  WRITE(KLUOUT,*) 'BUILD_PRONOSLIST: Affichage resultats'
  CURRENT=>HEAD
  DO WHILE(ASSOCIATED(CURRENT))
    WRITE(KLUOUT,*) 'BUILD_PRONOSLIST: Species ',TRIM(CNAMES(CURRENT%NAMINDEX)),' (index ',&
            CURRENT%NAMINDEX,' in CNAMES)'  
    DO JI=1,CURRENT%NBCOEFF
      WRITE(KLUOUT,*) CURRENT%XCOEFF(JI),TPEMISS(CURRENT%NEFINDEX(JI))%CNAME
    END DO
    CURRENT=>CURRENT%NEXT
  END DO
  WRITE(KLUOUT,*) '******** END SUBROUTINE (CHIMIE): BUILD_PRONOSLIST     ********'
END IF

IF (LHOOK) CALL DR_HOOK('BUILD_PRONOSLIST_N',1,ZHOOK_HANDLE)
CONTAINS 
!!
!!    ###########################
      SUBROUTINE TAB2SPACE(HTEXT)
!!    ###########################
!!
!!*** *TAB2SPACE*
!!
!!    PURPOSE
!!    -------
!!     Convert 'tab' character to 'space' character in the string HTEXT
!!
!!**  METHOD
!!    ------
!!
!!    AUTHOR
!!    ------
!!    D. Gazen
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 01/02/2000
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!
CHARACTER(len=*),INTENT(INOUT) :: HTEXT
!
!*       0.2  declaration of local variables
!
CHARACTER, PARAMETER :: YPTAB = CHAR(9) ! TAB character is ASCII : 9
INTEGER              :: JI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
!*    EXECUTABLE STATEMENTS
!     ---------------------
!
IF (LHOOK) CALL DR_HOOK('TAB2SPACE',0,ZHOOK_HANDLE)
DO JI=1,LEN_TRIM(HTEXT)
  IF (HTEXT(JI:JI) == YPTAB) HTEXT(JI:JI) = ' '
END DO
IF (LHOOK) CALL DR_HOOK('TAB2SPACE',1,ZHOOK_HANDLE)
END SUBROUTINE TAB2SPACE

END SUBROUTINE BUILD_PRONOSLIST_n
