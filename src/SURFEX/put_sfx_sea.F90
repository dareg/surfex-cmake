!     #########
      SUBROUTINE PUT_SFX_SEA(KLUOUT,OCPL_SEAICE,OWATER,PSEA_SST,PSEA_UCU,        &
                             PSEA_VCU,PSEAICE_SIT,PSEAICE_CVR,PSEAICE_ALB )  
!     ####################################################
!
!!****  *PUT_SFX_SEA* - routine to put some variables from
!!                       an oceanic general circulation model
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
!!	B. Decharme      *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2009
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_CSTS,       ONLY : XTT, XTTS, XICEC
!
USE MODD_SURF_ATM_n, ONLY : NR_SEA, NSIZE_SEA, XSEA, &
                            NR_WATER, NSIZE_WATER,   &
                            XWATER, CWATER
!
USE MODI_PACK_SAME_RANK
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,           INTENT(IN)  :: KLUOUT
LOGICAL,           INTENT(IN)  :: OCPL_SEAICE
LOGICAL,           INTENT(IN)  :: OWATER
!
REAL, DIMENSION(:), INTENT(IN) :: PSEA_SST
REAL, DIMENSION(:), INTENT(IN) :: PSEA_UCU
REAL, DIMENSION(:), INTENT(IN) :: PSEA_VCU
REAL, DIMENSION(:), INTENT(IN) :: PSEAICE_SIT
REAL, DIMENSION(:), INTENT(IN) :: PSEAICE_CVR
REAL, DIMENSION(:), INTENT(IN) :: PSEAICE_ALB
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
CHARACTER(LEN=50)     :: YCOMMENT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_SEA',0,ZHOOK_HANDLE)
!
!*       1.0   Initialization
!              --------------
!
!
!*       2.0   Get variable over sea
!              ---------------------
!
IF(NSIZE_SEA>0)THEN
! 
  CALL TREAT_SEA(NSIZE_SEA)
!
ENDIF
!
!*       3.0   Get variable over water without flake
!              -------------------------------------
!
IF(OWATER.AND.NSIZE_WATER>0)THEN
!
  CALL TREAT_WATER(NSIZE_WATER)
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_SEA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE TREAT_SEA(KLU)
!
USE MODD_SEAFLUX_n,  ONLY : XSST,XTICE,XDIR_ALB,XSCA_ALB, &
                            XICE_ALB,XUMER,XVMER
USE MODI_PACK_SAME_RANK
!
IMPLICIT NONE
!
INTEGER,     INTENT(IN) :: KLU
!
REAL,    DIMENSION(KLU) :: ZSST     ! sea surface temperature
REAL,    DIMENSION(KLU) :: ZICE_FRAC! ice fraction
REAL                    :: ZTMIN    ! Minimum temperature over this proc
REAL                    :: ZTMAX    ! Maximum temperature over this proc
CHARACTER(LEN=50)       :: YCOMMENT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_SEA:TREAT_SEA',0,ZHOOK_HANDLE)
!
YCOMMENT='Sea surface temperature'
CALL PACK_SAME_RANK(NR_SEA(:),PSEA_SST(:),ZSST(:))
WHERE (ZSST(:)/=0.0) XSST(:)=ZSST(:)
CALL CHECK_SEA(YCOMMENT,XSST(:))
!
ZTMIN=MINVAL(XSST(:))
ZTMAX=MAXVAL(XSST(:))
!
IF(ZTMIN<=0.0.OR.ZTMAX>500.)THEN
  WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(KLUOUT,*)'SST from ocean model not define or not physic'
  WRITE(KLUOUT,*)'SST MIN =',ZTMIN,'SST MAX =',ZTMAX
  WRITE(KLUOUT,*)'There is certainly a problem between         '
  WRITE(KLUOUT,*)'SURFEX and OASIS sea/land mask               '
  WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  CALL ABOR1_SFX('PUT_SFX_SEA: SST from ocean model not define or not physic')
ENDIF
!
YCOMMENT='Sea u-current stress'
CALL PACK_SAME_RANK(NR_SEA(:),PSEA_UCU(:),XUMER(:))
CALL CHECK_SEA(YCOMMENT,XUMER(:))
!
YCOMMENT='Sea v-current stress'
CALL PACK_SAME_RANK(NR_SEA(:),PSEA_VCU(:),XVMER(:))
CALL CHECK_SEA(YCOMMENT,XVMER(:))
!
IF(OCPL_SEAICE)THEN
!
  YCOMMENT='Sea-ice Temperature'
  CALL PACK_SAME_RANK(NR_SEA(:),PSEAICE_SIT(:),XTICE(:))
  CALL CHECK_SEA(YCOMMENT,XTICE(:))
!
  YCOMMENT='Sea-ice cover'
  CALL PACK_SAME_RANK(NR_SEA(:),PSEAICE_CVR(:),ZICE_FRAC(:))
  CALL CHECK_SEA(YCOMMENT,ZICE_FRAC(:))
!
  WHERE(ZICE_FRAC(:)>=XICEC)
    XSST(:) = MIN(XSST(:),XTTS-0.01)
  ELSEWHERE
    XSST(:) = MAX(XSST(:),XTTS)
  ENDWHERE
!
  YCOMMENT='Sea-ice albedo'
  CALL PACK_SAME_RANK(NR_SEA(:),PSEAICE_ALB(:),XICE_ALB(:))
  CALL CHECK_SEA(YCOMMENT,XICE_ALB(:))
!
  XDIR_ALB(:)=XICE_ALB(:)
  XSCA_ALB(:)=XICE_ALB(:)
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_SEA:TREAT_SEA',1,ZHOOK_HANDLE)
!
END SUBROUTINE TREAT_SEA
!
!-------------------------------------------------------------------------------
!
SUBROUTINE TREAT_WATER(KLU)
!
USE MODD_WATFLUX_n,  ONLY : XTS,XTICE,XDIR_ALB,XSCA_ALB,XICE_ALB
!
IMPLICIT NONE
!
INTEGER,     INTENT(IN) :: KLU
!
REAL,    DIMENSION(KLU) :: ZICE_FRAC! ice fraction
REAL                    :: ZTMIN    ! Minimum temperature over this proc
REAL                    :: ZTMAX    ! Maximum temperature over this proc
CHARACTER(LEN=50)       :: YCOMMENT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_SEA:TREAT_WATER',0,ZHOOK_HANDLE)
!
YCOMMENT='Water surface temperature'
CALL PACK_SAME_RANK(NR_SEA(:),PSEA_SST(:),XTS(:))
CALL CHECK_SEA(YCOMMENT,XTS(:))
!
ZTMIN=MINVAL(XTS(:))
ZTMAX=MAXVAL(XTS(:))
!
IF(ZTMIN<=0.0.OR.ZTMAX>500.)THEN
  WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(KLUOUT,*)'TS_WATER from ocean model not define or not physic'
  WRITE(KLUOUT,*)'TS_WATER MIN =',ZTMIN,'TS_WATER MAX =',ZTMAX
  WRITE(KLUOUT,*)'There is certainly a problem between              '
  WRITE(KLUOUT,*)'SURFEX and OASIS sea/land mask                    '
  WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  CALL ABOR1_SFX('PUT_SFX_SEA: SST from ocean model not define or not physic')
ENDIF
!
YCOMMENT='Water-ice Temperature'
CALL PACK_SAME_RANK(NR_SEA(:),PSEAICE_SIT(:),XTICE(:))
CALL CHECK_SEA(YCOMMENT,XTICE(:))
!
YCOMMENT='Water-ice cover'
CALL PACK_SAME_RANK(NR_SEA(:),PSEAICE_CVR(:),ZICE_FRAC(:))
CALL CHECK_SEA(YCOMMENT,ZICE_FRAC(:))
!
WHERE(ZICE_FRAC(:)>=XICEC)
  XTS(:) = MIN(XTS(:),XTT-0.01)
ELSEWHERE
  XTS(:) = MAX(XTS(:),XTT)
ENDWHERE
!
YCOMMENT='Water-ice albedo'
CALL PACK_SAME_RANK(NR_SEA(:),PSEAICE_ALB(:),XICE_ALB(:))
CALL CHECK_SEA(YCOMMENT,XICE_ALB(:))
!
XDIR_ALB(:)=XICE_ALB(:)
XSCA_ALB(:)=XICE_ALB(:)
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_SEA:TREAT_WATER',1,ZHOOK_HANDLE)
!
END SUBROUTINE TREAT_WATER
!
!-------------------------------------------------------------------------------
!
SUBROUTINE CHECK_SEA(HCOMMENT,PFIELD)
!
IMPLICIT NONE
!
CHARACTER(LEN=*),   INTENT(IN) :: HCOMMENT
REAL, DIMENSION(:), INTENT(IN) :: PFIELD
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_SEA:CHECK_SEA',0,ZHOOK_HANDLE)
!
IF(ANY(PFIELD(:)>=XUNDEF))THEN
  WRITE(KLUOUT,*)'PUT_SFX_SEA: problem after get '//TRIM(HCOMMENT)//' from OASIS'
  WRITE(KLUOUT,*)'PUT_SFX_SEA: some points not defined = ',COUNT(PFIELD(:)>=XUNDEF)
  CALL ABOR1_SFX('PUT_SFX_SEA: problem after get '//TRIM(HCOMMENT)//' from OASIS')          
ENDIF
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_SEA:CHECK_SEA',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHECK_SEA
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PUT_SFX_SEA
