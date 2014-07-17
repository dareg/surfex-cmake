!     ###############################################################################
SUBROUTINE ASSIM_INLAND_WATER_n(HPROGRAM,KI,PTS_IN,PITM,HTEST, &
                                OLKEEPEXTZONE,OD_MASKEXT,PLON_IN,PLAT_IN)

!     ###############################################################################
!
!!****  *ASSIM_INLAND_WATER_n * - Chooses the surface assimilation schemes for INLAND_WATER parts  
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     T. Aspelien
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2012
!!      Trygve Aspelien, Separating IO  06/2013
!!--------------------------------------------------------------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_ASSIM,          ONLY : NPRINTLEV,LEXTRAP_WATER,LWATERTG2
!
USE MODD_SURF_ATM_n,     ONLY : CWATER,NR_WATER,NR_NATURE,NSIZE_NATURE
USE MODD_ISBA_n,         ONLY : XTG
USE MODD_WATFLUX_n,      ONLY : XTS,XZS
USE MODD_WATFLUX_GRID_n, ONLY : XLAT, XLON
USE YOMHOOK,             ONLY : LHOOK,DR_HOOK
USE PARKIND1,            ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_OI_HOR_EXTRAPOL_SURF
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN) :: HPROGRAM  ! program calling surf. schemes
INTEGER,            INTENT(IN) :: KI
REAL,DIMENSION(KI), INTENT(IN) :: PTS_IN
REAL,DIMENSION(KI), INTENT(IN) :: PITM
CHARACTER(LEN=2),   INTENT(IN) :: HTEST ! must be equal to 'OK'
LOGICAL, INTENT(IN) :: OLKEEPEXTZONE
LOGICAL, DIMENSION(KI), INTENT(IN) :: OD_MASKEXT
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLON_IN
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLAT_IN
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
REAL, DIMENSION(KI)              :: ZLST
REAL, DIMENSION(KI)              :: ZLST0
REAL, DIMENSION(KI)              :: ZLSTINC
REAL, DIMENSION(:), ALLOCATABLE  :: ZLST01, ZLST1, ZLON1, ZLAT1, ZALT1 
!
LOGICAL,DIMENSION(KI) :: GINTERP_LST
LOGICAL, DIMENSION(:), ALLOCATABLE :: GINTERP_LST1
INTEGER  :: IRESP,I,J,IS1,J1
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ASSIM_INLAND_WATER_N',0,ZHOOK_HANDLE)

IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_INLAND_WATER_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF

WRITE(*,*) 'UPDATING LST FOR INLAND_WATER: ',TRIM(CWATER)
!
!*     ZLST updated!
!
ZLST(:) = XUNDEF
IF (.NOT.LWATERTG2 ) THEN
  !*     ZLST updated from from CANARI analysis
  DO I=1,KI
    IF ( PITM(I)<0.5 ) ZLST(I) = PTS_IN(I)
  ENDDO
  !
ELSE
  ! Set TG2 from global array
  DO I=1,KI
    IF ( PITM(I)>0.5 ) THEN
      !*     ZLST updated from LAND values of climatological TS
      DO J=1,NSIZE_NATURE
        IF ( NR_WATER(I)==NR_NATURE(J) ) THEN
          ZLST(I) = XTG(J,2,1)
          EXIT
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  !
ENDIF
!
! Set local array from global
GINTERP_LST(:) = .FALSE.
DO I=1,KI
  IF ( ZLST(I)/=XUNDEF ) THEN
    ZLST0(I) = ZLST(I)
  ELSEIF ( LEXTRAP_WATER ) THEN
    ! Keep ZLST or do extrapolation from neighbour points
    ZLST0(I) = XUNDEF
    GINTERP_LST(I) = .TRUE.  
  ELSE
    ZLST0(I) = XTS(I)
  ENDIF
ENDDO
!
IF ( LEXTRAP_WATER ) THEN
  !
  IF (OLKEEPEXTZONE) THEN
    !     
    ZLST(:) = ZLST0(:)
    WHERE ( OD_MASKEXT(:) ) ZLST0(:) = XUNDEF
    CALL OI_HOR_EXTRAPOL_SURF(KI,PLAT_IN,PLON_IN,ZLST0,PLAT_IN,PLON_IN,ZLST,GINTERP_LST,XZS)
    !
  ELSE
    !
    IS1 = COUNT (.NOT.OD_MASKEXT)
    ALLOCATE (ZLST1(IS1), ZLST01(IS1), ZLAT1(IS1), ZLON1(IS1), ZALT1(IS1), GINTERP_LST1(IS1))
    !
    ! remove extension zone
    J = 1
    DO J1 = 1, KI
      IF ( .NOT.OD_MASKEXT(J1) )  THEN
        ZLST01(J) = ZLST0(J1)
        ZLAT1 (J) = PLAT_IN (J1)
        ZLON1 (J) = PLON_IN (J1)
        ZALT1 (J) = XZS  (J1)
        GINTERP_LST1(J) = GINTERP_LST(J1)
        J = J + 1
      ENDIF
    ENDDO
       
    ZLST1(:) = ZLST01(:)
    CALL OI_HOR_EXTRAPOL_SURF(IS1,ZLAT1,ZLON1,ZLST01,ZLAT1,ZLON1,ZLST1,GINTERP_LST1,ZALT1)
    !
    ! copy back
    J = 1
    DO J1 = 1, KI
      IF ( .NOT.OD_MASKEXT(J1) ) THEN
        ZLST(J1) = ZLST1(J)
        J = J + 1
      ENDIF
    ENDDO
    !
    DEALLOCATE (ZLST01, ZLST1, ZLAT1, ZLON1, ZALT1, GINTERP_LST1)
    !
  ENDIF
  !
ENDIF
!
!*     Print values produced by OI_HO_EXTRAPOL_SURF
IF ( NPRINTLEV > 2 ) THEN
  DO I=1,KI
    IF (GINTERP_LST(I)) THEN
      PRINT *,'Lake surface temperature set to ',ZLST(I),'from nearest neighbour at I=',NR_WATER(I)
    ENDIF
  ENDDO
ENDIF
!
! Sum the increments
ZLSTINC(:) = ZLST(:) - XTS(:)
WRITE(*,*) 'Mean LST increments over inland water   ',SUM(ZLSTINC)/KI
!
! Setting modified variables
XTS(:) = ZLST(:)
!
IF (LHOOK) CALL DR_HOOK('ASSIM_INLAND_WATER_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_INLAND_WATER_n
