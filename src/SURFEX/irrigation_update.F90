!     #########
      SUBROUTINE IRRIGATION_UPDATE (AG, IMI, PTSTEP, KMONTH, KDAY, PTIME) 
!     ####################################################################
!
!!****  *IRRIGATION_UPDATE* - routine to update irrigation fields
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
!!      P. Le Moigne  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2006
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_IRRIG_t
USE MODD_AGRI_n, ONLY : AGRI_t
!
USE MODD_AGRI,   ONLY   : JPSTAGE, XTHRESHOLD
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE(ISBA_PARAM_IRRIG_t), INTENT(INOUT) :: IMI
TYPE(AGRI_t), INTENT(INOUT) :: AG
!
REAL,    INTENT(IN)  :: PTSTEP, PTIME
INTEGER, INTENT(IN)  :: KMONTH, KDAY
!
INTEGER              :: IL, JL                        
LOGICAL              :: GMASK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.1   Declarations of arguments
!-------------------------------------------------------------------------------
!
! Mask to realize update only once a day
!
IF (LHOOK) CALL DR_HOOK('MODI_IRRIGATION_UPDATE:IRRIGATION_UPDATE',0,ZHOOK_HANDLE)
GMASK = ( PTIME - PTSTEP < 0. ) .AND. ( PTIME >= 0. )
!
IF (GMASK) THEN

  WHERE( (IMI%XIRRIG(:,:).GT.0.).AND.(AG%LIRRIDAY(:,:)) .AND.(AG%NIRRINUM(:,:).LT.JPSTAGE))
    AG%NIRRINUM (:,:) = AG%NIRRINUM(:,:) + 1
    AG%LIRRIDAY (:,:) = .FALSE.
  ENDWHERE
!   
  DO IL=1,SIZE(IMI%XIRRIG,1)
    DO JL=1,SIZE(IMI%XIRRIG,2)
      AG%XTHRESHOLDSPT(IL,JL)=XTHRESHOLD(AG%NIRRINUM(IL,JL))
    ENDDO
  ENDDO
!
END IF
!
! Reinitialization of irrigation stage (necessary for runs from August to August)
!
IF((KMONTH==1).AND.(KDAY==1)) THEN
  AG%NIRRINUM(:,:) = 1
ENDIF
!
AG%LIRRIGATE(:,:) = .FALSE.
DO IL=1,SIZE(IMI%XIRRIG,1)
  DO JL=1,SIZE(IMI%XIRRIG,2)
    !
    ! Activate irrigation after seeding date
    !
    IF (KMONTH == IMI%TSEED(IL,JL)%TDATE%MONTH .AND. KDAY .GE. IMI%TSEED(IL,JL)%TDATE%DAY) THEN
      AG%LIRRIGATE(IL,JL) = .TRUE.
    END IF
    IF (KMONTH > IMI%TSEED(IL,JL)%TDATE%MONTH) THEN
      AG%LIRRIGATE(IL,JL) = .TRUE.
    END IF
    !
    ! Stop irrigation after reaping date
    !
    IF (KMONTH == IMI%TREAP(IL,JL)%TDATE%MONTH .AND. KDAY .GT. IMI%TREAP(IL,JL)%TDATE%DAY) THEN
      AG%LIRRIGATE(IL,JL) = .FALSE.
    END IF
    IF (KMONTH > IMI%TREAP(IL,JL)%TDATE%MONTH) THEN
      AG%LIRRIGATE(IL,JL) = .FALSE.
    END IF
  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('MODI_IRRIGATION_UPDATE:IRRIGATION_UPDATE',1,ZHOOK_HANDLE)
!
END SUBROUTINE IRRIGATION_UPDATE
