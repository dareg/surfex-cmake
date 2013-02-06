!##################
MODULE MODI_READ_SURF
!##################
!
!!****  *READX0* - routine to read a real scalar
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READX0 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
  INTERFACE READ_SURF
!
     SUBROUTINE READ_SURFX0(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC     ! name of the article to be read
REAL, INTENT(INOUT) :: PFIELD            ! real scalar to be read  
INTEGER,INTENT(OUT) :: KRESP             ! KRESP  : return-code if a problem appears 
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT  ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
END SUBROUTINE READ_SURFX0
!
     SUBROUTINE READ_SURFX1(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC     ! name of the article to be read
REAL, DIMENSION(:), INTENT(OUT) ::PFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP            ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
END SUBROUTINE READ_SURFX1
!
     SUBROUTINE READ_SURFX2(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM    ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC        ! name of the article to be read
REAL, DIMENSION(:,:), INTENT(OUT) :: PFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP               ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
END SUBROUTINE READ_SURFX2
!
      SUBROUTINE READ_SURFX2COV(HPROGRAM,HREC,PFIELD,OFLAG,KRESP,HCOMMENT,HDIR)
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM    ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC        ! name of the article to be read
REAL, DIMENSION(:,:), INTENT(OUT) :: PFIELD ! array containing the data field
LOGICAL,DIMENSION(:), INTENT(IN)  ::OFLAG   ! mask for array filling
INTEGER, INTENT(OUT) :: KRESP               ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
END SUBROUTINE READ_SURFX2COV
!
     SUBROUTINE READ_SURFX3(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM      ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC          ! name of the article to be read
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP                 ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
END SUBROUTINE READ_SURFX3
!
     SUBROUTINE READ_SURFN0(HPROGRAM,HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC     ! name of the article to be read
INTEGER, INTENT(OUT) :: KFIELD           ! integer to be read  
INTEGER, INTENT(OUT) :: KRESP            ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
END SUBROUTINE READ_SURFN0
!
     SUBROUTINE READ_SURFN1(HPROGRAM,HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM     ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC         ! name of the article to be read
INTEGER, DIMENSION(:), INTENT(OUT) :: KFIELD ! integer to be read  
INTEGER, INTENT(OUT) :: KRESP                ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
END SUBROUTINE READ_SURFN1
!
     SUBROUTINE READ_SURFC0(HPROGRAM,HREC,HFIELD,KRESP,HCOMMENT,HDIR)
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM   ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC       ! name of the article to be read
CHARACTER(LEN=*), INTENT(OUT) :: HFIELD    ! caracter to be read  
INTEGER, INTENT(OUT) :: KRESP              ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
END SUBROUTINE READ_SURFC0
!
      SUBROUTINE READ_SURFL0(HPROGRAM,HREC,OFIELD,KRESP,HCOMMENT,HDIR)
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC     ! name of the article to be read
LOGICAL, INTENT(OUT)         :: OFIELD   ! array containing the data field
INTEGER, INTENT(OUT) :: KRESP            ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
END SUBROUTINE READ_SURFL0
!
      SUBROUTINE READ_SURFL1(HPROGRAM,HREC,OFIELD,KRESP,HCOMMENT,HDIR)
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM     ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC         ! name of the article to be read
LOGICAL, DIMENSION(:), INTENT(OUT) :: OFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP                ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
END SUBROUTINE READ_SURFL1
!
      SUBROUTINE READ_SURFT0(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT,HDIR)
!
USE MODD_TYPE_DATE_SURF
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM  ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC      ! name of the article to be read
TYPE(DATE_TIME), INTENT(INOUT) ::TFIELD   ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP             ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
END SUBROUTINE READ_SURFT0
!
      SUBROUTINE READ_SURFT1(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT,HDIR)
!
USE MODD_TYPE_DATE_SURF
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC     ! name of the article to be read
TYPE (DATE_TIME), DIMENSION(:), INTENT(INOUT) :: TFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
END SUBROUTINE READ_SURFT1
!
      SUBROUTINE READ_SURFT2(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT,HDIR)
!
USE MODD_TYPE_DATE_SURF
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN) ::HREC      ! name of the article to be read
TYPE (DATE_TIME), DIMENSION(:,:), INTENT(INOUT)::TFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
END SUBROUTINE READ_SURFT2
!
END INTERFACE
!
END MODULE MODI_READ_SURF
!
!     #############################################################
      SUBROUTINE READ_SURFX0(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC, XTIME_NPIO_READ, XTIME_COMM_READ
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURF0_OL
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURF0_LFI
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURF0_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURF0_FA
#endif
#ifdef MNH
USE MODI_READ_SURFX0_MNH
#endif
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN)  :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN)  :: HREC     ! name of the article to be read
REAL, INTENT(OUT) :: PFIELD               ! the real scalar to be read  
INTEGER, INTENT(OUT) :: KRESP             ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
DOUBLE PRECISION   :: XTIME0
INTEGER            :: INFOMPI
REAL               :: ZWORK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX0',0,ZHOOK_HANDLE)
!
KRESP = 0
YCOMMENT = ""
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH
  CALL READ_SURFX0_MNH(YREC,ZWORK,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFX0_ARO(YREC,ZWORK,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='OFFLIN' .OR. HPROGRAM=='ASCII ' .OR. &
    HPROGRAM=='FA    ' .OR. HPROGRAM=='LFI   ' ) THEN 
  !
  IF (NRANK==NPIO) THEN
    ! 
#ifndef NOMPI    
    XTIME0 = MPI_WTIME()
#endif    
    !
!$OMP SINGLE
    !    
    IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL
      CALL READ_SURF0_OL(YREC,ZWORK,KRESP,YCOMMENT)
#endif
    ELSEIF (HPROGRAM=='LFI   ') THEN
#ifdef LFI
      CALL READ_SURF0_LFI(YREC,ZWORK,KRESP,YCOMMENT)
#endif
    ELSEIF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
      CALL READ_SURF0_ASC(YREC,ZWORK,KRESP,YCOMMENT)
#endif
    ELSEIF (HPROGRAM=='FA    ') THEN
#ifdef FA 
      CALL READ_SURF0_FA(YREC,ZWORK,KRESP,YCOMMENT)
#endif
    ENDIF
    !
!$OMP END SINGLE COPYPRIVATE(ZWORK,KRESP,YCOMMENT)
    !
#ifndef NOMPI    
    XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif    
    !    
    IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
    !    
  ENDIF
  !
#ifndef NOMPI
  IF (YDIR/='A' .AND. NPROC>1) THEN
    XTIME0 = MPI_WTIME()
!$OMP SINGLE   
    CALL MPI_BCAST(ZWORK,KIND(ZWORK)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
!$OMP END SINGLE COPYPRIVATE(ZWORK)
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
  ENDIF
#endif
  !
ENDIF
!
PFIELD=ZWORK
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX0',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX0
!
!     #############################################################
      SUBROUTINE READ_SURFX1(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFN_LFI
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFN_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFN_FA
#endif
#ifdef MNH 
USE MODI_READ_SURFX1_MNH
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM  ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC      ! name of the article to be read
REAL, DIMENSION(:), INTENT(OUT) :: PFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP             ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: IL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX1',0,ZHOOK_HANDLE)
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL = SIZE(PFIELD)
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFX1_MNH(YREC,IL,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='AROME ') THEN
#ifdef ARO
  CALL READ_SURFX1_ARO(YREC,IL,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL 
  CALL READ_SURFN_OL(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI
  CALL READ_SURFN_LFI(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURFN_ASC(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
  CALL READ_SURFN_FA(YREC,IL,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX1',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX1
!
!     #############################################################
      SUBROUTINE READ_SURFX2(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFN_LFI
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFN_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFN_FA
#endif
#ifdef MNH 
USE MODI_READ_SURFX2_MNH
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM    ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC        ! name of the article to be read
REAL, DIMENSION(:,:), INTENT(OUT) :: PFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP               ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: IL1, IL2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX2',0,ZHOOK_HANDLE)
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL1 = SIZE(PFIELD,1)
IL2 = SIZE(PFIELD,2)
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFX2_MNH(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFX2_ARO(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL
  CALL READ_SURFN_OL(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
  CALL READ_SURFN_LFI(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC
  CALL READ_SURFN_ASC(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
  CALL READ_SURFN_FA(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX2',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX2
!
!     #############################################################
      SUBROUTINE READ_SURFX2COV(HPROGRAM,HREC,PFIELD,OFLAG,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFN_LFI
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFN_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFN_FA
#endif
#ifdef MNH        
USE MODI_READ_SURFX2COV_MNH
#endif
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM    ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC        ! name of the article to be read
REAL, DIMENSION(:,:), INTENT(OUT) :: PFIELD ! array containing the data field
LOGICAL,DIMENSION(:), INTENT(IN) :: OFLAG   ! mask for array filling
INTEGER, INTENT(OUT) :: KRESP               ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: JCOVER
INTEGER            :: IL1, IL2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX2COV',0,ZHOOK_HANDLE)
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL1 = SIZE(PFIELD,1)
IL2 = SIZE(PFIELD,2)
!
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH        
    CALL READ_SURFX2COV_MNH(YREC,IL1,IL2,PFIELD,OFLAG,KRESP,YCOMMENT,YDIR)
#endif
ELSE
  !
  DO JCOVER=1,IL2
    !
    WRITE(YREC,'(A5,I3.3)') 'COVER',JCOVER
    YCOMMENT='X_Y_'//YREC
    PFIELD(:,JCOVER)=0.
    IF (.NOT. OFLAG(JCOVER)) CYCLE
    !
    IF (HPROGRAM=='AROME ') THEN
#ifdef ARO
      CALL READ_SURFX1_ARO(YREC,IL1,PFIELD(:,JCOVER),KRESP,YCOMMENT,YDIR)
#endif
    ENDIF
    !    
    IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL
      CALL READ_SURFN_OL(YREC,PFIELD(:,JCOVER),KRESP,YCOMMENT,YDIR)
#endif
    ENDIF
    !
    IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC        
      CALL READ_SURFN_ASC(YREC,PFIELD(:,JCOVER),KRESP,YCOMMENT,YDIR)
#endif
    ENDIF
    !
    IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
      CALL READ_SURFN_FA(YREC,IL1,PFIELD(:,JCOVER),KRESP,YCOMMENT,YDIR)
#endif
    ENDIF
    !
    IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
      CALL READ_SURFN_LFI(YREC,PFIELD(:,JCOVER),KRESP,YCOMMENT,YDIR)
#endif
    ENDIF
  END DO
  !
ENDIF
!
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX2COV',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX2COV
!
!     #############################################################
      SUBROUTINE READ_SURFX3(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM      ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC          ! name of the article to be read
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP                 ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: IL1, IL2, IL3
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX3',0,ZHOOK_HANDLE)
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL1 = SIZE(PFIELD,1)
IL2 = SIZE(PFIELD,2)
IL3 = SIZE(PFIELD,3)
!
!plmIF (HPROGRAM=='MESONH') THEN
!plm  CALL READ_SURFX3_MNH(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
!plmENDIF
!
IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL
  CALL READ_SURFN_OL(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
!plmIF (HPROGRAM=='ASCII ') THEN
!plm  CALL READ_SURFX3_ASC(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
!plmENDIF
!
!plmIF (HPROGRAM=='AROME ') THEN
!plm  CALL READ_SURFX3_ARO(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
!plmENDIF
!
!plmIF (HPROGRAM=='FA    ') THEN
!plm  CALL READ_SURFX3_FA(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
!plmENDIF
!
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX3',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX3
!
!     #############################################################
      SUBROUTINE READ_SURFN0(HPROGRAM,HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_SURFEX_OMP, ONLY : NBLOCK
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC, XTIME_NPIO_READ, XTIME_COMM_READ
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURF0_OL
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURF0_LFI
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURF0_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURF0_FA
#endif
#ifdef MNH
USE MODI_READ_SURFN0_MNH
#endif
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC     ! name of the article to be read
INTEGER, INTENT(OUT) :: KFIELD           ! the integer to be read  
INTEGER, INTENT(OUT) :: KRESP            ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
DOUBLE PRECISION   :: XTIME0
INTEGER            :: IWORK
INTEGER            :: INFOMPI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFN0',0,ZHOOK_HANDLE)
!
KRESP = 0
YCOMMENT = ""
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IWORK = 0
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH
  CALL READ_SURFN0_MNH(YREC,IWORK,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFN0_ARO(YREC,IWORK,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='OFFLIN' .OR. HPROGRAM=='ASCII ' .OR. &
    HPROGRAM=='FA    ' .OR. HPROGRAM=='LFI   ' ) THEN 
  !
  IF (NRANK==NPIO) THEN
    !
#ifndef NOMPI    
    XTIME0 = MPI_WTIME()
#endif
    !
!$OMP SINGLE
    !    
    IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL
      CALL READ_SURF0_OL(YREC,IWORK,KRESP,YCOMMENT)
#endif
    ENDIF
    !
    IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
      CALL READ_SURF0_LFI(YREC,IWORK,KRESP,YCOMMENT)
#endif
    ENDIF
    !
    IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
      CALL READ_SURF0_ASC(YREC,IWORK,KRESP,YCOMMENT)
#endif
    ENDIF
    !
    IF (HPROGRAM=='FA    ') THEN
#ifdef FA
      CALL READ_SURF0_FA(YREC,IWORK,KRESP,YCOMMENT)
#endif
    ENDIF
    !
!$OMP END SINGLE COPYPRIVATE(IWORK,KRESP,YCOMMENT)  
    !    
#ifndef NOMPI    
    XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
    !
    IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT  
    !
  ENDIF
  !
#ifndef NOMPI
  IF (YDIR/='A' .AND. NPROC>1) THEN          
    XTIME0 = MPI_WTIME()  
!$OMP SINGLE
    CALL MPI_BCAST(IWORK,KIND(IWORK)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
!$OMP END SINGLE COPYPRIVATE(IWORK)    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)    
  ENDIF    
#endif
  !
ENDIF
!
KFIELD=IWORK
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFN0',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFN0
!
!     #############################################################
      SUBROUTINE READ_SURFN1(HPROGRAM,HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFN_LFI
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFN_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFN_FA
#endif
#ifdef MNH 
USE MODI_READ_SURFN1_MNH
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM     ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC         ! name of the article to be read
INTEGER, DIMENSION(:), INTENT(OUT) :: KFIELD ! the integer to be read  
INTEGER, INTENT(OUT) :: KRESP                ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: IL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFN1',0,ZHOOK_HANDLE)
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL = SIZE(KFIELD,1)
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFN1_MNH(YREC,IL,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFN1_ARO(YREC,IL,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL
  CALL READ_SURFN_OL(YREC,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
  CALL READ_SURFN_LFI(YREC,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURFN_ASC(YREC,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='FA    ') THEN
#ifdef FA
  CALL READ_SURFN_FA(YREC,IL,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFN1',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFN1

!     #############################################################
      SUBROUTINE READ_SURFC0(HPROGRAM,HREC,HFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC, XTIME_NPIO_READ, XTIME_COMM_READ
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURF0_OL
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURF0_LFI
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURF0_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURF0_FA
#endif
#ifdef MNH 
USE MODI_READ_SURFC0_MNH
#endif
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN)  :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN)  :: HREC     ! name of the article to be read
CHARACTER(LEN=*), INTENT(OUT) :: HFIELD   ! the integer to be read  
INTEGER, INTENT(OUT) :: KRESP             ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL,INTENT(IN)  :: HDIR
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=40)  :: YFIELD
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
DOUBLE PRECISION   :: XTIME0
INTEGER            :: INFOMPI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFC0',0,ZHOOK_HANDLE)
!
KRESP = 0
YCOMMENT = ""
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFC0_MNH(YREC,YFIELD,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFC0_ARO(YREC,YFIELD,KRESP,YCOMMENT)
#endif
ELSEIF (HPROGRAM=='OFFLIN' .OR. HPROGRAM=='ASCII ' .OR. &
    HPROGRAM=='FA    ' .OR. HPROGRAM=='LFI   ' ) THEN 
  !
  IF (NRANK==NPIO) THEN
    !
#ifndef NOMPI    
    XTIME0 = MPI_WTIME()
#endif
    !  
!$OMP SINGLE    
    !    
    IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL
      CALL READ_SURF0_OL(YREC,YFIELD,KRESP,YCOMMENT)
#endif
    ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
      CALL READ_SURF0_LFI(YREC,YFIELD,KRESP,YCOMMENT)
#endif
    ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
      CALL READ_SURF0_ASC(YREC,YFIELD,KRESP,YCOMMENT)
#endif
    ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef FA
      CALL READ_SURF0_FA(YREC,YFIELD,KRESP,YCOMMENT)
#endif
    ENDIF
    !  
!$OMP END SINGLE COPYPRIVATE(YFIELD,KRESP,YCOMMENT)
    !
#ifndef NOMPI    
    XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)  
#endif    
    !    
    IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
    !    
  ENDIF
  !
#ifndef NOMPI
  IF (YDIR/='A' .AND. NPROC>1) THEN
    XTIME0 = MPI_WTIME()            
!$OMP SINGLE
    CALL MPI_BCAST(YFIELD,LEN(YFIELD),MPI_CHARACTER,NPIO,NCOMM,INFOMPI)
!$OMP END SINGLE COPYPRIVATE(YFIELD)    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
  ENDIF    
#endif
  !
ENDIF
!
HFIELD = YFIELD(1:LEN(HFIELD))
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFC0',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFC0
!
!     #############################################################
      SUBROUTINE READ_SURFL0(HPROGRAM,HREC,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC, XTIME_NPIO_READ, XTIME_COMM_READ
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURF0_OL
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURF0_LFI
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURF0_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURF0_FA
#endif
#ifdef MNH 
USE MODI_READ_SURFL0_MNH
#endif
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC     ! name of the article to be read
LOGICAL, INTENT(OUT) :: OFIELD           ! array containing the data field
INTEGER, INTENT(OUT) :: KRESP           ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT  ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
DOUBLE PRECISION   :: XTIME0
LOGICAL            :: GWORK
INTEGER            :: INFOMPI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFL0',0,ZHOOK_HANDLE)
!
KRESP = 0
YCOMMENT = ""
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFL0_MNH(YREC,GWORK,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFL0_ARO(YREC,GWORK,KRESP,YCOMMENT)
#endif
ELSEIF (HPROGRAM=='OFFLIN' .OR. HPROGRAM=='ASCII ' .OR. &
    HPROGRAM=='FA    ' .OR. HPROGRAM=='LFI   ' ) THEN 
  !  
  IF (NRANK==NPIO) THEN
    !
#ifndef NOMPI    
    XTIME0 = MPI_WTIME()
#endif    
    ! 
!$OMP SINGLE
    !    
    IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL
      CALL READ_SURF0_OL(YREC,GWORK,KRESP,YCOMMENT)
#endif
    ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
      CALL READ_SURF0_LFI(YREC,GWORK,KRESP,YCOMMENT)
#endif
    ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
      CALL READ_SURF0_ASC(YREC,GWORK,KRESP,YCOMMENT)
#endif
    ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
      CALL READ_SURF0_FA(YREC,GWORK,KRESP,YCOMMENT)
#endif
    ENDIF
    !
!$OMP END SINGLE COPYPRIVATE(GWORK,KRESP,YCOMMENT)
    ! 
#ifndef NOMPI    
    XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
    !    
    IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
    !    
  ENDIF
  !
#ifndef NOMPI
  IF (YDIR/='A' .AND. NPROC>1) THEN
    XTIME0 = MPI_WTIME()       
!$OMP SINGLE    
    CALL MPI_BCAST(GWORK,1,MPI_LOGICAL,NPIO,NCOMM,INFOMPI)
!$OMP END SINGLE COPYPRIVATE(GWORK)    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)    
  ENDIF
#endif 
  !
ENDIF
!
OFIELD=GWORK
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFL0',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFL0
!
!     #############################################################
      SUBROUTINE READ_SURFL1(HPROGRAM,HREC,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFN_LFI
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFN_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFN_FA
#endif
#ifdef MNH 
USE MODI_READ_SURFL1_MNH
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM     ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC         ! name of the article to be read
LOGICAL, DIMENSION(:), INTENT(OUT) :: OFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP                ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: IL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFL1',0,ZHOOK_HANDLE)
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL = SIZE(OFIELD)
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFL1_MNH(YREC,IL,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFL1_ARO(YREC,IL,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ELSE IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL
  CALL READ_SURFN_OL(YREC,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
  CALL READ_SURFN_LFI(YREC,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURFN_ASC(YREC,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
  CALL READ_SURFN_FA(YREC,IL,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFL1',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFL1
!
!     #############################################################
      SUBROUTINE READ_SURFT0(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_TYPE_DATE_SURF
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC, XTIME_NPIO_READ, XTIME_COMM_READ
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFT_OL
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFT_LFI
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFT_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFT_FA
#endif
#ifdef MNH 
USE MODI_READ_SURFT0_MNH
#endif
!
USE MODI_GET_LUOUT
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC     ! name of the article to be read
TYPE(DATE_TIME), INTENT(OUT) :: TFIELD   ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP            ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
!
DOUBLE PRECISION   :: XTIME0
REAL    :: ZTIME
INTEGER :: IDAY
INTEGER :: IMONTH
INTEGER :: IYEAR
INTEGER :: ILUOUT
INTEGER :: INFOMPI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT0',0,ZHOOK_HANDLE)
!
KRESP = 0
YCOMMENT = ""
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFT0_MNH(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFT0_ARO(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSEIF (HPROGRAM=='OFFLIN' .OR. HPROGRAM=='ASCII ' .OR. &
    HPROGRAM=='FA    ' .OR. HPROGRAM=='LFI   ' ) THEN 
  !  
  IF (NRANK==NPIO) THEN
    !
#ifndef NOMPI    
    XTIME0 = MPI_WTIME()
#endif    
    !
!$OMP SINGLE
    !    
    IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL
      CALL READ_SURFT_OL(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
    ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
      CALL READ_SURFT_LFI(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
    ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
      CALL READ_SURFT_ASC(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
    ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
      CALL READ_SURFT_FA(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
    ENDIF
    !
!$OMP END SINGLE COPYPRIVATE(IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)    
    ! 
#ifndef NOMPI    
    XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif    
    !    
    IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
    !    
  ENDIF
  !
#ifndef NOMPI
  IF (YDIR/='A' .AND. NPROC>1) THEN
    XTIME0 = MPI_WTIME() 
!$OMP SINGLE    
    CALL MPI_BCAST(IYEAR,KIND(IYEAR)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(IMONTH,KIND(IMONTH)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(IDAY,KIND(IDAY)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(ZTIME,KIND(ZTIME)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
!$OMP END SINGLE COPYPRIVATE(IYEAR,IMONTH,IDAY,ZTIME)    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)    
  ENDIF
#endif
  !
ENDIF
!
IF (KRESP==-2) THEN
  CALL GET_LUOUT(HPROGRAM,ILUOUT)
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) 'Date is not present file'
  WRITE(ILUOUT,*) 'Forcing value is kept'
  WRITE(ILUOUT,*) ' '
ELSE
  TFIELD%TDATE%YEAR = IYEAR
  TFIELD%TDATE%MONTH = IMONTH
  TFIELD%TDATE%DAY = IDAY
  TFIELD%TIME = ZTIME
END IF
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT0',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT0
!
!     #############################################################
      SUBROUTINE READ_SURFT1(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_TYPE_DATE_SURF
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC, XTIME_NPIO_READ, XTIME_COMM_READ
!
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFT_ASC
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFT_LFI
#endif
#ifdef MNH 
USE MODI_READ_SURFT1_MNH
#endif
!
USE MODI_ABOR1_SFX
USE MODI_GET_LUOUT
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM   ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC       ! name of the article to be read
TYPE(DATE_TIME), DIMENSION(:), INTENT(INOUT)::TFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP              ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: ILUOUT
INTEGER            :: INFOMPI
!
DOUBLE PRECISION   :: XTIME0
REAL , DIMENSION(SIZE(TFIELD,1))   :: ZTIME
INTEGER, DIMENSION(SIZE(TFIELD,1)) :: IDAY
INTEGER, DIMENSION(SIZE(TFIELD,1)) :: IMONTH
INTEGER, DIMENSION(SIZE(TFIELD,1)) :: IYEAR
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT1',0,ZHOOK_HANDLE)
!
KRESP = 0
YCOMMENT = ""
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFT1_MNH(YREC,IL1,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef ARO        
  CALL READ_SURFT1_ARO(YREC,IL1,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSEIF (HPROGRAM=='OFFLIN' .OR. HPROGRAM=='ASCII ' .OR. &
    HPROGRAM=='FA    ' .OR. HPROGRAM=='LFI   ' ) THEN 
  !
  IF (NRANK==NPIO) THEN
    !
#ifndef NOMPI    
    XTIME0 = MPI_WTIME()
#endif    
    !
!$OMP SINGLE
    !    
    IF (HPROGRAM=='OFFLIN') THEN
      CALL ABOR1_SFX('READ_SURFT1: NOT AVAILABLE FOR OFFLIN')
    ELSE IF (HPROGRAM=='FA    ') THEN
      CALL ABOR1_SFX('READ_SURFT1: NOT AVAILABLE FOR FA')      
    ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
      CALL READ_SURFT_ASC(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif 
    ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
      CALL READ_SURFT_LFI(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif  
    ENDIF
    !
#ifndef NOMPI    
    XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif    
    !      
!$OMP END SINGLE COPYPRIVATE(IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
    !  
  ENDIF
  !
#ifndef NOMPI
  IF (YDIR/='A' .AND. NPROC>1) THEN
    XTIME0 = MPI_WTIME()         
!$OMP SINGLE    
    CALL MPI_BCAST(IYEAR,SIZE(IYEAR)*KIND(IYEAR)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(IMONTH,SIZE(IMONTH)*KIND(IMONTH)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(IDAY,SIZE(IDAY)*KIND(IDAY)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(ZTIME,SIZE(ZTIME)*KIND(ZTIME)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
!$OMP END SINGLE COPYPRIVATE(IYEAR,IMONTH,IDAY,ZTIME)    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)    
  ENDIF
#endif
  !
ENDIF  
! 
IF (KRESP==-2) THEN
  CALL GET_LUOUT(HPROGRAM,ILUOUT)
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) 'Date is not present file'
  WRITE(ILUOUT,*) 'Forcing value is kept'
  WRITE(ILUOUT,*) ' '
ELSE
  TFIELD(:)%TDATE%YEAR  = IYEAR (:)
  TFIELD(:)%TDATE%MONTH = IMONTH(:)
  TFIELD(:)%TDATE%DAY   = IDAY  (:)  
  TFIELD(:)%TIME        = ZTIME (:)    
END IF
!
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT1',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT1
!
!     #############################################################
      SUBROUTINE READ_SURFT2(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC, XTIME_NPIO_READ, XTIME_COMM_READ
!
USE MODD_TYPE_DATE_SURF
!
USE MODI_ABOR1_SFX
USE MODI_GET_LUOUT
!
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFT_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFT_FA
#endif
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! calling program
CHARACTER(LEN=*), INTENT(IN) :: HREC     ! name of the article to be read
TYPE(DATE_TIME), DIMENSION(:,:), INTENT(INOUT) :: TFIELD ! array containing the data field  
INTEGER, INTENT(OUT) :: KRESP            ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=12)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: ILUOUT
INTEGER            :: INFOMPI
!
INTEGER :: IL1, IL2
DOUBLE PRECISION   :: XTIME0
REAL , DIMENSION(SIZE(TFIELD,1),SIZE(TFIELD,2))   :: ZTIME
INTEGER, DIMENSION(SIZE(TFIELD,1),SIZE(TFIELD,2)) :: IDAY
INTEGER, DIMENSION(SIZE(TFIELD,1),SIZE(TFIELD,2)) :: IMONTH
INTEGER, DIMENSION(SIZE(TFIELD,1),SIZE(TFIELD,2)) :: IYEAR
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT2',0,ZHOOK_HANDLE)
!
KRESP = 0
YCOMMENT = ""
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL1 = SIZE(TFIELD,1)
IL2 = SIZE(TFIELD,2)
!
IF (HPROGRAM=='MESONH') THEN
  CALL ABOR1_SFX('READ_SURFT2: NOT AVAILABLE FOR MESONH')
ELSE IF (HPROGRAM=='AROME ') THEN
  CALL ABOR1_SFX('READ_SURFT2: NOT AVAILABLE FOR AROME')  
ELSEIF (HPROGRAM=='OFFLIN' .OR. HPROGRAM=='ASCII ' .OR. &
    HPROGRAM=='FA    ' .OR. HPROGRAM=='LFI   ' ) THEN 
  !  
  IF (NRANK==NPIO) THEN
    !
#ifndef NOMPI   
    XTIME0 = MPI_WTIME()
#endif    
    !
!$OMP SINGLE
    !     
    IF (HPROGRAM=='OFFLIN') THEN
      CALL ABOR1_SFX('READ_SURFT2: NOT AVAILABLE FOR OFFLIN')
    ELSE IF (HPROGRAM=='LFI   ') THEN
      CALL ABOR1_SFX('READ_SURFT2: NOT AVAILABLE FOR LFI')      
    ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
      CALL READ_SURFT_ASC(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
    ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
      CALL READ_SURFT_FA(YREC,IL1,IL2,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
    ENDIF
    !
#ifndef NOMPI
    XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
    !    
!$OMP END SINGLE COPYPRIVATE(IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
    !    
  ENDIF
  !
#ifndef NOMPI
  IF (YDIR/='A' .AND. NPROC>1) THEN
    XTIME0 = MPI_WTIME()   
!$OMP SINGLE
    CALL MPI_BCAST(IYEAR,SIZE(IYEAR)*KIND(IYEAR)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(IMONTH,SIZE(IMONTH)*KIND(IMONTH)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(IDAY,SIZE(IDAY)*KIND(IDAY)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(ZTIME,SIZE(ZTIME)*KIND(ZTIME)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
!$OMP END SINGLE COPYPRIVATE(IYEAR,IMONTH,IDAY,ZTIME)    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
  ENDIF
#endif
  !
ENDIF
!
IF (KRESP==-2) THEN
  CALL GET_LUOUT(HPROGRAM,ILUOUT)
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) 'Date is not present file'
  WRITE(ILUOUT,*) 'Forcing value is kept'
  WRITE(ILUOUT,*) ' '
ELSE
  TFIELD(:,:)%TDATE%YEAR  = IYEAR (:,:)
  TFIELD(:,:)%TDATE%MONTH = IMONTH(:,:)
  TFIELD(:,:)%TDATE%DAY   = IDAY  (:,:)  
  TFIELD(:,:)%TIME        = ZTIME (:,:)  
END IF
!
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT2',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT2
