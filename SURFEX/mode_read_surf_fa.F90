!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
MODULE MODE_READ_SURF_FA
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURF_FA is
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
!!      J.Escobar      10/06/2013: replace DOUBLE PRECISION by REAL to handle problem for promotion of real on IBM SP
!----------------------------------------------------------------------------
!

USE MODN_IO_OFFLINE

INTERFACE READ_SURF0_FA
        MODULE PROCEDURE READ_SURFX0_FA
        MODULE PROCEDURE READ_SURFN0_FA
        MODULE PROCEDURE READ_SURFL0_FA
        MODULE PROCEDURE READ_SURFC0_FA
END INTERFACE
INTERFACE READ_SURFX_FA
        MODULE PROCEDURE READ_SURFX1_FA
        MODULE PROCEDURE READ_SURFX2_FA
END INTERFACE
INTERFACE READ_SURFN_FA
        MODULE PROCEDURE READ_SURFN1_FA
        MODULE PROCEDURE READ_SURFL1_FA
END INTERFACE
INTERFACE READ_SURFT_FA
        MODULE PROCEDURE READ_SURFT0_FA
        MODULE PROCEDURE READ_SURFT1_FA
        MODULE PROCEDURE READ_SURFT2_FA
END INTERFACE
!
!----------------------------------------------------------------------------
!
CONTAINS

      SUBROUTINE FACILE_NOEXTZ (KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, LDCOSP)
      IMPLICIT NONE
      INTEGER :: KREP, KNUMER, KNIVAU
      CHARACTER (LEN=*) :: CDPREF, CDSUFF
      LOGICAL :: LDCOSP
      REAL, TARGET :: PCHAMP (:)
      REAL, POINTER :: ZCHAMP (:)

      INTEGER :: ITYPTR, ITRONC, INLATI, INXLON, INIVER
      INTEGER, PARAMETER              :: MAX_LATS=5000
      INTEGER, DIMENSION (MAX_LATS)   :: INLOPA
      INTEGER, DIMENSION ((1+MAX_LATS)/2) :: INOZPA
      REAL, DIMENSION (MAX_LATS)      :: ZSINLA
      REAL, DIMENSION (2)     :: ZAHYBR, ZBHYBR
      REAL :: ZSLAPO, ZCLOPO, ZSLOPO, ZCODIL, ZREFER
      LOGICAL :: LLNOMM, LLERFA, LLIMST, LLMLAM, LLGARD
      CHARACTER (LEN=256) :: CLNOMF
      CHARACTER (LEN=16) :: CLSTTU, CLNOMC
      INTEGER :: INIMES
      INTEGER :: INXLON_NOEXTZ, INLATI_NOEXTZ
      INTEGER :: JLAT, JLON

      CALL FAIOPT (KREP,   KNUMER, LLNOMM, CLNOMF, CLSTTU,  &
                 & LLERFA, LLIMST, INIMES, CLNOMC)

      CALL FACIES (CLNOMC, ITYPTR, ZSLAPO, ZCLOPO, ZSLOPO, &
                 & ZCODIL, ITRONC, INLATI, INXLON, INLOPA, &
                 & INOZPA, ZSINLA, INIVER, ZREFER, ZAHYBR, &
                 & ZBHYBR, LLGARD)
      IF ( INLATI > MAX_LATS ) THEN
          CALL ABOR1('YOU MUST INCREASE MAX_LATS')
      ENDIF
      LLMLAM = ITYPTR < 0

      IF (LLMLAM) THEN
        ALLOCATE (ZCHAMP (INXLON * INLATI))
      ELSE
        ZCHAMP => PCHAMP
      ENDIF

      CALL FACILE (KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, ZCHAMP, LDCOSP)

      IF (LLMLAM) THEN

        INXLON_NOEXTZ = INLOPA (4)
        INLATI_NOEXTZ = INLOPA (6)

        DO JLON = 1, INXLON_NOEXTZ
          DO JLAT = 1, INLATI_NOEXTZ
            PCHAMP ((JLAT-1) * INXLON_NOEXTZ + JLON) = ZCHAMP ((JLAT-1) * INXLON + JLON)
          ENDDO
        ENDDO

        DEALLOCATE (ZCHAMP)

      ENDIF

      END SUBROUTINE FACILE_NOEXTZ
!
!----------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE SFX_FA_VERSION(ONEW)
!     #############################################################
!
!!****  *SFX_FA_VERSION* - routine to find which fa version 
!                          (convergence with GMAP var name)
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, CPREFIX1D
!
USE MODI_ERROR_READ_SURF_FA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
LOGICAL, INTENT(OUT) :: ONEW
!
LOGICAL :: GOLD, GWORK
INTEGER :: INGRIB, INBITS, ISTRON, IPUILA, IRESP
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:SFX_FA_VERSION',0,ZHOOK_HANDLE)
!
IRESP=0
!
CALL FANION(IRESP,NUNIT_FA,CPREFIX1D,0,'VERSION',ONEW,GWORK,INGRIB,INBITS,ISTRON,IPUILA)
IF (IRESP/=0) CALL ERROR_READ_SURF_FA('FULLVERSION',IRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:SFX_FA_VERSION',1,ZHOOK_HANDLE)
!
END SUBROUTINE SFX_FA_VERSION
!
!----------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE READ_SURFX0_FA (&
                                 HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READX0* - routine to read a real scalar
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NLUOUT, CMASK, CPREFIX1D
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF
USE MODI_ERROR_READ_SURF_FA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
!
!
CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
REAL,              INTENT(OUT) :: PFIELD   ! the real scalar to be read
INTEGER,           INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=50) :: YCOMMENT
CHARACTER(LEN=6)  :: YMASK
CHARACTER(LEN=18) :: YNAME ! Field Name
LOGICAL           :: GV8, GFOUND
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX0_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
CALL IO_BUFF(HREC,'R',GFOUND)
IF (.NOT. LFAGMAP) THEN
!
  CALL SFX_FA_VERSION(GV8)
  IF(GV8)THEN
    YNAME=CPREFIX1D//TRIM(HREC)
  ELSE
    YMASK=CMASK
    IF (GFOUND) YMASK='FULL  '
    YNAME=TRIM(YMASK)//TRIM(HREC)
  ENDIF
!
  CALL FALIT_R(KRESP,NUNIT_FA,YNAME,PFIELD)
  IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
ELSE
  CALL LFILEC (KRESP, NUNIT_FA, CPREFIX1D//TRIM(HREC), PFIELD, 1)
  IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
ENDIF
!
YCOMMENT = TRIM(YNAME)
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX0_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX0_FA
!
!----------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE READ_SURFX1_FA(HREC,KL,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX1* - routine to fill a real 1D array for the externalised surface 
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ, &
                            WLOG_MPI, NREQ, LSFX_MPI
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NLUOUT, NMASK, NFULL, NFULL_EXT, &
                            NDGL, NDLON, NDGUX, NDLUX, CPREFIX1D
!
USE MODE_FASURFEX
!
USE MODI_PACK_SAME_RANK
USE MODI_ERROR_READ_SURF_FA
USE MODI_READ_AND_SEND_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),    INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,             INTENT(IN)  :: KL       ! number of points
REAL, DIMENSION(:),  INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,             INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),  INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),    INTENT(IN)  :: HDIR     ! type of field :
                                             ! 'H' : field with
                                             !       horizontal spatial dim.
                                             ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=4)           :: YPREFIX
CHARACTER(LEN=3)           :: YPREF
LOGICAL                    :: GV8
!
INTEGER ::  I, J, INFOMPI
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE,NPROC-1) :: ISTATUS
#endif
!
REAL, DIMENSION(:), ALLOCATABLE :: ZWORK2, ZWORK
REAL   :: XTIME0
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX1_FA',0,ZHOOK_HANDLE)
!
!
KRESP=0
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
  ALLOCATE(ZWORK(NFULL))
  !
  YPREF=HREC(1:3)
  !
  IF (YPREF=='CLS' .OR. YPREF=='SUR' .OR. YPREF=='PRO' .OR. YPREF=='ATM') THEN
    ALLOCATE(ZWORK2(NFULL_EXT))
    CALL FACILE(KRESP,NUNIT_FA,HREC(1:4),0,HREC(5:16),ZWORK2,.FALSE.)
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
    DO J=1,NDGUX
      DO I=1,NDLUX
        ZWORK((J-1)*NDLUX + I) = ZWORK2((J-1)*NDLON + I)
      ENDDO
    ENDDO
    DEALLOCATE(ZWORK2)
    HCOMMENT = TRIM(HREC)
  ELSEIF (.NOT. LFAGMAP) THEN
    CALL SFX_FA_VERSION(GV8)
    IF(GV8)THEN
      YPREFIX=CPREFIX1D
    ELSE
      YPREFIX='S1D_'
    ENDIF
    CALL FACILE(KRESP,NUNIT_FA,YPREFIX,0,HREC,ZWORK,.FALSE.)
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)  
    HCOMMENT = YPREFIX//TRIM(HREC)
  ELSE
    CALL FACILE_NOEXTZ (KRESP, NUNIT_FA, CPREFIX1D, 0, HREC, ZWORK, .FALSE.)
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)  
    HCOMMENT = YPREFIX//TRIM(HREC)
  ENDIF
  !
ELSEIF (HDIR/='-') THEN
  ALLOCATE(ZWORK(0))
ENDIF
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
!
IF (HDIR=='E') THEN
  IF ( NRANK==NPIO ) THEN
    CALL PACK_SAME_RANK(NMASK,ZWORK(:),PFIELD(:))
  ENDIF
ELSEIF (HDIR=='A') THEN  ! no distribution on other tasks
  IF ( NRANK==NPIO ) THEN
#ifdef SFX_MPI
    IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
    PFIELD(:) = ZWORK(1:KL)
#ifdef SFX_MPI
    IF (LSFX_MPI) XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif
  ENDIF
ELSEIF (HDIR=='-') THEN ! distribution of the total field on other tasks
#ifdef SFX_MPI
  IF (LSFX_MPI .AND. NPROC>1) THEN
    XTIME0 = MPI_WTIME()
    CALL MPI_BCAST(NFULL,KIND(NFULL)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    IF ( NRANK/=NPIO ) ALLOCATE(ZWORK(NFULL))
    CALL MPI_BCAST(ZWORK(1:KL),KL*KIND(ZWORK)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
  ENDIF
#endif
  PFIELD(:) = ZWORK(1:KL)
ELSE
  CALL READ_AND_SEND_MPI(ZWORK,PFIELD,NMASK)
  !IF (NRANK==NPIO) THEN
  !  CALL MPI_WAITALL(NPROC-1,NREQ,ISTATUS,INFOMPI)
  !ENDIF  
ENDIF
!
DEALLOCATE(ZWORK)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX1_FA',1,ZHOOK_HANDLE)
! 
END SUBROUTINE READ_SURFX1_FA
!
!----------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE READ_SURFX2_FA(HREC,KL1,KL2,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX2* - routine to fill a real 2D array for the externalised surface 
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ, &
                            WLOG_MPI, NREQ, LSFX_MPI
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NLUOUT, NMASK, NFULL, CPREFIX2D
!
USE MODE_FASURFEX
!
USE MODI_PACK_SAME_RANK
USE MODI_ERROR_READ_SURF_FA
USE MODI_READ_AND_SEND_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),         INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                  INTENT(IN)  :: KL1      ! number of points
INTEGER,                  INTENT(IN)  :: KL2      ! 2nd dimension
REAL, DIMENSION(:,:),     INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),         INTENT(IN)  :: HDIR     ! type of field :
                                                  ! 'H' : field with
                                                  !       horizontal spatial dim.
                                                  ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=4)           :: YPREFIX
CHARACTER(LEN=2)           :: YPATCH
CHARACTER(LEN=3)           :: YNUM
LOGICAL                    :: GV8
!
INTEGER :: JL, I, INFOMPI ! loop counter
!
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE,NPROC-1) :: ISTATUS
#endif
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK2   ! work array read in the file
REAL:: XTIME0
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX2_FA',0,ZHOOK_HANDLE)
!
!
KRESP=0
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
  ALLOCATE(ZWORK2(NFULL,KL2))
  IF (.NOT. LFAGMAP) THEN 
  !  
    CALL SFX_FA_VERSION(GV8)
  !
    DO JL=1,KL2
      IF(GV8)THEN
        WRITE(YNUM,'(I3.3)')JL
        YPREFIX=CPREFIX2D//YNUM
      ELSE
        WRITE(YPATCH,'(I2.2)')JL
        YPREFIX='S'//YPATCH//'_'
      ENDIF
      CALL FACILE(KRESP,NUNIT_FA,YPREFIX,JL,HREC,ZWORK2(:,JL),.FALSE.)
      IF (KRESP/=0) THEN
        HCOMMENT = YPREFIX//TRIM(HREC)
        CALL ERROR_READ_SURF_FA(HCOMMENT,KRESP)
      ENDIF  
    END DO
  ELSE
    DO JL = 1, KL2
      CALL FACILE_NOEXTZ (KRESP, NUNIT_FA, CPREFIX2D, JL, HREC, ZWORK2(:,JL), .FALSE.)
      IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)  
    ENDDO
  ENDIF
  !
  HCOMMENT = 'PATCH_'//TRIM(HREC)
  !
ELSEIF (HDIR/='-') THEN
  ALLOCATE(ZWORK2(0,0))
ENDIF
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
!
IF (HDIR=='E') THEN
  IF ( NRANK==NPIO ) THEN
    CALL PACK_SAME_RANK(NMASK,ZWORK2(:,:),PFIELD(:,:))
  ENDIF
ELSEIF (HDIR=='A') THEN  ! no distribution on other tasks
  IF ( NRANK==NPIO ) THEN
#ifdef SFX_MPI
    IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
    PFIELD(:,:) = ZWORK2(1:KL1,1:KL2)
#ifdef SFX_MPI
    IF (LSFX_MPI) XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif
  ENDIF
ELSEIF (HDIR=='-') THEN ! distribution of the total field on other tasks
#ifdef SFX_MPI
  IF (LSFX_MPI) THEN
    IF (NPROC>1) THEN
      XTIME0 = MPI_WTIME()
      CALL MPI_BCAST(NFULL,KIND(NFULL)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
      IF ( NRANK/=NPIO ) ALLOCATE(ZWORK2(NFULL,KL2))    
      CALL MPI_BCAST(ZWORK2(1:KL1,1:KL2),KL1*KL2*KIND(ZWORK2)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
      XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
    ENDIF
    PFIELD(:,:) = ZWORK2(1:KL1,1:KL2)
  ENDIF
#endif
ELSE
  CALL READ_AND_SEND_MPI(ZWORK2,PFIELD,NMASK)
  !IF (NRANK==NPIO) THEN
  !  CALL MPI_WAITALL(NPROC-1,NREQ,ISTATUS,INFOMPI)
  !ENDIF  
ENDIF
!
DEALLOCATE(ZWORK2)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX2_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX2_FA
!
!----------------------------------------------------------------------------
!
!     #############################################################
      RECURSIVE SUBROUTINE READ_SURFN0_FA (HREC,KFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NLUOUT, NMASK, CMASK, CPREFIX1D
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF
USE MODI_ERROR_READ_SURF_FA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_IO_SURF_LFI, ONLY : LMNH_COMPATIBLE, NIU, NIB, NIE, NJU, NJB, NJE
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
!
!
CHARACTER(LEN=*),   INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,            INTENT(OUT) :: KFIELD   ! the integer to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100):: YCOMMENT
CHARACTER(LEN=6)  :: YMASK
CHARACTER(LEN=18) :: YNAME ! Field Name
LOGICAL           :: GV8, GFOUND
INTEGER(KIND=8)   :: IWORK8
CHARACTER(LEN=40) :: YGRID
INTEGER           :: IIMAX, IJMAX
!
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN0_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
CALL IO_BUFF(HREC,'R',GFOUND)
IF (.NOT. LFAGMAP) THEN
!
  CALL SFX_FA_VERSION(GV8)
  IF(GV8)THEN
    YNAME=CPREFIX1D//TRIM(HREC)
  ELSE
    YMASK=CMASK
    IF (GFOUND) YMASK='FULL  '
    YNAME=TRIM(YMASK)//TRIM(HREC)
  ENDIF
!
  CALL FALIT_I(KRESP,NUNIT_FA,YNAME,KFIELD)
  IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
ELSE
  CALL LFILEC (KRESP, NUNIT_FA, CPREFIX1D//TRIM(HREC), IWORK8, 1)
  IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
  KFIELD = IWORK8
  IF (HREC == 'DIM_FULL') THEN
    CALL READ_SURFC0_FA ("GRID_TYPE", YGRID, KRESP, HCOMMENT)
    LMNH_COMPATIBLE = (YGRID=="CARTESIAN " .OR. YGRID=="CONF PROJ ")
!
    IF (LMNH_COMPATIBLE) THEN
      CALL READ_SURFN0_FA ("IMAX", IIMAX, KRESP, HCOMMENT)
      NIU = IIMAX+2
      NIB = 2
      NIE = IIMAX+1
      CALL READ_SURFN0_FA ("JMAX", IJMAX, KRESP, HCOMMENT)
      NJU = IJMAX+2
      NJB = 2
      NJE = IJMAX+1
    END IF
  ENDIF
ENDIF
!
YCOMMENT = YNAME
HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN0_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFN0_FA
!
!----------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE READ_SURFN1_FA (&
                                 HREC,KL,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!
!
!
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ, & 
                            WLOG_MPI, LSFX_MPI
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NLUOUT, NMASK, NFULL, CMASK, CPREFIX1D
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF
USE MODI_ERROR_READ_SURF_FA
USE MODI_READ_AND_SEND_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
!
!
CHARACTER(LEN=*),       INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                INTENT(IN)  :: KL       ! number of points
INTEGER, DIMENSION(:),  INTENT(OUT) :: KFIELD   ! the integer to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
                                                ! 'H' : field with
                                                !       horizontal spatial dim.
                                                ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=6)  :: YMASK
CHARACTER(LEN=18) :: YNAME ! Field Name
LOGICAL           :: GV8, GFOUND
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
INTEGER ::  I, INFOMPI
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
!
REAL   :: XTIME0
INTEGER (KIND=8), ALLOCATABLE :: IWORK8 (:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN1_FA',0,ZHOOK_HANDLE)

!
KRESP = 0
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  ! 
  IF (.NOT. LFAGMAP) THEN
    CALL IO_BUFF(HREC,'R',GFOUND)
    !
    CALL SFX_FA_VERSION(GV8)
    IF(GV8)THEN
      YNAME=CPREFIX1D//TRIM(HREC)
    ELSE
      YMASK=CMASK
      IF (GFOUND) YMASK='FULL  '
      YNAME=TRIM(YMASK)//TRIM(HREC)
    ENDIF
    !
    IF (HDIR=='A') THEN
      ALLOCATE(IWORK(KL))
    ELSEIF (HDIR/='-') THEN
      ALLOCATE(IWORK(NFULL))
    END IF    
    !
    CALL FALIT_I_D(KRESP,NUNIT_FA,YNAME,SIZE(IWORK),IWORK)
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
    !
  ELSE
    ALLOCATE (IWORK(KL))
    ALLOCATE (IWORK8(KL))
    CALL LFILEC (KRESP, NUNIT_FA, CPREFIX1D//TRIM(HREC), IWORK8, KL)
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
    IWORK = IWORK8
  ENDIF
  !
  HCOMMENT = YNAME
  !
ELSEIF (HDIR/='-') THEN
  ALLOCATE(IWORK(0))
ELSE
  ALLOCATE(IWORK(KL))
  IF (LFAGMAP) ALLOCATE(IWORK8(KL))
ENDIF
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
!
IF (HDIR=='A') THEN  ! no distribution on other tasks
  IF ( NRANK==NPIO ) THEN
#ifdef SFX_MPI
    IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
    KFIELD(:) = IWORK(1:KL)
#ifdef SFX_MPI
    IF (LSFX_MPI) XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif
  ENDIF
ELSEIF (HDIR=='-') THEN ! distribution of the total field on other tasks
#ifdef SFX_MPI
  IF (LSFX_MPI) THEN
    IF (NPROC>1) THEN
      XTIME0 = MPI_WTIME()
      CALL MPI_BCAST(NFULL,KIND(NFULL)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI) 
      CALL MPI_BCAST(IWORK(1:KL),KL*KIND(IWORK)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
      XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
    ENDIF
  ENDIF
#endif
  KFIELD(:) = IWORK(1:KL)
ELSE
  CALL READ_AND_SEND_MPI(IWORK,KFIELD,NMASK)
ENDIF
!
DEALLOCATE(IWORK)
IF (LFAGMAP) DEALLOCATE(IWORK8)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN1_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFN1_FA
!
!----------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE READ_SURFC0_FA (&
                                 HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READC0* - routine to read a character
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NLUOUT, CMASK, CPREFIX1D
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF
USE MODI_ERROR_READ_SURF_FA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
!
!
CHARACTER(LEN=*),   INTENT(IN)  :: HREC      ! name of the article to be read
CHARACTER(LEN=40),  INTENT(OUT) :: HFIELD    ! the integer to be read
INTEGER,            INTENT(OUT) :: KRESP     ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT  ! comment
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=50)       :: YCOMMENT
CHARACTER(LEN=6)        :: YMASK
CHARACTER(LEN=18)       :: YNAME ! Field Name
CHARACTER,DIMENSION(40) :: YFIELD
LOGICAL                 :: GV8, GFOUND
INTEGER (KIND=8), ALLOCATABLE :: IWORK8 (:)
INTEGER                 :: J, N
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFC0_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
CALL IO_BUFF(HREC,'R',GFOUND)
!
IF (.NOT. LFAGMAP) THEN
  CALL SFX_FA_VERSION(GV8)
  IF(GV8)THEN
    YNAME=CPREFIX1D//TRIM(HREC)
  ELSE
    YMASK=CMASK
    IF (GFOUND) YMASK='FULL  '
    YNAME=TRIM(YMASK)//TRIM(HREC)
  ENDIF
  !
  CALL FALIT_C(KRESP,NUNIT_FA,YNAME,40,YFIELD)
  IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
ELSE
  N = 40
  ALLOCATE (IWORK8(N))
  CALL LFILEC (KRESP, NUNIT_FA, CPREFIX1D//TRIM(HREC), IWORK8, SIZE(IWORK8))
  IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
  DO J = 1, N
    YFIELD (J:J) = CHAR (IWORK8(J))
  ENDDO
  DEALLOCATE(IWORK8)
ENDIF
WRITE(HFIELD,'(40A1)') YFIELD(:)
!
YCOMMENT = YNAME
HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFC0_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFC0_FA
!
!
!     #############################################################
      SUBROUTINE READ_SURFL0_FA (&
                                 HREC,OFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READL0* - routine to read a logical
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NLUOUT, CMASK, CPREFIX1D
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF
USE MODI_ERROR_READ_SURF_FA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
!
!
CHARACTER(LEN=*),   INTENT(IN)  :: HREC     ! name of the article to be read
LOGICAL,            INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=50) :: YCOMMENT
CHARACTER(LEN=6)  :: YMASK
CHARACTER(LEN=18) :: YNAME ! Field Name
LOGICAL           :: GV8, GFOUND
INTEGER (KIND=8)  :: IWORK8
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL0_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
CALL IO_BUFF(HREC,'R',GFOUND)
IF (.NOT. LFAGMAP) THEN
!
  CALL SFX_FA_VERSION(GV8)
  IF(GV8)THEN
    YNAME=CPREFIX1D//TRIM(HREC)
  ELSE
    YMASK=CMASK
    IF (GFOUND) YMASK='FULL  '
    YNAME=TRIM(YMASK)//TRIM(HREC)
  ENDIF
!
  CALL FALIT_L(KRESP,NUNIT_FA,YNAME,OFIELD)
  IF (KRESP/=0)CALL ERROR_READ_SURF_FA(HREC,KRESP)
ELSE
  CALL LFILEC (KRESP, NUNIT_FA, CPREFIX1D//TRIM(HREC), IWORK8, 1)
  IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
  OFIELD = IWORK8 /= 0
ENDIF
!
YCOMMENT = YNAME
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL0_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFL0_FA
!
!
!     #############################################################
      SUBROUTINE READ_SURFL1_FA (&
                                 HREC,KL,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READL1* - routine to read a logical array
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ, &
                            WLOG_MPI, LSFX_MPI
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NLUOUT, CMASK, CPREFIX1D
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF
USE MODI_ERROR_READ_SURF_FA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
!
!
CHARACTER(LEN=*),       INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                INTENT(IN)  :: KL       ! number of points
LOGICAL, DIMENSION(:),  INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
                                                ! 'H' : field with
                                                !       horizontal spatial dim.
                                                ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=6)  :: YMASK
CHARACTER(LEN=18) :: YNAME ! Field Name
LOGICAL           :: GV8, GFOUND
!
INTEGER           :: INFOMPI
INTEGER (KIND=8), ALLOCATABLE :: IWORK8 (:)
REAL  :: XTIME0
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL1_FA',0,ZHOOK_HANDLE)
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
!
KRESP = 0
!
IF (NRANK==NPIO) THEN
  !
  ! 
  CALL IO_BUFF(HREC,'R',GFOUND)
  !
  IF (.NOT. LFAGMAP) THEN
    CALL SFX_FA_VERSION(GV8)
    IF(GV8)THEN
      YNAME=CPREFIX1D//TRIM(HREC)
    ELSE
      YMASK=CMASK
      IF (GFOUND) YMASK='FULL  '
      YNAME=TRIM(YMASK)//TRIM(HREC)
    ENDIF
    !
    CALL FALIT_L_D(KRESP,NUNIT_FA,YNAME,KL,OFIELD)
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
  ELSE
    ALLOCATE (IWORK8(KL))
    CALL LFILEC (KRESP, NUNIT_FA, CPREFIX1D//TRIM(HREC), IWORK8, KL)
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
    OFIELD = IWORK8 /= 0
    DEALLOCATE(IWORK8)
  ENDIF
  !
  HCOMMENT = YNAME
  !
  !
ENDIF
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
!
#ifdef SFX_MPI
IF (LSFX_MPI .AND. NPROC>1 .AND. HDIR/='A') THEN
  XTIME0 = MPI_WTIME()
  CALL MPI_BCAST(OFIELD,KL,MPI_LOGICAL,NPIO,NCOMM,INFOMPI)
  XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
ENDIF
#endif
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL1_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFL1_FA
!
!----------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE READ_SURFT0_FA (&
                                 HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a date
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NLUOUT, CMASK, CPREFIX1D
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF
USE MODI_ERROR_READ_SURF_FA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB, JPRD
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
!
!
CHARACTER(LEN=*),  INTENT(IN)   :: HREC     ! name of the article to be read
INTEGER,            INTENT(OUT) :: KYEAR    ! year
INTEGER,            INTENT(OUT) :: KMONTH   ! month
INTEGER,            INTENT(OUT) :: KDAY     ! day
REAL,               INTENT(OUT) :: PTIME    ! year
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment

!*      0.2   Declarations of local variables
!
CHARACTER(LEN=50) :: YCOMMENT
CHARACTER(LEN=6)  :: YMASK
CHARACTER(LEN=18) :: YNAME ! Field Name
LOGICAL           :: GV8, GFOUND
!
INTEGER, DIMENSION(3) :: ITDATE
INTEGER (KIND=8) :: IWORK8(4)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT0_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
CALL IO_BUFF(HREC,'R',GFOUND)
IF (.NOT. LFAGMAP) THEN
!
  CALL SFX_FA_VERSION(GV8)
  IF(GV8)THEN
    YNAME=CPREFIX1D//TRIM(HREC)//'%TDATE'
  ELSE
    YMASK=CMASK
    IF (GFOUND) YMASK='FULL  '
    YNAME=TRIM(YMASK)//TRIM(HREC)//'%TDATE'
  ENDIF
!
  CALL FALIT_I_D(KRESP,NUNIT_FA,YNAME,3,ITDATE)
  IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
!
  KYEAR  = ITDATE(1)
  KMONTH = ITDATE(2)
  KDAY   = ITDATE(3)
!
  CALL SFX_FA_VERSION(GV8)
  IF(GV8)THEN
    YNAME=CPREFIX1D//TRIM(HREC)//'%TIME'
  ELSE
    YMASK=CMASK
    IF (GFOUND) YMASK='FULL  '
    YNAME=TRIM(YMASK)//TRIM(HREC)//'%TIME'
  ENDIF
!
  CALL FALIT_R(KRESP,NUNIT_FA,YNAME,PTIME)
  IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
ELSE
  CALL LFILEC (KRESP, NUNIT_FA, TRIM(CPREFIX1D)//HREC, IWORK8, SIZE (IWORK8))
  IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
  KYEAR  = INT (IWORK8(1))
  KMONTH = INT (IWORK8(2))
  KDAY   = INT (IWORK8(3))
  PTIME  = REAL(TRANSFER(IWORK8(4),1.0_JPRD),JPRB)
ENDIF
!
YCOMMENT = TRIM(HREC)
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT0_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT0_FA
!----------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE READ_SURFT1_FA (&
                                 HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READN0* - routine to read a 1D array of dates
!
!
!
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, XTIME_NPIO_READ, LSFX_MPI
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NLUOUT, CPREFIX1D
!
!USE MODE_FASURFEX
!
!USE MODI_IO_BUFF
USE MODI_ERROR_READ_SURF_FA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
!
!
CHARACTER(LEN=*),       INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:),  INTENT(OUT) :: KYEAR    ! year
INTEGER, DIMENSION(:),  INTENT(OUT) :: KMONTH   ! month
INTEGER, DIMENSION(:),  INTENT(OUT) :: KDAY     ! day
REAL,    DIMENSION(:),  INTENT(OUT) :: PTIME    ! year
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=50) :: YCOMMENT
CHARACTER(LEN=18) :: YNAME ! Field Name
!
REAL   :: XTIME0
INTEGER :: J
INTEGER (KIND=8), ALLOCATABLE :: IWORK (:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN1_FA',0,ZHOOK_HANDLE)

!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
  IF (.NOT. LFAGMAP) THEN
    WRITE(NLUOUT,*) ' READ_SURFT1_FA : time in 1 dimension with LFAGMAP=.T. not yet implemented'
    CALL ABOR1_SFX('MODE_READ_SURF_FA:READ_SURFT1_FA: time in 1 dimension with LFAGMAP=.T. not yet implemented')
  ELSE
    ALLOCATE (IWORK (SIZE(KYEAR)*4))
    CALL LFILEC (KRESP, NUNIT_FA, CPREFIX1D//TRIM(HREC), IWORK, SIZE(IWORK))
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
    DO J = 1, SIZE(KYEAR)
      KYEAR(J) = IWORK(1+(J-1)*4)
      KMONTH(J) = IWORK(2+(J-1)*4)
      KDAY(J) = IWORK(3+(J-1)*4)
      PTIME(J) = TRANSFER(IWORK(4+(J-1)*4), PTIME(J))
    ENDDO
  ENDIF
  !
ENDIF
!
YCOMMENT = TRIM(HREC)
HCOMMENT = YCOMMENT
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
!
!$OMP BARRIER
!
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT1_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT1_FA
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE READ_SURFT2_FA (&
                                 HREC,KL1,KL2,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT2* - routine to read a date
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NLUOUT, CMASK, CPREFIX1D
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF
USE MODI_ABOR1_SFX
USE MODI_ERROR_READ_SURF_FA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
!
!
CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER                                  :: KL1, KL2
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KYEAR    ! year
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KMONTH   ! month
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KDAY     ! day
REAL,    DIMENSION(:,:), INTENT(OUT) :: PTIME    ! year
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment

!*      0.2   Declarations of local variables
!
CHARACTER(LEN=50) :: YCOMMENT
CHARACTER(LEN=6)  :: YMASK
CHARACTER(LEN=18) :: YNAME ! Field Name
LOGICAL           :: GV8, GFOUND
!
INTEGER, DIMENSION(3,SIZE(KYEAR,1),SIZE(KYEAR,2)) :: ITDATE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT2_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
KYEAR=0
KMONTH=0
KDAY=0
PTIME=0.
!
HCOMMENT=""
!
CALL IO_BUFF( HREC,'R',GFOUND)
!
CALL SFX_FA_VERSION(GV8)
IF(GV8)THEN
  YNAME=CPREFIX1D//TRIM(HREC)
ELSE
  YMASK=CMASK
  IF (GFOUND) YMASK='FULL  '
  YNAME=TRIM(YMASK)//TRIM(HREC)
ENDIF
!
WRITE(NLUOUT,*) ' READ_SURFT2_FA : time in 2 dimensions not yet implemented : YNAME=',YNAME
 CALL ABOR1_SFX('MODE_READ_SURF_FA:READ_SURFT2_FA: time in 2 dimensions not yet implemented')
!
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT2_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT2_FA
!
END MODULE MODE_READ_SURF_FA

