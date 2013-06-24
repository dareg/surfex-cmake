MODULE MODE_WRITE_SURF_NC
!
INTERFACE WRITE_SURF0_NC
        MODULE PROCEDURE WRITE_SURFX0_NC
        MODULE PROCEDURE WRITE_SURFN0_NC
        MODULE PROCEDURE WRITE_SURFC0_NC
        MODULE PROCEDURE WRITE_SURFL0_NC
END INTERFACE
INTERFACE WRITE_SURFN_NC
        MODULE PROCEDURE WRITE_SURFX1_NC
        MODULE PROCEDURE WRITE_SURFN1_NC
        MODULE PROCEDURE WRITE_SURFN2_NC
        MODULE PROCEDURE WRITE_SURFL1_NC
        MODULE PROCEDURE WRITE_SURFX2_NC
END INTERFACE
INTERFACE WRITE_SURFT_NC
        MODULE PROCEDURE WRITE_SURFT0_NC
        MODULE PROCEDURE WRITE_SURFT1_NC
        MODULE PROCEDURE WRITE_SURFT2_NC
END INTERFACE
!
CONTAINS
!
!     #############################################################
      SUBROUTINE WRITE_SURFX0_NC(HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *WRITEX0* - routine to read a real scalar
!
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_DEF_VAR_NETCDF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=12),  INTENT(IN)  :: HREC     ! name of the article to be read
REAL,               INTENT(IN)  :: PFIELD   ! the real scalar to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(IN)  :: HCOMMENT
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
INTEGER, DIMENSION(0) :: IDIMS
INTEGER :: IRET
INTEGER :: IVAR_ID,JRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFX0_NC',0,ZHOOK_HANDLE)
!
KRESP=0
!
YATT_TITLE(1) = "comment"
YATT(1) = HCOMMENT
!
IF (NID_NC /= 0) THEN        
  ! 1. Define the variable
  !----------------------------
  CALL DEF_VAR_NETCDF(NID_NC, HREC, HREC, IDIMS, YATT_TITLE, YATT,  IVAR_ID, NF_DOUBLE)
  
  ! 2. Put variable
  !----------------------------
  IRET = NF_PUT_VAR_DOUBLE (NID_NC,IVAR_ID,PFIELD)
ENDIF
!
! 3. Check for errors
!--------------------
IF (NID_NC==0 .OR. IRET.NE.NF_NOERR) KRESP=1
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFX0_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_SURFX0_NC
!
!     #############################################################
      SUBROUTINE WRITE_SURFN0_NC(HREC,KFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *WRITEN0* - routine to read an integer
!
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_DEF_VAR_NETCDF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=12),  INTENT(IN) :: HREC     ! name of the article to be read
INTEGER,            INTENT(IN)  :: KFIELD   ! the integer scalar to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(IN) :: HCOMMENT
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
INTEGER, DIMENSION(0) :: IDIMS
INTEGER              :: IVAR_ID, JRET
INTEGER :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFN0_NC',0,ZHOOK_HANDLE)
!
KRESP=0
!
YATT_TITLE(1) = "comment"
YATT(1) = HCOMMENT
!
! 0. find filename
! -----------------
!
IF (NID_NC /= 0) THEN    
  ! 1. Find id of the variable
  !----------------------------
  CALL DEF_VAR_NETCDF(NID_NC, HREC, HREC, IDIMS, YATT_TITLE, YATT, IVAR_ID, NF_INT)
  ! 
  ! 2. Get variable
  !----------------------------
  IRET=NF_PUT_VAR_INT(NID_NC,IVAR_ID,KFIELD)
ENDIF
!
! 3. Check for errors
!--------------------
IF (NID_NC==0 .OR. IRET.NE.NF_NOERR) KRESP=1
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFN0_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_SURFN0_NC
!
!     #############################################################
      SUBROUTINE WRITE_SURFC0_NC(HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *WRITEC0* - routine to read a STRING
!
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_DEF_VAR_NETCDF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=12),  INTENT(IN)  :: HREC     ! name of the article to be read
 CHARACTER(LEN=40),  INTENT(IN)  :: HFIELD   ! the integer scalar to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(IN)  :: HCOMMENT
!
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
INTEGER, DIMENSION(1) :: IDIMS
 CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: YFIELD
INTEGER :: IVAR_ID, JRET
INTEGER :: IRET, J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFC0_NC',0,ZHOOK_HANDLE)
!
KRESP=0
!
YATT_TITLE(1) = "comment"
YATT(1) = HCOMMENT
!
!
IF (NID_NC /= 0) THEN 
  ! 1. Find id of the variable
  !----------------------------
  IRET = NF_INQ_DIMID(NID_NC,'char_len',IDIMS(1))
  !
  CALL DEF_VAR_NETCDF(NID_NC, HREC, HREC, IDIMS, YATT_TITLE, YATT, IVAR_ID, NF_CHAR,LEN_TRIM(HFIELD))
  !
  ! 2. Get variable
  !----------------------------
  ALLOCATE(YFIELD(LEN(HFIELD)))
  DO J=1,LEN(HFIELD)
    YFIELD(J) = HFIELD(J:J)
  ENDDO
  IRET=NF_PUT_VAR_TEXT(NID_NC,IVAR_ID,YFIELD)
ENDIF
!
! 3. Check for errors
!--------------------
IF (NID_NC==0 .OR. IRET.NE.NF_NOERR) KRESP=1
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFC0_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_SURFC0_NC
!
!     #############################################################
      SUBROUTINE WRITE_SURFL0_NC(HREC,OFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *WRITEL0* - routine to read a logical
!
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_DEF_VAR_NETCDF
!
USE MODI_HANDLE_ERR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=12),   INTENT(IN) :: HREC     ! name of the article to be read
LOGICAL,             INTENT(IN) :: OFIELD   ! array containing the data field
INTEGER,             INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),  INTENT(IN) :: HCOMMENT
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
INTEGER, DIMENSION(0) :: IDIMS
CHARACTER(LEN=1)    :: YFIELD   ! work array read in the file
INTEGER :: IFIELD
INTEGER              :: IVAR_ID, JRET
INTEGER :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFL0_NC',0,ZHOOK_HANDLE)
!
KRESP=0
!
YATT_TITLE(1) = "comment"
YATT(1) = HCOMMENT
!
IF (NID_NC /= 0) THEN        
  ! 1. Find id of the variable
  !----------------------------
  IF (OFIELD) THEN
    YFIELD = 'T'
    IFIELD = 1
  ELSE
    YFIELD = 'F'
    IFIELD = 0
  ENDIF
  !
  ! 2. Put variable
  !----------------------------
  CALL DEF_VAR_NETCDF(NID_NC, HREC, HREC, IDIMS, YATT_TITLE, YATT, IVAR_ID, NF_INT)
  !
  IRET=NF_PUT_VAR_INT(NID_NC,IVAR_ID,IFIELD)
  CALL HANDLE_ERR(IRET,HREC)
ENDIF
!
! 3. Check for errors
!--------------------
IF (NID_NC==0 .OR. IRET.NE.NF_NOERR) KRESP=1
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFL0_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_SURFL0_NC
!
!
!     #############################################################
      SUBROUTINE WRITE_SURFX1_NC(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *WRITEX1* - routine to fill a real 1D array for the externalised surface 
! 
USE MODD_SURFEX_OMP, ONLY : LWORK0
!
USE MODD_IO_SURF_NC, ONLY : NID_NC, LMASK, NMASK, NMASK_IGN, CFILEOUT_NC
!
USE MODI_DEF_VAR_NETCDF
!
USE MODI_IO_BUFF_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=12),   INTENT(IN) :: HREC     ! name of the article to be read
REAL, DIMENSION(:),  INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,             INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),  INTENT(IN) :: HCOMMENT
 CHARACTER(LEN=1),    INTENT(IN) :: HDIR     ! type of field :
                                            ! 'H' : field with
                                            !       horizontal spatial dim.
                                            ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
INTEGER, DIMENSION(2) :: IDIMIDS
INTEGER, DIMENSION(2) :: IDIMLEN
 CHARACTER(LEN=100)    :: YNAME
INTEGER               :: IVAR_ID, JDIM, INDIMS
INTEGER               :: JRET
INTEGER,DIMENSION(5)  :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFX1_NC',0,ZHOOK_HANDLE)
!
KRESP = 0
!
IRET(:) = 0
!
YATT_TITLE(1) = "comment"
YATT(1) = HCOMMENT
!
 CALL IO_BUFF_n(HREC,'W',LWORK0)
!
IF (LWORK0 .AND. LHOOK) CALL DR_HOOK("WRITE_SURF_NC:WRITE_SURFX1_NC",1,ZHOOK_HANDLE)
IF (LWORK0) RETURN
!
IF (NID_NC /= 0) THEN 
  !    
  ! 0. find filename
  ! -----------------
  !
  IRET(1) = NF_INQ_NDIMS(NID_NC,INDIMS)
  IRET(2) = NF_INQ_DIMID(NID_NC,'Number_of_points',IDIMIDS(1))
  IF (IRET(2)/=0) THEN
    IRET(2) = NF_INQ_DIMID(NID_NC,'lon',IDIMIDS(1))
    IF (IRET(2)/=0) THEN
      IRET(2) = NF_INQ_DIMID(NID_NC,'xx',IDIMIDS(1))
      IRET(3) = NF_INQ_DIMID(NID_NC,'yy',IDIMIDS(2))
    ELSE
      IRET(3) = NF_INQ_DIMID(NID_NC,'lat',IDIMIDS(2))
    ENDIF
  ENDIF
  DO JDIM=1,SIZE(IDIMIDS)
    JRET=NF_INQ_DIMLEN(NID_NC,IDIMIDS(JDIM),IDIMLEN(JDIM))
  ENDDO
  !
  IRET(4)=NF_INQ_DIMNAME(NID_NC,IDIMIDS(1),YNAME)
  !
  IF (YNAME .EQ. 'Number_of_points') THEN
    CALL DEF_VAR_NETCDF(NID_NC, HREC, HREC, IDIMIDS(1:1), YATT_TITLE, YATT, IVAR_ID, NF_DOUBLE)
    CALL WRITE_DATAX1_NC(IDIMLEN(1),INDIMS)
  ELSE
   CALL DEF_VAR_NETCDF(NID_NC, HREC, HREC, IDIMIDS(1:2), YATT_TITLE, YATT, IVAR_ID, NF_DOUBLE)
   CALL WRITE_DATAX1_NC(IDIMLEN(1)*IDIMLEN(2),INDIMS)
  ENDIF
  !
ENDIF
!
DO JRET=1,5
 IF (NID_NC==0 .OR. IRET(JRET).NE.NF_NOERR) KRESP = 1
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFX1_NC',1,ZHOOK_HANDLE)
!
CONTAINS
!
SUBROUTINE WRITE_DATAX1_NC(KDIM,KNDIMS)
!
USE MODI_GATHER_AND_WRITE_MPI
USE MODI_UNPACK_SAME_RANK
USE MODI_HANDLE_ERR
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KDIM
INTEGER, INTENT(IN) :: KNDIMS
!
REAL, DIMENSION(KDIM) :: ZTAB1D
REAL, DIMENSION(KDIM) :: ZWORK_IGN
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFX1_NC:WRITE_DATAX1_NC',0,ZHOOK_HANDLE)
!
IF(.NOT.ALLOCATED(NMASK_IGN))THEN
  IF (LMASK) THEN
    CALL GATHER_AND_WRITE_MPI(PFIELD,ZTAB1D,NMASK)
  ELSE 
    CALL GATHER_AND_WRITE_MPI(PFIELD,ZTAB1D)
  ENDIF
ELSE
  !ign grid 
  IF (LMASK) THEN
    CALL GATHER_AND_WRITE_MPI(PFIELD,ZWORK_IGN(1:SIZE(NMASK_IGN)),NMASK)
  ELSE 
    CALL GATHER_AND_WRITE_MPI(PFIELD,ZWORK_IGN(1:SIZE(NMASK_IGN)))
  ENDIF
  CALL UNPACK_SAME_RANK(NMASK_IGN,ZWORK_IGN(1:SIZE(NMASK_IGN)),ZTAB1D)
ENDIF
!
IRET(5)=NF_PUT_VAR_DOUBLE(NID_NC,IVAR_ID,ZTAB1D)
!
CALL HANDLE_ERR(IRET(5),HREC)
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFX1_NC:WRITE_DATAX1_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_DATAX1_NC
!
END SUBROUTINE WRITE_SURFX1_NC
!
!     #############################################################
      SUBROUTINE WRITE_SURFX2_NC(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *WRITEX2* - routine to fill a real 2D array for the externalised surface 
!
USE MODD_IO_SURF_NC, ONLY : NID_NC, LMASK, NMASK, NMASK_IGN
!
USE MODI_DEF_VAR_NETCDF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=12),    INTENT(IN) :: HREC     ! name of the article to be read
REAL, DIMENSION(:,:), INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,              INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),   INTENT(IN) :: HCOMMENT
 CHARACTER(LEN=1),     INTENT(IN) :: HDIR     ! type of field :
                                             ! 'H' : field with
                                             !       horizontal spatial dim.
                                             ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
INTEGER, DIMENSION(3) :: IDIMIDS
INTEGER, DIMENSION(3) :: IDIMLEN
 CHARACTER(LEN=100)    :: YNAME
INTEGER               :: IVAR_ID, JDIM, INDIMS
INTEGER               :: JRET
INTEGER, DIMENSION(5) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFX2_NC',0,ZHOOK_HANDLE)
!
KRESP = 0
!
IRET(:) = 0
!
YATT_TITLE(1) = "comment"
YATT(1) = HCOMMENT
!
IF (NID_NC /= 0) THEN 
  !    
  ! 0. find filename
  ! -----------------
  !
  IRET(1) = NF_INQ_NDIMS(NID_NC,INDIMS)
  IRET(2) = NF_INQ_DIMID(NID_NC,'Number_of_points',IDIMIDS(1))
  IF (IRET(2)==0) THEN
    IRET(3) = NF_INQ_DIMID(NID_NC,'Number_of_Tile',IDIMIDS(2))
  ELSE
    IRET(2) = NF_INQ_DIMID(NID_NC,'lon',IDIMIDS(1))
    IF (IRET(2)/=0) THEN
      IRET(2) = NF_INQ_DIMID(NID_NC,'xx',IDIMIDS(1))
      IRET(3) = NF_INQ_DIMID(NID_NC,'yy',IDIMIDS(2))
    ELSE
      IRET(3) = NF_INQ_DIMID(NID_NC,'lat',IDIMIDS(2))
    ENDIF
    IRET(4) = NF_INQ_DIMID(NID_NC,'Number_of_Tile',IDIMIDS(3))
  ENDIF
  DO JDIM=1,SIZE(IDIMLEN)
    JRET=NF_INQ_DIMLEN(NID_NC,IDIMIDS(JDIM),IDIMLEN(JDIM))
  ENDDO
  !
  IRET(5)=NF_INQ_DIMNAME(NID_NC,IDIMIDS(1),YNAME)
  !
  IF (YNAME .EQ. 'Number_of_points') THEN
    CALL DEF_VAR_NETCDF(NID_NC, HREC, HREC, IDIMIDS(1:2), YATT_TITLE, YATT, IVAR_ID, NF_DOUBLE)
    CALL WRITE_DATAX2_NC(IDIMLEN(1),IDIMLEN(2),INDIMS)
  ELSE
    CALL DEF_VAR_NETCDF(NID_NC, HREC, HREC, IDIMIDS(1:3), YATT_TITLE, YATT, IVAR_ID, NF_DOUBLE)
    CALL WRITE_DATAX2_NC(IDIMLEN(1)*IDIMLEN(2),IDIMLEN(3),INDIMS)
  ENDIF
  !
ENDIF
!
DO JRET=1,5
 IF (NID_NC==0 .OR. IRET(JRET).NE.NF_NOERR) KRESP = 1
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFX2_NC',1,ZHOOK_HANDLE)
!
CONTAINS
!
SUBROUTINE WRITE_DATAX2_NC(KDIM1,KDIM2,KNDIMS)
!
USE MODI_GATHER_AND_WRITE_MPI
USE MODI_UNPACK_SAME_RANK
USE MODI_HANDLE_ERR
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KDIM1
INTEGER, INTENT(IN) :: KDIM2
INTEGER, INTENT(IN) :: KNDIMS
!
REAL, DIMENSION(KDIM1,KDIM2) :: ZTAB2D    ! work array read in the file
REAL, DIMENSION(KDIM1,KDIM2) :: ZWORK_IGN ! work array read in the file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFX2_NC:WRITE_DATAX2_NC',0,ZHOOK_HANDLE)
!
IF(.NOT.ALLOCATED(NMASK_IGN))THEN
  IF (LMASK) THEN
    CALL GATHER_AND_WRITE_MPI(PFIELD,ZTAB2D,NMASK)
  ELSE 
    CALL GATHER_AND_WRITE_MPI(PFIELD,ZTAB2D)
  ENDIF
ELSE
  !ign grid 
  IF (LMASK) THEN
    CALL GATHER_AND_WRITE_MPI(PFIELD,ZWORK_IGN(1:SIZE(NMASK_IGN),:),NMASK)
  ELSE 
    CALL GATHER_AND_WRITE_MPI(PFIELD,ZWORK_IGN(1:SIZE(NMASK_IGN),:))
  ENDIF
  CALL UNPACK_SAME_RANK(NMASK_IGN,ZWORK_IGN(1:SIZE(NMASK_IGN),:),ZTAB2D)
ENDIF
!
IRET(5)=NF_PUT_VAR_DOUBLE(NID_NC,IVAR_ID,ZTAB2D)
!    
CALL HANDLE_ERR(IRET(5),HREC)  
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFX2_NC:WRITE_DATAX2_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_DATAX2_NC
!
END SUBROUTINE WRITE_SURFX2_NC

!     #############################################################
      SUBROUTINE WRITE_SURFN1_NC(HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *WRITEN0* - routine to read an integer
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=12),      INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:),  INTENT(IN)  :: KFIELD   ! the integer scalar to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),     INTENT(IN)  :: HCOMMENT
 CHARACTER(LEN=1),       INTENT(IN) :: HDIR     ! type of field :
                                               ! 'H' : field with
                                               !       horizontal spatial dim.
                                               ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(KFIELD)) :: ZFIELD 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFN1_NC',0,ZHOOK_HANDLE)
!
ZFIELD=FLOAT(KFIELD)
 CALL WRITE_SURFX1_NC(HREC,ZFIELD,KRESP,HCOMMENT,HDIR)
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFN1_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_SURFN1_NC
!

!     #############################################################
      SUBROUTINE WRITE_SURFN2_NC(HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *WRITEN0* - routine to read an integer
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=12),      INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:,:),  INTENT(IN)  :: KFIELD   ! the integer scalar to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),     INTENT(IN)  :: HCOMMENT
 CHARACTER(LEN=1),       INTENT(IN) :: HDIR     ! type of field :
                                               ! 'H' : field with
                                               !       horizontal spatial dim.
                                               ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(KFIELD,1),SIZE(KFIELD,2)) :: ZFIELD 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFN2_NC',0,ZHOOK_HANDLE)
!
ZFIELD=FLOAT(KFIELD)
 CALL WRITE_SURFX2_NC(HREC,ZFIELD,KRESP,HCOMMENT,HDIR)
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFN2_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_SURFN2_NC
!
!     #############################################################
      SUBROUTINE WRITE_SURFL1_NC(HREC,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *WRITEL1* - routine to read a logical array
!    
USE MODD_IO_SURF_NC, ONLY : NID_NC, LMASK, NMASK, NMASK_IGN
!
USE MODI_DEF_VAR_NETCDF
!
USE MODI_HANDLE_ERR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),      INTENT(IN) :: HREC     ! name of the article to be read
LOGICAL, DIMENSION(:), INTENT(IN) :: OFIELD   ! array containing the data field
INTEGER,               INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),    INTENT(IN) :: HCOMMENT
 CHARACTER(LEN=1),      INTENT(IN) :: HDIR     ! type of field :
                                              ! 'H' : field with
                                              !       horizontal spatial dim.
                                              ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
INTEGER, DIMENSION(1) :: IDIMIDS
INTEGER, DIMENSION(1) :: IDIMLEN
INTEGER               :: IVAR_ID, JDIM, INDIMS
INTEGER               :: JRET
INTEGER, DIMENSION(3) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFL1_NC',0,ZHOOK_HANDLE)
!
IRET(:) = 0
!
YATT_TITLE(1) = "comment"
YATT(1) = HCOMMENT
!
KRESP=0
!
IF (NID_NC /= 0) THEN 
  !    
  ! 0. find filename
  ! -----------------
  !
  IRET(1) = NF_INQ_DIMID(NID_NC,'Number_of_covers',IDIMIDS(1))
  IRET(2) = NF_INQ_DIMLEN(NID_NC,IDIMIDS(1),IDIMLEN(1))
  !
  CALL DEF_VAR_NETCDF(NID_NC, HREC, HREC, IDIMIDS(1:1), YATT_TITLE, YATT, IVAR_ID, NF_INT)
  CALL WRITE_DATAL1_NC(IDIMLEN(1))
  !
ENDIF
!
DO JRET=1,3
  IF (NID_NC==0 .OR. IRET(JRET).NE.NF_NOERR) KRESP=1
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFL1_NC',1,ZHOOK_HANDLE)
!
CONTAINS
!
SUBROUTINE WRITE_DATAL1_NC(KDIM)
!
INTEGER, INTENT(IN) :: KDIM
!
 CHARACTER(LEN=1), DIMENSION(KDIM) :: YTAB1D  ! work array read in the file
INTEGER, DIMENSION(KDIM) :: ITAB1D
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFL1_NC:WRITE_DATAL1_NC',0,ZHOOK_HANDLE)
!
YTAB1D(:) = ""
ITAB1D(:) = 0
!
DO JRET=1,SIZE(OFIELD)
  IF (OFIELD(JRET)) THEN
    YTAB1D(JRET) ='T'
    ITAB1D(JRET) = 1
  ELSE
    YTAB1D(JRET) ='F'
    ITAB1D(JRET) = 0
  ENDIF
ENDDO  
!
! 2. Put variable
!-----------------
IRET(3)=NF_PUT_VAR_INT(NID_NC,IVAR_ID,ITAB1D)
!
 CALL HANDLE_ERR(IRET(3),HREC)
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFL1_NC:WRITE_DATAL1_NC',1,ZHOOK_HANDLE)
END SUBROUTINE WRITE_DATAL1_NC
!
END SUBROUTINE WRITE_SURFL1_NC
!
!
!     #############################################################
      SUBROUTINE WRITE_SURFT0_NC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *WRITET0* - routine to read a NETCDF  date_time scalar
!
USE MODD_IO_SURF_NC, ONLY : NID_NC, LMASK, NMASK, NMASK_IGN
!
USE MODI_DEF_VAR_NETCDF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=12),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,            INTENT(IN)  :: KYEAR    ! year
INTEGER,            INTENT(IN)  :: KMONTH   ! month
INTEGER,            INTENT(IN)  :: KDAY     ! day
REAL,               INTENT(IN)  :: PTIME    ! time
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(IN)  :: HCOMMENT
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
INTEGER, DIMENSION(0) :: IDIMIDS
 CHARACTER(LEN=12) :: YRECFM    ! Name of the article to be written
INTEGER :: IVAR_ID, JRET, JWRK
INTEGER :: JLEN
INTEGER,DIMENSION(4) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFT0_NC',0,ZHOOK_HANDLE)
!
KRESP=0
!
YATT_TITLE(1) = "comment"
YATT(1) = HCOMMENT
!
DO JWRK=1,4
  !
  IF (JWRK == 1) THEN 
    YRECFM = TRIM(HREC)//'-YEAR'
  ELSEIF (JWRK == 2) THEN
    YRECFM = TRIM(HREC)//'-MONTH'
  ELSEIF (JWRK == 3) THEN
    YRECFM = TRIM(HREC)//'-DAY'
  ELSEIF (JWRK == 4) THEN
    YRECFM=TRIM(HREC)//'-TIME'
  ENDIF
  ! 0. find filename
  ! -----------------
  !
  IF (NID_NC /= 0) THEN
    !
    !
    IF (JWRK==1) THEN
      CALL DEF_VAR_NETCDF(NID_NC, YRECFM, YRECFM, IDIMIDS, YATT_TITLE, YATT, IVAR_ID,NF_INT)
      IRET(JWRK)=NF_PUT_VAR_INT(NID_NC,IVAR_ID,KYEAR)
    ELSEIF (JWRK==2) THEN
      CALL DEF_VAR_NETCDF(NID_NC, YRECFM, YRECFM, IDIMIDS, YATT_TITLE, YATT, IVAR_ID,NF_INT)
      IRET(JWRK)=NF_PUT_VAR_INT(NID_NC,IVAR_ID,KMONTH)
    ELSEIF (JWRK==3) THEN
      CALL DEF_VAR_NETCDF(NID_NC, YRECFM, YRECFM, IDIMIDS, YATT_TITLE, YATT, IVAR_ID,NF_INT)
      IRET(JWRK)=NF_PUT_VAR_INT(NID_NC,IVAR_ID,KDAY)
    ELSEIF (JWRK==4) THEN
      CALL DEF_VAR_NETCDF(NID_NC, YRECFM, YRECFM, IDIMIDS, YATT_TITLE, YATT, IVAR_ID,NF_DOUBLE)
      IRET(JWRK)=NF_PUT_VAR_DOUBLE(NID_NC,IVAR_ID,PTIME)
    ENDIF
  ENDIF
ENDDO
!
! 3. Check for errors
!--------------------
DO JRET=1,4
  IF (NID_NC==0.OR.IRET(JRET).NE.NF_NOERR) KRESP=1
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFT0_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_SURFT0_NC
!
!     #############################################################
      SUBROUTINE    WRITE_SURFT1_NC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *WRITET0* - routine to read a NETCDF  date_time scalar
!
USE MODD_IO_SURF_NC, ONLY : NID_NC, LMASK, NMASK, NMASK_IGN
!
USE MODI_DEF_VAR_NETCDF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=12),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:), INTENT(IN)  :: KYEAR    ! year
INTEGER, DIMENSION(:), INTENT(IN)  :: KMONTH   ! month
INTEGER, DIMENSION(:), INTENT(IN)  :: KDAY     ! day
REAL, DIMENSION(:), INTENT(IN)  :: PTIME    ! time
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(IN)  :: HCOMMENT
!
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=100)    :: YNAME
 CHARACTER(LEN=12) :: YRECFM    ! Name of the article to be written
INTEGER :: JRET, JWRK
INTEGER,DIMENSION(4) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFT1_NC',0,ZHOOK_HANDLE)
!
KRESP=0
!
DO JWRK=1,4
  !
  IF (JWRK == 1) THEN 
    YRECFM = TRIM(HREC)//'-YEAR'
  ELSEIF (JWRK == 2) THEN
    YRECFM = TRIM(HREC)//'-MONTH'
  ELSEIF (JWRK == 3) THEN
    YRECFM = TRIM(HREC)//'-DAY'
  ELSEIF (JWRK == 4) THEN
    YRECFM=TRIM(HREC)//'-TIME'
  ENDIF
  !
  IF (JWRK==1) THEN
    CALL WRITE_SURFN1_NC(YRECFM,KYEAR,IRET(JWRK),HCOMMENT,'-')
  ELSEIF (JWRK==2) THEN
    CALL WRITE_SURFN1_NC(YRECFM,KMONTH,IRET(JWRK),HCOMMENT,'-')
  ELSEIF (JWRK==3) THEN
    CALL WRITE_SURFN1_NC(YRECFM,KDAY,IRET(JWRK),HCOMMENT,'-')
  ELSEIF (JWRK==4) THEN
    CALL WRITE_SURFX1_NC(YRECFM,PTIME,IRET(JWRK),HCOMMENT,'-')
  ENDIF
  !
ENDDO
!
! 3. Check for errors
!--------------------
DO JRET=1,4
  IF (NID_NC==0.OR.IRET(JRET).NE.NF_NOERR) KRESP=1
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFT1_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_SURFT1_NC
!
!     #############################################################
      SUBROUTINE    WRITE_SURFT2_NC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *WRITET0* - routine to read a NETCDF  date_time scalar
!
USE MODD_IO_SURF_NC, ONLY : NID_NC, LMASK, NMASK, NMASK_IGN
!
USE MODI_DEF_VAR_NETCDF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=12),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:,:), INTENT(IN)  :: KYEAR    ! year
INTEGER, DIMENSION(:,:), INTENT(IN)  :: KMONTH   ! month
INTEGER, DIMENSION(:,:), INTENT(IN)  :: KDAY     ! day
REAL, DIMENSION(:,:),    INTENT(IN)  :: PTIME    ! time
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(IN)  :: HCOMMENT
!
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=100)    :: YNAME
 CHARACTER(LEN=12) :: YRECFM    ! Name of the article to be written
INTEGER :: JRET, JWRK
INTEGER,DIMENSION(4) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFT2_NC',0,ZHOOK_HANDLE)
!
KRESP=0
!
DO JWRK=1,4
  !
  IF (JWRK == 1) THEN 
    YRECFM = TRIM(HREC)//'-YEAR'
  ELSEIF (JWRK == 2) THEN
    YRECFM = TRIM(HREC)//'-MONTH'
  ELSEIF (JWRK == 3) THEN
    YRECFM = TRIM(HREC)//'-DAY'
  ELSEIF (JWRK == 4) THEN
    YRECFM=TRIM(HREC)//'-TIME'
  ENDIF
  !
  IF (JWRK==1) THEN
    CALL WRITE_SURFN2_NC(YRECFM,KYEAR,IRET(JWRK),HCOMMENT,'-')
  ELSEIF (JWRK==2) THEN
    CALL WRITE_SURFN2_NC(YRECFM,KMONTH,IRET(JWRK),HCOMMENT,'-')
  ELSEIF (JWRK==3) THEN
    CALL WRITE_SURFN2_NC(YRECFM,KDAY,IRET(JWRK),HCOMMENT,'-')
  ELSEIF (JWRK==4) THEN
    CALL WRITE_SURFX2_NC(YRECFM,PTIME,IRET(JWRK),HCOMMENT,'-')
  ENDIF
  !
ENDDO
!
! 3. Check for errors
!--------------------
DO JRET=1,4
  IF (NID_NC==0.OR.IRET(JRET).NE.NF_NOERR) KRESP=1
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_NC:WRITE_SURFT2_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_SURFT2_NC
!
END MODULE MODE_WRITE_SURF_NC
