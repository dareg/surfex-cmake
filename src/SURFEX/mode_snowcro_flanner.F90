MODULE MODE_SNOWCRO_FLANNER

!!****  SNOWCRO_FLANNER - read "drdt_bst_fit_60.nc" file, which containes parameters from Flanner and Zender, 2006
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
!!     C. Carmagnola
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2013

USE MODI_ABOR1_SFX

USE MODE_READ_CDF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), POINTER :: XDRDT0,XTAU,XKAPPA   ! field read
!
CONTAINS
!
!------------------------------------------------------------------
!
SUBROUTINE READ_FZ06(HFILE)
!
IMPLICIT NONE
!
INCLUDE 'netcdf.inc'
!
!*      1.    declarations of arguments
!
 CHARACTER(LEN=18),  INTENT(IN)  :: HFILE     ! name of file
 CHARACTER(LEN=5),DIMENSION(3),PARAMETER :: HVARNAME=(/'drdt0','tau  ','kappa'/)
!
!*      2.    declarations of local variables
!
INTEGER :: IERROR !error status
!
INTEGER :: ID_FILE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO_FLANNER',0,ZHOOK_HANDLE)
!
!*      3.     Reading of field
!
! Open netcdf file
!
IERROR = NF_OPEN(HFILE,NF_NOWRITE,ID_FILE)
!
CALL HANDLE_ERR_CDF(IERROR,"can't open file "//TRIM(HFILE))
!
CALL READ_VAR_FZ06(ID_FILE,HVARNAME(1),XDRDT0)
CALL READ_VAR_FZ06(ID_FILE,HVARNAME(2),XTAU)
CALL READ_VAR_FZ06(ID_FILE,HVARNAME(3),XKAPPA)
!
! Close netcdf file
IERROR=NF_CLOSE(ID_FILE)
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO_FLANNER',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_FZ06
!------------------------------------------------------------------
SUBROUTINE READ_VAR_FZ06(ID_FILE,HSURF,PVAR)
!
IMPLICIT NONE
!
INCLUDE 'netcdf.inc'
!
INTEGER,INTENT(IN) :: ID_FILE
 CHARACTER(LEN=5),INTENT(IN) :: HSURF
REAL, DIMENSION(:,:,:), POINTER :: PVAR
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IVARDIMSID
!
INTEGER :: INVARDIMS !number of dimensions of netcdf input variable
INTEGER :: ILENDIM,ILENDIM1,ILENDIM2,ILENDIM3
INTEGER :: ID_VAR ! Netcdf IDs for  variable
INTEGER :: IERROR !error status
!
! Look for variable ID
IERROR = NF_INQ_VARID(ID_FILE,TRIM(HSURF),ID_VAR)
CALL HANDLE_ERR_CDF(IERROR,"can't find variable "//TRIM(HSURF))
!
! Number of dimensions
IERROR = NF_INQ_VARNDIMS(ID_FILE,ID_VAR,INVARDIMS)
IF ( IERROR/=NF_NOERR ) CALL HANDLE_ERR_CDF(IERROR,"can't get variable dimensions number")

! Id of dimensions
ALLOCATE(IVARDIMSID(INVARDIMS))
IERROR = NF_INQ_VARDIMID(ID_FILE,ID_VAR,IVARDIMSID)
IF ( IERROR/=NF_NOERR ) CALL HANDLE_ERR_CDF(IERROR,"can't get variable dimensions ids")
!
SELECT CASE (INVARDIMS)
  !
  CASE (3)
    IERROR = NF_INQ_DIMLEN(ID_FILE,IVARDIMSID(1),ILENDIM1)
    IF ( IERROR/=NF_NOERR ) CALL HANDLE_ERR_CDF(IERROR,"can't get variable dimensions lengths")
    IERROR = NF_INQ_DIMLEN(ID_FILE,IVARDIMSID(2),ILENDIM2)
    IF ( IERROR/=NF_NOERR ) CALL HANDLE_ERR_CDF(IERROR,"can't get variable dimensions lengths")
    IERROR = NF_INQ_DIMLEN(ID_FILE,IVARDIMSID(3),ILENDIM3)
    IF ( IERROR/=NF_NOERR ) CALL HANDLE_ERR_CDF(IERROR,"can't get variable dimensions lengths")
    ALLOCATE(PVAR(ILENDIM1,ILENDIM2,ILENDIM3))
  !
  CASE DEFAULT
    CALL ABOR1_SFX('SNOWCRO_FLANNER: incorrect number of dimensions for variable '//TRIM(HSURF))
  !
END SELECT
!
! Read 3D variable
IERROR = NF_GET_VAR_DOUBLE(ID_FILE,ID_VAR,PVAR)
 CALL HANDLE_ERR_CDF(IERROR,"can't read variable "//TRIM(HSURF))
!
END SUBROUTINE READ_VAR_FZ06
!------------------------------------------------------------------
END MODULE MODE_SNOWCRO_FLANNER
