!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
       SUBROUTINE INIT_OUTFN_SURF_ATM_n (DGO, UG, U, CHE, CHU, SV, HSELECT, HPROGRAM, KLUOUT)
!     ###############################
!
!
!!****  *INIT_OUTFN_SURF_ATM_n* -  create output files and defines variables
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
!!      P. LeMoigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05-04
!!      Modified 06-10 by S. Faroux
!!      modified   06-13  B. Decharme  : Add some key
!!                                       Add diag (Qs, Evap, Subl)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CH_EMIS_FIELD_n, ONLY : CH_EMIS_FIELD_t
USE MODD_CH_SURF_n, ONLY : CH_SURF_t
USE MODD_SV_n, ONLY : SV_t
!
USE MODD_DIAG_n, ONLY : DIAG_OPTIONS_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_OL_FILEID,       ONLY : XNETCDF_FILEID_OUT, XNETCDF_FILENAME_OUT
!
!
USE MODN_IO_OFFLINE,      ONLY : XTSTEP_OUTPUT
!
USE MODI_GET_DIM_FULL_n
USE MODI_OL_DEFINE_DIM
USE MODI_GET_DATE_OL
USE MODI_CREATE_FILE
USE MODI_OL_WRITE_COORD
USE MODI_OL_WRITE_PROJ
USE MODD_TYPE_SNOW,ONLY : SURF_SNOW
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
USE NETCDF
!
IMPLICIT NONE
!
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DGO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
TYPE(CH_EMIS_FIELD_t), INTENT(INOUT) :: CHE
TYPE(CH_SURF_t), INTENT(INOUT) :: CHU
TYPE(SV_t), INTENT(INOUT) :: SV
!
 CHARACTER(LEN=*),DIMENSION(:),INTENT(IN) :: HSELECT
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM
INTEGER, INTENT(IN) :: KLUOUT
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
 CHARACTER(LEN=100), DIMENSION(:), POINTER :: YNAME_DIM
 CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT 
 CHARACTER(LEN=40),DIMENSION(1)   :: YDATE
 CHARACTER(LEN=13),DIMENSION(1)   :: YUNIT1, YUNIT2
 CHARACTER(LEN=100)               :: YCOMMENT 
 CHARACTER(LEN=50)                :: YFILE
!
REAL,DIMENSION(:), POINTER    :: ZX, ZY
REAL,DIMENSION(:), POINTER    :: ZLAT,ZLON
TYPE(SURF_SNOW)::TSNOW
!
INTEGER, DIMENSION(:), POINTER   :: IDIMS, IDDIM
INTEGER                          :: IFILE_ID, IVAR_ID, IDIMID
INTEGER                          :: IDIM1, IDIM2, INDIMS
INTEGER                          :: INI, JFILE
INTEGER                          :: JRET
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------

! 1. Compute output lenght dimension
!-----------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_SURF_ATM_N',0,ZHOOK_HANDLE)
!
!
 CALL GET_DIM_FULL_n(U%NDIM_FULL, INI)
!
NULLIFY(ZX)
NULLIFY(ZY)
!
 CALL OL_DEFINE_DIM(UG, U%NSIZE_FULL, HPROGRAM, KLUOUT, INI,TSNOW, &
		   						IDIM1, YUNIT1, YUNIT2, ZX, ZY, IDIMS, IDDIM, YNAME_DIM,   &
                  .FALSE., PLAT=ZLAT,PLON=ZLON)
 CALL GET_DATE_OL(U%TTIME,XTSTEP_OUTPUT,YDATE(1))
!
INDIMS = SIZE(IDDIM)
!
! 4. Create output file for fraction of tiles
!--------------------------------------------
!
YATT_TITLE(1)='units'
!
YFILE='SURF_ATM.OUT.nc'
 CALL CREATE_FILE(YFILE,IDIMS(1:INDIMS-1),YNAME_DIM(1:INDIMS-1),IFILE_ID,IDDIM(1:INDIMS-1))
JRET=NF90_REDEF(IFILE_ID)
!
 CALL OL_WRITE_PROJ(HSELECT,IFILE_ID,UG)
!
DO JFILE = 1,SIZE(XNETCDF_FILENAME_OUT) 
  IF (TRIM(YFILE)==TRIM(XNETCDF_FILENAME_OUT(JFILE))) THEN
    XNETCDF_FILEID_OUT(JFILE) = IFILE_ID
    EXIT
  ENDIF
ENDDO
!
 CALL OL_WRITE_COORD(DGO%CSELECT,YFILE,IFILE_ID,IDDIM(1:INDIMS-1),YATT_TITLE,YNAME_DIM(1:INDIMS-1),&
                     YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY,ZLON,ZLAT)
!
!
! 5. Create output file for diagnostic variables
!-----------------------------------------------
YFILE='SURF_ATM_DIAGNOSTICS.OUT.nc'
 CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
JRET=NF90_REDEF(IFILE_ID)
!
 CALL OL_WRITE_PROJ(HSELECT,IFILE_ID,UG)
!
DO JFILE = 1,SIZE(XNETCDF_FILENAME_OUT) 
  IF (TRIM(YFILE)==TRIM(XNETCDF_FILENAME_OUT(JFILE))) THEN
    XNETCDF_FILEID_OUT(JFILE) = IFILE_ID
    EXIT
  ENDIF
ENDDO
!
IF (CHU%LCH_EMIS .AND. SV%NBEQ>0 .AND. CHU%LCH_SURF_EMIS) THEN
  !
  IF (CHU%CCH_EMIS=='AGGR') JRET = NF90_DEF_DIM(IFILE_ID,"Temporal_emiss",CHE%NTIME_MAX,IDIMID)
  !
ENDIF
!
 CALL OL_WRITE_COORD(DGO%CSELECT,YFILE,IFILE_ID,IDDIM,YATT_TITLE,&
                    YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY,ZLON,ZLAT)
!
IF (ASSOCIATED(ZX)) DEALLOCATE(ZX,ZY)
DEALLOCATE(ZLON,ZLAT)
!
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_SURF_ATM_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_OUTFN_SURF_ATM_n
