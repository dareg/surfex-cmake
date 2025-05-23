!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_BINLLV (UG, U, USS, &
                              HPROGRAM,HSUBROUTINE,HFILENAME)
!     ##############################################################
!
!!**** *READ_BINLLV* reads a binary latlonvalue file and call treatment 
!!                   subroutine
!!
!!    PURPOSE
!!    -------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson          Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    12/09/95
!!                03/2004  externalization (V. Masson)
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SSO_n, ONLY : SSO_t
!
USE MODD_PGD_GRID,   ONLY : LLATLONMASK
!
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
USE MODI_PT_BY_PT_TREATMENT
USE MODI_GET_LUOUT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SSO_t), INTENT(INOUT) :: USS
!
 CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM      ! Type of program
 CHARACTER(LEN=6),  INTENT(IN) :: HSUBROUTINE   ! Name of the subroutine to call
 CHARACTER(LEN=28), INTENT(IN) :: HFILENAME     ! Name of the field file.
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER      :: IGLB                       ! logical unit
!
INTEGER      :: JLAT, JLON                 ! indexes of OLATLONMASK array
!
INTEGER, PARAMETER :: ILONG=1000
!
REAL*4         :: ZVALUER
REAL, DIMENSION(ILONG) :: ZVALUE          ! values of a data point
REAL*4         :: ZLATR
REAL, DIMENSION(ILONG) :: ZLAT              ! latitude of data point
REAL*4         :: ZLONR, ZLONR2
REAL, DIMENSION(ILONG) :: ZLON              ! longitude of data point
!
INTEGER :: ICPT, ISTAT
INTEGER      :: ILUOUT                     ! output listing
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
!*    1.      Open the file
!             -------------
!
IF (LHOOK) CALL DR_HOOK('READ_BINLLV',0,ZHOOK_HANDLE)
 CALL OPEN_FILE(HPROGRAM,IGLB,HFILENAME,'UNFORMATTED',HACTION='READ')
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
ICPT = 0
!
ZLAT(:) = 0
ZLON(:) = 0
ZVALUE(:) = 0
!
!----------------------------------------------------------------------------
DO
!----------------------------------------------------------------------------
!
!*    3.     Reading of a data point
!            -----------------------
!
  READ(IGLB,IOSTAT=ISTAT) ZLATR,ZLONR,ZVALUER
!
!----------------------------------------------------------------------------
!
!*    4.     Test if point is in the domain
!            ------------------------------
!
  IF (ISTAT==0) THEN
    !
    ZLONR2=ZLONR+NINT((180.-ZLONR)/360.)*360.
    !
    JLAT = 1 + INT( ( ZLATR + 90. ) * 2. )
    JLAT = MIN(JLAT,360)
    JLON = 1 + INT( ( ZLONR2      ) * 2. )
    JLON = MIN(JLON,720)
    !
    IF (.NOT. LLATLONMASK(JLON,JLAT)) CYCLE
    !
    ICPT = ICPT + 1
    !
    IF (ICPT<=ILONG) THEN
      !
      ZLAT  (ICPT) = ZLATR
      ZLON  (ICPT) = ZLONR
      ZVALUE(ICPT) = ZVALUER
      !
    ENDIF
    !
  ENDIF
    !
  IF (ISTAT==-1 .OR. ICPT==ILONG) THEN
    !
    !-------------------------------------------------------------------------------
    !
    !*    5.     Call to the adequate subroutine (point by point treatment)
    !            ----------------------------------------------------------
    !     
    CALL PT_BY_PT_TREATMENT(UG, U, USS, ILUOUT, &
            ZLAT(1:ICPT), ZLON(1:ICPT), ZVALUE(1:ICPT), HSUBROUTINE    )  
    !
    ICPT = 0
    ZLAT  (:) = 0.
    ZLON  (:) = 0.
    ZVALUE(:) = 0.
    !
  ENDIF
  !
  IF (ISTAT==-1) EXIT
  !
  !-------------------------------------------------------------------------------
ENDDO
!
!----------------------------------------------------------------------------
!
!*    8.    Closing of the data file
!           ------------------------
!
 CALL CLOSE_FILE (HPROGRAM,IGLB)
IF (LHOOK) CALL DR_HOOK('READ_BINLLV',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_BINLLV
