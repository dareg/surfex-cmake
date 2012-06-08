!     ##################
      MODULE MODN_PREP_TEB
!     ##################
!
!!****  *MODN_PREP_TEB* - declaration of namelist NAM_PREP_TEB
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_PREP_TEB
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!       
!!    AUTHOR
!!    ------
!!	V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004                    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PREP_TEB, ONLY : CFILE_TEB, CTYPE, CFILEPGD_TEB, CTYPEPGD,          &
                            CFILE_WS, CTYPE_WS, CFILE_TS, CTYPE_TS,          &
                            XWS_ROOF, XWS_ROAD,                              &
                            XTS_ROOF, XTS_ROAD, XTS_WALL, XTI_BLD, XTI_ROAD  

!
IMPLICIT NONE
!
INTEGER           :: NYEAR        ! YEAR for surface
INTEGER           :: NMONTH       ! MONTH for surface
INTEGER           :: NDAY         ! DAY for surface
REAL              :: XTIME        ! TIME for surface
LOGICAL           :: LTEB_CANOPY  ! flag to use air layers inside the canopy
!
NAMELIST/NAM_PREP_TEB/CFILE_TEB, CTYPE, CFILEPGD_TEB, CTYPEPGD,  &
                      CFILE_WS, CTYPE_WS, XWS_ROOF, XWS_ROAD,    &
                      CFILE_TS, CTYPE_TS, XTS_ROOF, XTS_ROAD,    &
                      XTS_WALL, XTI_BLD, XTI_ROAD,               &
                      NYEAR, NMONTH, NDAY, XTIME, LTEB_CANOPY  
!
END MODULE MODN_PREP_TEB
