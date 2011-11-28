
!     ##################
      MODULE MODN_PREP_SEAFLUX
!     ##################
!
!!****  *MODN_PREP_SEAFLUX* - declaration of namelist NAM_PREP_SEAFLUX
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_PREP_SEAFLUX
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
!!	S.Malardel    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2003                    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PREP_SEAFLUX, ONLY : CFILE_SEAFLX, CTYPE, XSST_UNIF

!
IMPLICIT NONE
!
INTEGER           :: NYEAR            ! YEAR for surface
INTEGER           :: NMONTH           ! MONTH for surface
INTEGER           :: NDAY             ! DAY for surface
REAL              :: XTIME            ! TIME for surface
LOGICAL           :: LSEA_SBL         ! flag to use air layers inside the SBL
LOGICAL           :: LOCEAN_MERCATOR  ! oceanic variables initialized from 
                                      !   MERCATOR if true
LOGICAL           :: LOCEAN_CURRENT   ! initial ocean state with current 
                                      !   (if false ucur=0, vcur=0)
!
NAMELIST/NAM_PREP_SEAFLUX/CFILE_SEAFLX, CTYPE, XSST_UNIF,       &
                            NYEAR, NMONTH, NDAY, XTIME, LSEA_SBL, &
                            LOCEAN_MERCATOR, LOCEAN_CURRENT  
!
END MODULE MODN_PREP_SEAFLUX
