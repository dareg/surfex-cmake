!     ##################
      MODULE MODN_IDEAL_FLUX
!     ##################
!
!!****  *MODN_IDEAL_FLUX - declaration of keys for 
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	S. Faroux   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/10/10
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_IDEAL_FLUX, ONLY : XSFTH, XSFTQ, XSFTS, XSFCO2, CUSTARTYPE, XUSTAR, &
                            XZ0, XALB, XEMIS, XTSRAD    
  
IMPLICIT NONE
!
CHARACTER(LEN=7) :: CSFTQ ! Unit for the evaporation flux :
                          !'kg/m2/s'
                          !'W/m2   '
!
!
NAMELIST/NAM_IDEAL_FLUX/XSFTH, CSFTQ, XSFTQ, XSFCO2, CUSTARTYPE, XUSTAR, &
                        XZ0, XALB, XEMIS, XTSRAD
!
END MODULE MODN_IDEAL_FLUX
