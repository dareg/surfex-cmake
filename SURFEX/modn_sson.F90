!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!##################
MODULE MODN_SSO_n
!##################
!
!!****  *MODN_SSO_n* - declaration of namelist NAM_SSO_n
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_SSO_n
!     which concern the roughness parameterization for orography.
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
!!      V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2010                    
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
CHARACTER(LEN=4)  :: CROUGH
REAL              :: XFRACZ0
REAL              :: XCOEFBE
REAL              :: XCOROT
REAL              :: XVOROT
REAL              :: XSOROT
LOGICAL           :: LDSV, LDSH, LDSL
!
NAMELIST/NAM_SSOn/CROUGH, XFRACZ0, XCOEFBE, XCOROT, XSOROT, XVOROT, &
                  LDSV, LDSH, LDSL
!
CONTAINS
!
SUBROUTINE INIT_NAM_SSOn (USS)
!
  USE MODD_SSO_n, ONLY : SSO_t
!
  IMPLICIT NONE
!
  TYPE(SSO_t), INTENT(INOUT) :: USS
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODN_SSO_N:INIT_NAM_SSON',0,ZHOOK_HANDLE)
  CROUGH = USS%CROUGH
  XFRACZ0= USS%XFRACZ0  
  XCOEFBE= USS%XCOEFBE  
  XCOROT= USS%XCOROT
  XSOROT= USS%XSOROT
  XVOROT= USS%XVOROT
  LDSV= USS%LDSV
  LDSH= USS%LDSH
  LDSL= USS%LDSL
  IF (LHOOK) CALL DR_HOOK('MODN_SSO_N:INIT_NAM_SSON',1,ZHOOK_HANDLE)
END SUBROUTINE INIT_NAM_SSOn

SUBROUTINE UPDATE_NAM_SSOn (USS)
!
  USE MODD_SSO_n, ONLY : SSO_t
!
  IMPLICIT NONE
!
  TYPE(SSO_t), INTENT(INOUT) :: USS
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('MODN_SSO_N:UPDATE_NAM_SSON',0,ZHOOK_HANDLE)
  USS%CROUGH = CROUGH
  USS%XFRACZ0= XFRACZ0  
  USS%XCOEFBE= XCOEFBE  
  USS%XCOROT = XCOROT
  USS%XSOROT = XSOROT
  USS%XVOROT = XVOROT
  USS%LDSV   = LDSV
  USS%LDSH   = LDSH
  USS%LDSL   = LDSL  
  IF (LHOOK) CALL DR_HOOK('MODN_SSO_N:UPDATE_NAM_SSON',1,ZHOOK_HANDLE)
END SUBROUTINE UPDATE_NAM_SSOn

END MODULE MODN_SSO_n
