!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE AVERAGE_PHY(PFRAC_TILE,             &
                   PTSURF_TILE, PZ0_TILE,            &
                   PZ0H_TILE, PQSURF_TILE,           &   
                   PUREF, PZREF,                     &
                   PTSURF, PZ0, PZ0H, PQSURF         )  
!     ######################################################################
!
!
!!****  *AVERAGE_PHY*  
!!
!!    PURPOSE
!!    -------
!      Average the physical properties from the land and water surfaces depending on the
!      fraction of each surface cover type in the mesh area.
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!      B. Decharme          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/04/2013 
!
!      B. Decharme 07/2015 - Modification to deal with E-zone points in Arome/Aladin
!      F. Svabik   08/2023 - Unapproximated roughness length averaging
!-----------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_SURF_ATM,   ONLY : XZ0_OFFSET
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN) :: PFRAC_TILE ! Fraction in a mesh-area of 
!
REAL, DIMENSION(:,:), INTENT(IN) :: PTSURF_TILE ! surface effective temperature (K)
REAL, DIMENSION(:,:), INTENT(IN) :: PZ0_TILE    ! roughness length for momentum (m)
REAL, DIMENSION(:,:), INTENT(IN) :: PZ0H_TILE   ! roughness length for heat     (m)
REAL, DIMENSION(:,:), INTENT(IN) :: PQSURF_TILE ! specific humidity at surface  (kg/kg)
!
REAL, DIMENSION(:),   INTENT(IN) :: PUREF     ! height of wind forcing                (m)
REAL, DIMENSION(:),   INTENT(IN) :: PZREF     ! height of T,q forcing                 (m)
REAL, DIMENSION(:),   INTENT(OUT):: PTSURF    ! surface effective temperature   (K)
REAL, DIMENSION(:),   INTENT(OUT):: PZ0       ! roughness length for momentum   (m)
REAL, DIMENSION(:),   INTENT(OUT):: PZ0H      ! roughness length for heat       (m)
REAL, DIMENSION(:),   INTENT(OUT):: PQSURF    ! specific humidity at surface    (kg/kg)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PUREF)) :: ZWORK_Z0       ! work array for roughness length for momentum 
REAL, DIMENSION(SIZE(PUREF)) :: ZWORK_Z0H      ! work array for roughness length for heat     
!
INTEGER :: INI, INP  ! dimenssion
INTEGER ::  JI, JP    ! loop counter
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!       0.     Initialization
!              --------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_PHY',0,ZHOOK_HANDLE)
!
INI  = SIZE(PFRAC_TILE,1)
INP  = SIZE(PFRAC_TILE,2)
!
PTSURF     (:)   = 0.
PZ0        (:)   = 0.
PZ0H       (:)   = 0.
PQSURF     (:)   = 0.
!
ZWORK_Z0   (:)   = 0.
ZWORK_Z0H  (:)   = 0.
!
!       1.     Grid-Box average 
!              ----------------
DO JP = 1, INP
  WHERE ( PFRAC_TILE(:,JP) > 0. )
!  
!   surface effective temperature
!
    PTSURF(:) = PTSURF(:) + PFRAC_TILE(:,JP) * PTSURF_TILE(:,JP)
!
!   specific humidity at surface
!
    PQSURF(:) = PQSURF(:) + PFRAC_TILE(:,JP) * PQSURF_TILE(:,JP)
!
!   roughness length for momentum and heat
!
    ZWORK_Z0(:) = ZWORK_Z0(:) + PFRAC_TILE(:,JP) * 1. / (LOG(XZ0_OFFSET + PUREF(:) / PZ0_TILE (:,JP)))**2
    ZWORK_Z0H(:) = ZWORK_Z0H(:) + PFRAC_TILE(:,JP) * 1. / (LOG(XZ0_OFFSET + PZREF(:) / PZ0H_TILE(:,JP)) * LOG(XZ0_OFFSET * (1. + PUREF(:) / PZ0_TILE(:,JP) - PZREF(:) / PZ0H_TILE(:,JP)) + PZREF(:) / PZ0H_TILE(:,JP)))
  END WHERE
END DO
!
WHERE ( ZWORK_Z0(:) /= 0. )
  PZ0(:) = PUREF(:) / (EXP(SQRT(1. / ZWORK_Z0(:))) - XZ0_OFFSET)
  PZ0H(:) = PZREF(:) / (EXP((XZ0_OFFSET * (SQRT(ZWORK_Z0(:) / ZWORK_Z0H(:)) - 1.) + 1.) / SQRT(ZWORK_Z0H(:))) - XZ0_OFFSET)
ELSEWHERE
  PZ0(:) = XUNDEF
  PZ0H(:) = XUNDEF
END WHERE
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_PHY',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE_PHY
