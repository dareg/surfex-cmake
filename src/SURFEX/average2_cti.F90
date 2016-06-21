!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      SUBROUTINE AVERAGE2_CTI
!     #######################
!
!!**** *AVERAGE2_CTI* computes the topo index stats
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    B. Decharme         Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    06/2009
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_PGDWORK,       ONLY : NSIZE, XSUMVAL,  &
                               XMEAN_WORK, XSTD_WORK, XSKEW_WORK, &
                               XMIN_WORK, XMAX_WORK 
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
REAL, DIMENSION(SIZE(NSIZE,1)) :: ZSIZE
!
integer :: I
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGE2_CTI',0,ZHOOK_HANDLE)
ZSIZE(:)=REAL(NSIZE(:,1))
!
WHERE (NSIZE(:,1)>=36)
!
!----------------------------------------------------------------------------
!
!*    1.     Mean CTI
!            --------------
!
  XMEAN_WORK(:) = XSUMVAL(:,1)/ZSIZE(:)
!
!-------------------------------------------------------------------------------
!
!*    2.     Standard deviation
!            ------------------
!
  WHERE (XMAX_WORK(:)-XMIN_WORK(:)>=1.0) 
    XSTD_WORK(:) = SQRT( MAX(0.,XSUMVAL(:,2)/NSIZE(:,1) - XMEAN_WORK(:)*XMEAN_WORK(:)) )
  ELSEWHERE
    XSTD_WORK(:) = 0.0
  END WHERE
!
!-------------------------------------------------------------------------------
!
!*    3.     Skewness
!            --------
!
  WHERE(XSTD_WORK(:)>0.0)
!          
        XSKEW_WORK(:) = XSUMVAL(:,3)-ZSIZE(:)*XMEAN_WORK(:)*XMEAN_WORK(:)*XMEAN_WORK(:) &
                       -3.0*ZSIZE(:)*XMEAN_WORK(:)*XSTD_WORK(:)*XSTD_WORK(:)  
!
        XSKEW_WORK(:) = XSKEW_WORK(:)/(ZSIZE(:)*XSTD_WORK(:)*XSTD_WORK(:)*XSTD_WORK(:))
!
  END WHERE
!
!----------------------------------------------------------------------------
!
END WHERE
!
WHERE (XMEAN_WORK(:)/=XUNDEF) 
  XMEAN_WORK(:) = AINT(XMEAN_WORK(:)) + &
              NINT((XMEAN_WORK(:)-AINT(XMEAN_WORK(:)))*100000000.)/100000000.
  XMIN_WORK(:) = AINT(XMIN_WORK(:)) + &
              NINT((XMIN_WORK(:)-AINT(XMIN_WORK(:)))*100000000.)/100000000.
  XMAX_WORK(:) = AINT(XMAX_WORK(:)) + &
              NINT((XMAX_WORK(:)-AINT(XMAX_WORK(:)))*100000000.)/100000000.
  XSTD_WORK(:) = AINT(XSTD_WORK(:)) + &
              NINT((XSTD_WORK(:)-AINT(XSTD_WORK(:)))*100000000.)/100000000.
  XSKEW_WORK(:) = AINT(XSKEW_WORK(:)) + &
              NINT((XSKEW_WORK(:)-AINT(XSKEW_WORK(:)))*100000000.)/100000000.
END WHERE
!
IF (LHOOK) CALL DR_HOOK('AVERAGE2_CTI',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE2_CTI
