!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE WRITESURF_SSO_n (HSELECT, USS, HPROGRAM)
!     #################################
!
!!****  *WRITESURF_SSO_n* - writes orographic fields
!!
!!    PURPOSE
!!    -------
!!       
!!
!!
!!**  METHOD
!!    ------
!!      
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_PGDWORK, ONLY : LORORAD
!
USE MODI_WRITE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
 CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT
!
TYPE(SSO_t), INTENT(INOUT) :: USS
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
INTEGER           :: JK             ! Loop index
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.     Orography :
!               ---------
!
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_SSO_N',0,ZHOOK_HANDLE)
YRECFM='AVG_ZS'
YCOMMENT='X_Y_AVG_ZS (M)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XAVG_ZS(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='SIL_ZS'
YCOMMENT='X_Y_SIL_ZS (M)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XSIL_ZS(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!*       2.     Subgrid Orography :
!               -----------------
!
YRECFM='SSO_STDEV'
YCOMMENT='X_Y_SSO_STDEV (M)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XSSO_STDEV(:),IRESP,HCOMMENT=YCOMMENT)
!
!
YRECFM='MIN_ZS'
YCOMMENT='X_Y_MIN_ZS (M)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XMIN_ZS(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='MAX_ZS'
YCOMMENT='X_Y_MAX_ZS (M)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XMAX_ZS(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='SSO_ANIS'
YCOMMENT='X_Y_SSO_ANIS (-)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XSSO_ANIS(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='SSO_DIR'
YCOMMENT='X_Y_SSO_DIR (DEG)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XSSO_DIR(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='SSO_SLOPE'
YCOMMENT='X_Y_SSO_SLOPE (-)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XSSO_SLOPE(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!*       3.     Subgrid Orography roughness:
!               ---------------------------
!
!
YRECFM='HO2IP'
YCOMMENT='X_Y_HO2IP (M)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XHO2IP(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='HO2JP'
YCOMMENT='X_Y_HO2JP (M)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XHO2JP(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='HO2IM'
YCOMMENT='X_Y_HO2IM (M)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XHO2IM(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='HO2JM'
YCOMMENT='X_Y_HO2JM (M)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XHO2JM(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='AOSIP'
YCOMMENT='X_Y_AOSIP (-)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XAOSIP(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='AOSJP'
YCOMMENT='X_Y_AOSJP (-)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XAOSJP(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='AOSIM'
YCOMMENT='X_Y_AOSIM (-)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XAOSIM(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='AOSJM'
YCOMMENT='X_Y_AOSJM (-)'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,USS%XAOSJM(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!*       4.     Orographic radiation parameters:
!               -------------------------------
!
YCOMMENT='(-)'
  CALL WRITE_SURF(HSELECT, &
                  HPROGRAM,'LORORAD',LORORAD,IRESP,YCOMMENT)
!
IF (USS%NSECTORS > 0)THEN
  YCOMMENT='(-)'
  CALL WRITE_SURF(HSELECT, &
                  HPROGRAM,'ASPSECTORS',USS%NSECTORS,IRESP,YCOMMENT)

  YRECFM='AVG_SLO'
  YCOMMENT='X_Y_AVG_SLO (RAD)'
  CALL WRITE_SURF(HSELECT, &
                  HPROGRAM,YRECFM,USS%XAVG_SLO(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='SLOPE'
  YCOMMENT='X_Y_SLOPE (RAD)'
  CALL WRITE_SURF(HSELECT, &
                  HPROGRAM,YRECFM,USS%XSLOPE(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='ASPECT'
  YCOMMENT='X_Y_ASPECT (RAD)'
  CALL WRITE_SURF(HSELECT, &
                  HPROGRAM,YRECFM,USS%XASPECT(:),IRESP,HCOMMENT=YCOMMENT)
  !
  DO JK=1,USS%NSECTORS
    WRITE(YRECFM,"(a9,i0.2)") 'SLOPE_DIR',JK
    WRITE(YCOMMENT,"(a13,i0.2,a6)") 'X_Y_SLOPE_DIR',JK,' (RAD)'
    CALL WRITE_SURF(HSELECT, &
                    HPROGRAM,YRECFM,USS%XSLOPE_DIR(:,JK),IRESP,HCOMMENT=YCOMMENT)
  END DO
  !
  DO JK=1,USS%NSECTORS
    WRITE(YRECFM,"(a8,i0.2)") 'FRAC_DIR',JK
    WRITE(YCOMMENT,"(a12,i0.2,a4)") 'X_Y_FRAC_DIR',JK,' (%)'
    CALL WRITE_SURF(HSELECT, &
                    HPROGRAM,YRECFM,USS%XFRAC_DIR(:,JK),IRESP,HCOMMENT=YCOMMENT)
  END DO
  !
  YRECFM='SVF'
  YCOMMENT='X_Y_SVF (-)'
  CALL WRITE_SURF(HSELECT, &
                  HPROGRAM,YRECFM,USS%XSVF(:),IRESP,HCOMMENT=YCOMMENT)
  !
  DO JK=1,USS%NSECTORS
    WRITE(YRECFM,"(a9,i0.2)") 'HMINS_DIR',JK
    WRITE(YCOMMENT,"(a13,i0.2,a4)") 'X_Y_HMINS_DIR',JK,' (-)'
    CALL WRITE_SURF(HSELECT, &
                    HPROGRAM,YRECFM,USS%XHMINS_DIR(:,JK),IRESP,HCOMMENT=YCOMMENT)
  END DO
  !
  DO JK=1,USS%NSECTORS
    WRITE(YRECFM,"(a9,i0.2)") 'HMAXS_DIR',JK
    WRITE(YCOMMENT,"(a13,i0.2,a4)") 'X_Y_HMAXS_DIR',JK,' (-)'
    CALL WRITE_SURF(HSELECT, &
                    HPROGRAM,YRECFM,USS%XHMAXS_DIR(:,JK),IRESP,HCOMMENT=YCOMMENT)
  END DO
END IF
IF (LHOOK) CALL DR_HOOK('WRITESURF_SSO_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_SSO_n
