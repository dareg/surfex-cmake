!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######spl
SUBROUTINE OL_TIME_INTERP_ATM (KSURF_STEP,KNB_ATM,                       &
                               PTA1,PTA2,PQA1,PQA2,PWIND1,PWIND2,        &
                               PDIR_SW1,PDIR_SW2,PSCA_SW1,PSCA_SW2,      &
                               PLW1,PLW2,PSNOW2,PRAIN2,                  &
                               PPS1,PPS2,PCO21,PCO22,PDIR1,PDIR2,        &
                               PO31,PO32,PAE1,PAE2,PIMPWET2,PIMPDRY2,    &
                               PZEN,PSUMZEN )  
!**************************************************************************
!
!!    PURPOSE
!!    -------
!        Time interpolation of the atmospheric forcing
!        So far, it is a simple linear interpolation.
!        More complex interpolation may be added, especially for the atmospheric
!        radiation (Option to use).
!        Output are in the module 
!!
!!**  METHOD
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
!!
!!    AUTHOR
!!    ------
!!      F. Habets   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2003
!
USE MODN_IO_OFFLINE, ONLY : LINTERP_SW
!
USE MODD_CSTS,       ONLY : XPI, XRD, XRV, XG
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODN_IO_OFFLINE,  ONLY : NIMPUROF, LFORCIMP, LFORCATMOTARTES
USE MODD_FORC_ATM,  ONLY: XTA         ,&! air temperature forcing               (K)
                            XQA       ,&! air specific humidity forcing         (kg/m3)
                            XRHOA     ,&! air density forcing                   (kg/m3)
                            XZS       ,&! orography                             (m)
                            XU        ,&! zonal wind                            (m/s)
                            XV        ,&! meridian wind                         (m/s)
                            XDIR_SW   ,&! direct  solar radiation (on horizontal surf.)
                            XSCA_SW   ,&! diffuse solar radiation (on horizontal surf.)
                            XLW       ,&! longwave radiation (on horizontal surf.)
                            XPS       ,&! pressure at atmospheric model surface (Pa)
                            XPA       ,&! pressure at forcing level             (Pa)
                            XRHOA     ,&! density at forcing level              (kg/m3)
                            XCO2      ,&! CO2 concentration in the air          (kg/kg)
                            XO3       ,&! Ozone
                            XAE       ,&! Aerosol optical depth
                            XIMPWET   ,&! wet deposit coefficient
                            XIMPDRY   ,&! dry deposit coefficient
                            XSNOW     ,&! snow precipitation                    (kg/m2/s)
                            XRAIN     ,&! liquid precipitation                  (kg/m2/s)
                            XZREF       ! height of T,q forcing                 (m)  
!
USE MODI_GET_LUOUT
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
!$ USE OMP_LIB
!
IMPLICIT NONE
!
!
! global variables
INTEGER,INTENT(IN) :: KSURF_STEP, KNB_ATM
REAL, DIMENSION(:),INTENT(IN) :: PTA1,PTA2,PQA1,PQA2,PWIND1,PWIND2
REAL, DIMENSION(:),INTENT(IN) :: PDIR_SW1,PDIR_SW2,PSCA_SW1,PSCA_SW2,PLW1,PLW2
REAL, DIMENSION(:),INTENT(IN) :: PSNOW2,PRAIN2,PPS1,PPS2,PCO21,PCO22,PDIR1,PDIR2,PO31,PO32,PAE1,PAE2
REAL, DIMENSION(:,:),INTENT(IN) :: PIMPWET2,PIMPDRY2
REAL, DIMENSION(:),INTENT(IN) :: PZEN,PSUMZEN 

! local variables
REAL :: ZDTA, ZDQA, ZDDIR_SW, ZDSCA_SW, ZDLW,  &
        ZDPS, ZDCO2, ZDU, ZDV, ZU1, ZV1, ZU2, ZV2,ZDO3,ZDAE  
REAL :: ZPI, ZNB_ATM, ZSURF_STEP, ZCOEF, ZCOEF2
INTEGER :: J,JIMP
INTEGER :: ILUOUT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_OMP
!========================================================================
!
IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_1',0,ZHOOK_HANDLE)
!
 CALL GET_LUOUT('OFFLIN',ILUOUT)
!
ZPI = XPI/180.
ZNB_ATM = KNB_ATM*1.
ZSURF_STEP = KSURF_STEP*1.-1.
ZCOEF = ZSURF_STEP / ZNB_ATM
!
IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_1',1,ZHOOK_HANDLE)
!
!$OMP PARALLEL PRIVATE(ZHOOK_HANDLE_OMP)
IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_2',0,ZHOOK_HANDLE_OMP)
!$OMP DO PRIVATE(J,ZU1,ZU2,ZV1,ZV2,ZDU,ZDV,ZDTA, &
!$OMP ZDQA,ZDLW,ZDPS,ZDCO2,ZDDIR_SW,ZDSCA_SW,ZCOEF2)
DO J = 1,SIZE(PTA1)
  !
  IF (PTA1(J)/=XUNDEF) THEN
    ! 
    ! Compute wind components
    !
    ! zonal wind
    ZU1 = PWIND1(J) * SIN(PDIR1(J)*ZPI)
    ZU2 = PWIND2(J) * SIN(PDIR2(J)*ZPI)
    ZDU = (ZU2-ZU1)*ZCOEF
    XU(J) = ZU1 + ZDU
    !
    ZV1 = PWIND1(J) * COS(PDIR1(J)*ZPI) 
    ZV2 = PWIND2(J) * COS(PDIR2(J)*ZPI)
    ZDV = (ZV2-ZV1)*ZCOEF
    XV(J) = ZV1 + ZDV
    !    
    ! Compute variation from atmospheric time step J and J+1
    !
    ZDTA     = (PTA2(J)-PTA1(J))*ZCOEF
    XTA(J) = PTA1(J) + ZDTA
    !
    ZDQA     = (PQA2(J)-PQA1(J))*ZCOEF
    XQA(J) = PQA1(J) + ZDQA
    !
    ZDLW     = (PLW2(J)-PLW1(J))*ZCOEF
    XLW(J) = PLW1(J) + ZDLW
    !
    ZDPS     = (PPS2(J)-PPS1(J))*ZCOEF
    XPS(J) = PPS1(J) + ZDPS
    !
    ZDCO2    = (PCO22(J)-PCO21(J))*ZCOEF
    XCO2(J) = PCO21(J) + ZDCO2
    !
    IF (LFORCATMOTARTES) THEN
      ZDO3    = (PO32(J)-PO31(J))*ZCOEF
      XO3(J) = PO31(J) + ZDO3
      
      ZDAE    = (PAE2(J)-PAE1(J))*ZCOEF
      XAE(J) = PAE1(J) + ZDAE
    ENDIF
    !
    IF (LFORCIMP) THEN
      DO JIMP=1,NIMPUROF
        XIMPWET(J,JIMP) = PIMPWET2(J,JIMP)
        
        XIMPDRY(J,JIMP) = PIMPDRY2(J,JIMP)
      
      ENDDO
    ENDIF
    !
    IF (LINTERP_SW) THEN
      !
      ZCOEF2=0.
      IF (PSUMZEN(J)>0.) ZCOEF2=MAX((COS(PZEN(J))/PSUMZEN(J)),0.)
      !
      XDIR_SW(J,1) = MIN(PDIR_SW2(J)*ZCOEF2,1300.0*MAX(COS(PZEN(J)),0.))
      !
      XSCA_SW(J,1) = MIN(PSCA_SW2(J)*ZCOEF2,1300.0*MAX(COS(PZEN(J)),0.))
      !
    ELSE
      !
      ZDDIR_SW = (PDIR_SW2(J)-PDIR_SW1(J))*ZCOEF
      XDIR_SW(J,1) = PDIR_SW1(J)+ZDDIR_SW
      !
      ZDSCA_SW = (PSCA_SW2(J)-PSCA_SW1(J))*ZCOEF
      XSCA_SW(J,1) = PSCA_SW1(J)+ZDSCA_SW
      !
    ENDIF
    !
    !
    XRAIN (J)= PRAIN2(J)
    XSNOW (J)= PSNOW2(J)
    !
    !
    XRHOA(J) = XPS(J) / ( XTA(J)*XRD * ( 1.+((XRV/XRD)-1.)*XQA(J) ) + XZREF(J)*XG )
    !
    ! humidity in kg/m3
    XQA(J) = XQA(J) * XRHOA(J)
    !
  ENDIF
  !
ENDDO
!$OMP END DO
IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_2',1,ZHOOK_HANDLE_OMP)
!$OMP END PARALLEL
!
IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_3',0,ZHOOK_HANDLE)
! air density
!
! Check No value data
!---------------------
! Error cases
!
IF ((MINVAL(XTA)  .EQ.XUNDEF).OR.(MINVAL(XQA).EQ.XUNDEF).OR.&
      (MINVAL(XU).EQ.XUNDEF).OR.(MINVAL(XRAIN).EQ.XUNDEF).OR.&
      (MINVAL(XSNOW).EQ.XUNDEF)) THEN  
    WRITE(ILUOUT,*)'MINVAL(XTA),MINVAL(XQA),MINVAL(XU),MINVAL(XRAIN),MINVAL(XSNOW)'
    WRITE(ILUOUT,*)MINVAL(XTA),MINVAL(XQA),MINVAL(XU),MINVAL(XRAIN),MINVAL(XSNOW)
    CALL ABOR1_SFX('OL_TIME_INTERP_ATM: UNDEFINED VALUE IN ATMOSPHERIC FORCING')
ENDIF
!
IF ((MINVAL(XDIR_SW).EQ.XUNDEF).AND.(MINVAL(XSCA_SW).EQ.XUNDEF)) THEN
    WRITE(ILUOUT,*)'MINVAL(XSCA_SW),MINVAL(XDIR_SW)'
    WRITE(ILUOUT,*)MINVAL(XSCA_SW),MINVAL(XDIR_SW)
    CALL ABOR1_SFX('OL_TIME_INTERP_ATM: UNDEFINED VALUE IN ATMOSPHERIC FORCING')
ENDIF
!
IF ((MINVAL(XPS).EQ.XUNDEF).AND.(MINVAL(XZS).EQ.XUNDEF)) THEN
    WRITE(ILUOUT,*)'MINVAL(XPS),MINVAL(XZS)'
    WRITE(ILUOUT,*)MINVAL(XPS),MINVAL(XZS)
    CALL ABOR1_SFX('OL_TIME_INTERP_ATM: UNDEFINED VALUE IN ATMOSPHERIC FORCING')
ENDIF
!
IF (MINVAL(XDIR_SW).EQ.XUNDEF) XDIR_SW(:,:)=0. ! No direct solar radiation
IF (MINVAL(XSCA_SW).EQ.XUNDEF) XSCA_SW(:,:)=0. ! No diffuse solar radiation
IF (MINVAL(XPS)    .EQ.XUNDEF) THEN            ! No surface Pressure 
   WRITE(ILUOUT,*)' OL_TIME_INTERP_ATM: SURFACE PRESSURE COMPUTED FROM ZS'
   XPS(:)  = 101325*(1-0.0065 * XZS(:)/288.15)**5.31
ENDIF
!
!* forcing level pressure from hydrostatism
WHERE(XPS(:)/=XUNDEF)
  XPA(:) = XPS(:) - XRHOA(:) * XZREF(:) * XG
ENDWHERE
!
IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_3',1,ZHOOK_HANDLE)
!
END SUBROUTINE OL_TIME_INTERP_ATM
