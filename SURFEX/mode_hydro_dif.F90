!     ######spl
      MODULE MODE_HYDRO_DIF 
!     ################
!
!!****  *MODE_HYDRO_DIF * - pedo-transfert functions
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!    
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
!!	B. Decharme       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        04/2010
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
INTERFACE W_FUNC
  MODULE PROCEDURE W_FUNC_2D
  MODULE PROCEDURE W_FUNC_1D
END INTERFACE
!
INTERFACE PSI_FUNC
  MODULE PROCEDURE PSI_FUNC_2D
  MODULE PROCEDURE PSI_FUNC_1D
END INTERFACE
!
INTERFACE DPSI_FUNC
  MODULE PROCEDURE DPSI_FUNC
END INTERFACE
!
INTERFACE K_FUNC
  MODULE PROCEDURE K_FUNC_2D
  MODULE PROCEDURE K_FUNC_1D
END INTERFACE
!
INTERFACE DK_FUNC
  MODULE PROCEDURE DK_FUNC
END INTERFACE
!
INTERFACE WDRAIN_FUNC
  MODULE PROCEDURE WDRAIN_FUNC
END INTERFACE
!
INTERFACE INFMAX_FUNC
  MODULE PROCEDURE INFMAX_FUNC
END INTERFACE
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!Soil moisture using Brook and Corey or Van Genuchten
!-------------------------------------------------------------------------------
!
FUNCTION W_FUNC_2D(HDIF,PPSI,PWSAT,PMPOTSAT,PBCOEF,PALPHA,PN,PM,PWR)
USE MODD_ISBA_PAR, ONLY : XWGMIN
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
CHARACTER(LEN=*),     INTENT(IN) :: HDIF            !BC or VG DIF function
REAL, DIMENSION(:,:), INTENT(IN) :: PPSI,PWSAT      !Always used
REAL, DIMENSION(:,:), INTENT(IN) :: PMPOTSAT,PBCOEF !BC parameters
REAL, DIMENSION(:,:), INTENT(IN) :: PALPHA,PN,PM,PWR!VG parameters
!
REAL, DIMENSION(SIZE(PPSI,1),SIZE(PPSI,2)) :: W_FUNC_2D
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:W_FUNC_2D',0,ZHOOK_HANDLE)
IF(HDIF=='BC')THEN
  W_FUNC_2D(:,:) = PWSAT(:,:)*(MAX(1.0,PPSI(:,:)/PMPOTSAT(:,:))**(-1./PBCOEF(:,:)))
ELSE
  W_FUNC_2D(:,:) = PWR(:,:)+((PWSAT(:,:)-PWR(:,:))*(1.0/(1.0+(-PALPHA(:,:)*PPSI(:,:))**PN(:,:)))**PM(:,:))
  W_FUNC_2D(:,:) = MAX(MIN(W_FUNC_2D(:,:),PWSAT(:,:)),XWGMIN)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:W_FUNC_2D',1,ZHOOK_HANDLE)
END FUNCTION W_FUNC_2D
!
!-------------------------------------------------------------------------------
!
FUNCTION W_FUNC_1D(HDIF,PPSI,PWSAT,PMPOTSAT,PBCOEF,PALPHA,PN,PM,PWR)
USE MODD_ISBA_PAR, ONLY : XWGMIN
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
CHARACTER(LEN=*),   INTENT(IN) :: HDIF            !BC or VG DIF function
REAL, DIMENSION(:), INTENT(IN) :: PPSI,PWSAT      !Always used
REAL, DIMENSION(:), INTENT(IN) :: PMPOTSAT,PBCOEF !BC parameters
REAL, DIMENSION(:), INTENT(IN) :: PALPHA,PN,PM,PWR!VG parameters
!
REAL, DIMENSION(SIZE(PPSI))    :: W_FUNC_1D
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:W_FUNC_1D',0,ZHOOK_HANDLE)
IF(HDIF=='BC')THEN
  W_FUNC_1D(:) = PWSAT(:)*(MAX(1.0,PPSI(:)/PMPOTSAT(:))**(-1./PBCOEF(:))) 
ELSE
  W_FUNC_1D(:) = PWR(:)+((PWSAT(:)-PWR(:))*(1.0/(1.0+(-PALPHA(:)*PPSI(:))**PN(:)))**PM(:))
  W_FUNC_1D(:) = MAX(MIN(W_FUNC_1D(:),PWSAT(:)),XWGMIN)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:W_FUNC_1D',1,ZHOOK_HANDLE)
END FUNCTION W_FUNC_1D
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!Matric potential using Brook and Corey or Van Genuchten
!-------------------------------------------------------------------------------
!
FUNCTION PSI_FUNC_2D(HDIF,PWG,PWSAT,PMPOTSAT,PBCOEF,PALPHA,PN,PM,PWR)
USE MODD_ISBA_PAR, ONLY : XWGMIN
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
CHARACTER(LEN=*),     INTENT(IN) :: HDIF            !BC or VG DIF function
REAL, DIMENSION(:,:), INTENT(IN) :: PWG,PWSAT       !Always used
REAL, DIMENSION(:,:), INTENT(IN) :: PMPOTSAT,PBCOEF !BC parameters
REAL, DIMENSION(:,:), INTENT(IN) :: PALPHA,PN,PM,PWR!VG parameters
!
REAL, DIMENSION(SIZE(PM ,1),SIZE(PM ,2)) :: ZS
REAL, DIMENSION(SIZE(PWG,1),SIZE(PWG,2)) :: PSI_FUNC_2D
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:PSI_FUNC_2D',0,ZHOOK_HANDLE)
IF(HDIF=='BC')THEN
  PSI_FUNC_2D(:,:) = PMPOTSAT(:,:)*(MIN(1.0,PWG(:,:)/PWSAT(:,:))**(-PBCOEF(:,:))) 
ELSE
  ZS         (:,:) = MAX(XWGMIN,PWG(:,:)-PWR(:,:))/(PWSAT(:,:)-PWR(:,:))
  PSI_FUNC_2D(:,:) = -(1.0/PALPHA(:,:))*(MIN(1.0,ZS(:,:))**(-1.0/PM(:,:))-1.0)**(1.0/PN(:,:))
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:PSI_FUNC_2D',1,ZHOOK_HANDLE)
END FUNCTION PSI_FUNC_2D
!
!-------------------------------------------------------------------------------
!
FUNCTION PSI_FUNC_1D(HDIF,PWG,PWSAT,PMPOTSAT,PBCOEF,PALPHA,PN,PM,PWR)
USE MODD_ISBA_PAR, ONLY : XWGMIN
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
CHARACTER(LEN=*),   INTENT(IN) :: HDIF            !BC or VG DIF function
REAL, DIMENSION(:), INTENT(IN) :: PWG,PWSAT       !Always used
REAL, DIMENSION(:), INTENT(IN) :: PMPOTSAT,PBCOEF !BC parameters
REAL, DIMENSION(:), INTENT(IN) :: PALPHA,PN,PM,PWR!VG parameters
!
REAL, DIMENSION(SIZE(PM ))     :: ZS
REAL, DIMENSION(SIZE(PWG))     :: PSI_FUNC_1D
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:PSI_FUNC_1D',0,ZHOOK_HANDLE)
IF(HDIF=='BC')THEN
  PSI_FUNC_1D(:) = PMPOTSAT(:)*(MIN(1.0,PWG(:)/PWSAT(:))**(-PBCOEF(:))) 
ELSE
  ZS         (:) = MAX(XWGMIN,PWG(:)-PWR(:))/(PWSAT(:)-PWR(:))
  PSI_FUNC_1D(:) =-(1.0/PALPHA(:))*(MIN(1.0,ZS(:))**(-1.0/PM(:))-1.0)**(1.0/PN(:))
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:PSI_FUNC_1D',1,ZHOOK_HANDLE)
END FUNCTION PSI_FUNC_1D
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!Derrivative of Matric potential using Brook and Corey or Van Genuchten
!-------------------------------------------------------------------------------
!
FUNCTION DPSI_FUNC(HDIF,PPSI,PWG,PWSAT,PBCOEF,PN,PM,PWR)
USE MODD_ISBA_PAR, ONLY : XWGMIN
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
CHARACTER(LEN=*),   INTENT(IN) :: HDIF            !BC or VG DIF function
REAL, DIMENSION(:), INTENT(IN) :: PWG,PWSAT,PPSI  !Always used
REAL, DIMENSION(:), INTENT(IN) :: PBCOEF          !BC parameters
REAL, DIMENSION(:), INTENT(IN) :: PN,PM,PWR       !VG parameters
!
REAL, DIMENSION(SIZE(PM ))     :: ZWG
REAL, DIMENSION(SIZE(PM ))     :: ZS
REAL, DIMENSION(SIZE(PM ))     :: ZCOEF
REAL, DIMENSION(SIZE(PWG))     :: DPSI_FUNC
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:DPSI_FUNC',0,ZHOOK_HANDLE)
IF(HDIF=='BC')THEN
!        
  DPSI_FUNC(:) = -PBCOEF(:)*PPSI(:)/PWG(:)
!  
ELSE
!        
  ZWG      (:) = MAX(XWGMIN,PWG(:)-PWR(:))
  ZS       (:) = MIN(1.0,ZWG(:)/(PWSAT(:)-PWR(:)))
  ZCOEF    (:) = MAX(1.E-6,ZS(:)**(-1.0/PM(:))-1.0)
  DPSI_FUNC(:) = (ZCOEF(:)+1.0)*PPSI(:)/((1.0-PN(:))*ZCOEF(:)*ZWG(:))
!  
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:DPSI_FUNC',1,ZHOOK_HANDLE)
END FUNCTION DPSI_FUNC
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!Hydraulic conductivity using Brook and Corey or Van Genuchten
!-------------------------------------------------------------------------------
!
FUNCTION K_FUNC_2D(HDIF,PPSI,PFRZ,PCONDSAT,PMPOTSAT,PBCOEF,PALPHA,PN,PM,PL)
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
CHARACTER(LEN=*),     INTENT(IN) :: HDIF               !BC or VG DIF function
REAL, DIMENSION(:,:), INTENT(IN) :: PPSI,PFRZ,PCONDSAT !Always used
REAL, DIMENSION(:,:), INTENT(IN) :: PMPOTSAT,PBCOEF    !BC parameters
REAL, DIMENSION(:,:), INTENT(IN) :: PALPHA,PN,PM,PL    !VG parameters
!
REAL, DIMENSION(SIZE(PN  ,1),SIZE(PN  ,2)) :: ZP
REAL, DIMENSION(SIZE(PPSI,1),SIZE(PPSI,2)) :: K_FUNC_2D
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:K_FUNC_2D',0,ZHOOK_HANDLE)
IF(HDIF=='BC')THEN
  K_FUNC_2D(:,:) = PFRZ(:,:)*PCONDSAT(:,:)*(MAX(1.0,PPSI(:,:)/PMPOTSAT(:,:)) **(-(2.*PBCOEF(:,:)+3.)/PBCOEF(:,:)))
ELSE
  ZP       (:,:) = -PALPHA(:,:)*PPSI(:,:)
  K_FUNC_2D(:,:) = PFRZ(:,:)*PCONDSAT(:,:) * (((1.0+ZP(:,:)**PN(:,:))**PM(:,:)-ZP(:,:)**(PN(:,:)-1.0) )**2)     &
                                           / (( 1.0+ZP(:,:)**PN(:,:))**(PM(:,:)*(PL(:,:)+2.0)))
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:K_FUNC_2D',1,ZHOOK_HANDLE)
END FUNCTION K_FUNC_2D
!
!-------------------------------------------------------------------------------
!
FUNCTION K_FUNC_1D(HDIF,PPSI,PFRZ,PCONDSAT,PMPOTSAT,PBCOEF,PALPHA,PN,PM,PL)
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
CHARACTER(LEN=*),   INTENT(IN) :: HDIF               !BC or VG DIF function
REAL, DIMENSION(:), INTENT(IN) :: PPSI,PFRZ,PCONDSAT !Always used
REAL, DIMENSION(:), INTENT(IN) :: PMPOTSAT,PBCOEF    !BC parameters
REAL, DIMENSION(:), INTENT(IN) :: PALPHA,PN,PM,PL    !VG parameters
!
REAL, DIMENSION(SIZE(PN  ))    :: ZP
REAL, DIMENSION(SIZE(PPSI))    :: K_FUNC_1D
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:K_FUNC_1D',0,ZHOOK_HANDLE)
IF(HDIF=='BC')THEN
  K_FUNC_1D(:) = PFRZ(:)*PCONDSAT(:)*(MAX(1.0,PPSI(:)/PMPOTSAT(:)) **(-(2.*PBCOEF(:)+3.)/PBCOEF(:)))
ELSE
  ZP       (:) = -PALPHA(:)*PPSI(:)
  K_FUNC_1D(:) = PFRZ(:)*PCONDSAT(:) * (((1.0+ZP(:)**PN(:))**PM(:)-ZP(:)**(PN(:)-1.0) )**2)     &
                                     / (( 1.0+ZP(:)**PN(:))**(PM(:)*(PL(:)+2.0)))
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:K_FUNC_1D',1,ZHOOK_HANDLE)
END FUNCTION K_FUNC_1D
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!Derrivative of Matric potential using Brook and Corey or Van Genuchten
!-------------------------------------------------------------------------------
!
FUNCTION DK_FUNC(HDIF,PWG,PWSAT,PK1,PK2,PBCOEF,PM,PL,PWR)
USE MODD_ISBA_PAR, ONLY : XWGMIN
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
CHARACTER(LEN=*),   INTENT(IN) :: HDIF               !BC or VG DIF function
REAL, DIMENSION(:), INTENT(IN) :: PWG,PWSAT,PK1,PK2  !Always used
REAL, DIMENSION(:), INTENT(IN) :: PBCOEF             !BC parameters
REAL, DIMENSION(:), INTENT(IN) :: PM,PL,PWR          !VG parameters
!
REAL, DIMENSION(SIZE(PM ))     :: ZWG
REAL, DIMENSION(SIZE(PM ))     :: ZS
REAL, DIMENSION(SIZE(PM ))     :: ZXCOEF
REAL, DIMENSION(SIZE(PM ))     :: ZYCOEF
REAL, DIMENSION(SIZE(PWG))     :: DK_FUNC
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:DK_FUNC',0,ZHOOK_HANDLE)
IF(HDIF=='BC')THEN
!        
  DK_FUNC(:) = 0.5*(2.*PBCOEF(:)+3.)*SQRT(PK1(:)*PK2(:))/PWG(:)
!  
ELSE
!        
  ZWG(:)=MAX(XWGMIN,PWG(:)-PWR(:))
  ZS (:)=MIN(1.0,ZWG(:)/(PWSAT(:)-PWR(:)))
!  
  ZXCOEF(:)=1.0-ZS(:)**(1.0/PM(:))
  ZYCOEF(:)=1.0-ZXCOEF(:)**PM(:)
!  
  WHERE(ZXCOEF(:)==0.0.OR.ZYCOEF(:)==0.0)
!
    DK_FUNC(:) = 0.5*PL(:)*SQRT(PK1(:)*PK2(:))/ZWG(:)
!          
  ELSEWHERE
!          
    ZXCOEF(:)=(ZXCOEF(:)**(PM(:)-1.0))/ZYCOEF(:)
    ZXCOEF(:)=ZXCOEF(:)*ZS(:)**(1.0/PM(:))
!    
    DK_FUNC(:) = 0.5*(2.0*ZXCOEF(:)+PL(:))*SQRT(PK1(:)*PK2(:))/ZWG(:)
!
  ENDWHERE
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:DK_FUNC',1,ZHOOK_HANDLE)
END FUNCTION DK_FUNC
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!Linear (in time) sub-grid drainage term
!-------------------------------------------------------------------------------
!
FUNCTION WDRAIN_FUNC(HDIF,PWG,PWFC,PWSAT,PWLIM,PWDRAIN,PFRZ,PCONDSAT_EXP,PBCOEF,PM,PL,PWR)
USE MODD_ISBA_PAR, ONLY : XWGMIN
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
CHARACTER(LEN=*),     INTENT(IN) :: HDIF                      !BC or VG DIF function
REAL, DIMENSION(:,:), INTENT(IN) :: PWG,PWFC,PWSAT,PWLIM,    &
                                    PWDRAIN,PFRZ,PCONDSAT_EXP !Always used
REAL, DIMENSION(:,:), INTENT(IN) :: PBCOEF                    !BC parameters
REAL, DIMENSION(:,:), INTENT(IN) :: PM,PL,PWR                 !VG parameters
!
REAL, DIMENSION(SIZE(PWG,1),SIZE(PWG,2)) :: ZS
REAL, DIMENSION(SIZE(PWG,1),SIZE(PWG,2)) :: ZWDRAIN, ZTHETA
REAL, DIMENSION(SIZE(PWG,1),SIZE(PWG,2)) :: WDRAIN_FUNC
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:WDRAIN_FUNC',0,ZHOOK_HANDLE)
WDRAIN_FUNC(:,:) = 0.0
!
IF (ALL(PWDRAIN(:,:)==0.0) .AND. LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:WDRAIN_FUNC',1,ZHOOK_HANDLE)
IF(ALL(PWDRAIN(:,:)==0.0))RETURN
!
ZWDRAIN(:,:) = 0.0
ZTHETA (:,:) = 0.0
!
IF(HDIF=='BC')THEN  
  WHERE(PWDRAIN(:,:) > 0.0)
        ZS     (:,:) = (PWFC(:,:)+PWDRAIN(:,:))/PWSAT(:,:)
        ZTHETA (:,:) = MIN(1.0,ZS(:,:))
        ZWDRAIN(:,:) = PFRZ(:,:)*PCONDSAT_EXP(:,:) * (ZTHETA(:,:)**(2.*PBCOEF(:,:)+3.0))
  ENDWHERE
ELSE
  WHERE(PWDRAIN(:,:) > 0.0)
        ZS     (:,:) = MAX(XWGMIN,PWFC(:,:)+PWDRAIN(:,:)-PWR(:,:))/(PWSAT(:,:)-PWR(:,:))
        ZTHETA (:,:) = MIN(1.0,ZS(:,:))
        ZWDRAIN(:,:) = PFRZ(:,:)*PCONDSAT_EXP(:,:) * (ZTHETA(:,:)**PL(:,:))                           &
                                                   * (1.0-(1.0-ZTHETA(:,:)**(1.0/PM(:,:)))**PM(:,:))**2
  ENDWHERE        
ENDIF
!
WDRAIN_FUNC(:,:) = ZWDRAIN(:,:) * MAX(0.0, MIN(PWFC(:,:),PWG(:,:))-PWLIM)/(PWFC(:,:)-PWLIM)
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:WDRAIN_FUNC',1,ZHOOK_HANDLE)
END FUNCTION WDRAIN_FUNC
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!Green-Ampt approximation for maximum infiltration over ~20cm
!-------------------------------------------------------------------------------
!
FUNCTION INFMAX_FUNC(HDIF,PPSI,PFRZ,PCONDSAT,PMPOTSAT,PDZG,PDG)
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
CHARACTER(LEN=*),     INTENT(IN) :: HDIF                        !BC or VG DIF function
REAL, DIMENSION(:,:), INTENT(IN) :: PPSI,PFRZ,PCONDSAT,PDZG,PDG !Always used
REAL, DIMENSION(:,:), INTENT(IN) :: PMPOTSAT                    !BC parameters
REAL, DIMENSION(SIZE(PPSI,1),SIZE(PPSI,2)) :: ZPSI0
REAL, DIMENSION(SIZE(PPSI,1)) :: ZGREEN_AMPT, ZDEPTH
INTEGER                       :: II,JJ,INI,INLVLD
REAL, DIMENSION(SIZE(PPSI,1)) :: INFMAX_FUNC
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:INFMAX_FUNC',0,ZHOOK_HANDLE)
INI   =SIZE(PPSI,1)
INLVLD=SIZE(PPSI,2)
!
ZGREEN_AMPT(:) = 0.0
ZDEPTH     (:) = 0.0
!
!Saturated surface properties
IF(HDIF=='BC')THEN
   ZPSI0(:,:)=PMPOTSAT(:,:)
ELSE
   ZPSI0(:,:)=0.0
ENDIF
!
DO II=1,INI
   DO JJ=1,INLVLD
      ZGREEN_AMPT(II)=ZGREEN_AMPT(II)+PDZG(II,JJ)*PFRZ(II,JJ)*PCONDSAT(II,JJ)*(1.0+(ZPSI0(II,JJ)-PPSI(II,JJ))/PDZG(II,JJ))
      ZDEPTH     (II)=PDG(II,JJ)
      IF(ZDEPTH(II)>=0.2)EXIT
   ENDDO
ENDDO
!
INFMAX_FUNC(:) = ZGREEN_AMPT(:)/ZDEPTH(:)
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:INFMAX_FUNC',1,ZHOOK_HANDLE)
END FUNCTION INFMAX_FUNC
!
!-------------------------------------------------------------------------------
!
END MODULE MODE_HYDRO_DIF
