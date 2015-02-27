      SUBROUTINE GWF_CPL_UPDATE(PTABGW_H,PTABGW_F,OMASK_GW,PTOPO_RIV,     &
                                PELEV,PHGROUND,PHG_OLD,PWTD,PFWTD,PWTDELEV)
!     ##########################################################################
!
!!****  *FLOOD_UPDATE*  
!!
!!    PURPOSE
!!    -------
!
!     update groundwater diag
!     
!!**  METHOD
!!    ------
!
!     Direct calculation
!
!!    EXTERNAL
!!    --------
!
!     None
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
!!      B. Decharme
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/11/06 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,    DIMENSION(:,:,:), INTENT(IN)    :: PTABGW_H
REAL,    DIMENSION(:,:,:), INTENT(IN)    :: PTABGW_F
LOGICAL, DIMENSION(:,:),   INTENT(IN)    :: OMASK_GW
REAL,    DIMENSION(:,:),   INTENT(IN)    :: PTOPO_RIV
REAL,    DIMENSION(:,:),   INTENT(IN)    :: PELEV
REAL,    DIMENSION(:,:),   INTENT(IN)    :: PHGROUND
!
REAL,    DIMENSION(:,:),   INTENT(INOUT) :: PHG_OLD
REAL,    DIMENSION(:,:),   INTENT(OUT)   :: PWTD
REAL,    DIMENSION(:,:),   INTENT(OUT)   :: PFWTD
REAL,    DIMENSION(:,:),   INTENT(OUT)   :: PWTDELEV
!
!*      0.2    declarations of local variables
!
INTEGER, DIMENSION(SIZE(PTABGW_H,1),SIZE(PTABGW_H,2)) :: ISUP
INTEGER, DIMENSION(SIZE(PTABGW_H,1),SIZE(PTABGW_H,2)) :: IINF
REAL,    DIMENSION(SIZE(PTABGW_H,1),SIZE(PTABGW_H,2)) :: ZFWTD
!
INTEGER         :: ILON, ILAT, JLON, JLAT, JFRAC, INFRAC
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('GWF_CPL_UPDATE',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
! * Init
!-------------------------------------------------------------------------------
!
ILON   = SIZE(PTABGW_H,1)
ILAT   = SIZE(PTABGW_H,2)
INFRAC = SIZE(PTABGW_H,3)
!
!-------------------------------------------------------------------------------
! * Evolution of water table depth
!-------------------------------------------------------------------------------
!
WHERE(OMASK_GW(:,:).AND.PHGROUND(:,:)/=PHG_OLD(:,:))
      PWTD    (:,:) = PHGROUND(:,:)-PTOPO_RIV(:,:)
      PWTDELEV(:,:) = PHGROUND(:,:)-PELEV(:,:)
      PFWTD   (:,:) = PTABGW_F(:,:,1)
ENDWHERE
!
!-------------------------------------------------------------------------------
! * Evolution of the fraction of water table to rise
!-------------------------------------------------------------------------------
!
ISUP (:,:)=0
IINF (:,:)=0
ZFWTD(:,:)=PFWTD(:,:)
!
DO JLAT=1,ILAT
   DO JLON=1,ILON
      IF(OMASK_GW(JLON,JLAT).AND.PHGROUND(JLON,JLAT)/=PHG_OLD(JLON,JLAT))THEN
        DO JFRAC=1,INFRAC-1
           IF(PHGROUND(JLON,JLAT)>=PTABGW_H(JLON,JLAT,JFRAC))THEN
             ISUP(JLON,JLAT)=JFRAC+1
             IINF(JLON,JLAT)=JFRAC
           ENDIF           
        ENDDO
        IF(IINF(JLON,JLAT)>0)THEN
          ZFWTD(JLON,JLAT) =  PTABGW_F(JLON,JLAT,IINF(JLON,JLAT))                                       &
                           + (PHGROUND(JLON,JLAT                )-PTABGW_H(JLON,JLAT,IINF(JLON,JLAT))) &
                           * (PTABGW_F(JLON,JLAT,ISUP(JLON,JLAT))-PTABGW_F(JLON,JLAT,IINF(JLON,JLAT))) &
                           / (PTABGW_H(JLON,JLAT,ISUP(JLON,JLAT))-PTABGW_H(JLON,JLAT,IINF(JLON,JLAT)))
        ENDIF
        PFWTD(JLON,JLAT)=MIN(1.0,ZFWTD(JLON,JLAT))        
      ENDIF
   ENDDO
ENDDO
!
PHG_OLD(:,:)=PHGROUND(:,:)
!
IF (LHOOK) CALL DR_HOOK('GWF_CPL_UPDATE',1,ZHOOK_HANDLE)
!
END SUBROUTINE GWF_CPL_UPDATE

