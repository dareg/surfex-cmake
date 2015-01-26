!     #########
SUBROUTINE TRIP_SURFACE_WATER (KLISTING,PTSTEP,KGRCN,KSEQ,KNEXTX,KNEXTY,KSEQMAX, &
                               OPRINT,OMASK_VEL,PLEN,PSLOPEBED,PWIDTH,PN,PRUNOFF,&
                               PSURF_STO,PSURF_STO2,PGOUT,PSIN,PSOUT,PVEL,PHS,   &
                               PAREA,PQFR,PQRF,                                  &
                               PSSTO_ALL,PSSTO2_ALL,PSIN_ALL,PDRUN_ALL,          &
                               PSOUT_ALL,PVEL_ALL,PHS_ALL                        ) 
!     ################################################################
!
!!****  *TRIP_SURFACE_WATER*  
!!
!!    PURPOSE
!!    -------
!
!     Calculate the river storage in the next time step based on the storage of current time step 
!     Where OMASK_VEL=true the Manning equation is used to compute a variable flow velocity.
!
!     
!!**  METHOD
!!    ------
!
!     RK Ordre 4 Rang 4
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
!!	B. Decharme     
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/02/09 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODN_TRIP,     ONLY : XCVEL
USE MODD_TRIP_PAR, ONLY : XUNDEF, XM, XVELMIN, &
                          XHSMIN, XRHOLW
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABORT_TRIP
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER, INTENT(IN)                  :: KLISTING
!
REAL, INTENT(IN)                     :: PTSTEP ! Trip timestep value (10800s)
!
INTEGER, DIMENSION(:,:),INTENT(IN)   :: KGRCN  ! Flow direction (1->8)
INTEGER, DIMENSION(:,:),INTENT(IN)   :: KSEQ   ! River sequence
INTEGER, DIMENSION(:,:),INTENT(IN)   :: KNEXTX ! returns x and y point
INTEGER, DIMENSION(:,:),INTENT(IN)   :: KNEXTY ! of destination grid:
!                                                                    8 1 2
!                                                                    7   3
!                                                                    6 5 4
!
INTEGER, INTENT(IN)                    :: KSEQMAX ! maximum down flow
LOGICAL, INTENT(IN)                    :: OPRINT   !Printable budget key
!
LOGICAL, DIMENSION(:,:), INTENT(IN)    :: OMASK_VEL  !Variable velocity mask
REAL,    DIMENSION(:,:), INTENT(IN)    :: PLEN       ! river length       [m] 
REAL,    DIMENSION(:,:), INTENT(IN)    :: PSLOPEBED  ! river bed slopes             [m/m]
REAL,    DIMENSION(:,:), INTENT(IN)    :: PWIDTH     ! river widths                 [m]
REAL,    DIMENSION(:,:), INTENT(IN)    :: PN         ! Manning roughness coeficient [-] (0.03 to 0.065)
REAL,    DIMENSION(:,:), INTENT(IN)    :: PAREA      ! Grid-cell area    [m²]
REAL,    DIMENSION(:,:), INTENT(IN)    :: PRUNOFF    ! Surface runoff from ISBA    [kg/s]
REAL,    DIMENSION(:,:), INTENT(IN)    :: PGOUT      ! ground water outflow        [kg/s]
REAL,    DIMENSION(:,:), INTENT(IN)    :: PQFR       ! Flood flow to river         [kg/s]
REAL,    DIMENSION(:,:), INTENT(IN)    :: PQRF       ! River flow to floodplain    [kg/s]
REAL,    DIMENSION(:,:), INTENT(IN)    :: PSURF_STO  ! river channel storage at t  [kg]
REAL,    DIMENSION(:,:), INTENT(INOUT) :: PSURF_STO2 ! river channel storage at t+1[kg]
REAL,    DIMENSION(:,:), INTENT(OUT)   :: PHS   ! river channel height [m]
REAL,    DIMENSION(:,:), INTENT(OUT)   :: PSIN  ! Inflow to the surface river reservoir [kg/s]
REAL,    DIMENSION(:,:), INTENT(OUT)   :: PSOUT ! Outflow from the surface river reservoir [kg/s]
REAL,    DIMENSION(:,:), INTENT(OUT)   :: PVEL  ! River channel velocity  [m/s]
!
REAL,                    INTENT(OUT)   :: PSSTO_ALL,PSSTO2_ALL,PSIN_ALL,    &
                                          PDRUN_ALL,PSOUT_ALL,PVEL_ALL,     &
                                          PHS_ALL 
!                                         Final budget variable
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PLEN,1),SIZE(PLEN,2)) :: ZQIN
REAL, DIMENSION(SIZE(PLEN,1),SIZE(PLEN,2)) :: ZRADIUS
REAL, DIMENSION(SIZE(PLEN,1),SIZE(PLEN,2)) :: ZHS
REAL, DIMENSION(SIZE(PLEN,1),SIZE(PLEN,2)) :: ZVEL
REAL, DIMENSION(SIZE(PLEN,1),SIZE(PLEN,2)) :: ZRC
REAL, DIMENSION(SIZE(PLEN,1),SIZE(PLEN,2)) :: ZQOUT
REAL, DIMENSION(SIZE(PLEN,1),SIZE(PLEN,2)) :: ZSTOMAX
!
REAL    :: ZAREA
!
INTEGER :: ILON, ILAT, JLON, JLAT, ISEQ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
! * Init
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIP_SURFACE_WATER',0,ZHOOK_HANDLE)
!
ILON = SIZE(PLEN,1)
ILAT = SIZE(PLEN,2)
!
PSURF_STO2 (:,:) = 0.0
PSIN       (:,:) = 0.0
PSOUT      (:,:) = 0.0
PVEL       (:,:) = 0.0
PHS        (:,:) = XUNDEF
!
ZQIN       (:,:) = 0.0
ZRADIUS    (:,:) = 0.0
ZVEL       (:,:) = 0.0
ZHS        (:,:) = 0.0
ZRC        (:,:) = 0.0
ZQOUT      (:,:) = 0.0
ZSTOMAX    (:,:) = 0.0
!
!-------------------------------------------------------------------------------
! * River channel velocity             
!-------------------------------------------------------------------------------
!           
WHERE(OMASK_VEL(:,:))
    PHS     (:,:)=PSURF_STO(:,:)/(XRHOLW*PLEN(:,:)*PWIDTH(:,:))
    ZHS     (:,:)=MAX(XHSMIN,PHS(:,:))
    ZRADIUS (:,:)=LOG(PWIDTH(:,:)*ZHS(:,:)/(PWIDTH(:,:)+2.0*ZHS(:,:)))
    ZVEL    (:,:)=MAX(XVELMIN,EXP(XM*ZRADIUS(:,:))*SQRT(PSLOPEBED(:,:))/PN(:,:))
    PVEL    (:,:)=MIN(ZVEL(:,:),PLEN(:,:)/PTSTEP)
    ZRC     (:,:)=PVEL(:,:)/PLEN(:,:)
ELSEWHERE(KSEQ(:,:)>0)
    ZRC     (:,:)=XCVEL/PLEN(:,:)
ENDWHERE
!
!-------------------------------------------------------------------------------
! * Sequence loop
!-------------------------------------------------------------------------------
!
DO ISEQ=1,KSEQMAX
   DO JLAT=1,ILAT
      DO JLON=1,ILON
!      
        IF(KSEQ(JLON,JLAT)==ISEQ)THEN
!
!        ---------------------------------------------------------------------
!        inflow calculation
!
         ZQIN(JLON,JLAT)=ZQIN(JLON,JLAT)+PRUNOFF(JLON,JLAT)+PGOUT(JLON,JLAT)+PQFR(JLON,JLAT)-PQRF(JLON,JLAT)
         PSIN(JLON,JLAT)=ZQIN(JLON,JLAT)
!
!        ------------------------------------------------------------------
!        river channel storage calculation
!
         PSURF_STO2(JLON,JLAT) = PSURF_STO(JLON,JLAT)*EXP(-(ZRC(JLON,JLAT)*PTSTEP)) &
                               + (1.0-EXP(-(ZRC(JLON,JLAT)*PTSTEP)))*ZQIN(JLON,JLAT)&
                               / ZRC(JLON,JLAT)
!
!        -------------------------------------------------------------------
!        supress numerical artifacs
!
         ZSTOMAX(JLON,JLAT)=ZQIN(JLON,JLAT)*PTSTEP+PSURF_STO(JLON,JLAT)
!      
         PSURF_STO2(JLON,JLAT)=MIN(ZSTOMAX(JLON,JLAT),PSURF_STO2(JLON,JLAT))
!
!        ------------------------------------------------------------------
!        river channel outflow calculation and supress numerical artifacs
!
         ZQOUT(JLON,JLAT) = (PSURF_STO(JLON,JLAT)-PSURF_STO2(JLON,JLAT))/PTSTEP+ZQIN(JLON,JLAT)
         PSOUT(JLON,JLAT) = MAX(ZQOUT(JLON,JLAT),0.0)
!             
         PSURF_STO2(JLON,JLAT) = PSURF_STO2(JLON,JLAT) + (PSOUT(JLON,JLAT)-ZQOUT(JLON,JLAT))*PTSTEP        
!
!        ------------------------------------------------------------------
         IF(KGRCN(JLON,JLAT)>=1.AND.KGRCN(JLON,JLAT)<=8)THEN
           ZQIN(KNEXTX(JLON,JLAT),KNEXTY(JLON,JLAT))=ZQIN(KNEXTX(JLON,JLAT),KNEXTY(JLON,JLAT))+PSOUT(JLON,JLAT)
         ENDIF
!
        ENDIF
!
      ENDDO
   ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
! * Budget calculation
!-------------------------------------------------------------------------------
!
IF(OPRINT)THEN
!
  PDRUN_ALL   = 0.0
  PSSTO_ALL   = 0.0
  PSSTO2_ALL  = 0.0
  PSIN_ALL    = 0.0
  PSOUT_ALL   = 0.0
  PVEL_ALL    = 0.0
  PHS_ALL     = 0.0
  ZAREA       = 0.0
!
  DO JLAT=1,ILAT
     DO JLON=1,ILON
        IF(KSEQ(JLON,JLAT)>0)THEN
           PDRUN_ALL  = PDRUN_ALL  + PRUNOFF(JLON,JLAT)+PGOUT(JLON,JLAT)+PQFR(JLON,JLAT)-PQRF(JLON,JLAT)
           PSSTO_ALL  = PSSTO_ALL  + PSURF_STO (JLON,JLAT) / PAREA(JLON,JLAT)
           PSSTO2_ALL = PSSTO2_ALL + PSURF_STO2(JLON,JLAT) / PAREA(JLON,JLAT)
           PSIN_ALL   = PSIN_ALL   + ZQIN      (JLON,JLAT) / PAREA(JLON,JLAT)
           PSOUT_ALL  = PSOUT_ALL  + PSOUT     (JLON,JLAT) / PAREA(JLON,JLAT)
        ENDIF
        IF(OMASK_VEL(JLON,JLAT))THEN
         PVEL_ALL   = PVEL_ALL   + PVEL (JLON,JLAT) * PAREA(JLON,JLAT)
         PHS_ALL    = PHS_ALL    + PHS  (JLON,JLAT) * PAREA(JLON,JLAT)              
         ZAREA      = ZAREA      + PAREA(JLON,JLAT)         
        ENDIF
     ENDDO
  ENDDO
!
  IF(ZAREA>0.0)THEN
    PVEL_ALL = PVEL_ALL / ZAREA
    PHS_ALL  = PHS_ALL  / ZAREA
  ENDIF
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_SURFACE_WATER',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE TRIP_SURFACE_WATER
