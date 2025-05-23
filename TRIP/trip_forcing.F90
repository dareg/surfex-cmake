!     ######
SUBROUTINE TRIP_FORCING (TPG, &
                         KLUOUT,KLON,KLAT,KNB_TSTEP_RUN, &
                        PDRAIN,PRUNOFF,PSRC_FLOOD       )
!######################################################################
!
!!****  *TRIP_FORCING* - prepare the forcing for running trip
!!
!!    PURPOSE
!!    -------
!!
!!    AUTHOR
!!    ------
!!      B. decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2008 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_TRIP_GRID, ONLY : TRIP_GRID_t
!
USE MODN_TRIP_RUN, ONLY : CFILE_FRC,CREADFRC,CDRAIN, &
                          CRUNOFF,LCUMFRC,CSRC_FLOOD
!
!
USE MODD_TRIP_PAR,ONLY : XUNDEF
!
USE MODE_RW_TRIP
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(TRIP_GRID_t), INTENT(INOUT) :: TPG
!
INTEGER,                INTENT(IN)  :: KLUOUT
INTEGER,                INTENT(IN)  :: KLON
INTEGER,                INTENT(IN)  :: KLAT
INTEGER,                INTENT(IN)  :: KNB_TSTEP_RUN
!
REAL, DIMENSION(KLON,KLAT,KNB_TSTEP_RUN), INTENT(OUT) :: PDRAIN         
REAL, DIMENSION(KLON,KLAT,KNB_TSTEP_RUN), INTENT(OUT) :: PRUNOFF
REAL, DIMENSION(KLON,KLAT,KNB_TSTEP_RUN), INTENT(OUT) :: PSRC_FLOOD
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(KLON,KLAT,KNB_TSTEP_RUN) :: ZREAD_LATLON_DRAIN         
REAL, DIMENSION(KLON,KLAT,KNB_TSTEP_RUN) :: ZREAD_LATLON_RUNOFF
REAL, DIMENSION(KLON,KLAT,KNB_TSTEP_RUN) :: ZREAD_LATLON_SRC_FLOOD
!
REAL, DIMENSION(KLON*KLAT,KNB_TSTEP_RUN) :: ZREAD_VECTOR_DRAIN         
REAL, DIMENSION(KLON*KLAT,KNB_TSTEP_RUN) :: ZREAD_VECTOR_RUNOFF
REAL, DIMENSION(KLON*KLAT,KNB_TSTEP_RUN) :: ZREAD_VECTOR_SRC_FLOOD
!
 CHARACTER(LEN=6) :: YVAR
INTEGER          :: JLON,JLAT,JSTEP,ICOUNT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TRIP_FORCING',0,ZHOOK_HANDLE)
!
!*       1.     Initialize :
!               ------------
!
PDRAIN    (:,:,:) = 0.0
PRUNOFF   (:,:,:) = 0.0
PSRC_FLOOD(:,:,:) = 0.0
!
!-------------------------------------------------------------------------------
!
!*       2.     Get fields :
!               ------------
!
IF(CREADFRC=='LATLON')THEN
!
! * Lat lon case
!
  CALL READ_TRIP(KLUOUT,CFILE_FRC,CDRAIN ,ZREAD_LATLON_DRAIN (:,:,:))
  CALL READ_TRIP(KLUOUT,CFILE_FRC,CRUNOFF,ZREAD_LATLON_RUNOFF(:,:,:))
!
  WHERE(ZREAD_LATLON_DRAIN (:,:,:)==XUNDEF)ZREAD_LATLON_DRAIN (:,:,:)=0.0
  WHERE(ZREAD_LATLON_RUNOFF(:,:,:)==XUNDEF)ZREAD_LATLON_RUNOFF(:,:,:)=0.0
! 
  IF(LEN_TRIM(CSRC_FLOOD)/=0)THEN
    CALL READ_TRIP(KLUOUT,CFILE_FRC,CSRC_FLOOD,ZREAD_LATLON_SRC_FLOOD(:,:,:))
    WHERE(ZREAD_LATLON_SRC_FLOOD(:,:,:)==XUNDEF)ZREAD_LATLON_SRC_FLOOD(:,:,:)=0.0
  ELSE
    ZREAD_LATLON_SRC_FLOOD(:,:,:) = 0.0
  ENDIF
!
  IF(LCUMFRC)THEN
    DO JSTEP = KNB_TSTEP_RUN,2,-1
       ZREAD_LATLON_DRAIN    (:,:,JSTEP) = ZREAD_LATLON_DRAIN    (:,:,JSTEP) - ZREAD_LATLON_DRAIN    (:,:,JSTEP-1)
       ZREAD_LATLON_RUNOFF   (:,:,JSTEP) = ZREAD_LATLON_RUNOFF   (:,:,JSTEP) - ZREAD_LATLON_RUNOFF   (:,:,JSTEP-1)
       ZREAD_LATLON_SRC_FLOOD(:,:,JSTEP) = ZREAD_LATLON_SRC_FLOOD(:,:,JSTEP) - ZREAD_LATLON_SRC_FLOOD(:,:,JSTEP-1)
    ENDDO
  ENDIF
!
  PDRAIN    (:,:,:) = ZREAD_LATLON_DRAIN    (:,:,:)
  PRUNOFF   (:,:,:) = ZREAD_LATLON_RUNOFF   (:,:,:)
  PSRC_FLOOD(:,:,:) = ZREAD_LATLON_SRC_FLOOD(:,:,:)
!
ELSE
!        
! * Vector case
!
  CALL READ_TRIP(KLUOUT,CFILE_FRC,CDRAIN ,ZREAD_VECTOR_DRAIN (:,:))
  CALL READ_TRIP(KLUOUT,CFILE_FRC,CRUNOFF,ZREAD_VECTOR_RUNOFF(:,:))
!
  WHERE(ZREAD_VECTOR_DRAIN (:,:)==XUNDEF)ZREAD_VECTOR_DRAIN (:,:)=0.0
  WHERE(ZREAD_VECTOR_RUNOFF(:,:)==XUNDEF)ZREAD_VECTOR_RUNOFF(:,:)=0.0
!
  IF(LEN_TRIM(CSRC_FLOOD)/=0)THEN
    CALL READ_TRIP(KLUOUT,CFILE_FRC,CSRC_FLOOD,ZREAD_VECTOR_SRC_FLOOD(:,:))
    WHERE(ZREAD_VECTOR_SRC_FLOOD(:,:)==XUNDEF)ZREAD_VECTOR_SRC_FLOOD(:,:)=0.0
  ELSE
    ZREAD_VECTOR_SRC_FLOOD(:,:) = 0.0
  ENDIF
!
  IF(LCUMFRC)THEN
    DO JSTEP = KNB_TSTEP_RUN,2,-1
       ZREAD_VECTOR_DRAIN    (:,JSTEP) = ZREAD_VECTOR_DRAIN    (:,JSTEP) - ZREAD_VECTOR_DRAIN    (:,JSTEP-1)
       ZREAD_VECTOR_RUNOFF   (:,JSTEP) = ZREAD_VECTOR_RUNOFF   (:,JSTEP) - ZREAD_VECTOR_RUNOFF   (:,JSTEP-1)
       ZREAD_VECTOR_SRC_FLOOD(:,JSTEP) = ZREAD_VECTOR_SRC_FLOOD(:,JSTEP) - ZREAD_VECTOR_SRC_FLOOD(:,JSTEP-1)
    ENDDO
  ENDIF
!  
  ICOUNT=0
  DO JLAT=1,KLAT
     DO JLON=1,KLON
        ICOUNT=ICOUNT+1
        PDRAIN    (JLON,JLAT,:)= ZREAD_VECTOR_DRAIN    (ICOUNT,:)
        PRUNOFF   (JLON,JLAT,:)= ZREAD_VECTOR_RUNOFF   (ICOUNT,:)
        PSRC_FLOOD(JLON,JLAT,:)= ZREAD_VECTOR_SRC_FLOOD(ICOUNT,:)
     ENDDO
  ENDDO
!
ENDIF
!
DO JSTEP=1,KNB_TSTEP_RUN
   PDRAIN    (:,:,JSTEP) = PDRAIN    (:,:,JSTEP)*TPG%XAREA(:,:)
   PRUNOFF   (:,:,JSTEP) = PRUNOFF   (:,:,JSTEP)*TPG%XAREA(:,:)
   PSRC_FLOOD(:,:,JSTEP) = PSRC_FLOOD(:,:,JSTEP)*TPG%XAREA(:,:)
ENDDO
!
IF (LHOOK) CALL DR_HOOK('TRIP_FORCING',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIP_FORCING
