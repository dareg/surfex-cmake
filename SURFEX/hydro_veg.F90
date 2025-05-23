!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE HYDRO_VEG(HRAIN, PTSTEP, PMUF, PRR, PLEV, PLETR,    &
                              PVEG, PPSNV, PWR, PWRMAX, PPG, PDRIP,    &
                              PRRVEG, PLVTT  )    
!     #####################################################################
!
!!****  *HYDRO_VEG*  
!!
!!    PURPOSE
!!    -------
!
!     Calculates the evolution of the liquid water retained in the vegetation 
!     canopy (Wr). Also determine the runoff from the canopy that reaches the
!     ground (Mahfouf et al. 1995). This routine take into account the spatially 
!     exponential distribution of precip introduced by Entekhabi and Eagleson (1989).
!         
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    USE MODD_CST
!!      
!!    REFERENCE
!!    ---------
!!
!!    Noilhan and Planton (1989)
!!    Belair (1995)
!!    Mahfouf et al. 1995
!!    Decharme and Douville (2006)
!!      
!!    AUTHOR
!!    ------
!!
!!      S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    14/03/95 
!!                  31/08/98 (V. Masson and F. Habets) add Dumenil et Todini
!!                           runoff scheme
!!                  31/08/98 (V. Masson and A. Boone) add the third soil-water
!!                           reservoir (WG3,D3)
!!                  31/05/04 (B. Decharme) add the rainfall spatial distribution
!!                      2008 (B. Decharme) add the dripping rate as new diag
!!                  11/2009  (S.Senesi) returns precipitation intercepted by  
!                               the vegetation 
!!                  07/2011  (B. Decharme) delete SGH for very fine precipitation
!!                  09/2012  (B. Decharme) Computation efficiency for HRAIN=='SGH'
!!                  10/2012  (B. Decharme) PPG intent(out)
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SGH_PAR, ONLY : X001
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRD
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=*),     INTENT(IN)   :: HRAIN   ! Rainfall spatial distribution
                                              ! 'DEF' = No rainfall spatial distribution
                                              ! 'SGH' = Rainfall exponential spatial distribution
                                              ! 
!
REAL, INTENT(IN)                    :: PTSTEP
!                                      timestep of the integration
!
REAL, DIMENSION(:), INTENT(IN)    :: PRR,  PLEV, PLETR, PMUF, PLVTT
!                                      PRR   = rain rate
!                                      PLEV = latent heat of evaporation over vegetation
!                                      PLETR = evapotranspiration of the vegetation
!                                      PMUF   = fraction of the grid cell reached by the precipitation
!                                      PLVTT  = latent heat of vaporization (J/kg)
!
REAL, DIMENSION(:), INTENT(IN)    :: PVEG, PWRMAX
!                                      PVEG   = fraction of vegetation
!                                      PWRMAX = maximum equivalent water content
!                                               in the vegetation canopy
!
REAL, DIMENSION(:), INTENT(IN)    :: PPSNV
!                                      PPSNV = vegetation covered by snow
!
REAL, DIMENSION(:), INTENT(INOUT) :: PWR
!                                      PWR = liquid water retained on the foliage
!                                             of the vegetation at time 't+dt'
!
REAL, DIMENSION(:), INTENT(OUT)   :: PPG,PDRIP
!                                      PPG   = total water reaching the ground
!                                      PDRIP = Dripping from the vegetation
REAL, DIMENSION(:), INTENT(OUT)   :: PRRVEG  
!                                      PRRVEG = Precip. intercepted by vegetation (kg/m2/s)
!
!
!*      0.2    declarations of local variables
!
#ifdef HIRLAM_SP_HACKS
REAL(KIND=JPRD), DIMENSION(SIZE(PVEG)) :: ZER
!                                  ZER = evaporation rate from the canopy
!
REAL(KIND=JPRD), DIMENSION(SIZE(PVEG)) :: ZWR ! for time stability scheme
!
REAL(KIND=JPRD), DIMENSION(SIZE(PVEG)) :: ZRUIR, ZRUIR2 ! dripping from the vegetation
!
REAL(KIND=JPRD)                        :: ZLIM
REAL(KIND=JPRD), DIMENSION(SIZE(PVEG)) :: PWR_DP, PPG_DP, PDRIP_DP, PRRVEG_DP
#else
REAL, DIMENSION(SIZE(PVEG)) :: ZER
!                                  ZER = evaporation rate from the canopy
!
REAL, DIMENSION(SIZE(PVEG)) :: ZWR ! for time stability scheme
!
REAL, DIMENSION(SIZE(PVEG)) :: ZRUIR, ZRUIR2 ! dripping from the vegetation
!
REAL                        :: ZLIM
REAL, DIMENSION(SIZE(PVEG)) :: PWR_DP, PPG_DP, PDRIP_DP, PRRVEG_DP
#endif
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('HYDRO_VEG',0,ZHOOK_HANDLE)
ZRUIR (:) = 0.
ZRUIR2(:) = 0.
PDRIP_DP (:) = 0.
ZWR   (:) = 0.
!
!*       1.     EVOLUTION OF THE EQUIVALENT WATER CONTENT Wr
!               --------------------------------------------
!
!evaporation rates
!
ZER(:)    = (PLEV(:)-PLETR(:))  / PLVTT(:)
!
!intercepted rainfall rate
!
PRRVEG_DP(:) = PVEG(:) * (1.-PPSNV(:)) * PRR(:)
!
!evolution of the intercepted water
!(if we don't consider the runoff)
!
PWR_DP(:)  = PWR(:) - PTSTEP * (ZER(:) - PRRVEG_DP(:))
!
!When Wr < 0, the direct evaporation
!(i.e., EV-ETR) removes too much
!liquid water from the vegetation
!reservoir.  This is considered as
!negative runoff, and it is stocked
!in ZRUIR2.
!
ZRUIR2(:) = MIN(0.,PWR_DP(:)/PTSTEP)
!
!Wr must be positive
!
PWR_DP(:)    = MAX(0., PWR_DP(:))
!
IF(HRAIN=='SGH')THEN
!        
!*       2.     SPATIALLY EXPONENTIAL DISTRIBUTION OF PRECIPITATION
!               ---------------------------------------------------
!
!
!  Subgrid dripping from Wr
!
   ZLIM=X001/PTSTEP
!
   WHERE(PRRVEG_DP(:)>ZLIM.AND.PWR_DP(:)>0.0)
        ZRUIR(:) = PRRVEG_DP(:)*EXP(PMUF(:)*(PWR_DP(:)-PWRMAX(:))/(PRRVEG_DP(:)*PTSTEP))
        ZRUIR(:) = MIN(ZRUIR(:),PWR_DP(:)/PTSTEP) 
   ENDWHERE
!
   IF(PTSTEP>300.)THEN
!
!    if the isba time step is coarser than 5min, the "prediction/correction" method is applied
!    to Wr using the predicted Wr* at the end of the time step for time numerical stability
!
     ZWR(:)   = PWR_DP(:)-PTSTEP*ZRUIR(:)
     ZRUIR(:) = 0.0
!
!    if the dripping is too big, the "prediction/correction" method is applied to Wr using
!    the predicted Wr* at the midle of the time step for time numerical stability 
!    (<=> Runge-Kutta order 1 rang 1)
!
     WHERE(PRRVEG_DP(:)>ZLIM.AND.ZWR(:)<=0.0)
           ZRUIR(:) = PRRVEG_DP(:)*EXP(PMUF(:)*(PWR_DP(:)-PWRMAX(:))/(PRRVEG_DP(:)*PTSTEP/2.))
           ZRUIR(:) = MIN(ZRUIR(:),PWR_DP(:)/(PTSTEP/2.))
           ZWR  (:) = PWR_DP(:)-PTSTEP*ZRUIR(:)/2.
           ZRUIR(:) = 0.0
     ENDWHERE
!
!    Calculate the corrected dripping from the predicted Wr*
!
     WHERE(PRRVEG_DP(:)>ZLIM.AND.ZWR(:)>0.0)
          ZRUIR(:) = PRRVEG_DP(:)*EXP(PMUF(:)*(ZWR(:)-PWRMAX(:))/(PRRVEG_DP(:)*PTSTEP))
          ZRUIR(:) = MIN(ZRUIR(:),PWR_DP(:)/PTSTEP) 
     ENDWHERE
!
   ENDIF
!
   PWR_DP(:)   = PWR_DP(:)-PTSTEP*ZRUIR(:)  
!
!  As previously Wr must be positive (numerical artefact)
!
   ZRUIR2(:) = ZRUIR2(:) + MIN(0.,PWR_DP(:)/PTSTEP)
   PWR_DP(:)    = MAX( 0., PWR_DP(:) )
!
!  Wr must be smaller then Wrmax
!  Then if Wr remain > Wrmax, there is runoff
!
   ZRUIR(:) = ZRUIR(:) + MAX(0., (PWR_DP(:) - PWRMAX(:)) / PTSTEP )
!   
ELSE
!
! if Wr > Wrmax, there is runoff
!
  ZRUIR(:) = MAX(0., (PWR_DP(:) - PWRMAX(:)) / PTSTEP )
!
ENDIF
!
!Wr must be smaller then Wrmax
!
PWR_DP(:)   = MIN(PWR_DP(:), PWRMAX(:))
!
!
!*       3.     LIQUID WATER REACHING THE GROUND Pg
!               -----------------------------------
!
!Thus, the rate of liquid water reaching the ground is the 
!precipitation plus the vegetation runoff (we also consider the
!negative runoff).
!
PPG_DP(:) = (1.-PVEG(:)*(1-PPSNV(:))) * PRR(:) + ZRUIR(:) + ZRUIR2(:)
!
PDRIP_DP(:) = ZRUIR(:) + ZRUIR2(:)
PWR = REAL(PWR_DP)
PPG = REAL(PPG_DP)
PDRIP = REAL(PDRIP_DP)
PRRVEG = REAL(PRRVEG_DP)
IF (LHOOK) CALL DR_HOOK('HYDRO_VEG',1,ZHOOK_HANDLE)

!
!-------------------------------------------------------------------------------
!
END SUBROUTINE HYDRO_VEG
