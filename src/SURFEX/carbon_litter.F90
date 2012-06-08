!     #########
SUBROUTINE CARBON_LITTER (PTSTEP, PTURNOVER, PLITTER, PLIGNIN_STRUC,          &
                           PCONTROL_TEMP, PCONTROL_MOIST,                       &
                           PRESP_HETERO_LITTER, PSOILCARBON_INPUT)  

!   ###############################################################
!!**  CARBON_LITTER 
!!
!!    PURPOSE
!!    -------
!!    Calculates litter evolution.
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    none
!!
!!    REFERENCE
!!    ---------
!!
!!      Parton et al., Biogeochemestry, 1988
!!      Krinner et al., Global Biochemical Cycles, 2005
!!      Gibelin et al. 2008, AFM
!!      
!!    AUTHOR
!!    ------
!!
!!	A.-L. Gibelin           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/06/09
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
    
!USE MODD_CARBON_PAR
USE MODD_CO2V_PAR,       ONLY : XLC, XTAU_LITTER, XFRAC_LITTER, XFRAC_SOILCARB
USE MODD_CSTS,           ONLY : XDAY, XTT

!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE


!*       0.1 input

! time step in s
REAL, INTENT(IN)                                                 :: PTSTEP
! Turnover rates (gC/m**2/s)
REAL, DIMENSION(:,:), INTENT(IN)                                 :: PTURNOVER
! temperature control of heterotrophic respiration, above and below
REAL, DIMENSION(:,:), INTENT(IN)                                 :: PCONTROL_TEMP
! moisture control of heterotrophic respiration
REAL, DIMENSION(:,:), INTENT(IN)                                 :: PCONTROL_MOIST


!*       0.2 modified fields

! metabolic and structural litter, above and below ground (gC/m**2)
REAL, DIMENSION(:,:,:), INTENT(INOUT)                            :: PLITTER
! ratio Lignin/Carbon in structural litter, above and below ground (gC/m**2)
REAL, DIMENSION(:,:), INTENT(INOUT)                              :: PLIGNIN_STRUC

!*       0.3 output

! litter heterotrophic respiration (in gC/m**2/day)
REAL, DIMENSION(:), INTENT(OUT)                                  :: PRESP_HETERO_LITTER
! quantity of carbon going into carbon pools from litter decomposition
!   (gC/m**2/day)
REAL, DIMENSION(:,:), INTENT(OUT)                                :: PSOILCARBON_INPUT

!*       0.4 local

! time step in days
REAL                                                             :: ZDT
! fraction of structural or metabolic litter decomposed
REAL, DIMENSION(SIZE(PLITTER,1))                                 :: ZFD
! quantity of structural or metabolic litter decomposed (gC/m**2)
REAL, DIMENSION(SIZE(PLITTER,1))                                 :: ZQD
! old structural litter, above and below (gC/m**2)
REAL, DIMENSION(SIZE(PLITTER,1),SIZE(PLITTER,3))                 :: ZOLD_STRUC
! increase of metabolic and structural litter, above and below ground (gC/m**2)
REAL, DIMENSION(SIZE(PLITTER,1),SIZE(PLITTER,2),SIZE(PLITTER,3)) :: ZLITTER_INC
! lignin increase in structural litter, above and below ground (gC/m**2)
REAL, DIMENSION(SIZE(PLITTER,1),SIZE(PLITTER,3))                 :: ZLIGNIN_STRUC_INC
! dimensions
INTEGER                                                          :: INLITTER,INLITTLEVS
! indices
INTEGER                                                          :: K,L
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! correspondence between array indices and biomass compartments
! LEAF = 1
! STRUCT_ACT = 2
! STRUCT_PAS = 3
! STRUCT_BELOW = 4
! WOOD_ABOVE = 5
! WOOD_BELOW = 6
! correspondence between array indices and litter type
! LT_METABOLIC = 1
! LT_STRUCTURAL = 2
! correspondence between array indices and litter levels
! LT_ABOVE = 1
! LT_BELOW = 2
! correspondence between array indices and soil carbon pools
! SL_ACTIVE = 1
! SL_SLOW = 2
! SL_PASSIVE = 3
!-------------------------------------------------------------------------------
!
!*       1 Initialisations
!

!
!*       1.1 dimensions
!
IF (LHOOK) CALL DR_HOOK('CARBON_LITTER',0,ZHOOK_HANDLE)
INLITTER = SIZE(PLITTER,2)
INLITTLEVS = SIZE(PLITTER,3)

!
!*       1.2 set output to zero
!

PRESP_HETERO_LITTER(:) = 0.0
PSOILCARBON_INPUT(:,:) = 0.0

!
!*       2 Add biomass to different litterpools
!

ZDT = PTSTEP/XDAY

!
!*       2.1 first, save old structural litter (needed for lignin fractions).
!            (above/below)
!

DO L = 1, INLITTLEVS
    ZOLD_STRUC(:,L) = PLITTER(:,2,L)
ENDDO

!
! *      2.2 update litter, and lignin content in structural litter
!

ZLITTER_INC(:,:,:) = 0.0
ZLIGNIN_STRUC_INC(:,:) = 0.0

DO K = 1, INLITTER    ! METABOLIC AND STRUCTURAL

  !*       2.2.1 calculate litter increase (per m**2 of ground).
  !              Litter increase for structural and metabolic, above/below


  ZLITTER_INC(:,K,1) = ( XFRAC_LITTER(1,K) * PTURNOVER(:,1) + &
      XFRAC_LITTER(2,K) * PTURNOVER(:,2) + XFRAC_LITTER(3,K) * PTURNOVER(:,3) + &
      XFRAC_LITTER(5,K) * PTURNOVER(:,5) ) * PTSTEP  

  ZLITTER_INC(:,K,2) = ( XFRAC_LITTER(4,K) * PTURNOVER(:,4) + &
      XFRAC_LITTER(6,K) * PTURNOVER(:,6 ) ) * PTSTEP  

  !*       2.2.2 lignin increase in structural litter

  IF ( K .EQ. 2 ) THEN

    ZLIGNIN_STRUC_INC(:,1) = ZLIGNIN_STRUC_INC(:,1) + &
        ( XLC(1) * PTURNOVER(:,1) + XLC(2) * PTURNOVER(:,2) + &
          XLC(3) * PTURNOVER(:,3) + XLC(5) * PTURNOVER(:,5 ) ) * PTSTEP  

    ZLIGNIN_STRUC_INC(:,2) = ZLIGNIN_STRUC_INC(:,2) + &
        ( XLC(4)*PTURNOVER(:,4) + XLC(6)*PTURNOVER(:,6 ) ) * PTSTEP  

  ENDIF

ENDDO

!*       2.2.3 add new litter (struct/met, above/below)

PLITTER(:,:,:) = PLITTER(:,:,:) + ZLITTER_INC(:,:,:)

!*       2.2.4 for security: can't add more lignin than structural litter
!              (above/below)

DO L = 1, INLITTLEVS
    ZLIGNIN_STRUC_INC(:,L) = MIN( ZLIGNIN_STRUC_INC(:,L), ZLITTER_INC(:,2,L) )
ENDDO

!*       2.2.5 new lignin content: add old lignin and lignin increase, divide by 
!              total structural litter (above/below)

WHERE ( PLITTER(:,2,:) .GT. 0.0 )

  PLIGNIN_STRUC(:,:) = &
             ( PLIGNIN_STRUC(:,:)*ZOLD_STRUC(:,:) + ZLIGNIN_STRUC_INC(:,:) ) / &
             PLITTER(:,2,:)  

ENDWHERE

!
!*       3 fluxes from litter to carbon pools and respiration
!

DO L = 1, INLITTLEVS

    !
    !*       3.1 structural litter: goes into active and slow carbon pools + respiration
    !

    !*       3.1.1 total quantity of structural litter which is decomposed

    ZFD(:) = PTSTEP/XTAU_LITTER(2) * &
              PCONTROL_TEMP(:,L) * PCONTROL_MOIST(:,L) * EXP( -3. * PLIGNIN_STRUC(:,L) )  

    ZQD(:) = PLITTER(:,2,L) * ZFD(:)

    PLITTER(:,2,L) = PLITTER(:,2,L) - ZQD(:)

    !*       3.1.2 non-lignin fraction of structural litter goes into
    !              active carbon pool + respiration

    PSOILCARBON_INPUT(:,1) = PSOILCARBON_INPUT(:,1) + &
          XFRAC_SOILCARB(2,1,L) * ZQD(:) * ( 1. - PLIGNIN_STRUC(:,L) ) / ZDT  

    PRESP_HETERO_LITTER(:) = PRESP_HETERO_LITTER(:) + &
          ( 1. - XFRAC_SOILCARB(2,1,L) ) * ZQD(:) * &
          ( 1. - PLIGNIN_STRUC(:,L) ) / ZDT  

    !*       3.1.3 lignin fraction of structural litter goes into
    !              slow carbon pool + respiration

    PSOILCARBON_INPUT(:,2) = PSOILCARBON_INPUT(:,2) + &
          XFRAC_SOILCARB(2,2,L) * ZQD(:) * PLIGNIN_STRUC(:,L) / ZDT  

    PRESP_HETERO_LITTER(:) = PRESP_HETERO_LITTER(:) + &
          ( 1. - XFRAC_SOILCARB(2,2,L) ) * ZQD(:) * PLIGNIN_STRUC(:,L) / ZDT  

    !
    !*       3.2 metabolic litter goes into active carbon pool + respiration
    !

    !*       3.2.1 total quantity of metabolic litter that is decomposed

    ZFD(:) = PTSTEP/XTAU_LITTER(1) * PCONTROL_TEMP(:,L) * PCONTROL_MOIST(:,L)

    ZQD(:) = PLITTER(:,1,L) * ZFD(:)

    PLITTER(:,1,L) = PLITTER(:,1,L) - ZQD(:)

    !*       3.2.2 put decomposed litter into carbon pool + respiration

    PSOILCARBON_INPUT(:,1) = PSOILCARBON_INPUT(:,1) + &
                            XFRAC_SOILCARB(1,1,L) * ZQD(:) / ZDT  

    PRESP_HETERO_LITTER(:) = PRESP_HETERO_LITTER(:) + &
                         ( 1. - XFRAC_SOILCARB(1,1,L) ) * ZQD(:) / ZDT  

ENDDO
IF (LHOOK) CALL DR_HOOK('CARBON_LITTER',1,ZHOOK_HANDLE)

!
END SUBROUTINE CARBON_LITTER
