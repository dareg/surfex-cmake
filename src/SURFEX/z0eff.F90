!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE Z0EFF (HSNOW_SCHEME, &
                      OMEB, PALFA, PZREF, PUREF, PZ0, PZ0REL, PPSN,             &
                      PPALPHAN,PZ0LITTER, PWSNOW, PRHOSNOW, ISS, PFF, PZ0_FLOOD,&
                      PZ0_O_Z0H, PZ0_WITH_SNOW, PZ0H_WITH_SNOW,PZ0EFF,          &
                      PZ0G_WITHOUT_SNOW,                                        &
                      PZ0_MEBV,PZ0H_MEBV,PZ0EFF_MEBV,                           &
                      PZ0_MEBN,PZ0H_MEBN,PZ0EFF_MEBN                            )

!   ############################################################################
!
!!****  *Z0EFF*  
!!
!!    PURPOSE
!!    -------
!
!     Calculates the z0eff for momentum fluxes according to wind direction.
!         
!     
!!**  METHOD
!!    ------
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    Mascart et al. (1995)
!!    Belair (1995)
!!      
!!    AUTHOR
!!    ------
!!
!!      S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    13/03/95 
!!      (J.Stein)   15/11/95  use the potential temperature to compute Ri
!!                            and PVMOD instead of ZVMOD
!!      (P.Lacarrere)15/03/96 replace * PEXNS by / PEXNS
!!      (V.Masson)   22/12/97 computation of z0eff after snow treatment
!!      (V.Masson)   05/10/98 clear routine
!!      (A.Boone)    11/26/98 Option for PDELTA: forested vs default surface
!!      (V Masson)   12/07/01 new formulation for aggregation with snow z0
!!      (P.LeMoigne) 09/02/06 computation of z0h in presence of snow
!!      (B. Decharme)    2008 floodplains
!!      (P. Samuelsson) 10/2014 MEB
!!      (P. LeMoigne)   12/2014 EBA scheme update
!!      (F. Svabik)     08/2023 unapproximated roughness length averaging;
!!                              quadratic averaging snow-nowsnow for EBA scheme;
!!                              snow roughness length taken from XZ0SN;
!!                              orographic contribution via key LZ0_EFF;
!!                              snow effect on roughness length via snow height
!!                              (key LZ0SNOWH_ARP)
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_SSO_n, ONLY : SSO_t
!
USE MODD_CSTS,     ONLY : XPI, XG
USE MODD_SNOW_PAR, ONLY : XZ0SN, XZ0HSN
!
USE MODI_SUBSCALE_Z0EFF
USE MODD_SURF_ATM, ONLY : LZ0_EFF, XZ0_OFFSET, LZ0SNOWH_ARP, XRZ0_TO_HEIGHT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
 CHARACTER(LEN=*), INTENT(IN) :: HSNOW_SCHEME
!
LOGICAL, INTENT(IN)             :: OMEB           ! True = patch with multi-energy balance 
!                                                 ! False = patch with classical ISBA
REAL, DIMENSION(:), INTENT(IN)  :: PALFA          ! wind direction from J axis (clockwise)
REAL, DIMENSION(:), INTENT(IN)  :: PZREF          ! height of atmospheric level
REAL, DIMENSION(:), INTENT(IN)  :: PUREF          ! reference height for wind
REAL, DIMENSION(:), INTENT(IN)  :: PZ0            ! vegetation roughness length
REAL, DIMENSION(:), INTENT(IN)  :: PZ0REL         ! 1d orographic roughness length
REAL, DIMENSION(:), INTENT(IN)  :: PPSN           ! fraction of snow
REAL, DIMENSION(:), INTENT(IN)  :: PPALPHAN       ! snow/canopy transition coefficient
TYPE(SSO_t), INTENT(INOUT) :: ISS
REAL, DIMENSION(:), INTENT(IN)  :: PZ0_O_Z0H      ! ratio between heat and momentum z0
!
REAL, DIMENSION(:), INTENT(IN)  :: PFF            ! fraction of flood
REAL, DIMENSION(:), INTENT(IN)  :: PZ0_FLOOD      ! floodplains roughness length
!
! For multi-energy balance
REAL, DIMENSION(:), INTENT(IN)  :: PZ0LITTER      ! ground litter roughness length for MEB
!
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW         ! equivalent snow water content
REAL, DIMENSION(:), INTENT(IN)  :: PRHOSNOW       ! equivalent snow water content
!
REAL, DIMENSION(:), INTENT(OUT) :: PZ0_WITH_SNOW  ! vegetation z0 modified by snow
REAL, DIMENSION(:), INTENT(OUT) :: PZ0H_WITH_SNOW ! vegetation z0h modified by snow
REAL, DIMENSION(:), INTENT(OUT) :: PZ0EFF         ! effective z0
!
! For multi-energy balance
REAL, DIMENSION(:), INTENT(OUT) :: PZ0G_WITHOUT_SNOW  ! roughness length for momentum at snow-free canopy floor
!
REAL, DIMENSION(:), INTENT(OUT) :: PZ0_MEBV           ! roughness length for momentum over MEB vegetation part of patch
REAL, DIMENSION(:), INTENT(OUT) :: PZ0H_MEBV          ! roughness length for heat over MEB vegetation part of path
REAL, DIMENSION(:), INTENT(OUT) :: PZ0EFF_MEBV        ! roughness length for momentum over MEB vegetation part of patch
!                                                     ! eventually including orograhic roughness
REAL, DIMENSION(:), INTENT(OUT) :: PZ0_MEBN           ! roughness length for momentum over MEB snow part of patch
REAL, DIMENSION(:), INTENT(OUT) :: PZ0H_MEBN          ! roughness length for heat over MEB snow part of path
REAL, DIMENSION(:), INTENT(OUT) :: PZ0EFF_MEBN        ! roughness length for momentum over MEB snow part of patch
!
!
!
!*      0.2    declarations of local variables
!
!
!
REAL, DIMENSION(SIZE(PZ0EFF)) :: ZWORK, ZWORKH, ZALFA, ZPFF, ZZ0_FLOOD, ZZ0H_FLOOD
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('Z0EFF',0,ZHOOK_HANDLE)
ZALFA(:) = PALFA(:)
WHERE(ZALFA(:)<=-XPI) ZALFA = ZALFA + 2.*XPI
WHERE(ZALFA(:)>  XPI) ZALFA = ZALFA - 2.*XPI
!
! Initialisation of MEB roughness lengths
PZ0G_WITHOUT_SNOW=0.
PZ0_MEBV=0.
PZ0H_MEBV=0.
PZ0EFF_MEBV=0.
PZ0_MEBN=0.
PZ0H_MEBN=0.
PZ0EFF_MEBN=0.
!
!*       1.     GRID-AVERAGED ROUGHNESS LENGTHS
!               -------------------------------
!       (considering the effect of snow-flood-covered surfaces and orography)
!
!*       1.1    for heat
!               --------
!
PZ0_WITH_SNOW(:)  = PZ0(:)
PZ0H_WITH_SNOW(:) = PZ0(:) / PZ0_O_Z0H(:)
!
IF(HSNOW_SCHEME=='EBA') THEN
       
  ! Flooding scheme not implemented with this option
  IF (LZ0SNOWH_ARP) THEN
    WHERE (PPSN(:) > 0..AND.PPSN(:) <= 1.)
      PZ0_WITH_SNOW (:)=MAX(PZ0(:)-XRZ0_TO_HEIGHT*PWSNOW(:)/PRHOSNOW(:),XZ0SN)
      PZ0H_WITH_SNOW(:)=PZ0_WITH_SNOW(:)/PZ0_O_Z0H(:)
    ENDWHERE 
  ELSE
    WHERE (PPSN(:) > 0..AND.PPSN(:) <= 1.)
      PZ0_WITH_SNOW (:)=SQRT((1.-PPSN(:))*PZ0(:)**2+PPSN(:)*XZ0SN**2)
      PZ0H_WITH_SNOW(:)=PZ0_WITH_SNOW(:)/PZ0_O_Z0H(:)
    ENDWHERE 
  ENDIF

ELSE

  ! PZ0_FLOOD secured from undef values
  ZZ0_FLOOD (:)=MIN(PZ0_FLOOD(:),1000.)
  ZZ0H_FLOOD(:)=ZZ0_FLOOD(:)/PZ0_O_Z0H(:)

  IF (LZ0SNOWH_ARP) THEN

    WHERE (PPSN(:) > 0..OR.PFF(:) > 0.)

      ! roughness of non-flooded part (including snow)
      PZ0_WITH_SNOW (:)=MAX(PZ0(:)-XRZ0_TO_HEIGHT*PWSNOW(:)/PRHOSNOW(:),XZ0SN)
      PZ0H_WITH_SNOW(:)=PZ0_WITH_SNOW(:)/PZ0_O_Z0H(:)

      ! averaging with flooded part
      ZWORK(:) =  PFF(:) /(LOG(XZ0_OFFSET+PUREF(:)/ZZ0_FLOOD    (:)))**2 &
        +(1.-PFF(:))/(LOG(XZ0_OFFSET+PUREF(:)/PZ0_WITH_SNOW(:)))**2
      ZWORKH(:)=  PFF(:) /(LOG(XZ0_OFFSET+PZREF(:)/ZZ0H_FLOOD    (:))*LOG(XZ0_OFFSET*(1. &
        +PUREF(:)/ZZ0_FLOOD    (:)-PZREF(:)/ZZ0H_FLOOD    (:))+PZREF(:)/ZZ0H_FLOOD    (:)))   &
        +(1.-PFF(:))/(LOG(XZ0_OFFSET+PZREF(:)/PZ0H_WITH_SNOW(:))*LOG(XZ0_OFFSET*(1. &
        +PUREF(:)/PZ0_WITH_SNOW(:)-PZREF(:)/PZ0H_WITH_SNOW(:))+PZREF(:)/PZ0H_WITH_SNOW(:)))

      PZ0_WITH_SNOW (:)=PUREF(:)/(EXP(SQRT(1./ZWORK(:)))-XZ0_OFFSET)
      PZ0H_WITH_SNOW(:)=PZREF(:)/(EXP((XZ0_OFFSET*(SQRT(ZWORK(:)/ZWORKH(:))-1.)+1.)/SQRT(ZWORKH(:)))-XZ0_OFFSET)
 
    ENDWHERE

  ELSE

    WHERE (PPSN(:) > 0..OR.PFF(:) > 0.)
        
      ZWORK(:) =          PPSN(:) /(LOG(XZ0_OFFSET + PUREF(:)/XZ0SN       ))**2 &
        +                 PFF (:) /(LOG(XZ0_OFFSET + PUREF(:)/ZZ0_FLOOD(:)))**2 &
        +(1.-PPSN(:)-PFF (:))/(LOG(XZ0_OFFSET + PUREF(:)/PZ0      (:)))**2
      ZWORKH(:)=         PPSN(:) /(LOG(XZ0_OFFSET + PZREF(:)/XZ0HSN             )*LOG(XZ0_OFFSET*(1. &
        +PUREF(:)/XZ0SN       -PZREF(:)/XZ0HSN             )+PZREF(:)/XZ0HSN             ))               &
        +                 PFF(:) /(LOG(XZ0_OFFSET + PZREF(:)/ZZ0H_FLOOD(:)      )*LOG(XZ0_OFFSET*(1. &
        +PUREF(:)/ZZ0_FLOOD(:)-PZREF(:)/ZZ0H_FLOOD(:)      )+PZREF(:)/ZZ0H_FLOOD(:)      ))               &
        +(1.-PPSN(:)-PFF(:))/(LOG(XZ0_OFFSET + PZREF(:)*PZ0_O_Z0H(:)/PZ0(:))*LOG(XZ0_OFFSET*(1. &
        +PUREF(:)/PZ0(:)      -PZREF(:)*PZ0_O_Z0H(:)/PZ0(:))+PZREF(:)*PZ0_O_Z0H(:)/PZ0(:)))

      PZ0_WITH_SNOW (:)=PUREF(:)/(EXP(SQRT(1./ZWORK(:)))-XZ0_OFFSET)
      PZ0H_WITH_SNOW(:)=PZREF(:)/(EXP((XZ0_OFFSET*(SQRT(ZWORK(:)/ZWORKH(:))-1.)+1.)/SQRT(ZWORKH(:)))-XZ0_OFFSET)

    ENDWHERE

  ENDIF
        
ENDIF

! For multi-energy balance (LSNOWH_ARP not implemented for MEB)
IF (OMEB) THEN

  ! roughness length for momentum at snow-free canopy floor
  PZ0G_WITHOUT_SNOW(:)=PZ0LITTER
  WHERE (PFF(:) > 0.)
    ZPFF(:)=PFF(:)/(1.-PPSN(:)+1.E-6)
    ZWORK(:)=(    ZPFF(:) *LOG(PZ0_FLOOD(:)) ) &
      +( (1.-ZPFF(:))*LOG(PZ0LITTER(:)) )
    PZ0G_WITHOUT_SNOW(:)=EXP(ZWORK(:))
  ENDWHERE

  ! roughness lengths for momentum and heat over MEB vegetation part of patch
  PZ0_MEBV (:)=PZ0(:)
  PZ0H_MEBV(:)=PZ0(:)/PZ0_O_Z0H(:)

  ! roughness lengths for momentum and heat over MEB snow part of patch
  ZWORK(:) =  PPALPHAN(:) /(LOG(XZ0_OFFSET+PUREF(:)/XZ0SN      ))**2 &
    +(1.-PPALPHAN(:))/(LOG(XZ0_OFFSET+PUREF(:)/PZ0_MEBV(:)))**2
  ZWORKH(:)=  PPALPHAN(:) /(LOG(XZ0_OFFSET+PZREF(:)/XZ0HSN      )*                                 &
    LOG(XZ0_OFFSET*(1. + PUREF(:)/XZ0SN       -PZREF(:)/XZ0HSN      )+PZREF(:)/XZ0HSN      )) &
    +(1.-PPALPHAN(:))/(LOG(XZ0_OFFSET+PZREF(:)/PZ0H_MEBV(:))*                                 &
    LOG(XZ0_OFFSET*(1. + PUREF(:)/PZ0_MEBV(:) -PZREF(:)/PZ0H_MEBV(:))+PZREF(:)/PZ0H_MEBV(:)))

  PZ0_MEBN (:)=PUREF(:)/(EXP(SQRT(1./ZWORK(:)))-XZ0_OFFSET)
  PZ0H_MEBN(:)=PZREF(:)/(EXP((XZ0_OFFSET*(SQRT(ZWORK(:)/ZWORKH(:))-1.)+1.)/SQRT(ZWORKH(:)))-XZ0_OFFSET)

  ! roughness lengths for momentum and heat over MEB total patch
  ZWORK(:) =  PPSN(:) /(LOG(XZ0_OFFSET+PUREF(:)/PZ0_MEBN(:)))**2 &
    +(1.-PPSN(:))/(LOG(XZ0_OFFSET+PUREF(:)/PZ0_MEBV(:)))**2
  ZWORKH(:)=  PPSN(:) /(LOG(XZ0_OFFSET+PZREF(:)/PZ0H_MEBN(:))*                                  &
    LOG(XZ0_OFFSET*(1.+PUREF(:)/PZ0_MEBN(:)-PZREF(:)/PZ0H_MEBN(:))+PZREF(:)/PZ0H_MEBN(:))) &
    +(1.-PPSN(:))/(LOG(XZ0_OFFSET+PZREF(:)/PZ0H_MEBV(:))*                                  &
    LOG(XZ0_OFFSET*(1.+PUREF(:)/PZ0_MEBV(:)-PZREF(:)/PZ0H_MEBV(:))+PZREF(:)/PZ0H_MEBV(:)))

  PZ0_WITH_SNOW (:) = PUREF(:)/(EXP(SQRT(1./ZWORK(:))) - XZ0_OFFSET)
  PZ0H_WITH_SNOW(:) = PZREF(:)/(EXP((XZ0_OFFSET*(SQRT(ZWORK(:)/ZWORKH(:)) - 1.) + 1.)/SQRT(ZWORKH(:))) - XZ0_OFFSET)

ENDIF
!
!*       1.2    for momentum
!               ------------
!
!                                     In this particular case, we now use
!                                     the roughness length due to the coupled
!                                     effect of vegetation and topography
!                                     Snow and Flood effects are yet taken
!                                     into account through ZZ0EFF
!
IF (LZ0_EFF) THEN
  ! add orographic component
  PZ0EFF       (:)=SQRT(PZ0_WITH_SNOW(:)**2+PZ0REL(:)**2)
  PZ0_WITH_SNOW(:)=PZ0EFF(:)  ! use z0eff consistently everywhere
  IF (OMEB) THEN
    PZ0EFF_MEBV(:)=SQRT(PZ0_MEBV(:)**2+PZ0REL(:)**2)
    PZ0EFF_MEBN(:)=SQRT(PZ0_MEBN(:)**2+PZ0REL(:)**2)
    PZ0_MEBV   (:)=PZ0EFF_MEBV(:)
    PZ0_MEBN   (:)=PZ0EFF_MEBN(:)
  ENDIF
ELSE
  ! do not add orographic component
  PZ0EFF(:) = PZ0_WITH_SNOW(:)
  IF(OMEB)THEN
    PZ0EFF_MEBV(:) = PZ0_MEBV(:)
    PZ0EFF_MEBN(:) = PZ0_MEBN(:)
  ENDIF
ENDIF
!
IF (LHOOK) CALL DR_HOOK('Z0EFF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE Z0EFF
