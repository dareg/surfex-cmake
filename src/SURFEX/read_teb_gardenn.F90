!     #########
      SUBROUTINE READ_TEB_GARDEN_n(HPROGRAM,HPATCH)
!     ##################################
!
!!****  *READ_TEB_GARDEN_n* - routine to initialise ISBA variables
!!                         
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!!
!!      READ_SURF for general reading : 08/2003 (S.Malardel)
!!      B. Decharme  2008    : Floodplains
!!      B. Decharme  01/2009 : Optional Arpege deep soil temperature read
!!      B. Decharme  09/2012 : suppress NWG_LAYER (parallelization problems)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_CO2V_PAR,       ONLY : XANFMINIT, XCONDCTMIN
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
!                                
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_SNOW_PAR,       ONLY : XZ0SN
!
USE MODI_READ_SURF
!
USE MODI_INIT_IO_SURF_n
USE MODI_SET_SURFEX_FILEIN
USE MODI_END_IO_SURF_n
USE MODI_TOWN_PRESENCE
USE MODI_ALLOCATE_GR_SNOW
USE MODI_READ_GR_SNOW
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_TYPE_DIM_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! calling program
 CHARACTER(LEN=3),  INTENT(IN)  :: HPATCH   ! current TEB patch identificator
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
LOGICAL           :: GTOWN          ! town variables written in the file
INTEGER           :: IVERSION, IBUGFIX
INTEGER           :: ILU            ! 1D physical dimension
INTEGER           :: IRESP          ! Error code after redding
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=4)  :: YLVL
REAL, DIMENSION(:),ALLOCATABLE  :: ZWORK      ! 2D array to write data in file
!
INTEGER :: IWORK   ! Work integer
!
INTEGER :: JLAYER, JNBIOMASS  ! loop counter on layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_TEB_GARDEN_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_TOWN'
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     'TOWN  ',ILU)
!
YRECFM='VERSION'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,IVERSION,IRESP)
!
YRECFM='BUG'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,IBUGFIX,IRESP)
!
!*       2.     Prognostic fields:
!               -----------------
!
ALLOCATE(ZWORK(ILU))
!* soil temperatures
!
IWORK=TGDO%NGROUND_LAYER
!
ALLOCATE(TGD%CUR%XTG(ILU,IWORK))
DO JLAYER=1,IWORK
  WRITE(YLVL,'(I2)') JLAYER
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
    YRECFM=HPATCH//'GD_TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ELSE
    YRECFM='TWN_TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ENDIF
  YRECFM=ADJUSTL(YRECFM)  
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,ZWORK(:),IRESP)
  TGD%CUR%XTG(:,JLAYER)=ZWORK
END DO
!
!
!* soil liquid water content
!
ALLOCATE(TGD%CUR%XWG(ILU,IWORK))
DO JLAYER=1,TGDO%NGROUND_LAYER
  WRITE(YLVL,'(I2)') JLAYER
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
    YRECFM=HPATCH//'GD_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ELSE
    YRECFM='TWN_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ENDIF  
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,ZWORK(:),IRESP)
  TGD%CUR%XWG(:,JLAYER)=ZWORK
END DO
!
!* soil ice water content
!
ALLOCATE(TGD%CUR%XWGI(ILU,IWORK))
DO JLAYER=1,TGDO%NGROUND_LAYER
  WRITE(YLVL,'(I2)') JLAYER
! ajouter ici un test pour lire les anciens fichiers
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
    YRECFM=HPATCH//'GD_WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ELSE
    YRECFM='TWN_WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ENDIF  
  YRECFM=ADJUSTL(YRECFM)  
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,ZWORK(:),IRESP)
  TGD%CUR%XWGI(:,JLAYER)=ZWORK
END DO
!
!* water intercepted on leaves
!
ALLOCATE(TGD%CUR%XWR(ILU))
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
  YRECFM=HPATCH//'GD_WR'
ELSE
  YRECFM='TWN_WR'
ENDIF
YRECFM=ADJUSTL(YRECFM)
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TGD%CUR%XWR(:),IRESP)
!
!* Leaf Area Index (if prognostic)
!
IF (TVG%CPHOTO=='LAI' .OR. TVG%CPHOTO=='LST' .OR. TVG%CPHOTO=='NIT' .OR. TVG%CPHOTO=='NCB') THEN
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
    YRECFM=HPATCH//'GD_LAI'
  ELSE
    YRECFM='TWN_LAI'
  ENDIF        
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TGDPE%CUR%XLAI(:),IRESP)        
END IF
!
!* snow mantel
!
 CALL END_IO_SURF_n(HPROGRAM)
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ')
 CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
!
 CALL TOWN_PRESENCE(IOB, &
                    HPROGRAM,GTOWN)
!
 CALL END_IO_SURF_n(HPROGRAM)
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP')
 CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
!
IF (.NOT. GTOWN) THEN
  TGD%CUR%TSNOW%SCHEME='1-L'
  CALL ALLOCATE_GR_SNOW(TGD%CUR%TSNOW,ILU,1)
ELSE
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
    CALL READ_GR_SNOW(IOB, &
                      HPROGRAM,'GD',HPATCH,ILU,1,TGD%CUR%TSNOW  )
  ELSE
    CALL READ_GR_SNOW(IOB, &
                      HPROGRAM,'GARD',HPATCH,ILU,1,TGD%CUR%TSNOW  )
  ENDIF
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.  Semi-prognostic variables
!            -------------------------
!
!* aerodynamical resistance
!
ALLOCATE(TGD%CUR%XRESA(ILU))
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
  YRECFM=HPATCH//'GD_RES'
ELSE
  YRECFM='TWN_RESA'
ENDIF
YRECFM=ADJUSTL(YRECFM)
TGD%CUR%XRESA(:) = 100.
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TGD%CUR%XRESA(:),IRESP)
!
ALLOCATE(TGD%CUR%XLE(ILU))
TGD%CUR%XLE(:) = XUNDEF
!
!* ISBA-AGS variables
!
IF (TVG%CPHOTO/='NON') THEN
  ALLOCATE(TGD%CUR%XAN   (ILU)) 
  ALLOCATE(TGD%CUR%XANDAY(ILU)) 
  ALLOCATE(TGD%CUR%XANFM (ILU))
  ALLOCATE(TGDP%XANF  (ILU))
  TGD%CUR%XAN(:)    = 0.
  TGD%CUR%XANDAY(:) = 0.
  TGD%CUR%XANFM(:)  = XANFMINIT
  TGD%CUR%XLE(:)    = 0.
ELSE
  ALLOCATE(TGD%CUR%XAN   (0)) 
  ALLOCATE(TGD%CUR%XANDAY(0)) 
  ALLOCATE(TGD%CUR%XANFM (0))
  ALLOCATE(TGDP%XANF  (0))
ENDIF
!
IF(TVG%CPHOTO/='NON') THEN
  ALLOCATE(TGD%CUR%XBIOMASS         (ILU,TVG%NNBIOMASS))
  ALLOCATE(TGD%CUR%XRESP_BIOMASS    (ILU,TVG%NNBIOMASS))
ELSE
  ALLOCATE(TGD%CUR%XBIOMASS         (0,0))
  ALLOCATE(TGD%CUR%XRESP_BIOMASS    (0,0))
END IF
!
IF (TVG%CPHOTO=='AGS' .OR. TVG%CPHOTO=='AST') THEN
  !
  TGD%CUR%XBIOMASS(:,:) = 0.
  TGD%CUR%XRESP_BIOMASS(:,:) = 0.
ELSEIF (TVG%CPHOTO=='LAI' .OR. TVG%CPHOTO=='LST') THEN
  !
  TGD%CUR%XBIOMASS(:,1) = TGDP%XBSLAI(:) * TGDPE%CUR%XLAI(:)
  TGD%CUR%XRESP_BIOMASS(:,:) = 0.
ELSEIF (TVG%CPHOTO=='NIT' .OR. TVG%CPHOTO=='NCB') THEN
  !
  TGD%CUR%XBIOMASS(:,:) = 0.
  DO JNBIOMASS=1,TVG%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
      YRECFM=HPATCH//'GD_BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='TWN_BIOMASS'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF
    CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TGD%CUR%XBIOMASS(:,JNBIOMASS),IRESP)
  END DO

  TGD%CUR%XRESP_BIOMASS(:,:) = 0.
  DO JNBIOMASS=2,TVG%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
      YRECFM=HPATCH//'GD_RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='TWN_RESP_BIOM'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF    
    CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TGD%CUR%XRESP_BIOMASS(:,JNBIOMASS),IRESP)
  END DO
  !
ENDIF
!
DEALLOCATE(ZWORK)
IF (LHOOK) CALL DR_HOOK('READ_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_TEB_GARDEN_n
