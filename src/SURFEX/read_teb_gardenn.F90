!     #########
      SUBROUTINE READ_TEB_GARDEN_n (DTCO, DGU, U, GDM, &
                                    HPROGRAM,HPATCH)
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
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURFEX_n, ONLY : TEB_GARDEN_MODEL_t
!
USE MODD_CO2V_PAR,       ONLY : XANFMINIT, XCONDCTMIN
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
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_t), INTENT(INOUT) :: DGU
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(TEB_GARDEN_MODEL_t), INTENT(INOUT) :: GDM
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
 CALL READ_SURF(&
                HPROGRAM,YRECFM,IVERSION,IRESP)
!
YRECFM='BUG'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,IBUGFIX,IRESP)
!
!*       2.     Prognostic fields:
!               -----------------
!
ALLOCATE(ZWORK(ILU))
!* soil temperatures
!
IWORK=GDM%TV%O%NGROUND_LAYER
!
ALLOCATE(GDM%TV%R%CUR%XTG(ILU,IWORK,1))
DO JLAYER=1,IWORK
  WRITE(YLVL,'(I2)') JLAYER
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
    YRECFM=HPATCH//'GD_TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ELSE
    YRECFM='TWN_TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ENDIF
  YRECFM=ADJUSTL(YRECFM)  
  CALL READ_SURF(&
                HPROGRAM,YRECFM,ZWORK(:),IRESP)
  GDM%TV%R%CUR%XTG(:,JLAYER,1)=ZWORK
END DO
!
!
!* soil liquid water content
!
ALLOCATE(GDM%TV%R%CUR%XWG(ILU,IWORK,1))
DO JLAYER=1,GDM%TV%O%NGROUND_LAYER
  WRITE(YLVL,'(I2)') JLAYER
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
    YRECFM=HPATCH//'GD_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ELSE
    YRECFM='TWN_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ENDIF  
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(&
                HPROGRAM,YRECFM,ZWORK(:),IRESP)
  GDM%TV%R%CUR%XWG(:,JLAYER,1)=ZWORK
END DO
!
!* soil ice water content
!
ALLOCATE(GDM%TV%R%CUR%XWGI(ILU,IWORK,1))
DO JLAYER=1,GDM%TV%O%NGROUND_LAYER
  WRITE(YLVL,'(I2)') JLAYER
! ajouter ici un test pour lire les anciens fichiers
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
    YRECFM=HPATCH//'GD_WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ELSE
    YRECFM='TWN_WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ENDIF  
  YRECFM=ADJUSTL(YRECFM)  
  CALL READ_SURF(&
                HPROGRAM,YRECFM,ZWORK(:),IRESP)
  GDM%TV%R%CUR%XWGI(:,JLAYER,1)=ZWORK
END DO
!
!* water intercepted on leaves
!
ALLOCATE(GDM%TV%R%CUR%XWR(ILU,1))
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
  YRECFM=HPATCH//'GD_WR'
ELSE
  YRECFM='TWN_WR'
ENDIF
YRECFM=ADJUSTL(YRECFM)
 CALL READ_SURF(&
                HPROGRAM,YRECFM,GDM%TV%R%CUR%XWR(:,1),IRESP)
!
!* Leaf Area Index (if prognostic)
!
IF (GDM%TV%O%CPHOTO=='LAI' .OR. GDM%TV%O%CPHOTO=='LST' .OR. &
                GDM%TV%O%CPHOTO=='NIT' .OR. GDM%TV%O%CPHOTO=='NCB') THEN
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
    YRECFM=HPATCH//'GD_LAI'
  ELSE
    YRECFM='TWN_LAI'
  ENDIF        
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(&
                HPROGRAM,YRECFM,GDM%TV%M%T%CUR%XLAI(:,1),IRESP)        
END IF
!
!* snow mantel
!
 CALL END_IO_SURF_n(HPROGRAM)
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ')
CALL INIT_IO_SURF_n(DTCO, DGU, U, &
                     HPROGRAM,'FULL  ','SURF  ','READ ')
!
 CALL TOWN_PRESENCE(&
                    HPROGRAM,GTOWN)
!
 CALL END_IO_SURF_n(HPROGRAM)
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP')
CALL INIT_IO_SURF_n(DTCO, DGU, U, &
                     HPROGRAM,'TOWN  ','TEB   ','READ ')
!
IF (.NOT. GTOWN) THEN
  GDM%TV%R%CUR%TSNOW%SCHEME='1-L'
  CALL ALLOCATE_GR_SNOW(GDM%TV%R%CUR%TSNOW,ILU,1)
ELSE
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
    CALL READ_GR_SNOW(&
                      HPROGRAM,'GD',HPATCH,ILU,1,GDM%TV%R%CUR%TSNOW  )
  ELSE
    CALL READ_GR_SNOW(&
                      HPROGRAM,'GARD',HPATCH,ILU,1,GDM%TV%R%CUR%TSNOW  )
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
ALLOCATE(GDM%TV%R%CUR%XRESA(ILU,1))
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
  YRECFM=HPATCH//'GD_RES'
ELSE
  YRECFM='TWN_RESA'
ENDIF
YRECFM=ADJUSTL(YRECFM)
GDM%TV%R%CUR%XRESA(:,:) = 100.
 CALL READ_SURF(&
                HPROGRAM,YRECFM,GDM%TV%R%CUR%XRESA(:,1),IRESP)
!
ALLOCATE(GDM%TV%R%CUR%XLE(ILU,1))
GDM%TV%R%CUR%XLE(:,1) = XUNDEF
!
!* ISBA-AGS variables
!
IF (GDM%TV%O%CPHOTO/='NON') THEN
  ALLOCATE(GDM%TV%R%CUR%XAN   (ILU,1)) 
  ALLOCATE(GDM%TV%R%CUR%XANDAY(ILU,1)) 
  ALLOCATE(GDM%TV%R%CUR%XANFM (ILU,1))
  GDM%TV%R%CUR%XAN(:,1)    = 0.
  GDM%TV%R%CUR%XANDAY(:,1) = 0.
  GDM%TV%R%CUR%XANFM(:,1)  = XANFMINIT
  GDM%TV%R%CUR%XLE(:,1)    = 0.
ELSE
  ALLOCATE(GDM%TV%R%CUR%XAN   (0,0)) 
  ALLOCATE(GDM%TV%R%CUR%XANDAY(0,0)) 
  ALLOCATE(GDM%TV%R%CUR%XANFM (0,0))
ENDIF
!
IF(GDM%TV%O%CPHOTO/='NON') THEN
  ALLOCATE(GDM%TV%R%CUR%XBIOMASS         (ILU,GDM%TV%O%NNBIOMASS,1))
  ALLOCATE(GDM%TV%R%CUR%XRESP_BIOMASS    (ILU,GDM%TV%O%NNBIOMASS,1))
ELSE
  ALLOCATE(GDM%TV%R%CUR%XBIOMASS         (0,0,0))
  ALLOCATE(GDM%TV%R%CUR%XRESP_BIOMASS    (0,0,0))
END IF
!
IF (GDM%TV%O%CPHOTO=='AGS' .OR. GDM%TV%O%CPHOTO=='AST') THEN
  !
  GDM%TV%R%CUR%XBIOMASS(:,:,:) = 0.
  GDM%TV%R%CUR%XRESP_BIOMASS(:,:,:) = 0.
ELSEIF (GDM%TV%O%CPHOTO=='LAI' .OR. GDM%TV%O%CPHOTO=='LST') THEN
  !
  GDM%TV%R%CUR%XBIOMASS(:,1,:) = GDM%TV%M%T%CUR%XBSLAI(:,:) * GDM%TV%M%T%CUR%XLAI(:,:)
  GDM%TV%R%CUR%XRESP_BIOMASS(:,:,:) = 0.
ELSEIF (GDM%TV%O%CPHOTO=='NIT' .OR. GDM%TV%O%CPHOTO=='NCB') THEN
  !
  GDM%TV%R%CUR%XBIOMASS(:,:,:) = 0.
  DO JNBIOMASS=1,GDM%TV%O%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
      YRECFM=HPATCH//'GD_BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='TWN_BIOMASS'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF
    YRECFM=ADJUSTL(YRECFM)
    CALL READ_SURF(&
                HPROGRAM,YRECFM,GDM%TV%R%CUR%XBIOMASS(:,JNBIOMASS,1),IRESP)
  END DO

  GDM%TV%R%CUR%XRESP_BIOMASS(:,:,:) = 0.
  DO JNBIOMASS=2,GDM%TV%O%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
      YRECFM=HPATCH//'GD_RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ELSE
      YRECFM='TWN_RESP_BIOM'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    ENDIF    
    YRECFM=ADJUSTL(YRECFM)
    CALL READ_SURF(&
                HPROGRAM,YRECFM,GDM%TV%R%CUR%XRESP_BIOMASS(:,JNBIOMASS,1),IRESP)
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
