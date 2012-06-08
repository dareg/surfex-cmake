!     ######spl
PROGRAM MAIN_WOOD_SPINUP
!
!
!!****  *MAIN_WOOD_SPINUP*  
!!
!!    PURPOSE
!!    -------
!!    Spinup of woody biomass.
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
!!      Gibelin et al. 2008, AFM
!!        Modelling energy and CO2 fluxes with an interactive vegetation land surface model -
!!        Evaluation at high and middle latitudes.
!!
!!    AUTHOR
!!    ------
!!	A.L. Gibelin           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/05/09
!!      B. Decharme : read pgd+prep
!!
!-------------------------------------------------------------------------------
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE, NVT_NO, NVT_ROCK, NVT_SNOW, &
  NVT_TREE, NVT_CONI, NVT_EVER, NVT_C3, NVT_C4, NVT_IRR, NVT_GRAS, NVT_TROG, NVT_PARK  
USE MODD_IO_SURF_ASC,ONLY : CFILEIN, CFILEIN_SAVE, CFILEOUT, CFILEPGD
USE MODD_IO_SURF_FA, ONLY : CFILEIN_FA, CFILEIN_FA_SAVE, CFILEOUT_FA, &
                            LFANOCOMPACT, CFILEPGD_FA
USE MODD_IO_SURF_LFI,ONLY : CLUOUT_LFI, CFILEIN_LFI, CFILEIN_LFI_SAVE,&
                            CFILEOUT_LFI, CFILEPGD_LFI
USE MODD_IO_SURF_OL,       ONLY : XSTART, XCOUNT, XSTRIDE,          &
                              LDEFINED_NATURE, LDEFINED_SEA,    &
                              LDEFINED_WATER,  LDEFINED_TOWN,   &
                              LDEFINED_SURF_ATM, LPARTW,        &
                              XSTARTW, XCOUNTW, NSTEP_OUTPUT   
USE MODD_SURF_ATM_n,     ONLY : NDIM_FULL
USE MODD_SURF_PAR,     ONLY : NUNDEF


USE MODE_POS_SURF

USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_OL_READ_ATM_CONF
USE MODI_INIT_IO_SURF_n
USE MODI_TEST_NAM_VAR_SURF
USE MODI_READ_SURF
USE MODI_SPINUP_WOOD_BIOMASS
USE MODI_WRITE_SURF
!
USE MODI_SET_SURFEX_FILEIN
!
USE MODN_IO_OFFLINE
!
! ------------------------------------------------------------------------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_GOTO_SURFEX
!
USE MODI_IO_BUFF_CLEAN_n
USE MODI_ALLOC_SURFEX
USE MODI_DEALLOC_SURFEX
USE MODI_END_IO_SURF_n
USE MODI_GET_LUOUT
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.    declarations of local variables
!
INTEGER                              :: NNBIOMASS    ! number of biomass pools
INTEGER                              :: NPATCH         ! maximum number of sub-tiles (patches)
!                                                      ! used at any grid point within a 
!
INTEGER                           :: ILUOUT              ! ascii output unit number
INTEGER                           :: ILUNAM              ! namelist unit number
CHARACTER(LEN=28)                 :: YLUOUT    = 'LISTING_WOOD                '
LOGICAL                           :: GFOUND              ! return logical when reading namelist
!
REAL                              :: ZDURATION           ! duration of run                     (s)
REAL                              :: ZTSTEP              ! atmospheric time-step            (s)
INTEGER                           :: INI                 ! grid dimension
INTEGER                           :: IYEAR               ! current year (UTC)
INTEGER                           :: IMONTH              ! current month (UTC)
INTEGER                           :: IDAY                ! current day (UTC)
REAL                              :: ZTIME               ! current time since start of the run (s)
REAL, DIMENSION(:), POINTER       :: ZLAT                ! latitude                         (rad)
REAL, DIMENSION(:), POINTER       :: ZLON                ! longitude                        (rad)
REAL, DIMENSION(:), POINTER       :: ZZS_FORC            ! orography                        (m)  
REAL, DIMENSION(:), POINTER       :: ZZREF               ! Forcing level for T
REAL, DIMENSION(:), POINTER       :: ZUREF               ! Forcing level for U
!
!REAL, DIMENSION(:,:), ALLOCATABLE :: ZTA                 ! air temperature forcing               (K)!
!
!
!
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4)  :: YLVL
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=14) :: YFORM          ! Writing format
INTEGER           :: IRESP          ! error return code
INTEGER           :: JSTEP, JNBIOMASS, JPATCH
!
CHARACTER(LEN=3)                      :: CPHOTO   ! type of photosynthesis
!
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZTAU_WOOD
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZVEGTYPE_PATCH
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZINCREASE_CURR
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZBIOMASS
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZBIOMASS_LAST
REAL(KIND=JPRB) :: ZHOOK_HANDLE


! --------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MAIN_WOOD_SPINUP',0,ZHOOK_HANDLE)
CALL ALLOC_SURFEX(1)
CALL GOTO_SURFEX(1,.TRUE.)

!
!*      1.    Initialisations
!             ---------------
!
! Vegetation type

NVEGTYPE = 12
!
NVT_NO = 1 ! no vegetation (smooth)
NVT_ROCK = 2 ! no vegetation (rocks)
NVT_SNOW = 3 ! permanent snow and ice
NVT_TREE = 4 ! forest and trees
NVT_CONI = 5 ! forest and trees (coniferous)
NVT_EVER = 6 ! forest and trees (broadleaf evergreen)
NVT_C3 = 7 ! C3 cultures types
NVT_C4 = 8 ! C4 cultures types
NVT_IRR = 9 ! irrigated crops
NVT_GRAS = 10 ! grassland
NVT_TROG = 11 ! tropical grassland
NVT_PARK = 12 ! peat bogs, parks and gardens (irrigated grass)

!
!*      2.    Reading of inputs
!             -----------------
!
CLUOUT_LFI =  ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
CALL GET_LUOUT('ASCII ',ILUOUT)
OPEN(UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YLUOUT)//'.txt'),FORM='FORMATTED',ACTION='WRITE')
!
CFILEPGD     = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')
CFILEPGD_LFI = CPGDFILE
CFILEPGD_FA  = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
!
CFILEIN     = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
CFILEIN_LFI = CPREPFILE
CFILEIN_FA  = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
!
CFILEIN_SAVE     = CFILEIN
CFILEIN_LFI_SAVE = CFILEIN_LFI
CFILEIN_FA_SAVE  = CFILEIN_FA
!
!*      2.1   Namelist
!
CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
!
CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND,ILUOUT)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
!
CALL TEST_NAM_VAR_SURF(ILUOUT,'CSURF_FILETYPE',CSURF_FILETYPE,'ASCII ','LFI   ','FA    ')
!IF (CSURF_FILETYPE/='ASCII ') THEN
!  WRITE(ILUOUT,*) '*****************************************'
!  WRITE(ILUOUT,*) '* Wood spinup is only implemented for CSURF_FILETYPE=''ASCII'' *'
!  WRITE(ILUOUT,*) '* and CSURF_FILETYPE = ',CSURF_FILETYPE,' in this run          *'
!  WRITE(ILUOUT,*) '*****************************************'
!  CALL ABOR1_SFX('MAIN_WOOD_SPINUP: WOOD SPINUP REQUIRES CSURF_FILETYPE=''ASCII''')
!ENDIF

CALL TEST_NAM_VAR_SURF(ILUOUT,'CTIMESERIES_FILETYPE',CTIMESERIES_FILETYPE, &
       'NETCDF','TEXTE ','BINARY','ASCII ','LFI   ','FA    ','NONE  ','OFFLIN')  
IF (CTIMESERIES_FILETYPE=='NETCDF') CTIMESERIES_FILETYPE='OFFLIN'

CALL CLOSE_NAMELIST('ASCII ',ILUNAM)

!
!*      2.2   Netcdf file handling
!
XSTART            = NUNDEF
XSTRIDE           = NUNDEF
XCOUNT            = NUNDEF
XSTARTW           = 0
XCOUNTW           = 1
LPARTW            = .TRUE.
LDEFINED_SURF_ATM = .TRUE.
LDEFINED_NATURE   = .TRUE.
LDEFINED_TOWN     = .TRUE.
LDEFINED_WATER    = .TRUE.
LDEFINED_SEA      = .TRUE.

!
!*      2.3   Dimensions
!

CALL OL_READ_ATM_CONF(CSURF_FILETYPE, CFORCING_FILETYPE,            &
                        ZDURATION, ZTSTEP, INI, IYEAR, IMONTH, IDAY,  &
                        ZTIME, ZLAT, ZLON, ZZS_FORC, ZZREF, ZUREF     )  

NSTEP_OUTPUT  = INT(ZDURATION / XTSTEP_OUTPUT)
!
!-----------------------------------------------------------------------------------------------------
! READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
!         Initialisation for IO
!
CALL SET_SURFEX_FILEIN(CSURF_FILETYPE,'PGD ') ! change input file name to pgd name
! Read NFULL before reading 2D arrays
CALL INIT_IO_SURF_n(CSURF_FILETYPE,'FULL  ','SURF  ','READ ')

CALL INIT_IO_SURF_n(CSURF_FILETYPE,'NATURE','SURF  ','READ ')

CALL READ_SURF(CSURF_FILETYPE,'PATCH_NUMBER',NPATCH,IRESP)
CALL READ_SURF(CSURF_FILETYPE,'NBIOMASS',NNBIOMASS,IRESP)

CALL READ_SURF(CSURF_FILETYPE,'PHOTO',CPHOTO,IRESP)
IF (CPHOTO/='NCB') THEN
  WRITE(ILUOUT,*) '**********************************************'
  WRITE(ILUOUT,*) '* Running wood spinup does not have interest *'
  WRITE(ILUOUT,*) '* with option CPHOTO = ',CPHOTO,'            *'
  WRITE(ILUOUT,*) '**********************************************'
  CALL ABOR1_SFX('MAIN_WOOD_SPINUP: DO NOT RUN WOOD SPINUP WITH CPHOTO='//CPHOTO)
ENDIF

ALLOCATE(ZVEGTYPE_PATCH(INI,NPATCH,NVEGTYPE))

! Fraction of each vegetation type for each patch

DO JPATCH=1,NPATCH

  WRITE(YLVL,'(I2)') JPATCH 
  YRECFM='VEGTYPE_P'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(CSURF_FILETYPE,YRECFM,ZVEGTYPE_PATCH(:,JPATCH,:),IRESP)

END DO

CALL END_IO_SURF_n(CSURF_FILETYPE)
CALL SET_SURFEX_FILEIN(CSURF_FILETYPE,'PREP') ! restore input file name
!
!-----------------------------------------------------------------------------------------------------
! END READ PGD FILE
!-----------------------------------------------------------------------------------------------------


      
!
!*      2.4   Model variables
!

ALLOCATE(ZTAU_WOOD(INI,NPATCH))
ALLOCATE(ZINCREASE_CURR(INI,NNBIOMASS,NPATCH,NSTEP_OUTPUT))
ALLOCATE(ZBIOMASS(INI,NNBIOMASS,NPATCH,NSTEP_OUTPUT))
ALLOCATE(ZBIOMASS_LAST(INI,NNBIOMASS,NPATCH))


! Biomass increments

DO JNBIOMASS=1,NNBIOMASS
    
  WRITE(YLVL,'(I2)') JNBIOMASS
  YRECFM='INCREASE'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))

  CALL READ_SURF(CTIMESERIES_FILETYPE,YRECFM,ZINCREASE_CURR(:,JNBIOMASS,:,:),IRESP)

END DO
!
YRECFM='TAU_WOOD'
CALL READ_SURF(CTIMESERIES_FILETYPE,YRECFM,ZTAU_WOOD(:,:),IRESP)
!
! Read NFULL before reading 2D arrays
CALL INIT_IO_SURF_n(CSURF_FILETYPE,'FULL  ','SURF  ','READ ')

CALL INIT_IO_SURF_n(CSURF_FILETYPE,'NATURE','SURF  ','READ ')

! Last values of biomass

DO JNBIOMASS=1,NNBIOMASS    
  WRITE(YLVL,'(I2)') JNBIOMASS
  YRECFM='BIOMASS'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(CSURF_FILETYPE,YRECFM,ZBIOMASS_LAST(:,JNBIOMASS,:),IRESP)
END DO
!

CALL END_IO_SURF_n(CSURF_FILETYPE)

! Initialise variables

DO JSTEP=1,NSTEP_OUTPUT
  ZBIOMASS(:,1:4,:,JSTEP) = ZBIOMASS_LAST(:,1:4,:)
ENDDO
        
ZBIOMASS(:,5:6,:,1) = ZBIOMASS_LAST(:,5:6,:)

!
!
!*      3.    Calculations
!             ------------
!

DO JSTEP=1,NSTEP_OUTPUT

  ! Evolution of woody biomass

  DO JPATCH=1,NPATCH
  
    CALL SPINUP_WOOD_BIOMASS( ZVEGTYPE_PATCH(:,JPATCH,:), &
                                ZBIOMASS(:,:,JPATCH,JSTEP), &
                                ZINCREASE_CURR(:,:,JPATCH,JSTEP), &
                                ZTAU_WOOD(:,JPATCH) )  

  ENDDO

  ! Prepare next time step

  IF ( JSTEP < NSTEP_OUTPUT ) THEN

      ZBIOMASS(:,5:6,:,JSTEP+1) = ZBIOMASS(:,5:6,:,JSTEP)
      
  ENDIF

ENDDO     
!
!*      4.    Writing of outputs
!             ------------------
! 
CFILEOUT    = ADJUSTL(ADJUSTR(CSURFFILE)//'.txt')
CFILEOUT_LFI= CFILEIN_LFI_SAVE
CFILEOUT_FA = CFILEIN_FA_SAVE
!
IF(CSURF_FILETYPE=='FA    ')LFANOCOMPACT=.TRUE.
!
CALL IO_BUFF_CLEAN_n
!
! Read NFULL before writing 2D arrays
CALL INIT_IO_SURF_n(CSURF_FILETYPE,'FULL  ','ISBA  ','READ ')
!
CALL INIT_IO_SURF_n(CSURF_FILETYPE,'NATURE','ISBA  ','WRITE')
!
DO JNBIOMASS=1,NNBIOMASS
  WRITE(YLVL,'(I4)') JNBIOMASS
  YRECFM='BIOMASS'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A11,I1.1,A8)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_BIOMASS',JNBIOMASS,' (kg/m2)'
  CALL WRITE_SURF(CSURF_FILETYPE,YRECFM,ZBIOMASS(:,JNBIOMASS,:,NSTEP_OUTPUT),IRESP,HCOMMENT=YCOMMENT)
END DO
!  
CALL END_IO_SURF_n(CSURF_FILETYPE)
!
DEALLOCATE(ZTAU_WOOD)
DEALLOCATE(ZVEGTYPE_PATCH)
DEALLOCATE(ZINCREASE_CURR) 
DEALLOCATE(ZBIOMASS) 
DEALLOCATE(ZBIOMASS_LAST) 

!
!
!*      5.    Close parallelized I/O
!             ----------------------
!
CLOSE(ILUOUT)
!
WRITE(*,*) ' '
WRITE(*,*) '    -----------------------------------'
WRITE(*,*) '    | MAIN_WOOD_SPINUP ENDS CORRECTLY |'
WRITE(*,*) '    -----------------------------------'
WRITE(*,*) ' '
CALL DEALLOC_SURFEX
IF (LHOOK) CALL DR_HOOK('MAIN_WOOD_SPINUP',1,ZHOOK_HANDLE)
! --------------------------------------------------------------------------------------
!
END PROGRAM MAIN_WOOD_SPINUP
