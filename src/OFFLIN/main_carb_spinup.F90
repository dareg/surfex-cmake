!     ######spl
PROGRAM MAIN_CARB_SPINUP
!
!
!!****  *MAIN_CARB_SPINUP*  
!!
!!    PURPOSE
!!    -------
!!    Spinup of soil carbon.
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
!!      Original    02/07/09
!!      S.Lafont : computation only on the vegetation patch : from 4 to 12
!!      S.Lafont :  add test to see if the PATCH is present.
!!      B. Decharme : read pgd+prep
!!
!-------------------------------------------------------------------------------
!
USE MODD_IO_SURF_ASC,ONLY : CFILEIN,CFILEIN_SAVE,CFILEOUT,CFILEPGD
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
!!!
!!!
USE MODE_POS_SURF
USE MODE_SOIL
!!!
USE MODI_CARBON_INIT
USE MODI_CONTROL_MOIST_FUNC
USE MODI_CONTROL_TEMP_FUNC
USE MODI_CARBON_LITTER
USE MODI_CARBON_SOIL
USE MODI_INI_CSTS
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_OL_READ_ATM_CONF
USE MODI_INIT_IO_SURF_n
USE MODI_TEST_NAM_VAR_SURF
USE MODI_READ_SURF
USE MODI_WRITE_SURF
USE MODI_CLOSE_FILEOUT_OL
!
USE MODI_SET_SURFEX_FILEIN
!
!!!
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
!
USE MODI_ALLOC_SURFEX
USE MODI_DEALLOC_SURFEX
USE MODI_END_IO_SURF_n
USE MODI_GET_LUOUT
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.    declarations of local variables
!
INTEGER                           :: NPATCH              ! maximum number of sub-tiles (patches)
!                                                        ! used at any grid point within a 
INTEGER                           :: NGROUND_LAYER       ! number of ground layers
INTEGER                           :: NNBIOMASS           ! number of biomass pools
INTEGER                           :: NNLITTER            ! number of litter pools
INTEGER                           :: NNLITTLEVS          ! number of litter levels
INTEGER                           :: NNSOILCARB          ! number of soil carbon pools
!
INTEGER                           :: ILUOUT              ! ascii output unit number
INTEGER                           :: ILUNAM              ! namelist unit number
CHARACTER(LEN=28)                 :: YLUOUT    = 'LISTING_CARB                '
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
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4)  :: YLVL
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=25) :: YFORM          ! Writing format
INTEGER           :: IRESP          ! error return code
INTEGER           :: JSTEP, JLAYER, JNBIOMASS, JPATCH, JNLITTER, JNLITTLEVS, JNSOILCARB, NVEGTYPE
!!!!
CHARACTER(LEN=4)                      :: CPEDOTF  ! type of pedo-transfert
CHARACTER(LEN=3)                      :: CPHOTO   ! type of photosynthesis
CHARACTER(LEN=3)                      :: CRESPSL  ! type of soil respiration
!!!!
REAL                                  :: ZDTSTEP
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSOILTEMP_CURR
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSOILMOIST_CURR
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZTURNOVER_CURR
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZCLAY, ZSAND, ZWWILT, ZWFC, ZWSAT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZLITTER, ZLITTER_LAST
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZLIGNIN_STRUC, ZLIGNIN_STRUC_LAST
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZSOILCARB, ZSOILCARB_LAST
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZCONTROL_MOIST, ZCONTROL_TEMP
REAL, DIMENSION(:),       ALLOCATABLE :: ZRESP_HETERO_DAY_LITTER, ZRESP_HETERO_DAY_SOIL
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZSOILCARBON_INPUT
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZVEGTYPE_PATCH
REAL                                  :: ZTEST
REAL(KIND=JPRB) :: ZHOOK_HANDLE


! --------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MAIN_CARB_SPINUP',0,ZHOOK_HANDLE)
CALL ALLOC_SURFEX(1)
CALL GOTO_SURFEX(1,.TRUE.)
NVEGTYPE = 12
!
!*      1.    Reading of inputs
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
!*      1.1   Namelist
!

CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
!
CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND,ILUOUT)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
!
CALL TEST_NAM_VAR_SURF(ILUOUT,'CSURF_FILETYPE',CSURF_FILETYPE,'ASCII ','LFI   ','FA    ')
!IF (CSURF_FILETYPE/='ASCII ') THEN
!  WRITE(ILUOUT,*) '*****************************************'
!  WRITE(ILUOUT,*) '* Carbon spinup is only implemented for CSURF_FILETYPE=''ASCII'' *'
!  WRITE(ILUOUT,*) '* and CSURF_FILETYPE = ',CSURF_FILETYPE,' in this run            *'
!  WRITE(ILUOUT,*) '*****************************************'
!  CALL ABOR1_SFX('MAIN_CARB_SPINUP: CARBON SPINUP REQUIRES CSURF_FILETYPE=''ASCII''')
!ENDIF

CALL TEST_NAM_VAR_SURF(ILUOUT,'CTIMESERIES_FILETYPE',CTIMESERIES_FILETYPE, &
       'NETCDF','TEXTE ','BINARY','ASCII ','LFI   ','FA    ','NONE  ','OFFLIN')  
IF (CTIMESERIES_FILETYPE=='NETCDF') CTIMESERIES_FILETYPE='OFFLIN'

CALL CLOSE_NAMELIST('ASCII ',ILUNAM)


!
!*      1.2   Netcdf file handling
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
!*      1.3   Dimensions
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
CALL READ_SURF(CSURF_FILETYPE,'GROUND_LAYER',NGROUND_LAYER,IRESP)
CALL READ_SURF(CSURF_FILETYPE,'NBIOMASS',NNBIOMASS,IRESP)
CALL READ_SURF(CSURF_FILETYPE,'NLITTER',NNLITTER,IRESP)
CALL READ_SURF(CSURF_FILETYPE,'NLITTLEVS',NNLITTLEVS,IRESP)
CALL READ_SURF(CSURF_FILETYPE,'NSOILCARB',NNSOILCARB,IRESP)
CALL READ_SURF(CSURF_FILETYPE,'PEDOTF',CPEDOTF,IRESP)
CALL READ_SURF(CSURF_FILETYPE,'PHOTO',CPHOTO,IRESP)

CALL READ_SURF(CSURF_FILETYPE,'RESPSL',CRESPSL,IRESP)
IF (CRESPSL/='CNT') THEN
  WRITE(ILUOUT,*) '*****************************************************'
  WRITE(ILUOUT,*) '* Running soil carbon spinup does not have interest *'
  WRITE(ILUOUT,*) '* with option CRESPSL = ',CRESPSL,'                 *'
  WRITE(ILUOUT,*) '*****************************************************'
  CALL ABOR1_SFX('MAIN_CARB_SPINUP: DO NOT RUN SOIL CARBON SPINUP WITH CRESPSL='//CRESPSL)
ENDIF
!
! Fraction of each vegetation type for each patch
ALLOCATE(ZVEGTYPE_PATCH(INI,NPATCH,NVEGTYPE))
!
DO JPATCH=1,NPATCH
  WRITE(YLVL,'(I2)') JPATCH 
  YRECFM='VEGTYPE_P'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(CSURF_FILETYPE,YRECFM,ZVEGTYPE_PATCH(:,JPATCH,:),IRESP)
END DO
!
! Clay and sand fractions
!
ALLOCATE(ZCLAY(INI,NGROUND_LAYER))
ALLOCATE(ZSAND(INI,NGROUND_LAYER))
!
YRECFM='CLAY'
CALL READ_SURF(CSURF_FILETYPE,YRECFM,ZCLAY(:,1),IRESP)
!
YRECFM='SAND'
CALL READ_SURF(CSURF_FILETYPE,YRECFM,ZSAND(:,1),IRESP)
!
DO JLAYER=2,NGROUND_LAYER
  ZCLAY(:,JLAYER)=ZCLAY(:,1)
  ZSAND(:,JLAYER)=ZSAND(:,1)
END DO
!
CALL END_IO_SURF_n(CSURF_FILETYPE)
CALL SET_SURFEX_FILEIN(CSURF_FILETYPE,'PREP') ! restore input file name
!
!-----------------------------------------------------------------------------------------------------
! END READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
!*      1.4   Model variables
!

ALLOCATE(ZSOILTEMP_CURR(INI,NGROUND_LAYER,NPATCH,NSTEP_OUTPUT))
ALLOCATE(ZSOILMOIST_CURR(INI,NGROUND_LAYER,NPATCH,NSTEP_OUTPUT))
ALLOCATE(ZTURNOVER_CURR(INI,NNBIOMASS,NPATCH,NSTEP_OUTPUT))
ALLOCATE(ZWWILT(INI,NGROUND_LAYER))
ALLOCATE(ZWFC(INI,NGROUND_LAYER))
ALLOCATE(ZWSAT(INI,NGROUND_LAYER))
ALLOCATE(ZLITTER(INI,NNLITTER,NNLITTLEVS,NPATCH))
ALLOCATE(ZLIGNIN_STRUC(INI,NNLITTLEVS,NPATCH))
ALLOCATE(ZSOILCARB(INI,NNSOILCARB,NPATCH))
ALLOCATE(ZLITTER_LAST(INI,NNLITTER,NNLITTLEVS,NPATCH))
ALLOCATE(ZLIGNIN_STRUC_LAST(INI,NNLITTLEVS,NPATCH))
ALLOCATE(ZSOILCARB_LAST(INI,NNSOILCARB,NPATCH))
ALLOCATE(ZCONTROL_MOIST(INI,NNLITTLEVS))
ALLOCATE(ZCONTROL_TEMP(INI,NNLITTLEVS))
ALLOCATE(ZRESP_HETERO_DAY_LITTER(INI))
ALLOCATE(ZRESP_HETERO_DAY_SOIL(INI))
ALLOCATE(ZSOILCARBON_INPUT(INI,NNSOILCARB))

! Soil temperature and soil liquid water content

DO JLAYER=1,NGROUND_LAYER

  WRITE(YLVL,'(I4)') JLAYER

  YRECFM='TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(CTIMESERIES_FILETYPE,YRECFM,ZSOILTEMP_CURR(:,JLAYER,:,:),IRESP)

  YRECFM='WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(CTIMESERIES_FILETYPE,YRECFM,ZSOILMOIST_CURR(:,JLAYER,:,:),IRESP)

END DO


! Turnover rates

DO JNBIOMASS=1,NNBIOMASS
    
  WRITE(YLVL,'(I4)') JNBIOMASS
  YRECFM='TURNOVER'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(CTIMESERIES_FILETYPE,YRECFM,ZTURNOVER_CURR(:,JNBIOMASS,:,:),IRESP)

END DO

DO JLAYER=1,NGROUND_LAYER
  ZWWILT(:,JLAYER) = WWILT_FUNC(ZCLAY(:,JLAYER),ZSAND(:,JLAYER),CPEDOTF)
  ZWFC  (:,JLAYER) = WFC_FUNC  (ZCLAY(:,JLAYER),ZSAND(:,JLAYER),CPEDOTF)
  ZWSAT (:,JLAYER) = WSAT_FUNC (ZCLAY(:,JLAYER),ZSAND(:,JLAYER),CPEDOTF)
END DO
!
!         Initialisation for IO
!
! Read NFULL before reading 2D arrays
CALL INIT_IO_SURF_n(CSURF_FILETYPE,'FULL  ','SURF  ','READ ')

CALL INIT_IO_SURF_n(CSURF_FILETYPE,'NATURE','SURF  ','READ ')

! Last values of soil carbon variables

DO JNLITTER=1,NNLITTER
  DO JNLITTLEVS=1,NNLITTLEVS

    WRITE(YLVL,'(I1,A1,I1)') JNLITTER,'_',JNLITTLEVS
    YRECFM='LITTER'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL READ_SURF(CSURF_FILETYPE,YRECFM,ZLITTER_LAST(:,JNLITTER,JNLITTLEVS,:),IRESP)

  END DO
END DO

DO JNSOILCARB=1,NNSOILCARB

  WRITE(YLVL,'(I4)') JNSOILCARB
  YRECFM='SOILCARB'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(CSURF_FILETYPE,YRECFM,ZSOILCARB_LAST(:,JNSOILCARB,:),IRESP)

END DO

DO JNLITTLEVS=1,NNLITTLEVS

  WRITE(YLVL,'(I4)') JNLITTLEVS
  YRECFM='LIGNIN_STR'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(CSURF_FILETYPE,YRECFM,ZLIGNIN_STRUC_LAST(:,JNLITTLEVS,:),IRESP)

END DO
!

CALL END_IO_SURF_n(CSURF_FILETYPE)


!
!*      2.    Initialisations
!             ---------------

! Initialise surface parameters

CALL INI_CSTS

! Initialise carbon parameters
IF (CRESPSL=='CNT' .AND. CPHOTO == 'NCB') THEN
   CALL CARBON_INIT(NNBIOMASS, NNLITTER, NNLITTLEVS, NNSOILCARB)
ENDIF

! Initialise variables

ZDTSTEP = XTSTEP_OUTPUT

ZLITTER(:,:,:,:) = ZLITTER_LAST(:,:,:,:)
ZLIGNIN_STRUC(:,:,:) = ZLIGNIN_STRUC_LAST(:,:,:)
ZSOILCARB(:,:,:) = ZSOILCARB_LAST(:,:,:)

ZCONTROL_TEMP(:,:) = 0.
ZCONTROL_MOIST(:,:) = 0.
ZRESP_HETERO_DAY_LITTER(:) = 0.
ZRESP_HETERO_DAY_SOIL(:) = 0.
ZSOILCARBON_INPUT(:,:) = 0.


!
!*      3.    Calculations
!             ------------
!

DO JSTEP=1,NSTEP_OUTPUT
  ! computation only on vegetation patch : start JPATCH start at 4
  DO JPATCH=1,NPATCH

    ZTEST=SUM(ZVEGTYPE_PATCH(:,JPATCH,JPATCH))
    IF (ZTEST>0.0) THEN
      ! Control of heterotrophic respiration by soil temperature and moisture

      ZCONTROL_TEMP(:,1) = CONTROL_TEMP_FUNC (ZSOILTEMP_CURR(:,1,JPATCH,JSTEP))
      ZCONTROL_TEMP(:,2) = CONTROL_TEMP_FUNC (ZSOILTEMP_CURR(:,2,JPATCH,JSTEP))
      ZCONTROL_MOIST(:,1) = &
          CONTROL_MOIST_FUNC (ZSOILMOIST_CURR(:,1,JPATCH,JSTEP),ZWWILT(:,1), &
                              ZWFC(:,1),ZWSAT(:,1))  
      ZCONTROL_MOIST(:,2) = &
          CONTROL_MOIST_FUNC (ZSOILMOIST_CURR(:,2,JPATCH,JSTEP),ZWWILT(:,2), &
                              ZWFC(:,2),ZWSAT(:,2))  
     
      ! Evolution of litter

      CALL CARBON_LITTER (ZDTSTEP,ZTURNOVER_CURR(:,:,JPATCH,JSTEP),              &
                          ZLITTER(:,:,:,JPATCH),ZLIGNIN_STRUC(:,:,JPATCH),       &
                          ZCONTROL_TEMP,ZCONTROL_MOIST,                          &
                          ZRESP_HETERO_DAY_LITTER,ZSOILCARBON_INPUT)  

      ! Evolution of soil carbon

      CALL CARBON_SOIL (ZDTSTEP,ZSAND(:,2),ZSOILCARBON_INPUT,                    &
                        ZCONTROL_TEMP,ZCONTROL_MOIST,                            &
                        ZSOILCARB(:,:,JPATCH),ZRESP_HETERO_DAY_SOIL)

    ENDIF  
  ENDDO

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
DO JNLITTER=1,NNLITTER
  DO JNLITTLEVS=1,NNLITTLEVS
    WRITE(YLVL,'(I1,A1,I1)') JNLITTER,'_',JNLITTLEVS
    YRECFM='LITTER'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YFORM='(A10,I1.1,A1,I1.1,A8)'
    WRITE(YCOMMENT,FMT=YFORM) 'X_Y_LITTER',JNLITTER,' ',JNLITTLEVS,' (gC/m2)'
    CALL WRITE_SURF(CSURF_FILETYPE,YRECFM,ZLITTER(:,JNLITTER,JNLITTLEVS,:),IRESP,HCOMMENT=YCOMMENT)
  END DO
END DO

DO JNSOILCARB=1,NNSOILCARB
  WRITE(YLVL,'(I4)') JNSOILCARB
  YRECFM='SOILCARB'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A8,I1.1,A8)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_SOILCARB',JNSOILCARB,' (gC/m2)'
  CALL WRITE_SURF(CSURF_FILETYPE,YRECFM,ZSOILCARB(:,JNSOILCARB,:),IRESP,HCOMMENT=YCOMMENT)
END DO

DO JNLITTLEVS=1,NNLITTLEVS
  WRITE(YLVL,'(I4)') JNLITTLEVS
  YRECFM='LIGNIN_STR'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A12,I1.1,A8)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_LIGNIN_STRUC',JNLITTLEVS,' (gC/m2)'
  CALL WRITE_SURF(CSURF_FILETYPE,YRECFM,ZLIGNIN_STRUC(:,JNLITTLEVS,:),IRESP,HCOMMENT=YCOMMENT)
END DO
!
CALL END_IO_SURF_n(CSURF_FILETYPE)


DEALLOCATE(ZSOILTEMP_CURR) 
DEALLOCATE(ZSOILMOIST_CURR) 
DEALLOCATE(ZTURNOVER_CURR) 
DEALLOCATE(ZCLAY)
DEALLOCATE(ZSAND)
DEALLOCATE(ZWWILT)
DEALLOCATE(ZWFC)
DEALLOCATE(ZWSAT)
DEALLOCATE(ZLITTER)
DEALLOCATE(ZLIGNIN_STRUC)
DEALLOCATE(ZSOILCARB)
DEALLOCATE(ZLITTER_LAST)
DEALLOCATE(ZLIGNIN_STRUC_LAST)
DEALLOCATE(ZSOILCARB_LAST)
DEALLOCATE(ZCONTROL_MOIST)
DEALLOCATE(ZCONTROL_TEMP)
DEALLOCATE(ZRESP_HETERO_DAY_LITTER)
DEALLOCATE(ZRESP_HETERO_DAY_SOIL)
DEALLOCATE(ZSOILCARBON_INPUT)
DEALLOCATE(ZVEGTYPE_PATCH)
!
!
!*      5.    Close parallelized I/O
!             ----------------------
!
IF (CTIMESERIES_FILETYPE=='OFFLIN') CALL CLOSE_FILEOUT_OL
!
CLOSE(ILUOUT)
!
WRITE(*,*) ' '
WRITE(*,*) '    -----------------------------------'
WRITE(*,*) '    | MAIN_CARB_SPINUP ENDS CORRECTLY |'
WRITE(*,*) '    -----------------------------------'
WRITE(*,*) ' '
CALL DEALLOC_SURFEX
IF (LHOOK) CALL DR_HOOK('MAIN_CARB_SPINUP',1,ZHOOK_HANDLE)
! --------------------------------------------------------------------------------------
!
END PROGRAM MAIN_CARB_SPINUP
