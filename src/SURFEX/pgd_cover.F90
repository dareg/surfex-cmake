
!     #########
      SUBROUTINE PGD_COVER(HPROGRAM)
!     ##############################################################
!
!!**** *PGD_COVER* monitor for averaging and interpolations of cover fractions
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    10/12/97
!!    B. Decharme  06/2008  limit of coast coverage under which the coast is replaced by sea or inland water
!!    B. Decharme  06/2009  remove lack and sea as the user want
!!    B. Decharme  07/2009  compatibility between Surfex and Orca (Nemo) grid (Earth Model)
!!    B. Decharme  07/2012  if sea or water imposed to 1 in a grid cell: no extrapolation
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_PGD_GRID,       ONLY : CGRID, NL, XGRID_PAR, NGRID_PAR
USE MODD_PGDWORK,        ONLY : XSUMCOVER, NSIZE
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER, NROCK, NSEA, NWATER, NPERMSNOW
USE MODD_DATA_COVER,     ONLY : XDATA_TOWN, XDATA_SEA, XDATA_NATURE, XDATA_WATER
USE MODD_SURF_ATM_n,      ONLY : CNATURE, CSEA, CTOWN, CWATER,            &
                                  XNATURE, XSEA, XTOWN, XWATER,           &
                                  XCOVER, LCOVER,                         &
                                  NSIZE_NATURE, NSIZE_SEA,                &
                                  NSIZE_TOWN, NSIZE_WATER,NSIZE_FULL,     &
                                  NDIM_NATURE, NDIM_SEA,                  &
                                  NDIM_TOWN,NDIM_WATER  
!
USE MODI_GET_LUOUT
USE MODE_GRIDTYPE_GAUSS

USE MODI_TREAT_FIELD
USE MODI_INTERPOL_FIELD2D
USE MODI_CONVERT_COVER_FRAC
!
USE MODI_READ_LCOVER
USE MODI_READ_SURF
USE MODI_SUM_ON_ALL_PROCS
!
USE MODI_READ_NAM_PGD_COVER
!
USE MODI_INIT_IO_SURF_n
USE MODI_END_IO_SURF_n
!
USE MODI_ABOR1_SFX
!
USE MODI_PGD_ECOCLIMAP2_DATA
!
#ifdef ASC
USE MODD_IO_SURF_ASC, ONLY : CFILEIN
#endif
#ifdef FA
USE MODD_IO_SURF_FA,  ONLY : CFILEIN_FA
#endif
#ifdef LFI
USE MODD_IO_SURF_LFI, ONLY : CFILEIN_LFI
#endif
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
 CHARACTER(LEN=6),    INTENT(IN)    :: HPROGRAM     ! Type of program
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
 CHARACTER(LEN=10)       :: YFIELD
 CHARACTER(LEN=28)       :: YCOVER      ! file name for cover types
 CHARACTER(LEN=6)        :: YFILETYPE   ! data file type
!
REAL                     :: XRM_COVER   ! limit of coverage under which the
                                        ! cover is removed. Default is 1.E-6
REAL                     :: XRM_COAST   ! limit of coast coverage under which
                                        ! the coast is replaced by sea or
                                        ! inland water. Default is 1.
REAL                     :: XRM_LAKE    ! limit of inland lake coverage under which
                                        ! the water is removed. Default is 0.0                     
REAL                     :: XRM_SEA     ! limit of sea coverage under which
                                        ! the sea is removed. Default is 0.0
REAL                     :: XLAT_ANT    ! Lattitude limit from Orca grid (Antartic)
!
REAL, DIMENSION(:), ALLOCATABLE :: ZDEF
REAL, DIMENSION(:), ALLOCATABLE :: ZLAT
REAL, DIMENSION(:), ALLOCATABLE :: XUNIF_COVER ! value of each cover (cover will be
!                                                uniform on the horizontal)
REAL, DIMENSION(:), ALLOCATABLE :: ZSEA   !to check compatibility between 
REAL, DIMENSION(:), ALLOCATABLE :: ZWATER !prescribed fractions and ECOCLIMAP
REAL, DIMENSION(:), ALLOCATABLE :: ZNATURE
REAL, DIMENSION(:), ALLOCATABLE :: ZTOWN
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCOVER_NATURE, ZCOVER_TOWN, ZCOVER_SEA, ZCOVER_WATER, ZCOVER
!
!
INTEGER               :: ILUOUT    ! output listing logical unit
INTEGER               :: IRESP     ! Error code after redding
INTEGER               :: JCOVER    ! loop counter on covers
INTEGER               :: JL        ! loop counter on horizontal points
INTEGER               :: ICOVER, ICOVERSUM, ICOVER_OLD, ICPT  ! 0 if cover is not present, >1 if present somewhere
INTEGER               :: IPERMSNOW, IECO2 
INTEGER               :: IC_NAT, IC_TWN, IC_WAT, IC_SEA
!
INTEGER, DIMENSION(1) :: IMAXCOVER ! index of maximum cover for the given point
INTEGER, DIMENSION(:), POINTER :: IMASK_COVER
INTEGER, DIMENSION(:), ALLOCATABLE :: IMASK_SEA, IMASK_WATER
!
LOGICAL                  :: LORCA_GRID  ! flag to compatibility between Surfex and Orca grid 
                                        ! (Earth Model over Antarctic)
LOGICAL                  :: LIMP_COVER  ! Imposed values for Cover from another PGD file
!
LOGICAL                  :: GPRESENT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!---------------------------------------------------------------
!
!*    1.      Initializations
!             ---------------
!
IF (LHOOK) CALL DR_HOOK('PGD_COVER',0,ZHOOK_HANDLE)
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
ALLOCATE(LCOVER     (JPCOVER))
ALLOCATE(XUNIF_COVER(JPCOVER))
!
LCOVER = .FALSE.
XUNIF_COVER = XUNDEF
!-------------------------------------------------------------------------------
!
!*    2.      Input file for cover types
!             --------------------------
!
 CALL READ_NAM_PGD_COVER(HPROGRAM, YCOVER, YFILETYPE, XUNIF_COVER,  &
                          XRM_COVER, XRM_COAST, XRM_LAKE, XRM_SEA,   &
                          LORCA_GRID, XLAT_ANT, LIMP_COVER           )  
!
!-------------------------------------------------------------------------------
!
!*    3.      Uniform field is prescribed
!             ---------------------------
!-------------------------------------------------------------------------------
!
IF (ANY(XUNIF_COVER/=0.)) THEN
!
!*    3.1     Verification of the total input cover fractions
!             -----------------------------------------------
!
  IF (ABS(SUM(XUNIF_COVER)-1.)>1.E-6) THEN
    WRITE(ILUOUT,*) ' '
    WRITE(ILUOUT,*) '***************************************************'
    WRITE(ILUOUT,*) '* Error in COVER fractions preparation            *'
    WRITE(ILUOUT,*) '* The prescribed covers does not fit              *'
    WRITE(ILUOUT,*) '* The sum of all cover must be equal to 1 exactly *'
    WRITE(ILUOUT,*) '***************************************************'
    WRITE(ILUOUT,*) ' '
    CALL ABOR1_SFX('PGD_COVER: SUM OF ALL COVER FRACTIONS MUST BE 1.')
!
!*    3.2     Use of the presribed cover fractions
!             ------------------------------------
!
  ELSE
    ICOVER = COUNT(XUNIF_COVER(:)/=0.)
    ALLOCATE(XCOVER(NL,ICOVER))
    ICPT = 0
    DO JCOVER=1,JPCOVER
      IF (XUNIF_COVER(JCOVER)/=0.) THEN
        LCOVER(JCOVER) = .TRUE.
        ICPT = ICPT + 1
        XCOVER(:,ICPT) = XUNIF_COVER(JCOVER)
      ENDIF
    END DO
    XCOVER(:,:)=XCOVER(:,:)/SPREAD(SUM(XCOVER(:,:),2),2,ICOVER)
  END IF
!
!*    3.3     No data
!             -------
!
ELSEIF (LEN_TRIM(YCOVER)==0) THEN
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) '***********************************************************'
  WRITE(ILUOUT,*) '* Error in COVER fractions preparation                    *'
  WRITE(ILUOUT,*) '* There is no prescribed cover fraction and no input file *'
  WRITE(ILUOUT,*) '***********************************************************'
  WRITE(ILUOUT,*) ' '
  CALL ABOR1_SFX('PGD_COVER: NO PRESCRIBED COVER NOR INPUT FILE')
!
!-------------------------------------------------------------------------------
ELSEIF(LIMP_COVER)THEN !LIMP_COVER (impose cover from input file at the same resolution)
!
  IF(YFILETYPE=='NETCDF')THEN
    CALL ABOR1_SFX('Use another format than netcdf for cover input file with LIMP_COVER')
  ELSE
#ifdef ASC
    CFILEIN     = ADJUSTL(ADJUSTR(YCOVER)//'.txt')
#endif
#ifdef FA
    CFILEIN_FA  = ADJUSTL(ADJUSTR(YCOVER)//'.fa')
#endif
#ifdef LFI
    CFILEIN_LFI = ADJUSTL(YCOVER)
#endif
    CALL INIT_IO_SURF_n(YFILETYPE,'FULL  ','SURF  ','READ ')
  ENDIF
!
  ALLOCATE(LCOVER(JPCOVER))
  CALL READ_LCOVER(YFILETYPE,LCOVER)
!
  CALL READ_SURF(YFILETYPE,'COVER',XCOVER(:,:),LCOVER,IRESP)
!
  CALL END_IO_SURF_n(YFILETYPE)
!
ELSE 
!-------------------------------------------------------------------------------
!
!*    3.      Averages the field
!             ------------------
!
  ALLOCATE(NSIZE     (NL)        )
  ALLOCATE(XSUMCOVER (NL,JPCOVER))
!
  NSIZE    (:)   = 0.
  XSUMCOVER(:,:) = 0.
  CALL TREAT_FIELD(HPROGRAM,'SURF  ',YFILETYPE,'A_COVR',YCOVER,  &
                     'COVER               '                      ) 

!
!*    4.      Interpolation if some points are not initialized (no data for these points) (same time)
!             ---------------------------------------------------------------------------------------
!
  WRITE(YFIELD,FMT='(A)') 'covers'
  CALL INTERPOL_FIELD2D(HPROGRAM,ILUOUT,NSIZE,XCOVER(:,:),YFIELD)
!
!-------------------------------------------------------------------------------
!
!*    5.      Coherence check
!             ---------------
!
  ICOVER = SIZE(XCOVER,2)
!
  XCOVER(:,:)=XCOVER(:,:)/SPREAD(SUM(XCOVER(:,:),2),2,ICOVER)
!
  DEALLOCATE(NSIZE    )
  DEALLOCATE(XSUMCOVER)
!
  CALL MAKE_MASK_COVER(IMASK_COVER,ICOVER)
!
  ALLOCATE(IMASK_SEA(SIZE(NSEA)))
  IMASK_SEA(:) = 0
  DO JL=1,SIZE(NSEA)
    DO JCOVER=1,ICOVER
      IF (IMASK_COVER(JCOVER)==NSEA(JL)) IMASK_SEA(JL) = JCOVER
    ENDDO
  ENDDO
  !
  ALLOCATE(IMASK_WATER(SIZE(NWATER)))
  IMASK_WATER(:) = 0
  DO JL=1,SIZE(NWATER)
    DO JCOVER=1,ICOVER
      IF (IMASK_COVER(JCOVER)==NWATER(JL)) IMASK_WATER(JL) = JCOVER
    ENDDO
  ENDDO
  !
  IPERMSNOW=0
  DO JCOVER=1,ICOVER
    IF (IMASK_COVER(JCOVER)==NPERMSNOW) IPERMSNOW = JCOVER
  ENDDO
  !
  IECO2 = 0
  DO JCOVER=1,ICOVER
    IF (IMASK_COVER(JCOVER)>300) THEN
      IECO2 = JCOVER
      EXIT
    ENDIF
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*    6.      Special treatments asked by user
!             --------------------------------
!
! * removes cover with very small coverage
  DO JL=1,SIZE(XCOVER,1)
    IMAXCOVER(:) = MAXLOC(XCOVER(JL,:))
    DO JCOVER=1,ICOVER
      IF (XCOVER(JL,JCOVER)/=0.) THEN
        IF (XCOVER(JL,JCOVER)<=XRM_COVER .AND. JCOVER /= IMAXCOVER(1)) THEN
          XCOVER(JL,JCOVER) = 0.
        END IF
      ENDIF
    END DO
  END DO
  !
  ! * removes lake as the user want
  IF(XRM_LAKE>0.0)THEN
     DO JL=1,SIZE(NWATER)
       IF (IMASK_WATER(JL)/=0) THEN
         WHERE(XCOVER(:,IMASK_WATER(JL))<=XRM_LAKE)
           XCOVER(:,IMASK_WATER(JL)) = 0.
         ENDWHERE
       ENDIF
     ENDDO          
  ENDIF
  !
  ! * removes sea as the user want
  IF(XRM_SEA>0.0)THEN
     DO JL=1,SIZE(NSEA)
       IF (IMASK_SEA(JL)/=0) THEN
         WHERE(XCOVER(:,IMASK_SEA(JL))<=XRM_SEA)
           XCOVER(:,IMASK_SEA(JL)) = 0.
         ENDWHERE
       ENDIF
     ENDDO          
  ENDIF
  !
  !
  ! * removes cover; replace by sea or inland water if sea or inland water > XRM_COAST
  DO JCOVER=1,ICOVER
    !
    DO JL=1,SIZE(NSEA)
      IF (IMASK_SEA(JL)/=0) THEN
        WHERE(XCOVER(:,IMASK_SEA(JL))>=XRM_COAST)
          XCOVER(:,JCOVER) = 0.
          XCOVER(:,IMASK_SEA(JL)) = 1.
        END WHERE 
      ENDIF
    ENDDO
    !
    DO JL=1,SIZE(NWATER)
      IF (IMASK_WATER(JL)/=0) THEN
        WHERE(XCOVER(:,IMASK_WATER(JL))>=XRM_COAST)
          XCOVER(:,JCOVER) = 0.
          XCOVER(:,IMASK_WATER(JL)) = 1.
        END WHERE
      ENDIF
    ENDDO 
    !    
  ENDDO
!
!
! * Compatibility between Surfex and Orca grid 
!   (Earth Model over water bodies and Antarctic)
!
  IF(LORCA_GRID.AND.CGRID=='GAUSS     ')THEN
!
!     No river or inland water bodies
    IF (IMASK_WATER(2)/=0) XCOVER(:,IMASK_WATER(2)) = 0.
    IF (IMASK_WATER(3)/=0) XCOVER(:,IMASK_WATER(3)) = 0.
!
    ALLOCATE(ZLAT(NL))
    CALL GET_GRIDTYPE_GAUSS(XGRID_PAR,PLAT=ZLAT)
!
      DO JL=1,SIZE(NSEA)
        IF (IMASK_SEA(JL)/=0) THEN
          WHERE(ZLAT(:)<XLAT_ANT.AND.XCOVER(:,IMASK_SEA(JL))>0.0)
            XCOVER(:,IPERMSNOW) = 1.0
            XCOVER(:,IMASK_SEA(JL))  = 0.0
          ENDWHERE 
        ENDIF
      ENDDO

      DO JL=1,SIZE(NWATER)
        IF (IMASK_WATER(JL)/=0) THEN
          WHERE(ZLAT(:)<XLAT_ANT.AND.XCOVER(:,IMASK_WATER(JL))>0.0)
            XCOVER(:,IPERMSNOW)  = 1.0
            XCOVER(:,IMASK_WATER(JL)) = 0.0
          ENDWHERE
        ENDIF
      ENDDO
!  
    DEALLOCATE(ZLAT)
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*    7.      Coherence check
!             ---------------
!
  XCOVER(:,:)=XCOVER(:,:)/SPREAD(SUM(XCOVER(:,:),2),2,ICOVER)
!
  DEALLOCATE(IMASK_SEA)
  DEALLOCATE(IMASK_WATER)
!
!*    8.      List of cover present
!             ---------------------
!
  ALLOCATE(ZCOVER(NL,ICOVER))
  ZCOVER(:,:) = 0.
!
  ICOVER_OLD = ICOVER
  ICOVER = 0
!
  IECO2 = 0
!
  LCOVER(:) = .FALSE.
  DO JCOVER=1,ICOVER_OLD
    ICOVERSUM = SUM_ON_ALL_PROCS(HPROGRAM,CGRID,XCOVER(:,JCOVER)/=0., 'COV')
    IF (ICOVERSUM>0) THEN
      LCOVER(IMASK_COVER(JCOVER))=.TRUE.
      ICOVER = ICOVER+1
      ZCOVER(:,ICOVER) = XCOVER(:,JCOVER)
      IF (IMASK_COVER(JCOVER)>300) IECO2 = ICOVER
    ENDIF
  END DO
!
  DEALLOCATE(XCOVER)
  ALLOCATE(XCOVER(NL,ICOVER))
  XCOVER(:,:) = ZCOVER(:,1:ICOVER)
!
  DEALLOCATE(ZCOVER)
  DEALLOCATE(IMASK_COVER)
!
!-------------------------------------------------------------------------------
END IF
!
DEALLOCATE(XUNIF_COVER)
!-------------------------------------------------------------------------------
!
!
IF(.NOT.LIMP_COVER)THEN
        
!*    8.      List of cover present
!             ---------------------
!
  IF (IECO2/=0) THEN
    IF ( SUM_ON_ALL_PROCS(HPROGRAM,CGRID,ANY(XCOVER(:,IECO2:)>0.,DIM=2),'COV' ) >0 ) &
      CALL PGD_ECOCLIMAP2_DATA(HPROGRAM)
  ENDIF
!
!-------------------------------------------------------------------------------
ENDIF
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*    9.      Land - sea fractions
!             --------------------
!
IF (.NOT.ASSOCIATED(XSEA)) THEN

  ALLOCATE(XSEA   (NL))
  ALLOCATE(XWATER (NL))
  ALLOCATE(XNATURE(NL))
  ALLOCATE(XTOWN  (NL))
  CALL CONVERT_COVER_FRAC(XCOVER,LCOVER,XSEA,XNATURE,XTOWN,XWATER)

ELSE
  !
  ICOVER = SIZE(XCOVER,2)
  !
  CALL MAKE_MASK_COVER(IMASK_COVER,ICOVER)
  !
!if fractions are prescribed, it has to be verified that the locations of
!ECOCLIMAP covers are compatible with the fractions of surface types
  ALLOCATE(ZSEA   (NL))
  ALLOCATE(ZWATER (NL))
  ALLOCATE(ZNATURE(NL))
  ALLOCATE(ZTOWN  (NL))
  CALL CONVERT_COVER_FRAC(XCOVER,LCOVER,ZSEA,ZNATURE,ZTOWN,ZWATER)
  !
  CALL FIT_COVERS(XDATA_NATURE,XNATURE,4,ICOVER,IC_NAT)
  CALL FIT_COVERS(XDATA_TOWN,XTOWN,7,ICOVER,IC_TWN)
  CALL FIT_COVERS(XDATA_WATER,XWATER,2,ICOVER,IC_WAT)
  CALL FIT_COVERS(XDATA_SEA,XSEA,1,ICOVER,IC_SEA)
  !
  ALLOCATE(ZCOVER_NATURE(NL,ICOVER))
  ALLOCATE(ZCOVER_TOWN  (NL,ICOVER))
  ALLOCATE(ZCOVER_SEA   (NL,ICOVER))
  ALLOCATE(ZCOVER_WATER (NL,ICOVER))
  !
  ZCOVER_NATURE(:,:) = XCOVER(:,:)
  ZCOVER_TOWN  (:,:) = XCOVER(:,:)
  ZCOVER_SEA   (:,:) = XCOVER(:,:)
  ZCOVER_WATER (:,:) = XCOVER(:,:)
  !
  ALLOCATE(NSIZE(NL))
  !
  ALLOCATE(ZDEF(ICOVER))
  !
  WRITE(ILUOUT,FMT=*) &
  '*********************************************************************'
  WRITE(ILUOUT,FMT=*) &
  '*  Coherence computation between covers and imposed nature fraction *'
  WRITE(ILUOUT,FMT=*) &
  '*********************************************************************'
  NSIZE(:) = 1
  WHERE (XNATURE(:).NE.0. .AND. ZNATURE(:).EQ.0.) NSIZE(:)=0
          
  DO JL=1,SIZE(XCOVER,1)
    IF (XNATURE(JL).EQ.0.) NSIZE(JL)=-1
  ENDDO
  ZDEF(:)=0.
  DO JCOVER=1,ICOVER
    IF (XDATA_NATURE(IMASK_COVER(JCOVER))/=0.) THEN
      ZDEF(JCOVER) = 1.
      EXIT
    ENDIF
  ENDDO
  CALL INTERPOL_FIELD2D(HPROGRAM,ILUOUT,NSIZE,ZCOVER_NATURE(:,:),YFIELD,ZDEF)  
!
  WRITE(ILUOUT,FMT=*) &
  '*********************************************************************'
  WRITE(ILUOUT,FMT=*) &
  '*  Coherence computation between covers and imposed town   fraction *'
  WRITE(ILUOUT,FMT=*) &
  '*********************************************************************'
  NSIZE(:) = 1
  WHERE (XTOWN(:).NE.0. .AND. ZTOWN(:).EQ.0.) NSIZE(:)=0
  DO JL=1,SIZE(XCOVER,1)
    IF (XTOWN(JL).EQ.0.) NSIZE(JL)=-1
  ENDDO
  ZDEF(:)=0.
  DO JCOVER=1,ICOVER
    IF (XDATA_TOWN(IMASK_COVER(JCOVER))/=0.) THEN
      ZDEF(JCOVER) = 1.
      EXIT
    ENDIF
  ENDDO  
  CALL INTERPOL_FIELD2D(HPROGRAM,ILUOUT,NSIZE,ZCOVER_TOWN (:,:),YFIELD,ZDEF) 

  WRITE(ILUOUT,FMT=*) &
  '*********************************************************************'
  WRITE(ILUOUT,FMT=*) &
  '*  Coherence computation between covers and imposed water  fraction *'
  WRITE(ILUOUT,FMT=*) &
  '*********************************************************************'
  NSIZE(:) = 1
  WHERE (XWATER(:).NE.0. .AND. ZWATER(:).EQ.0.) NSIZE(:)=0
! if water imposed to 1 in a grid cell: no extrapolation
  DO JL=1,SIZE(XCOVER,1)
     IF(XWATER(JL)==1.0)THEN
        ZCOVER_WATER(JL,:)=0.0             
        ZCOVER_WATER(JL,IC_WAT)=1.0
        NSIZE(JL)=1
     ELSEIF(XWATER(JL)==0.0)THEN
        NSIZE(JL)=-1
     ENDIF
  ENDDO
  ZDEF(:)=0.
  DO JCOVER=1,ICOVER
    IF (XDATA_WATER(IMASK_COVER(JCOVER))/=0.) THEN
      ZDEF(JCOVER) = 1.
      EXIT
    ENDIF
  ENDDO    
  CALL INTERPOL_FIELD2D(HPROGRAM,ILUOUT,NSIZE,ZCOVER_WATER (:,:),YFIELD,PDEF=ZDEF)
  WRITE(ILUOUT,FMT=*) &
  '*********************************************************************'
  WRITE(ILUOUT,FMT=*) &
  '*  Coherence computation between covers and imposed sea    fraction *'
  WRITE(ILUOUT,FMT=*) &
  '*********************************************************************'
  NSIZE(:) = 1
  WHERE (XSEA(:).NE.0. .AND. ZSEA(:).EQ.0.) NSIZE(:)=0
! if sea imposed to 1 in a grid cell: no extrapolation          
  DO JL=1,SIZE(XCOVER,1)
     IF(XSEA(JL)==1.0)THEN
        ZCOVER_SEA(JL,:)=0.0             
        ZCOVER_SEA(JL,IC_SEA)=1.0
        NSIZE(JL)=1
     ELSEIF(XSEA(JL)==0.0)THEN
        NSIZE(JL)=-1
     ENDIF
  ENDDO
  ZDEF(:)=0.
  DO JCOVER=1,ICOVER
    IF (XDATA_SEA(IMASK_COVER(JCOVER))/=0.) THEN
      ZDEF(JCOVER) = 1.
      EXIT
    ENDIF
  ENDDO    
  CALL INTERPOL_FIELD2D(HPROGRAM,ILUOUT,NSIZE,ZCOVER_SEA (:,:),YFIELD,PDEF=ZDEF)
  !
  XCOVER(:,:) = XCOVER(:,:) + 0.001 * ( ZCOVER_NATURE(:,:) + ZCOVER_TOWN(:,:) + &
                                        ZCOVER_WATER (:,:) + ZCOVER_SEA (:,:) )
  !
  XCOVER(:,:)=XCOVER(:,:)/SPREAD(SUM(XCOVER(:,:),2),2,ICOVER)
  !
  DEALLOCATE(ZCOVER_NATURE)
  DEALLOCATE(ZCOVER_TOWN  )
  DEALLOCATE(ZCOVER_WATER )
  DEALLOCATE(ZCOVER_SEA   )
  !
  DEALLOCATE(NSIZE    )
  DEALLOCATE(ZSEA     )
  DEALLOCATE(ZWATER   )
  DEALLOCATE(ZNATURE  )
  DEALLOCATE(ZTOWN    )
  !
  DEALLOCATE(ZDEF)
  DEALLOCATE(IMASK_COVER)
  !
ENDIF
!
NSIZE_NATURE    = COUNT(XNATURE(:) > 0.0)
NSIZE_WATER     = COUNT(XWATER (:) > 0.0)
NSIZE_SEA       = COUNT(XSEA   (:) > 0.0)
NSIZE_TOWN      = COUNT(XTOWN  (:) > 0.0)
NSIZE_FULL      = NL
!
NDIM_NATURE    = SUM_ON_ALL_PROCS(HPROGRAM,CGRID,XNATURE(:) > 0., 'DIM')
NDIM_WATER     = SUM_ON_ALL_PROCS(HPROGRAM,CGRID,XWATER (:) > 0., 'DIM')
NDIM_SEA       = SUM_ON_ALL_PROCS(HPROGRAM,CGRID,XSEA   (:) > 0., 'DIM')
NDIM_TOWN      = SUM_ON_ALL_PROCS(HPROGRAM,CGRID,XTOWN  (:) > 0., 'DIM')
!
IF (LHOOK) CALL DR_HOOK('PGD_COVER',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
CONTAINS
!
SUBROUTINE FIT_COVERS(PDATA_SURF,PSURF,KSURF,KCOVER,KC_SURF)
!
REAL, DIMENSION(:), INTENT(IN) :: PDATA_SURF
REAL, DIMENSION(:), INTENT(IN) :: PSURF
INTEGER, INTENT(IN) :: KSURF
INTEGER, INTENT(INOUT) :: KCOVER
INTEGER, INTENT(OUT) :: KC_SURF
!
LOGICAL :: GPRESENT
!
GPRESENT = .FALSE.
DO JCOVER=1,KCOVER
  IF (PDATA_SURF(IMASK_COVER(JCOVER))/=0.) THEN
    GPRESENT = .TRUE.
    EXIT
  ENDIF
ENDDO
!
IF (.NOT.GPRESENT .AND. ANY(PSURF(:)/=0.)) THEN
  !
  LCOVER(KSURF) = .TRUE.
  KCOVER = KCOVER + 1
  ALLOCATE(ZCOVER(NL,KCOVER))
  DO JCOVER = 1,KCOVER
    IF (JCOVER<KCOVER) THEN
      IF (IMASK_COVER(JCOVER)<KSURF) CYCLE
    ENDIF
    KC_SURF = JCOVER
    IF (JCOVER>1) ZCOVER(:,1:JCOVER-1) = XCOVER(:,1:JCOVER-1)
    ZCOVER(:,JCOVER) = 0.
    IF (JCOVER<KCOVER) ZCOVER(:,JCOVER+1:KCOVER) = XCOVER(:,JCOVER:KCOVER-1)
    EXIT
  ENDDO
  DEALLOCATE(XCOVER)
  ALLOCATE(XCOVER(NL,KCOVER))
  XCOVER(:,:) = ZCOVER(:,:)
  DEALLOCATE(ZCOVER)
  !
  CALL MAKE_MASK_COVER(IMASK_COVER,KCOVER)
  !
ENDIF
!
END SUBROUTINE FIT_COVERS
!
!------------------------------------------------------
!
SUBROUTINE MAKE_MASK_COVER(KMASK_COVER,KCOVER)
!
INTEGER, DIMENSION(:), POINTER :: KMASK_COVER
INTEGER, INTENT(IN) :: KCOVER
!
INTEGER :: ICPT
!
IF (ASSOCIATED(KMASK_COVER)) DEALLOCATE(KMASK_COVER)
ALLOCATE(KMASK_COVER(KCOVER))
ICPT = 0
DO JCOVER=1,JPCOVER
  IF (LCOVER(JCOVER)) THEN
    ICPT = ICPT + 1
    KMASK_COVER(ICPT) = JCOVER
  ENDIF
ENDDO
!
END SUBROUTINE MAKE_MASK_COVER
!
END SUBROUTINE PGD_COVER
