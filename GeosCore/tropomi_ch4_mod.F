!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tropomi_ch4_mod.F
!
! !DESCRIPTION: Module TROPOMI\_CH4\_MOD
!\\
!\\
! !INTERFACE: 
!
      MODULE TROPOMI_CH4_MOD
!
! !USES:
!
      USE m_netcdf_io_open       ! netCDF open
      USE m_netcdf_io_get_dimlen ! netCDF dimension queries
      USE m_netcdf_io_read       ! netCDF data reads
      USE m_netcdf_io_close      ! netCDF close
      USE PRECISION_MOD          ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE 
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CALC_TROPOMI_CH4_FORCE
!
! !REVISION HISTORY:
!  20 Jun 2019 - M. Sulprizio- Initial version based on gosat_ch4_mod.F 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Parameters
      INTEGER,  PARAMETER  :: MaxLev = 13
      INTEGER,  PARAMETER  :: MaxObs = 1e7
      INTEGER,  PARAMETER  :: NT_FD = 335 !1055 ! Index for the FD test
      REAL(fp), PARAMETER  :: GC_XCH4_ERROR = 0e+0_fp !ppb
      REAL(fp), PARAMETER  :: PMAG = 1e+0_fp ! Perturbation magnitude (%)
      LOGICAL,  PARAMETER  :: LDCH4SAT    = .TRUE.
!
! !MODULE VARIABLES
!
      ! Record to store data from each TROP obs
      TYPE TROP_CH4_OBS 
         INTEGER           :: Lev(1)
         REAL(fp)          :: Lat(1)
         REAL(fp)          :: Lon(1)
         INTEGER           :: Year(1)
         INTEGER           :: Month(1)
         INTEGER           :: Day(1)
         INTEGER           :: Hour(1)
         INTEGER           :: Minute(1)
         INTEGER           :: Second(1)
         INTEGER           :: Time(1)
         REAL(fp)          :: CH4(1)
         REAL(fp)          :: CH4_Error(1)
         REAL(fp)          :: QFlag(1)
         REAL(fp)          :: Pres(MaxLev)
         REAL(fp)          :: Prior(MaxLev)
         REAL(fp)          :: Avg_Kernel(MaxLev)
      ENDTYPE TROP_CH4_OBS  

      TYPE(TROP_CH4_OBS)    :: TROP(MaxObs)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_trop_ch4_obs
!
! !DESCRIPTION: Subroutine READ\_TROP\_CH4\_OBS reads the file and passes back
!  info contained therein. (dkh, 10/12/10) 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_TROP_CH4_OBS( YYYYMMDD, nTrop )
!
! !USES:
!
!
      USE TIME_MOD,        ONLY : EXPAND_DATE
      USE ERROR_MOD,       ONLY : ALLOC_ERR
      USE ERROR_MOD,       ONLY : GEOS_CHEM_STOP
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)      :: YYYYMMDD  ! Current date
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT)     :: nTrop      ! Number of TROPOMI retrievals
!     
! !REVISION HISTORY:
!  20 Jun 2019 - M. Sulprizio- Initial version based on gosat_ch4_mod.F 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: L, M, N, AS
      INTEGER              :: nLev, nScan, nPixel
      LOGICAL              :: EXIST_VAR
      INTEGER, ALLOCATABLE :: time(:)
      REAL*8,  ALLOCATABLE :: delta_time(:,:)
      REAL*4,  ALLOCATABLE :: lat(:,:)
      REAL*4,  ALLOCATABLE :: lon(:,:)
      INTEGER, ALLOCATABLE :: qflag(:,:)
      REAL*4,  ALLOCATABLE :: xch4(:,:)
      REAL*4,  ALLOCATABLE :: xch4_error(:,:)
      REAL*4,  ALLOCATABLE :: psurf(:,:)
      REAL*4,  ALLOCATABLE :: dp(:,:)
      REAL*4,  ALLOCATABLE :: prior(:,:,:)
      REAL*4,  ALLOCATABLE :: ak(:,:,:)

      ! For reading netCDF file
      INTEGER            :: fId              ! netCDF file ID
      INTEGER            :: st1d(1), ct1d(1) ! Start + count for 1D arrays
      INTEGER            :: st2d(2), ct2d(2) ! Start + count for 2D arrays
      INTEGER            :: st3d(3), ct3d(3) ! Start + count for 3D arrays
      INTEGER            :: st4d(4), ct4d(4) ! Start + count for 4D arrays
      CHARACTER(LEN=255) :: nc_file          ! netCDF file name
      CHARACTER(LEN=255) :: v_name           ! netCDF variable name
      CHARACTER(LEN=255) :: dir              ! Data directory path
      CHARACTER(LEN=255) :: prefix
      INTEGER                         :: I, nFiles, Count
      CHARACTER(LEN=255)              :: cmd
      CHARACTER(LEN=255), ALLOCATABLE :: Files(:) 

      !=================================================================
      ! READ_TROP_CH4_OBS begins here!
      !=================================================================

      ! Filename
      dir    = '/n/scratchlfs/jacob_lab/msulprizio/CH4/TROPOMI/' //
     &         'processed_201908/'
      prefix = 'S5P_RPRO_L2__CH4____YYYYMMDD'
      CALL Expand_Date( prefix, YYYYMMDD, 9999 ) 

      ! Determine what files to read
      CALL System('rm -rf filelist.txt')
      cmd = 'ls -1 ' // TRIM(dir) // TRIM(prefix) // '* > filelist.txt'
      CALL System(cmd)
      OPEN(31,FILE='filelist.txt',action="read")
      nFiles = 0
      DO
         READ(31,FMT='(a)',iostat=AS)
         IF (AS/=0) EXIT
         nFiles = nFiles+1
      ENDDO
      Print*, " # TROPOMI overpasses on this day = " , nFiles
      ALLOCATE(Files(nFiles))
      REWIND(31)
      DO I = 1,nFiles
         READ(31,'(a)') Files(I)
      ENDDO
      CLOSE(31)
      !CALL System('rm -rf filelist.txt')

      ! Loop over files for this day
      DO I = 1,nFiles

         IF ( I == 1 ) THEN
            ! Initialize observation counter for this day
            Count = 0
         ENDIF

         ! Initialize dimensions to zero
         nScan  = 0
         nPixel = 0
         nLev   = 0
         nTrop  = 0

         ! Construct complete file path
         nc_file = TRIM( Files(I) )
         WRITE( 6, 10 ) TRIM( nc_file )
 10      FORMAT( '     - Reading ', a)

         ! Make sure the file exists (ajt, 03/31/2013)
         INQUIRE( FILE=TRIM( nc_file ), EXIST=EXIST_VAR )
         IF ( .not. EXIST_VAR ) THEN
            Print*, 'File does not exist'
            nTrop = -1
            RETURN
         ENDIF

         ! Open netCDF file
         CALL NcOp_Rd( fId, TRIM( nc_file ) )

         ! Read the dimensions from the netCDF file
         CALL NcGet_DimLen( fId, 'layer',        nLev   )
         CALL NcGet_DimLen( fId, 'scanline',     nScan  )
         CALL NcGet_DimLen( fId, 'ground_pixel', nPixel )

         IF ( nLev > MaxLev ) THEN
            print*, '  # levels for this overpass = ', nLev
            print*, 'WARNING: nLev > MaxLev. Need to increase'
            print*, ' MaxLev in tropomi_ch4_mod.f.'
            CALL GEOS_CHEM_STOP
         ENDIF

         print*, '  # scan lines for this overpass    = ', nScan
         print*, '  # ground pixels for this overpass = ', nPixel
         print*, '  # levels for this overpass        = ', nLev

         !--------------------------------
         ! Read Data
         !--------------------------------

         ALLOCATE( time(                 1     ), STAT=AS )
         ALLOCATE( delta_time(   nPixel, nScan ), STAT=AS )
         ALLOCATE( lat(          nPixel, nScan ), STAT=AS )
         ALLOCATE( lon(          nPixel, nScan ), STAT=AS )
         ALLOCATE( qflag(        nPixel, nScan ), STAT=AS )
         ALLOCATE( xch4(         nPixel, nScan ), STAT=AS )
         ALLOCATE( xch4_error(   nPixel, nScan ), STAT=AS )
         ALLOCATE( psurf(        nPixel, nScan ), STAT=AS )
         ALLOCATE( dp(           nPixel, nScan ), STAT=AS )
         ALLOCATE( prior(  nLev, nPixel, nScan ), STAT=AS )
         ALLOCATE( ak(     nLev, nPixel, nScan ), STAT=AS )

         ! Reference time for measurements [seconds since 2010-01-01 00:00:00]
         v_name = 'time'
         st1d   = (/ 1 /)
         ct1d   = (/ 1 /)
         CALL NcRd( time, fId, TRIM(v_name), st1d, ct1d )

         ! Offset from reference start time of measurement [seconds]
         v_name = 'delta_time'
         st3d   = (/ 1,      1,     1 /)
         ct3d   = (/ nPixel, nScan, 1 /)
         CALL NcRd( delta_time, fId, TRIM(v_name), st3d, ct3d )

         ! Latitude
         v_name = 'latitude'
         st3d   = (/ 1,      1,     1 /)
         ct3d   = (/ nPixel, nScan, 1 /)
         CALL NcRd( lat, fId, TRIM(v_name), st3d, ct3d )

         ! Longitude
         v_name = 'longitude'
         st3d   = (/ 1,      1,     1 /)
         ct3d   = (/ nPixel, nScan, 1 /)
         CALL NcRd( lon, fId, TRIM(v_name), st3d, ct3d )

         ! Data quality value [0 (no data) - 1 (full quality data)]
         ! Recommend to ignore data with qa_value < 0.5
         v_name = 'qa_value'
         st3d   = (/ 1,      1,     1 /)
         ct3d   = (/ nPixel, nScan, 1 /)
         CALL NcRd( qflag, fId, TRIM(v_name), st3d, ct3d )

         ! Methane mixing ratio [ppb]
         v_name = 'methane_mixing_ratio_bias_corrected'
         st3d   = (/ 1,      1,     1 /)
         ct3d   = (/ nPixel, nScan, 1 /)
         CALL NcRd( xch4, fId, TRIM(v_name), st3d, ct3d )

         ! Methane mixing ratio precision [ppb]
         v_name = 'methane_mixing_ratio_precision'
         st3d   = (/ 1,      1,     1 /)
         ct3d   = (/ nPixel, nScan, 1 /)
         CALL NcRd( xch4_error, fId, TRIM(v_name), st3d, ct3d )

         ! Surface pressure [Pa]
         v_name = 'surface_pressure'
         st3d   = (/ 1,      1,     1 /)
         ct3d   = (/ nPixel, nScan, 1 /)
         CALL NcRd( psurf, fId, TRIM(v_name), st3d, ct3d )

         ! Pressure interval [Pa]
         v_name = 'pressure_interval'
         st3d   = (/ 1,      1,     1 /)
         ct3d   = (/ nPixel, nScan, 1 /)
         CALL NcRd( dp, fId, TRIM(v_name), st3d, ct3d )

         ! APRIORI [mol/m2]
         v_name = 'methane_profile_apriori'
         st4d   = (/ 1,    1,      1,     1 /)
         ct4d   = (/ nLev, nPixel, nScan, 1 /)
         CALL NcRd( prior, fId, TRIM(v_name), st4d, ct4d )

         ! Averaging Kernel [unitless]
         v_name = 'column_averaging_kernel'
         st4d   = (/ 1,    1,      1,     1 /)
         ct4d   = (/ nLev, nPixel, nScan, 1 /)
         CALL NcRd( ak, fId, TRIM(v_name), st4d, ct4d )

         !-------------------------------- 
         ! Place data into TROP structure
         !-------------------------------- 
         DO M = 1, nPixel
         DO N = 1, nScan

            ! Observation number on this day
            nTrop = Count + ( M * N )

            ! TIME is the same value for all observations on the same day;
            ! need to add Delta_Time to get hour, min, sec from start of day
            TROP(nTrop)%Time(1)            = time(1) + delta_time(M,N)

            TROP(nTrop)%Lev(1)             = nLev
            TROP(nTrop)%Lat(1)             = lat(M,N)
            TROP(nTrop)%Lon(1)             = lon(M,N)
            TROP(nTrop)%CH4(1)             = xch4(M,N)
            TROP(nTrop)%CH4_Error(1)       = xch4_error(M,N)
            TROP(nTrop)%QFlag(1)           = qflag(M,N)
            TROP(nTrop)%Prior(1:nLev)      = prior(1:nLev,M,N)
            TROP(nTrop)%Avg_Kernel(1:nLev) = ak(1:nLev,M,N)

            DO L = 1, nLev
               ! Convert Pa to hPa
               TROP(nTrop)%Pres(L) = (psurf(M,N) + dp(M,N)*(L-1))*1e-2
            ENDDO

         ENDDO
         ENDDO
         
         ! Update the number of observations on this day
         Count = nTrop

         ! Deallocate arrays
         IF ( ALLOCATED( time       ) ) DEALLOCATE( time       )
         IF ( ALLOCATED( delta_time ) ) DEALLOCATE( delta_time )
         IF ( ALLOCATED( lat        ) ) DEALLOCATE( lat        )
         IF ( ALLOCATED( lon        ) ) DEALLOCATE( lon        )
         IF ( ALLOCATED( qflag      ) ) DEALLOCATE( qflag      )
         IF ( ALLOCATED( xch4       ) ) DEALLOCATE( xch4       )
         IF ( ALLOCATED( xch4_error ) ) DEALLOCATE( xch4_error )
         IF ( ALLOCATED( psurf      ) ) DEALLOCATE( psurf      )
         IF ( ALLOCATED( dp         ) ) DEALLOCATE( dp         )
         IF ( ALLOCATED( prior      ) ) DEALLOCATE( prior      )
         IF ( ALLOCATED( ak         ) ) DEALLOCATE( ak         )

         !-------------------------------- 
         ! Close netCDF file
         !-------------------------------- 

         ! Echo info
         Print*, 'Number of TROPOMI observations on this day so far =',
     &            Count

         ! Close netCDF file
         CALL NcCl( fId )

      ENDDO ! End loop over files for this day

      ! Echo info
      WRITE( 6, 20 ) YYYYMMDD
 20      FORMAT('     - Finished reading TROPOMI CH4 observations for ',
     &          i8)
         Print*, 'Total number of TROPOMI observations on this day  =',
     &            Count

      END SUBROUTINE READ_TROP_CH4_OBS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_tropomi_ch4_force
!
! !DESCRIPTION: Subroutine CALC\_TROP\_CH4\_FORCE calculates the adjoint forcing
!  from the TROPOMI CH4 observations and updates the cost function.
!  (dkh, 10/12/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_TROPOMI_CH4_FORCE( Input_Opt,  State_Chm,
     &                                   State_Grid, State_Met,
     &                                   COST_FUNC )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
      USE GC_GRID_MOD,        ONLY : GET_IJ
      USE Input_Opt_Mod,      ONLY : OptInput
      USE JULDAY_MOD,         ONLY : CALDATE,     JULDAY
      USE TIME_MOD
      USE PhysConstants,      ONLY : XNUMOLAIR, AIRMW
      USE State_Chm_Mod,      ONLY : ChmState, Ind_
      USE State_Grid_Mod,     ONLY : GrdState
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
      TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
      TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      REAL(fp), INTENT(INOUT)    :: COST_FUNC   ! Cost function [unitless]
!
! !REVISION HISTORY:
!  20 Jun 2019 - M. Sulprizio- Initial version based on gosat_ch4_mod.F 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: id_CH4
      INTEGER            :: NTSTART, NTSTOP, NT
      INTEGER            :: YYYYMM,  YYYYMMDD, HHMMSS
      INTEGER            :: YEAR,    MON,      DAY
      INTEGER            :: HOUR,    MIN,      SEC
      INTEGER            :: IIJJ(2), I,      J,     N
      INTEGER            :: L,       LL,     LTROP
      INTEGER            :: JLOOP,   NOBS,   IND
      INTEGER            :: INDS(MAXOBS)
      REAL(fp)           :: REF_DATE, TIME
      REAL(fp)           :: GC_PRES(State_Grid%NZ)
      REAL(fp)           :: GC_PEDGE(State_Grid%NZ+1)
      REAL(fp)           :: GC_CH4_NATIVE(State_Grid%NZ)
      REAL(fp)           :: GC_CH4(MAXLEV)
      REAL(fp)           :: GC_CH4_cm2(MAXLEV)
      REAL(fp)           :: GC_PSURF
      REAL(fp)           :: MAP(State_Grid%NZ,MAXLEV)
      REAL(fp)           :: CH4_HAT(MAXLEV)
      REAL(fp)           :: NEW_COST(MAXOBS)
      REAL(fp)           :: OLD_COST
      REAL(fp), SAVE     :: TIME_FRAC(MAXOBS)
      INTEGER,  SAVE     :: NTROP
      REAL(fp)           :: frac, RHS, LHS
      REAL(fp)           :: CH4_PRIOR(MAXLEV)
      REAL(fp)           :: CH4_PRIOR_cm2(MAXLEV)
      REAL(fp)           :: molecontrop(MAXLEV)
      REAL(fp)           :: TROP_XCH4, GC_XCH4, TROP_XCH4_ERROR
      REAL(fp)           :: p(MAXLEV), h(MAXLEV)
      REAL(fp)           :: ak(MAXLEV), prior(MAXLEV)
      REAL(fp)           :: pres_w(MAXLEV)
      REAL(fp)           :: XCH4m, XCH4a
      REAL(fp)           :: SATELLITE_BIAS(3) ! Hardcode for now (mps,6/16/17)
      REAL(fp)           :: MEAN_MODEL_BIAS   ! Hardcode for now (mps,6/16/17)
      REAL(fp)           :: FORCE
      REAL(fp)           :: DIFF
      REAL(fp)           :: S_OBS

      ! For miscellaneous
      LOGICAL, SAVE      :: FIRST = .TRUE. 
      INTEGER            :: IOS
      INTEGER, SAVE      :: TotalObs = 0
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! CALC_TROP_CH4_FORCE begins here!
      !=================================================================

      print*, '     - CALC_TROP_CH4_FORCE '

      ! Initialize species ID flag
      id_CH4     = Ind_('CH4'       )

      ! Reset
      NEW_COST(:) = 0.0_fp

      ! Hardcode bias values for now (mps, 6/16/17)
      SATELLITE_BIAS(1) = 0.0e+0_fp
      SATELLITE_BIAS(2) = 0.0e+0_fp
      SATELLITE_BIAS(3) = 0.0e+0_fp
      MEAN_MODEL_BIAS   = 0.0e+0_fp

      ! Read Observations at first call and at end of the day
      IF ( FIRST .OR. ITS_A_NEW_DAY() ) THEN

         ! Set YYYYMMDD
         YYYYMMDD = 1d4*GET_YEAR() + 1d2*GET_MONTH() + GET_DAY()

         ! For recording model and observation values
         ! Create new file for each day
         FILENAME = 'sat_obs.tropomi.YYYYMMDD.m'
         CALL Expand_Date( FILENAME, YYYYMMDD, 9999 )
         FILENAME = TRIM( Input_Opt%RUN_DIR ) //  TRIM( FILENAME )
         OPEN( 117,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         ! Write header of sat_obs.NN.m
         WRITE( 117, 281 ) '       NNN',
     &                     '   I', '   J', '     LON','     LAT','YYYY',
     &                     'MM', 'DD', 'hh', 'mm', 'ss',
     &                     '         TAU', '       TROPOMI', 
     &                     '       model', '       S_OBS', 
     &                     '    COST_FUN'
 281     FORMAT( A10,2x,A4,2x,A4,2x,A8,2x,A8,2x,A4,2x,A2,2x,A2,2x,A2,2x,
     &           A2,2x,A2,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12)

         ! Set Total Observations = 0
         TotalObs = 0

         ! Save a value of the cost function first
         OLD_COST = COST_FUNC

         ! Read the TROPOMI CH4 file for this day
         CALL READ_TROP_CH4_OBS( YYYYMMDD, NTROP )

         IF ( FIRST ) FIRST = .FALSE.

         ! Make sure there are observations on this day
         IF ( NTROP .EQ. -1 ) RETURN

         ! TIME is seconds since 2010-01-01 00:00:00
         ! Convert to calendar date by expressing as days since REF_DATE
         ! and adding Julian day value for REF_DATE (mpayer 9/12/12)
         REF_DATE = JULDAY( 2010, 1, 1.d0 )
         DO N = 1, NTROP
            TIME = ( TROP(N)%TIME(1) / 86400.d0 ) + REF_DATE
            CALL CALDATE ( TIME, YYYYMMDD, HHMMSS )
            CALL YMD_EXTRACT( YYYYMMDD, YEAR, MON, DAY )
            CALL YMD_EXTRACT( HHMMSS,   HOUR, MIN, SEC )
            TROP(N)%YEAR(1)   = YEAR
            TROP(N)%MONTH(1)  = MON
            TROP(N)%DAY(1)    = DAY
            TROP(N)%HOUR(1)   = HOUR
            TROP(N)%MINUTE(1) = MIN
            TROP(N)%SECOND(1) = SEC
         ENDDO

      ENDIF 

      ! Get indices of TROP observations in the current hour
      !   At the start of each hour, assimilate observations that 
      !   were made in the previous 60 minutes.
      !   For example, at time 18:00, assimilate observations 
      !   made from 18:00 - 18:59
      INDS(:) = 0
      NOBS    = 0
!      print*,'Looking for observations at MONTH, DAY, HOUR = ',
!     &          GET_MONTH(), GET_DAY(), GET_HOUR()
      DO NT = 1, NTROP
         IF ( TROP(NT)%MONTH(1) .EQ. GET_MONTH() .AND.
     &        TROP(NT)%DAY(1)   .EQ. GET_DAY()   .AND.
     &        TROP(NT)%HOUR(1)  .EQ. GET_HOUR()  ) THEN
            NOBS = NOBS + 1
            INDS(NOBS) = NT
            !print*,'Found a good observation! NT = ', NT
         ENDIF
      ENDDO

      IF ( NOBS == 0 ) THEN 
         print*, ' No matching TROPOMI CH4 obs for this hour'
         RETURN
      ENDIF 

      print*, ' for hour range: ', GET_HOUR(), GET_HOUR()+1
      print*, ' found # TROPOMI observations: ', NOBS
      

!! need to update this in order to do i/o with this loop parallel 
!!      ! Now do a parallel loop for analyzing data 
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( PRIVATE )
!!!$OMP+PRIVATE( IND, NT, MAP, LTROP, IIJJ,  I, J,  L,   LL, JLOOP )
!!!$OMP+PRIVATE( GC_CH4, FORCE, CH4_PRIOR, GC_PRES, FILENAME       )
!!!$OMP+PRIVATE( GC_PEDGE, GC_PSURF, GC_CH4_NATIVE, TROP_XCH4      )
!!!$OMP+PRIVATE( TROP_XCH4_ERROR, S_OBS, h, p, XCH4a, XCH4m        )
!!!$OMP+PRIVATE( GC_XCH4, DIFF,                                    )
!!!$OMP+PRIVATE( TotalObs                                          )
      DO IND = 1, NOBS

         NT = INDS(IND)
         !print*, '     - CALC_TROP_CH4_FORCE: analyzing record ', NT 

         ! Data quality value [0 (no data) - 1 (full quality data)]
         ! Recommend to ignore data with qa_value < 0.5
         IF ( TROP(NT)%QFLAG(1) <= 0.5 ) THEN 
            ! Comment out to avoid excessive lines in log file
            !print*, ' Bad QFLAG, skipping record         ', NT
            CYCLE
         ENDIF

         ! Check for NaN in data
         IF ( IT_IS_NAN( TROP(NT)%CH4(1) ) ) THEN 
            print*, ' XCH4 is NaN, skipping record       ', NT
            CYCLE
         ENDIF
         IF ( IT_IS_NAN( TROP(NT)%PRIOR(1) ) ) THEN 
            print*, ' PRIOR is NaN, skipping record      ', NT
            CYCLE
         ENDIF
         IF ( IT_IS_NAN( TROP(NT)%AVG_KERNEL(1) ) ) THEN 
            print*, ' AVG_KERNEL is NaN, skipping record ', NT
            CYCLE
         ENDIF

         ! Check for infinity in data
         IF ( .not. IT_IS_FINITE( TROP(NT)%CH4(1) ) ) THEN
            print*, ' XCH4 is infinity, skipping record  ', NT
            CYCLE
         ENDIF

         ! Check for negative/zero data
         IF ( TROP(NT)%CH4(1) <= 0d0 ) THEN 
            print*, ' XCH4 is <= 0, skipping record      ', NT
            CYCLE
         ENDIF

         ! Check for negative/zero error values
         ! Skip zero error values for now to avoid returning S_OBS=0
         ! which will cause COST_FUN=infinity (mps, 6/6/16)
         IF ( TROP(NT)%CH4_ERROR(1) <= 0d0 ) THEN 
            print*, ' XCH4 ERROR is <= 0, skipping record ', NT
            CYCLE
         ENDIF

         ! skip Observations outside the nested domain
         ! Only skip if we're running nested (ajt, 1/26/13)
#if   defined( NESTED_NA )
         IF ( TROP(NT)%LAT(1) <   10e0   .OR. 
     &        TROP(NT)%LAT(1) >   70e0   .OR. 
     &        TROP(NT)%LON(1) < -140e0   .OR.
     &        TROP(NT)%LON(1) >  -40e0 ) THEN
            print*, ' Outside nested domain, skipping record ', NT
            CYCLE
         ENDIF
#endif

         ! Skip data at high lat (ajt, 7/23/13)
         !JDMIF ( TROP(NT)%LAT(1) >   60d0 ) THEN
         !JDMprint*, ' Skip all data above 60N ', NT
         !JDMCYCLE
         !JDMENDIF
	 !Turn of this for now to reproduce plots made by AJT
	 !Should probably turn this back on with SH once we start
         !To do the actual adjoint run

         ! Get grid box of current record
         IIJJ = GET_IJ(REAL(TROP(NT)%LON(1),4),REAL(TROP(NT)%LAT(1),4),
     &                 State_Grid)
         I    = IIJJ(1)
         J    = IIJJ(2)

         ! skip observations where the TROPOMI surface pressure is 
         ! less than the model
         IF ( (TROP(NT)%PRES(1) - State_Met%PEDGE(I,J,1)) > 50e0 ) THEN
            print*, ' Psurf threshold not met, skipping record ', NT
            CYCLE
         ENDIF

         !------------------------------
         ! Begin good observations
         !------------------------------
         print*,' Begin assimilating good observation. NT = ', NT

         ! For safety, initialize these up to LTROP
         LTROP            = 0
         GC_CH4(:)       = 0.0_fp 
         MAP(:,:)        = 0.0_fp 
         FORCE           = 0.0_fp 
         CH4_PRIOR(:)    = 0.0_fp

         ! Copy variable names to make coding a bit cleaner
         LTROP = TROP(NT)%Lev(1)
         DO L = 1, LTROP
            ! Convert [mol/m2] --> [v/v]
            CH4_PRIOR(L) = TROP(NT)%PRIOR(L) * State_Grid%Area_M2(I,J)*
     &       AIRMW / ( State_Met%AD(I,J,L) * 1000.0_fp )
         ENDDO

         ! Get GC pressure levels (mbar) 
         DO L = 1, State_Grid%NZ
            GC_PRES(L) = State_Met%PMID(I,J,L)
         ENDDO

         ! Get GC pressure edges (mbar) 
         DO L = 1, State_Grid%NZ+1
            GC_PEDGE(L) = State_Met%PEDGE(I,J,L)
         ENDDO

         ! Get GC surface pressure (mbar) 
         GC_PSURF = State_Met%PEDGE(I,J,1) 

         ! Calculate the interpolation weight matrix 
         MAP(:,:) = 0.0_fp
         CALL GET_INTMAP( State_Grid, GC_PRES, GC_PSURF, TROP(NT)%PRES,
     &                    LTROP, MAP )

         ! Get CH4 values at native model resolution
         GC_CH4_NATIVE(:) = 0.0_fp

         ! Get species concentrations
         ! Convert from [kg/box] --> [v/v]
         GC_CH4_NATIVE(:) = State_Chm%Species(I,J,:,id_CH4) * ( AIRMW / 
     &                      State_Chm%SpcData(id_CH4)%Info%emMW_g )

 
         ! Interpolate GC CH4 column to TROPOMI grid
         DO LL = 1, LTROP
            GC_CH4(LL) = 0.0_fp
            DO L = 1, State_Grid%NZ 
               GC_CH4(LL) = GC_CH4(LL) 
     &                    + MAP(L,LL) * GC_CH4_NATIVE(L) 
            ENDDO
         ENDDO

         ! TROPOMI Proxy XCH4 [ppb] --> [v/v]
         TROP_XCH4       = TROP(NT)%CH4(1)       * 1e-9_fp
         TROP_XCH4_ERROR = TROP(NT)%CH4_ERROR(1) * 1e-9_fp

         ! Remove any TROPOMI bias
         TROP_XCH4 = TROP_XCH4 + 1e-9 * ( SATELLITE_BIAS(1) 
     &                       + SATELLITE_BIAS(2)*(TROP(NT)%LAT(1))
     &                       + SATELLITE_BIAS(3)*(TROP(NT)%LAT(1))**2 )

         ! Get the S_obs, assume stddev adds in quadrature, variance
         ! adds linearly.  (ajt, 03/27/2013) 
         S_OBS = TROP_XCH4_ERROR**2 + (GC_XCH4_ERROR * 1e-9_fp)**2

         !--------------------------------------------------------------
         ! Apply TROPOMI observation operator
         !
         !   Xch4_m = Xch4_a + SUM_j( h_j * a_j * (x_m - x_a) )
         !  
         !   Xch4_m  - model XCH4
         !   Xch4_a  - apriori XCH4 = h^T * x_a
         !   h       - pressure weighting function
         !   a       - column averaging kernel
         !   x_m     - model CH4 [v/v]
         !   x_a     - apriori CH4 [v/v]
         !
         !   The pressure weighting function is defined in Connor et al. 2008 
         !     and the OCO-2 ATBD
         !--------------------------------------------------------------

         ! Pressure weighting array
         h(:)          = 0.0_fp
         p(1:LTROP)     = TROP(NT)%PRES(1:LTROP)
         ak(1:LTROP)    = TROP(NT)%AVG_KERNEL(1:LTROP)
         DO L = 1, LTROP
            ! Convert [mol/m2] --> [v/v]
            prior(L) = TROP(NT)%PRIOR(L) * State_Grid%Area_M2(I,J)*
     &       AIRMW / ( State_Met%AD(I,J,L) * 1000.0_fp )
         ENDDO

         ! Need to integrate from the toa to surface (ajt, 05/21/13)
         IF (LTROP .GT. 1) THEN
            IF(TROP(NT)%PRES(2) .LT. TROP(NT)%PRES(1)) THEN
               p(1:LTROP) = p(LTROP:1:-1)
            ENDIF
         ENDIF

         L = 1
         h(L) = 1./TROP(NT)%PRES(1) * ABS( 
     &         ( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) )
         L = LTROP
         h(L) = 1./TROP(NT)%PRES(1) * ABS( 
     &         (  p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) ) )
         DO L=2,LTROP-1
            h(L) = 1./TROP(NT)%PRES(1) * ABS( 
     &         ( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) +  
     &         (      p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) )   )
         ENDDO

         ! Now return to the orientation of the other variables
         IF (LTROP .GT. 1) THEN 
            IF(TROP(NT)%PRES(2) .LT. TROP(NT)%PRES(1)) THEN 
               h(1:LTROP) = h(LTROP:1:-1)
               p(1:LTROP) = p(LTROP:1:-1)
            ENDIF
         ENDIF

         ! Calculate Xch4_a
         XCH4a = 0.0_fp
         DO L = 1,LTROP
            XCH4a = XCH4a + h(L) * prior(L)
         ENDDO

         ! Calculate Xch4_m
         XCH4m = 0.0_fp
         XCH4m = XCH4a
         DO L = 1, LTROP
            XCH4m = XCH4m + ( h(L) * ak(L) * ( GC_CH4(L) - prior(L) ) )
         ENDDO
         GC_XCH4 = 0.0_fp
         GC_XCH4 = XCH4m

         ! Remove mean bias from GEOS-Chem XCH4
         GC_XCH4 = GC_XCH4 - ( 1e-9_fp * MEAN_MODEL_BIAS )

         !--------------------------------------------------------------
         ! Calculate cost function, given S is error in vmr
         ! J = 1/2 [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
         !--------------------------------------------------------------

         ! Calculate difference between modeled and observed profile
         DIFF = GC_XCH4 - TROP_XCH4

         ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF
         ! Need to account for the model error (ajt, 03/27/2013)
         FORCE        = ( 1.0_fp / (S_OBS) ) * DIFF
         NEW_COST(NT) = NEW_COST(NT) + 0.5e0 * DIFF * FORCE

         TotalObs = TotalObs + 1

         ! Record information for satellite diagnostics
         IF ( LDCH4SAT ) THEN 
            WRITE( 117, 282 ) TotalObs, I, J, TROP(NT)%LON(1),
     &           TROP(NT)%LAT(1),TROP(NT)%YEAR(1), 
     &           TROP(NT)%MONTH(1),TROP(NT)%DAY(1), TROP(NT)%HOUR(1), 
     &           TROP(NT)%MINUTE(1), TROP(NT)%SECOND(1),
     &           GET_TAU(), TROP_XCH4, GC_XCH4, S_OBS, 0.5d0*FORCE*DIFF
         ENDIF



      ENDDO  ! NT
!!$OMP END PARALLEL DO

      ! Number of observations processed in total
      !TotalObs = TotalObs + NOBS

      ! Update cost function 
      COST_FUNC = COST_FUNC + SUM(NEW_COST(:))

 282  FORMAT( I10,2x,I4,2x,I4,2x,F8.3,2x,F8.4,2x,I4,2x,I2,2x,I2,2x,I2,
     &        2x,I2,2x,I2,2x,F12.3,2x,E12.6,2x,E12.6,2x,E12.6,2x,E12.6)

      print*, ' Updated value of COST_FUNC = ', COST_FUNC 
      print*, ' TROP contribution           = ', COST_FUNC - OLD_COST  
      print*, ' Number of observations this hour = ', NOBS
      print*, ' Number of observations total     = ', TotalObs

      END SUBROUTINE CALC_TROPOMI_CH4_FORCE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_intmap
!
! !DESCRIPTION: Function GET\_INTMAP linearly interpolates column quatities
!   based upon the centered (average) pressue levels.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_INTMAP( State_Grid, GCPCEN, GCPSURF, TROPPEDGE,
     &                       nlev, INTMAP )
!
! !USES:
!
      USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
      TYPE(GrdState), INTENT(IN)  :: State_Grid       ! Grid State object
      REAL(fp),       INTENT(IN)  :: GCPCEN(State_Grid%NZ)
      REAL(fp),       INTENT(IN)  :: GCPSURF
      REAL(fp),       INTENT(IN)  :: TROPPEDGE(MAXLEV) 
      INTEGER,        INTENT(IN)  :: nlev
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),       INTENT(OUT) :: INTMAP(State_Grid%NZ,MAXLEV)
!
! !REVISION HISTORY:
!  16 Jun 2017 - M. Sulprizio- Initial version based on TROPOMI CH4 observation
!                              operator from GC Adjoint v35j
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: LGC, LTM
      REAL(fp) :: DIFF, DELTA_SURFP
      REAL(fp) :: LOW, HI

      !=================================================================
      ! GET_INTMAP begins here!
      !=================================================================

      ! Initialize
      INTMAP(:,:) = 0e+0_fp

      ! Loop over each pressure level of TROP grid
      DO LTM = 1, nlev

         ! Find the levels from GC that bracket level LTM
         DO LGC = 1, State_Grid%NZ-1

            LOW = GCPCEN(LGC+1)
            HI  = GCPCEN(LGC)

            ! Match GEOS-Chem level to TROP level
            IF ( TROPPEDGE(LTM) <= HI .and. 
     &           TROPPEDGE(LTM)  > LOW) THEN 

               DIFF             = HI - LOW  
               INTMAP(LGC+1,LTM) = ( HI - TROPPEDGE(LTM)  ) / DIFF
               INTMAP(LGC  ,LTM) = ( TROPPEDGE(LTM) - LOW ) / DIFF

            ENDIF

          ENDDO

       ENDDO

       ! Correct for case where TROPOMI pressure is higher than the
       ! highest GC pressure center.  In this case, just 1:1 map. 
       DO LTM = 1, nlev
          IF ( TROPPEDGE(LTM) > GCPCEN(1) ) THEN
             INTMAP(:,LTM) = 0e+0_fp
             INTMAP(1,LTM) = 1e+0_fp
          ENDIF
       ENDDO

      END SUBROUTINE GET_INTMAP
!EOC
      END MODULE TROPOMI_CH4_MOD
