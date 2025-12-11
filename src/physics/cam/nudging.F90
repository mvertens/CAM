module nudging
!=====================================================================
!
! Purpose: Implement Nudging of the model state of U,V,T,Q, and/or PS
!          toward specified values from analyses.
!
! Authors: Patrick Callaghan (original)
!          Mariana Vertenstein (2025) refactored for CDEPS capability
!
! Description:
!
!    This module assumes that the user has {U,V,T,Q,PS} values from analyses
!    which have been preprocessed onto the current model grid and adjusted
!    for differences in topography. It is also assumed that these resulting
!    values and are stored in individual files which are indexed with respect
!    to year, month, day, and second of the day. When the model is inbetween
!    the given begining and ending times, a relaxation forcing is added to
!    nudge the model toward the analyses values determined from the forcing
!    option specified. After the model passes the ending analyses time, the
!    forcing discontinues.
!
!    Some analyses products can have gaps in the available data, where values
!    are missing for some interval of time. When files are missing, the nudging
!    force is switched off for that interval of time, so we effectively 'coast'
!    thru the gap.
!
!    Currently, the nudging module is set up to accomodate nudging of PS
!    values, however that functionality requires forcing that is applied in
!    the selected dycore and is not yet implemented.
!
!    The nudging of the model toward the analyses data is controlled by
!    the 'nudging_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which nudging is applied, the strength of the nudging
!    tendencies, and its spatial distribution.
!
!    FORCING:
!    --------
!    Nudging tendencies are applied as a relaxation force between the current
!    model state values and target state values derived from the avalilable
!    analyses. The form of the target values is selected by the 'Nudge_Force_Opt'
!    option, the timescale of the forcing is determined from the given
!    'Nudge_TimeScale_Opt', and the nudging strength Alpha=[0.,1.] for each
!    variable is specified by the 'Nudge_Xcoef' values. Where X={U,V,T,Q,PS}
!
!           F_nudge = Alpha*((Target-Model(t_curr))/TimeScale
!
!
!    WINDOWING:
!    ----------
!    The region of applied nudging can be limited using Horizontal/Vertical
!    window functions that are constructed using a parameterization of the
!    Heaviside step function.
!
!    The Heaviside window function is the product of separate horizonal and vertical
!    windows that are controled via 12 parameters:
!
!        Nudge_Hwin_lat0:     Specify the horizontal center of the window in degrees.
!        Nudge_Hwin_lon0:     The longitude must be in the range [0,360] and the
!                             latitude should be [-90,+90].
!        Nudge_Hwin_latWidth: Specify the lat and lon widths of the window as positive
!        Nudge_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999)
!                             renders the window a constant in that direction.
!        Nudge_Hwin_latDelta: Controls the sharpness of the window transition with a
!        Nudge_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step
!                             function while a large value yeilds a smoother transition.
!        Nudge_Hwin_Invert  : A logical flag used to invert the horizontal window function
!                             to get its compliment.(e.g. to nudge outside a given window).
!
!        Nudge_Vwin_Lindex:   In the vertical, the window is specified in terms of model
!        Nudge_Vwin_Ldelta:   level indcies. The High and Low transition levels should
!        Nudge_Vwin_Hindex:   range from [0,(NLEV+1)]. The transition lengths are also
!        Nudge_Vwin_Hdelta:   specified in terms of model indices. For a window function
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NLEV+1), and the transition
!                             lengths should be set to 0.001
!        Nudge_Vwin_Invert  : A logical flag used to invert the vertical window function
!                             to get its compliment.
!
!        EXAMPLE: For a channel window function centered at the equator and independent
!                 of the vertical (30 levels):
!                        Nudge_Hwin_lat0     = 0.         Nudge_Vwin_Lindex = 0.
!                        Nudge_Hwin_latWidth = 30.        Nudge_Vwin_Ldelta = 0.001
!                        Nudge_Hwin_latDelta = 5.0        Nudge_Vwin_Hindex = 31.
!                        Nudge_Hwin_lon0     = 180.       Nudge_Vwin_Hdelta = 0.001
!                        Nudge_Hwin_lonWidth = 999.       Nudge_Vwin_Invert = .false.
!                        Nudge_Hwin_lonDelta = 1.0
!                        Nudge_Hwin_Invert   = .false.
!
!                 If on the other hand one wanted to apply nudging at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Nudge_Hwin_Invert = .true.
!
!    A user can preview the window resulting from a given set of namelist values before
!    running the model. Lookat_NudgeWindow.ncl is a script avalable in the tools directory
!    which will read in the values for a given namelist and display the resulting window.
!
!    The module is currently configured for only 1 window function. It can readily be
!    extended for multiple windows if the need arises.
!
!
! Input/Output Values:
!    Forcing contributions are available for history file output by
!    the names:    {'Nudge_U','Nudge_V','Nudge_T',and 'Nudge_Q'}
!    The target values that the model state is nudged toward are available for history
!    file output via the variables:  {'Target_U','Target_V','Target_T',and 'Target_Q'}
!
!    &nudging_nl
!      Nudge_Model         - LOGICAL toggle to activate nudging.
!                              TRUE  -> Nudging is on.
!                              FALSE -> Nudging is off.                            [DEFAULT]
!
!      Nudge_Path          - CHAR path to the analyses files.
!                              (e.g. '/glade/scratch/USER/inputdata/nudging/ERAI-Data/')
!
!      Nudge_Filenames     - CHAR array of analysis files
!
!      Nudge_Times_Per_Day - INT Number of times to update the model state (used for nudging)
!                                each day. The value is restricted to be longer than the
!                                current model timestep. As this number is increased, the nudging
!                                force has the form of newtonian cooling.
!                              48 --> 1800 Second timestep.
!                              96 -->  900 Second timestep.
!
!      Nudge_Beg_Year      - INT nudging begining year.  [1979- ]
!      Nudge_Beg_Month     - INT nudging begining month. [1-12]
!      Nudge_Beg_Day       - INT nudging begining day.   [1-31]
!
!      Nudge_End_Year      - INT nudging ending year.    [1979-]
!      Nudge_End_Month     - INT nudging ending month.   [1-12]
!      Nudge_End_Day       - INT nudging ending day.     [1-31]
!
!      Nudge_Align_Year    - INT simulation year corresponding to NUDGE_BEG_YEAR.
!                            A common usage is to set this to the first year of the model run
!                            (corresponding to the xml variable RUN_STARTDATE). With this setting,
!                            the forcing in the first year of the run will be the forcing of year
!                            yearFirst.
!                            Another usage is to align the calendar of transient forcing with the
!                            model calendar. For example, setting yearAlign = yearFirst will lead
!                            to the forcing calendar being the same as the model calendar. The
!                            forcing for a given model year would be the forcing of the same
!                            year. This would be appropriate in transient runs where the model
!                            calendar is setup to span the same year range as the forcing data.
!                            If Nudge_Align_Year is not set - then it is set to NUDGE_BEG_YEAR.
!
!      Nudge_Force_Opt     - INT Index to select the nudging Target for a relaxation forcing of the form:
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -> NEXT-OBS: Target=Anal(t'_next)    [DEFAULT]
!                              1 -> LINEAR:   Target=(F*Anal(t'_curr) +(1-F)*Anal(t'_next))
!                                                 F =(t'_next - t_curr )/Tdlt_Anal
!
!      Nudge_TimeScale_Opt - INT Index to select the timescale for nudging.
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -->  TimeScale = 1/Tdlt_Anal [DEFAULT]
!                              1 -->  TimeScale = 1/(t'_next - t_curr )
!
!      Nudge_Uprof         - INT index of profile structure to use for U.  [0,1,2]
!      Nudge_Vprof         - INT index of profile structure to use for V.  [0,1,2]
!      Nudge_Tprof         - INT index of profile structure to use for T.  [0,1,2]
!      Nudge_Qprof         - INT index of profile structure to use for Q.  [0,1,2]
!      Nudge_PSprof        - INT index of profile structure to use for PS. [0,N/A]
!
!                                The spatial distribution is specified with a profile index.
!                                 Where:  0 == OFF      (No Nudging of this variable)
!                                         1 == CONSTANT (Spatially Uniform Nudging)
!                                         2 == HEAVISIDE WINDOW FUNCTION
!
!      Nudge_Ucoef         - REAL fractional nudging coeffcient for U.
!      Nudge_Vcoef         - REAL fractional nudging coeffcient for V.
!      Nudge_Tcoef         - REAL fractional nudging coeffcient for T.
!      Nudge_Qcoef         - REAL fractional nudging coeffcient for Q.
!      Nudge_PScoef        - REAL fractional nudging coeffcient for PS.
!
!                                 The strength of the nudging is specified as a fractional
!                                 coeffcient between [0,1].
!
!      Nudge_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Nudge_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Nudge_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Nudge_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Nudge_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Nudge_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Nudge_Hwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!
!      Nudge_Vwin_Lindex   - REAL LO model index of transition
!      Nudge_Vwin_Hindex   - REAL HI model index of transition
!      Nudge_Vwin_Ldelta   - REAL LO transition length
!      Nudge_Vwin_Hdelta   - REAL HI transition length
!      Nudge_Vwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!    /
!
!================
!
! TO DO:
! -----------
!    ** Implement Ps Nudging????
!
!=====================================================================
  ! Useful modules
  !------------------
  use ESMF
  use shr_kind_mod      , only : r8=>SHR_KIND_R8, cs=>SHR_KIND_CS, cl=>SHR_KIND_CL
  use time_manager      , only : get_curr_date, get_step_size
  use cam_abortutils    , only : endrun, handle_allocate_error
  use cam_logfile       , only : iulog
  use spmd_utils        , only : masterproc, masterprocid, mpicom, mpi_success, iam
  use spmd_utils        , only : mpi_integer, mpi_real8, mpi_logical, mpi_character
  use zonal_mean_mod    , only : ZonalMean_t
  use nuopc_shr_methods , only : chkerr
  use dshr_strdata_mod  , only : shr_strdata_type
  use cam_esmf_mod      , only : model_clock, model_mesh
  use string_utils      , only : int2str

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public  :: Nudge_Model
  public  :: nudging_readnl
  public  :: nudging_init
  public  :: nudging_timestep_init
  public  :: nudging_timestep_tend
  public  :: nudging_final

  private :: nudging_set_PSprofile
  private :: nudging_set_profile
  private :: calc_DryStaticEnergy
  private :: nudging_stream_init   ! position datasets for dynamic nudging
  private :: nudging_stream_interp ! interpolates between two years of nudging file data

  integer, parameter :: maxfiles = 100

  logical, public :: Nudge_On = .false.

  ! Nudging Parameters
  !--------------------
  logical                 :: Nudge_Model       =.false.
  logical                 :: Nudge_Initialized =.false.
  character(len=cl)       :: Nudge_Meshfile
  character(len=cl)       :: Nudge_Filenames(maxfiles)
  character(len=cl)       :: Nudge_Datapath

  integer                 :: Nudge_Beg_year
  integer                 :: Nudge_Beg_month
  integer                 :: Nudge_Beg_day
  integer                 :: Nudge_Beg_sec
  type(ESMF_Time)         :: Nudge_Beg_time

  integer                 :: Nudge_End_year
  integer                 :: Nudge_End_month
  integer                 :: Nudge_End_day
  integer                 :: Nudge_End_sec
  type(ESMF_Time)         :: Nudge_End_time

  integer                 :: Nudge_Align_year
  character(len=cs)       :: Nudge_mapalgo      ! [bilinear, consf, nn]
  character(len=cs)       :: Nudge_tintalgo     ! [linear, upper]
  character(len=cs)       :: Nudge_taxmode      ! [limit, extend]
  character(len=cs)       :: Nudge_levname

  integer                 :: Model_Update_Times_Per_Day
  type(ESMF_TimeInterval) :: Model_Update_Interval
  type(ESMF_Time)         :: Model_Update_Next_Time

  integer                 :: Nudge_Force_Opt
  integer                 :: Nudge_TimeScale_Opt
  integer                 :: Nudge_TSmode

  real(r8)                :: Nudge_Ucoef,Nudge_Vcoef
  integer                 :: Nudge_Uprof,Nudge_Vprof
  real(r8)                :: Nudge_Qcoef,Nudge_Tcoef
  integer                 :: Nudge_Qprof,Nudge_Tprof
  real(r8)                :: Nudge_PScoef
  integer                 :: Nudge_PSprof

  real(r8)                :: Nudge_Hwin_lat0
  real(r8)                :: Nudge_Hwin_latWidth
  real(r8)                :: Nudge_Hwin_latDelta
  real(r8)                :: Nudge_Hwin_lon0
  real(r8)                :: Nudge_Hwin_lonWidth
  real(r8)                :: Nudge_Hwin_lonDelta
  logical                 :: Nudge_Hwin_Invert = .false.
  real(r8)                :: Nudge_Hwin_lo
  real(r8)                :: Nudge_Hwin_hi

  real(r8)                :: Nudge_Vwin_Hindex
  real(r8)                :: Nudge_Vwin_Hdelta
  real(r8)                :: Nudge_Vwin_Lindex
  real(r8)                :: Nudge_Vwin_Ldelta
  logical                 :: Nudge_Vwin_Invert =.false.
  real(r8)                :: Nudge_Vwin_lo
  real(r8)                :: Nudge_Vwin_hi

  real(r8)                :: Nudge_Hwin_latWidthH
  real(r8)                :: Nudge_Hwin_lonWidthH
  real(r8)                :: Nudge_Hwin_max
  real(r8)                :: Nudge_Hwin_min

  ! Nudging Zonal Filter variables
  !---------------------------------
  logical             :: Nudge_ZonalFilter =.false.
  integer             :: Nudge_ZonalNbasis = -1
  type(ZonalMean_t)   :: ZM
  real(r8),allocatable:: Zonal_Bamp2d(:)
  real(r8),allocatable:: Zonal_Bamp3d(:,:)

  ! Nudging State Arrays
  !-----------------------
  real(r8),allocatable:: Nudge_Utau0 (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Vtau0 (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Stau0 (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Qtau0 (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_PStau0(:,:)    !(pcols,begchunk:endchunk)

  real(r8),allocatable:: Nudge_Ustep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Vstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Sstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Qstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_PSstep(:,:)    !(pcols,begchunk:endchunk)

  ! Stream functionality
  !-----------------------
  type(shr_strdata_type) :: sdat_nudging_multi, sdat_nudging_singl
  character(len=2)       :: nudge_varlist_multi(4) = (/'U ', 'V ','T ','Q '/)
  character(len=2)       :: nudge_varlist_singl(1) = (/'PS'/)

  integer :: iunset = -999
  logical :: stream_initialized = .false.

  character(*),parameter :: u_FILE_u = __FILE__

contains

  !================================================================
  subroutine nudging_readnl(nlfile)
   !
   ! NUDGING_READNL: Initialize default values controlling the Nudging
   !                 process. Then read namelist values to override
   !                 them.
   !===============================================================
   use ppgrid,         only: pver
   use namelist_utils, only: find_group_name
   !
   ! Arguments
   !-------------
   character(len=*), intent(in) :: nlfile
   !
   ! Local Values
   !---------------
   integer :: ierr, unitn
   integer :: nfile

   character(len=*), parameter :: subname = 'nudging_readnl: '

   namelist /nudging_nl/ Nudge_Model, Nudge_datapath, Nudge_Filenames, Nudge_Meshfile, &
                         Nudge_Force_Opt, Nudge_TimeScale_Opt,                 &
                         Nudge_Beg_Year, Nudge_Beg_Month, Nudge_Beg_Day,       &
                         Nudge_End_Year, Nudge_End_Month, Nudge_End_Day,       &
                         Nudge_Align_Year, Nudge_Mapalgo, Nudge_Taxmode,       &
                         Nudge_Levname,                                        &
                         Model_Update_Times_Per_Day,                           &
                         Nudge_Ucoef , Nudge_Uprof,                            &
                         Nudge_Vcoef , Nudge_Vprof,                            &
                         Nudge_Qcoef , Nudge_Qprof,                            &
                         Nudge_Tcoef , Nudge_Tprof,                            &
                         Nudge_PScoef, Nudge_PSprof,                           &
                         Nudge_Hwin_lat0, Nudge_Hwin_lon0,                     &
                         Nudge_Hwin_latWidth, Nudge_Hwin_lonWidth,             &
                         Nudge_Hwin_latDelta, Nudge_Hwin_lonDelta,             &
                         Nudge_Hwin_Invert,                                    &
                         Nudge_Vwin_Lindex, Nudge_Vwin_Hindex,                 &
                         Nudge_Vwin_Ldelta, Nudge_Vwin_Hdelta,                 &
                         Nudge_Vwin_Invert

   ! For Zonal Mean Filtering
   namelist /nudging_nl/ Nudge_ZonalFilter, Nudge_ZonalNbasis
   ! ----------------------------------------------------------------------------

   ! Nudging is NOT initialized yet, For now
   ! Nudging will always begin/end at midnight.
   !--------------------------------------------
   Nudge_Initialized =.false.
   Nudge_End_Sec     = 0
   Nudge_Beg_Sec     = 0

   ! Set Default Namelist values
   !-----------------------------
   Nudge_Model              = .false.
   Model_Update_Times_Per_Day = 4

   Nudge_Datapath           = 'unset'
   Nudge_Filenames(:)       = 'unset'
   Nudge_Meshfile           = 'unset'

   Nudge_Beg_Year           = iunset
   Nudge_Beg_Month          = iunset
   Nudge_Beg_Day            = iunset
   Nudge_End_Year           = iunset
   Nudge_End_Month          = iunset
   Nudge_End_Day            = iunset
   Nudge_Align_Year         = iunset

   Nudge_Mapalgo            = 'bilinear'
   Nudge_taxmode            = 'limit'
   Nudge_levname            = 'lev'

   Nudge_Force_Opt          = 0
   Nudge_TimeScale_Opt      = 0
   Nudge_TSmode             = 0

   Nudge_Ucoef              = 0._r8
   Nudge_Vcoef              = 0._r8
   Nudge_Qcoef              = 0._r8
   Nudge_Tcoef              = 0._r8
   Nudge_PScoef             = 0._r8

   Nudge_Uprof              = 0
   Nudge_Vprof              = 0
   Nudge_Qprof              = 0
   Nudge_Tprof              = 0
   Nudge_PSprof             = 0

   Nudge_Hwin_lat0          = 0._r8
   Nudge_Hwin_latWidth      = 9999._r8
   Nudge_Hwin_latDelta      = 1.0_r8
   Nudge_Hwin_lon0          = 180._r8
   Nudge_Hwin_lonWidth      = 9999._r8
   Nudge_Hwin_lonDelta      = 1.0_r8
   Nudge_Hwin_Invert        = .false.

   Nudge_Vwin_Hindex        = float(pver+1)
   Nudge_Vwin_Hdelta        = 0.001_r8
   Nudge_Vwin_Lindex        = 0.0_r8
   Nudge_Vwin_Ldelta        = 0.001_r8
   Nudge_Vwin_Invert        = .false.
   Nudge_Vwin_lo            = 0.0_r8
   Nudge_Vwin_hi            = 1.0_r8

   ! Read in namelist values
   !------------------------
   if(masterproc) then
      open(newunit=unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'nudging_nl', status=ierr)
      if(ierr == 0) then
         read(unitn,nudging_nl,iostat=ierr)
         if(ierr /= 0) then
            call endrun('nudging_readnl:: ERROR reading namelist')
         end if
      end if
      close(unitn)
   end if

   ! Broadcast namelist variables
   !------------------------------
   call MPI_bcast(Nudge_Model, 1, mpi_logical, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Model')

   call MPI_bcast(Nudge_Datapath, len(Nudge_Datapath), mpi_character, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Datapath '//trim(Nudge_Datapath))
   call MPI_bcast(Nudge_Filenames(:), len(Nudge_Filenames(1))*maxfiles, mpi_character, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Filenames')
   call MPI_bcast(Nudge_Meshfile, len(Nudge_Meshfile), mpi_character, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Meshfile '//trim(Nudge_Meshfile))

   call MPI_bcast(Nudge_Levname, len(Nudge_Taxmode),  mpi_character, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Taxmode '//trim(Nudge_Taxmode))

   call MPI_bcast(Model_Update_Times_Per_Day, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Model_Update_Times_Per_Day '//&
        int2str(Model_Update_Times_Per_Day))

   call MPI_bcast(Nudge_Beg_Year, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Beg_Year '//int2str(Nudge_Beg_Year))
   call MPI_bcast(Nudge_Beg_Month, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Beg_Month'//int2str(Nudge_Beg_Month))
   call MPI_bcast(Nudge_Beg_Day, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Beg_Day '//int2str(Nudge_Beg_Day))
   call MPI_bcast(Nudge_Beg_Sec, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Beg_Sec '//int2str(Nudge_Beg_Sec))

   call MPI_bcast(Nudge_End_Year, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_End_Year '//int2str(Nudge_End_Year))
   call MPI_bcast(Nudge_End_Month, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_End_Month '//int2str(Nudge_End_Month))
   call MPI_bcast(Nudge_End_Day, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_End_Day '//int2str(Nudge_End_Day))
   call MPI_bcast(Nudge_End_Sec, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_End_Sec '//int2str(Nudge_End_Sec))

   call MPI_bcast(Nudge_Align_Year, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Align_Year '//int2str(Nudge_Align_Year))

   call MPI_bcast(Nudge_Mapalgo, len(Nudge_Mapalgo),  mpi_character, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Mapalgo '//trim(Nudge_Mapalgo))
   call MPI_bcast(Nudge_Taxmode, len(Nudge_Taxmode),  mpi_character, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Taxmode '//trim(Nudge_TaxMode))

   call MPI_bcast(Nudge_Initialized, 1, mpi_logical, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Initialized')

   call MPI_bcast(Nudge_Force_Opt, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Force_Opt '//int2str(Nudge_Force_Opt))
   call MPI_bcast(Nudge_TimeScale_Opt, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_TimeScale_Opt '//int2str(Nudge_TimeScale_Opt))
   call MPI_bcast(Nudge_TSmode, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_TSmode '//int2str(Nudge_TSmode))

   call MPI_bcast(Nudge_Ucoef, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Ucoef')
   call MPI_bcast(Nudge_Vcoef, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Vcoef')
   call MPI_bcast(Nudge_Tcoef, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Tcoef')
   call MPI_bcast(Nudge_Qcoef, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Qcoef')
   call MPI_bcast(Nudge_PScoef, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_PScoef')

   call MPI_bcast(Nudge_Uprof, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Uprof '//int2str(Nudge_Uprof))
   call MPI_bcast(Nudge_Vprof, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Vprof '//int2str(Nudge_Vprof))
   call MPI_bcast(Nudge_Tprof, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Tprof '//int2str(Nudge_Tprof))
   call MPI_bcast(Nudge_Qprof, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Qprof '//int2str(Nudge_Qprof))
   call MPI_bcast(Nudge_PSprof, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_PSprof '//int2str(Nudge_PSprof))

   call MPI_bcast(Nudge_Hwin_lat0, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Hwin_lat0')
   call MPI_bcast(Nudge_Hwin_latWidth, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Hwin_latWidth')
   call MPI_bcast(Nudge_Hwin_latDelta, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Hwin_latDelta')
   call MPI_bcast(Nudge_Hwin_lon0, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Hwin_lon0')
   call MPI_bcast(Nudge_Hwin_lonWidth, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Hwin_lonWidth')
   call MPI_bcast(Nudge_Hwin_lonDelta, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Hwin_lonDelta')
   call MPI_bcast(Nudge_Hwin_Invert,   1, mpi_logical, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Hwin_Invert')

   call MPI_bcast(Nudge_Vwin_Hindex, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Vwin_Hindex')
   call MPI_bcast(Nudge_Vwin_Hdelta, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Vwin_Hdelta')
   call MPI_bcast(Nudge_Vwin_Lindex, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Vwin_Lindex')
   call MPI_bcast(Nudge_Vwin_Ldelta, 1, mpi_real8,  masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Vwin_Ldelta')
   call MPI_bcast(Nudge_Vwin_Invert,   1, mpi_logical, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_Vwin_Invert')

   call MPI_bcast(Nudge_ZonalFilter, 1, mpi_logical, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_ZonalFilter ')
   call MPI_bcast(Nudge_ZonalNbasis, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(subname//'FATAL: mpi_bcast: Nudge_ZonalNbasis '//int2str(Nudge_ZonalNbasis))

   ! Set hi/lo values according to the given '_Invert' parameters
   !--------------------------------------------------------------
   if(Nudge_Hwin_Invert) then
      Nudge_Hwin_lo = 1.0_r8
      Nudge_Hwin_hi = 0.0_r8
   else
      Nudge_Hwin_lo = 0.0_r8
      Nudge_Hwin_hi = 1.0_r8
   end if

   if(Nudge_Vwin_Invert) then
      Nudge_Vwin_lo = 1.0_r8
      Nudge_Vwin_hi = 0.0_r8
   else
      Nudge_Vwin_lo = 0.0_r8
      Nudge_Vwin_hi = 1.0_r8
   end if

   ! Check for valid namelist values
   !----------------------------------
   if (Nudge_Beg_year == iunset) then
      call endrun('nudging_readnl:: Nudge_Beg_year '//int2str(Nudge_Beg_year)//' is not a valid value')
   end if
   if (Nudge_end_year == iunset) then
      call endrun('nudging_readnl:: Nudge_end_year '//int2str(Nudge_end_year)//' is not a valid value')
   end if
   if (Nudge_Beg_year > Nudge_end_year) then
      call endrun('nudging_readnl:: Nudge_Beg_year '//int2str(Nudge_Beg_year)//&
           'cannot be greater than Nudge_end_year '//int2str(Nudge_end_year))
   end if

   ! Determine nudge_align_year if not set
   if (Nudge_align_year == iunset) then
      Nudge_align_year = Nudge_Beg_year
   end if

   ! Determine Nudge_tintalgo
   if (Nudge_Force_Opt == 0) then
      Nudge_tintalgo = 'upper'
   elseif(Nudge_Force_Opt == 1) then
      Nudge_tintalgo = 'linear'
   else
      call endrun('nudging_timestep_init:: ERROR unknown Nudge_Force_Opt '//int2str(Nudge_Force_Opt))
   endif

   if((Nudge_Hwin_lat0 < -90._r8) .or. (Nudge_Hwin_lat0 > +90._r8)) then
     if (masterproc) then
        write(iulog,*) 'NUDGING: Window lat0 must be in [-90,+90]'
        write(iulog,*) 'NUDGING:  Nudge_Hwin_lat0=',Nudge_Hwin_lat0
     end if
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_lon0 < 0._r8) .or. (Nudge_Hwin_lon0 >= 360._r8)) then
     if (masterproc) then
        write(iulog,*) 'NUDGING: Window lon0 must be in [0,+360)'
        write(iulog,*) 'NUDGING:  Nudge_Hwin_lon0=',Nudge_Hwin_lon0
     end if
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Vwin_Lindex > Nudge_Vwin_Hindex)                          .or.   &
      (Nudge_Vwin_Hindex > float(pver+1)) .or. (Nudge_Vwin_Hindex < 0._r8) .or.  &
      (Nudge_Vwin_Lindex > float(pver+1)) .or. (Nudge_Vwin_Lindex < 0._r8)   ) then
     if (masterproc) then
        write(iulog,*) 'NUDGING: Window Lindex must be in [0,pver+1]'
        write(iulog,*) 'NUDGING: Window Hindex must be in [0,pver+1]'
        write(iulog,*) 'NUDGING: Lindex must be LE than Hindex'
        write(iulog,*) 'NUDGING:  Nudge_Vwin_Lindex=',Nudge_Vwin_Lindex
        write(iulog,*) 'NUDGING:  Nudge_Vwin_Hindex=',Nudge_Vwin_Hindex
     end if
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_latDelta <= 0._r8) .or. (Nudge_Hwin_lonDelta <= 0._r8) .or. &
      (Nudge_Vwin_Hdelta <= 0._r8) .or. (Nudge_Vwin_Ldelta <= 0._r8)    ) then
     if (masterproc) then
        write(iulog,*) 'NUDGING: Window Deltas must be positive'
        write(iulog,*) 'NUDGING:  Nudge_Hwin_latDelta=',Nudge_Hwin_latDelta
        write(iulog,*) 'NUDGING:  Nudge_Hwin_lonDelta=',Nudge_Hwin_lonDelta
        write(iulog,*) 'NUDGING:  Nudge_Vwin_Hdelta=',Nudge_Vwin_Hdelta
        write(iulog,*) 'NUDGING:  Nudge_Vwin_Ldelta=',Nudge_Vwin_Ldelta
     end if
     call endrun('nudging_readnl:: ERROR in namelist')

   endif

   if ((Nudge_Hwin_latWidth <= 0._r8) .or. (Nudge_Hwin_lonWidth <= 0._r8)) then
      if (masterproc) then
         write(iulog,*) 'NUDGING: Window widths must be positive'
         write(iulog,*) 'NUDGING:  Nudge_Hwin_latWidth=',Nudge_Hwin_latWidth
         write(iulog,*) 'NUDGING:  Nudge_Hwin_lonWidth=',Nudge_Hwin_lonWidth
      end if
      call endrun('nudging_readnl:: ERROR in namelist for Nudge_Hwin_LonWidgth')
   endif
   ! End Routine
   !------------

  end subroutine nudging_readnl
  !================================================================


  !================================================================
  subroutine nudging_init
   !
   ! NUDGING_INIT: Allocate space and initialize Nudging values
   !===============================================================
   use ppgrid        ,only: pver,pcols,begchunk,endchunk
   use error_messages,only: alloc_err
   use dycore        ,only: dycore_is
   use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
   use cam_history   ,only: addfld
   use shr_const_mod ,only: SHR_CONST_PI

   ! Local values
   !----------------
   type(ESMF_Time) :: curr_time
   type(ESMF_Time) :: Model_Update_Current_Time
   integer         :: Year,Month,Day,Sec
   logical         :: After_Beg,Before_End
   integer         :: Model_Update_Step
   integer         :: lchnk,ncol,icol,ilev
   integer         :: istat, ierr, rc
   integer         :: dtime
   real(r8)        :: rlat,rlon
   real(r8)        :: Wprof(pver)
   real(r8)        :: lonp,lon0,lonn,latp,lat0,latn
   real(r8)        :: Val1_p,Val2_p,Val3_p,Val4_p
   real(r8)        :: Val1_0,Val2_0,Val3_0,Val4_0
   real(r8)        :: Val1_n,Val2_n,Val3_n,Val4_n
   integer         :: nn, nf
   integer         :: size2d, size3d
   character(len=*), parameter :: subname = "nudging_init: "
   !----------------------------------------------------------

   ! Allocate Space
   !-----------------------------------------

   size3d = pcols*pver*((endchunk-begchunk)+1)
   size2d = pcols*((endchunk-begchunk)+1)

   allocate(Nudge_Utau0(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,subname,'Nudge_Utau',size3d)
   allocate(Nudge_Vtau0(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,subname,'Nudge_Vtau',size3d)
   allocate(Nudge_Stau0(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,subname,'Nudge_Stau',size3d)
   allocate(Nudge_Qtau0(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,subname,'Nudge_Qtau',size3d)
   allocate(Nudge_PStau0(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,subname,'Nudge_PStau',size2d)

   allocate(Nudge_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,subname,'Nudge_Ustep',size3d)
   allocate(Nudge_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,subname,'Nudge_Vstep',size3d)
   allocate(Nudge_Sstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,subname,'Nudge_Sstep',size3d)
   allocate(Nudge_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,subname,'Nudge_Qstep',size3d)
   allocate(Nudge_PSstep(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,subname,'Nudge_PSstep',size2d)

   ! Register output fields with the cam history module
   !-----------------------------------------------------
   call addfld( 'Nudge_U',(/ 'lev' /),'A','m/s/s'  ,'U Nudging Tendency')
   call addfld( 'Nudge_V',(/ 'lev' /),'A','m/s/s'  ,'V Nudging Tendency')
   call addfld( 'Nudge_T',(/ 'lev' /),'A','K/s'    ,'T Nudging Tendency')
   call addfld( 'Nudge_Q',(/ 'lev' /),'A','kg/kg/s','Q Nudging Tendency')
   call addfld('Target_U',(/ 'lev' /),'A','m/s'    ,'U Nudging Target'  )
   call addfld('Target_V',(/ 'lev' /),'A','m/s'    ,'V Nudging Target'  )
   call addfld('Target_T',(/ 'lev' /),'A','K'      ,'T Nudging Target'  )
   call addfld('Target_Q',(/ 'lev' /),'A','kg/kg'  ,'Q Nudging Target  ')

   ! Set the Stepping intervals for Model and Nudging values
   !--------------------------------------------------------

   ! Get the CAM time step size
   dtime = get_step_size()
   Model_Update_Step = 86400/Model_Update_Times_Per_Day

   if(Model_Update_Step < dtime) then
      write(iulog,*) ' '
      write(iulog,*) 'NUDGING: Model_Update_Step cannot be less than a model timestep'
      write(iulog,*) 'NUDGING:  Setting Model_Update_Step=dtime , dtime=',dtime
      write(iulog,*) ' '
      Model_Update_Step = dtime
   endif

   ! Set module time and time interval variables
   !------------------------------------------------

   call get_curr_date(Year, Month, Day, Sec)
   call ESMF_TimeSet(curr_time, yy=Year, mm=Month, dd=Day, s=Sec, rc=rc)
   call chkrc(rc,__LINE__,u_FILE_u)

   call ESMF_TimeSet(Nudge_beg_time, &
        yy=Nudge_Beg_Year, mm=Nudge_Beg_Month, dd=Nudge_Beg_Day, s=Nudge_Beg_Sec, rc=rc)
   call chkrc(rc,__LINE__,u_FILE_u)

   call ESMF_TimeSet(Nudge_end_time, &
        yy=Nudge_End_Year, mm=Nudge_End_Month, dd=Nudge_End_Day, s=Nudge_End_Sec, rc=rc)
   call chkrc(rc,__LINE__,u_FILE_u)

   call ESMF_TimeIntervalSet(Model_Update_Interval, s=Model_Update_Step, rc=rc)
   call chkrc(rc,__LINE__,u_FILE_u)

   ! Initialize the time relative to the nudging window
   !------------------------------------------------

   After_Beg  = (curr_time >= Nudge_beg_time)
   Before_End = (curr_time <= Nudge_end_time)

   if ((After_Beg) .and. (Before_End)) then

      ! Set Time indicies so that the next call to timestep_init will initialize the Model_Update_Next_time
      call ESMF_TimeSet(Model_Update_next_time, &
           yy=Year, mm=Month, dd=Day, s=(Sec/Model_Update_Step)*Model_Update_Step, rc=rc)
      call chkrc(rc,__LINE__,u_FILE_u)

   elseif (.not.After_Beg) then

      ! Set Time indicies to Nudging start so next call to timestep_init will initialize the Model_Update_Next_time
      call ESMF_TimeSet(Model_Update_next_time, &
           yy=Nudge_Beg_Year, mm=Nudge_Beg_Month, dd=Nudge_Beg_Day, s=Nudge_Beg_Sec, rc=rc)
      call chkrc(rc,__LINE__,u_FILE_u)
      ! Still need to have nudge on so that streams can be initialized - but then it will be turned off
      ! in nudging_timestep_init

   elseif (.not.Before_End) then

      ! Nudging will never occur, so switch it off
      Nudge_Model = .false.
      write(iulog,*) ' '
      write(iulog,*) 'NUDGING: WARNING - Nudging has been requested by it will'
      write(iulog,*) 'NUDGING:           never occur for the given time values'
      write(iulog,*) ' '

   endif

   ! Initialize values for window function
   !----------------------------------------
   lonp =  180._r8
   lon0 =    0._r8
   lonn = -180._r8
   latp =   90._r8-Nudge_Hwin_lat0
   lat0 =    0._r8
   latn =  -90._r8-Nudge_Hwin_lat0

   Nudge_Hwin_lonWidthH = Nudge_Hwin_lonWidth/2._r8
   Nudge_Hwin_latWidthH = Nudge_Hwin_latWidth/2._r8

   Val1_p = (1._r8+tanh((Nudge_Hwin_lonWidthH+lonp)/Nudge_Hwin_lonDelta))/2._r8
   Val2_p = (1._r8+tanh((Nudge_Hwin_lonWidthH-lonp)/Nudge_Hwin_lonDelta))/2._r8
   Val3_p = (1._r8+tanh((Nudge_Hwin_latWidthH+latp)/Nudge_Hwin_latDelta))/2._r8
   Val4_p = (1._r8+tanh((Nudge_Hwin_latWidthH-latp)/Nudge_Hwin_latDelta))/2_r8

   Val1_0 = (1._r8+tanh((Nudge_Hwin_lonWidthH+lon0)/Nudge_Hwin_lonDelta))/2._r8
   Val2_0 = (1._r8+tanh((Nudge_Hwin_lonWidthH-lon0)/Nudge_Hwin_lonDelta))/2._r8
   Val3_0 = (1._r8+tanh((Nudge_Hwin_latWidthH+lat0)/Nudge_Hwin_latDelta))/2._r8
   Val4_0 = (1._r8+tanh((Nudge_Hwin_latWidthH-lat0)/Nudge_Hwin_latDelta))/2._r8

   Val1_n = (1._r8+tanh((Nudge_Hwin_lonWidthH+lonn)/Nudge_Hwin_lonDelta))/2._r8
   Val2_n = (1._r8+tanh((Nudge_Hwin_lonWidthH-lonn)/Nudge_Hwin_lonDelta))/2._r8
   Val3_n = (1._r8+tanh((Nudge_Hwin_latWidthH+latn)/Nudge_Hwin_latDelta))/2._r8
   Val4_n = (1._r8+tanh((Nudge_Hwin_latWidthH-latn)/Nudge_Hwin_latDelta))/2._r8

   Nudge_Hwin_max =      Val1_0*Val2_0*Val3_0*Val4_0
   Nudge_Hwin_min = min((Val1_p*Val2_p*Val3_n*Val4_n), &
                        (Val1_p*Val2_p*Val3_p*Val4_p), &
                        (Val1_n*Val2_n*Val3_n*Val4_n), &
                        (Val1_n*Val2_n*Val3_p*Val4_p))

   ! Initialization is done,
   !--------------------------
   Nudge_Initialized = .true.

   if (masterproc) then

     ! Informational Output
     !---------------------------
     write(iulog,*) ' '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) '  MODEL NUDGING INITIALIZED WITH THE FOLLOWING SETTINGS: '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,'(a,l4)') 'NUDGING: Nudge_Model                = ',Nudge_Model
     write(iulog,'(2a)'  ) 'NUDGING: Nudge_Datapath             = ',trim(Nudge_Datapath)
     write(iulog,'(2a)'  ) 'NUDGING: Nudge_Meshfile             = ',trim(Nudge_Meshfile)
     write(iulog,'(2a)'  ) 'NUDGING: Nudge_Levname              = ',trim(Nudge_Levname)
     do nf = 1,maxfiles
        if (trim(Nudge_Filenames(nf)) /= ' ') then
           write(iulog,'(a,a)')'NUDGING: Nudge_Datapath             = ',len_trim(Nudge_Datapath)
        end if
     end do
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_Beg_Year             = ',Nudge_Beg_Year
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_Beg_Month            = ',Nudge_Beg_Month
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_Beg_Day              = ',Nudge_Beg_Day
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_End_Year             = ',Nudge_End_Year
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_End_Month            = ',Nudge_End_Month
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_End_Day              = ',Nudge_End_Day
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_Align_Year           = ',Nudge_Align_Year
     write(iulog,'(2a)'      ) 'NUDGING: Nudge_Mapalgo              = ',trim(Nudge_Mapalgo)
     write(iulog,'(2a)'      ) 'NUDGING: Nudge_Tintalgo             = ',trim(Nudge_Tintalgo)
     write(iulog,'(2a)'      ) 'NUDGING: Nudge_Taxmode              = ',trim(Nudge_Taxmode)
     write(iulog,'(a,i8)'    ) 'NUDGING: Model_Update_Times_Per_Day = ',Model_Update_Times_Per_Day
     write(iulog,'(a,i8)'    ) 'NUDGING: Model_Update_Step          = ',Model_Update_Step
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_PSprof               = ',Nudge_PSprof
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_Force_Opt            = ',Nudge_Force_Opt
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_TimeScale_Opt        = ',Nudge_TimeScale_Opt
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_TSmode               = ',Nudge_TSmode
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_ZonalFilter          = ',Nudge_ZonalFilter
     write(iulog,'(a,i8)'    ) 'NUDGING: Nudge_ZonalNbasis          = ',Nudge_ZonalNbasis
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Ucoef                = ',Nudge_Ucoef
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Vcoef                = ',Nudge_Vcoef
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Qcoef                = ',Nudge_Qcoef
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Tcoef                = ',Nudge_Tcoef
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_PScoef               = ',Nudge_PScoef
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Uprof                = ',Nudge_Uprof
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Vprof                = ',Nudge_Vprof
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Qprof                = ',Nudge_Qprof
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Tprof                = ',Nudge_Tprof
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_lat0            = ',Nudge_Hwin_lat0
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_latWidth        = ',Nudge_Hwin_latWidth
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_latDelta        = ',Nudge_Hwin_latDelta
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_lon0            = ',Nudge_Hwin_lon0
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_lonWidth        = ',Nudge_Hwin_lonWidth
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_lonDelta        = ',Nudge_Hwin_lonDelta
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_Invert          = ',Nudge_Hwin_Invert
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_lo              = ',Nudge_Hwin_lo
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_hi              = ',Nudge_Hwin_hi
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Vwin_Hindex          = ',Nudge_Vwin_Hindex
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Vwin_Hdelta          = ',Nudge_Vwin_Hdelta
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Vwin_Lindex          = ',Nudge_Vwin_Lindex
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Vwin_Ldelta          = ',Nudge_Vwin_Ldelta
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Vwin_Invert          = ',Nudge_Vwin_Invert
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Vwin_lo              = ',Nudge_Vwin_lo
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Vwin_hi              = ',Nudge_Vwin_hi
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_latWidthH       = ',Nudge_Hwin_latWidthH
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_lonWidthH       = ',Nudge_Hwin_lonWidthH
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_max             = ',Nudge_Hwin_max
     write(iulog,'(a,f13.5)' ) 'NUDGING: Nudge_Hwin_min             = ',Nudge_Hwin_min
     write(iulog,'(a,l4)'    ) 'NUDGING: Nudge_Initialized          = ',Nudge_Initialized
     write(iulog,*) ' '

   endif ! (masterproc) then

   ! Initialize the Zonal Mean type if needed
   !------------------------------------------
   if (Nudge_ZonalFilter) then
     call ZM%init(Nudge_ZonalNbasis)

     allocate(Zonal_Bamp2d(Nudge_ZonalNbasis),stat=istat)
     call alloc_err(istat,subname,'Zonal_Bamp2d',Nudge_ZonalNbasis)

     allocate(Zonal_Bamp3d(Nudge_ZonalNbasis,pver),stat=istat)
     call alloc_err(istat,subname,'Zonal_Bamp3d',Nudge_ZonalNbasis*pver)
   endif

   ! Initialize Nudging Coeffcient profiles in local arrays
   ! Load zeros into nudging arrays
   !------------------------------------------------------
   do lchnk = begchunk,endchunk

     ncol = get_ncols_p(lchnk)
     do icol = 1,ncol
       rlat = get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
       rlon = get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI

       call nudging_set_profile(rlat, rlon, Nudge_Uprof, Wprof, pver)
       Nudge_Utau0(icol,:,lchnk) = Wprof(:)
       call nudging_set_profile(rlat, rlon, Nudge_Vprof, Wprof, pver)
       Nudge_Vtau0(icol,:,lchnk) = Wprof(:)
       call nudging_set_profile(rlat, rlon, Nudge_Tprof, Wprof, pver)
       Nudge_Stau0(icol,:,lchnk) = Wprof(:)
       call nudging_set_profile(rlat, rlon, Nudge_Qprof, Wprof, pver)
       Nudge_Qtau0(icol,:,lchnk) = Wprof(:)
       Nudge_PStau0(icol,lchnk) = nudging_set_PSprofile(rlat, rlon, Nudge_PSprof)
     end do
   end do

   ! End Routine
   !------------

  end subroutine nudging_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_init(phys_state)
   !
   ! NUDGING_TIMESTEP_INIT:
   !   Check the current time and update Model/Nudging
   !   arrays when necessary. Toggle the Nudging flag
   !   when the time is withing the nudging window.
   !===============================================================

   use physconst    ,only: cpair
   use physics_types,only: physics_state
   use constituents ,only: cnst_get_ind
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use phys_grid    ,only: get_ncols_p
   use cam_history  ,only: outfld
   use shr_cal_mod  ,only: shr_cal_timeSet

   ! Arguments
   !-----------
   type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

   ! Local values
   !----------------
   integer                 :: Year, Month, Day, Sec
   logical                 :: Update_Model, Sync_Error
   logical                 :: Update_Nudge
   logical                 :: After_Beg, Before_End
   integer                 :: lchnk,ncol,icol,indw
   type(ESMF_Time)         :: model_curr_time
   type(ESMF_Time)         :: curr_time
   type(ESMF_TimeInterval) :: date_diff
   type(ESMF_Time)         :: time_data_LB  ! data lb time
   type(ESMF_Time)         :: time_data_ub  ! data ub  time
   type(ESMF_Time)         :: time_model    ! will have same calendar as input data
   type(ESMF_TimeInterval) :: timeint_file  ! time_data_ub - time_data_lb
   type(ESMF_TimeInterval) :: timeint_nudge ! time_data_ub - time_model
   real(r8)                :: Model_U(pcols,pver,begchunk:endchunk)
   real(r8)                :: Model_V(pcols,pver,begchunk:endchunk)
   real(r8)                :: Model_T(pcols,pver,begchunk:endchunk)
   real(r8)                :: Model_S(pcols,pver,begchunk:endchunk)
   real(r8)                :: Model_Q(pcols,pver,begchunk:endchunk)
   real(r8)                :: Model_PS(pcols,begchunk:endchunk)
   real(r8)                :: Nudge_Utau(pcols,pver,begchunk:endchunk)
   real(r8)                :: Nudge_Vtau(pcols,pver,begchunk:endchunk)
   real(r8)                :: Nudge_Stau(pcols,pver,begchunk:endchunk)
   real(r8)                :: Nudge_Qtau(pcols,pver,begchunk:endchunk)
   real(r8)                :: Nudge_PStau(pcols,begchunk:endchunk)
   real(r8)                :: Target_U(pcols,pver,begchunk:endchunk)
   real(r8)                :: Target_V(pcols,pver,begchunk:endchunk)
   real(r8)                :: Target_T(pcols,pver,begchunk:endchunk)
   real(r8)                :: Target_S(pcols,pver,begchunk:endchunk)
   real(r8)                :: Target_Q(pcols,pver,begchunk:endchunk)
   real(r8)                :: Target_PS(pcols,begchunk:endchunk)
   character(CS)           :: calendar      ! calendar name
   integer                 :: mcdate        ! current model date (yyyymmdd)
   integer                 :: Nudge_File_Step
   integer                 :: DeltaT
   real(r8)                :: Tscale
   integer                 :: rc
   logical                 :: first_call = .true.
   character(len=*), parameter :: sub = "(nudging_timestep_init) "
   !--------------------------------------------------------------

   ! Check if Nudging is initialized
   !---------------------------------
   if (.not.Nudge_Initialized) then
     call endrun('nudging_timestep_init:: Nudging NOT Initialized')
   endif

   !-------------------------------------------------------
   ! Determine if the current CAM time is AFTER the begining nudging time
   ! and if it is BEFORE the ending nudging time.
   ! Toggle Nudging flag when the time interval is between
   ! beginning and ending times, and all of the analyses files exist.
   ! When past the NEXT nudge time, update model
   !----------------------------------------------------------------

   ! Get Current CAM time
   call get_curr_date(Year,Month,Day,Sec)
   mcdate = year*10000 + month*100 + day
   call ESMF_TimeSet(curr_time, yy=Year, mm=Month, dd=Day, s=Sec, rc=rc)
   call chkrc(rc,__LINE__,u_FILE_u)

   After_Beg  = (curr_time >= Nudge_beg_time)
   Before_End = (curr_time <= Nudge_end_time)
   Nudge_On = (After_Beg .and. Before_End)
   Update_Model = (Nudge_on .and. (curr_time >= Model_Update_Next_Time))
   if (masterproc) then
      write(iulog,'(a,4(i6,2x),l8)')' Nudge Status: year, month, day, sec, update_model = ',&
           year, month, day, sec, update_model
   end if

   if (Update_Model) then

     ! Initialize nudging stream data type
     ! NOTE: this must be done once the ESMF mesh for the model is
     ! actually created - so it cannot be called out of nudging_init
     ! since that occurs before the creation of the model mesh
     !----------------------------------------------------------
     if (.not. stream_initialized) then
        call nudging_stream_init()
     end if

     ! Load values at Current into the Model arrays
     !-----------------------------------------------
     call cnst_get_ind('Q',indw)
     do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       Model_U(:ncol,:pver,lchnk) = phys_state(lchnk)%u(:ncol,:pver)
       Model_V(:ncol,:pver,lchnk) = phys_state(lchnk)%v(:ncol,:pver)
       Model_T(:ncol,:pver,lchnk) = phys_state(lchnk)%t(:ncol,:pver)
       Model_Q(:ncol,:pver,lchnk) = phys_state(lchnk)%q(:ncol,:pver,indw)
       Model_PS(:ncol,lchnk) = phys_state(lchnk)%ps(:ncol)
     end do

     ! Load Dry Static Energy values for Model
     !-----------------------------------------
     if(Nudge_TSmode == 0) then
       ! Calculate DSE from Temperature only
       do lchnk = begchunk,endchunk
         ncol = phys_state(lchnk)%ncol
         Model_S(:ncol,:pver,lchnk) = cpair*Model_T(:ncol,:pver,lchnk)
       end do
     elseif(Nudge_TSmode == 1) then
       ! Calculate DSE from Temperature, Water Vapor, and Surface Pressure
       do lchnk = begchunk,endchunk
         ncol = phys_state(lchnk)%ncol
         call calc_DryStaticEnergy(Model_T(:,:,lchnk)  , Model_Q(:,:,lchnk), &
              phys_state(lchnk)%phis,  Model_PS(:,lchnk), Model_S(:,:,lchnk), ncol)
       end do
     endif

     ! Optionally: Apply Zonal Filtering to Model state data
     !-------------------------------------------------------
     if (Nudge_ZonalFilter) then
        call ZM%calc_amps(Model_U,Zonal_Bamp3d)
        call ZM%eval_grid(Zonal_Bamp3d,Model_U)

        call ZM%calc_amps(Model_V,Zonal_Bamp3d)
        call ZM%eval_grid(Zonal_Bamp3d,Model_V)

        call ZM%calc_amps(Model_T,Zonal_Bamp3d)
        call ZM%eval_grid(Zonal_Bamp3d,Model_T)

        call ZM%calc_amps(Model_S,Zonal_Bamp3d)
        call ZM%eval_grid(Zonal_Bamp3d,Model_S)

        call ZM%calc_amps(Model_Q,Zonal_Bamp3d)
        call ZM%eval_grid(Zonal_Bamp3d,Model_Q)

        call ZM%calc_amps(Model_PS,Zonal_Bamp2d)
        call ZM%eval_grid(Zonal_Bamp2d,Model_PS)
     endif

     !-------------------------------------------------------
     ! HERE Implement time dependence of Nudging Coefs HERE
     !-------------------------------------------------------

     ! Using CDEPS:
     ! Read new nudging data and interpolate to model grid and Model_Update_Time
     !---------------------------------------------------
     call nudging_stream_interp(Target_U, Target_V, Target_T, Target_Q, Target_PS)

     do lchnk = begchunk,endchunk
        call outfld('Target_U',Target_U(:,:,lchnk),pcols,lchnk)
        call outfld('Target_V',Target_V(:,:,lchnk),pcols,lchnk)
        call outfld('Target_T',Target_T(:,:,lchnk),pcols,lchnk)
        call outfld('Target_Q',Target_Q(:,:,lchnk),pcols,lchnk)
     end do

     ! Now load Dry Static Energy values for Target
     !---------------------------------------------
     if (Nudge_TSmode == 0) then
       ! Calculate DSE from Temperature only
       do lchnk = begchunk,endchunk
         ncol = phys_state(lchnk)%ncol
         Target_S(:ncol,:pver,lchnk) = cpair*Target_T(:ncol,:pver,lchnk)
       end do
     else if(Nudge_TSmode == 1) then
        ! Calculate DSE from Temperature, Water Vapor, and Surface Pressure
        do lchnk = begchunk,endchunk
           ncol = phys_state(lchnk)%ncol
           call calc_DryStaticEnergy(Target_T(:,:,lchnk), Target_Q(:,:,lchnk), &
                phys_state(lchnk)%phis, Target_PS(:,lchnk), Target_S(:,:,lchnk), ncol)
        end do
     endif

     ! Determine Nudge_File_Step
     call get_calendar(sdat_nudging_multi, year, month, day, calendar)
     call shr_cal_timeSet(time_data_lb, &
          sdat_nudging_multi%pstrm(1)%ymdLB, sdat_nudging_multi%pstrm(1)%todLB, calendar, rc=rc)
     call chkrc(rc,__LINE__,u_FILE_u)
     call shr_cal_timeSet(time_data_ub, &
          sdat_nudging_multi%pstrm(1)%ymdUB, sdat_nudging_multi%pstrm(1)%todUB, calendar, rc=rc)
     call chkrc(rc,__LINE__,u_FILE_u)
     timeint_file = time_data_ub - time_data_lb
     call ESMF_TimeIntervalGet(timeint_file, s=Nudge_File_Step)
     call chkrc(rc,__LINE__,u_FILE_u)

     ! Determine deltaT
     call shr_cal_timeset(time_model, mcdate, sec, calendar, rc=rc)
     call chkrc(rc,__LINE__,u_FILE_u)
     timeint_nudge = time_data_ub - time_model
     call ESMF_TimeIntervalGet(timeint_nudge, s=DeltaT)
     call chkrc(rc,__LINE__,u_FILE_u)

     if (masterproc) then
        write(iulog,*)'Nudging: sdat%ymdLB, sdat%todLB ',&
             sdat_nudging_multi%pstrm(1)%ymdLB,sdat_nudging_multi%pstrm(1)%todLB
        write(iulog,*)'Nudging: sdat%ymdUB, sdat%todUB ',&
             sdat_nudging_multi%pstrm(1)%ymdUB,sdat_nudging_multi%pstrm(1)%todUB
     end if

     ! Set Tscale for the specified Forcing Option
     !-----------------------------------------------
     if(Nudge_TimeScale_Opt == 0) then
       Tscale=1._r8
     elseif (Nudge_TimeScale_Opt == 1) then
       Tscale = float(Nudge_File_Step)/float(DeltaT)
     else
       if (masterproc) then
          write(iulog,*) 'NUDGING: Unknown Nudge_TimeScale_Opt=',Nudge_TimeScale_Opt
       end if
       call endrun('nudging_timestep_init:: ERROR unknown Nudging_TimeScale_Opt')
     endif

     ! Update the nudging tendencies
     !--------------------------------
     do lchnk=begchunk,endchunk
        ncol = phys_state(lchnk)%ncol

        Nudge_Utau(:ncol,:pver,lchnk) = Nudge_Utau0(:ncol,:pver,lchnk) * Nudge_Ucoef/float(Nudge_File_Step)
        Nudge_Vtau(:ncol,:pver,lchnk) = Nudge_Vtau0(:ncol,:pver,lchnk) * Nudge_Vcoef/float(Nudge_File_Step)
        Nudge_Stau(:ncol,:pver,lchnk) = Nudge_Stau0(:ncol,:pver,lchnk) * Nudge_Tcoef/float(Nudge_File_Step)
        Nudge_Qtau(:ncol,:pver,lchnk) = Nudge_Qtau0(:ncol,:pver,lchnk) * Nudge_Qcoef/float(Nudge_File_Step)
        Nudge_PStau(:ncol,lchnk)      = Nudge_PStau0(:ncol,lchnk)      * Nudge_PScoef/float(Nudge_File_Step)

        Nudge_Ustep(:ncol,:pver,lchnk) = &
             (Target_U(:ncol,:pver,lchnk) - Model_U(:ncol,:pver,lchnk))*Tscale*Nudge_Utau(:ncol,:pver,lchnk)
        Nudge_Vstep(:ncol,:pver,lchnk) = &
             (Target_V(:ncol,:pver,lchnk) - Model_V(:ncol,:pver,lchnk))*Tscale*Nudge_Vtau(:ncol,:pver,lchnk)
        Nudge_Sstep(:ncol,:pver,lchnk) = &
             (Target_S(:ncol,:pver,lchnk) - Model_S(:ncol,:pver,lchnk))*Tscale*Nudge_Stau(:ncol,:pver,lchnk)
        Nudge_Qstep(:ncol,:pver,lchnk) = &
             (Target_Q(:ncol,:pver,lchnk) - Model_Q(:ncol,:pver,lchnk))*Tscale*Nudge_Qtau(:ncol,:pver,lchnk)
        Nudge_PSstep(:ncol,lchnk) = &
             (Target_PS(:ncol,lchnk) - Model_PS(:ncol,lchnk))*Tscale*Nudge_PStau(:ncol,lchnk)
     end do

     ! Increment the Model times by the current interval
     Model_Update_next_time = model_update_next_time + Model_Update_Interval

     ! Check for Sync Error where NEXT model time after the update
     ! is before the current time. If so, reset the next model
     ! time to a Model_Update_Step after the current time.
     Sync_Error = (curr_time >= Model_Update_next_time)
     if (Sync_Error) then
        Model_Update_next_time = curr_time + Model_Update_Interval
        write(iulog,*) 'NUDGING: WARNING - Model_Update_Time Sync ERROR... CORRECTED'
     endif

   endif ! (Update_Model)

   ! End Routine
   !------------

  end subroutine nudging_timestep_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_tend(phys_state,phys_tend)
   !
   ! NUDGING_TIMESTEP_TEND:
   !   If Nudging is ON, return the Nudging contributions
   !   to forcing using the current contents of the Nudge
   !   arrays. Send output to the cam history module as well.
   !===============================================================

   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(physics_ptend), intent(out):: phys_tend

   ! Local variables
   !--------------------
   integer :: indw,ncol,lchnk
   logical :: lq(pcnst)

   call cnst_get_ind('Q',indw)
   lq(:)    = .false.
   lq(indw) = .true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'nudging',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   if (Nudge_ON) then
      lchnk = phys_state%lchnk
      ncol  = phys_state%ncol
      Phys_tend%u(:ncol,:pver)      = Nudge_Ustep(:ncol,:pver,lchnk)
      phys_tend%v(:ncol,:pver)      = Nudge_Vstep(:ncol,:pver,lchnk)
      phys_tend%s(:ncol,:pver)      = Nudge_Sstep(:ncol,:pver,lchnk)
      phys_tend%q(:ncol,:pver,indw) = Nudge_Qstep(:ncol,:pver,lchnk)

      call outfld( 'Nudge_U',phys_tend%u          ,pcols,lchnk)
      call outfld( 'Nudge_V',phys_tend%v          ,pcols,lchnk)
      call outfld( 'Nudge_T',phys_tend%s/cpair    ,pcols,lchnk)
      call outfld( 'Nudge_Q',phys_tend%q(1,1,indw),pcols,lchnk)
   end if

   ! End Routine
   !------------

  end subroutine nudging_timestep_tend
  !================================================================


  !================================================================
  subroutine nudging_set_profile(rlat, rlon, Nudge_prof, Wprof, nlev)
   !
   ! NUDGING_SET_PROFILE:
   !   for the given lat,lon, and Nudging_prof, set the verical profile
   !   of window coeffcients.  Values range from 0. to 1. to affect
   !   spatial variations on nudging strength.
   ! ===============================================================

   ! Arguments
   !--------------
   integer  :: nlev,Nudge_prof
   real(r8) :: rlat,rlon
   real(r8) :: Wprof(nlev)

   ! Local variables
   !----------------
   integer  :: ilev
   real(r8) :: Hcoef,latx,lonx,Vmax,Vmin
   real(r8) :: lon_lo,lon_hi,lat_lo,lat_hi,lev_lo,lev_hi

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_prof == 0) then

     ! No Nudging
     Wprof(:)=0.0_r8

   elseif(Nudge_prof == 1) then

     ! Uniform Nudging
     Wprof(:)=1.0_r8

   elseif(Nudge_prof == 2) then

     ! Localized Nudging with specified Heaviside window function
     if(Nudge_Hwin_max <= Nudge_Hwin_min) then

       ! For a constant Horizontal window function,
       ! just set Hcoef to the maximum of Hlo/Hhi.
       Hcoef=max(Nudge_Hwin_lo,Nudge_Hwin_hi)

     else

       ! get lat/lon relative to window center
       latx=rlat-Nudge_Hwin_lat0
       lonx=rlon-Nudge_Hwin_lon0
       if(lonx > 180._r8) lonx=lonx-360._r8
       if(lonx <= -180._r8) lonx=lonx+360._r8

       ! Calcualte RAW window value
       lon_lo=(Nudge_Hwin_lonWidthH+lonx)/Nudge_Hwin_lonDelta
       lon_hi=(Nudge_Hwin_lonWidthH-lonx)/Nudge_Hwin_lonDelta
       lat_lo=(Nudge_Hwin_latWidthH+latx)/Nudge_Hwin_latDelta
       lat_hi=(Nudge_Hwin_latWidthH-latx)/Nudge_Hwin_latDelta
       Hcoef=((1._r8+tanh(lon_lo))/2._r8)*((1._r8+tanh(lon_hi))/2._r8) &
            *((1._r8+tanh(lat_lo))/2._r8)*((1._r8+tanh(lat_hi))/2._r8)

       ! Scale the horizontal window coef for specfied range of values.
       Hcoef=(Hcoef-Nudge_Hwin_min)/(Nudge_Hwin_max-Nudge_Hwin_min)
       Hcoef=(1._r8-Hcoef)*Nudge_Hwin_lo + Hcoef*Nudge_Hwin_hi

     endif

     ! Load the RAW vertical window
     do ilev=1,nlev
       lev_lo=(float(ilev)-Nudge_Vwin_Lindex)/Nudge_Vwin_Ldelta
       lev_hi=(Nudge_Vwin_Hindex-float(ilev))/Nudge_Vwin_Hdelta
       Wprof(ilev)=((1._r8+tanh(lev_lo))/2._r8)*((1._r8+tanh(lev_hi))/2._r8)
     end do

     ! Scale the Window function to span the values between Vlo and Vhi:
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if((Vmax <= Vmin) .or. ((Nudge_Vwin_Hindex >= (nlev+1)) .and.  &
                             (Nudge_Vwin_Lindex <= 0      )     )) then

       ! For a constant Vertical window function,
       ! load maximum of Vlo/Vhi into Wprof()
       Vmax=max(Nudge_Vwin_lo,Nudge_Vwin_hi)
       Wprof(:)=Vmax

     else

       ! Scale the RAW vertical window for specfied range of values.
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=Nudge_Vwin_lo + Wprof(:)*(Nudge_Vwin_hi-Nudge_Vwin_lo)

     endif

     ! The desired result is the product of the vertical profile
     ! and the horizontal window coeffcient.
     Wprof(:)=Hcoef*Wprof(:)

   else
     call endrun('nudging_set_profile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------

  end subroutine nudging_set_profile
  !================================================================

  !================================================================
  subroutine nudging_final

    if (allocated(Zonal_Bamp2d)) deallocate(Zonal_Bamp2d)
    if (allocated(Zonal_Bamp3d)) deallocate(Zonal_Bamp3d)

    if (allocated(Nudge_Utau0)) deallocate(Nudge_Utau0)
    if (allocated(Nudge_Vtau0)) deallocate(Nudge_Vtau0)
    if (allocated(Nudge_Stau0)) deallocate(Nudge_Stau0)
    if (allocated(Nudge_Qtau0)) deallocate(Nudge_Qtau0)
    if (allocated(Nudge_PStau0)) deallocate(Nudge_PStau0)

    if (allocated(Nudge_Ustep)) deallocate(Nudge_Ustep)
    if (allocated(Nudge_Vstep)) deallocate(Nudge_Vstep)
    if (allocated(Nudge_Sstep)) deallocate(Nudge_Sstep)
    if (allocated(Nudge_Qstep)) deallocate(Nudge_Qstep)
    if (allocated(Nudge_PSstep)) deallocate(Nudge_PSstep)

    call ZM%final()

  end subroutine nudging_final
  !================================================================

  !================================================================
  real(r8) function nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)
   !
   ! NUDGING_SET_PSPROFILE: for the given lat and lon set the surface
   !                      pressure profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) :: rlat,rlon
   integer  :: Nudge_PSprof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_PSprof == 0) then
     ! No Nudging
     !-------------
     nudging_set_PSprofile=0.0_r8
   elseif(Nudge_PSprof == 1) then
     ! Uniform Nudging
     !-----------------
     nudging_set_PSprofile=1.0_r8
   else
     call endrun('nudging_set_PSprofile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------

  end function nudging_set_PSprofile
  !================================================================


  !================================================================
  subroutine calc_DryStaticEnergy(t, q, phis, ps, dse, ncol)
   !
   ! calc_DryStaticEnergy: Given the temperature, specific humidity, surface pressure,
   !                       and surface geopotential for a chunk containing 'ncol' columns,
   !                       calculate and return the corresponding dry static energy values.
   !--------------------------------------------------------------------------------------
   use ppgrid,       only: pver, pverp
   use dycore,       only: dycore_is
   use hycoef,       only: hyai, hybi, ps0, hyam, hybm
   use physconst,    only: zvir, gravit, cpair, rair
   !
   ! Input/Output arguments
   !-----------------------
   integer , intent(in) :: ncol      ! Number of columns in chunk
   real(r8), intent(in) :: t(:,:)    ! (pcols,pver) - temperature
   real(r8), intent(in) :: q(:,:)    ! (pcols,pver) - specific humidity
   real(r8), intent(in) :: ps(:)     ! (pcols)      - surface pressure
   real(r8), intent(in) :: phis(:)   ! (pcols)      - surface geopotential
   real(r8), intent(out):: dse(:,:)  ! (pcols,pver)  - dry static energy
   !
   ! Local variables
   !------------------
   integer  :: ii,kk                 ! Lon, level, level indices
   real(r8) :: tvfac                 ! Virtual temperature factor
   real(r8) :: hkk(ncol)             ! diagonal element of hydrostatic matrix
   real(r8) :: hkl(ncol)             ! off-diagonal element
   real(r8) :: pint(ncol,pverp)      ! Interface pressures
   real(r8) :: pmid(ncol,pver )      ! Midpoint pressures
   real(r8) :: zi(ncol,pverp)        ! Height above surface at interfaces
   real(r8) :: zm(ncol,pver )        ! Geopotential height at mid level

   ! Load Pressure values and midpoint pressures
   !----------------------------------------------
   do kk=1,pverp
     do ii=1,ncol
       pint(ii,kk)=(hyai(kk)*ps0)+(hybi(kk)*ps(ii))
     end do
   end do
   do kk=1,pver
     do ii=1,ncol
       pmid(ii,kk)=(hyam(kk)*ps0)+(hybm(kk)*ps(ii))
     end do
   end do

   ! The surface height is zero by definition.
   !-------------------------------------------
   do ii = 1,ncol
     zi(ii,pverp) = 0.0_r8
   end do

   ! Compute the dry static energy, zi, zm from bottom up
   ! Note, zi(i,k) is the interface above zm(i,k)
   !---------------------------------------------------------
   do kk=pver,1,-1

     ! First set hydrostatic elements consistent with dynamics
     !--------------------------------------------------------
     if (dycore_is ('LR')) then
       do ii=1,ncol
         hkl(ii)=log(pint(ii,kk+1))-log(pint(ii,kk))
         hkk(ii)=1._r8-(hkl(ii)*pint(ii,kk)/(pint(ii,kk+1)-pint(ii,kk)))
       end do
     else
       do ii=1,ncol
         hkl(ii)=(pint(ii,kk+1)-pint(ii,kk))/pmid(ii,kk)
         hkk(ii)=0.5_r8*hkl(ii)
       end do
     endif

     ! Now compute zm, zi, and dse  (WACCM-X vars rairv/zairv/cpairv not used!)
     !------------------------------------------------------------------------
     do ii=1,ncol
       tvfac=t(ii,kk)*rair*(1._r8+(zvir*q(ii,kk)))/gravit
       zm (ii,kk)=zi(ii,kk+1) + (tvfac*hkk(ii))
       zi (ii,kk)=zi(ii,kk+1) + (tvfac*hkl(ii))
       dse(ii,kk)=(t(ii,kk)*cpair) + phis(ii) + (gravit*zm(ii,kk))
     end do

   end do ! kk=pver,1,-1

   ! End Routine
   !-----------

  end subroutine calc_DryStaticEnergy
  !================================================================


  !================================================================
  subroutine nudging_stream_init()

    use dshr_strdata_mod, only: shr_strdata_init_from_inline

    ! local variables
    integer                 :: rc
    integer                 :: nfile
    integer                 :: nudge_year_first
    integer                 :: nudge_year_last
    integer                 :: nudge_year_align
    character(*), parameter :: sub = "('nudging_stream_init')"
    !----------------------------------------------------------------

    ! Write output log info

    call ESMF_TimeGet(nudge_beg_time, yy=nudge_year_first, rc=rc)
    call chkrc(rc,__LINE__,u_FILE_u)
    call ESMF_TimeGet(nudge_end_time, yy=nudge_year_last, rc=rc)
    call chkrc(rc,__LINE__,u_FILE_u)
    nudge_year_align = nudge_align_year

    if (masterproc) then
       write(iulog,'(a)'   ) ' '
       write(iulog,'(a,i8)')  'stream nudging settings:'
       write(iulog,'(a,a,a)') '  nudge varlist    = ','U,V,T,Q,PS'
       write(iulog,'(a,i8)')  '  nudge year first = ',nudge_year_first
       write(iulog,'(a,i8)')  '  nudge year last  = ',nudge_year_last
       write(iulog,'(a,i8)')  '  nudge year align = ',nudge_year_align
       write(iulog,'(a,a)')   '  nudge tintalgo   = ',trim(nudge_tintalgo)
       write(iulog,'(a,a)' )  '  nudge meshfile   = ',trim(nudge_meshfile)
       write(iulog,'(a,a)' )  '  nudge datapath   = ',trim(nudge_datapath)
       do nfile = 1,size(nudge_filenames)
          if (trim(nudge_filenames(nfile)) /= ' ') then
             write(iulog,'(a,i8,2x,a)' )  '  nudge files = ',nfile,trim(nudge_filenames(nfile))
          end if
       end do
       write(iulog,'(a)'   )  ' '
    endif

    do nfile = 1,size(nudge_filenames)
       if (trim(nudge_filenames(nfile)) /= ' ') then
          nudge_filenames(nfile) = trim(nudge_datapath)//'/'//trim(nudge_filenames(nfile))
       end if
    end do

    ! Create module stream data type sdat_nudging

    call shr_strdata_init_from_inline(sdat_nudging_multi, &
         my_task             = iam,                       &
         logunit             = iulog,                     &
         compname            = 'ATM',                     &
         model_clock         = model_clock,               &
         model_mesh          = model_mesh,                &
         stream_meshfile     = trim(nudge_meshfile),      &
         stream_filenames    = nudge_filenames,           &
         stream_yearFirst    = nudge_year_first,          &
         stream_yearLast     = nudge_year_last,           &
         stream_yearAlign    = nudge_year_align,          &
         stream_fldlistFile  = nudge_varlist_multi,       &
         stream_fldListModel = nudge_varlist_multi,       &
         stream_lev_dimname  = trim(nudge_levname),       &
         stream_mapalgo      = trim(nudge_mapalgo),       &
         stream_offset       = 0,                         &
         stream_taxmode      = trim(nudge_taxmode),       &
         stream_dtlimit      = 1.0e30_r8,                 &
         stream_tintalgo     = nudge_tintalgo,            &
         stream_name         = 'NUDGING forcing data ',   &
         rc                  = rc)
    call chkrc(rc,__LINE__,u_FILE_u)

    call shr_strdata_init_from_inline(sdat_nudging_singl, &
         my_task             = iam,                       &
         logunit             = iulog,                     &
         compname            = 'ATM',                     &
         model_clock         = model_clock,               &
         model_mesh          = model_mesh,                &
         stream_meshfile     = trim(nudge_meshfile),      &
         stream_filenames    = nudge_filenames,           &
         stream_yearFirst    = nudge_year_first,          &
         stream_yearLast     = nudge_year_last,           &
         stream_yearAlign    = nudge_year_first,          &
         stream_fldlistFile  = nudge_varlist_singl,       &
         stream_fldListModel = nudge_varlist_singl,       &
         stream_lev_dimname  = 'null',                    &
         stream_mapalgo      = trim(nudge_mapalgo),       &
         stream_offset       = 0,                         &
         stream_taxmode      = trim(nudge_taxmode),       &
         stream_dtlimit      = 1.0e30_r8,                 &
         stream_tintalgo     = nudge_tintalgo,            &
         stream_name         = 'NUDGING forcing data ',   &
         rc                  = rc)
    call chkrc(rc,__LINE__,u_FILE_u)

    stream_initialized = .true.

  end subroutine nudging_stream_init
  !================================================================


  !================================================================
  subroutine nudging_stream_interp(Target_U, Target_V, Target_T, Target_Q, Target_PS)

    ! Caculate Target_T, Target_U, Target_V, Target_Q and Target_PS

    use dshr_strdata_mod , only : shr_strdata_advance
    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    use ppgrid           , only : pcols, pver, begchunk, endchunk
    use phys_grid        , only : get_ncols_p

    ! Arguments
    real(r8), intent(out) :: Target_U(pcols,pver,begchunk:endchunk)
    real(r8), intent(out) :: Target_V(pcols,pver,begchunk:endchunk)
    real(r8), intent(out) :: Target_T(pcols,pver,begchunk:endchunk)
    real(r8), intent(out) :: Target_Q(pcols,pver,begchunk:endchunk)
    real(r8), intent(out) :: Target_PS(pcols,begchunk:endchunk)

    ! Local variables
    integer :: rc     ! ESMF error return
    integer :: istat  ! allocate return
    integer :: nvar   ! variable index
    integer :: klev   ! level index
    integer :: icol   ! column index
    integer :: ncol   ! number of columns in chunk
    integer :: lchnk  ! chunk index
    integer :: g      ! counter index
    integer :: year   ! year (0, ...) for nstep+1
    integer :: mon    ! month (1, ..., 12) for nstep+1
    integer :: day    ! day of month (1, ..., 31) for nstep+1
    integer :: sec    ! seconds into current date for nstep+1
    integer :: mcdate ! current model date (yyyymmdd)
    real(r8), pointer :: dataptr2d(:,:) ! first dimension is level, second is data on that level
    real(r8), pointer :: dataptr1d(:)
    real(r8)          :: Tmp3D(pcols,pver,begchunk:endchunk)
    real(r8)          :: Tmp2D(pcols,begchunk:endchunk)
    character(len=*), parameter :: sub = "(nudging_stream_interp) "
    !-----------------------------------------------------------------------

    ! Extract YMD from model_update_next_time
    call ESMF_TimeGet(Model_Update_Next_Time, yy=year, mm=mon, dd=day, s=sec, rc=rc)
    call chkrc(rc,__LINE__,u_FILE_u)
    mcdate = year*10000 + mon*100 + day
    if (masterproc) then
       write(iulog,'(a,4(i6,2x))')' nudging_stream_interp: interpolating nudge to ',year,mon,day,sec
    end if

    ! Advance sdat streams
    call shr_strdata_advance(sdat_nudging_multi, ymd=mcdate, tod=sec, logunit=iulog, istr='nudging', rc=rc)
    call chkrc(rc,__LINE__,u_FILE_u)

    call shr_strdata_advance(sdat_nudging_singl, ymd=mcdate, tod=sec, logunit=iulog, istr='nudging', rc=rc)
    call chkrc(rc,__LINE__,u_FILE_u)

    ! Get pointer for stream data that is time and spatially interpolated to model time and grid
    ! Obtain Target_U, Target_V, Target_T and Target_Q

    do nvar = 1,4
       if ( trim(nudge_varlist_multi(nvar)) == 'U' .or. &
            trim(nudge_varlist_multi(nvar)) == 'V' .or. &
            trim(nudge_varlist_multi(nvar)) == 'T' .or. &
            trim(nudge_varlist_multi(nvar)) == 'Q' )  then

          call dshr_fldbun_getFldPtr(sdat_nudging_multi%pstrm(1)%fldbun_model, &
               nudge_varlist_multi(nvar), fldptr2=dataptr2d, rc=rc)
          call chkrc(rc,__LINE__,u_FILE_u)

          ! Obtain TMP3d
          do klev = 1, pver
             g = 1
             do lchnk = begchunk,endchunk
                ncol = get_ncols_p(lchnk)
                do icol = 1,ncol
                   Tmp3d(icol,klev,lchnk) = dataptr2d(klev,g)
                   g = g + 1
                end do
             end do
          end do

          ! Apply zonal mean filtering
          if (Nudge_ZonalFilter) then
             call ZM%calc_amps(Tmp3D, Zonal_Bamp3d)
             call ZM%eval_grid(Zonal_Bamp3d, Tmp3D)
          endif

          ! Determine output variables
          if (trim(nudge_varlist_multi(nvar)) == 'U') then
             do lchnk = begchunk,endchunk
                ncol = get_ncols_p(lchnk)
                Target_U(:ncol,:pver,lchnk) = Tmp3d(:ncol,:pver,lchnk)
             end do
          else if (trim(nudge_varlist_multi(nvar)) == 'V') then
             do lchnk = begchunk,endchunk
                ncol = get_ncols_p(lchnk)
                Target_V(:ncol,:pver,lchnk) = Tmp3d(:ncol,:pver,lchnk)
             end do
          else if (trim(nudge_varlist_multi(nvar)) == 'T') then
             do lchnk = begchunk,endchunk
                ncol = get_ncols_p(lchnk)
                Target_T(:ncol,:pver,lchnk) = Tmp3d(:ncol,:pver,lchnk)
             end do
          else if (trim(nudge_varlist_multi(nvar)) == 'Q') then
             do lchnk = begchunk,endchunk
                ncol = get_ncols_p(lchnk)
                Target_Q(:ncol,:pver,lchnk) = Tmp3d(:ncol,:pver,lchnk)
             end do
          end if
       end if
    end do

    ! Obtain Target_PS

    call dshr_fldbun_getFldPtr(sdat_nudging_singl%pstrm(1)%fldbun_model, 'PS', fldptr1=dataptr1d, rc=rc)
    call chkrc(rc,__LINE__,u_FILE_u)

    g = 1
    do lchnk = begchunk,endchunk
       ncol = get_ncols_p(lchnk)
       do icol = 1,ncol
          Tmp2d(icol,lchnk) = dataptr1d(g)
          g = g + 1
       end do
    end do

    if (Nudge_ZonalFilter) then
       call ZM%calc_amps(Tmp2D,Zonal_Bamp2d)
       call ZM%eval_grid(Zonal_Bamp2d,Tmp2D)
    endif

    do lchnk=begchunk,endchunk
       ncol = get_ncols_p(lchnk)
       Target_PS(:ncol,lchnk)= Tmp2d(:ncol,lchnk)
    end do

  end subroutine nudging_stream_interp
  !================================================================


  !================================================================
  subroutine get_calendar(sdat, model_year, model_month, model_day, calendar)

     use shr_cal_mod, only : shr_cal_noleap, shr_cal_gregorian
     use shr_cal_mod, only : shr_cal_date2ymd, shr_cal_leapyear

     ! Arguments
     type(shr_strdata_type), intent(in)  :: sdat
     integer,                intent(in)  :: model_year  ! model year
     integer,                intent(in)  :: model_month ! model month
     integer,                intent(in)  :: model_day   ! model day
     character(len=*),       intent(out) :: calendar

     ! Local Variables
     integer :: data_year, data_month, data_day ! data date year month day
     !-----------------------------------------------------------------------

     call shr_cal_date2ymd(sdat%pstrm(1)%ymdUB, data_year, data_month, data_day)

     calendar = trim(sdat%stream(1)%calendar)
     if (trim(sdat%model_calendar) /= trim(sdat%stream(1)%calendar)) then
        if (( trim(sdat%model_calendar) == trim(shr_cal_gregorian)) .and. &
             (trim(sdat%stream(1)%calendar) == trim(shr_cal_noleap))) then
           ! set feb 29 = feb 28
           if (model_month == 2 .and. model_day == 29) then
              calendar = shr_cal_noleap
           end if
        elseif ((trim(sdat%model_calendar) == trim(shr_cal_noleap)) .and. &
             (trim(sdat%stream(1)%calendar) == trim(shr_cal_gregorian))) then
           ! feb 29 input data will be skipped automatically
           if (data_month==3 .and. data_day==1 .and. model_month==2 .and. model_day==28) then
              calendar = shr_cal_noleap
           endif
        endif
     else ! calendars are the same
        if (trim(sdat%model_calendar) == trim(shr_cal_gregorian)) then
           ! Both are in gregorian - but it's possible that there is a mismatch
           ! such that the model is in leapyear but the data is not
           if (model_month == 2 .and. model_day >= 28) then
              if (shr_cal_leapyear(model_year) .and. .not. shr_cal_leapyear(data_year)) then
                 ! model is in leap year but data is not
                 calendar = shr_cal_noleap
              endif
           else
              calendar = sdat%model_calendar
           endif
        else
           calendar = sdat%model_calendar
        endif
     endif

  end subroutine get_calendar
  !================================================================


  !================================================================
  subroutine chkrc(rc, line, file)
     integer          , intent(in) :: rc
     integer          , intent(in) :: line
     character(len=*) , intent(in) :: file

     if ( rc /= ESMF_SUCCESS ) then
        call ESMF_LogWrite('ERROR:', ESMF_LOGMSG_ERROR, line=line, file=file)
        call endrun('chkrc')
     end if
  end subroutine chkrc

end module nudging
