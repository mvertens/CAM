module co2_data_flux

   !-------------------------------------------------------------------------------
   ! read and interpolate co2 fossil fuel surface flux
   !-------------------------------------------------------------------------------

   use shr_kind_mod,     only: r8=>shr_kind_r8, cl=>shr_kind_cl, cs=> shr_kind_cs
   use ESMF,             only: ESMF_Mesh, ESMF_Finalize, ESMF_LogFoundError
   use ESMF,             only: ESMF_END_ABORT, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
   use cam_logfile,      only: iulog
   use cam_abortutils,   only: endrun
   use spmd_utils,       only: iam, masterproc
   use dshr_strdata_mod, only: shr_strdata_type

   implicit none
   private

   ! Public interfaces
   public  :: co2_data_flux_readnl
   public  :: co2_data_flux_advance

   ! Private interfaces
   private :: co2_data_flux_init

   type :: co2_data_flux_type
      character(len=cs) :: varname = "CO2_flux"
      real(r8), pointer :: co2flx(:,:)  ! Interpolated output (pcols,begchunk:endchunk)
   end type co2_data_flux_type
   type(co2_data_flux_type), public :: data_flux_fuel

   type(shr_strdata_type)   :: sdat_co2

   character(len=cl)  :: co2flux_fuel_datafile = 'unset' ! co2 flux from fossil fuel
   character(len=cl)  :: co2flux_fuel_meshfile = 'unset' ! ESMF mesh corresponding to co2flux_fuel_datafile
   integer            :: co2flux_fuel_year_first = -999  ! first year in stream to use
   integer            :: co2flux_fuel_year_last = -999   ! last year in stream to use
   integer            :: co2flux_fuel_year_align = -999  ! align stream_year_first
   character(len=cs)  :: co2flux_fuel_tintalgo = 'unset' ! time interpolation [linear, lower, upper]
   character(len=cs)  :: co2flux_fuel_taxmode = 'unset'  ! time extraploation [cycle, extend or limit]

   logical :: debug = .false.

   character(*),parameter :: u_FILE_u = __FILE__

!===============================================================================
contains
!===============================================================================

   subroutine co2_data_flux_readnl(nlfile)

      !--------------------------------------------
      ! Purpose: Read co2_ffuel_nl namelist group.
      !--------------------------------------------
      use namelist_utils,  only: find_group_name
      use spmd_utils,      only: masterproc, mpicom, masterprocid
      use spmd_utils,      only: mpi_logical, mpi_character, mpi_integer
      use cam_logfile,     only: iulog
      use cam_abortutils,  only: endrun

      ! Arguments
      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer            :: unitn, ierr
      character(len=256) :: msg
      character(len=*), parameter :: subname = 'co2_cycle_readnl'

      namelist /co2_ffuel_nl/       &
           co2flux_fuel_datafile,   & ! input fuel dataset
           co2flux_fuel_meshfile,   & ! ESMF mesh file for input dataset
           co2flux_fuel_year_first, & ! first year in stream to use
           co2flux_fuel_year_last,  & ! last year in stream to use
           co2flux_fuel_year_align, & ! align stream_year_first
           co2flux_fuel_tintalgo,   & ! time extraploation [linear, lower, upper]
           co2flux_fuel_taxmode       ! time extraploation [cycle, extend or limit]
      !--------------------------------------------

      if (masterproc) then
         open( newunit=unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'co2_ffuel_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, co2_ffuel_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading co2_ffuel_nl namelist')
            end if
         end if
         close(unitn)
      end if

      call mpi_bcast(co2flux_fuel_datafile, len(co2flux_fuel_datafile), mpi_character, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2flux_fuel_datafile")
      call mpi_bcast(co2flux_fuel_meshfile, len(co2flux_fuel_meshfile), mpi_character, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2flux_fuel_meshfile")
      call mpi_bcast(co2flux_fuel_year_first, 1, mpi_integer, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2flux_fuel_year_first")
      call mpi_bcast(co2flux_fuel_year_last, 1, mpi_integer, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2flux_fuel_year_last")
      call mpi_bcast(co2flux_fuel_year_align, 1, mpi_integer, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2flux_fuel_year_align")
      call mpi_bcast(co2flux_fuel_tintalgo, len(co2flux_fuel_tintalgo), mpi_character, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2flux_fuel_tintalgo")
      call mpi_bcast(co2flux_fuel_taxmode, len(co2flux_fuel_taxmode), mpi_character, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2flux_fuel_taxmode")

   end subroutine co2_data_flux_readnl

   !===============================================================================
   subroutine co2_data_flux_init ()

      !--------------------------------------------
      ! Initialize co2_data_flux_type instance
      !   including initial read of input and interpolation to the current timestep
      !--------------------------------------------

      use ppgrid,           only: begchunk, endchunk, pcols
      use cam_esmf_mod,     only: model_mesh, model_clock
      use error_messages,   only: alloc_err
      use dshr_strdata_mod, only: shr_strdata_init_from_inline

      ! Local variables
      integer           :: istat
      integer           :: rc
      character(len=*), parameter :: subname = 'co2_data_flux_init'
      !--------------------------------------------

      ! Allocate data_flux_fuel%co2flx
      allocate( data_flux_fuel%co2flx(pcols,begchunk:endchunk), stat=istat)
      call alloc_err(istat, subname, 'data_flux_fuel%co2flx', pcols*(endchunk-begchunk+1))

      ! Initialize sdat_co2
      call shr_strdata_init_from_inline(sdat_co2,             &
           my_task             = iam,                         &
           logunit             = iulog,                       &
           compname            = 'ATM',                       &
           model_clock         = model_clock,                 &
           model_mesh          = model_mesh,                  &
           stream_meshfile     = trim(co2flux_fuel_meshfile), &
           stream_filenames    = (/co2flux_fuel_datafile/),   &
           stream_yearFirst    = co2flux_fuel_year_first,     &
           stream_yearLast     = co2flux_fuel_year_last,      &
           stream_yearAlign    = co2flux_fuel_year_align,     &
           stream_fldlistFile  = (/data_flux_fuel%varname/),  &
           stream_fldListModel = (/data_flux_fuel%varname/),  &
           stream_lev_dimname  = 'null',                      &
           stream_mapalgo      = 'consf',                     &
           stream_offset       = 0,                           &
           stream_taxmode      = trim(co2flux_fuel_taxmode),  &
           stream_dtlimit      = 1.0e30_r8,                   &
           stream_tintalgo     = trim(co2flux_fuel_tintalgo), &
           stream_name         = 'CO2 forcing data ',         &
           rc                  = rc)
      call chkrc(rc,__LINE__,u_FILE_u)

   end subroutine co2_data_flux_init

   !===============================================================================
   subroutine co2_data_flux_advance()

      !-------------------------------------------------------------------------------
      ! Advance the contents of a co2_data_flux_type sdat (map and interpolate in time)
      !-------------------------------------------------------------------------------

      use dshr_methods_mod , only : dshr_fldbun_getfldptr
      use dshr_strdata_mod , only : shr_strdata_advance
      use ppgrid           , only : begchunk, endchunk
      use phys_grid        , only : get_ncols_p
      use time_manager     , only : get_curr_date
      use cam_esmf_mod     , only : cam_esmf_global_sum

      ! Local variables
      integer :: icol,lchnk,g
      integer :: year    ! year (0, ...) for nstep+1
      integer :: mon     ! month (1, ..., 12) for nstep+1
      integer :: day     ! day of month (1, ..., 31) for nstep+1
      integer :: sec     ! seconds into current date for nstep+1
      integer :: mcdate  ! Current model date (yyyymmdd)
      integer :: rc
      logical :: first_time = .true.
      real(r8) :: global_sum_model, global_sum_mesh
      real(r8), pointer :: dataptr1d(:)
      character(len=*), parameter :: subname = 'co2_data_flux_advance'
      !----------------------------------------------------------------------------

      ! Initialize stream data type for fossil fuel read
      if (first_time) then
         call co2_data_flux_init()
         first_time = .false.
      end if

      ! Advance sdat stream
      call get_curr_date(year, mon, day, sec)
      mcdate = year*10000 + mon*100 + day
      call shr_strdata_advance(sdat_co2, ymd=mcdate, tod=sec, logunit=iulog, istr='co2_advance', rc=rc)
      call chkrc(rc,__LINE__,u_FILE_u)

      ! Get pointer for stream data that is time and spatially interpolated to model time and grid
      call dshr_fldbun_getFldPtr(sdat_co2%pstrm(1)%fldbun_model, data_flux_fuel%varname, fldptr1=dataptr1d, rc=rc)
      call chkrc(rc,__LINE__,u_FILE_u)

      if (debug) then
         if (masterproc) then
            write(iulog,*)
            write(iulog,'(a)')'Calling cam_esmf_global_sum from co2_data_flux'
         end if
         call cam_esmf_global_sum(trim(data_flux_fuel%varname), dataptr1d, &
              global_sum_model, global_sum_mesh, rc)
         call chkrc(rc,__LINE__,u_FILE_u)
         write(iulog,'(a)') 'Global sum for forcing field '//trim(data_flux_fuel%varname)
         write(iulog,'(a,d20.10)') ' global sum with model areas = ',global_sum_model
         write(iulog,'(a,d20.10)') ' global sum with mesh areas  = ',global_sum_mesh
      end if

      g = 1
      do lchnk = begchunk,endchunk
         do icol = 1,get_ncols_p(lchnk)
            data_flux_fuel%co2flx(icol,lchnk) = dataptr1d(g)
            g = g + 1
         end do
      end do

   end subroutine co2_data_flux_advance

   !================================================================
   subroutine chkrc(rc, line, file)
      use ESMF, only: ESMF_LOGMSG_ERROR, ESMF_SUCCESS, ESMF_LogWrite

      ! Arguments
      integer          , intent(in) :: rc
      integer          , intent(in) :: line
      character(len=*) , intent(in) :: file

      if ( rc /= ESMF_SUCCESS ) then
         call ESMF_LogWrite('ERROR:', ESMF_LOGMSG_ERROR, line=line, file=file)
         call endrun('chkrc')
      end if
   end subroutine chkrc

end module co2_data_flux
