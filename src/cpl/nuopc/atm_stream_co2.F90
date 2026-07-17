module atm_stream_co2

  !-----------------------------------------------------------------------
  ! Contains methods for reading in co2_surface_source deposition data file
  ! Also includes functions for dynamic co2_surface_source file handling and
  ! interpolation.
  !-----------------------------------------------------------------------
  !
  use ESMF              , only : ESMF_SUCCESS 
  use nuopc_shr_methods , only : chkerr
  use dshr_strdata_mod  , only : shr_strdata_type
  use shr_kind_mod      , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use spmd_utils        , only : mpicom, masterproc, iam, masterprocid
  use spmd_utils        , only : mpi_character, mpi_integer, mpi_logical
  use cam_logfile       , only : iulog
  use cam_abortutils    , only : endrun
  use cam_esmf_mod      , only : model_clock, model_mesh

  implicit none
  private

  public :: stream_co2_surface_source_readnl    ! read runtime options
  public :: stream_co2_surface_source_init      ! position datasets for dynamic co2_surface_source
  public :: stream_co2_surface_source_interp    ! interpolates between two years of co2_surface_source file data

  type(shr_strdata_type) :: sdat_co2_surface_source     ! input data stream

  ! namelist variables
  character(len=CL) :: stream_co2_surface_source_mesh_filename
  character(len=CL) :: stream_co2_surface_source_data_filename
  character(len=CL) :: stream_co2_surface_source_data_varname ! variable name for co2_surface_source on stream file(s)
  character(len=CS) :: stream_co2_surface_source_taxmode      ! 'cycle' or 'extend' or 'limit'
  integer           :: stream_co2_surface_source_year_first   ! first year in stream to use
  integer           :: stream_co2_surface_source_year_last    ! last year in stream to use
  integer           :: stream_co2_surface_source_year_align   ! align stream_year_first

  logical, public, protected :: co2_surface_source
  logical, public, protected :: stream_co2_surface_source_is_initialized = .false.

  character(len=*), parameter :: u_FILE_u = __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine stream_co2_surface_source_readnl(nlfile)

    ! Uses:
    use shr_nl_mod,     only: shr_nl_find_group_name

    ! input/output variables
    character(len=*), intent(in) :: nlfile

    ! local variables
    integer :: nu_nml    ! unit for namelist file
    integer :: nml_error ! namelist i/o error flag
    integer :: ierr      ! error status
    integer :: nf        ! field counter
    character(*), parameter :: subName = "('stream_co2_surface_source_readnl')"
    !-----------------------------------------------------------------------

    namelist /co2_surface_source_stream_nl/       &
         co2_surface_source,                      &
         stream_co2_surface_source_mesh_filename, &
         stream_co2_surface_source_data_filename, &
         stream_co2_surface_source_data_varname,  &
         stream_co2_surface_source_taxmode,       &
         stream_co2_surface_source_year_first,    &
         stream_co2_surface_source_year_last,     &
         stream_co2_surface_source_year_align

    ! Default values for namelist
    co2_surface_source = .false.
    stream_co2_surface_source_data_filename = ' '
    stream_co2_surface_source_mesh_filename = ' '
    stream_co2_surface_source_data_varname  = ' '
    stream_co2_surface_source_taxmode       = 'unset'
    stream_co2_surface_source_year_first    = -999 ! first year in stream to use
    stream_co2_surface_source_year_last     = -999 ! last  year in stream to use
    stream_co2_surface_source_year_align    = -999 ! align stream_co2_surface_source_year_first with this model year

    ! Read co2_surface_source_stream namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(nlfile), status='old', iostat=nml_error )
       if (nml_error /= 0) then
          call endrun(subName//': ERROR opening '//trim(nlfile)//errMsg(u_FILE_u, __LINE__))
       end if
       call shr_nl_find_group_name(nu_nml, 'co2_surface_source_stream_nl', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=co2_surface_source_stream_nl, iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(' ERROR reading co2_surface_source_stream_nl namelist'//errMsg(u_FILE_u, __LINE__))
          end if
       end if
       close(nu_nml)

       ! Error check
       if ( trim(stream_co2_surface_source_taxmode) /= 'cycle'  .and. &
            trim(stream_co2_surface_source_taxmode) /= 'extend' .and. &
            trim(stream_co2_surface_source_taxmode) /= 'limit') then
          call endrun(subName//': ERROR stream_co2_surface_source_taxmode '&
               //trim(stream_co2_surface_source_taxmode)&
               //' must be either cycle, extend or limit')
       end if
    endif

    call mpi_bcast(co2_surface_source, 1, mpi_logical, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: co2_surface_source")

    if (co2_surface_source) then
       call mpi_bcast(stream_co2_surface_source_mesh_filename, &
            len(stream_co2_surface_source_mesh_filename), mpi_character, masterprocid, mpicom, ierr)
       if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_co2_surface_source_mesh_filename")
       call mpi_bcast(stream_co2_surface_source_data_filename, &
            len(stream_co2_surface_source_data_filename), mpi_character, masterprocid, mpicom, ierr)
       if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_co2_surface_source_data_filename")
       call mpi_bcast(stream_co2_surface_source_data_varname, &
            len(stream_co2_surface_source_data_varname), mpi_character, masterprocid, mpicom, ierr)
       if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_co2_surface_source_data_varname")
       call mpi_bcast(stream_co2_surface_source_taxmode, &
            len(stream_co2_surface_source_taxmode), mpi_character, masterprocid, mpicom, ierr)
       if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_co2_surface_source_taxmode")
       call mpi_bcast(stream_co2_surface_source_year_first, &
            1, mpi_integer, masterprocid, mpicom, ierr)
       if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_co2_surface_source_year_first")
       call mpi_bcast(stream_co2_surface_source_year_last, &
            1, mpi_integer, masterprocid, mpicom, ierr)
       if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_co2_surface_source_year_last")
       call mpi_bcast(stream_co2_surface_source_year_align, &
            1, mpi_integer, masterprocid, mpicom, ierr)
    end if

    if (masterproc) then
       write(iulog,'(a)')  ' '
       if (co2_surface_source) then
          write(iulog,'(2a)')    subname,' co2 surface source override settings::'
          write(iulog,'(3a)')    subname,'  stream_co2_surface_source_data_filename = ',&
               trim(stream_co2_surface_source_data_filename)
          write(iulog,'(3a)')    subname,'  stream_co2_surface_source_mesh_filename = ',&
               trim(stream_co2_surface_source_mesh_filename)
          write(iulog,'(3a)')    subname,'  stream_co2_surface_source_data_varname  = ',&
               trim(stream_co2_surface_source_data_varname)
          write(iulog,'(3a)')    subname,'  stream_co2_surface_source_taxmode       = ',&
               trim(stream_co2_surface_source_taxmode)
          write(iulog,'(2a,i0)') subname,'  stream_co2_surface_source_year_first    = ',&
               stream_co2_surface_source_year_first
          write(iulog,'(2a,i0)') subname,'  stream_co2_surface_source_year_last     = ',&
               stream_co2_surface_source_year_last
          write(iulog,'(2a,i0)') subname,'  stream_co2_surface_source_year_align    = ',&
               stream_co2_surface_source_year_align
          write(iulog,'(a)') ' '
       else
          write(iulog, '(2a)') subname, 'co2 surface source will not be overwritten'
       end if
    endif

  end subroutine stream_co2_surface_source_readnl

  !================================================================
  subroutine stream_co2_surface_source_init(rc)
    use dshr_strdata_mod, only: shr_strdata_init_from_inline

    ! input/output variables
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: subName = "('stream_co2_surface_source_init')"
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Initialize the cdeps data type sdat_co2_surface_source
    call shr_strdata_init_from_inline(sdat_co2_surface_source,                    &
         my_task             = iam,                                               &
         logunit             = iulog,                                             &
         compname            = 'ATM',                                             &
         model_clock         = model_clock,                                       &
         model_mesh          = model_mesh,                                        &
         stream_meshfile     = trim(stream_co2_surface_source_mesh_filename),     &
         stream_filenames    = (/trim(stream_co2_surface_source_data_filename)/), &
         stream_yearFirst    = stream_co2_surface_source_year_first,              &
         stream_yearLast     = stream_co2_surface_source_year_last,               &
         stream_yearAlign    = stream_co2_surface_source_year_align,              &
         stream_fldlistFile  = (/stream_co2_surface_source_data_varname/),        &
         stream_fldListModel = (/stream_co2_surface_source_data_varname/),        &
         stream_lev_dimname  = 'null',                                            &
         stream_mapalgo      = 'bilinear',                                        &
         stream_offset       = 0,                                                 &
         stream_taxmode      = trim(stream_co2_surface_source_taxmode),           &
         stream_dtlimit      = 1.0e30_r8,                                         &
         stream_tintalgo     = 'linear',                                          &
         stream_name         = 'CO2_SURFACE_SOURCE data ',                        &
         rc                  = rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    stream_co2_surface_source_is_initialized = .true.

  end subroutine stream_co2_surface_source_init

  !================================================================
  subroutine stream_co2_surface_source_interp(cam_out, rc)

    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod , only : shr_strdata_advance
    use camsrfexch       , only : cam_out_t
    use ppgrid           , only : begchunk, endchunk
    use time_manager     , only : get_curr_date
    use phys_grid        , only : get_ncols_p

    ! input/output variables
    type(cam_out_t) , intent(inout)  :: cam_out(begchunk:endchunk)
    integer         , intent(out)    :: rc

    ! local variables
    integer  :: ig,icol,lchnk
    integer  :: year    ! year (0, ...) for nstep+1
    integer  :: mon     ! month (1, ..., 12) for nstep+1
    integer  :: day     ! day of month (1, ..., 31) for nstep+1
    integer  :: sec     ! seconds into current date for nstep+1
    integer  :: mcdate  ! Current model date (yyyymmdd)
    real(r8), pointer :: dataptr1d(:)
    character(*), parameter :: subName = "('stream_co2_surface_source_interp')"
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Advance co2_surface_source sdat stream
    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day
    call shr_strdata_advance(sdat_co2_surface_source, ymd=mcdate, tod=sec, logunit=iulog, &
         istr='co2_surface_source_diag', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get pointer for stream data that is time and spatially interpolated to model time and grid
    call dshr_fldbun_getFldPtr(sdat_co2_surface_source%pstrm(1)%fldbun_model, &
         stream_co2_surface_source_data_varname, fldptr1=dataptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set output diagnostic co2_surface_source
    ! Input data is in [mol/mol] but land and ocean expect to receive [ppm].
    ! Unit conversion is to multiply by 1e6.
    ig = 1
    do lchnk = begchunk,endchunk
       do icol = 1,get_ncols_p(lchnk)
          cam_out(lchnk)%co2diag(icol) = dataptr1d(ig) * 1.0e6_r8
          ig = ig + 1
       end do
    end do

  end subroutine stream_co2_surface_source_interp

end module atm_stream_co2
