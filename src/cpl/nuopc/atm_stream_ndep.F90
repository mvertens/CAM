module atm_stream_ndep

  !-----------------------------------------------------------------------
  ! Contains methods for reading in nitrogen deposition data file
  ! Also includes functions for dynamic ndep file handling and
  ! interpolation.
  !-----------------------------------------------------------------------
  !
  use ESMF              , only : ESMF_Clock, ESMF_Mesh
  use ESMF              , only : ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
  use ESMF              , only : ESMF_Finalize, ESMF_LogFoundError
  use nuopc_shr_methods , only : chkerr
  use dshr_strdata_mod  , only : shr_strdata_type
  use shr_kind_mod      , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use spmd_utils        , only : mpicom, masterproc, iam
  use spmd_utils        , only : mpi_character, mpi_integer
  use cam_logfile       , only : iulog
  use cam_abortutils    , only : endrun
  use cam_esmf_mod      , only : model_clock, model_mesh

  implicit none
  private

  public :: stream_ndep_readnl    ! read runtime options
  public :: stream_ndep_init      ! position datasets for dynamic ndep
  public :: stream_ndep_interp    ! interpolates between two years of ndep file data

  private :: stream_ndep_check_units   ! Check the units and make sure they can be used

  ! The ndep stream is not needed for aquaplanet or simple model configurations.  It
  ! is disabled by setting the namelist variable stream_ndep_data_filename to 'UNSET' or empty string.
  logical, public, protected :: ndep_stream_active = .false.
  logical, public, protected :: stream_ndep_is_initialized = .false.

  type(shr_strdata_type) :: sdat_ndep     ! input data stream

  character(len=*), parameter :: sourcefile = __FILE__

  ! namelist variables
  character(len=CL) :: stream_ndep_data_filename
  character(len=CL) :: stream_ndep_mesh_filename
  character(len=CL) :: stream_ndep_varlist    ! colon delimited string of ndep field names
  integer           :: stream_ndep_year_first ! first year in stream to use
  integer           :: stream_ndep_year_last  ! last year in stream to use
  integer           :: stream_ndep_year_align ! align stream_year_firstndep with

  character(len=CS), allocatable :: stream_ndep_varnames(:) ! array of ndep field names

!==============================================================================
contains
!==============================================================================

  subroutine stream_ndep_readnl(nlfile)

    ! Uses:
    use shr_nl_mod,     only: shr_nl_find_group_name
    use shr_string_mod, only: shr_string_listGetNum, shr_string_listGetName
    use error_messages, only: alloc_err

    ! input/output variables
    character(len=*), intent(in) :: nlfile

    ! local variables
    integer :: nu_nml    ! unit for namelist file
    integer :: nml_error ! namelist i/o error flag
    integer :: ierr      ! error status
    integer :: nf        ! field counter
    integer :: numflds   ! number of fields in stream_ndep_varlist
    character(*), parameter :: subName = "('stream_ndep_readnl')"
    !-----------------------------------------------------------------------

    namelist /ndep_stream_nl/       &
         stream_ndep_data_filename, &
         stream_ndep_mesh_filename, &
         stream_ndep_year_first,    &
         stream_ndep_year_last,     &
         stream_ndep_year_align,    &
         stream_ndep_varlist

    ! Default values for namelist
    stream_ndep_data_filename = ' '
    stream_ndep_mesh_filename = ' '
    stream_ndep_varlist       = ' '
    stream_ndep_year_first    = 1 ! first year in stream to use
    stream_ndep_year_last     = 1 ! last  year in stream to use
    stream_ndep_year_align    = 1 ! align stream_ndep_year_first with this model year

    ! Read ndep_stream namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(nlfile), status='old', iostat=nml_error )
       if (nml_error /= 0) then
          call endrun(subName//': ERROR opening '//trim(nlfile)//errMsg(sourcefile, __LINE__))
       end if
       call shr_nl_find_group_name(nu_nml, 'ndep_stream_nl', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=ndep_stream_nl, iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(' ERROR reading ndep_stream_nl namelist'//errMsg(sourcefile, __LINE__))
          end if
       end if
       close(nu_nml)
    endif
    call mpi_bcast(stream_ndep_mesh_filename, len(stream_ndep_mesh_filename), mpi_character, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_ndep_mesh_filename")
    call mpi_bcast(stream_ndep_data_filename, len(stream_ndep_data_filename), mpi_character, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_ndep_data_filename")
    call mpi_bcast(stream_ndep_varlist, len(stream_ndep_varlist), mpi_character, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_ndep_varlist")
    call mpi_bcast(stream_ndep_year_first, 1, mpi_integer, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_ndep_year_first")
    call mpi_bcast(stream_ndep_year_last, 1, mpi_integer, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_ndep_year_last")
    call mpi_bcast(stream_ndep_year_align, 1, mpi_integer, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//": FATAL: mpi_bcast: stream_ndep_year_align")

    ! Determine if ndep stream is active, and if not return
    ndep_stream_active = (len_trim(stream_ndep_data_filename)>0 .and. stream_ndep_data_filename/='UNSET')

    ! Check whether the stream is being used.
    if (.not. ndep_stream_active) then
       if (masterproc) then
          write(iulog,'(a)') ' '
          write(iulog,'(2a)') trim(subname),' NDEP STREAM IS NOT USED.'
          write(iulog,'(a)') ' '
       endif
       RETURN
    endif

    ! Create array of variable names on ndep forcing file - needed to initialize sdat
    numflds = shr_string_listGetNum(stream_ndep_varlist)
    allocate(stream_ndep_varnames(numflds), stat=ierr)
    call alloc_err(ierr, subname, 'varnames', numflds)
    do nf = 1,numflds
       call shr_string_listGetName(stream_ndep_varlist, nf, stream_ndep_varnames(nf))
    end do

    if (masterproc) then
       write(iulog,'(a)'   ) ' '
       write(iulog,'(2a,i0)') trim(subname),' stream ndep settings:'
       write(iulog,'(3a)')    trim(subname),'  stream_ndep_data_filename = ',trim(stream_ndep_data_filename)
       write(iulog,'(3a)')    trim(subname),'  stream_ndep_mesh_filename = ',trim(stream_ndep_mesh_filename)
       write(iulog,'(3a)')    trim(subname),'  stream_ndep_varlist       = ',trim(stream_ndep_varlist)
       write(iulog,'(2a,i0)') trim(subname),'  stream_ndep_year_first    = ',stream_ndep_year_first
       write(iulog,'(2a,i0)') trim(subname),'  stream_ndep_year_last     = ',stream_ndep_year_last
       write(iulog,'(2a,i0)') trim(subname),'  stream_ndep_year_align    = ',stream_ndep_year_align
       write(iulog,'(a)'   )  ' '
    endif

  end subroutine stream_ndep_readnl

  !================================================================
  subroutine stream_ndep_init(rc)
    use dshr_strdata_mod, only: shr_strdata_init_from_inline

    ! input/output variables
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: subName = "('stream_ndep_init')"
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (.not.ndep_stream_active) then
       RETURN
    end if

    ! Read in units
    call stream_ndep_check_units(stream_ndep_data_filename)

    ! Initialize the cdeps data type sdat_ndep
    call shr_strdata_init_from_inline(sdat_ndep,                    &
         my_task             = iam,                                 &
         logunit             = iulog,                               &
         compname            = 'ATM',                               &
         model_clock         = model_clock,                         &
         model_mesh          = model_mesh,                          &
         stream_meshfile     = trim(stream_ndep_mesh_filename),     &
         stream_filenames    = (/trim(stream_ndep_data_filename)/), &
         stream_yearFirst    = stream_ndep_year_first,              &
         stream_yearLast     = stream_ndep_year_last,               &
         stream_yearAlign    = stream_ndep_year_align,              &
         stream_fldlistFile  = stream_ndep_varnames,                &
         stream_fldListModel = stream_ndep_varnames,                &
         stream_lev_dimname  = 'null',                              &
         stream_mapalgo      = 'bilinear',                          &
         stream_offset       = 0,                                   &
         stream_taxmode      = 'cycle',                             &
         stream_dtlimit      = 1.0e30_r8,                           &
         stream_tintalgo     = 'linear',                            &
         stream_name         = 'Nitrogen deposition data ',         &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    stream_ndep_is_initialized = .true.

  end subroutine stream_ndep_init

  !================================================================
  subroutine stream_ndep_check_units( stream_fldFileName_ndep)

    !--------------------------------------------------------
    ! Check that units are correct on the file and if need any conversion
    !--------------------------------------------------------

     use cam_pio_utils , only : cam_pio_createfile, cam_pio_openfile, cam_pio_closefile, pio_subsystem
     use pio           , only : file_desc_t, io_desc_t, var_desc_t, pio_double, pio_def_dim
     use pio           , only : pio_bcast_error, pio_seterrorhandling, pio_inq_varid, pio_get_att
     use pio           , only : PIO_NOERR, PIO_NOWRITE

    ! Arguments
    character(len=*), intent(in)  :: stream_fldFileName_ndep  ! ndep filename
    !
    ! Local variables
    type(file_desc_t) :: fileid       ! NetCDF filehandle for ndep file
    type(var_desc_t)  :: vardesc      ! variable descriptor
    integer           :: ierr         ! error status
    integer           :: err_handling ! temporary
    character(len=CS) :: ndepunits    ! ndep units
    character(*), parameter :: subName = "('stream_ndep_check_units')"
    !-----------------------------------------------------------------------

    call cam_pio_openfile( fileid, trim(stream_fldFileName_ndep), PIO_NOWRITE)
    call pio_seterrorhandling(fileid, PIO_BCAST_ERROR, err_handling)
    ierr = pio_inq_varid(fileid, stream_ndep_varnames(1), vardesc)
    if (ierr /= PIO_NOERR) then
       call endrun(' ERROR finding variable: '//trim(stream_ndep_varnames(1))//" in file: "// &
            trim(stream_fldFileName_ndep)//errMsg(sourcefile, __LINE__))
    else
       ierr = PIO_get_att(fileid, vardesc, "units", ndepunits)
    end if
    call pio_seterrorhandling(fileid, err_handling)
    call cam_pio_closefile(fileid)

    select case (trim(stream_ndep_varlist))
    case ('NDEP_NHx_month:NDEP_NOy_month')
       ! Now check to make sure they are correct
       if (.not. trim(ndepunits) == "g(N)/m2/s"  )then
          call endrun(' ERROR in units for nitrogen deposition equal to: '//trim(ndepunits)//" not units expected"// &
               errMsg(sourcefile, __LINE__))
       end if
    case ('drynhx:wetnhx:drynoy:wetnoy')
       if (.not. trim(ndepunits) == "kg m-2 s-1")then
          call endrun(' ERROR in units for nitrogen deposition equal to: '//trim(ndepunits)//" not units expected"// &
               errMsg(sourcefile, __LINE__))
       end if
    case default
       call endrun(trim(subname)//'stream_ndep_varlist '//trim(stream_ndep_varlist)//' is not supported')
    end select

  end subroutine stream_ndep_check_units

  !================================================================
  subroutine stream_ndep_interp(cam_out, rc)

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
    real(r8) :: scale_ndep
    real(r8), pointer :: dataptr1d_nhx(:)
    real(r8), pointer :: dataptr1d_noy(:)
    real(r8), pointer :: dataptr1d_nhx_dry(:)
    real(r8), pointer :: dataptr1d_nhx_wet(:)
    real(r8), pointer :: dataptr1d_noy_dry(:)
    real(r8), pointer :: dataptr1d_noy_wet(:)
    character(*), parameter :: subName = "('stream_ndep_interp')"
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (.not.ndep_stream_active) then
       RETURN
    end if

    ! Advance sdat stream
    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day
    call shr_strdata_advance(sdat_ndep, ymd=mcdate, tod=sec, logunit=iulog, istr='ndepdyn', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Get pointer for stream data that is time and spatially interpolated to model time and grid
    select case (trim(stream_ndep_varlist))
    case ('NDEP_NHx_month:NDEP_NOy_month')

       call dshr_fldbun_getFldPtr(sdat_ndep%pstrm(1)%fldbun_model, 'NDEP_NHX_month', fldptr1=dataptr1d_nhx, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       call dshr_fldbun_getFldPtr(sdat_ndep%pstrm(1)%fldbun_model, 'NDEP_NOy_month', fldptr1=dataptr1d_noy, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if

       ! NDEP read from forcing is in units of gN/m2/sec - but the mediator
       ! expects units of kgN/m2/sec
       scale_ndep = .001_r8
       ig = 1
       do lchnk = begchunk,endchunk
          do icol = 1,get_ncols_p(lchnk)
             cam_out(lchnk)%nhx_nitrogen_flx(icol) = dataptr1d_nhx(ig) * scale_ndep
             cam_out(lchnk)%noy_nitrogen_flx(icol) = dataptr1d_noy(ig) * scale_ndep
             ig = ig + 1
          end do
       end do

    case ('drynhx:wetnhx:drynoy:wetnoy')

       call dshr_fldbun_getFldPtr(sdat_ndep%pstrm(1)%fldbun_model, 'drynhx', fldptr1=dataptr1d_nhx_dry, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       call dshr_fldbun_getFldPtr(sdat_ndep%pstrm(1)%fldbun_model, 'wetnhx', fldptr1=dataptr1d_nhx_wet, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       call dshr_fldbun_getFldPtr(sdat_ndep%pstrm(1)%fldbun_model, 'drynoy', fldptr1=dataptr1d_noy_dry, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       call dshr_fldbun_getFldPtr(sdat_ndep%pstrm(1)%fldbun_model, 'wetnoy', fldptr1=dataptr1d_noy_wet, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if

       ! NDEP read from forcing is in units of kgN/m2/sec
       ig = 1
       do lchnk = begchunk,endchunk
          do icol = 1,get_ncols_p(lchnk)
             cam_out(lchnk)%nhx_nitrogen_flx(icol) = dataptr1d_nhx_dry(ig) + dataptr1d_nhx_wet(ig)
             cam_out(lchnk)%noy_nitrogen_flx(icol) = dataptr1d_noy_dry(ig) + dataptr1d_noy_wet(ig)
             ig = ig + 1
          end do
       end do

    end select

  end subroutine stream_ndep_interp

end module atm_stream_ndep
