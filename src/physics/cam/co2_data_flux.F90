module co2_data_flux

!-------------------------------------------------------------------------------
! utilities for reading and interpolating co2 surface fluxes
!-------------------------------------------------------------------------------

   use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
   use cam_abortutils,   only: endrun

   implicit none
   private

   ! Public interfaces
   public co2_data_flux_type

   public co2_data_flux_init
   public co2_data_flux_advance

   type :: co2_data_flux_type
      type(shr_strdata_type) :: sdat_co2
      real(r8), pointer      :: co2flx(:,:)  ! Interpolated output (pcols,begchunk:endchunk)
   end type co2_data_flux_type

   ! dimension names for physics grid (physgrid)
   logical           :: dimnames_set = .false.
   character(len=8)  :: dim1name, dim2name

!===============================================================================
contains
!===============================================================================

  subroutine co2_data_flux_init (input_file, input_mesh, &
       varname, year_first, year_last, year_align, tintalgo, taxmode, xin)

!-------------------------------------------------------------------------------
! Initialize co2_data_flux_type instance
!   including initial read of input and interpolation to the current timestep
!-------------------------------------------------------------------------------

   use ESMF,             only: ESMF_Mesh
   use ppgrid,           only: begchunk, endchunk, pcols
   use cam_grid_support, only: cam_grid_id, cam_grid_check
   use cam_grid_support, only: cam_grid_get_dim_names
   use atm_shr         , only: model_mesh, model_clock

   ! Arguments
   character(len=*),         intent(in)    :: input_file
   type(ESMF_Mesh) ,         intent(in)    :: input_mesh
   character(len=*),         intent(in)    :: varname
   integer,                  intent(in)    :: year_first      
   integer,                  intent(in)    :: year_last
   integer,                  intent(in)    :: year_align
   character(len=*),         intent(in)    :: tintalgo
   character(len=*),         intent(in)    :: taxalgo
   type(co2_data_flux_type), intent(inout) :: xin

   ! Local variables
   integer  :: grid_id
   real(r8) :: dtime
   character(len=*), parameter :: subname = 'co2_data_flux_init'
   !----------------------------------------------------------------------------

   if (.not. dimnames_set) then
      grid_id = cam_grid_id('physgrid')
      if (.not. cam_grid_check(grid_id)) then
         call endrun(subname // ': ERROR: no "physgrid" grid')
      endif
      call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
      dimnames_set = .true.
   end if

   call shr_strdata_init_from_inline(xin%sdat_co2,  &
         my_task             = iam,                 &
         logunit             = iulog,               &
         compname            = 'ATM',               &
         model_clock         = model_clock,         &
         model_mesh          = model_mesh,          &
         stream_meshfile     = trim(input_mesh),    &
         stream_filenames    = (/input_file/),      &
         stream_yearFirst    = year_first,          &
         stream_yearLast     = year_last,           &
         stream_yearAlign    = year_first,          &
         stream_fldlistFile  = (/varname/),         &
         stream_fldListModel = (/varname/),         &
         stream_lev_dimname  = 'null',              &
         stream_mapalgo      = 'bilinear',          &
         stream_offset       = 0,                   &
         stream_taxmode      = trim(taxmode),       &
         stream_dtlimit      = 1.0e30_r8,           &
         stream_tintalgo     = trim(tintalgo),      &
         stream_name         = 'CO2 forcing data ', &
         rc                  = rc)
   call chkrc(rc, sub//': error return from shr_strdata_init_from_inline')

   allocate( xin%co2flx(pcols,begchunk:endchunk) )

   call co2_data_flux_advance(xin)

end subroutine co2_data_flux_init

!===============================================================================
subroutine co2_data_flux_advance (xin)

!-------------------------------------------------------------------------------
! Advance the contents of a co2_data_flux_type sdat (map and interpolate in time) 
!-------------------------------------------------------------------------------

    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod , only : shr_strdata_advance
    use ppgrid           , only : begchunk, endchunk
    use phys_grid        , only : get_ncols_p
    use time_manager     , only : get_curr_date

   ! Arguments
   type(co2_data_flux_type),  intent(inout) :: xin

   ! Local variables
   integer :: icol,lchnk,g
   integer :: year    ! year (0, ...) for nstep+1
   integer :: mon     ! month (1, ..., 12) for nstep+1
   integer :: day     ! day of month (1, ..., 31) for nstep+1
   integer :: sec     ! seconds into current date for nstep+1
   integer :: mcdate  ! Current model date (yyyymmdd)
   real(r8), pointer :: dataptr1d(:)
   character(len=*), parameter :: subname = 'co2_data_flux_advance'
   !----------------------------------------------------------------------------


    ! Advance sdat stream
    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day
    call shr_strdata_advance(sdat_ndep, ymd=mcdate, tod=sec, logunit=iulog, istr='ndepdyn', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Get pointer for stream data that is time and spatially interpolated to model time and grid
    call dshr_fldbun_getFldPtr(sdat_ndep%pstrm(1)%fldbun_model, stream_varlist_ndep(1), fldptr1=dataptr1d_nhx, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(sdat_ndep%pstrm(1)%fldbun_model, stream_varlist_ndep(2), fldptr1=dataptr1d_noy, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    g = 1
    do lchnk = begchunk,endchunk
       do icol = 1,get_ncols_p(lchnk)
          xin%co2flx(icol,lchnk) = dataptr1d(g)
          g = g + 1
       end do
    end do

  end subroutine co2_data_flux_advance

!===============================================================================

end module co2_data_flux
