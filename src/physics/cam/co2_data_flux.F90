module co2_data_flux

   !-------------------------------------------------------------------------------
   ! read and interpolate co2 surface fluxes
   !-------------------------------------------------------------------------------

   use shr_kind_mod,     only: r8=>shr_kind_r8, cl=>shr_kind_cl, cs=> shr_kind_cs
   use ESMF,             only: ESMF_Mesh, ESMF_Finalize, ESMF_LogFoundError
   use ESMF,             only: ESMF_END_ABORT, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
   use cam_logfile,      only: iulog
   use cam_abortutils,   only: endrun
   use spmd_utils,       only: iam
   use dshr_strdata_mod, only: shr_strdata_type

   implicit none
   private

   ! Public interfaces
   public co2_data_flux_type
   public co2_data_flux_init
   public co2_data_flux_advance

   type :: co2_data_flux_type
      type(shr_strdata_type) :: sdat_co2
      character(len=cs)      :: varname
      real(r8), pointer      :: co2flx(:,:)  ! Interpolated output (pcols,begchunk:endchunk)
   end type co2_data_flux_type

!===============================================================================
contains
!===============================================================================

   subroutine co2_data_flux_init (input_file, input_meshfile, &
        varname, year_first, year_last, year_align, tintalgo, taxmode, data_flux)

      !-------------------------------------------------------------------------------
      ! Initialize co2_data_flux_type instance
      !   including initial read of input and interpolation to the current timestep
      !-------------------------------------------------------------------------------

      use ppgrid,           only: begchunk, endchunk, pcols
      use atm_shr,          only: model_mesh, model_clock
      use dshr_strdata_mod, only: shr_strdata_init_from_inline

      ! Arguments
      character(len=*),         intent(in)    :: input_file ! assumes only one input file
      character(len=*),         intent(in)    :: input_meshfile
      character(len=*),         intent(in)    :: varname    ! assume only one varname for sdat
      integer,                  intent(in)    :: year_first
      integer,                  intent(in)    :: year_last
      integer,                  intent(in)    :: year_align
      character(len=*),         intent(in)    :: tintalgo
      character(len=*),         intent(in)    :: taxmode
      type(co2_data_flux_type), intent(inout) :: data_flux

      ! Local variables
      integer :: rc
      character(len=*), parameter :: subname = 'co2_data_flux_init'
      !----------------------------------------------------------------------------

      ! Initialize data_flux%sdat_co2
      call shr_strdata_init_from_inline(data_flux%sdat_co2, &
           my_task             = iam,                       &
           logunit             = iulog,                     &
           compname            = 'ATM',                     &
           model_clock         = model_clock,               &
           model_mesh          = model_mesh,                &
           stream_meshfile     = trim(input_meshfile),      &
           stream_filenames    = (/input_file/),            &
           stream_yearFirst    = year_first,                &
           stream_yearLast     = year_last,                 &
           stream_yearAlign    = year_align,                &
           stream_fldlistFile  = (/varname/),               &
           stream_fldListModel = (/varname/),               &
           stream_lev_dimname  = 'null',                    &
           stream_mapalgo      = 'bilinear',                &
           stream_offset       = 0,                         &
           stream_taxmode      = trim(taxmode),             &
           stream_dtlimit      = 1.0e30_r8,                 &
           stream_tintalgo     = trim(tintalgo),            &
           stream_name         = 'CO2 forcing data ',       &
           rc                  = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      ! Initialize data_flux%varname
      data_flux%varname = trim(varname)

      ! Initialize data_flux%co2flx
      allocate( data_flux%co2flx(pcols,begchunk:endchunk) )

   end subroutine co2_data_flux_init

   !===============================================================================
   subroutine co2_data_flux_advance (data_flux)

      !-------------------------------------------------------------------------------
      ! Advance the contents of a co2_data_flux_type sdat (map and interpolate in time)
      !-------------------------------------------------------------------------------

      use dshr_methods_mod , only : dshr_fldbun_getfldptr
      use dshr_strdata_mod , only : shr_strdata_advance
      use ppgrid           , only : begchunk, endchunk
      use phys_grid        , only : get_ncols_p
      use time_manager     , only : get_curr_date

      ! Arguments
      type(co2_data_flux_type),  intent(inout) :: data_flux

      ! Local variables
      integer :: icol,lchnk,g
      integer :: year    ! year (0, ...) for nstep+1
      integer :: mon     ! month (1, ..., 12) for nstep+1
      integer :: day     ! day of month (1, ..., 31) for nstep+1
      integer :: sec     ! seconds into current date for nstep+1
      integer :: mcdate  ! Current model date (yyyymmdd)
      integer :: rc
      real(r8), pointer :: dataptr1d(:)
      character(len=*), parameter :: subname = 'co2_data_flux_advance'
      !----------------------------------------------------------------------------

      ! Advance sdat stream
      call get_curr_date(year, mon, day, sec)
      mcdate = year*10000 + mon*100 + day
      call shr_strdata_advance(data_flux%sdat_co2, ymd=mcdate, tod=sec, logunit=iulog, istr='co2_advance', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      ! Get pointer for stream data that is time and spatially interpolated to model time and grid
      call dshr_fldbun_getFldPtr(data_flux%sdat_co2%pstrm(1)%fldbun_model, data_flux%varname, fldptr1=dataptr1d, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      g = 1
      do lchnk = begchunk,endchunk
         do icol = 1,get_ncols_p(lchnk)
            data_flux%co2flx(icol,lchnk) = dataptr1d(g)
            g = g + 1
         end do
      end do

   end subroutine co2_data_flux_advance

end module co2_data_flux
