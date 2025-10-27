module atm_stream_nudging

  !-----------------------------------------------------------------------
  ! Contains methods for reading in nitrogen deposition data file
  ! Also includes functions for dynamic nudging file handling and
  ! interpolation.
  !-----------------------------------------------------------------------
  !
  use ESMF              , only : ESMF_Clock, ESMF_Mesh
  use ESMF              , only : ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
  use ESMF              , only : ESMF_Finalize, ESMF_LogFoundError
  use ESMF              , only : ESMF_Time, ESMF_TimeInterval
  use ESMF              , only : ESMF_TimeGet, ESMF_TimeIntervalGet, ESMF_TimeIntervalSet
  use nuopc_shr_methods , only : chkerr
  use dshr_strdata_mod  , only : shr_strdata_type
  use shr_kind_mod      , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use spmd_utils        , only : masterproc, iam
  use cam_logfile       , only : iulog
  use cam_abortutils    , only : endrun
  use atm_shr           , only : model_mesh

  implicit none
  private

  public :: stream_nudging_init      ! position datasets for dynamic nudging
  public :: stream_nudging_interp    ! interpolates between two years of nudging file data

  type(shr_strdata_type) :: sdat_nudging

  character(len=2)       :: nudging_varlist(5) = (/'U ', 'V ','T ','Q ','PS'/)

  character(*),parameter :: u_FILE_u = __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine stream_nudging_init(nudge_path, nudge_files, nudge_mesh, &
       nudge_beg_time, nudge_end_time, model_update_interval, nudge_force_opt)

    use dshr_strdata_mod, only: shr_strdata_init_from_inline

    ! input/output arguments
    character(len=*)        , intent(in) :: nudge_path
    character(len=*)        , intent(in) :: nudge_files(:)
    character(len=*)        , intent(in) :: nudge_mesh
    type(ESMF_Time)         , intent(in) :: nudge_beg_time
    type(ESMF_Time)         , intent(in) :: nudge_end_time
    type(ESMF_TimeInterval) , intent(in  :: model_udpate_interval
    integer                 , intent(in) :: nudge_force_opt

    ! local variables
    integer                 :: rc
    integer                 :: nfile
    integer                 :: nudge_year_first
    integer                 :: nudge_year_last
    type(ESMF_Clock)        :: nudging_clock
    character(*), parameter :: sub = "('stream_nudging_init')"
    !----------------------------------------------------------------

    ! Create a nudging_clock for nudging - this is different than the CAM clock - it's time step is from the input
    ! nudging information

    call ESMF_TimeGet(nudge_beg_time, year=nudge_year_first, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    call ESMF_TimeGet(nudge_end_time, year=nudge_year_last, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! TODO: should this be initialized with a gregorian calendar
    ! the only use of the model clock in CDEPS is to extract the calendar

    nudging_clock = ESMF_ClockCreate(name="Nudging Model Clock", &
         model_update_interval, nudge_beg_time, stop_time=nudge_end_time, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Output info

    if (masterproc) then
       write(iulog,'(a)'   ) ' '
       write(iulog,'(a,i8)')  'stream nudging settings:'
       write(iulog,'(a,a)' )  '  nudge_mesh       = ',trim(nudge_mesh)
       write(iulog,'(a,a,a)') '  nudge_varlist    = ','U,V,T,Q,PS'
       write(iulog,'(a,i8)')  '  nudge_year_first = ',nudge_year_first
       write(iulog,'(a,i8)')  '  nudge_year_last  = ',nudge_year_last
       write(iulog,'(a,i8)')  '  nudge_year_align = ',nudge_year_first
       do nfile = 1,size(nudge_files)
          write(iulog,'(a,i8,a)' )  '  nudge_files = ',nfile,trim(nudge_files(nfile))
       end do
       write(iulog,'(a)'   )  ' '
    endif

    ! Create stream data type sdat_nudging

    if (Nudge_Force_Opt == 0) then
       tintalgo = 'upper'
    elseif(Nudge_Force_Opt == 1) then
       tintalgo = 'linear'
    else
       write(iulog,*) 'NUDGING: Unknown Nudge_Force_Opt=',Nudge_Force_Opt
       call endrun('nudging_timestep_init:: ERROR unknown Nudge_Force_Opt')
    endif

    ! Initialize the cdeps data type sdat_nudging
    call shr_strdata_init_from_inline(sdat_nudging,     &
         my_task             = iam,                     &
         logunit             = iulog,                   &
         compname            = 'ATM',                   &
         model_clock         = nudge_clock,             &
         model_mesh          = model_mesh,              &
         stream_meshfile     = trim(nudge_mesh),        &
         stream_filenames    = nudge_files,             &
         stream_yearFirst    = nudge_year_first,        &
         stream_yearLast     = nudge_year_last,         &
         stream_yearAlign    = nudge_year_align,        &
         stream_fldlistFile  = nudge_varlist,           &
         stream_fldListModel = nudge_varlist,           &
         stream_lev_dimname  = 'null',                  &
         stream_mapalgo      = 'bilinear',              &
         stream_offset       = 0,                       &
         stream_taxmode      = 'limit',                 &
         stream_dtlimit      = 1.0e30_r8,               & ! change dtlimit to be twice the step size
         stream_tintalgo     = tintalgo,                &
         stream_name         = 'NUDGING forcing data ', &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

  end subroutine stream_nudging_init

  !================================================================

  subroutine stream_nudging_interp(Model_Update_Time, &
       Nudge_ZonalFilter, ZM, Zonal_Bamp2d, Zonal_Bamp3d, &
       Target_U, Target_V, Target_T, Target_Q, Target_PS)

    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod , only : shr_strdata_advance
    use ppgrid           , only : pcols, pver, begchunk,endchunk
    use ppgrid           , only : begchunk, endchunk
    use time_manager     , only : get_curr_date
    use phys_grid        , only : get_ncols_p
    use zonal_mean_mod   , only : ZonalMean_t
    use cam_abortutils   , only : endrun, handle_allocate_error

    ! input/output variables
    type(ESMF_Time)   , intent(in)  :: Model_Update_Time
    logical           , intent(in)  :: Nudge_ZonalFilter
    type(ZonalMean_t) , intent(in)  :: ZM
    real(r8)          , intent(in)  :: Zonal_Bamp2d(:)
    real(r8)          , intent(in)  :: Zonal_Bamp3d(:,:)
    real(r8)          , intent(out) :: Target_U(pcols,pver,begchunk:endchunk)
    real(r8)          , intent(out) :: Target_V(pcols,pver,begchunk:endchunk)
    real(r8)          , intent(out) :: Target_T(pcols,pver,begchunk:endchunk)
    real(r8)          , intent(out) :: Target_Q(pcols,pver,begchunk:endchunk)

    ! Local variables
    integer :: rc     ! ESMF error return
    integer :: istat  ! allocate return
    integer :: nvar   ! variable index
    integer :: ilev   ! level index
    integer :: icol   ! column index
    integer :: ichnk  ! chunk index
    integer :: g      ! counter index
    integer :: year   ! year (0, ...) for nstep+1
    integer :: mon    ! month (1, ..., 12) for nstep+1
    integer :: day    ! day of month (1, ..., 31) for nstep+1
    integer :: sec    ! seconds into current date for nstep+1
    integer :: mcdate ! current model date (yyyymmdd)
    real(r8), pointer    :: dataptr2d(:,:) ! first dimension is level, second is data on that level
    real(r8), pointer    :: dataptr1d(:)
    real(r8),allocatable :: Tmp3D(:,:,:)
    real(r8),allocatable :: Tmp2D(:,:)
    character(len=*), parameter :: sub = "(stream_nudging_interp) "
    !-----------------------------------------------------------------------

    ! Extract YMD from model_nudge_time
    call ESMF_TimeGet(Model_Update_Time, year=year, month=month, day=day, sec=sec, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    mcdate = year*10000 + mon*100 + day

    ! Advance sdat stream
    call shr_strdata_advance(sdat_nudging, ymd=mcdate, tod=sec, logunit=iulog, istr='nudging', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Get pointer for stream data that is time and spatially interpolated to model time and grid
    allocate(Tmp3D(pcols,pver,begchunk:endchunk), stat=istta)
    call handle_allocate_error(istat, sub, 'TMP3d')

    allocate(Tmp2D(pcols,begchunk:endchunk))
    call handle_allocate_error(istat, sub, 'TM23d')

    ! Determine 3d nudging fields
    do nvar = 1,4

       if ( trim(stream_varlist_nudging(nvar)) == 'U' .or. &
            trim(stream_varlist_nudging(nvar)) == 'V' .or. &
            trim(stream_varlist_nudging(nvar)) == 'T' .or. &
            trim(stream_varlist_nudging(nvar)) == 'Q' )  then

          call dshr_fldbun_getFldPtr(sdat_nudging%pstrm(1)%fldbun_model, stream_varlist_nudging(nvar), fldptr2=dataptr2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if

          ! Obtain TMP3d
          g = 1
          do ichnk = begchunk,endchunk
             do ilev = 1, plev
                do icol = 1,get_ncols_p(c)
                   Tmp3d(icol,ilev,ichnk) = dataptr2d(ilev,g)
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
          if (trim(stream_varlist_nudging(nvar) == 'U')) then
             do lchnk = begchunk,endchunk
                ncol = phys_state(lchnk)%ncol
                Target_U(:ncol,:pver,lchnk) = Tmp3d(:ncol,:pver,lchnk)
             end do
          else if (trim(stream_varlist_nudging(nvar) == 'V')) then
             do lchnk = begchunk,endchunk
                ncol = phys_state(lchnk)%ncol
                Target_V(:ncol,:pver,lchnk) = Tmp3d(:ncol,:pver,lchnk)
             end do
          else if (trim(stream_varlist_nudging(nvar) == 'T')) then
             do lchnk = begchunk,endchunk
                ncol = phys_state(lchnk)%ncol
                Target_T(:ncol,:pver,lchnk) = Tmp3d(:ncol,:pver,lchnk)
             end do
          else if (trim(stream_varlist_nudging(nvar) == 'Q')) then
             do lchnk = begchunk,endchunk
                ncol = phys_state(lchnk)%ncol
                Target_Q(:ncol,:pver,lchnk) = Tmp3d(:ncol,:pver,lchnk)
             end do
          end if

       else if (trim(stream_varlist_nudging(nvar)) == 'PS') then

          call dshr_fldbun_getFldPtr(sdat_nudging%pstrm(1)%fldbun_model, stream_varlist_nudging(nvar), fldptr2=dataptr1d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if

          g = 1
          do ichnk = begchunk,endchunk
             do icol = 1,get_ncols_p(c)
                Tmp2d(icol, ichnk) = dataptr2d(g)
                g = g + 1
             end do
          end do

          if (Nudge_ZonalFilter) then
             call ZM%calc_amps(Tmp2D,Zonal_Bamp2d)
             call ZM%eval_grid(Zonal_Bamp2d,Tmp2D)
          endif

          do lchnk=begchunk,endchunk
             ncol=phys_state(lchnk)%ncol
             Target_PS(:ncol,lchnk)= Tmp3d(:ncol,lchnk)
          end do

       end if !

    end do

  end subroutine stream_nudging_interp

end module atm_stream_nudging
