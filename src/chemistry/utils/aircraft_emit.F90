module aircraft_emit
   !-----------------------------------------------------------------------
   !
   ! Purpose:
   ! Manages reading and interpolation of aircraft aerosols
   !
   ! Authors:
   !   Chih-Chieh (Jack) Chen and Cheryl Craig -- February 2010
   !   Mariana Vertenstein (Refactored using CDEPS in-line functionality) -- November 2025
   !
   !-----------------------------------------------------------------------

   use perf_mod,         only : t_startf, t_stopf
   use shr_kind_mod,     only : r8 => shr_kind_r8, cl=>shr_kind_cl, cs=>shr_kind_cs
   use cam_abortutils,   only : endrun
   use cam_logfile,      only : iulog
   use spmd_utils,       only : masterproc, iam
   use dshr_strdata_mod, only : shr_strdata_type

   implicit none
   private

   public :: aircraft_emit_init
   public :: aircraft_emit_adv
   public :: aircraft_emit_register
   public :: aircraft_emit_readnl
   public :: get_aircraft

   private :: get_vertical_dimension
   private :: interpz_conserve

   type :: forcing_type
      type(shr_strdata_type) :: sdat
      character(len=cs)      :: fldname     = 'unset '
      character(len=cs)      :: fldunits    = 'unset'
      character(len=cl)      :: datafile    = 'unset'
      character(len=cl)      :: meshfile    = 'unset'
      character(len=cs)      :: mapalgo     = 'consf'
      character(len=cs)      :: tintalgo    = 'lower'
      character(len=cs)      :: taxmode     = 'unset'
      integer                :: year_first  = -999
      integer                :: year_last   = -999
      integer                :: year_align  = -999
      integer                :: nilev       = -1
      integer                :: nlev        = -1
      integer                :: pbuf_index  = -1
      real(r8), pointer      :: altitude_int(:)
      real(r8), pointer      :: altitude_lev(:)
   end type forcing_type

   integer, parameter  :: N_AERO = 3
   type(forcing_type)  :: forcing(N_AERO)
   character(len=3)    :: mixtype(N_AERO) = 'wet'
   real(r8), parameter :: molmass(N_AERO) = 1._r8

   character(len=*),parameter :: u_FILE_u = __FILE__

!============================================================================
contains
!============================================================================

   subroutine aircraft_emit_readnl(nlfile)

      !-------------------------------------------------------------------
      ! **** Read in the aircraft_emit namelist *****
      !-------------------------------------------------------------------

      use namelist_utils, only: find_group_name
      use spmd_utils,     only: mpicom, masterprocid
      use spmd_utils,     only: mpi_integer, mpi_logical, mpi_character
      use co2_cycle,      only: co2_readflux_aircraft
      use cam_pio_utils,  only: cam_pio_openfile
      use string_utils,   only: int2str
      use pio,            only: PIO_BCAST_ERROR, PIO_NOERR, PIO_NOWRITE
      use pio,            only: file_desc_t, pio_seterrorhandling, pio_inq_varid
      use pio,            only: pio_closefile

      ! Arguments
      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer           :: nf, ni
      integer           :: unitn, ierr
      type(file_desc_t) :: fileid
      integer           :: err_handling
      integer           :: varid
      logical           :: use_time_bnds

      character(len=cs) :: aircraft_co2_fldname          = 'ac_CO2'
      character(len=cl) :: aircraft_co2_datafile         = 'unset'
      character(len=cl) :: aircraft_co2_meshfile         = 'unset'
      character(len=cs) :: aircraft_co2_taxmode          = 'unset'
      character(len=cs) :: aircraft_co2_tintalgo         = 'unset'
      integer           :: aircraft_co2_year_first       = -999
      integer           :: aircraft_co2_year_last        = -999
      integer           :: aircraft_co2_year_align       = -999

      character(len=cs) :: aircraft_h2o_fldname          = 'ac_H2O'
      character(len=cl) :: aircraft_h2o_datafile         = 'unset'
      character(len=cl) :: aircraft_h2o_meshfile         = 'unset'
      character(len=cs) :: aircraft_h2o_taxmode          = 'unset'
      character(len=cs) :: aircraft_h2o_tintalgo         = 'unset'
      integer           :: aircraft_h2o_year_first       = -999
      integer           :: aircraft_h2o_year_last        = -999
      integer           :: aircraft_h2o_year_align       = -999

      character(len=cs) :: aircraft_slant_dist_fldname   = 'ac_SLANT_DIST'
      character(len=cl) :: aircraft_slant_dist_datafile  = 'unset'
      character(len=cl) :: aircraft_slant_dist_meshfile  = 'unset'
      character(len=cs) :: aircraft_slant_dist_tintalgo  = 'unset'
      character(len=cs) :: aircraft_slant_dist_taxmode   = 'unset'
      integer           :: aircraft_slant_dist_year_first= -999
      integer           :: aircraft_slant_dist_year_last = -999
      integer           :: aircraft_slant_dist_year_align= -999

      character(len=*), parameter :: subname = 'aircraft_emit_readnl'

      namelist /aircraft_emit_nl/  &
           aircraft_co2_datafile, aircraft_co2_meshfile, &
           aircraft_co2_year_first, aircraft_co2_year_last, aircraft_co2_year_align, &
           aircraft_co2_taxmode, aircraft_co2_tintalgo, &
           aircraft_h2o_datafile, aircraft_h2o_meshfile, &
           aircraft_h2o_year_first, aircraft_h2o_year_last, aircraft_h2o_year_align, &
           aircraft_h2o_taxmode, aircraft_h2o_tintalgo, &
           aircraft_slant_dist_datafile, aircraft_slant_dist_meshfile, &
           aircraft_slant_dist_year_first, aircraft_slant_dist_year_last, aircraft_slant_dist_year_align, &
           aircraft_slant_dist_taxmode, aircraft_slant_dist_tintalgo
      !-----------------------------------------------------------------------------

      ! Read namelist
      if (masterproc) then

         open( newunit=unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'aircraft_emit_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, aircraft_emit_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)

         ! Note - the following call assumes that co2_readflux_aircraft is
         ! set in co2_cycle_readnl and this is called before this routine in
         ! runtime_opts.F90. If co2_readflux_aircraft is .false. then, the
         ! forcing(nf)%datafile = 'unset' and this logic will be triggered
         ! in the other routines in this module
         if (co2_readflux_aircraft) then
            if (trim(aircraft_co2_datafile) /= 'unset') then
               nf = 1
               forcing(nf)%fldname    = aircraft_co2_fldname
               forcing(nf)%datafile   = aircraft_co2_datafile
               forcing(nf)%meshfile   = aircraft_co2_meshfile
               forcing(nf)%year_first = aircraft_co2_year_first
               forcing(nf)%year_last  = aircraft_co2_year_last
               forcing(nf)%year_align = aircraft_co2_year_align
               forcing(nf)%taxmode    = aircraft_co2_taxmode
               forcing(nf)%tintalgo   = aircraft_co2_tintalgo
            end if
         end if
         if (trim(aircraft_h2o_datafile) /= 'unset') then
            nf = 2
            forcing(nf)%datafile   = aircraft_h2o_datafile
            forcing(nf)%fldname    = aircraft_h2o_fldname
            forcing(nf)%meshfile   = aircraft_h2o_meshfile
            forcing(nf)%year_first = aircraft_h2o_year_first
            forcing(nf)%year_last  = aircraft_h2o_year_last
            forcing(nf)%year_align = aircraft_h2o_year_align
            forcing(nf)%taxmode    = aircraft_h2o_taxmode
            forcing(nf)%tintalgo   = aircraft_h2o_tintalgo
         end if
         if (trim(aircraft_slant_dist_datafile) /= 'unset') then
            nf = 3
            forcing(nf)%datafile   = aircraft_slant_dist_datafile
            forcing(nf)%fldname    = aircraft_slant_dist_fldname
            forcing(nf)%meshfile   = aircraft_slant_dist_meshfile
            forcing(nf)%year_first = aircraft_slant_dist_year_first
            forcing(nf)%year_last  = aircraft_slant_dist_year_last
            forcing(nf)%year_align = aircraft_slant_dist_year_align
            forcing(nf)%taxmode    = aircraft_slant_dist_taxmode
            forcing(nf)%tintalgo   = aircraft_slant_dist_tintalgo
         end if

      end if

      n_aero_loop: do nf = 1,N_AERO

         ! Broadcast namelist variables
         call mpi_bcast(forcing(nf)%datafile, len(forcing(nf)%datafile), mpi_character, masterprocid, mpicom, ierr)
         if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: forcing("//int2str(nf)//"%datapath")
         call mpi_bcast(forcing(nf)%fldname,len(forcing(nf)%fldname), mpi_character, masterprocid, mpicom, ierr)
         if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: forcing("//int2str(nf)//"%fldname")
         call mpi_bcast(forcing(nf)%meshfile, len(forcing(nf)%meshfile), mpi_character, masterprocid, mpicom, ierr)
         if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: forcing("//int2str(nf)//"%meshfile")
         call mpi_bcast(forcing(nf)%year_first, 1, mpi_integer, masterprocid, mpicom, ierr)
         if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: forcing("//int2str(nf)//"%year_first")
         call mpi_bcast(forcing(nf)%year_last, 1, mpi_integer, masterprocid, mpicom, ierr)
         if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: forcing("//int2str(nf)//"%year_last")
         call mpi_bcast(forcing(nf)%year_align, 1, mpi_integer, masterprocid, mpicom, ierr)
         if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: forcing("//int2str(nf)//"%year_align")
         call mpi_bcast(forcing(nf)%tintalgo, len(forcing(nf)%tintalgo), mpi_character, masterprocid, mpicom, ierr)
         if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: forcing("//int2str(nf)//"%year_tintalgo")
         call mpi_bcast(forcing(nf)%taxmode, len(forcing(nf)%taxmode), mpi_character, masterprocid, mpicom, ierr)
         if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: forcing("//int2str(nf)//"%year_taxmode")

         datafile_isnot_unset: if (trim(forcing(nf)%datafile) /= 'unset') then
            ! overwrite mapalgo for ac_SLANT_DIST
            if ( trim(forcing(nf)%fldname) == 'ac_SLANT_DIST') then
               forcing(nf)%mapalgo = 'nn'
            end if

            ! Overwrite forcing(nf)%tintalgo if it is set to 'unset'
            ! Check if the data file has a time_bnds variable and if so set the time interpolation
            ! type to 'nearest' otherwise set it to 'linear'

            if (trim(forcing(nf)%tintalgo) == 'unset') then
               call cam_pio_openfile( fileid, forcing(nf)%datafile, PIO_NOWRITE )
               call pio_seterrorhandling( fileid, PIO_BCAST_ERROR, oldmethod=err_handling )
               ierr = pio_inq_varid( fileid, 'time_bnds', varid )
               call pio_seterrorhandling( fileid, err_handling)
               use_time_bnds = (ierr == PIO_NOERR)
               if (use_time_bnds) then
                  forcing(nf)%tintalgo = 'nearest'
               else
                  forcing(nf)%tintalgo = 'linear'
               end if
               call pio_closefile( fileid )
            end if

            !  diagnostics
            if (masterproc) then
               write(iulog,*) ' '
               write(iulog,'(2a)' ) ' aircraft init settings for: ',trim(forcing(nf)%fldname)
               write(iulog,'(2a)' ) '   aircraft datafile   = ',trim(forcing(nf)%datafile)
               write(iulog,'(2a)' ) '   aircraft meshfile   = ',trim(forcing(nf)%meshfile)
               write(iulog,'(2a)' ) '   aircraft mapalgo    = ',trim(forcing(nf)%mapalgo)
               write(iulog,'(2a)' ) '   aircraft tintalgo   = ',trim(forcing(nf)%tintalgo)
               write(iulog,'(2a)' ) '   aircraft taxmode    = ',trim(forcing(nf)%taxmode)
               write(iulog,'(a,i0)')'   aircraft year_first = ',forcing(nf)%year_first
               write(iulog,'(a,i0)')'   aircraft year_last  = ',forcing(nf)%year_last
               write(iulog,'(a,i0)')'   aircraft year_align = ',forcing(nf)%year_align
               write(iulog,*) ' '
            end if
         end if datafile_isnot_unset

      end do n_aero_loop

   end subroutine aircraft_emit_readnl

   !=========================================================================
   subroutine aircraft_emit_register()

      !------------------------------------------------------------------
      ! **** Add the aircraft aerosol data to the physics buffer ****
      !------------------------------------------------------------------
      use ppgrid,         only: pver, pcols
      use physics_buffer, only: pbuf_add_field, dtype_r8

      ! Local variables
      integer           :: nf
      !--------------------------------------------

      do nf = 1,N_AERO
         if (trim(forcing(nf)%datafile) /= 'unset') then
            ! Add fldname to pbuf and obtain pbuf_index
            call pbuf_add_field(forcing(nf)%fldname, 'physpkg', dtype_r8, (/pcols,pver/), &
                 forcing(nf)%pbuf_index)
         end if
      end do

   end subroutine aircraft_emit_register

   !=========================================================================
   subroutine aircraft_emit_init()

      !-------------------------------------------------------------------
      ! **** Initialize the aircraft aerosol data handling ****
      !-------------------------------------------------------------------
      use cam_history,    only: addfld, add_default
      use phys_control,   only: phys_getopts
      use physics_buffer, only: pbuf_get_chunk, pbuf_get_index
      use cam_pio_utils,  only: cam_pio_openfile, cam_pio_closefile
      use pio,            only: file_desc_t, var_desc_t
      use pio,            only: pio_inq_varid, pio_get_att
      use pio,            only: PIO_NOERR, PIO_NOWRITE

      ! Local variables
      type(file_desc_t) :: pioid
      type(var_desc_t)  :: varid
      integer           :: ierr, rc
      integer           :: klev
      integer           :: nf
      logical           :: history_chemistry
      character(len=*), parameter :: subname = 'aircraft_emit_init'
      !-----------------------------------------------

      call phys_getopts(history_chemistry_out=history_chemistry)

      loop_n_aero: do nf = 1,N_AERO
         if (trim(forcing(nf)%datafile) /= 'unset') then

            ! Open file
            call cam_pio_openfile( pioid, forcing(nf)%datafile, PIO_NOWRITE)

            ! Determine units
            ierr = pio_inq_varid( pioid, forcing(nf)%fldname, varid )
            if (ierr/=pio_noerr) then
               call endrun(trim(subname)//' Cannot find variable '//trim(forcing(nf)%fldname)// &
                    ' in file '//trim(forcing(nf)%datafile))
            endif
            ierr = pio_get_att( pioid, varid, 'units', forcing(nf)%fldunits)
            if (ierr/=pio_noerr) then
               call endrun(trim(subname)//' Cannot get attribute units from '//trim(forcing(nf)%fldname)// &
                    ' in file '//trim(forcing(nf)%datafile))
            endif

            ! Determine vertical levels - altitude_int and altitude_lev
            call get_vertical_dimension(fid=pioid, dname='altitude_int', dsize=forcing(nf)%nilev, &
                 data=forcing(nf)%altitude_int)
            call get_vertical_dimension(fid=pioid, dname='altitude'    , dsize=forcing(nf)%nlev , &
                 data=forcing(nf)%altitude_lev)

            ! Write out info to log file
            if (masterproc) then
               write(iulog,'(a)') trim(subname)// ' file: ',trim(forcing(nf)%datafile)
               write(iulog,'(a)')'  variable '//trim(forcing(nf)%fldname)
               write(iulog,'(a)')'  units '//trim(forcing(nf)%fldunits)
               write(iulog,'(a)')'  altitude levels '
               do klev=1,size(forcing(nf)%altitude_lev)
                  write(iulog,*)'  ',klev,forcing(nf)%altitude_lev(klev)
               end do
               write(iulog,'(a)')'  altitude interfaces '
               do klev=1,size(forcing(nf)%altitude_int)
                  write(iulog,*)'  ',klev,forcing(nf)%altitude_int(klev)
               end do
            end if

            ! Close file
            call cam_pio_closefile(pioid)

            ! Add field to cam history output
            call addfld(trim(forcing(nf)%fldname), (/ 'lev' /), 'A', trim(forcing(nf)%fldunits), &
                 'aircraft emission '//trim(forcing(nf)%fldname))
            if (history_chemistry) then
               call add_default( trim(forcing(nf)%fldname), 1, ' ' )
            end if

         end if
      end do loop_n_aero

   end subroutine aircraft_emit_init

   !=========================================================================
   subroutine aircraft_emit_adv( state, pbuf2d )

      !-------------------------------------------------------------------
      ! **** Advance to the aircraft data ****
      !-------------------------------------------------------------------

      use dshr_methods_mod , only : dshr_fldbun_getfldptr
      use dshr_strdata_mod , only : shr_strdata_init_from_inline, shr_strdata_advance
      use cam_esmf_mod,      only : model_mesh, model_clock
      use physics_types,     only : physics_state
      use ppgrid,            only : begchunk, endchunk, pcols, pver, pverp
      use string_utils,      only : to_lower, GLC
      use cam_history,       only : outfld
      use physconst,         only : mwdry ! molecular weight dry air ~ kg/kmole
      use physconst,         only : boltz ! J/K/molecule
      use phys_grid,         only : get_wght_all_p, get_ncols_p
      use physics_buffer,    only : physics_buffer_desc, pbuf_get_field
      use physics_buffer,    only : pbuf_get_chunk, pbuf_get_index
      use time_manager,      only : get_curr_date

      ! Arguments
      type(physics_state), intent(in)    :: state(begchunk:endchunk)
      type(physics_buffer_desc), pointer :: pbuf2d(:,:)

      ! Local variables
      integer               :: gcell, ind, nf
      integer               :: lchnk, icol, klev, ncol
      integer               :: caseid
      integer               :: year, mon, day, sec
      integer               :: mcdate
      real(r8)              :: to_mmr(pcols,pver)
      real(r8)              :: wght(pcols)
      real(r8), pointer     :: tmpptr(:,:)
      real(r8), pointer     :: data_out(:,:)
      real(r8), pointer     :: dataptr2d(:,:)
      real(r8)              :: datain3d(pcols,pver,begchunk:endchunk)
      real(r8)              :: data_col(pver)
      real(r8)              :: model_z(pverp)
      character(len=cs)     :: units
      integer               :: rc
      logical               :: first_time = .true.
      type(physics_buffer_desc), pointer :: pbuf_chnk(:)
      real(r8), parameter :: m2km  = 1.e-3_r8
      character(len=*), parameter :: subname = 'aircraft_emit_adv'
      !------------------------------------------------------------------

      call t_startf('All_aircraft_emit_adv')

      !------------------------------------------------------------------
      ! The stream initialization must be called after the cam initailization
      !------------------------------------------------------------------
      n_aero_loop: do nf = 1,N_AERO
         unset_file: if (trim(forcing(nf)%datafile) /= 'unset') then

            first_call: if (first_time) then
               ! Initialize forcing%sdat
               call shr_strdata_init_from_inline(forcing(nf)%sdat,    &
                    my_task             = iam,                        &
                    logunit             = iulog,                      &
                    compname            = 'ATM',                      &
                    model_clock         = model_clock,                &
                    model_mesh          = model_mesh,                 &
                    stream_meshfile     = trim(forcing(nf)%meshfile), &
                    stream_filenames    = (/forcing(nf)%datafile/),   &
                    stream_yearFirst    = forcing(nf)%year_first,     &
                    stream_yearLast     = forcing(nf)%year_last,      &
                    stream_yearAlign    = forcing(nf)%year_align,     &
                    stream_fldlistFile  = (/forcing(nf)%fldname/),    &
                    stream_fldListModel = (/forcing(nf)%fldname/),    &
                    stream_lev_dimname  = 'altitude',                 &
                    stream_mapalgo      = trim(forcing(nf)%mapalgo),  &
                    stream_offset       = 0,                          &
                    stream_taxmode      = trim(forcing(nf)%taxmode),  &
                    stream_dtlimit      = 1.0e30_r8,                  &
                    stream_tintalgo     = trim(forcing(nf)%tintalgo), &
                    stream_name         = 'Aircraft forcing data ',   &
                    rc                  = rc)
               call chkrc(rc,__LINE__,u_FILE_u)

               first_time = .false.
            end if first_call

            !-------------------------------------------------------------------
            !  For each field, interpolate data in time and to model horizontal grid
            !-------------------------------------------------------------------

            ! Extract YMD from model_update_next_time
            call get_curr_date(year, mon, day, sec)
            mcdate = year*10000 + mon*100 + day

            ! Advance sdat streams
            call shr_strdata_advance(forcing(nf)%sdat, ymd=mcdate, tod=sec, logunit=iulog, &
                 istr='aircraft_stream', rc=rc)
            call chkrc(rc,__LINE__,u_FILE_u)

            ! Get pointer to horizontally interpolated data
            call dshr_fldbun_getFldPtr(forcing(nf)%sdat%pstrm(1)%fldbun_model, trim(forcing(nf)%fldname), &
                 fldptr2=dataptr2d, rc=rc)
            call chkrc(rc,__LINE__,u_FILE_u)

            ! Obtain datain on model horizontal grid but the same vertical levels as the forcing dataset
            do klev = 1, forcing(nf)%nlev  !nlev is the number of levels in the forcing data
               gcell = 1
               do lchnk = begchunk,endchunk
                  ncol = get_ncols_p(lchnk)
                  do icol = 1,ncol
                     datain3d(icol,klev,lchnk) = dataptr2d(klev,gcell)
                     gcell = gcell + 1
                  end do
               end do
            end do

            ! Do vertical interpolation - aircraft data is vertically
            ! interpolated to conserve the total column
            do lchnk = begchunk,endchunk
               call pbuf_get_field(pbuf2d, lchnk, forcing(nf)%pbuf_index, data_out)
               ncol = get_ncols_p(lchnk)
               do icol = 1,ncol
                  model_z(1:pverp) = m2km * state(lchnk)%zi(icol,pverp:1:-1)
                  call interpz_conserve( forcing(nf)%nlev, pver, forcing(nf)%altitude_int, model_z, &
                       datain3d(icol,:,lchnk), data_col(:) )
                  data_out(icol,:) = data_col(pver:1:-1)
               end do
            end do

            !-------------------------------------------------------------------
            ! set the tracer fields with the correct units
            !-------------------------------------------------------------------

            ! GLC IS position of last significant character in string.
            units = to_lower(trim(forcing(nf)%fldunits(:GLC(forcing(nf)%fldunits))))
            select case (trim(units))
            case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
               caseid = 1
            case ('kg/kg','mmr')
               caseid = 2
            case ('mol/mol','mole/mole','vmr','fraction')
               caseid = 3
            case ('kg/kg/sec')
               caseid = 4
            case ('kg m-2 s-1')
               caseid = 5
            case ('m/sec' )
               caseid = 6
            case default
               if (masterproc) then
                  write(iulog,*)'aircraft_emit_adv: units = '//trim(units)//' are not recognized'
               end if
               call endrun('aircraft_emit_adv: units are not recognized')
            end select

            !$OMP PARALLEL DO PRIVATE (lchnk, ncol, to_mmr, tmpptr, pbuf_chnk, wght)
            do lchnk = begchunk,endchunk
               ncol = state(lchnk)%ncol

               ! Turn emission data to mixing ratio
               call get_wght_all_p(lchnk, ncol, wght(:ncol))

               if (caseid == 1) then
                  to_mmr(:ncol,:) = (molmass(nf)*1.e6_r8*boltz*state(lchnk)%t(:ncol,:)) &
                                   /(mwdry*state(lchnk)%pmiddry(:ncol,:))
               elseif (caseid == 2) then
                  to_mmr(:ncol,:) = 1._r8
               elseif (caseid == 4) then
                  to_mmr(:ncol,:) = 1.0_r8
               elseif (caseid == 5) then
                  to_mmr(:ncol,:) = 1.0_r8
               elseif (caseid == 6) then
                  to_mmr(:ncol,:) = 1.0_r8
               else
                  to_mmr(:ncol,:) = molmass(nf)/mwdry
               endif

               pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
               call pbuf_get_field(pbuf_chnk, forcing(nf)%pbuf_index, tmpptr)
               tmpptr(:ncol,:) = tmpptr(:ncol,:)*to_mmr(:ncol,:)
               call outfld( forcing(nf)%fldname, tmpptr(:ncol,:), ncol, state(lchnk)%lchnk )
            enddo

         end if unset_file
      end do n_aero_loop

      call t_stopf('All_aircraft_emit_adv')

   end subroutine aircraft_emit_adv

   !=========================================================================
   subroutine interpz_conserve( nsrc, ndst, src_x, dst_x, src, dst)

      ! Arguments
      integer, intent(in)   :: nsrc                  ! dimension source array
      integer, intent(in)   :: ndst                  ! dimension target array
      real(r8), intent(in)  :: src_x(nsrc+1)         ! source coordinates
      real(r8), intent(in)  :: dst_x(ndst+1)         ! target coordinates
      real(r8), intent(in)  :: src(nsrc)             ! source array
      real(r8), intent(out) :: dst(ndst)             ! target array

      ! local variables
      integer  :: i, j
      integer  :: isrc
      real(r8) :: y
      real(r8) :: bot, top
      !---------------------------------------------------------------

      do i = 1, ndst
         if ( (dst_x(i)<src_x(nsrc+1)) .and. (dst_x(i+1)>src_x(1)) ) then
            do isrc = 1,nsrc
               if ( (dst_x(i)-src_x(isrc))*(dst_x(i)-src_x(isrc+1))<=0.0_r8 ) then
                  exit
               end if
            end do

            if ( dst_x(i)<src_x(1) ) isrc = 1

            y = 0.0_r8
            bot = max(dst_x(i),src_x(1))
            top = dst_x(i+1)
            do j = isrc, nsrc
               if ( top>src_x(j+1) ) then
                  y = y+(src_x(j+1)-bot)*src(j)/(src_x(j+1)-src_x(j))
                  bot = src_x(j+1)
               else
                  y = y+(top-bot)*src(j)/(src_x(j+1)-src_x(j))
                  exit
               endif
            enddo
            dst(i) = y
         else
            dst(i) = 0.0_r8
         end if
      end do

      if ( dst_x(1)>src_x(1) ) then
         top = dst_x(1)
         bot = src_x(1)
         y = 0.0_r8
         do j = 1, nsrc
            if ( top>src_x(j+1) ) then
               y = y+(src_x(j+1)-bot)*src(j)/(src_x(j+1)-src_x(j))
               bot = src_x(j+1)
            else
               y = y+(top-bot)*src(j)/(src_x(j+1)-src_x(j))
               exit
            endif
         end do
         dst(1) = dst(1)+y
      end if

   end subroutine interpz_conserve

   !=========================================================================
   subroutine get_aircraft(cnt, spc_name_list_out)

      ! Arguments
      integer,          intent(out) :: cnt
      character(len=*), intent(out) :: spc_name_list_out(:)

      ! Local variables
      integer :: nf
      !------------------------------------------------------------------

      cnt = 0
      spc_name_list_out(:) = ''

      do nf = 1,N_AERO
         if (trim(forcing(nf)%datafile) /= 'unset') then
            cnt = cnt + 1
            spc_name_list_out(nf) = trim(forcing(nf)%fldname)
         end if
      end do

   end subroutine get_aircraft

   !=========================================================================
   subroutine get_vertical_dimension( fid, dname, dsize, data )

      use pio, only : file_desc_t, pio_seterrorhandling
      use pio, only : pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, pio_get_var
      use pio, only : PIO_BCAST_ERROR, PIO_NOERR

      ! Arguments
      type(file_desc_t), intent(inout) :: fid
      character(*),      intent(in)    :: dname
      integer,           intent(out)   :: dsize
      real(r8),          pointer       :: data(:)

      ! Local variables
      integer :: vid, ierr, id
      integer :: err_handling
      character(len=*), parameter :: subname = 'get_vertical_dimension'
      !------------------------------------------------------------------

      call pio_seterrorhandling(fid, PIO_BCAST_ERROR, oldmethod=err_handling)
      ierr = pio_inq_dimid( fid, dname, id )
      if ( ierr == PIO_NOERR ) then
         ierr = pio_inq_dimlen( fid, id, dsize )
         if (ierr /= PIO_NOERR) then
            call endrun(trim(subname)//': failed on pio_inq_dimid')
         end if
         allocate( data(dsize), stat=ierr )
         if ( ierr /= 0 ) then
            call endrun(trim(subname)//': failed to allocate data array')
         end if
         ierr = pio_inq_varid( fid, dname, vid )
         if (ierr /= PIO_NOERR) then
            call endrun(trim(subname)//': failed on pio_inq_varid')
         end if
         ierr = pio_get_var( fid, vid, data )
         if (ierr /= PIO_NOERR) then
            call endrun(trim(subname)//': failed on pio_get_var')
         end if
      endif
      call pio_seterrorhandling(fid, err_handling)

   end subroutine get_vertical_dimension

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

end module aircraft_emit
