module cam_esmf_mod

  use shr_kind_mod      , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use ESMF              , only : ESMF_Mesh, ESMF_Clock
  use ESMF              , only : ESMF_VM, ESMF_VMAllreduce, ESMF_VMGetCurrent
  use ESMF              , only : ESMF_SUCCESS, ESMF_REDUCE_SUM
  use shr_sys_mod       , only : shr_sys_abort
  use cam_abortutils    , only : endrun
  use nuopc_shr_methods , only : chkerr
  use error_messages    , only : alloc_err
  use cam_logfile       , only : iulog

  implicit none
  private

  public :: cam_esmf_set_clock
  public :: cam_esmf_set_mesh
  public :: cam_esmf_set_areas
  public :: cam_esmf_global_sum

  private :: cam_set_mesh_for_single_column

  type(ESMF_Mesh) , public, protected :: model_mesh     ! model mesh
  type(ESMF_Clock), public, protected :: model_clock    ! model clock

  real(r8), allocatable, public, protected :: model_areas(:)
  real(r8), allocatable, public, protected :: mesh_areas(:)

  logical :: model_clock_initialized = .false.
  logical :: model_mesh_initialized = .false.

  character(*), parameter :: u_FILE_u = &
       __FILE__

!=====================================================================
contains
!=====================================================================

   subroutine cam_esmf_set_clock(clock_in, rc)
      use ESMF, only : ESMF_Clock

      type(ESMF_Clock), intent(in) :: clock_in

      rc = ESMF_SUCCESS

      model_clock = ESMF_ClockCreate(clock_in, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      if (model_clock_initialized) then
         call shr_sys_abort('initialize_model_clock: model clock already initialized')
      else
         model_clock = model_clock_in
         model_clock_initialized = .true.
      end if

   end subroutine cam_esmf_set_clock

   !=====================================================================
   subroutine cam_esmf_set_mesh(model_mesh_in)
      type(ESMF_Mesh) , intent(in) :: model_mesh_in

      if (model_mesh_initialized) then
         call shr_sys_abort('initialize_model_mesh: model mesh already initialized')
      else
         model_mesh  = model_mesh_in
         model_mesh_initialized = .true.
      end if
   end subroutine cam_esmf_set_mesh

   !=====================================================================
   subroutine cam_esmf_set_areas(model_areas_in, mesh_areas_in, rc)

      ! Arguments
      real(r8), intent(in)  :: model_areas_in(:)
      real(r8), intent(in)  :: mesh_areas_in(:)
      integer , intent(out) :: rc

      ! Local variables
      type(ESMF_VM) :: vm
      integer       :: locsize
      integer       :: ng
      integer       :: istat
      real(r8)      :: global_mesh_area(1)
      real(r8)      :: global_model_area(1)
      real(r8)      :: local_mesh_area(1)
      real(r8)      :: local_model_area(1)
      !---------------------------------------

      rc = ESMF_SUCCESS

      locsize = size(model_areas_in)

      allocate(model_areas(locsize), stat=istat)
      call alloc_err(istat,'cam_esmf_set_areas','model_areas',locsize)
      model_areas(:) = model_areas_in(:)

      allocate(mesh_areas(locsize), stat=istat)
      call alloc_err(istat,'cam_esmf_set_areas','mesh_areas',locsize)
      mesh_areas(:) = mesh_areas_in(:)

      ! Compare global sum of model_areas and mesh_areas
      local_model_area(1) = 0._r8
      local_mesh_area(1) = 0._r8
      do ng = 1,locsize
         local_model_area(1) = local_model_area(1) + model_areas(ng)
         local_mesh_area(1) = local_mesh_area(1) + mesh_areas(ng)
      end do

      call ESMF_VMGetCurrent(vm, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_VMAllreduce(vm, senddata=local_model_area, recvdata=global_model_area, &
           count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_VMAllreduce(vm, senddata=local_mesh_area, recvdata=global_mesh_area, &
           count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      write(iulog,'(a,d13.5)') ' global mesh area  = ',global_mesh_area(1)
      write(iulog,'(a,d13.5)') ' global model area = ',global_model_area(1)

   end subroutine cam_esmf_set_areas

   !=====================================================================
   subroutine cam_esmf_global_sum(fldname, flddata, rc)

      ! Arguments
      character(len=*), intent(in)  :: fldname
      real(r8),         intent(in)  :: flddata(:)
      integer ,         intent(out) :: rc

      ! local variables
      type(ESMF_VM) :: vm
      integer       :: ng
      real(r8)      :: local_sum_model(1)
      real(r8)      :: global_sum_model(1)
      real(r8)      :: local_sum_mesh(1)
      real(r8)      :: global_sum_mesh(1)
      !---------------------------------------

      rc = ESMF_SUCCESS

      local_sum_model(1) = 0._r8
      local_sum_mesh(1) = 0._r8
      do ng=1,size(flddata)
         local_sum_model(1) = local_sum_model(1) + flddata(ng) * model_areas(ng)
         local_sum_mesh(1) = local_sum_mesh(1) + flddata(ng) * mesh_areas(ng)
      end do

      call ESMF_VMGetCurrent(vm, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_VMAllreduce(vm, senddata=local_sum_model, recvdata=global_sum_model, &
           count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_VMAllreduce(vm, senddata=local_sum_mesh, recvdata=global_sum_mesh, &
           count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      write(iulog,'(a)') 'Global sum for forcing field '//trim(fldname)
      write(iulog,'(a,d13.5)') ' global sum with model areas = ',global_sum_model(1)
      write(iulog,'(a,d13.5)') ' global sum with mesh areas  = ',global_sum_mesh(1)

   end subroutine cam_esmf_global_sum

end module cam_esmf_mod
