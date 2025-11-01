module cam_esmf_mod

  use ESMF, only : ESMF_Mesh, ESMF_Clock

  implicit none
  public

  type(ESMF_Mesh) , protected :: model_mesh     ! model mesh
  type(ESMF_Clock), protected :: model_clock    ! model clock

contains

   subroutine cam_esmf_set_mesh_and_clock(model_mesh_in, model_clock_in)
      type(ESMF_Mesh) , intent(in) :: model_mesh_in
      type(ESMF_Clock), intent(in) :: model_clock_in

      model_mesh  = model_mesh_in
      model_clock = model_clock_in
   end subroutine cam_esmf_set_mesh_and_clock

end module cam_esmf_mod
