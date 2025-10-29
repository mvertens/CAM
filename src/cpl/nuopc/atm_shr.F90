module atm_shr

  use ESMF, only : ESMF_Mesh, ESMF_Clock

  implicit none
  public

  type(ESMF_Mesh)  :: model_mesh     ! model mesh
  type(ESMF_Clock) :: model_clock    ! model clock

end module atm_shr
