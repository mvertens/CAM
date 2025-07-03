      module mo_lin_matrix
      private
      public :: linmat
      contains
      subroutine linmat01( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)
         mat(1) = -( rxt(9) + rxt(10) + rxt(12) + het_rates(3) )
         mat(5) = -( rxt(11) + het_rates(1) )
         mat(2) = rxt(9) + rxt(10) + .750_r8*rxt(12)
         mat(7) = -( rxt(1) + rxt(8) + het_rates(4) )
         mat(9) = -( het_rates(5) )
         mat(10) = -( het_rates(6) )
         mat(11) = -( het_rates(7) )
         mat(12) = -( het_rates(8) )
         mat(13) = -( het_rates(9) )
         mat(14) = -( het_rates(10) )
         mat(15) = -( het_rates(11) )
         mat(16) = -( het_rates(12) )
         mat(17) = -( het_rates(13) )
         mat(18) = -( het_rates(14) )
         mat(19) = -( het_rates(15) )
         mat(20) = -( het_rates(16) )
         mat(21) = -( het_rates(17) )
         mat(22) = -( het_rates(18) )
         mat(23) = -( het_rates(19) )
         mat(24) = -( het_rates(20) )
         mat(25) = -( het_rates(21) )
         mat(26) = -( het_rates(22) )
         mat(27) = -( het_rates(23) )
         mat(28) = -( het_rates(2) )
         mat(6) = rxt(11)
         mat(29) = -( het_rates(24) )
         mat(30) = -( het_rates(25) )
         mat(31) = -( het_rates(26) )
         mat(3) = .029_r8*rxt(12)
         mat(33) = .150_r8*rxt(13)
         mat(32) = -( het_rates(27) )
         mat(4) = .114_r8*rxt(12)
         mat(34) = .150_r8*rxt(14) + .150_r8*rxt(15)
         mat(36) = .025_r8*rxt(16) + .025_r8*rxt(17) + .025_r8*rxt(18)
         mat(35) = -( rxt(13) + rxt(14) + rxt(15) + het_rates(28) )
         mat(37) = -( rxt(16) + rxt(17) + rxt(18) + het_rates(29) )
         mat(38) = -( rxt(4) + het_rates(30) )
         mat(39) = -( rxt(5) + het_rates(31) )
         mat(40) = -( rxt(6) + het_rates(32) )
         mat(42) = -( rxt(3) + het_rates(34) )
         mat(43) = -( rxt(2) + het_rates(35) )
         mat(41) = 2.000_r8*rxt(6)
         mat(8) = rxt(8)
      end subroutine linmat01
      subroutine linmat( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)
      call linmat01( mat, y, rxt, het_rates )
      end subroutine linmat
      end module mo_lin_matrix
