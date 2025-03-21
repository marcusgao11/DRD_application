


function ch_bulk_energy(phi, psi, phi_mode)
  implicit none
  real, intent(in) :: phi, psi
  integer, intent(in) :: phi_mode
  real :: ch_bulk_energy
  real aa

  if (phi_mode == 1) then
    !如果需要使用注释掉的公式，取消下一行的注释
    !ch_bulk_energy = 0.5 * phi * (1.0 - phi) * (1.0 - 2.0 * phi)

    ! 1 - \psi
    !如果需要使用注释掉的公式，取消下一行的注释
    ch_bulk_energy = 0.5 * phi * (1.0 - phi) * (1.0 - 2.0 * phi) * (1-psi)

    ! !localize
    ! aa = phi**2.0 * (1.0 - phi) **2.0
    ! aa = 1.0 - sqrt( aa / (aa + 1.0E-2) ) * psi;
    ! ch_bulk_energy = 0.5 * phi * (1.0 - phi) * (1.0 - 2.0 * phi) * aa  !- phi * (1 - phi - psi) * psi ! yibao li's work

    ! ! yibao li's work
    ! ch_bulk_energy = phi * (phi -0.5)*(phi-1.0)   &
    !        & - phi*(1.0-phi-psi)*psi
  elseif (phi_mode == 2) then 

      ! \phi = -1, 1
      ch_bulk_energy = (phi**3.0 - phi) * (1-psi)

    ! !localize
    ! aa = (1.0 - phi**2) **2.0
    ! aa = 1.0 - sqrt( aa / (aa + 1.0E-2) ) * psi;
    ! ch_bulk_energy = (phi**3.0 - phi) * aa  !- phi * (1 - phi - psi) * psi ! yibao li's work

  endif

end function ch_bulk_energy

