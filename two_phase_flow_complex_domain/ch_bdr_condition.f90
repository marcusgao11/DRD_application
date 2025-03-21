
! boundary condition for phi
function ch_bdr_condition(angle, phi_val, pi, phi_mode)
  implicit none
  real, intent(in) :: angle, phi_val, pi
  integer, intent(in):: phi_mode
  real :: ch_bdr_condition
  real cf

  if (phi_mode == 1) then
    ! cf = sqrt(2.0)*cos(angle*pi/180.0)/3.0*pi/2.
    ! ch_bdr_condition = cf*cos(phi_val*pi/2)

    cf = cos(angle*pi/180.0)/sqrt(2.0)
    ch_bdr_condition = cf*(1.0 - phi_val)*phi_val
    
  elseif (phi_mode == 2) then 
    ! cf = sqrt(2.0)*cos(angle*pi/180.0)/3.0*pi/2.
    ! ch_bdr_condition = cf*cos(phi_val*pi/2)

    cf = cos(angle*pi/180.0)/sqrt(2.0)
    ch_bdr_condition = cf*(1.0 - phi_val**2.0)
  endif



end function ch_bdr_condition

