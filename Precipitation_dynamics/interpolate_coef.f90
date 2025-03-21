!
! this function is used to interpolate related parameters,
! such as eta, density, slip length
subroutine interpolate_coef(indi)
use global
implicit none
real indi(-1:nr,-1:nz)
integer i,j

  do j=-1,nz
      do i=-1,nr
         if(indi(i,j) <= -1.0) then
            rho(i,j)=1.0
            eta(i,j)=1.0
         elseif(indi(i,j) >= 1.0) then
            rho(i,j)=lambda_rho
            eta(i,j)=lambda_eta
         else
            rho(i,j) = 0.5*(1.0-indi(i,j)) + 0.5*lambda_rho*(1.0+indi(i,j))
            eta(i,j) = 0.5*(1.0-indi(i,j)) + 0.5*lambda_eta*(1.0+indi(i,j))
         endif
      enddo
  enddo

end subroutine interpolate_coef
