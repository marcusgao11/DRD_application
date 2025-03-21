!
! this function is used to interpolate related parameters,
! such as eta, density, slip length
subroutine interpolate_coef() !indi)
use global
use particles
implicit none
!real indi(-1:nr,-1:nz)
integer i,k,i_num
real u_aver,v_aver,sum_val, dd
real eta_tmp(-1:nr,-1:nz)
real,external::dist

!  do j=-1,nz
!      do i=-1,nr
!         if(indi(i,j) <= -1.0) then
!            rho(i,j)=1.0
!            eta(i,j)=1.0
!         elseif(indi(i,j) >= 1.0) then
!            rho(i,j)=lambda_rho
!            eta(i,j)=lambda_eta
!         else
!            rho(i,j) = 0.5*(1.0-indi(i,j)) + 0.5*lambda_rho*(1.0+indi(i,j))
!            eta(i,j) = 0.5*(1.0-indi(i,j)) + 0.5*lambda_eta*(1.0+indi(i,j))
!         endif
!      enddo
!  enddo


!    do i_num = 1,fpd_num
!      u_aver = 0.
!      v_aver = 0.
!      sum_val = 0.0
!
!      do k=0,nz-1
!        do i=0,nr-1
!          dd = sqrt((r(i,k)-fpd_x_pos(i_num))**2 + (z(i,k)-fpd_y_pos(i_num))**2) - fpd_radius(i_num)
!          dd = (tanh( -1.0*dd /sqrt(2.0)/I_thickness ) + 1.0)/2.0
!
!          u_aver = u_aver + (u(i,k)+u(i,k+1))/2.0*dd
!          v_aver = v_aver + (v(i,k)+v(i,k+1))/2.0*dd
!          sum_val = sum_val + dd
!        enddo
!      enddo
!      u_aver = u_aver / sum_val
!      v_aver = v_aver / sum_val
!
!      ! update position
!      fpd_x_pos(i_num) = fpd_x_pos(i_num) + u_aver*dt
!      fpd_y_pos(i_num) = fpd_y_pos(i_num) + v_aver*dt
!      write(*,*) "i_num: ", i_num, ", u_aver: ", u_aver, ", v_aver: ", v_aver
!    enddo

    ! update phi
    !    call random_seed()
    do k=0,nz-1
      do i=0,nr-1
        fpd_ph(i,k) = tanh( -1.0*dist(r(i,k), z(i,k)) /sqrt(2.0)/I_thickness )
         !dist(r(i,k), z(i,k)) !tanh( -1.0*dist(r(i,k), z(i,k)) /sqrt(2.0)/I_thickness )
      enddo
    enddo

    ! re-scale to [0.0,1.0]
    fpd_ph = (fpd_ph + 1.0)/2.0


    ! incorporate boundary phi
    fpd_ph = fpd_ph + bdr_ph

    ! re-scale to [0.0,1.0] need: since already done in dist.f90
    !ph = (ph + 1.0)/2.0
    do k=0,nz-1
      fpd_ph(-1,k) = fpd_ph(0,k)
      fpd_ph(nr,k) = fpd_ph(nr-1,k)
    enddo

    do i=-1,nr
      fpd_ph(i,-1) = fpd_ph(i,0)
      fpd_ph(i,nz) = fpd_ph(i,nz-1)
    end do



  ! Calculte rho and eta
   do k=-1,nz
        do i=-1,nr
            rho(i,k) = 1.0 + (lambda_rho-1.0)*fpd_ph(i,k)
            eta(i,k) = 1.0 + (lambda_eta-1.0)*fpd_ph(i,k)
        enddo
  enddo
  
  ! !TODO: two-phase viscosity/density ratio
  ! do k=-1,nz
  !   do i=-1,nr
  !      if(ph(i,k) <= -1.0) then
  !         rho(i,k)=1.0
  !         eta_tmp(i,k) = 1.0
  !      elseif(ph(i,k) >= 1.0) then
  !         rho(i,k)=lambda_rho
  !         eta_tmp(i,k) = lambda_eta
  !       else
  !        rho(i,k) = 1.0 + (lambda_rho-1.0)*ph(i,k)
  !        eta_tmp(i,k) = 1.0 + (lambda_eta-1.0)*ph(i,k)
  !      endif
  !   enddo
  ! enddo
  ! eta = max(eta, eta_tmp)


! rho = 1.0 + (lambda_rho-1.0)*ph

! eta = 1.0 + (lambda_eta-1.0)*ph


end subroutine interpolate_coef
