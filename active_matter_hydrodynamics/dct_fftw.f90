  ! rhs: will be changed
  ! fftw for solving poisson equation
  subroutine dct_fftw(nr,nz,plan,plan_back, rhs, psi, pre_cof, coef)
  implicit none
  integer*8 plan,plan_back
  integer nr,nz
  real psi(-1:nr,-1:nz), pre_cof(0:nr,0:nz)
  real rhs(0:nr-1,0:nz-1), rhs_wk(0:nr-1,0:nz-1)
  real coef

  integer i,k

  rhs_wk = rhs

! cannot adjust at will
!  ! shift to zero
!  rhs = rhs - SUM(SUM(rhs,1))/nr/nz


   call dfftw_execute_r2r(plan, rhs_wk, rhs_wk)

   do k=0,nz-1
     do i=0,nr-1
        rhs_wk(i,k) = rhs_wk(i,k)/pre_cof(i,k)
     enddo
   enddo

   !rhs_wk(0,0) = 0.0

   call dfftw_execute_r2r(plan_back, rhs_wk, rhs_wk)

   do k=0,nz-1
       do i=0,nr-1
          psi(i,k) = rhs_wk(i,k)/coef
       enddo
    enddo


!!!!!!!!!!!!!!!!!  ghost points for pressure
! here assume Neumann boundary condition
! and reflective point only over one grid distance, not two grids
! since we use staggered grid here
    do k=0,nz-1
     psi(-1,k)=psi(0,k)
     psi(nr,k)=psi(nr-1,k)
    enddo

    do i=-1,nr
     psi(i,-1) = psi(i,0)
     psi(i,nz) = psi(i,nz-1)
    end do


  end subroutine


