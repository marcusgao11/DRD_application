  subroutine PreSolve_2(plan,plan_back,psi)
  use global
  implicit none
  integer*8 plan,plan_back
  real psi(-1:nr,-1:nz)

  !real uin(-1:nr+1,-1:nz+1)
  !real w1(-1:nr+1,-1:nz+1),w2(-1:nr+1,-1:nz+1)
  !real pnb(0:nr),pnt(0:nr),pnl(0:nz),pnr(0:nz)

  integer i,k
  real uin(0:nr-1,0:nz-1)


  uin = 0.0
  do k=0,nz-1
     do i=0,nr-1
      uin(i,k) = ( u(i+1,k)-u(i,k) )/dr + ( v(i,k+1)-v(i,k) )/dz
     enddo
  enddo


!TODO: check which is faster
!TODO: solve the linear system
   !call DCT(nr,nz,uin,psi,Pre_Lcc,wsavexc,wsaveyc)
   !psi(1:nr, 1:nz) = uin

   call dfftw_execute_r2r(plan, uin, uin)

   do k=0,nz-1
     do i=0,nr-1
        uin(i,k) = uin(i,k)/Pre_Lcc(i,k)
     enddo
   enddo

   call dfftw_execute_r2r(plan_back, uin, uin)

   do k=0,nz-1
       do i=0,nr-1
          psi(i,k) = uin(i,k)/(4.0*nr*nz)
       enddo
    enddo

!        call output(nr,nz,r,z,w1,w2,p,(w1+w2)*3*re/2/dt,0,'testing')
!        print *,'output test'
!	pause

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


! update velocity if needed
    ! TODO: update velocity
    do k=0,nz-1 ! Correct the u-velocity
       do i=0,nr
        u(i,k) = u(i,k) - (psi(i,k)-psi(i-1,k))/dr;
       enddo
    enddo

    do k=1,nz-1 ! Correct the v-velocity
      do i=0,nr-1
        v(i,k) = v(i,k) - (psi(i,k)-psi(i,k-1))/dz;
      enddo
    enddo

    ! bdr
    do k=0,nz-1
       u(-1,k) = 2.0*uL(k) - u(1,k)
    enddo

    do k=0,nz-1
       u(nr+1,k) = 2.0*uR(k) - u(nr-1,k)
    enddo

    do i=-1,nr+1
       u(i,-1) = u(i,nz-1)
       u(i,nz) = u(i,0)
    enddo

    !TODO: need extend v to [0, nz-1] for periodic bdrs
    ! Here, we don't extend boundary, and instead approximate the periodic
    ! bdrs by interpolate values at boundary. in this case, we don't need change
    ! the size of v in iteration process
    do i=-1,nr
       v(i,0) = ( v(i,nz-1) + v(i, 1) )/2.0
       v(i,nz) = ( v(i,nz-1) + v(i,1) ) /2.0
    enddo

    do k=0,nz
       v(-1,k) = v(0,k)
       v(nr,k) = v(nr-1,k)
    enddo


! verify


  end subroutine

