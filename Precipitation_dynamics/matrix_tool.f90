

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  transform two matrices to a single vector in natrual sense
  subroutine Atob_uv(nr,nz,f1,f2,b)
  implicit none
  integer nr,nz
  real f1(-1:nr+1,-1:nz), f2(-1:nr,0:nz)
  real b((nr+1)*nz + nr*(nz-1))
  integer i,j,mn

    do j=0,nz-1
       do i=0,nr
          b(i+1+(nr+1)*j) = f1(i,j)
       enddo
    enddo

    mn = (nr+1)*nz

    do j=1,nz-1
       do i=0,nr-1
          b(i+1+nr*(j-1) + mn) = f2(i,j)
       enddo
    enddo

  end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  transform two matrices to a single vector in natrual sense
  subroutine btoA_uv(nr,nz,f1,f2,b)
  implicit none
  integer nr,nz
  real f1(-1:nr+1,-1:nz), f2(-1:nr,0:nz)
  real b((nr+1)*nz + nr*(nz-1))
  integer i,j,mn

    do j=0,nz-1
       do i=0,nr
          f1(i,j) = b(i+1+(nr+1)*j)
       enddo
    enddo

    mn = (nr+1)*nz

    do j=1,nz-1
       do i=0,nr-1
          f2(i,j) = b(i+1+nr*(j-1) + mn)
       enddo
    enddo


    !!reflective boundary condition to enforce zero velocity
    ! !zero boundary condition
!    do i=0,nr
!       f1(i,-1) = -f1(i,0)
!       f1(i,nz) = -f1(i,nz-1)
!    enddo
    do i=0,nr
       f1(i,-1) = 0.0
       f1(i,nz) = 0.0
    enddo

    !zero boundary condition
    do i=0,nr-1
       f2(i,0) = 0
       f2(i,nz) = 0
    enddo
  end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  transform two matrices to a single vector in natrual sense
  subroutine Atob_ph_mu(nr,nz,f1,f2,b)
  implicit none
  integer nr,nz
  real f1(-1:nr,-1:nz), f2(-1:nr,-1:nz)
  real b(nr*nz*2)
  integer i,j,mn

  mn = nr*nz

  do j=0,nz-1
     do i=0,nr-1
        b(i+1+nr*j) = f1(i,j)
        b(i+1+nr*j + mn) = f2(i,j)
     enddo
  enddo

  end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  transform two matrices to a single vector in natrual sense
  subroutine btoA_ph_mu(nr,nz,f1,f2,b)
  implicit none
  integer nr,nz
  real f1(-1:nr,-1:nz), f2(-1:nr,-1:nz)
  real b(nr*nz*2)
  integer i,j,mn

  mn = nr*nz

  do j=0,nz-1
     do i=0,nr-1
        f1(i,j) = b(i+1+nr*j)
        f2(i,j) = b(i+1+nr*j + mn)
     enddo
  enddo

  do i=0,nr-1
     f1(i,-1)=f1(i,0)
     f1(i,nz)=f1(i,nz-1)

     f2(i,-1)=f2(i,0)
     f2(i,nz)=f2(i,nz-1)
  enddo

  do j=-1,nz
     f1(-1,j)=f1(0,j)
     f1(nr,j)=f1(nr-1,j)

     f2(-1,j)=f2(0,j)
     f2(nr,j)=f2(nr-1,j)
  enddo

  end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  transform matrix to a single vector in natrual sense
! the transformed data range is (0, m)x(0, n), and input data range is (-1, m+1)x(-1, n+1),
  subroutine Atob(m,n,f1,b)
  implicit none
  integer m,n
  real f1(-1:m+1,-1:n+1)
  real b((m+1)*(n+1))
  integer i,j

  do j=0,n
     do i=0,m
        b(i+1+(m+1)*j)=f1(i,j)
     enddo
  enddo

  end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! decompose a single vector to two matrices, inverse process of the above subroutine Atob
  subroutine btoA(m,n,f1,b)
  implicit none
  integer m,n
  real f1(-1:m+1,-1:n+1)
  real b((m+1)*(n+1))
  integer i,j

  do j=0,n
     do i=0,m
        f1(i,j)=b(i+1+(m+1)*j)
     enddo
  enddo

  end subroutine


!  subroutine Atob(nr,nz,f1,b)
!  implicit none
!  integer nr,nz
!  real f1(-1:nr+1,-1:nz+1)
!  real b((nr+1)*(nz+1))
!  integer i,j
!
!  do j=0,nz
!     do i=0,nr
!        b(i+1+(nr+1)*j)=f1(i,j)
!     enddo
!  enddo
!
!  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! decompose a single vector to two matrices, inverse process of the above subroutine Atob
!  subroutine btoA(nr,nz,f1,b)
!  implicit none
!  integer nr,nz
!  real f1(-1:nr+1,-1:nz+1)
!  real b((nr+1)*(nz+1))
!  integer i,j
!
!  do j=0,nz
!     do i=0,nr
!        f1(i,j)=b(i+1+(nr+1)*j)
!     enddo
!  enddo
!
!  end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  transform two matrices to a single vector in natrual sense
!
!  subroutine Atob_2(nr,nz,f1,b)
!  implicit none
!  integer nr,nz
!  real f1(-1:nr+1,-1:nz+1)
!  real b((nr+1)*(nz+1))
!  integer i,j
!
!  do j=0,nz
!     do i=0,nr
!        b(i+1+(nr+1)*j)=f1(i,j)
!     enddo
!  enddo
!
!  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! decompose a single vector to two matrices, inverse process of the above subroutine Atob
!
!  subroutine btoA_2(nr,nz,f1,b)
!  implicit none
!  integer nr,nz
!  real f1(-1:nr+1,-1:nz+1)
!  real b((nr+1)*(nz+1))
!  integer i,j
!
!  do j=0,nz
!     do i=0,nr
!        f1(i,j)=b(i+1+(nr+1)*j)
!     enddo
!  enddo
!
!  end subroutine

! finite difference operators

! laplace periodic in x, homogeneous in z
  subroutine  laplace_pn(a,b,nx,ny,dx2,dy2)
  implicit none
  !real, intent(in):: a(-1:nx,-1:ny)
  real:: a(-1:nx,-1:ny),b(-1:nx,-1:ny)
  integer,intent(in):: nx,ny
  real,intent(in)::dx2,dy2
  integer i,j

  do i=0,nx-1
     a(i,-1) = a(i,0)
     a(i,ny) = a(i,ny-1)
  enddo

  do j=-1,ny
     a(-1,j) = a(nx-1,j)
     a(nx,j) = a(0,j)
  enddo

  do j=0,ny-1
     do i=0,nx-1
        b(i,j) = (a(i+1,j)+a(i-1,j)-2*a(i,j))/dx2 + (a(i,j+1)+a(i,j-1)-2*a(i,j))/dy2
      enddo
  enddo

  end subroutine laplace_pn


  subroutine  laplace_nn(a,b,nx,ny,dx2,dy2)
  implicit none
  !real, intent(in):: a(-1:nx,-1:ny)
  real:: a(-1:nx,-1:ny),b(-1:nx,-1:ny)
  integer,intent(in):: nx,ny
  real,intent(in)::dx2,dy2
  integer i,j

  do i=0,nx-1
     a(i,-1) = a(i,0)
     a(i,ny) = a(i,ny-1)
  enddo

  do j=-1,ny
     a(-1,j) = a(0,j)
     a(nx,j) = a(nx-1,j)
  enddo

  do j=0,ny-1
     do i=0,nx-1
        b(i,j) = (a(i+1,j)+a(i-1,j)-2*a(i,j))/dx2 + (a(i,j+1)+a(i,j-1)-2*a(i,j))/dy2
      enddo
  enddo

  end subroutine laplace_nn


! laplace homogeneous in x, direchlet in z
  subroutine  laplace_nd(a,b,nx,ny,dx2,dy2, c_bot, c_top)
  implicit none
  !real, intent(in):: a(-1:nx,-1:ny)
  real:: a(-1:nx,-1:ny),b(-1:nx,-1:ny)
  integer,intent(in):: nx,ny
  real,intent(in)::dx2,dy2
  real c_bot, c_top
  integer i,j

    do i=0,nx-1
       a(i,-1) = 2.0*c_bot - a(i,0)
       a(i,ny) = 2.0*c_top - a(i,ny-1)
    enddo

    do j=-1,ny
       a(-1,j) = a(0,j)
       a(nx,j) = a(nx-1,j)
    enddo

    do j=0,ny-1
       do i=0,nx-1
          b(i,j) = (a(i+1,j)+a(i-1,j)-2*a(i,j))/dx2 + (a(i,j+1)+a(i,j-1)-2*a(i,j))/dy2
        enddo
    enddo

  end subroutine laplace_nd


  subroutine  laplace_np(a,b,nx,ny,dx2,dy2)
  implicit none
  !real, intent(in):: a(-1:nx,-1:ny)
  real:: a(-1:nx,-1:ny),b(-1:nx,-1:ny)
  integer,intent(in):: nx,ny
  real,intent(in)::dx2,dy2
  integer i,j

  do i=0,nx-1
     a(i,-1) = a(i,ny-1)
     a(i,ny) = a(i,0)
  enddo

  do j=-1,ny
     a(-1,j) = a(0,j)
     a(nx,j) = a(nx-1,j)
  enddo

  do j=0,ny-1
     do i=0,nx-1
        b(i,j) = (a(i+1,j)+a(i-1,j)-2*a(i,j))/dx2 + (a(i,j+1)+a(i,j-1)-2*a(i,j))/dy2
      enddo
  enddo

  end subroutine laplace_np


! laplace homogeneous in x, direchlet in z, with variable mobility
  subroutine  laplace_nd_var_coef(a,b,nx,ny,dx2,dy2, mo)
  implicit none
  !real, intent(in):: a(-1:nx,-1:ny)
  real:: a(-1:nx,-1:ny),b(-1:nx,-1:ny), mo(-1:nx,-1:ny)
  integer,intent(in):: nx,ny
  real,intent(in)::dx2,dy2
  integer i,j
  real moL, moR

    do i=0,nx-1
       a(i,-1) = - a(i,0)
       a(i,ny) = - a(i,ny-1)
    enddo

    do j=-1,ny
       a(-1,j) = a(0,j)
       a(nx,j) = a(nx-1,j)
    enddo

    do j=0,ny-1
       do i=0,nx-1
          moL = (mo(i-1,j) + mo(i,j))/2.0
          moR = (mo(i,j) + mo(i+1,j))/2.0
          b(i,j) = ( moL*a(i-1,j) + moR*a(i+1,j) - (moL+moR)*a(i,j) )/dx2

          moL = (mo(i,j-1) + mo(i,j))/2.0
          moR = (mo(i,j) + mo(i,j+1))/2.0

          b(i,j) = b(i,j) + ( moL*a(i,j-1) + moR*a(i,j+1) - (moL+moR)*a(i,j) )/dy2
        enddo
    enddo

  end subroutine laplace_nd_var_coef



! laplace direchlet at left, homo at right and periodic at top and bot
  subroutine  laplace_dp_var_coef(a,b,nx,ny,dx2,dy2, mo)
  implicit none
  !real, intent(in):: a(-1:nx,-1:ny)
  real:: a(-1:nx,-1:ny),b(-1:nx,-1:ny), mo(-1:nx,-1:ny)
  integer,intent(in):: nx,ny
  real,intent(in)::dx2,dy2
  integer i,j
  real moL, moR

    ! periodic
    do i=0,nx-1
       a(i,-1) = a(i,ny-1)
       a(i,ny) = a(i,0)
    enddo

    ! Direchelt at left
    do j=-1,ny
       a(-1,j) = -a(0,j)
    enddo

    ! Homo at right
    do j=-1,ny
       a(nx,j) = a(nx-1,j)
    enddo


    do j=0,ny-1
       do i=0,nx-1
          moL = (mo(i-1,j) + mo(i,j))/2.0
          moR = (mo(i,j) + mo(i+1,j))/2.0
          b(i,j) = ( moL*a(i-1,j) + moR*a(i+1,j) - (moL+moR)*a(i,j) )/dx2

          moL = (mo(i,j-1) + mo(i,j))/2.0
          moR = (mo(i,j) + mo(i,j+1))/2.0

          b(i,j) = b(i,j) + ( moL*a(i,j-1) + moR*a(i,j+1) - (moL+moR)*a(i,j) )/dy2
        enddo
    enddo

  end subroutine laplace_dp_var_coef



