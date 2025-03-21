

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


    !reflective boundary condition to enforce zero velocity
    do i=0,nr
       f1(i,-1) = -f1(i,0)
       f1(i,nz) = -f1(i,nz-1)
    enddo
    ! Nuemann bcd
    do j=-1,nz
       f1(-1,j) = f1(1,j)
       f1(nr+1,j) = f1(nr-1,j)
    enddo

    !zero boundary condition
    do i=0,nr-1
       f2(i,0) = 0
       f2(i,nz) = 0
    enddo
    ! Nuemann bcd, half grid
    do j=0,nz
       f2(-1,j) = f2(0,j)
       f2(nr,j) = f2(nr-1,j)
    enddo
  end



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  transform two matrices to a single vector in natrual sense
  subroutine Atob_uv_periodic(nr,nz,f1,f2,b)
   implicit none
   integer nr,nz
   real f1(-1:nr,-1:nz), f2(-1:nr,-1:nz)
   real b(nr*nz + nr*nz)
   integer i,j,mn
 
     do j=0,nz-1
        do i=0,nr-1
           b(i+1+nr*j) = f1(i,j)
        enddo
     enddo
 
     mn = nr*nz
 
     do j=0,nz-1
        do i=0,nr-1
           b(i+1+nr*j + mn) = f2(i,j)
        enddo
     enddo
 
   end subroutine


   subroutine btoA_uv_periodic(nr,nz,f1,f2,b)
      implicit none
      integer nr,nz
      real f1(-1:nr,-1:nz), f2(-1:nr,-1:nz)
      real b(nr*nz + nr*nz)
      integer i,j,mn
    
        do j=0,nz-1
           do i=0,nr-1
              f1(i,j) = b(i+1+nr*j)
           enddo
        enddo
    
        mn = nr*nz
    
        do j=0,nz-1
           do i=0,nr-1
              f2(i,j) = b(i+1+nr*j + mn)
           enddo
        enddo
    
        !periodic boundary condition
        f1(nr,:) = f1(0, :)
        f1(-1,:) = f1(nr-1, :)

        f1(:, nz) = f1(:, 0)
        f1(:, -1) = f1(:, nz-1)

        f2(-1,:) = f2(nr-1, :)
        f2(nr, :) = f2(0, :)

        f2(:, nz) = f2(:, 0)
        f2(:, -1) = f2(:, nz-1)

   end subroutine !btoA_uv_periodic(nr,nz,f1,f2,b)


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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  transform matrix to a single vector in natrual sense
! the transformed data range is (0, m)x(0, n), and input data range is (-1, m+1)x(-1, n+1),
  subroutine Atob_Periodic(m,n,f1,b)
   implicit none
   integer m,n
   real f1(-1:m,-1:n)
   real b(m*n)
   integer i,j
 
   do j=0,n-1
      do i=0,m-1
         b(i+1+m*j)=f1(i,j)
      enddo
   enddo

   end subroutine
 
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! decompose a single vector to two matrices, inverse process of the above subroutine Atob
   subroutine btoA_Periodic(m,n,f1,b)
   implicit none
   integer m,n
   real f1(-1:m,-1:n)
   real b(m*n)
   integer i,j
 
   do j=0,n-1
      do i=0,m-1
         f1(i,j)=b(i+1+m*j)
      enddo
   enddo

   ! periodic bdrs
   f1(-1,:) = f1(m-1,:)
   f1(m,:) = f1(0,:)

   f1(:,-1) = f1(:,n-1)
   f1(:,n) = f1(:,0)

   ! ! Nuem bdrs
   ! f1(-1,:) = f1(0,:)
   ! f1(m,:) = f1(m-1,:)

   ! f1(:,-1) = f1(:,0)
   ! f1(:,n) = f1(:,n-1)

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


