  subroutine PreSolve_2(psi, coef)
  use global
  implicit none
  real psi(-1:nr,-1:nz)
  real coef

!   interface
!   subroutine  MatrixVector_Lap(f1,f2)
!      use global
!      implicit none
!      real f1(-1:nr,-1:nz),f2(-1:nr,0:nz)
!   end subroutine MatrixVector_Lap
!    end interface

  !real uin(-1:nr+1,-1:nz+1)
  !real w1(-1:nr+1,-1:nz+1),w2(-1:nr+1,-1:nz+1)
  !real pnb(0:nr),pnt(0:nr),pnl(0:nz),pnr(0:nz)

  integer i,k
  real uin(0:nr-1,0:nz-1), uout(0:nr-1,0:nz-1)
  real fO(-1:nr+1,-1:nz+1),fN(-1:nr+1,-1:nz+1)

  real temp_uv(-1:nr,-1:nz)

!   ! for iteration method
!   integer, parameter::mn = nr*nz

!   real rhs(mn),ComputedSolution(mn)
!   real work(mn,7)

!   external MatrixVector_Pre_CG, PSolve_Pre_CG
!   external MatrixVector_Pre_gmres,PSolve_Pre_gmres
!   integer ITER,INFO
!   real RESID



  uin = 0.0
  do k=0,nz-1
     do i=0,nr-1
      uin(i,k) = ( u(i+1,k)-u(i,k) )/dr + ( v(i,k+1)-v(i,k) )/dz
     enddo
  enddo

  !print *,'output SUM SUM', SUM(SUM(uin,1))/nr/nz


!TODO: check which is faster
!TODO: solve the linear system
   !call DCT(nr,nz,uin,psi,Pre_Lcc,wsavexc,wsaveyc)
   !psi(1:nr, 1:nz) = uin


!  ! test dct solver
!  do k=-1,nz
!     do i=-1,nr
!      p(i,k) = cos(20* pi*(i+0.5)*dr/ slen )*cos(10* pi*(k+0.5)*dz/shig)! + 1.0E-3
!     enddo
!  enddo
!
!  uin = 0.0
!  do k=0,nz-1
!     do i=0,nr-1
!      uin(i,k) = ( p(i+1,k) + p(i-1,k) - 2.0*p(i,k) )/dr2 +&
!       & ( p(i,k-1) + p(i,k+1) - 2.0*p(i,k))/dz2
!     enddo
!  enddo

  !call dct_fftw(nr,nz, plan,plan_back, uin, psi, Pre_Lcc, 4.0*nr*nz)


!   ! BICGSTAB: not converge
!   uin_tmp(0:nr-1,0:nz-1) = uin
!   call Atob_Periodic(nr,nz,uin_tmp,rhs)
!   call Atob_Periodic(nr,nz,psi,ComputedSolution) !initial guess

!    ! ITER = 800
!    ! RESID = tol
!    ! call BICGSTAB(nr*nz,rhs,ComputedSolution,work,mn,&
!    !    &ITER,RESID,MatrixVector_Pre_CG,PSolve_Pre_CG,INFO)
!    ! write(*,*) 'Pressure 1 BICGSTAB no slip',ITER
!    ! if(INFO .ne. 0) then
!    !    print *,'Pressure 1 no slip: error',ITER,INFO
!    !    !		stop
!    ! endif

!    call fgmres(nr*nz,150,ComputedSolution,rhs,tol,MatrixVector_Pre_gmres,PSolve_Pre_gmres)

!    call btoA_Periodic(nr,nz,psi,ComputedSolution)


   ! ! test fft solver
   ! do k=-1,nz
   !    do i=-1,nr
   !      p(i,k) = sin(20* pi*(i+0.5)*dr/ slen )*cos(10* pi*(k+0.5)*dz/shig)! + 1.0E-3
   !    enddo
   ! enddo

   ! uin = 0.0
   ! do k=0,nz-1
   !    do i=0,nr-1
   !    uin(i,k) = ( p(i+1,k) + p(i-1,k) - 2.0*p(i,k) )/dr2 +&
   !       & ( p(i,k-1) + p(i,k+1) - 2.0*p(i,k))/dz2
   !    enddo
   ! enddo


  fO(0:nr-1,0:nz-1) = uin
  call FFT(nr,nz,fO,fN,Pre_Lpp,wsavexf,wsaveyf)
  psi(0:nr-1,0:nz-1) = fN(0:nr-1,0:nz-1)

!   call FFT_regular(nr,nz,uin,uout,Pre_Lpp,wsavexf,wsaveyf)
!   psi(0:nr-1,0:nz-1) = uout(0:nr-1,0:nz-1)


  ! periodic bdrs
  psi(-1,:) = psi(nr-1,:);  psi(nr,:) = psi(0,:)
  psi(:,-1) = psi(:,nz-1);  psi(:,nz) = psi(:,0)


!   temp_uv = 0.0
!   temp_uv(0:nr-1,0:nz-1) = uin
! !   ! ! test fft solver
!     !call output(nr,nz,r,z,u,v,p,psi,0.0,'testin2')
!     call output2(nr,nz,r,z,u,v,p, psi, temp_uv, psi, psi, 0.0, 'testin3')
!     print *,'output testin3'

!    call MatrixVector_Lap(psi, temp_uv)
!    ! p = psi
!    ! temp_uv = 0.0
!    ! do k=0,nz-1
!    !    do i=0,nr-1
!    !       temp_uv(i,k) = ( p(i+1,k) + p(i-1,k) - 2.0*p(i,k) )/dr2 +&
!    !       & ( p(i,k-1) + p(i,k+1) - 2.0*p(i,k))/dz2
!    !    enddo
!    ! enddo


!    print *, '======================================='
!     print *, MAXVAL(MAXVAL(ABS(temp_uv(0:nr-1,0:nz-1) - uin(0:nr-1,0:nz-1) ),1))
! !  print *, MINmakeVAL(MINVAL(ABS(p(0:nr-1,0:nz-1) - psi(0:nr-1,0:nz-1) ),1))
! !    pause
!     call output2(nr,nz,r,z,u,v,p, psi, temp_uv, psi, psi, 0.0, 'testin4')
!     print *,'output testin4' 



!  ! for test
!  ph = 0.0
!  do k=0,nz-1
!     do i=0,nr-1
!      ph(i,k) = ( psi(i+1,k) + psi(i-1,k) - 2.0*psi(i,k) )/dr2 +&
!       & ( psi(i,k-1) + psi(i,k+1) - 2.0*psi(i,k))/dz2
!     enddo
!  enddo
!
!  print *, MAXVAL(MAXVAL(ABS(uin(0:nr-1,0:nz-1) - ph(0:nr-1,0:nz-1) ),1))
!  print *, MINVAL(MINVAL(ABS(uin(0:nr-1,0:nz-1) - ph(0:nr-1,0:nz-1) ),1))
!
!  call output(nr,nz,r,z,u,v,psi,ph,0,'testin2')
!  print *,'output testin2'


! ! comment: cannot do this here, since it's not project method
!   ! update velocity if needed
!   ! TODO: update velocity
!   do k=0,nz-1 ! Correct the u-velocity
!      do i=0,nr
!       u(i,k) = u(i,k) - (psi(i,k)-psi(i-1,k))/dr;
!      enddo
!   enddo

!   do k=1,nz-1 ! Correct the v-velocity
!     do i=0,nr-1
!       v(i,k) = v(i,k) - (psi(i,k)-psi(i,k-1))/dz;
!     enddo
!   enddo

!     ! apply boundary conditions
!     ! bdr
!     do k=0,nz-1
!        u(-1,k) = u(nr-1,k)
!        u(nr+1,k) = u(1,k)
!     enddo

!     do i=-1,nr+1
!        u(i,-1) = u(i,nz-1)
!        u(i,nz) = u(i,0)
!     enddo

!     !TODO: need extend v to [0, nz-1] for periodic bdrs
!     ! Here, we don't extend boundary, and instead approximate the periodic
!     ! bdrs by interpolate values at boundary. in this case, we don't need change
!     ! the size of v in iteration process
!     do k=1,nz-1
!        v(-1,k) = v(nr-1,k)
!        v(nr,k) = v(0,k)
!     enddo

!     do i=-1,nr
!        v(i,0) = ( v(i,nz-1) + v(i, 1) )/2.0
!        v(i,nz) = ( v(i,nz-1) + v(i,1) ) /2.0
!     enddo
  

!  ! for test
!  ph = 0.0
!  do k=0,nz-1
!   do i=0,nr-1
!      ph(i,k) = ( u(i+1,k)-u(i,k) )/dr  + ( v(i,k+1)-v(i,k) )/dz
!   enddo
!  enddo

!  print *, "div", MAXVAL(MAXVAL(ABS(ph(0:nr-1,0:nz-1)),1))
!  print *,"div", MINVAL(MINVAL(ABS(ph(0:nr-1,0:nz-1)),1))
!
!  call output(nr,nz,r,z,u,v,psi,ph,0,'testin3')
!  print *,'output testin3'
!	pause

  ! re-scale of pressure
  if (vel_solver_type == 1) then
     psi = psi*coef
     p = psi
  else if (vel_solver_type == 2) then

      ! update velocity if needed
      ! TODO: update velocity
      do k=0,nz-1 ! Correct the u-velocity
         do i=0,nr-1
            u(i,k) = u(i,k) - (psi(i,k)-psi(i-1,k))/dr;
         enddo
      enddo

      do k=0,nz-1 ! Correct the v-velocity
         do i=0,nr-1
            v(i,k) = v(i,k) - (psi(i,k)-psi(i,k-1))/dz;
         enddo
      enddo

         ! apply boundary conditions
         ! bdr
         !periodic boundary condition
         u(nr,:) = u(0, :)
         u(-1,:) = u(nr-1, :)
 
         u(:, nz) = u(:, 0)
         u(:, -1) = u(:, nz-1)
 
         v(-1,:) = v(nr-1, :)
         v(nr, :) = v(0, :)
 
         v(:, nz) = v(:, 0)
         v(:, -1) = v(:, nz-1)

         !TODO: need extend v to [0, nz-1] for periodic bdrs
         ! Here, we don't extend boundary, and instead approximate the periodic
         ! bdrs by interpolate values at boundary. in this case, we don't need change
         ! the size of v in iteration process


         ! do i=-1,nr
         !    v(i,0) = ( v(i,nz-1) + v(i, 1) )/2.0
         !    v(i,nz) = ( v(i,nz-1) + v(i,1) ) /2.0
         ! enddo  

         psi = psi*coef
         p = psi + p

      endif

  end subroutine


  subroutine  MatrixVector_Pre_CG(alpha,xx,beta,yy)
   use global
   implicit none
   real xx(nr*nz), yy(nr*nz)
   real zz(nr*nz)
   real alpha,beta
 
   real f1(-1:nr,-1:nz),f2(-1:nr,0:nz)
 
   integer i,k
 
 
   call btoA_Periodic(nr,nz,f1,xx)

   f2 = 0.0
   do k=0,nz-1
      do i=0,nr-1
         f2(i,k) = ( f1(i+1,k) + f1(i-1,k) - 2.0*f1(i,k) )/dr2 &
         & + ( f1(i,k+1) + f1(i,k-1) - 2.0*f1(i,k) )/dz2
      enddo
   enddo
 
   call Atob_Periodic(nr,nz,f2,zz)
   
   yy = alpha*zz + beta*yy
   
end subroutine
 


subroutine  MatrixVector_Pre_gmres(xx,yy)
   use global
   implicit none
   real xx(nr*nz), yy(nr*nz) 
   real f1(-1:nr,-1:nz),f2(-1:nr,0:nz)
 
   integer i,k
 
   call btoA_Periodic(nr,nz,f1,xx)

   f2 = 0.0
   do k=0,nz-1
      do i=0,nr-1
         f2(i,k) = ( f1(i+1,k) + f1(i-1,k) - 2.0*f1(i,k) )/dr2 &
         & + ( f1(i,k+1) + f1(i,k-1) - 2.0*f1(i,k) )/dz2
      enddo
   enddo
 
   call Atob_Periodic(nr,nz,f2,yy)   
end subroutine


subroutine  MatrixVector_Lap(f1,f2)
   use global
   implicit none
   real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)
 
   integer i,k
 
   f2 = 0.0
   do k=0,nz-1
      do i=0,nr-1
         f2(i,k) = ( f1(i+1,k) + f1(i-1,k) - 2.0*f1(i,k) )/dr2 &
         & + ( f1(i,k+1) + f1(i,k-1) - 2.0*f1(i,k) )/dz2
      enddo
   enddo
 
end subroutine




subroutine PSolve_Pre_gmres(xx, yy)
   use global
   implicit none
   real xx(nr*nz), yy(nr*nz) 
   real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)
 
   call btoA_Periodic(nr,nz,f1,xx)
  
   call dct_fftw(nr,nz, plan,plan_back, f1, f2, Pre_Lcc, 4.0*nr*nz)

   call Atob_Periodic(nr,nz,f2,yy)

!   yy = xx
!   return
 
   end subroutine

 
   subroutine PSolve_Pre_CG(yy,xx)
   use global
   implicit none
   real xx(nr*nz), yy(nr*nz) 
   real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)
 
   call btoA_Periodic(nr,nz,f1,xx)
  
   call dct_fftw(nr,nz, plan,plan_back, f1, f2, Pre_Lcc, 4.0*nr*nz)

   call Atob_Periodic(nr,nz,f2,yy)

!   yy = xx
!   return
 
   end subroutine



   subroutine FFT_regular(nr,nz,fin,fout,cof,wsavex,wsavey)

      implicit none
      integer nr,nz
      real fin(0:nr-1,0:nz-1),fout(0:nr-1,0:nz-1),cof(0:nr,0:nz)
    
      real wsavex(nr+15),wsavey(nz+15)

      integer i,k
      real f0(0:nr-1,0:nz-1)

      real tarrx(1:nr,1:nz),tarry(1:nz,1:nr)

      real temp(0:nz-1,0:nr-1)!,tempp(-1:nr+1,-1:nz+1)
    
    
      f0 = fin

      call VRFFTF(nr,nz,f0,tarrx,nr,wsavey)

      do k=0,nr-1
         do i=0,nz-1
            temp(i,k) = f0(k,i)
         enddo
      enddo
    
      call VRFFTF(nz,nr,temp,tarry,nz,wsavex)
    
      do k=0,nr-1
         do i=0,nz-1
            temp(i,k) = temp(i,k)/cof(k,i)
         enddo
      enddo
    
      call VRFFTB(nz,nr,temp,tarry,nz,wsavex)
    
       do k=0,nz-1
          do i=0,nr-1
             fout(i,k) = temp(k,i)
          enddo
       enddo
    
    
      call VRFFTB(nr,nz,fout,tarrx,nr,wsavey)
    
    return 
   
   end subroutine  ! FFT_regular