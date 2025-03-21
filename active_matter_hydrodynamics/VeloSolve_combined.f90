!#define FULL_IMPLICIT
!#define SEMI_IMPLICIT

  subroutine VeloSolve_combined(psi,u_after,v_after,rmu_after,ph_after,t,&
    & particle_force1, particle_force2)
  use global
  implicit none
  integer i,k
  real t
  real psi(-1:nr,-1:nz)

  real u_after(-1:nr,-1:nz)
  real v_after(-1:nr,-1:nz)

  real rmu_after(-1:nr,-1:nz)
  real ph_after(-1:nr,-1:nz)
  real particle_force1(-1:nr,-1:nz), particle_force2(-1:nr,-1:nz)

  real,external::force1,force2


  real fu(-1:nr,-1:nz),fv(-1:nr,-1:nz)

  real gT(0:nr),gB(0:nr),cf
  real cf1,cf2

  real temp_p(-1:nr,-1:nz) ! used for pressure, related to psi1, psi2

  ! ! store diff of eta on u and v locations
  ! real etay
  ! real tmp_hig

  integer, parameter::mn_u = nr*nz
  integer, parameter::mn_v = nr*nz

  real rhs_uv(mn_u + mn_v),ComputedSolution_uv(mn_u + mn_v)
  real work_uv(mn_u + mn_v,7)

  external MatrixVector_uv_CG, PSolve_uv_CG

  integer ITER,INFO
  real RESID
  logical is_debug

  is_debug = .false.

    !!pause
    if (is_debug) then
      call output_regular(nr,nz,r,z,fu,fv,ph_after,rmu_after,0.0,'testin1')
      print *,'output test 1'
    endif

!===================================================================================
!  first step to update velocity without enforce incompressbility
!==================================================================================

      ! ca*\mu*grad\phi
      do k=0,nz-1
         do i=0,nr-1
            fu(i,k) = ca*(rmu_after(i,k) + rmu_after(i-1,k))/2.0 * (ph_after(i,k) - ph_after(i-1,k)) /dr
            !fv(i,k) = ca*(rmu_after(i,k) + rmu_after(i,k-1))/2.0 * (ph_after(i,k) - ph_after(i,k-1)) /dz
            ! \phi is located at corners
            !fu(i,k) = ca*(rmu_after(i,k) + rmu_after(i,k+1))/2.0 * &
            !  & ( ph_after(i+1,k+1) - ph_after(i-1,k+1) + ph_after(i+1,k) - ph_after(i-1,k))/2.0 /2.0/dr
            !fv(i,k) = ca*(rmu_after(i,k) + rmu_after(i+1,k))/2.0 * &
            !  & ( ph_after(i,k+1) - ph_after(i,k-1) + ph_after(i+1,k+1) - ph_after(i+1,k-1))/2.0 /2.0/dz
         enddo
      enddo

      do k=0,nz-1
         do i=0,nr-1
            !fu(i,k) = ca*(rmu_after(i,k) + rmu_after(i-1,k))/2.0 * (ph_after(i,k) - ph_after(i-1,k)) /dr
            fv(i,k) = ca*(rmu_after(i,k) + rmu_after(i,k-1))/2.0 * (ph_after(i,k) - ph_after(i,k-1)) /dz
            ! \phi is located at corners
            !fu(i,k) = ca*(rmu_after(i,k) + rmu_after(i,k+1))/2.0 * &
            !  & ( ph_after(i+1,k+1) - ph_after(i-1,k+1) + ph_after(i+1,k) - ph_after(i-1,k))/2.0 /2.0/dr
            !fv(i,k) = ca*(rmu_after(i,k) + rmu_after(i+1,k))/2.0 * &
            !  & ( ph_after(i,k+1) - ph_after(i,k-1) + ph_after(i+1,k+1) - ph_after(i+1,k-1))/2.0 /2.0/dz
         enddo
      enddo


    if (is_debug) then
      call output_regular(nr,nz,r,z,fu,fv,rmu_after,ph_after,0.0,'testin2')
      print *,'output test 2'
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! re\rho*u\cdot\grad u
    do k=0,nz-1
          do i=0,nr-1
            fu(i,k) = fu(i,k) - re*(rho(i,k) + rho(i-1,k))/2.0 * ( u(i,k)* (u(i+1,k)-u(i-1,k) )/2.0/dr &
              & + (v(i,k) + v(i-1,k) + v(i,k+1) + v(i-1,k+1))/4.0 * ( u(i,k+1) - u(i,k-1) )/2.0/dz )
          end do
    end do

    do k=0,nz-1
          do i=0,nr-1
            fv(i,k) = fv(i,k) - re*(rho(i,k) + rho(i,k-1))/2.0 * &
            & ( (u(i,k) + u(i+1,k) + u(i,k-1) + u(i+1,k-1))/4.0 * ( v(i+1,k)-v(i-1,k) )/2.0/dr &
            & + v(i,k)* (v(i,k+1)-v(i,k-1))/2.0/dz )
          end do
    end do


    if (is_debug) then
      call output(nr,nz,r,z,fu,fv,p,ph,0.0,'testin3')
      print *,'output test 3'
    endif

!    call output(nr,nz,r,z,fu,fv,u,ph,0.0,'testin5')
!        print *,'output test 5'
!        pause

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef SEMI_IMPLICIT
!
!    temp = 2*u-u_pre
!    tempp = 2*v-v_pre
!
!    call  dif1x(nr,nz,xa1,xb1,xc1,temp,w1)
!    call  dif1y(nr,nz,ya1,yb1,yc1,temp,w2)
!    call  dif1x(nr,nz,xa1,xb1,xc1,tempp,w3)
!    call  dif1y(nr,nz,ya1,yb1,yc1,tempp,w4)
!
!    ! at locations of u
!    ! 2\eta_x * u_x + \eta_y * (u_y + v_x)
!    do k=0,nz
!      do i=0,nr
!        eta_x(i,k) = (eta(i,k) - eta(i-1,k))/dr  ! \partial_x\eta
!        eta_y(i,k) = (eta(i,k+1) - eta(i,k-1) + eta(i-1,k+1) - eta(i-1,k-1))/2.0/2.0/dz  ! \partial_y\eta
!      enddo
!    enddo
!
!    do k=0,nz
!      do i=0,nr
!        fu(i,k) = fu(i,k) + 2.*eta_x(i,k) *w1(i,k) &
!        & + eta_y(i,k)* (w2(i,k) + (tempp(i,k)-tempp(i-1,k) + tempp(i,k+1)-tempp(i-1,k+1))/2.0/dr)!\partial_x v
!      enddo
!    enddo
!
!
!    ! at locations of v
!    ! \eta_x * (u_y + v_x) + 2\eta_y * v_y
!    do k=0,nz
!      do i=0,nr
!        eta_x(i,k) = (eta(i+1,k)-eta(i-1,k) + eta(i+1,k-1)-eta(i-1,k-1))/2.0/2.0/dr ! \partial_x\eta
!        eta_y(i,k) = (eta(i,k)-eta(i,k-1))/dz  ! \partial_y\eta
!      enddo
!    enddo
!
!    do k=0,nz
!      do i=0,nr
!        fv(i,k) = fv(i,k) + 2.*eta_y(i,k)*w4(i,k) &
!        & + eta_x(i,k)*((temp(i+1,k)-temp(i+1,k-1) + temp(i,k)-temp(i,k-1))/2.0/dz  + w3(i,k))
!      enddo
!    enddo
!#endif


    if (is_debug) then
      call output(nr,nz,r,z,fu,fv,p,ph,0.0,'testin4')
      print *,'output test 4'
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! pressure related terms
!    call dif1x(nr,nz,xa1,xb1,xc1,p_pre,w1)
!    call dif1y(nr,nz,ya1,yb1,yc1,p_pre,w2)

    !call dif1x(nr,nz,xa1,xb1,xc1,p+4*psi2/3.-psi1/3.,w3)
    !call dif1y(nr,nz,ya1,yb1,yc1,p+4*psi2/3.-psi1/3.,w4)

    !TODO: check
    temp_p = p !+ psi
!    do k=0,nz
!      do i=0,nr
!        w3(i,k) = (temp_p(i+1, k) - temp_p(i, k) + temp_p(i+1, k+1) - temp_p(i, k+1))/2.0/dr
!        w4(i,k) = (temp_p(i, k+1) - temp_p(i, k) + temp_p(i+1, k+1) - temp_p(i+1, k))/2.0/dz
!      enddo
!    enddo

    ! p_x, p_y, re\rho*u/dt + ggx*\rho
    do k=0,nz-1
       do i=0,nr-1
          !fu(i,k)=fu(i,k) - 2*w3(i,k)*dt + w1(i,k)*dt + 2*u(i,k)*re - u_pre(i,k)*re/2.
          !fv(i,k)=fv(i,k) - 2*w4(i,k)*dt + w2(i,k)*dt + 2*v(i,k)*re - v_pre(i,k)*re/2.
          fu(i,k)=fu(i,k) - (temp_p(i, k) - temp_p(i-1, k))/dr &
           &+ (rho(i,k)+rho(i-1,k))/2.0*( u(i,k)*re/dt + ggx )
       end do
    enddo

    do k=0,nz-1
       do i=0,nr-1
          fv(i,k)=fv(i,k) - (temp_p(i, k) - temp_p(i, k-1))/dz &
           & + (rho(i,k)+rho(i,k-1))/2.0*( v(i,k)*re/dt + ggz )
       end do
    enddo


    if (is_debug) then
      call output(nr,nz,r,z,fu,fv,p,ph,0.0,'testin5')
      print *,'output test 5'
    endif

    ! external force

    !   do k=0,nz
    !      do i=0,nr
    !         fu(i,k) = fu(i,k) + (rho(i,k)+rho(i-1,k))/2.0*force1( r(i,k),z(i,k),t )
    !         fv(i,k) = fv(i,k) + (rho(i,k)+rho(i,k-1))/2.0*force2( r(i,k),z(i,k),t )
    !      enddo
    !   enddo


     ! TODO: force with density? and coefficient for scalinga
    do k=0,nz-1
       do i=0,nr-1
         fu(i,k) = fu(i,k) + (particle_force1(i,k) + particle_force1(i-1,k))/2.0
       end do
    enddo

    do k=0,nz-1
       do i=0,nr-1
         fv(i,k) = fv(i,k) + (particle_force2(i,k) + particle_force2(i,k-1))/2.0
       end do
    enddo



    ! !for test
    ! particle_force1(-1:nr,-1:nz) = fu(-1:nr,-1:nz)
    ! particle_force2(-1:nr,0:nz) = fv(-1:nr,0:nz)


    do k=0,nz-1
       do i=0,nr-1
         fu(i,k) = fu(i,k)*dt !/ 2.0*(eta(i,k) + eta(i-1,k))
       end do
    enddo

    do k=0,nz-1
       do i=0,nr-1
         fv(i,k) = fv(i,k)*dt !/ 2.0*(eta(i,k) + eta(i,k-1))
       end do
    enddo

    if (is_debug) then
      call output(nr,nz,r,z,fu,fv,p,ph,0.0,'testin6')
      print *,'output test 6'
      !pause
    endif


   !! Dirichlet BDs on L  BDr
   ! TODO
   !cf = 1.0 - exp(-t/twall)
   !cf1 = 1-shig/12/(xlsL+shig/4)  !scaling coef at left BDr
   !cf2 = 1-shig/12/(xlsG+shig/4) !scaling coef at right BDr

!   i = 0
!   do k=0,nz-1
!       !uL(k) = UU*( -(z(0,k)-shig/2)**2/shig/(xlsL+shig/4) +1 )*cf/cf1
!       !uR(k) = UU*( -(z(nr,k)-shig/2)**2/shig/(xlsG+shig/4) +1 )*cf/cf2
!       !fu(1,k) = fu(1,k) + uL(k)*dt/dr2
!       !fu(nr-1,k) = fu(nr-1,k) + uR(k)*dt/dr2
!
!       ! here bdr_ph is simply used to modify inlet velocity
!       ! when bdr_ph = 1.0 means boundary and will multiply by zero
!       !uL(k) = UU*( -(z(0,k)-shig/2)**2/shig/shig*4 + 1 )*cf * (1.0 - bdr_ph(0,k))
!
!       !if (z(0,k) > bdr_low .AND. z(0,k) < bdr_upp) then
!       uL(k) = UU  !*( -(z(0,k)-shig/2)**2/shig/shig*4 + 1 ) !*cf
!       !else
!       !  uL(k) = 0.0
!       !endif
!
!
!       !cr = eta(i,k)
!       !cl = eta(i-1,k)
!
!       fu(i,k) = fu(i,k) + 2.0*uL(k)*eta(i,k)*2.0*dt/dr2
!   enddo
!
!   i = nr
!   do k=0,nz-1
!       ! here simply use uL as same velocity to satisfy integral of div = 0
!       uR(k) = uL(k) !UU*( -(z(0,k)-shig/2)**2/shig/shig*4 + 1 )*cf * (1.0 - bdr_ph(0,k))
!
!       fu(i,k) = fu(i,k) + 2.0*uR(k)*eta(i,k)*2.0*dt/dr2
!   enddo


!!TODO
if(xlsB >= 1.E-8) then

!!!!! the value of u at boundary Robin BCs
!
!      do i=0,nr
!         cf=sqrt(2.0)*cos(angleB(i)*pi/180.0)/3.0*pi/2.
!         alphaB(i) = 0.5*(1.0+ph_after(i,0))/xlsB+0.5*(1.0-ph_after(i,0))/(xlsB*slipratio)
!
!         gB(i) =  uwB*alphaB(i) + ca*( I_thickness*(ph_after(i,-1)-ph_after(i,1))/(2.*dz)&
!& - 2*cf*cos(pi*ph(i,0)/2.) + cf*cos(pi*ph_pre(i,0)/2.) + alpha_s*(ph_after(i,0)-2*ph(i,0)+ph_pre(i,0)) )&
!&*(ph_after(i+1,0)-ph_after(i-1,0))/(2.*dr)/eta(i,0)
!      end do
!
!      do i=0,nr
!         cf=sqrt(2.0)*cos(angleT(i)*pi/180.0)/3.0*pi/2.
!         alphaT(i) = 0.5*(1.0+ph_after(i,nz))/xlsT+0.5*(1.0-ph_after(i,nz))/(xlsT*slipratio)
!
!         gT(i) = uwT*alphaT(i) + ca*( I_thickness*(ph_after(i,nz+1)-ph_after(i,nz-1))/(2.*dz)&
!& - 2*cf*cos(pi*ph(i,nz)/2.) + cf*cos(pi*ph_pre(i,nz)/2.) + alpha_s*(ph_after(i,nz)-2*ph(i,nz)+ph_pre(i,nz)) )&
!&*(ph_after(i+1,nz)-ph_after(i-1,nz))/(2.*dr)/eta(i,nz)
!      end do
!
!     do i=0,nr
!        fu(i,0) = fu(i,0) + 2.*dt*gB(i)/dz
!        fu(i,nz) = fu(i,nz) + 2.*dt*gT(i)/dz
!     enddo
!
!
!!     call ESOLP_x(nr,nz,ua_2,ub_2,uc_2,fu,u_after,wsavexf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!    transform the RHS to a single vector by calling AtoB
!!!!         initial guess for the solution using the value at last time step
!
!    !call Atob(nr,nz,fu,rhs_s)
!
!    !call Atob(nr,nz,u,ComputedSolution_s)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  call FFT based preconditioned GMRES method to solve linear system of PH
!!!  here 150 is the restart number of GMRES Iteration
!
!       ! call output(nr,nz,r,z,u,v,fu,fv,0.0,'testin5')
!       ! print *,'output testin5'
!!	pause
!
!!!	call fgmres((nr+1)*(nz+1),150,ComputedSolution_s,rhs_s,tol,MatrixVector_u_2_s,PSolve_u_2_s)
!!ITER = 500
!!RESID = tol
!!    call BICGSTAB((nr+1)*(nz+1),rhs_s,ComputedSolution_s,work,(nr+1)*(nz+1),ITER,RESID,MatrixVector_u_2_s_CG,PSolve_u_2_s_CG,INFO)
!!	write(*,*) 'u 2 BICGSTAB, slip',ITER
!!	if(INFO .ne. 0) then
!!		print *,'error'
!!!		stop
!!	endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! transform the solution back to matrix form
!!
!!    call btoA(nr,nz,u_after,ComputedSolution_s)
!
!!!!! ghost points for u_after
!   do k=0,nz
!       u_after(-1,k) = u_after(1,k)
!       u_after(nr+1,k) = u_after(nr-1,k)
!   enddo
!
!     do i=0,nr
!        u_after(i,-1) = u_after(i,1) + 2.*dz*(gB(i)-u_after(i,0)*alphaB(i))
!        u_after(i,nz+1) = u_after(i,nz-1) + 2.*dz*(gT(i)-u_after(i,nz)*alphaT(i))
!     enddo

else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!    transform the RHS to a single vector by calling AtoB
!!!         initial guess for the solution using the value at last time step

!    !!!! the value of u at boundary BCs
!    !TODO: checked  TODO
!     do i=0,nr
!        k = 0
!        etay = ( eta(i,k) + eta(i-1,k) + eta(i,k-1) + eta(i-1,k-1) ) /4.0
!        fu(i,k) = fu(i,k) + dt*(2.0*uwB/dz2 *etay)
!
!        k = nz - 1
!        etay = ( eta(i,k+1) + eta(i-1,k+1) + eta(i,k) + eta(i-1,k) ) /4.0
!        fu(i,k) = fu(i,k) + dt*(2.0*uwT/dz2 *etay)
!     enddo


    call Atob_uv_periodic(nr,nz,fu,fv,rhs_uv)
    call Atob_uv_periodic(nr,nz,u,v,ComputedSolution_uv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  call FFT based preconditioned GMRES method to solve linear system of PH
!!  here 150 is the restart number of GMRES Iteration

    !call fgmres_u_2((nr+1)*nz,150,ComputedSolution_u,rhs_u,tol)
    !call fgmres_u_2((nr-1)*(nz-1),150,ComputedSolution,rhs,tol)

ITER = 800
RESID = tol
    call BICGSTAB(nr*nz + nr*nz,rhs_uv,ComputedSolution_uv,work_uv,mn_u + mn_v,&
      &ITER,RESID,MatrixVector_uv_CG,PSolve_uv_CG,INFO)
	write(*,*) 'uv 1 BICGSTAB no slip',ITER
	if(INFO .ne. 0) then
		print *,'uv 1 BICGSTAB no slip: error',ITER,INFO
!		stop
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! transform the solution back to matrix form

    call btoA_uv_periodic(nr,nz,u_after,v_after,ComputedSolution_uv)

    !TODO: need extend v to [0, nz-1] for periodic bdrs
    ! Here, we don't extend boundary, and instead approximate the periodic
    ! bdrs by interpolate values at boundary. in this case, we don't need change
    ! the size of v in iteration process

!    call output_regular(nr,nz,r,z,u_after,v_after,u,v,0.0,'testin2')
!    print *,'output test2'
!    pause


!    !----------------------------------------------------------------
!    ! modify velicity to enforce zero velocity of boundary
!    do k=0,nz-1
!       do i=0,nr
!        tmp = (ph_after(i,k) + ph_after(i,k-1))/2.0
!        if ( tmp > 0.5 ) then
!          u_after(i,k) = 0.0
!        else
!          u_after(i,k) = u_after(i,k) * (1.0 - 2.0*tmp)
!        endif
!       end do
!    enddo
!
!    do k=1,nz-1
!       do i=0,nr-1
!        tmp = (ph_after(i,k) + ph_after(i,k-1))/2.0
!        if ( tmp > 0.5 ) then
!          v_after(i,k) = 0.0
!        else
!          v_after(i,k) = v_after(i,k) * (1.0 - 2.0*tmp)
!        endif
!       end do
!    enddo
!
!    ! bdr
!    do k=0,nz-1
!       u_after(-1,k) = 2.0*uL(k) - u_after(1,k)
!    enddo
!
!    do k=0,nz-1
!       u_after(nr+1,k) = 2.0*uR(k) - u_after(nr-1,k)
!    enddo
!
!    do i=-1,nr+1
!       u_after(i,-1) = u_after(i,nz-1)
!       u_after(i,nz) = u_after(i,0)
!    enddo
!
!    !TODO: need extend v to [0, nz-1] for periodic bdrs
!    ! Here, we don't extend boundary, and instead approximate the periodic
!    ! bdrs by interpolate values at boundary. in this case, we don't need change
!    ! the size of v in iteration process
!    do i=-1,nr
!       v_after(i,0) = ( v_after(i,nz-1) + v_after(i, 1) )/2.0
!       v_after(i,nz) = ( v_after(i,nz-1) + v_after(i,1) ) /2.0
!    enddo
!
!    do k=0,nz
!       v_after(-1,k) = v_after(0,k)
!       v_after(nr,k) = v_after(nr-1,k)
!    enddo


endif

    !call output(nr,nz,r,z,fu,fv,p,ph,0.0,'testing')
    !print *,'output test'
    !pause

!!!!!  discrete transform for velocity v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    no flux of incomressipility condition \partial_n (\grad \cdot \bf u)=0 at BD
!u_after(-1,nz+1)=u_after(nr-1,nz+1)
!u_after(-1,-1)=u_after(nr-1,-1)
!u_after(nr+1,nz+1)=u_after(1,nz+1)
!u_after(nr+1,-1)=u_after(1,-1)
!
!if(xlsB >= 1.E-8) then
!  do i=0,nr
!     k=0
!     v_after(i,k-1) = 2*v_after(i,k) - v_after(i,k+1) - dz/2*( (u_after(i+1,k+1)-u_after(i-1,k+1))/dr/2 &
!  & - (u_after(i+1,k-1)-u_after(i-1,k-1))/2/dr )
!     k=nz
!     v_after(i,k+1) = 2*v_after(i,k) - v_after(i,k-1) - dz/2*( (u_after(i+1,k+1)-u_after(i-1,k+1))/dr/2&
!  & - (u_after(i+1,k-1)-u_after(i-1,k-1))/2/dr )
!  enddo
!else !added
!  do i=0,nr
!     k=0
!     v_after(i,k-1) = 2*v_after(i,k) - v_after(i,k+1) - dz*( 4*(u_after(i+1,k+1)-u_after(i-1,k+1)) - &
!&(u_after(i+1,k+2)-u_after(i-1,k+2)) - 3*(u_after(i+1,k)-u_after(i-1,k)) )/2/dr/2
!     k=nz
!     v_after(i,k+1) = 2*v_after(i,k) - v_after(i,k-1) - dz*( (u_after(i+1,k-2)-u_after(i-1,k-2)) - &
!&4*(u_after(i+1,k-1)-u_after(i-1,k-1)) + 3*(u_after(i+1,k)-u_after(i-1,k))  )/2/dr/2
!  enddo
!endif



!!!!!  discrete transform for velocity v

!do i=0,nr
!   gB(i) = (u_after(i+1,0)-u_after(i-1,0))/2./dr
!   fv(i,0) = fv(i,0) + 2.*dt*gB(i)/dz
!
!   gT(i) = -(u_after(i+1,nz)-u_after(i-1,nz))/2./dr
!   fv(i,nz) = fv(i,nz) + 2.*dt*gT(i)/dz
!enddo

!     call DCT(nr,nz,fv,v_after,Ve_Lc_2,wsavexc,wsaveyc)

!!!!!!!!!!!!!!!!!  left and right bcd
!   do k=0,nz
!!       v_after(-1,k) = v_after(1,k)
!       v_after(nr+1,k) = v_after(nr-1,k)
!   enddo

!!!!!!! Top BD&Bottom BD
!    do i=0,nr
!     v_after(i,-1)=v_after(i,1) + 2.*dz*gB(i)
!     v_after(i,nz+1)=v_after(i,nz-1) + 2.*dz*gT(i)
!    enddo
!    !! Top BD&Bottom BD
!
!    do i=0,nr
!       v_after_after(i,-1)=v_after(i,1) + 2.*dz*gB(i)
!       v_after_after(i,nz+1)=v_after(i,nz-1) + 2.*dz*gT(i)
!    enddo
  end subroutine




  subroutine  MatrixVector_uv_CG(alpha,xx,beta,yy)
  use global
  implicit none
  real xx(nr*nz + nr*nz), yy(nr*nz + nr*nz)
  real zz(nr*nz + nr*nz)
  real alpha,beta

  real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)
  real f11(-1:nr,-1:nz),f22(-1:nr,-1:nz)

  ! store diff of eta on u and v locations
  !real eta_x(0:nr,0:nz),eta_y(0:nr,0:nz)
  !real w1(-1:nr+1,-1:nz+1),w2(-1:nr+1,-1:nz+1),w3(-1:nr+1,-1:nz+1),w4(-1:nr+1,-1:nz+1)

  integer i,k
  real cl, cr ! temporay coefficients, cl: for coef left, cr: coef right


  call btoA_uv_periodic(nr,nz,f1,f2,xx)

!!#ifdef SEMI_IMPLICIT
!    do k=0,nz-1
!       do i=0,nr
!          f11(i,k) = f1(i,k)*1.5*re*(rho(i,k) + rho(i-1,k))/2.0 &
!            & - dt*( (f1(i+1,k)+f1(i-1,k)-2*f1(i,k))/dr2 + (f1(i,k+1)+f1(i,k-1)-2*f1(i,k))/dz2 )&
!            & *(eta(i,k) + eta(i-1,k))/2.0
!        enddo
!    enddo
!
!    do k=1,nz-1
!       do i=0,nr-1
!          f22(i,k) = f2(i,k)*1.5*re* (rho(i,k) + rho(i,k-1))/2.0 &
!            & - dt*( (f2(i+1,k)+f2(i-1,k)-2*f2(i,k))/dr2 + (f2(i,k+1)+f2(i,k-1)-2*f2(i,k))/dz2 )&
!            & *(eta(i,k) + eta(i,k-1))/2.0
!        enddo
!    enddo
!!#endif


!#ifdef FULL_IMPLICIT
    ! u
    do k=0,nz-1
       do i=0,nr-1
          f11(i,k) = f1(i,k)*re*(rho(i,k) + rho(i-1,k))/2.0
            !& - dt*( (f1(i+1,k)+f1(i-1,k)-2*f1(i,k))/dr2 + (f1(i,k+1)+f1(i,k-1)-2*f1(i,k))/dz2 ) &
            !& *(eta(i,k) + eta(i-1,k))/2.0
        enddo
    enddo

    ! 2\partial_x(\eta\partial_x u)
    do k=0,nz-1
       do i=0,nr-1
        cr = eta(i,k)
        cl = eta(i-1,k)
        f11(i,k) = f11(i,k) &
           & - 2.0*dt*( cr*f1(i+1,k) + cl*f1(i-1,k) - (cl+cr)*f1(i,k) )/dr2
        enddo
    enddo

    ! \partial_y( \eta (\partial_y u + \partial_x v) )
    do k=0,nz-1
       do i=0,nr-1
        cr = ( eta(i,k+1) + eta(i-1,k+1) + eta(i,k) + eta(i-1,k) ) /4.0
        cl = ( eta(i,k) + eta(i-1,k) + eta(i,k-1) + eta(i-1,k-1) ) /4.0
        f11(i,k) = f11(i,k) &
           & - dt*( cr*f1(i,k+1) + cl*f1(i,k-1) - (cl+cr)*f1(i,k) )/dz2

        f11(i,k) = f11(i,k) &
           & - dt*( cr*( f2(i,k+1) - f2(i-1,k+1) ) - cl*( f2(i,k) - f2(i-1,k) ) )/dr/dz
        enddo
    enddo

    ! v
    do k=0,nz-1
       do i=0,nr-1
          f22(i,k) = f2(i,k)*re* (rho(i,k) + rho(i,k-1))/2.0
            !& - dt*( (f2(i+1,k)+f2(i-1,k)-2*f2(i,k))/dr2 + (f2(i,k+1)+f2(i,k-1)-2*f2(i,k))/dz2 ) &
            !& *(eta(i,k) + eta(i,k-1))/2.0
        enddo
    enddo

    ! 2\partial_y(\eta\partial_y v)
    do k=0,nz-1
       do i=0,nr-1
         cr = eta(i,k)
         cl = eta(i,k-1)
         f22(i,k) = f22(i,k) &
           & - 2.0*dt*( cr*f2(i,k+1) + cl*f2(i,k-1) - (cl+cr)*f2(i,k) )/dz2
        enddo
    enddo

    ! \partial_x( \eta (\partial_x v + \partial_y u) )
    do k=0,nz-1
       do i=0,nr-1
         cr = ( eta(i+1,k) + eta(i,k) + eta(i+1,k-1) + eta(i,k-1) )/4.0
         cl = ( eta(i,k) + eta(i-1,k) + eta(i,k-1) + eta(i-1,k-1) )/4.0
         f22(i,k) = f22(i,k) &
           & - dt*( cr*f2(i+1,k) + cl*f2(i-1,k) - (cl+cr)*f2(i,k) )/dr2

         f22(i,k) = f22(i,k) &
           & - dt*( cr*( f1(i+1,k) - f1(i+1,k-1) ) - cl*( f1(i,k) - f1(i,k-1) ) )/dr/dz
        enddo
    enddo

!#endif

    call Atob_uv_periodic(nr,nz,f11,f22,zz)
    yy = alpha*zz + beta*yy
  end subroutine


  subroutine PSolve_uv_CG(yy,xx)
  use global
  implicit none
  real xx(nr*nz + nr*nz), yy(nr*nz + nr*nz)

  real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)
  real f11(-1:nr,-1:nz),f22(-1:nr,-1:nz)

  integer i,j
  real tmp, cr,cl


  ! yy = xx
  ! return

   ! Jacobi preconditioner
  call btoA_uv_periodic(nr,nz,f1,f2,xx)

    do j=0,nz-1
       do i=0,nr-1
         tmp = re*(rho(i,j) + rho(i-1,j))/2.0

         cr = eta(i,j)
         cl = eta(i-1,j)
         tmp = tmp + 2.0*dt*(cl+cr)/dr2

          cr = ( eta(i,j+1) + eta(i-1,j+1) + eta(i,j) + eta(i-1,j) ) /4.0
          cl = ( eta(i,j) + eta(i-1,j) + eta(i,j-1) + eta(i-1,j-1) ) /4.0

          tmp = tmp + dt*(cl+cr)/dz2

          f11(i,j) = f1(i,j) / tmp
        enddo
    enddo

    do j=0,nz-1
       do i=0,nr-1
        tmp = re*(rho(i,j) + rho(i-1,j))/2.0

        cr = eta(i,j)
        cl = eta(i,j-1)
        tmp = tmp + 2.0*dt*(cl+cr)/dz2

        cr = ( eta(i+1,j) + eta(i,j) + eta(i+1,j-1) + eta(i,j-1) )/4.0
        cl = ( eta(i,j) + eta(i-1,j) + eta(i,j-1) + eta(i-1,j-1) )/4.0

        tmp = tmp + dt*(cl+cr)/dr2

        f22(i,j) = f2(i,j) / tmp
      enddo
    enddo

    call Atob_uv_periodic(nr,nz,f11,f22,yy)


!    ! FFT based preconditioner,needed to be improved
!    call btoA_uv(nr,nz,f1,f2,xx)
!
!     do k=0,nz-2
!        do i=0,nr-1
!           fvv(i,k) = f1(i,k+1)
!        enddo
!     enddo
!
!     call ESOLP_x(nr-1,nz-2,va_2,vb_2,vc_2,fvv,vv,wsavexc)
!
!     do k=1,nz-1
!        do i=0,nr-1
!           f11(i,k) = vv(i,k-1)
!        enddo
!     enddo
!
!
!     do k=0,nz-2
!        do i=0,nr-1
!           fvv(i,k) = f2(i,k+1)
!        enddo
!     enddo
!
!     call ESOLP_x(nr-1,nz-2,va_2,vb_2,vc_2,fvv,vv,wsavexc)
!
!     do k=1,nz-1
!        do i=0,nr-1
!           f22(i,k) = vv(i,k-1)
!        enddo
!     enddo
!
!    call Atob_uv(nr,nz,f11,f22,yy)


!  yy = xx
!  return

!    call btoA_uv(nr,nz,f1,f2,xx)
!
!    call ESOLP_x(nr,nz-1,ua,ub,uc,f1,f1,wsavexc)
!
!    !TODO: here whether f2(-1:nr,0:nz) will be maped to b(-1:nr,-1:nz-1)
!    call ESOLP_x(nr-1,nz-2,va,vb,vc,f2,f2,wsavexc_v)
!
!    call Atob_uv(nr,nz,f1,f2,yy)


  end subroutine






