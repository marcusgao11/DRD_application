    program main
    use global
    use particles
    use particles_move
    implicit none

    character(len=7)::fname
    character,dimension(0:9)::ch=(/'0','1','2','3','4','5','6','7','8','9'/)

    integer istep,iflag,i,j,k1,k2,k3,kcount,k
    real tf,t,time_step_change

    real gamma1
    real tmp

    real time_begin,time_end
    real dissip1,dissip2,mass1,mass2,mass0

    !integer,external:: OMP_GET_NUM_THREADS


    real u_pre(-1:nr+1,-1:nz),v_pre(-1:nr,0:nz)
    real u_after(-1:nr+1,-1:nz),v_after(-1:nr,0:nz)

    real ph_pre(-1:nr,-1:nz), ph_after(-1:nr,-1:nz)
    real rmu_after(-1:nr,-1:nz),rmu_pre(-1:nr,-1:nz)

    real psi1(-1:nr,-1:nz),psi2(-1:nr,-1:nz)
    real p_tmp(-1:nr,-1:nz)

    real particle_force1(-1:nr,-1:nz), particle_force2(-1:nr,-1:nz)

    !--------------------------------------------------------------------
    ! some open and close
    logical,parameter:: is_shift = .true.

    real particle_shift_dist !1.0 ! when particle distance to boundary smaller than
    ! this distance, all variables shift left
    integer shift_n, inter_n


    !real w1(-1:nr+1,-1:nz+1),w2(-1:nr+1,-1:nz+1)
    !real w3(-1:nr+1,-1:nz+1),w4(-1:nr+1,-1:nz+1)

    logical,parameter::is_fluid = .true.
    logical,parameter::is_modify_u = .FALSE. ! whether modify velocity


    ! initialization for pressure
    integer*8 plan
    integer*8 plan_back

    !-------------------------------------------------------
    !-  some important parameters
!    bdr_upp = 1.35
!    bdr_low = 0.35
!    particle_shift_dist = 1.0

    particle_shift_dist = 0.5

    !-------------------------------------------------------

    call PreSolve_init(nr,nz,plan,plan_back)
    print *,plan,plan_back
    write(*,*)

!parameters
        open(3,file='input')
        read(3,*) xld
        read(3,*) re
        read(3,*) ca
        read(3,*) vs
        read(3,*) xlsT
        read(3,*) xlsB
        read(3,*) xlsL
        read(3,*) xlsR
        read(3,*) uwT     ! wall speed top
        read(3,*) uwB     ! wall speed bottom
        read(3,*) slen    ! length of x directon
        read(3,*) shig      ! length of z directon
        read(3,*) angle  ! static contact angle
        read(3,*) slipratio
        read(3,*) iflag
        read(3,*) istep     !step to print
        read(3,*) tf              !! final time
        read(3,*) I_thickness           !! interface thickness ratio
        read(3,*) s_implicit           !! parameter to inforce satbility
        read(3,*) tol           !! parameter to inforce satbility
        read(3,*) ggz           !! external force along z direction
        read(3,*) ggx           !! external force along z direction
        read(3,*) time_step_change           !! change times step
        read(3,*) lambda_rho
        read(3,*) lambda_eta
        read(3,*) UU
        read(3,*) lr_bd_type
        read(3,*) top_type_ph
        read(3,*) bot_type_ph
        read(3,*) top_type_u
        read(3,*) bot_type_u

        close(3)

        write(6,100) 'xld= ', xld
        write(6,100) 're= ', re
        write(6,100) 'ca= ', ca
        write(6,100) 'vs= ', vs
        write(6,100) 'xlsT= ', xlsT
        write(6,100) 'xlsB= ', xlsB
        write(6,100) 'xlsL= ', xlsL
        write(6,100) 'xlsR= ', xlsR
        write(6,100) 'uwT= ', uwT
        write(6,100) 'uwB= ', uwB
        write(6,100) 'slen= ', slen
        write(6,100) 'shig= ', shig
        write(6,200) 'nr= ', nr
        write(6,200) 'nz= ', nz
        write(6,100) 'angle= ', angle
        write(6,100) 'slipratio= ', slipratio
        write(6,200) 'iflag= ', iflag
        write(6,100) 'tf= ', tf
        write(6,100) 'interface_thickness ratio',I_thickness           !! interface thickness ratio
        write(6,100) 's_implicit',s_implicit
        write(6,*) ''
        write(6,200) 'ouput data step',istep
        write(6,100) 'ggz',ggz
        write(6,100) 'ggx',ggx
        write(6,100) 'time_step_change',time_step_change
        write(6,100) 'lambda_rho',lambda_rho
        write(6,100) 'lambda_eta',lambda_eta
        write(6,100) 'UU', UU
        write(6,*) 'tolerence',tol
        write(*,*) 'lr_bd_type: ', lr_bd_type
        write(*,*) 'top_type_ph: ', top_type_ph
        write(*,*) 'bot_type_ph: ', bot_type_ph
        write(*,*) 'top_type_u: ', top_type_u
        write(*,*) 'bot_type_u: ', bot_type_u
        write(*,*) ''

!100 format(A30,10X,f10.5)
100 format(A30,10X,f10.5)
200 format(A30,10X,I4)

!!!!!! initialize some parameters
    kcount=0
    t=0.0
    mass0 = 0.5

    psi1 = 0.
    gamma1 = min(1.0,lambda_rho)  ! gamma in the scheme

    twall = 0.1 ! twall is used to smooth the velocity

    ! initialize all variables
    u = 0.0; v = 0.0; ph = 0.0; rmu = 0.0; p = 0.0; rho = 0.0; eta = 0.0;


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     call initial function: generate girds, time step, alpha_s,
!!                            initial condition for u v \phi \rmu, including ghost points

write(6,*) 're= finit', re
    call finit(time_step_change)

        print *,''
        write(6,*) 'dr= ', dr
        write(6,*) 'dz= ', dz
        write(6,*) 'mass0= ', mass0

    print *,''
    write(6,*) 'dt= ', dt
    write(6,*) 'alpha_s= ', alpha_s
    write(*,*) ''


    !=======================================================================
    ! particle initialize
    !call particle_initialize(slen/2.0, shig/2.0, slen, shig, is_fluid)
    call particle_initialize(slen/2.0, (0.17+0.67)/2.0, slen, shig, is_fluid)
    call interpolate_coef()

!  !TODO: how to choose suitable dt
!  write(*, *) "TODO: how to choose suitable dt"
!  dt = 1.0E-4


  ! !---------------------------------------------------------------
	! !  initialization for OPEN-OMP
  !  	call OMP_SET_NUM_THREADS(4)

  ! !$OMP PARALLEL
  ! !$OMP SINGLE
  ! print *,'numeber of OMP threads', OMP_GET_NUM_THREADS()
  ! !$OMP END SINGLE
  ! !$OMP END PARALLEL


  !pause


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  setup for the scheme
    call setup()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  output the intial data000 or input from data from data/flag

    if(iflag == 0) then
       print *,'output data000'
       call output(nr,nz,r,z,u,v,p,fpd_ph,ph,t,'data000')

       print *,'output databdr'
       call output(nr,nz,r,z,u,v,p,bdr_ph,ph,t,'databdr')

    else

       k1 = iflag/100
       k2 = (iflag - k1*100)/10
       k3 = iflag - k1*100 - k2*10
       fname = 'data'//ch(k1)//ch(k2)//ch(k3)

       print *,'input from ',fname
       call input(nr,nz,r,z,u,v,p,ph,t,fname)

       call bcd(t)
    endif


!	fname='testin1'
!        call output(nr,nz,r,z,u,v,p,u_after,t,fname)
!        print *,'output ',fname
!	pause


     call dissipation(ph_pre,dissip1,mass1,t)
     write(22,*) kcount,dt*dissip1,mass1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first order scheme is used as setup
    call PhSolve(rmu_after,ph_after,t+dt)
     !call ph_mu_solve(ph_after, rmu_after, t)

    if (is_fluid) then
      ! update particles
      !call force_caculation()
      !call update_position(dt)
      call force_caculation_vel(nr, nz, dr, dz, r, z, slen, shig,&
        & I_thickness, bdr_ph_dist, particle_force1, particle_force2, dt,t)
      call update_position_vel(nr,nz,r,z,slen, shig, u,v,dt,I_thickness)

      ! Calculte rho and eta
      call interpolate_coef()

    !   fname='testin1'
    !   call output(nr,nz,r,z,u,v,p,particle_force1,particle_force2,t,fname)
    !   print *,'output ',fname
	!     !pause

    endif

!!!!!!!!! test of cpu time per timestep
!    call cpu_time(time_begin)

! update without incompressibility constain
    if (is_fluid) then
      !TODO: to be modified
      call VeloSolve_combined(psi1,u_after, v_after, rmu_after,ph_after,&
      & t+dt, particle_force1, particle_force2)

      ! update velocity
      u_pre = u
      v_pre = v
      u = u_after
      v = v_after


      !!!!  test of cpu time per timestep
      !    call cpu_time(time_end)
      !    write(*,"('cpu time per time step  ',f10.6)"),(time_end-time_begin)/istep
      !    pause

      !project velocity onto divegence free space
      call PreSolve_2(plan,plan_back,psi2, re*gamma1/dt)

      !psi2 = psi2*re*gamma1/dt

      p = psi2 + p !- w3 - w4

    endif

!   !TODO: no need
!   do i=0,nr
!    p(i,-1) = 2*p(i,0) - p(i,1)
!    p(i,nz+1) = 2*p(i,nz) - p(i,nz-1)
!   enddo

!        fname='testin1'
!        call output(nr,nz,r,z,u,v,ph,rmu,t,fname)
!        print *,'output ',fname

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update phi
   ph_pre = ph
   rmu_pre = rmu

   rmu = rmu_after
   ph = ph_after

  !------------------------------------------------------------------------
  ! modify velicity to enforce zero velocity of boundary
  !------------------------------------------------------------------------
   if (is_modify_u) then
    do k=0,nz-1
      do i=0,nr
        tmp = (bdr_ph(i,k) + bdr_ph(i,k-1))/2.0
        !u(i,k) = u(i,k) * (1.0 - tmp)
        if ( tmp > 0.5 ) then
          u(i,k) = 0.0
        else
          u(i,k) = u(i,k) * (1.0 - 2.0*tmp)
        endif
      end do
    enddo

    do k=1,nz-1
      do i=0,nr-1
        tmp = (bdr_ph(i,k) + bdr_ph(i,k-1))/2.0
        !v(i,k) = v(i,k) * (1.0 - tmp)
        if ( tmp > 0.5 ) then
          v(i,k) = 0.0
        else
          v(i,k) = v(i,k) * (1.0 - 2.0*tmp)
        endif
      end do
    enddo
  endif



    t = t+dt
    kcount = kcount + 1

!     call bcd(t)

     if(kcount/istep*istep == kcount) then

       iflag=iflag+1
       k1 = iflag/100
       k2 = (iflag - k1*100)/10
       k3 = iflag - k1*100 - k2*10
       fname = 'data'//ch(k1)//ch(k2)//ch(k3)

        call output(nr,nz,r,z,u,v,p,fpd_ph,ph,t,fname)
        print *,'output ',fname
         if(iflag .GE. 999) then
           iflag = -1
           print *,'no engouh files to ouput'
           print *,'if you want to continue,it will input from the beginning'
         endif

       call dissipation(ph_pre,dissip2,mass2,t)
       write(22,*) kcount,dt*dissip2,mass2

       print *,'(dissip2-dissip1)/dissip1',(dissip2-dissip1)/dissip1
       print *,'mass 0 (mass2-mass1)/mass0',mass0,(mass2-mass1)/mass0
       print *,''

       dissip1 = dissip2
       mass1 = mass2

    endif

!!!!!!!!!!!!!!  time marching 2nd order
!
!!!!!!!!! test of cpu time per timestep
    call cpu_time(time_begin)


do while(t <= tf)


!call VeloSolve(psi1,rmu_after,ph_after,t+dt)

	!!  Calculate the ph, doesn't change ph and velocity. Store the new ph in phi_after
    !call PhSolve_2(u_pre,v_pre,ph_pre,ph_after,rmu_after,t+dt)
    !call ph_mu_solve(ph_after, rmu_after, t+dt)

    ! here we consider first order first
    call PhSolve(rmu_after,ph_after,t+dt)


  ! map \phi and \rum to centers

    if (is_fluid) then
      ! update particles
      !call force_caculation()
      !call update_position(dt)
      call force_caculation_vel(nr, nz, dr, dz, r, z, slen, shig,&
        & I_thickness, bdr_ph_dist, particle_force1, particle_force2, dt,t)
      call update_position_vel(nr,nz,r,z,slen, shig, u,v,dt,I_thickness)

      ! Calculte rho and eta
      call interpolate_coef()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  update the velocity without incompressibility constrain
          ! TODO: when shift, should deal with pressure
          ! since pressure condition is not easy shift, we simply
          ! do not use high order pressure, just no use pressure
          ! psi1 = 0.0; psi2 = 0.0; p = 0.0;

        !  fname='testin2'
        !  call output(nr,nz,r,z,u,v,eta,rho,ph,t,fname)
        !  print *,'output ',fname

        call VeloSolve_2_combined(psi1,psi2,u_pre,u_after,v_pre,v_after,ph_pre,ph_after,&
        rmu_pre,rmu_after,t+dt, particle_force1, particle_force2)

        ! !----------------------------------------------------------------
        ! ! modify velicity to enforce zero velocity of boundary
        ! if (is_modify_u) then
        !   do k=0,nz-1
        !     do i=0,nr
        !       tmp = (bdr_ph(i,k) + bdr_ph(i,k-1))/2.0
        !       if ( tmp > 0.5 ) then
        !         u_after(i,k) = 0.0
        !       else
        !         u_after(i,k) = u_after(i,k) * (1.0 - 2.0*tmp)
        !       endif
        !     end do
        !   enddo

        !   do k=1,nz-1
        !     do i=0,nr-1
        !       tmp = (bdr_ph(i,k) + bdr_ph(i,k-1))/2.0
        !       if ( tmp > 0.5 ) then
        !         v_after(i,k) = 0.0
        !       else
        !         v_after(i,k) = v_after(i,k) * (1.0 - 2.0*tmp)
        !       endif
        !     end do
        !   enddo
        ! endif

      !----------------------------------------------------------------
      ! modify velocity to enforce constant velocity inner particle
      ! implement in update_position_vel
        u_pre = u
        v_pre = v


        u = u_after
        v = v_after

    endif

    ! udpate phi
    ph_pre = ph
    rmu_pre = rmu

    ph = ph_after
    rmu = rmu_after


    !====================================================================
    !  project veolicity on the divegence free space
    !====================================================================
    if (is_fluid) then
      psi1 = psi2

      call PreSolve_2(plan,plan_back,psi2,3*re*gamma1/2./dt)

      !psi2 = psi2*3*re*gamma1/2./dt


      !  ! divergence of velocity,, if this is used, Flag: pressure ghost is needed
      !  do k=0,nz-1
      !     do i=0,nr-1
      !      w3(i,k) = ( u(i+1,k)-u(i,k) )/dr
      !      w4(i,k) = ( v(i,k+1)-v(i,k) )/dz
      !     enddo
      !  enddo


      ! here whether we need w3*eta - w4*eta
      p = psi2 + p !- w3*eta - w4*eta
    endif 

    !  Flag: pressure ghost, when use correction(w3*eta - w4*eta), we need ghost points for p
    !   do i=0,nr
    !    p(i,-1) = 2*p(i,0) - p(i,1)
    !    p(i,nz+1) = 2*p(i,nz) - p(i,nz-1)
    !   enddo

    ! !----------------------------------------------------------------
    ! ! modify velicity to enforce zero velocity of boundary
 
    if (is_modify_u) then
      do k=0,nz-1
        do i=0,nr
          tmp = (bdr_ph(i,k) + bdr_ph(i,k-1))/2.0
          !u(i,k) = u(i,k) * (1.0 - tmp)
          if ( tmp > 0.5 ) then
            u(i,k) = 0.0
          else
            u(i,k) = u(i,k) * (1.0 - 2.0*tmp)
          endif
        end do
      enddo
  
      do k=1,nz-1
        do i=0,nr-1
          tmp = (bdr_ph(i,k) + bdr_ph(i,k-1))/2.0
          !v(i,k) = v(i,k) * (1.0 - tmp)
          if ( tmp > 0.5 ) then
            v(i,k) = 0.0
          else
            v(i,k) = v(i,k) * (1.0 - 2.0*tmp)
          endif
        end do
      enddo
    endif

!     !----------------------------------------------------------------
!     !  shift particle to origins
!     !----------------------------------------------------------------
!     if (is_shift .AND. Np > 0 .AND. fpd_pos(0, 0) > slen -  particle_shift_dist) then

! !      print *, "shift-----------------------------------"
! !      inter_n = int( 2.0*particle_shift_dist/dr )
! !      shift_n = nr - inter_n
! !      print *, shift_n
! !
! !      u_pre(-1:-1+inter_n, :) = u_pre(-1+shift_n:nr-1, :)
! !      do i = inter_n,nr+1
! !        u_pre(i, :) = u_pre(nr, :)
! !      enddo
! !
! !      v_pre(-1:-1+inter_n,0:nz) = v_pre(-1+shift_n:nr-1,0:nz)
! !      do i = inter_n,nr
! !        v_pre(i, :) = v_pre(nr,:)
! !      enddo
! !
! !      u(-1:-1+inter_n, :) = u(-1+shift_n:nr-1, :)
! !      do i = inter_n,nr+1
! !        u(i, :) = u(nr, :)
! !      enddo
! !
! !      v(-1:-1+inter_n,0:nz) = v(-1+shift_n:nr-1,0:nz)
! !      do i = inter_n,nr
! !        v(i, :) = v(nr,:)
! !      enddo
! !
! !
! !      p(-1:-1+inter_n,-1:nz) = p(-1+shift_n:nr-1,-1:nz)
! !      do i = inter_n,nr
! !        p(i, :) = p(nr,:)
! !      enddo
! !
! !
! !      fpd_pos(0, 0) = fpd_pos(0, 0) - (slen - 2.0*particle_shift_dist)


!       ! way 2: shift half domain

!       !inter_n = int( 2.0*particle_shift_dist/dr )
!       shift_n = nr/2
!       print *, shift_n

!       u_after = u_pre ! for temp use
!       u_pre(0:shift_n-1, :) = u_after(shift_n:nr-1,:)
!       u_pre(shift_n:nr, :) = u_after(0:shift_n,:)
!       ! boundary....
!       u_pre(-1,:) = u_pre(0,:)
!       u_pre(nr+1,:) = u_pre(nr-1,:)

!       v_after = v_pre ! for temp use
!       v_pre(0:shift_n-1, :) = v_after(shift_n:nr-1,:)
!       v_pre(shift_n:nr, :) = v_after(0:shift_n,:)
!       ! boundary....
!       v_pre(-1,:) = v_pre(0,:)


!       u_after = u ! for temp use
!       u(0:shift_n-1, :) = u_after(shift_n:nr-1,:)
!       u(shift_n:nr, :) = u_after(0:shift_n,:)
!       ! boundary....
!       u(-1,:) = u(0,:)
!       u(nr+1,:) = u(nr-1,:)

!       v_after = v ! for temp use
!       v(0:shift_n-1, :) = v_after(shift_n:nr-1,:)
!       v(shift_n:nr, :) = v_after(0:shift_n,:)
!       ! boundary....
!       v(-1,:) = v(0,:)


!       p_tmp = p ! for temp use
!       p(0:shift_n-1, :) = p_tmp(shift_n:nr-1,:)
!       ! care needs be taken care, since pressure only continuous in sense of gradient
!       p(shift_n:nr, :) = p_tmp(0:shift_n,:) - p_tmp(0,0) + p_tmp(nr,0)
!       ! boundary....
!       p(-1,:) = p(0,:)

!       fpd_pos(0, 0) = fpd_pos(0, 0) - slen/2.0

!       ! no need to shift ph, since it's determined by particle position

!     endif

!      shift_n = int( (slen - 2.0*particle_shift_dist)/dr )
!      print *, shift_n



     t = t+dt
     kcount = kcount + 1

    !     call bcd(t)
    !!!!!!!!!!!!!!!!!!!!!!!!!!! end of forward Euler


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  output the data and carry out the stopping test by stability of dissipation

    !       call dissipation(ph_pre,dissip2,mass2,t)
    !       write(22,*) kcount,dt*dissip2,mass2

     if(kcount/istep*istep == kcount) then

!!!!  test of cpu time per timestep
    call cpu_time(time_end)
    write(*,"('cpu time per time step  ',f10.6)"),(time_end-time_begin)/istep
    !pause

       iflag=iflag+1
       k1 = iflag/100
       k2 = (iflag - k1*100)/10
       k3 = iflag - k1*100 - k2*10
       fname = 'data'//ch(k1)//ch(k2)//ch(k3)


        call output(nr,nz,r,z,u,v,p,fpd_ph,ph,t,fname)

!       fname = 'datb'//ch(k1)//ch(k2)//ch(k3)
!       call output_particles(fname)

!    call  dif1x(nr,nz,xa1,xb1,xc1,u,w3)
!    call  dif1y(nr,nz,ya1,yb1,yc1,v,w4)
! 	fname = 'divg'//ch(k1)//ch(k2)//ch(k3)
!        call output(nr,nz,w3+w4,z,u,v,p,ph,t,fname)

	print *,'output ',fname
         if(iflag .GE. 999) then
           iflag = -1
           print *,'no engouh files to ouput'
           print *,'if you want to continue,it will input from the beginning'
         endif

       call dissipation(ph_pre,dissip2,mass2,t)
       write(22,*) kcount,dt*dissip2,mass2

       print *,'(dissip2-dissip1)/dissip1',(dissip2-dissip1)/dissip1
       print *,'mass 0 (mass2-mass1)/mass0',mass0,(mass2-mass1)/mass0
       print *,''

       if( abs(dissip2-dissip1)/dissip2 <= 1.0E-4 ) then
          print *,'stable by dissipation'
!          stop
       endif

if(iflag .ge. 400) stop

!       if(k3+k2*10 .eq. 7) stop

       dissip1 = dissip2
       mass1 = mass2
    endif

enddo

 end program
