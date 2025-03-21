    program main
    use global
    implicit none

    character(len=7)::fname
    character,dimension(0:9)::ch=(/'0','1','2','3','4','5','6','7','8','9'/)

    integer istep,iflag,i,j,k, k1,k2,k3,kcount
    real tf,t,time_step_change

    real gamma1, tmp

    real time_begin,time_end
    real dissip1,dissip2,mass1,mass2,mass0


    real u_pre(-1:nr+1,-1:nz),v_pre(-1:nr,0:nz)
    real u_after(-1:nr+1,-1:nz),v_after(-1:nr,0:nz)

    real ph_pre(-1:nr,-1:nz), ph_after(-1:nr,-1:nz)
    real rmu_after(-1:nr,-1:nz),rmu_pre(-1:nr,-1:nz)

    real c_after(-1:nr,-1:nz)

    real psi1(-1:nr,-1:nz),psi2(-1:nr,-1:nz)

    logical,parameter::is_modify_u = .FALSE. ! whether modify velocity


    !real w1(-1:nr+1,-1:nz+1),w2(-1:nr+1,-1:nz+1)
    !real w3(-1:nr+1,-1:nz+1),w4(-1:nr+1,-1:nz+1)

    ! initialization for pressure
    integer*8 plan
    integer*8 plan_back

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
        read(3,*) uwT     ! wall speed top
        read(3,*) uwB     ! wall speed bottom
        read(3,*) slen    ! length of x directon
        read(3,*) shig      ! length of z directon
        read(3,*) angle_top  ! static contact angle
        read(3,*) angle_bot  ! static contact angle
        read(3,*) slipratio
        read(3,*) iflag
        read(3,*) istep     !step to print
        read(3,*) tf              !! final time
        read(3,*) I_thickness           !! interface thickness ratio
        read(3,*) s_implicit           !! parameter to inforce satbility
        read(3,*) alpha_s
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
        read(3,*) only_cahn_hillarid

        close(3)

        write(6,100) 'xld= ', xld
        write(6,100) 're= ', re
        write(6,100) 'ca= ', ca
        write(6,100) 'vs= ', vs
        write(6,100) 'xlsT= ', xlsT
        write(6,100) 'xlsB= ', xlsB
        write(6,100) 'uwT= ', uwT
        write(6,100) 'uwB= ', uwB
        write(6,100) 'slen= ', slen
        write(6,100) 'shig= ', shig
        write(6,200) 'nr= ', nr
        write(6,200) 'nz= ', nz
        write(6,100) 'angle_top ', angle_top
        write(6,100) 'angle_bot ', angle_bot
        write(6,100) 'slipratio= ', slipratio
        write(6,200) 'iflag= ', iflag
        write(6,100) 'tf= ', tf
        write(6,100) 'interface_thickness ratio',I_thickness           !! interface thickness ratio
        write(6,100) 's_implicit',s_implicit
        write(6,100) 'alpha_s', alpha_s
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

  !pause


! TODO:
C_0 = 0.8;  ! equal conc

D_0 = 1.0 !10.0; !1.0;  ! mobility

!c_bot = 1.0;  ! ~~
!c_top = 0.9;

c_left = 0.9

! initialize

!conc = 0.9  !0.85
!conc = c_top + (c_bot - c_top) * (1.0 - ph)/2.0

!conc = 0.9
!  do k=0,nz-1
!    do i=0,nr-1
!        if (ph(i,k) < 0 ) then
!          conc(i,k) = 1.0
!        else
!          conc(i,k) = 0.9
!        end if
!      end do
!  end do

conc = 0.0 !0.0  !0.85
 do k=0,nz
   do i=0,nr
       if (ph(i,k) < 0 ) then
         conc(i,k) = 0.85
       else
         conc(i,k) = 1.0
       end if
     end do
 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  setup for the scheme

    call setup()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  output the intial data000 or input from data from data/flag

    if(iflag == 0) then
       print *,'output data000'
       call output(nr,nz,r,z,u,v,p,ph,conc,t,'data000')

    else

       k1 = iflag/100
       k2 = (iflag - k1*100)/10
       k3 = iflag - k1*100 - k2*10
       fname = 'data'//ch(k1)//ch(k2)//ch(k3)

       print *,'input from ',fname
       call input(nr,nz,r,z,u,v,p,ph,conc, t,fname)

       call bcd(t)
    endif


!    ! for test
!    u = 0.2; v = 0.0

    ! initialize velocity and ph
    u_pre = u
    v_pre = v
    ph_pre = ph
    rmu_pre = rmu

    ph_after = ph ! when PhSolve is commented, this means ph not change

!	fname='testin1'
!        call output(nr,nz,r,z,u,v,p,u_after,t,fname)
!        print *,'output ',fname
!	pause


     call dissipation(ph_pre,dissip1,mass1,t)
     write(22,*) kcount,dt*dissip1,mass1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first order scheme is used as setup

! update ph_grad
! ph_grad for further use
  do k=0,nz-1
    do i=0,nr-1
      !TODO: other boundary conditions
      ! Neumann bcd
        if (i .eq. 0) then
            !ph_grad(i,k) = ((ph(i+1,k) - ph(i-1+nr,k))**2) / (4.0 * dr2)
            ph_grad(i,k) = 0.0  !((ph(i+1,k) - ph(i-1+nr,k))**2) / (4.0 * dr2)
        elseif (i .eq. nr-1) then
            ph_grad(i,k) = 0.0  !((ph(i+1-nr,k) - ph(i-1,k))**2) / (4.0 * dr2)
        else
            ph_grad(i,k) = ((ph(i+1,k) - ph(i-1,k))**2) / (4.0 * dr2)
        end if

        if (k .eq. 0) then
            ph_grad(i,k) = ph_grad(i,k) + 0.0  !((ph(i,k+1) - ph(i,k-1+nz))**2) / (4.0 * dz2)
        elseif (k .eq. nz-1) then
            ph_grad(i,k) = ph_grad(i,k) + 0.0  !((ph(i,k+1-nz) - ph(i,k-1))**2) / (4.0 * dz2)
        else
            ph_grad(i,k) = ph_grad(i,k) + ((ph(i,k+1) - ph(i,k-1))**2) / (4.0 * dz2)
        end if
      end do
  end do

  ph_grad = sqrt(ph_grad)*0.5



     call PhSolve(rmu_after,ph_after,t+dt)

     !call ph_mu_solve(ph_after, rmu_after, t)

     call diffusion_eqn(c_after,t)

    ! Calculte rho and eta
    call interpolate_coef(ph_after)


!!!!!!!!! test of cpu time per timestep
!    call cpu_time(time_begin)

! update without incompressibility constain

     !TODO: to be modified
     if (only_cahn_hillarid .ne. 'y') then
       call VeloSolve_combined(psi1,u_after, v_after, rmu_after,ph_after,t+dt)


        u_pre = u
        v_pre = v

        u = u_after
        v = v_after


!!!!  test of cpu time per timestep
!    call cpu_time(time_end)
!    write(*,"('cpu time per time step  ',f10.6)"),(time_end-time_begin)/istep
!    pause

!project velocity onto divegence free space
       call PreSolve_2(plan,plan_back,psi2)

       psi2 = psi2*re*gamma1/dt

       p = psi2 + p !- w3 - w4


      if (is_modify_u) then
        !----------------------------------------------------------------
        ! modify velocity to enforce zero velocity of boundary
!        do k=0,nz-1
!           do i=0,nr
!            tmp = (1.0 + (ph_after(i,k) + ph_after(i,k-1))/2.0) /2.0
!            if ( tmp > 0.5 ) then
!              u(i,k) = 0.0
!            else
!              u(i,k) = u(i,k) * (1.0 - 2.0*tmp)
!            endif
!           end do
!        enddo
!
!        do k=1,nz-1
!           do i=0,nr-1
!            tmp = (1.0 + (ph_after(i,k) + ph_after(i,k-1))/2.0) / 2.0
!            if ( tmp > 0.5 ) then
!              v(i,k) = 0.0
!            else
!              v(i,k) = v(i,k) * (1.0 - 2.0*tmp)
!            endif
!           end do
!        enddo

        do k=0,nz-1
           do i=0,nr
            tmp = (1.0 + (ph_after(i,k) + ph_after(i-1,k))/2.0) /2.0
            if ( tmp > 0.95 ) then
              u(i,k) = 0.0
            else
              u(i,k) = u(i,k) ! * (1.0 - 2.0*tmp)
            endif
           end do
        enddo

        do k=1,nz-1
           do i=0,nr-1
            tmp = (1.0 + (ph_after(i,k) + ph_after(i,k-1))/2.0) / 2.0
            if ( tmp > 0.95 ) then
              v(i,k) = 0.0
            else
              v(i,k) = v(i,k) !* (1.0 - 2.0*tmp)
            endif
           end do
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
      endif ! end of is_modify_u


     endif ! end of only_cahn_hillarid


!   !TODO: no need
!   do i=0,nr
!    p(i,-1) = 2*p(i,0) - p(i,1)
!    p(i,nz+1) = 2*p(i,nz) - p(i,nz-1)
!   enddo

!        fname='testin1'
!        call output(nr,nz,r,z,u,v,ph,rmu,t,fname)
!        print *,'output ',fname

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ph_pre = ph
    rmu_pre = rmu

    rmu = rmu_after
    ph = ph_after
    conc = c_after

    t = t+dt
    kcount = kcount + 1

!     call bcd(t)

     if(kcount/istep*istep == kcount) then
       iflag=iflag+1
       k1 = iflag/100
       k2 = (iflag - k1*100)/10
       k3 = iflag - k1*100 - k2*10
       fname = 'data'//ch(k1)//ch(k2)//ch(k3)

        call output(nr,nz,r,z,u,v,p,ph,conc,t,fname)
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! update ph_grad
! ph_grad for further use
  do k=0,nz-1
    do i=0,nr-1
      !TODO: other boundary conditions
      ! Neumann bcd
        if (i .eq. 0) then
            !ph_grad(i,k) = ((ph(i+1,k) - ph(i-1+nr,k))**2) / (4.0 * dr2)
            ph_grad(i,k) = 0.0  !((ph(i+1,k) - ph(i-1+nr,k))**2) / (4.0 * dr2)
        elseif (i .eq. nr-1) then
            ph_grad(i,k) = 0.0  !((ph(i+1-nr,k) - ph(i-1,k))**2) / (4.0 * dr2)
        else
            ph_grad(i,k) = ((ph(i+1,k) - ph(i-1,k))**2) / (4.0 * dr2)
        end if

        if (k .eq. 0) then
            ph_grad(i,k) = ph_grad(i,k) + 0.0  !((ph(i,k+1) - ph(i,k-1+nz))**2) / (4.0 * dz2)
        elseif (k .eq. nz-1) then
            ph_grad(i,k) = ph_grad(i,k) + 0.0  !((ph(i,k+1-nz) - ph(i,k-1))**2) / (4.0 * dz2)
        else
            ph_grad(i,k) = ph_grad(i,k) + ((ph(i,k+1) - ph(i,k-1))**2) / (4.0 * dz2)
        end if
      end do
  end do

  ph_grad = sqrt(ph_grad)*0.5



    ! first order scheme is used as setup
     call PhSolve(rmu_after,ph_after,t+dt)

     call diffusion_eqn(c_after,t)


!	  !!  Calculate the ph, doesn't change ph and velocity. Store the new ph in phi_after
!    call PhSolve_2(u_pre,v_pre,ph_pre,ph_after,rmu_after,t+dt)
!    !call ph_mu_solve(ph_after, rmu_after, t+dt)

  ! map \phi and \rum to centers


  ! Calculte rho and eta
  call interpolate_coef(ph_after)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  update the velocity without incompressibility constrain
     if (only_cahn_hillarid .ne. 'y') then

        !TODO: here only first order is used
       call VeloSolve_combined(psi1,u_after, v_after, rmu_after,ph_after,t+dt)

        u_pre = u
        v_pre = v

        u = u_after
        v = v_after

       !project velocity onto divegence free space
       call PreSolve_2(plan,plan_back,psi2)

       psi2 = psi2*re*gamma1/dt

       p = psi2 + p !- w3 - w4

      if (is_modify_u) then
        !----------------------------------------------------------------
        ! modify velocity to enforce zero velocity of boundary
!        do k=0,nz-1
!           do i=0,nr
!            tmp = (1.0 + (ph_after(i,k) + ph_after(i,k-1))/2.0) /2.0
!            if ( tmp > 0.5 ) then
!              u(i,k) = 0.0
!            else
!              u(i,k) = u(i,k) * (1.0 - 2.0*tmp)
!            endif
!           end do
!        enddo
!
!        do k=1,nz-1
!           do i=0,nr-1
!            tmp = (1.0 + (ph_after(i,k) + ph_after(i,k-1))/2.0) / 2.0
!            if ( tmp > 0.5 ) then
!              v(i,k) = 0.0
!            else
!              v(i,k) = v(i,k) * (1.0 - 2.0*tmp)
!            endif
!           end do
!        enddo

        do k=0,nz-1
           do i=0,nr
            tmp = (1.0 + (ph_after(i,k) + ph_after(i-1,k))/2.0) /2.0
            if ( tmp > 0.95 ) then
              u(i,k) = 0.0
            else
              u(i,k) = u(i,k) ! * (1.0 - 2.0*tmp)
            endif
           end do
        enddo

        do k=1,nz-1
           do i=0,nr-1
            tmp = (1.0 + (ph_after(i,k) + ph_after(i,k-1))/2.0) / 2.0
            if ( tmp > 0.95 ) then
              v(i,k) = 0.0
            else
              v(i,k) = v(i,k) !* (1.0 - 2.0*tmp)
            endif
           end do
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
      endif ! end of is_modify_u



       !TODO: not correct velocity? in PreSolve_2


!       call VeloSolve_2_combined(psi1,psi2,u_pre,u_after,v_pre,v_after,ph_pre,ph_after,rmu_pre,rmu_after,t+dt)
!        u_pre = u
!        v_pre = v
!
!        u = u_after
!        v = v_after

     endif ! end of only_cahn_hillarid


    ph_pre = ph
    rmu_pre = rmu

    ph = ph_after
    rmu = rmu_after
    conc = c_after


! second order's part
!    !====================================================================
!    !  project veolicity on the divegence free space
!    !====================================================================
!
!    psi1 = psi2
!
!     if (only_cahn_hillarid .ne. 'y') then
!      call PreSolve_2(plan,plan_back,psi2)
!
!      psi2 = psi2*3*re*gamma1/2./dt
!
!
!      !  ! divergence of velocity,, if this is used, Flag: pressure ghost is needed
!      !  do k=0,nz-1
!      !     do i=0,nr-1
!      !      w3(i,k) = ( u(i+1,k)-u(i,k) )/dr
!      !      w4(i,k) = ( v(i,k+1)-v(i,k) )/dz
!      !     enddo
!      !  enddo
!
!
!      ! here whether we need w3*eta - w4*eta
!      p = psi2 + p !- w3*eta - w4*eta
!
!      !  Flag: pressure ghost, when use correction(w3*eta - w4*eta), we need ghost points for p
!      !   do i=0,nr
!      !    p(i,-1) = 2*p(i,0) - p(i,1)
!      !    p(i,nz+1) = 2*p(i,nz) - p(i,nz-1)
!      !   enddo
!    endif


     t = t+dt
     kcount = kcount + 1

    !     call bcd(t)
    !!!!!!!!!!!!!!!!!!!!!!!!!!! end of forward Euler


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  output the data and carry out the stopping test by stability of dissipation

    !       call dissipation(ph_pre,dissip2,mass2,t)
    !       write(22,*) kcount,dt*dissip2,mass2

     if(kcount/istep*istep == kcount) then

       print *,'maximum:  ', MAXVAL(MAXVAL(ph,1))
       print *,''


!!!!  test of cpu time per timestep
    call cpu_time(time_end)
    write(*,"('cpu time per time step  ',f10.6)"),(time_end-time_begin)/istep
    !pause

       iflag=iflag+1
       k1 = iflag/100
       k2 = (iflag - k1*100)/10
       k3 = iflag - k1*100 - k2*10
       fname = 'data'//ch(k1)//ch(k2)//ch(k3)


        call output(nr,nz,r,z,u,v,p,ph,conc,t,fname)

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

if(iflag .ge. 100) stop

!       if(k3+k2*10 .eq. 7) stop

       dissip1 = dissip2
       mass1 = mass2
    endif

enddo

 end program
