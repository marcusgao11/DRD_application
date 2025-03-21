! initialize the code: grids, initial conditions for u v p \phi \rmu (including ghost points)
!                      alpha_s, contact angle, time step

    subroutine finit(time_step_change)
    use global
    implicit none
    integer i,j
    real t1,dt1,dt2,dt3,dt4
    real time_step_change
    real,external::dist
    real cal_alpha_s

    dr=slen/nr
    dz=shig/nz
    dr2=dr**2
    dz2=dz**2



!    call random_seed()
    do j=0,nz
      do i=0,nr
         r(i,j)=(i+0.5)*dr  ! postion of x at cell center
         z(i,j)=(j+0.5)*dz
         u(i,j)=0.0
         v(i,j)=0.0
         p(i,j)=0.0
!         ph(i,j) = 1.0
!          ph(i,j)= tanh( (r(i,j)-0.5*slen)/(sqrt(2.0)*I_thickness) )
!          ph(i,j) = tanh((sqrt((r(i,j)-0.5*slen)**2/1.5 + (z(i,j)-shig*0.5)**2) - 0.3 )/I_thickness )
!         ph(i,j) = tanh((sqrt((r(i,j)-0.5*slen)**2 + (z(i,j)-shig*0.5)**2) - 0.3 )/(sqrt(2.0)*I_thickness) )
         !ph(i,j) = tanh((sqrt((r(i,j)-0.5*slen)**2 + z(i,j)**2) - 0.2 )/(sqrt(2.0)*I_thickness) )
         ph(i,j) = tanh( -1.0*dist(r(i,j), z(i,j)) /sqrt(2.0)/I_thickness )
!         ph(i,j) = -tanh((sqrt((r(i,j)-0.5*slen)**2 + (z(i,j)-shig*0.1)**2) - 0.1 )/(sqrt(2.0)*I_thickness) )
         !ph(i,j) = -tanh((sqrt((r(i,j)-0.5*slen)**2/2. + (z(i,j))**2) - 0.04 )/I_thickness )
 	!ph(i,j)= 0.01*sin(2*pi*r(i,j)/slen)*cos(pi*z(i,j)/shig)
	!ph(i,j)=1.0E-6
        !  call random_number(dt2)
        !  ph(i,j) =  (dt2*2 - 1)*10**(-6.)
      enddo
    enddo


!! for test
!    do j=0,nz
!      do i=0,nr
!        if ( z(i,j) > shig/2.0 ) then
!          ph(i,j) = 1.0 ! large value
!        else
!          ph(i,j) = -1.0
!        endif
!      enddo
!    enddo



!!!!!!!!!!!!! contact angle
     do i=0,nr
        angleB(i) = angle_bot
        !angleT(i) = 180 - angle
        angleT(i) = angle_top
     enddo


   call bcd(0.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  alpha_s is stability coefficient at boundary
     cal_alpha_s = sqrt(2.0)*abs(cos(angle*pi/180.0))/3.0*(pi/2.)**2/2.
     write(*, *) 'cal_alpha_s, ', cal_alpha_s
     write(*,*)
     !alpha_s = 0.0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dt1 = min(dr,dz)**2.*re/15.
	print *,'Re',dt1
    dt2 = min(dr,dz)**4./xld/(I_thickness**2)/25
	print *,'xld',dt2
    dt3 = min(dr,dz)/60.
	print *,'CFL',dt3
	print *,''

!    dt1 = min(dr,dz)**2.*re/4.   !!!     time step size
!	print *,'Re',dt1
!    dt2 = min(dr,dz)**4./xld/(I_thickness**2)/16.   !!!     time step size
!	print *,'xld',dt2
    dt3 = min(dr,dz)!/0.2 !!!     time step size
	print *,'CFL',dt3
!    dt4 = min(dr,dz)**2/ca/I_thickness**2  !!!     time step size
!	print *,'\mu\grad\phi/(I_thickness**2)',dt4
!	print *,''

    dt = dt3*time_step_change !/10.

    !dt = 0.035*min(dr,dz)**2/xld/4./2./2./2.*100

    print *,'times large than explicit scheme',dt/(min(dt1,dt2,dt3))
    print *,'time step of explicit',min(dt1,dt2,dt3)

    end subroutine
