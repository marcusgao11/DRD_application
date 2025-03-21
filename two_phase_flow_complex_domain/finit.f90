! initialize the code: grids, initial conditions for u v p \phi \rmu (including ghost points)
!                      alpha_s, contact angle, time step

    subroutine finit(time_step_change)
    use global
    use particles_move
    implicit none
    integer i,j
    real t1,dt1,dt2,dt3,dt4
    real time_step_change
    real,external::dist
    real h0,r0
    real temp_h, tmp_d
    real ep_s ! used for boundary interface thickness

    ! for initialization for phi
    !real tmp_phi(1:nz,1:nr)
    integer i_dt
    character lr_bd_type_bp, top_type_ph_bp, bot_type_ph_bp
    real tmp_phi(1:nr,1:nz)
    !real tmp_phi(1:3,1:2)
    real ph_after(-1:nr,-1:nz)
    real rmu_after(-1:nr,-1:nz)

    ! for square pattern
    real square_a, square_b, square_d
    real square_start_x, square_start_y
    integer square_index, max_square_num
    real dd, x_pos
  
    ! for sine wave pattern
    real sine_start_x, sine_b, sine_a


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
         !ph(i,j) = tanh( -1.0*dist(r(i,j), z(i,j)) /sqrt(2.0)/I_thickness )
!         ph(i,j) = -tanh((sqrt((r(i,j)-0.5*slen)**2 + (z(i,j)-shig*0.1)**2) - 0.1 )/(sqrt(2.0)*I_thickness) )
         !ph(i,j) = -tanh((sqrt((r(i,j)-0.5*slen)**2/2. + (z(i,j))**2) - 0.04 )/I_thickness )
 	!ph(i,j)= 0.01*sin(2*pi*r(i,j)/slen)*cos(pi*z(i,j)/shig)
	!ph(i,j)=1.0E-6
        !  call random_number(dt2)
        !  ph(i,j) =  (dt2*2 - 1)*10**(-6.)
      enddo
    enddo


!!--------------------------------------------------------------------------
! initialize boundary phi
! boundary phi: used for boundary and not updated with velocity,
! which can be defined or input from files

    ! TODO
! when initialize bdr_ph, needs to use ep_s
    ep_s = I_thickness/4.0 !*2 !/5.0

!!--------------------------------------------------------------------------
! type 1:
!! sin wave
!    do j=0,nz
!      do i=0,nr
!
!        temp_h = shig/3.0
!
!        if (z(i,j) > (shig - temp_h)) then
!          tmp_d = -(z(i,j) - (shig - temp_h))
!        else if (z(i,j) > shig/2.0) then
!          tmp_d = (shig - temp_h) - z(i,j)
!        else
!          tmp_d = (z(i,j) - temp_h)
!        endif
!
!        bdr_ph(i,j) = tanh( -1.0*tmp_d /sqrt(2.0)/I_thickness )
!      enddo
!    enddo
!! re-scale to [0.0,1.0]
!bdr_ph = (bdr_ph + 1.0)/2.0


!--------------------------------------------------------------------------
! type 1.1: read from files and no need to initialize
!tmp_phi = 1.0
! bdr_upp = shig; bdr_low = 0.0 ! should be modified
!call output_var(nz,nr,tmp_phi,'ttt.bin')
!


! !--------------------------------------------------------------------------
! ! type 2: input from files

!     call read_var(nr,nz,tmp_phi,'domain_snake_420_180_dist.bin')
!     do j=0,nz-1
!       do i=0,nr-1
!          bdr_ph_dist(i,j) = tmp_phi(i+1,j+1)
!       enddo
!     enddo

!     ! when reading from files, boundary range should be specified
!     bdr_upp = 0.6725
!     bdr_low = 0.1725

! !call read_var(nr,nz,tmp_phi,'phi.bin')
! !call read_var(nr,nz,tmp_phi,'domain_snake_2.bin')
! !call read_var(nr,nz,tmp_phi,'domain_snake_200_100_2.bin')
! !call read_var(nr,nz,tmp_phi,'domain_snake_400_200_2.bin')

! call read_var(nr,nz,tmp_phi,"domain_snake_420_180.bin")

! !write(*,*) "ttt", MAXVAL(MAXVAL(tmp_phi, 1))
! !pause

! !call read_var(3,2,tmp_phi,'phi.bin')
! !write(*,*) "ttt", MAXVAL(MAXVAL(tmp_phi, 1))
! !print *, tmp_phi
! !pause

! do j=0,nz-1
!   do i=0,nr-1
!      bdr_ph(i,j) = tmp_phi(i+1,j+1)
!   enddo
! enddo


! ! initialize by solve CH equation for several steps for input structure
! ph = bdr_ph
! do j=0,nz-1
!  ph(-1,j) = ph(0,j)
!  ph(nr,j) = ph(nr-1,j)
! enddo

! do i=-1,nr
!  ph(i,-1) = ph(i,0)
!  ph(i,nz) = ph(i,nz-1)
! end do

! ! store boundary conditions type and restore it after initialization
! lr_bd_type_bp = lr_bd_type
! top_type_ph_bp = top_type_ph
! bot_type_ph_bp = bot_type_ph

! lr_bd_type = 'n'
! top_type_ph = 'n'
! bot_type_ph = 'n'
! dt = 1.0E-4


! !=================================================================
! ! need give value for initialization
! ph_M_0 = 1.0
! ph_M_inner = 1.0 ! mobility in CH equation


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  setup for the scheme
! call setup()

! do j=0,nz-1
!    do i=0,nr-1
!       rmu(i,j) = -I_thickness*( (ph(i+1,j)+ph(i-1,j)-2*ph(i,j))/dr2 &
!         + (ph(i,j+1)+ph(i,j-1)-2*ph(i,j))/dz2 )&
!        & - ph(i,j)/I_thickness + ph(i,j)**3/I_thickness
!    enddo
! enddo

! write(*,*) 'initialize boundary phi by CH equation'
! do i_dt=0,50
!   call PhSolve_ori(rmu_after,ph_after,0.0)
!   rmu = rmu_after
!   ph = ph_after
! enddo

! bdr_ph = (ph + 1.0)/2.0

! ! truncate
! do j=0,nz-1
!    do i=0,nr-1
!       if (bdr_ph(i,j) >= 1.0) then
!         bdr_ph(i,j) = 1.0
!       else if (bdr_ph(i,j) <= 0.0) then
!         bdr_ph(i,j) = 0.0
!       endif
!    enddo
! enddo

! ! store boundary conditions type and restore it after initialization
! lr_bd_type = lr_bd_type_bp
! top_type_ph = top_type_ph_bp
! bot_type_ph = bot_type_ph_bp
! ! end of initialize boundary phi
! !--------------------------------------------------------------------------


! !--------------------------------------------------------------------------
! ! type 2: input from files

!     call read_var(nr,nz,tmp_phi,"water_tank_400_300.bin")
!     do j=0,nz-1
!       do i=0,nr-1
!          bdr_ph_dist(i,j) = tmp_phi(i+1,j+1)
!       enddo
!     enddo

!     ! when reading from files, boundary range should be specified
!     bdr_upp = 3.0
!     bdr_low = 0.5

! !write(*,*) "ttt", MAXVAL(MAXVAL(tmp_phi, 1))
! !pause

! !call read_var(3,2,tmp_phi,'phi.bin')
! !write(*,*) "ttt", MAXVAL(MAXVAL(tmp_phi, 1))
! !print *, tmp_phi
! !pause

! do j=0,nz-1
!   do i=0,nr-1
!      bdr_ph(i,j) = tmp_phi(i+1,j+1)
!   enddo
! enddo

! call reinint_by_PhSolve()



! ! !--------------------------------------------------------------------------
! ! type 3: no boundary
! bdr_upp = shig; bdr_low = 0.0
! bdr_ph = 0.0

! do j=-1,nz
!   do i=-1,nr
!     if (bdr_ph(i,j) > 0.5) then
!       bdr_ph_penalty(i,j) = 1.0
!     else
!       bdr_ph_penalty(i,j) = 0.0
!     endif
!   enddo
! enddo 


!--------------------------------------------------------------------------
!type 4:  half plane
bdr_upp = shig; bdr_low = shig/5.0;
   do j=0,nz
     do i=0,nr
       bdr_ph(i,j) = tanh( -1.0*(z(i,j)-bdr_low)/sqrt(2.0)/ep_s )
     enddo
   enddo

  bdr_ph = (bdr_ph + 1.0)/2.0
  do j=0,nz-1
    bdr_ph(-1,j) = bdr_ph(0,j)
    bdr_ph(nr,j) = bdr_ph(nr-1,j)
  enddo

  do i=-1,nr
    bdr_ph(i,-1) = bdr_ph(i,0)
    bdr_ph(i,nz) = bdr_ph(i,nz-1)
  end do

  do j=-1,nz
    do i=-1,nr
      if (bdr_ph(i,j) > 0.5) then
        bdr_ph_penalty(i,j) = 1.0
      else
        bdr_ph_penalty(i,j) = 0.0
      endif
    enddo
  enddo 


! !--------------------------------------------------------------------------
!   ! type 5: sandwitch
!    bdr_upp = shig*4/5.0; bdr_low = shig*1.0/5.0;
!    !bdr_upp = shig*3.0/4.0; bdr_low = shig/4.0;
!    do j=0,nz
!      do i=0,nr
!        bdr_ph(i,j) = max( tanh( -1.0*(z(i,j)-bdr_low)/sqrt(2.0)/ep_s ), &
!         & tanh( (z(i,j)-bdr_upp)/sqrt(2.0)/ep_s ) )
!      enddo
!    enddo


!   bdr_ph = (bdr_ph + 1.0)/2.0
!   do j=0,nz-1
!     bdr_ph(-1,j) = bdr_ph(0,j)
!     bdr_ph(nr,j) = bdr_ph(nr-1,j)
!   enddo

!   do i=-1,nr
!     bdr_ph(i,-1) = bdr_ph(i,0)
!     bdr_ph(i,nz) = bdr_ph(i,nz-1)
!   end do


!   do j=-1,nz
!     do i=-1,nr
!       if (bdr_ph(i,j) > 0.5) then
!         bdr_ph_penalty(i,j) = 1.0
!       else
!         bdr_ph_penalty(i,j) = 0.0
!       endif
!     enddo
!   enddo


  ! !--------------------------------------------------------------------------
  ! ! type 6: square pattern
  ! ! square
  ! square_a = 0.1*slen; square_b = 0.1*slen; square_d = 0.2*shig;
  ! square_start_x = 0.1*slen + square_b/2.0;  square_start_y = 0.1*shig; 

  ! bdr_upp = shig; bdr_low = square_start_y
  ! max_square_num = 4

  ! do j=0,nz
  !   do i=0,nr
  !     square_index = floor( (r(i,j)-square_start_x)/(square_a+square_b) ) + 1
  !     if (square_index > max_square_num) then
  !       square_index = max_square_num + 1
  !     elseif ( square_index < 1) then
  !       square_index = 0
  !     endif

  !     !print *, r(i,j), square_index

  !     x_pos = square_start_x + (square_a+square_b)/2.0*(2.0*square_index-1)
  !     if (square_index > max_square_num) then
  !       dd = max( -(r(i,j)-x_pos) - square_a/2.0 ,  (z(i,j)-square_d/2.0) - square_d/2.0)
  !     elseif ( square_index < 1) then
  !       dd = max( (r(i,j)-x_pos) - square_a/2.0,  (z(i,j)-square_d/2.0) - square_d/2.0)
  !     else
  !       dd = max( abs(r(i,j)-x_pos) - square_a/2.0,  (z(i,j)-square_d/2.0) - square_d/2.0)
  !     endif

  !     bdr_ph(i,j) = tanh( - dd/sqrt(2.0)/ep_s )
  !   enddo
  ! enddo

  ! bdr_ph = (bdr_ph + 1.0)/2.0
  ! do j=0,nz-1
  !   bdr_ph(-1,j) = bdr_ph(0,j)
  !   bdr_ph(nr,j) = bdr_ph(nr-1,j)
  ! enddo

  ! do i=-1,nr
  !   bdr_ph(i,-1) = bdr_ph(i,0)
  !   bdr_ph(i,nz) = bdr_ph(i,nz-1)
  ! end do


  ! do j=-1,nz
  !   do i=-1,nr
  !     if (bdr_ph(i,j) > 0.5) then
  !       bdr_ph_penalty(i,j) = 1.0
  !     else
  !       bdr_ph_penalty(i,j) = 0.0
  !     endif
  !   enddo
  ! enddo


  ! !--------------------------------------------------------------------------
  ! ! type 6: sine wave pattern
  ! sine_start_x = 0.1*slen;
  ! sine_b = 0.1*shig;
  ! sine_a = 0.1*slen;

  ! bdr_upp = shig; bdr_low = sine_b;

  ! do j=0,nz
  !   do i=0,nr
  !     if (r(i,j) <= sine_start_x .OR. r(i,j) >= slen - sine_start_x) then
  !       if (z(i,j) <= sine_b) then
  !         bdr_ph(i,j) = 1.0
  !       else
  !         bdr_ph(i,j) = -1.0
  !       endif
  !     elseif ( z(i,j) <= sine_b*( sin(PI/sine_a * (r(i,j)-sine_start_x) ) + 1.0 ) ) then
  !       bdr_ph(i,j) = 1.0
  !     else
  !       bdr_ph(i,j) = -1.0
  !     endif     
  !   enddo
  ! enddo

  ! call reinint_by_PhSolve()


  ! do j=-1,nz
  !   do i=-1,nr
  !     if (bdr_ph(i,j) > 0.5) then
  !       bdr_ph_penalty(i,j) = 1.0
  !     else
  !       bdr_ph_penalty(i,j) = 0.0
  !     endif
  !   enddo
  ! enddo


!--------------------------------------------------------------------------
! initial for particle phi
do j=0,nz-1
  do i=0,nr-1
    fpd_ph(i,j) = tanh( -1.0*dist(r(i,j), z(i,j)) /sqrt(2.0)/I_thickness )
  enddo
enddo


! re-scale to [0.0,1.0]
fpd_ph = (fpd_ph + 1.0)/2.0

! incorporate boundary phi
fpd_ph = fpd_ph + bdr_ph


!-----------------------------------------------------------------------------------
! ! -------------------------------------
! ! initial for phi
! ! cannot put this above, since the variale ph is temp used
!   do j=0,nz
!     do i=0,nr
!         !ph(i,j) = tanh((sqrt((r(i,j)-0.5*slen)**2/2.0 + (z(i,j)-shig*0.5)**2) - 0.15 )/(sqrt(2.0)*I_thickness) )
!         ph(i,j) = -tanh((sqrt((r(i,j)-0.5*slen)**2 + (z(i,j)-shig*0.5)**2) - 0.2 )/(sqrt(2.0)*I_thickness) )
!         !ph(i,j) = (tanh(-(sqrt((r(i,j)-0.5*slen)**2 + (z(i,j)-0.1)**2) - 0.6 )/(sqrt(2.0)*I_thickness) ) + 1.0)/2.0
!         !ph(i,j) = (tanh(-(sqrt((r(i,j)-0.5*slen)**2/2.0 + (z(i,j)-0.7)**2) - 0.3 )/(sqrt(2.0)*I_thickness) ) + 1.0)/2.0
!         !ph(i,j) = tanh(-(sqrt((r(i,j)-0.5*slen)**2 + (z(i,j)-bdr_low)**2) - 0.25 )/(sqrt(2.0)*I_thickness) )
!     enddo
!   enddo

  ! ! sandiwtich
  ! do j=0,nz
  !   do i=0,nr
  !       !  ph(i,j) = min( tanh( ( r(i,j)-0.1*slen )/(sqrt(2.0)*I_thickness) ), &
  !       !    &  tanh( - ( r(i,j)-0.2*slen )/(sqrt(2.0)*I_thickness) ) )
  !          ph(i,j) = min( tanh( ( r(i,j)-0.4*slen )/(sqrt(2.0)*I_thickness) ), &
  !          &  tanh( - ( r(i,j)-0.6*slen )/(sqrt(2.0)*I_thickness) ) )
  !     enddo 
  ! enddo

  ! square
  do j=0,nz
    do i=0,nr
      ! ph(i,j) = tanh( - max( abs(r(i,j)-slen/2.0) - 1.5, abs(z(i,j)- (bdr_low + 0.2/2) ) - 0.2 ) &
      !   & / (sqrt(2.0)*I_thickness) )
      ! ph(i,j) = tanh( -max( abs(r(i,j)-slen/2.0) - 0.6, abs(z(i,j)- (bdr_low + 0.2/2) ) - 0.05 ) &
      !   & / (sqrt(2.0)*I_thickness) )

      ! ph(i,j) = tanh(-(sqrt((r(i,j)-0.5*slen)**2/2.0 + (z(i,j)-0.5*shig)**2) - 0.1 )&
      !  &/(sqrt(2.0)*I_thickness) )

       ph(i,j) = tanh(-(sqrt((r(i,j)-0.5*slen)**2 + (z(i,j)-bdr_low)**2) - 0.2 )&
         &/(sqrt(2.0)*I_thickness) )

      ! ph(i,j) = tanh(-(sqrt((r(i,j)-0.5*slen)**2 + (z(i,j)-0.8)**2) - 0.4 )&
      !    &/(sqrt(2.0)*I_thickness) )

      !ph(i,j) = tanh( -(r(i,j)-0.5*slen) /(sqrt(2.0)*I_thickness) )

    enddo 
  enddo

  if (phi_mode == 1) then
    ph = (ph + 1.0)/2.0
  endif

  
  ! substract bdr
  !ph = min(1.0 - bdr_ph, ph)
  ph = (1.0 - bdr_ph) * ph

  ! ! sandwitch
  !  do j=0,nz
  !    do i=0,nrbdr_ph
  !     tmp_d = min(-(r(i,j) - 0.5),  (r(i,j) - 0.2))

  !      ph(i,j) = tanh( 1.0*tmp_d /sqrt(2.0)/I_thickness )
  !    enddo
  !  enddo
  !  ! re-scale to [0.0,1.0]
  !  ph = (ph + 1.0)/2.0


  !TODO: only support bdrs
  do j=0,nz-1
    ph(-1,j) = ph(0,j)
    ph(nr,j) = ph(nr-1,j)
  enddo

  do i=-1,nr
    ph(i,-1) = ph(i,0)
    ph(i,nz) = ph(i,nz-1)
  end do


  do j=0,nz-1
    do i=0,nr-1
        rmu(i,j) = -I_thickness*( (ph(i+1,j)+ph(i-1,j)-2*ph(i,j))/dr2 &
            + (ph(i,j+1)+ph(i,j-1)-2*ph(i,j))/dz2 )&
        & - ph(i,j)/I_thickness + ph(i,j)**3/I_thickness
    enddo
  enddo
  do i=0,nr-1
    rmu(i,-1)=rmu(i,0)
    rmu(i,nz)=rmu(i,nz-1)
 enddo
   do j=-1,nz
    rmu(-1,j) = rmu(0,j)
    rmu(nr,j) = rmu(nr-1,j)
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


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! set mobiliy and bdr_ph_grd for enforcing bdrs
   bdr_ph_grad = 0.0
   do j=0,nz-1
    do i=0,nr-1
        bdr_ph_grad(i,j) = sqrt( ( (bdr_ph(i+1,j)-bdr_ph(i-1,j))/2.0/dr )**2.0 &
          & + ( (bdr_ph(i,j+1)-bdr_ph(i,j-1))/2.0/dz )**2.0 )
     enddo
   enddo

    !ph_M_0 = 1.0
  !  ph_M_inner = 1.0 - bdr_ph ! mobility in CH equation

    !ph_M_0 = max(1.0 - bdr_ph, 1.0E-3) ! 1.0 - bdr_ph 
    
    ph_M_0 = 1.0 - bdr_ph  ! 1.0 - bdr_ph 

      ! do j=-1,nz
      !   do i=-1,nr
      !     if (bdr_ph(i,j) > 0.5) then
      !       ph_M_0(i,j) = 0.0
      !     else
      !       ph_M_0(i,j) = 1.0
      !     endif
      !   enddo
      ! enddo

    ph_M_inner = 1.0 - bdr_ph ! mobility in CH equation

    bdr_ph_used = 1.0 - bdr_ph
    bdr_ph_used_adjust = 1.0E-6

    ! bdr_ph_used = 1.0
    ! bdr_ph_used_adjust = 0.0


!!!!!!!!!!!!! contact angle
     do i=0,nr
        angleB(i) = angle
        !angleT(i) = 180 - angle
        angleT(i) = angle
     enddo


   call bcd(0.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  alpha_s is stability coefficient at boundary
     alpha_s = sqrt(2.0)*abs(cos(angle*pi/180.0))/3.0*(pi/2.)**2/2.
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

    !dt = dt3*time_step_change !/10.

    dt = 1.0E-4/2.0/8 !5.0E-4

    !dt = 0.035*min(dr,dz)**2/xld/4./2./2./2.*100

    print *,'times large than explicit scheme',dt/(min(dt1,dt2,dt3))
    print *,'time step of explicit',min(dt1,dt2,dt3)

    end subroutine
