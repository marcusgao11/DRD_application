module particles_move

  contains

  ! subroutine Particles_Force_Velocity

  ! end subroutine

  subroutine particle_initialize(xc, yc, lx, ly, is_fluid) !N, l0)
    use particles_params
    use particles
    implicit none
    integer i, i_num
    real angle, xc, yc
    logical is_fluid
    real lx,ly
    real tmp, tmp_pos(0:1)
    real di_x,di_y,d
    integer i_count
    logical is_skip

    ! active force
    active_force   = 15*3/1.7 !5 !20.0 !1.0  !2.1250 !10 !0.
    if (.not. is_fluid) then
      active_force = 1.0
    endif

    active_force_ang   = 0.0 !PI/2

    call random_seed()

    !=============================================================
    ! fiber
    if (Np > 2) then
      i_num = floor(Np/2.0)
      points(i_num, 0) = lx/2.0  ! i_num = 7
      points(i_num, 1) = ly/2.0

      do i = 0, i_num-1
        points(i, 0) = lx/2.0 
        points(i, 1) = ly/2.0 + (i_num-i) * l0 !* 1.2

        points(Np-i-1, 0) = lx/2.0
        points(Np-i-1, 1) = ly/2.0 - (i_num-i) * l0 !* 1.2

        ! points(i, 0) = lx/2.0 
        ! points(i, 1) = ly/2.0 - (i_num-i) * l0

        ! points(Np-i-1, 0) = lx/2.0
        ! points(Np-i-1, 1) = ly/2.0 + (i_num-i) * l0
      enddo
    endif


!=============================================================
! swimmers

!    do i = 0, Np-1
!        angle = 2.0 * PI * i / Np
!        points(i, 0) = xc + 0.8 * cos(angle)
!        points(i, 1) = yc + 1.35 * sin(angle)
!    enddo


!    ! test case 1
!    points(:, 0) = 0.5
!    points(4, 1) = 0.5
!
!    ! Update particle positions
!    points(0, 1) = points(4, 1) + 8.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(1, 1) = points(4, 1) + 6.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(2, 1) = points(4, 1) + 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(3, 1) = points(4, 1) + 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(4, 1) = 0.5
!    points(5, 1) = points(4, 1) - 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(6, 1) = points(4, 1) - 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(7, 1) = points(4, 1) - 6.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(8, 1) = points(4, 1) - 8.0 * 2.0 * (fpd_radius + fpd_inter_dist)


!    ! test case 2
!    points(4, 0) = 0.5
!    points(0, 0) = points(4, 0) + 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(1, 0) = points(4, 0) + 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(2, 0) = points(4, 0) + 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(3, 0) = points(4, 0) + 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(4, 0) = 0.5
!    points(5, 0) = points(4, 0) - 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(6, 0) = points(4, 0) - 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(7, 0) = points(4, 0) - 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(8, 0) = points(4, 0) - 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!
!    points(4, 1) = 0.5
!    points(0, 1) = points(4, 1) + 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(1, 1) = points(4, 1) + 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(2, 1) = points(4, 1) + 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(3, 1) = points(4, 1) + 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(4, 1) = 0.5
!    points(5, 1) = points(4, 1) - 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(6, 1) = points(4, 1) - 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(7, 1) = points(4, 1) - 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!    points(8, 1) = points(4, 1) - 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)


!  !!string bending test 1
!  points(4, 0) = 0.5
!  points(0, 0) = points(4, 0) - 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(1, 0) = points(4, 0) - 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(2, 0) = points(4, 0) - 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(3, 0) = points(4, 0) - 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(4, 0) = 0.5
!  points(5, 0) = points(4, 0) - 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(6, 0) = points(4, 0) - 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(7, 0) = points(4, 0) - 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(8, 0) = points(4, 0) - 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!
!  points(4, 1) = 0.5
!  points(0, 1) = points(4, 1) + 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(1, 1) = points(4, 1) + 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(2, 1) = points(4, 1) + 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(3, 1) = points(4, 1) + 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(4, 1) = 0.5
!  points(5, 1) = points(4, 1) - 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(6, 1) = points(4, 1) - 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(7, 1) = points(4, 1) - 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
!  points(8, 1) = points(4, 1) - 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)


  ! !!!string bending test 2  Np=9
  ! points(4, 0) = 0.5
  ! points(0, 0) = points(4, 0) + 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(1, 0) = points(4, 0) + 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(2, 0) = points(4, 0) + 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(3, 0) = points(4, 0) + 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(4, 0) = 0.5
  ! points(5, 0) = points(4, 0) + 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(6, 0) = points(4, 0) + 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(7, 0) = points(4, 0) + 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(8, 0) = points(4, 0) + 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)

  ! points(4, 1) = 0.5
  ! points(0, 1) = points(4, 1) + 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(1, 1) = points(4, 1) + 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(2, 1) = points(4, 1) + 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(3, 1) = points(4, 1) + 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(4, 1) = 0.5
  ! points(5, 1) = points(4, 1) - 1.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(6, 1) = points(4, 1) - 2.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(7, 1) = points(4, 1) - 3.0 * 2.0 * (fpd_radius + fpd_inter_dist)
  ! points(8, 1) = points(4, 1) - 4.0 * 2.0 * (fpd_radius + fpd_inter_dist)


  !   do i = 0, Np_swm-1
  !       pos_swm(i, 0) = 0.3 !xc !1.6 + 0.15 !0.1*i
  !       pos_swm(i, 1) = yc
  !   enddo


    !!===================================
!    ! for test fpd_dist_periodic
!    do i = 0, Np_swm-1
!!        pos_swm(i, 0) = 0.18 !xc + 0.15 !0.1*i
!!        pos_swm(i, 1) = yc
!
!!        pos_swm(i, 0) = 0.05 !xc + 0.15 !0.1*i
!!        pos_swm(i, 1) = yc
!
!!        pos_swm(i, 0) = 0.05 !xc + 0.15 !0.1*i
!!        pos_swm(i, 1) = 0.05
!
!!        pos_swm(i, 0) = 2.0 - 0.05 !xc + 0.15 !0.1*i
!!        pos_swm(i, 1) = yc
!
!!        pos_swm(i, 0) = 2.0 - 0.05 !xc + 0.15 !0.1*i
!!        pos_swm(i, 1) = 0.05
!
!!        pos_swm(i, 0) = 2.0 - 0.05 !xc + 0.15 !0.1*i
!!        pos_swm(i, 1) = 1.0 - 0.05
!
!        pos_swm(i, 0) = xc + 0.15 !0.1*i
!        pos_swm(i, 1) = yc
!    enddo
!
!
!  points(0, 0) = 2.0 - 0.05 !xc + 0.15 !0.1*i
!  points(1, 0) = 1.0 - 0.05

    ! !===================================
    !test case 0:
    !i = 0
    !pos_swm(i, 0) = 0.7 !lx - 0.2 !/2.0 ! - 0.1 !/ 2.0 ! - 0.2 
    !pos_swm(i, 1) = ly/2 !+ 0.1 !/ 2.0
    !active_F_swm_ang(i) = 0
    

     !i = 0
     !pos_swm(i, 0) = lx/2 !lx - 0.2 !/2.0 ! - 0.1 !/ 2.0 ! - 0.2 
     !pos_swm(i, 1) = 0.7 !ly/2 !+ 0.1 !/ 2.0
     !active_F_swm_ang(i) = PI/2


    ! two paralel 
    !i = 0
      !pos_swm(i, 0) = lx - 0.2 !/2.0 ! - 0.1 !/ 2.0 ! - 0.2 
      !pos_swm(i, 1) = ly/2 + 0.1 !/ 2.0
    !active_F_swm_ang(i) = PI

    !i = 1
    !pos_swm(i, 0) = 0 + 0.2 !/2.0 ! - 0.1 !/ 2.0 ! - 0.2 
    !pos_swm(i, 1) = ly/2 - 0.1 !/ 2.0
    !active_F_swm_ang(i) = 0

    ! !test case 1
    ! !swimmer
    ! do i = 0, Np_swm-1
    !   call random_number(tmp)
    !   pos_swm(i, 0) = lx * tmp 
      
    !   call random_number(tmp)
    !   pos_swm(i, 1) = ly * tmp
    ! enddo

	! case:
    ! generate swimmers considering not interact exisiting fibers and swimmers
	i_count = 0
	do while (i_count < Np_swm)
	  is_skip = .false.
	  call random_number(tmp); tmp_pos(0) = lx * tmp 
	  call random_number(tmp); tmp_pos(1) = ly * tmp

	  ! check dist to fibers
	  do i = 0, Np-1
		di_x = MirrorDistance(tmp_pos(0) - points(i, 0), lx)
		di_y = MirrorDistance(tmp_pos(1) - points(i, 1), ly)
		d = sqrt( di_x**2 + di_y**2 )
		if (d < fpd_swm_dist) then ! regenerate
		  is_skip = .true.
		  exit
		endif
	  enddo
	  if (is_skip) then; CYCLE; endif
	  
	  ! check dist to existing swimmers
	  do i = 0, i_count-1
		di_x = MirrorDistance(tmp_pos(0) - pos_swm(i, 0), lx)
		di_y = MirrorDistance(tmp_pos(1) - pos_swm(i, 1), ly)
		d = sqrt( di_x**2 + di_y**2 )
		if (d < swm_swm_dist) then ! regenerate
		  is_skip = .true.
		  exit
		endif
	  enddo     
	  if (is_skip) then; CYCLE; endif

	  pos_swm(i_count,:) = tmp_pos
	  i_count = i_count + 1    
	enddo
    
    ! case:
    ! generate random swimmers
    ! do i = 0, Np_swm-1
    !   call random_number(tmp)
    !   pos_swm(i, 0) = lx * tmp 
      
    !   call random_number(tmp)
    !   pos_swm(i, 1) = ly * tmp
    ! enddo


    ! ! swimmer test
    ! ! two swimmers at left
    ! if (Np_swm > 0) then
    !   i = 0
    !   pos_swm(i, 0) = 0.8 !lx/2.0 - 0.1 !0.03
    !   pos_swm(i, 1) = ly/2.0 - 0.1
    ! endif

    ! if (Np_swm > 1) then
    !   i = 1
    !   pos_swm(i, 0) = 0.3 !lx/2.0 + 0.2
    !   pos_swm(i, 1) = 0.8 !ly/2.0 + 
    ! endif

    ! ! ! two swimmers at right
    ! ! i = 0
    ! ! pos_swm(i, 0) = lx/2.0 + 0.03
    ! ! pos_swm(i, 1) = ly/2.0

    ! ! i = 1
    ! ! pos_swm(i, 0) = lx/2.0 + 0.06
    ! ! pos_swm(i, 1) = ly/2.0 + 0.05*2.5


    write(*,*) "Np:", Np
    write(*,*) "Np_swm:", Np_swm

    ! initialize


    if (is_Tumbles_on) then
      swm_tP = 0.0            
    ! else
    !   do i = 0, Np_swm-1
    !     active_F_swm_ang(i) = active_force_ang
    !   enddo
    endif

    return
  end subroutine


  subroutine output_particles(fname) !N, l0)
    use particles_params
    use particles
    implicit none
    character fname*7
    integer i,k

    open(11,file=fname, status='unknown',form='unformatted')
    write(11)((points(i,k),i=0,Np-1),k=0,1)
    close(11)
    return
  end subroutine


  function fpd_dist_periodic(x,y, lx, ly, r_x, r_y, r)
  use particles
  real x,y, lx,ly
  real r_x, r_y, r

  real x_dist, y_dist

    ! for Phantoms, the position may be outside, needs to rescale
    if (r_x < 0.0) then
      r_x = r_x + lx
    elseif (r_x > lx) then
      r_x = r_x - lx
    endif

    if (r_y < 0.0) then
      r_y = r_y + ly
    elseif (r_y > ly) then
      r_y = r_y - ly
    endif


    ! justify where the the particle
    x_dist = 0.0; y_dist = 0.0;
    if (r_x < r) then
      x_dist = -lx
    else if (lx-r_x < r) then
      x_dist = lx
    endif

    if (r_y < r) then
      y_dist = -ly
    else if (ly-r_y < r) then
      y_dist = ly
    endif

    ! compare with four dist: TODO: here we can compare less than 4 points by conditions
    ! now, for simply, we calculate 4
    fpd_dist_periodic = sqrt((x-r_x)**2 + (y-r_y)**2) - r

    fpd_dist_periodic = min(fpd_dist_periodic,  sqrt((x+x_dist-r_x)**2 + (y-r_y)**2) - r)

    fpd_dist_periodic = min(fpd_dist_periodic,  sqrt((x-r_x)**2 + (y+y_dist-r_y)**2) - r)

    fpd_dist_periodic = min(fpd_dist_periodic,  sqrt((x+x_dist-r_x)**2 + (y+y_dist-r_y)**2) - r)

  end function fpd_dist_periodic



  function fpd_dist(x,y, lx, ly)
  use particles
  implicit none
  real x,y, fpd_dist
  real lx, ly
  integer i

  fpd_dist = 1.0E10
  do i = 0,Np-1
      fpd_dist = min(fpd_dist, fpd_dist_periodic(x,y, lx, ly, points(i, 0), points(i, 1), fpd_radius))
      !min(fpd_dist, sqrt((x-points(i, 0))**2 + (y-points(i, 1))**2) - fpd_radius)
  enddo

  ! swimmers
  do i = 0,Np_swm-1
      fpd_dist = min(fpd_dist, fpd_dist_periodic(x,y, lx, ly, pos_swm(i, 0), pos_swm(i, 1), swm_radius))
      !min(fpd_dist, sqrt((x-pos_swm(i, 0))**2 + (y-pos_swm(i, 1))**2) - swm_radius)
  enddo

  ! ! Phantoms: not considered as large viscosity part
  ! ! Phantoms
  ! do i = 0,Np_swm-1
  !     fpd_dist = min(fpd_dist, fpd_dist_periodic(x,y, lx, ly, &
  !         & pos_swm(i, 0) - swm_bp_dist*cos(active_F_swm_ang(i)), pos_swm(i, 1) - swm_bp_dist*sin(active_F_swm_ang(i)), swm_radius))
  !     !min(fpd_dist, sqrt( (x - (pos_swm(i, 0) - swm_bp_dist*cos(active_F_swm_ang(i))) )**2 &
  !     !      & + (y - (pos_swm(i, 1) - swm_bp_dist*sin(active_F_swm_ang(i))) )**2 ) - swm_radius )
  ! enddo

  end function fpd_dist


  function fpd_dist_swm(x,y, lx, ly)
    use particles
    implicit none
    real x,y, fpd_dist_swm
    real lx, ly
    integer i
  
    fpd_dist_swm = 1.0E10  
    ! swimmers
    do i = 0,Np_swm-1
      fpd_dist_swm = min(fpd_dist_swm, fpd_dist_periodic(x,y, lx, ly, pos_swm(i, 0), pos_swm(i, 1), swm_radius))
        !min(fpd_dist, sqrt((x-pos_swm(i, 0))**2 + (y-pos_swm(i, 1))**2) - swm_radius)
    enddo
  
    ! Phantoms: not considered as large viscosity part
  !  ! Phantoms
  !  do i = 0,Np_swm-1
  !      fpd_dist = min(fpd_dist, fpd_dist_periodic(x,y, lx, ly, &
  !          & pos_swm(i, 0) - swm_bp_dist*cos(active_F_swm_ang(i)), pos_swm(i, 1) - swm_bp_dist*sin(active_F_swm_ang(i)), swm_radius))
  !      !min(fpd_dist, sqrt( (x - (pos_swm(i, 0) - swm_bp_dist*cos(active_F_swm_ang(i))) )**2 &
  !      !      & + (y - (pos_swm(i, 1) - swm_bp_dist*sin(active_F_swm_ang(i))) )**2 ) - swm_radius )
  !  enddo
  
  
    end function fpd_dist_swm



  subroutine force_caculation(dt,t, lx, ly)
    use particles_params
    use particles
    implicit none
    real dt,t, lx,  ly
    logical,parameter::is_debug = .FALSE.

    !real :: tot_F(0:Np-1, 0:1),
    real :: Fs(0:Np-1, 0:1), Fb(0:Np-1, 0:1)
    real::rdiff(0:Np-1), tvect(0:Np-1, 0:1), nvect(0:Np-1, 0:1)
    real :: theta(0:Np-1), alpha(0:Np-1), G(0:Np-1)
    integer :: i,j
    real f(2)
    real disX, disY

    ! Initialize arrays
    ! forces from fibers
    tot_F_swm = 0.0

    tot_F = 0.0
    Fs = 0.0
    Fb = 0.0
    rdiff = 0.0
    tvect = 0.0
    nvect = 0.0
    theta = 0.0
    alpha = 0.0
    G = 0.0

    do i = 0, Np-1
        if (i < Np-1) then
          disX = MirrorDistance(points(i+1, 0)-points(i, 0), lx)
          disY = MirrorDistance(points(i+1, 1)-points(i, 1), ly)
          rdiff(i) = sqrt( disX**2 + disY**2 )
          tvect(i, :) = [disX/rdiff(i), disY/rdiff(i)]
          nvect(i, :) = [-disY/rdiff(i), disX/rdiff(i)]
          theta(i) = atan2(tvect(i, 1), tvect(i, 0))

          if (theta(i) < 0)  then
            theta(i) = theta(i) + 2*PI;
          endif

        end if
    end do

!    write(*,*) "rdiff: ", rdiff
!    write(*, *) "tvect: ", tvect
!
!    write(*,*) "nvect: ", nvect
!    write(*, *) "theta: ", theta
!    pause


    do i = 0, Np-1
        if (i < Np-2) then
            alpha(i) = theta(i+1) - theta(i)
            ! if (abs(alpha(i)) >= PI/2) then
            !     alpha(i) = alpha(i) + 2.0 * PI
            ! end if
            alpha(i) = MirrorDistance(alpha(i), 2*PI);

            G(i) = alpha(i) - Theta0
        end if
    end do

    do i = 0, Np-1
        if (i == 0) then
            !Fs(i, :) = ks * (rdiff(i) - l0) * tvect(i, :) - ks * (rdiff(Np-1) - l0) * tvect(Np-1, :)
            !Fb(i, :) = kb / l0 * ((G(Np-1) - G(Np-2)) / rdiff(Np-1) * nvect(Np-1, :) - (G(i) - G(Np-1)) / rdiff(i) * nvect(i, :))

            Fs(i, :) = ks/(l0*l0) * (rdiff(i) - l0) * tvect(i, :)
            Fb(i, :) = -kb * G(i) / rdiff(i) * nvect(i, :)
            cycle
        end if

        if (i == 1) then
            !Fs(i, :) = ks * (rdiff(i) - l0) * tvect(i, :) - ks * (rdiff(i-1) - l0) * tvect(i-1, :)
            !Fb(i, :) = kb / l0 * ((G(i-1) - G(Np-1)) / rdiff(i-1) * nvect(i-1, :) - (G(i) - G(i-1)) / rdiff(i) * nvect(i, :))
            Fs(i, :) = ks/(l0*l0) * ( (rdiff(i) - l0) * tvect(i, :) - (rdiff(i-1) - l0) * tvect(i-1, :) )
            Fb(i, :) = kb * ( G(i-1) / rdiff(i-1) * nvect(i-1, :) - (G(i) - G(i-1)) / rdiff(i) * nvect(i-1, :) )
            cycle
        end if

        if (i == Np-1) then
            Fs(i, :) = -ks/(l0*l0) *(rdiff(i-1)-l0)*tvect(i-1, :)
            Fb(i, :) = - kb *(G(i-2)/rdiff(i-1)*nvect(i-1,:))
            cycle
        endif

        if (i == Np - 2) then
            Fs(i, :) = ks/(l0*l0) * ( (rdiff(i) - l0) * tvect(i, :) - (rdiff(i-1) - l0) * tvect(i-1, :) )
            Fb(i, :) = kb * ((G(i-1) - G(i-2)) / rdiff(i-1) * nvect(i-1, :) + G(i-1) / rdiff(i) * nvect(i-1, :))
            cycle 
        end if

        Fs(i, :) = ks/(l0*l0) * ( (rdiff(i) - l0) * tvect(i, :) - (rdiff(i-1) - l0) * tvect(i-1, :) )
        Fb(i, :) = kb * ((G(i-1) - G(i-2)) / rdiff(i-1) * nvect(i-1, :) - (G(i) - G(i-1)) / rdiff(i) * nvect(i, :))
    end do

    tot_F = Fs + Fb

!    write(*,*) "Fs: ", Fs(:,0)
!    write(*, *) "Fs: ", Fs(:,1)
!
!    write(*,*) "Fb: ", Fb(:,0)
!    write(*, *) "Fb: ", Fb(:,1)
!
!    pause


    ! fibers' particle force from fibers' particle force 
    ! TODO: make it faster
    !TODO:  fpd_radius + swm_radius
    do i = 0, Np-1
      do j = i+1,Np-1
        call interaction_force(points(i,:), points(j,:), fpd_fpd_dist, k0_inter, f,lx,ly)
        tot_F(i,:) = tot_F(i,:) + f

        tot_F(j, :) = tot_F(j, :) - f
      enddo
    enddo

    ! fibers' particle force from swimmers, and swimmers from fibers' particle force
    ! TODO: make it faster
    !TODO:  fpd_radius + swm_radius
    do i = 0, Np-1
      do j = 0,Np_swm-1
        call interaction_force(points(i,:), pos_swm(j,:), fpd_swm_dist, k0_inter, f,lx,ly)
        tot_F(i,:) = tot_F(i,:) + f

        tot_F_swm(j, :) = tot_F_swm(j, :) - f
      enddo
    enddo

    ! ! swimmers from fibers' particle force: interaction force
    ! Fs = 0.0 ! temp use Fs
    ! do i = 0, Np-1
    !   do j = 0,Np_swm-1
    !     call interaction_force(points(i,:), pos_swm(j,:), fpd_radius+swm_radius, k0_inter, f)
    !     Fs(i,:) = Fs(i,:) + f

    !     tot_F_swm(j, :) = tot_F_swm(j, :) - f
    !   enddo
    ! enddo

    ! swimmers: interaction force among swimmers
    do i = 0, Np_swm-1
      do j = i+1, Np_swm-1
        call interaction_force(pos_swm(i,:), pos_swm(j,:), swm_swm_dist, k0_inter, f,lx,ly)
        tot_F_swm(i, :) = tot_F_swm(i, :) + f

        tot_F_swm(j, :) = tot_F_swm(j, :) - f
      enddo
    enddo


    ! TODO: active force for swimmers
    if (is_Tumbles_on) then
      call swm_tumbles(dt,t)
    else
      ! do i = 0, Np_swm-1
      !     active_F_swm_ang(i) = active_force_ang
      !     active_F_swm(i, 0) = active_force * cos(active_F_swm_ang(i))
      !     active_F_swm(i, 1) = active_force * sin(active_F_swm_ang(i))
      ! enddo
      if (Np_swm >= 1) then
        i = 0
        !active_F_swm_ang(i) = active_force_ang
        active_F_swm(i, 0) = active_force * cos(active_F_swm_ang(i))
        active_F_swm(i, 1) = active_force * sin(active_F_swm_ang(i))
      endif

      if (Np_swm >= 2) then
        i = 1
        !active_F_swm_ang(i) = active_force_ang + PI/2.0
        active_F_swm(i, 0) = active_force * cos(active_F_swm_ang(i))
        active_F_swm(i, 1) = active_force * sin(active_F_swm_ang(i))
      endif

    endif

    tot_F_swm = tot_F_swm + active_F_swm

    if (is_debug) then
      write(*,*) "tot_F_swm: ", tot_F_swm(:,0)
      write(*, *) "tot_F_swm: ", tot_F_swm(:,1)
      write(*, *) "tot_F_swm: end.."

      write(*,*) "tot_F: ", tot_F(:,0)
      write(*, *) "tot_F: ", tot_F(:,1)
      write(*, *) "tot_F: end.."
    endif


  return
  end subroutine

  subroutine swm_tumbles(dt,t)
  use particles_params
  use particles
  implicit none
  real dt,t
  integer i
  real tmp
    
  ! TODO: remember initialize  random_seed
  do i = 0, Np_swm-1
    if (swm_tP(i) <= 0.0) then ! when time smaller than zero, tumble


      ! ! write(*,*) "swm_tumbles", i
      ! ! for test 
      ! call random_number(tmp)
      ! swm_tP(i) = 1.0; 

      ! active_F_swm_ang(i) = 0.0;
      ! active_F_swm(i, 0) = active_force * cos(active_F_swm_ang(i))
      ! active_F_swm(i, 1) = active_force * sin(active_F_swm_ang(i))


      call random_number(tmp)
      swm_tP(i) = -log(1 - tmp) / fpd_Alpha; !Uniform2Exponential(fpd_Alpha);		

      	! double UniformPRN = genrand64_real2();
	      ! return -log(1 - UniformPRN) / lambd;

      call random_number(tmp)
      active_F_swm_ang(i) = 2.0*PI*tmp;
      active_F_swm(i, 0) = active_force * cos(active_F_swm_ang(i))
      active_F_swm(i, 1) = active_force * sin(active_F_swm_ang(i))

    else
      swm_tP(i) = swm_tP(i) - dt
    endif
  enddo  


  end subroutine


  subroutine update_position(lx,ly,dt)
    use particles_params
    use particles
    implicit none

    real :: lx,ly,dt
    integer :: j,i_num


    ! update position
    ! Assign appropriate values to r, tot_force, leed, dt, and bc
    points = points + dt * tot_F

    do i_num = 0,Np-1
      ! reset if particles move outside
      if ( points(i_num, 0) > lx) then
        points(i_num, 0) = points(i_num, 0) - lx
      else if ( points(i_num, 0) < 0) then
        points(i_num, 0) = points(i_num, 0) + lx
      endif

      if ( points(i_num, 1) > ly) then
        points(i_num, 1) = points(i_num, 1) - ly
      else if ( points(i_num, 1) < 0) then
        points(i_num, 1) = points(i_num, 1) + ly
      endif
    enddo

    !write(*,*) "i_num: ", i_num, ", u_aver: ", u_aver, ", v_aver: ", v_aver
    !write(*,*) "u_aver: ", i_num, sum_val, points(i_num, 0), points(i_num, 1)

    ! update position
    pos_swm = pos_swm + dt*tot_F_swm

    ! write(*,*) "dt*tot_F_swm: ", dt*tot_F_swm(:,0)
    ! write(*,*) "dt*tot_F_swm: ", dt*tot_F_swm(:,1)

    do i_num = 0,Np_swm-1
      ! reset if particles move outside
      if ( pos_swm(i_num, 0) > lx) then
        pos_swm(i_num, 0) = pos_swm(i_num, 0) - lx
      else if ( pos_swm(i_num, 0) < 0) then
        pos_swm(i_num, 0) = pos_swm(i_num, 0) + lx
      endif

      if ( pos_swm(i_num, 1) > ly) then
        pos_swm(i_num, 1) = pos_swm(i_num, 1) - ly
      else if ( pos_swm(i_num, 1) < 0) then
        pos_swm(i_num, 1) = pos_swm(i_num, 1) + ly
      endif
    enddo


!    write(*,*) "points: ", points(:,0)
!    write(*, *) "points: ", points(:,1)
!    pause


!!TODO: why here not used
!    if (bc == 1) then
!        do j = 0, Np-1
!            if (points(j, 1) <= 0) then
!                points(j, 1) = 0
!            end if
!        end do
!    end if
!
!    if (bc == 2) then
!        points(1, 1) = 0
!        points(Np-1, 1) = 0
!    end if
!
!    if (bc == 3) then
!        points(0, 1) = 0
!        points(1, 1) = 0
!        points(Np-1, 1) = 0
!        points(Np-2, 1) = 0
!    end if


  return
  end subroutine


  ! distribute total force on particle to body force density through phase variable
  ! particle_force1, particle_force2 are used in calculate velocity of fluid
  subroutine force_caculation_vel(nr, nz, dr, dz, r, z, lx, ly, &
      & I_thickness, particle_force1, particle_force2,dt,t)
    use particles_params
    use particles
    implicit none

    integer nr,nz
    real::r(0:nr,0:nz),z(0:nr,0:nz), dr, dz, I_thickness
    real lx, ly
    real particle_force1(-1:nr,-1:nz), particle_force2(-1:nr,-1:nz)
    real dt,t


    integer :: i,k, i_num
    real sum_val, dd
    real force1(-1:nr,-1:nz), force2(-1:nr,-1:nz)


    ! calculate total force: active_F_swm, tot_F_swm, tot_F
    call force_caculation(dt,t, lx, ly)


    ! distribute force on particles on fiber
    particle_force1 = 0.0
    particle_force2 = 0.0

    do i_num = 0,Np-1
      sum_val = 0.0
      force1 = 0.0;  force2 = 0.0;

      do k=0,nz-1
        do i=0,nr-1
          ! no need to use fpd_dist_periodic, since sphere here -> still need
          dd = fpd_dist_periodic(r(i,k),z(i,k), lx, ly, points(i_num, 0), points(i_num, 1), fpd_radius)
          !sqrt( (r(i,k)- points(i_num, 0))**2 + (z(i,k)- points(i_num, 1))**2 ) - fpd_radius
          dd = (tanh( -1.0*dd /sqrt(2.0)/I_thickness ) + 1.0)/2.0

          force1(i,k) = force1(i,k) + tot_F(i_num,0)*dd
          force2(i,k) = force2(i,k) + tot_F(i_num,1)*dd

          sum_val = sum_val + dd
        enddo
      enddo
      sum_val = sum_val*dr*dz

      particle_force1 = particle_force1 + force1 / sum_val
      particle_force2 = particle_force2 + force2 / sum_val

      !write(*,*) "tot_F: ", tot_F(i_num,0), tot_F(i_num,1), sum_val, MAXVAL(MAXVAL(ABS(particle_force1), 1)), MAXVAL(MAXVAL(ABS(particle_force2), 1))

    enddo


    ! distribute force on swimmers
    do i_num = 0,Np_swm-1
      sum_val = 0.0
      force1 = 0.0;  force2 = 0.0;

      do k=0,nz-1
        do i=0,nr-1
          dd = fpd_dist_periodic(r(i,k),z(i,k), lx, ly, pos_swm(i_num, 0), pos_swm(i_num, 1), swm_radius)
          !sqrt( (r(i,k)- pos_swm(i_num, 0))**2 + (z(i,k)- pos_swm(i_num, 1))**2 ) - swm_radius
          dd = (tanh( -1.0*dd /sqrt(2.0)/I_thickness ) + 1.0)/2.0

          force1(i,k) = force1(i,k) + tot_F_swm(i_num,0)*dd
          force2(i,k) = force2(i,k) + tot_F_swm(i_num,1)*dd

          sum_val = sum_val + dd
        enddo
      enddo
      sum_val = sum_val*dr*dz

      particle_force1 = particle_force1 + force1 / sum_val
      particle_force2 = particle_force2 + force2 / sum_val
    enddo

    ! shut down on off of Phantoms
   !TODO:  distribute force on Phantoms
   ! force along opposite directions
   if (is_Phantoms_on) then
      do i_num = 0,Np_swm-1
        sum_val = 0.0
        force1 = 0.0;  force2 = 0.0;

        do k=0,nz-1
          do i=0,nr-1
            dd = fpd_dist_periodic(r(i,k),z(i,k), lx, ly, &
              & (pos_swm(i_num, 0) - swm_bp_dist*cos(active_F_swm_ang(i_num))), &
              & (pos_swm(i_num, 1) - swm_bp_dist*sin(active_F_swm_ang(i_num))), swm_radius)
            !sqrt( (r(i,k)- (pos_swm(i_num, 0) - swm_bp_dist*cos(active_F_swm_ang(i_num))) )**2 &
            !  & + (z(i,k)- (pos_swm(i_num, 1) - swm_bp_dist*sin(active_F_swm_ang(i_num))) )**2 ) - swm_radius
            dd = (tanh( -1.0*dd /sqrt(2.0)/I_thickness ) + 1.0)/2.0

            force1(i,k) = force1(i,k) - active_F_swm(i_num,0)*dd
            force2(i,k) = force2(i,k) - active_F_swm(i_num,1)*dd

            sum_val = sum_val + dd
          enddo
        enddo
        sum_val = sum_val*dr*dz

        particle_force1 = particle_force1 + force1 / sum_val
        particle_force2 = particle_force2 + force2 / sum_val
      enddo
    endif ! end of is_Phantoms_on

    ! ghost points for force, needed in velocity calc
    particle_force1(nr,:) = particle_force1(0, :)
    particle_force1(-1,:) = particle_force1(nr-1, :)

    particle_force1(:, nz) = particle_force1(:, 0)
    particle_force1(:, -1) = particle_force1(:, nz-1)

    particle_force2(nr,:) = particle_force2(0, :)
    particle_force2(-1,:) = particle_force2(nr-1, :)

    particle_force2(:, nz) = particle_force2(:, 0)
    particle_force2(:, -1) = particle_force2(:, nz-1)


  return
  end subroutine


  ! update particles by fluid
  subroutine update_position_vel(nr,nz,r,z,lx,ly, u,v,dt,I_thickness)
    use particles_params
    use particles
    implicit none
    integer nr,nz
    real::r(0:nr,0:nz),z(0:nr,0:nz), I_thickness
    real lx, ly
    real::u(-1:nr,-1:nz),v(-1:nr,-1:nz)

    real :: dt
    integer :: i,k, i_num
    real u_aver, v_aver, sum_val, dd

    ! update fibers' particles' position
    do i_num = 0,Np-1
      u_aver = 0.
      v_aver = 0.
      sum_val = 0.0

      do k=0,nz-1
        do i=0,nr-1
          !dd = sqrt((r(i,k)- points(i_num, 0))**2 + (z(i,k)- points(i_num, 1))**2) - fpd_radius
          dd = fpd_dist_periodic(r(i,k),z(i,k), lx, ly, points(i_num, 0), points(i_num, 1), fpd_radius)

          dd = (tanh( -1.0*dd /sqrt(2.0)/I_thickness ) + 1.0)/2.0

          u_aver = u_aver + (u(i,k)+u(i,k+1))/2.0*dd
          v_aver = v_aver + (v(i,k)+v(i,k+1))/2.0*dd
          sum_val = sum_val + dd
        enddo
      enddo
      u_aver = u_aver / sum_val
      v_aver = v_aver / sum_val

      ! update position
      ! update position
      points(i_num, 0) = points(i_num, 0) + u_aver*dt
      points(i_num, 1) = points(i_num, 1) + v_aver*dt

      ! reset if particles move outside
      if ( points(i_num, 0) > lx) then
        points(i_num, 0) = points(i_num, 0) - lx
      else if ( points(i_num, 0) < 0) then
        points(i_num, 0) = points(i_num, 0) + lx
      endif

      if ( points(i_num, 1) > ly) then
        points(i_num, 1) = points(i_num, 1) - ly
      else if ( points(i_num, 1) < 0) then
        points(i_num, 1) = points(i_num, 1) + ly
      endif

      !write(*,*) "i_num: ", i_num, ", u_aver: ", u_aver, ", v_aver: ", v_aver
      !write(*,*) "u_aver: ", i_num, sum_val, points(i_num, 0), points(i_num, 1)

    enddo

    ! update swimmer's position
    do i_num = 0,Np_swm-1
      u_aver = 0.
      v_aver = 0.
      sum_val = 0.0

      do k=0,nz-1
        do i=0,nr-1
          dd = fpd_dist_periodic(r(i,k),z(i,k), lx, ly, pos_swm(i_num, 0), pos_swm(i_num, 1), swm_radius)
          !dd = sqrt((r(i,k)- pos_swm(i_num, 0))**2 + (z(i,k)- pos_swm(i_num, 1))**2) - swm_radius

          dd = (tanh( -1.0*dd /sqrt(2.0)/I_thickness ) + 1.0)/2.0

          u_aver = u_aver + (u(i,k)+u(i,k+1))/2.0*dd
          v_aver = v_aver + (v(i,k)+v(i,k+1))/2.0*dd
          sum_val = sum_val + dd
        enddo
      enddo
      u_aver = u_aver / sum_val
      v_aver = v_aver / sum_val

      ! update position
      ! update position
      pos_swm(i_num, 0) = pos_swm(i_num, 0) + u_aver*dt
      pos_swm(i_num, 1) = pos_swm(i_num, 1) + v_aver*dt

      ! reset if particles move outside
      if ( pos_swm(i_num, 0) > lx) then
        pos_swm(i_num, 0) = pos_swm(i_num, 0) - lx
      else if ( pos_swm(i_num, 0) < 0) then
        pos_swm(i_num, 0) = pos_swm(i_num, 0) + lx
      endif

      if ( pos_swm(i_num, 1) > ly) then
        pos_swm(i_num, 1) = pos_swm(i_num, 1) - ly
      else if ( pos_swm(i_num, 1) < 0) then
        pos_swm(i_num, 1) = pos_swm(i_num, 1) + ly
      endif


      ! write(*,*) "i_num: ", i_num, ", u_aver: ", u_aver, ", v_aver: ", v_aver
      ! write(*,*) "points: ", i_num, sum_val, points(i_num, 0), points(i_num, 1)
      ! write(*,*) "pos_swm: ", pos_swm(i_num, 0), pos_swm(i_num, 1)

    enddo


!    ! Assign appropriate values to r, tot_force, leed, dt, and bc
!    points = points + dt * tot_F


  return
  end subroutine

  ! ! calculate interaction dist
  ! subroutine inter_dist_periodic(lx, ly, pos_1, pos_2, di_x, di_y, d)
  ! implicit none
  ! real pos_1(0:1), pos_2(0:1)
  ! real lx, ly, di_x, di_y, d

  !   di_x = (pos_1(0) - pos_2(0))
  !   di_y = (pos_1(1) - pos_2(1))

  !   if (di_x > 0.5*lx) then 
  !     di_x = di_x - lx
  !   elseif (di_x < -0.5*lx) then
  !     di_x = di_x + lx
  !   endif

  !   if (di_y > 0.5*ly) then
  !     di_y = di_y - ly
  !   elseif (di_y < -0.5*ly) then
  !     di_y = di_y + ly
  !   endif

  !   !d = sqrt( (pos_1(0) - pos_2(0))**2 + (pos_1(1) - pos_2(1))**2 )
  !   d = sqrt( di_x**2 + di_y**2 )

  ! end subroutine  ! inter_dist_periodic


  ! calculate interaction force on particle 1 from particle 2
  subroutine interaction_force(pos_1, pos_2, eq_dist, k_inter, force,lx,ly)
  !use particles
  implicit none
  real pos_1(0:1), pos_2(0:1), eq_dist, k_inter
  real force(0:1), lx, ly

  real d, di_x, di_y

    di_x = MirrorDistance(pos_1(0) - pos_2(0), lx)
    di_y = MirrorDistance(pos_1(1) - pos_2(1), ly)
    d = sqrt( di_x**2 + di_y**2 )
    
    ! ! DOTO: consider periodic part
    ! !TODO: consider fiber
    ! di_x = (pos_1(0) - pos_2(0))
    ! di_y = (pos_1(1) - pos_2(1))

    ! if (di_x > 0.5*lx) then 
    !   di_x = di_x - lx
    ! elseif (di_x < -0.5*lx) then
    !   di_x = di_x + lx
    ! endif

    ! if (di_y > 0.5*ly) then
    !   di_y = di_y - ly
    ! elseif (di_y < -0.5*ly) then
    !   di_y = di_y + ly
    ! endif

    ! !d = sqrt( (pos_1(0) - pos_2(0))**2 + (pos_1(1) - pos_2(1))**2 )
    ! d = sqrt( di_x**2 + di_y**2 )

    if ( d < eq_dist) then
        di_x = di_x/d
        di_y = di_y/d

        force(0) = - k_inter/(eq_dist*eq_dist)* ( d - eq_dist ) * di_x
        force(1) = - k_inter/(eq_dist*eq_dist)* ( d - eq_dist ) * di_y
        !write(*,*) di_x, di_y, d, force(0), force(1)
    else
      force = 0.0
    endif

  end subroutine ! interaction_force

  function MirrorDistance(d, l)
  implicit none
  real d, l, MirrorDistance
      if (d < -0.5 * l) then
        MirrorDistance = (d + l);
      else if (d > 0.5 * l) then
        MirrorDistance = (d - l);
      else
        MirrorDistance = d
      endif
  end function MirrorDistance



  ! ! try to implement Squirmer Model
  ! subroutine calc_swm_vel(x,y, vel_x, vel_y, dist1, dist2,is_update)
  !   use particles_params
  !   use particles
  !   implicit none
  !   real x,y
  !   real vel_x, vel_y
  !   real dist1, dist2
  !   real r_x, r_y, d
  !   real cos_theta
  !   logical is_update
  
  !   !TODO: temp assume only one swimmer
  !   r_x = x - pos_swm(0,0)
  !   r_y = y - pos_swm(0,1)
  !   d = sqrt(r_x**2.0 + r_y**2.0)
  !   if (d < swm_TOL) then
  !     vel_x = 0.0; vel_y = 0.0;
  !     return
  !   endif
  
  !   if ( d <= dist1 .OR. d > dist2) then
  !     vel_x = 0.0; vel_y = 0.0;
  !     return
  !   endif
  
  !   is_update = .true.
  
  !   r_x = r_x/d; r_y = r_y/d;
  
  !   cos_theta = r_x * e_x + r_y * e_y
  
  !   vel_x = U_tilde*e_x - 1.0/3.0 * swm_B_1 * e_x + swm_B_1 * cos_theta * r_x
  !   vel_y = U_tilde*e_y - 1.0/3.0 * swm_B_1 * e_y + swm_B_1 * cos_theta * r_y
  !   return
  ! end subroutine
  
  
  ! subroutine correct_vel()
  !   use global
  !   use particles_move
  !   implicit none
  !   real vel_x, vel_y
  !   integer i,j
  !   real tmp
  !   logical is_update
  
  
  !   do j=0,nz
  !     do i=0,nr
  !        tmp = (ph_swm(i,j) + ph_swm(i-1,j))/2.0 
  !        tmp = 0.1
  !        if ( tmp <0.2 .AND. tmp > 0.05  ) then
  !         is_update = .false.
  !         call calc_swm_vel(i*dr, (j+0.5)*dz, vel_x, vel_y, 0.07, 0.08, is_update)
  !         if (is_update) then
  !           u(i,j) = vel_x
  !         endif
  !        endif
  
  !        tmp = (ph_swm(i,j) + ph_swm(i,j-1))/2.0 
  !        tmp = 0.1
  !        if ( tmp <0.2 .AND. tmp > 0.05  ) then
  !           is_update = .false.
  !           call calc_swm_vel((i+0.5)*dr, j*dz, vel_x, vel_y, 0.07, 0.08, is_update)
  !           if(is_update) then
  !             v(i,j) = vel_y
  !           endif
  
  !        endif
  
  !     enddo
  
  !   enddo

end module
