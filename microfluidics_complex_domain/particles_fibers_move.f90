module particles_fibers_move

  contains

  subroutine particle_initialize() !N, l0)
    use particles_params
    use particles
    implicit none
    integer i

    do i = 0, Np-1
        points(i, 0) = l0 * i
        points(i, 1) = 0.0
    end do

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

  function fpd_dist(x,y)
  use particles
  implicit none
  real x,y, fpd_dist
  integer i

  fpd_dist = 1.0E10
  do i = 0,Np-1
      fpd_dist = min(fpd_dist, sqrt((x-points(i, 0))**2 + (y-points(i, 1))**2) - fpd_radius)
  enddo

  end function fpd_dist


  subroutine force_caculation()
    use particles_params
    use particles
    implicit none

    !real :: tot_F(0:Np-1, 0:1),
    real :: Fs(0:Np-1, 0:1), Fb(0:Np-1, 0:1)
    real::rdiff(0:Np-1), tvect(0:Np-1, 0:1), nvect(0:Np-1, 0:1)
    real :: theta(0:Np-1), alpha(0:Np-1), G(0:Np-1)
    integer :: i

    ! Initialize arrays
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
        rdiff(i) = sqrt(sum((points(i+1,:) - points(i,:))**2))
        tvect(i,:) = (points(i+1,:) - points(i,:))/rdiff(i)
        nvect(i,:) = [-(points(i+1,1) - points(i,1)), (points(i+1,0) - points(i,0))]/rdiff(i)
        theta(i) = atan2(tvect(i,1), tvect(i,0))
      end if
    end do

    do i = 0, Np-2
      if (i < Np-2) then
        alpha(i) = theta(i+1) - theta(i)
        if (abs(alpha(i)) >= PI/2) alpha(i) = alpha(i) + 2*PI
        G(i) = alpha(i) - Theta0
      end if
    end do

    do i = 0, Np-1
      if (i == 0) then
        Fs(i,:) = ks*(rdiff(i) - l0)*tvect(i,:)
        Fb(i,:) = -kb/l0*G(i)/rdiff(i)*nvect(i,:)
      else if (i == Np-1) then
        Fs(i,:) = -ks*(rdiff(i-1) - l0)*tvect(i-1,:)
        Fb(i,:) = -kb/l0*(G(i-2)/rdiff(i-1)*nvect(i-1,:))
      else if (i == 1) then
        Fs(i,:) = ks*(rdiff(i) - l0)*tvect(i,:) - ks*(rdiff(i-1) - l0)*tvect(i-1,:)
        Fb(i,:) = kb/l0*(G(i-1)/rdiff(i-1)*nvect(i-1,:) - (G(i) - G(i-1))/rdiff(i)*nvect(i,:))
      else if (i == Np-2) then
        Fs(i,:) = ks*(rdiff(i) - l0)*tvect(i,:) - ks*(rdiff(i-1) - l0)*tvect(i-1,:)
        Fb(i,:) = kb/l0*((G(i-1) - G(i-2))/rdiff(i-1)*nvect(i-1,:) + G(i-1)/rdiff(i)*nvect(i,:))
      else
        Fs(i,:) = ks*(rdiff(i) - l0)*tvect(i,:) - ks*(rdiff(i-1) - l0)*tvect(i-1,:)
        Fb(i,:) = kb/l0*((G(i-1) - G(i-2))/rdiff(i-1)*nvect(i-1,:) - (G(i) - G(i-1))/rdiff(i)*nvect(i,:))
      end if
    end do

    tot_F = Fs + Fb
  return
  end subroutine


  subroutine update_position(dt)
    use particles_params
    use particles
    implicit none

    real :: dt
    integer :: j


    ! Assign appropriate values to r, tot_force, leed, dt, and bc
    points = points + dt * tot_F

    if (bc == 1) then
        do j = 0, Np-1
            if (points(j, 1) <= 0) then
                points(j, 1) = 0
            end if
        end do
    end if

    if (bc == 2) then
        points(1, 1) = 0
        points(Np-1, 1) = 0
    end if

    if (bc == 3) then
        points(0, 1) = 0
        points(1, 1) = 0
        points(Np-1, 1) = 0
        points(Np-2, 1) = 0
    end if


  return
  end subroutine



end module
