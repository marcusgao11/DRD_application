! calculate distance to interface

    function dist(x,y)
    use global
    use particles_move
    implicit none
    real x,y, dist

    real temp_h, x_pos

    integer square_index, i

    !dist = 1.0

    !dist = sqrt((x-0.5*slen)**2 + y**2) - 0.2
    !dist = sqrt((x-0.5*slen)**2 + (y-0.5*shig)**2) - 0.2

    !dist = - (x-0.5*slen)

!    temp_h = shig/3.0
!
!    if (y > (shig - temp_h)) then
!      dist = -(y - (shig - temp_h))
!    else if (y > shig/2.0) then
!      dist = (shig - temp_h) - y
!    else
!      dist = (y - temp_h)
!    endif


!! sin wave
!    temp_h = shig/5.0
!    if ( y > (temp_h + 0.1*cos(2.0*pi*x/slen*5.0) ) ) then
!      dist = 1.0 ! large value
!    else
!      dist = -1.0
!    endif


!    temp_h = shig/2.0
!    if ( y > temp_h ) then
!      dist = 1.0 ! large value
!    else
!      dist = -1.0
!    endif


!! square
!square_index = floor( (x-0.1)/0.4 + 0.5 )
!if (square_index > 4) then
!  square_index = 4
!endif
!!print *, x, square_index
!
!x_pos = 0.4*square_index + 0.1
!dist = max( abs(x-x_pos)-0.1, abs(y-0.05)-0.1)

!  ! two circles
!  dist = 1.0E-10
!  do i = 0,fpd_num
!      dist = min(dist, sqrt((x-fpd_x_pos(i))**2 + (y-fpd_y_pos(i))**2) - fpd_radius(i))
!  enddo
  dist = fpd_dist(x,y)

!  dist = 1.0E10
!  do i = 0,Np-1
!      dist = min(dist, sqrt((x-points(i, 0))**2 + (y-points(i, 1))**2) - fpd_radius)
!  enddo

  !! phi value
  !dist = tanh( -1.0*dist /sqrt(2.0)/I_thickness );
  !dist = (dist + 1.0)/2.0;


    end function dist
