! calculate distance to interface

    function dist(x,y)
    use global
    implicit none
    real x,y, dist

    real temp_h, x_pos

    integer square_index, i

    real radius, width

     !dist = sqrt((x-0.8*slen)**2 + y**2) - 0.5
     !dist = sqrt((x-1.45)**2 + y**2) - 0.5

    ! positive outside bubble

    !dist = -(y - 0.1)

    ! ! semi-circle
    ! dist = sqrt((x-0.5*slen)**2 + y**2) - 0.2


    dist = 1.0E10

    dist = min(dist, sqrt((x-0.4)**2 + (y-0.3)**2) - 0.1)

    dist = min(dist, sqrt((x-0.4)**2 + (y-0.7)**2) - 0.1)

    dist = min(dist, sqrt((x-0.7)**2 + (y-0.5)**2) - 0.1)

!    if ((x > 0.0).and.(x < 0.55).and.(y < 0.5)) then
!        dist = ( sqrt((x-0.4)**2 + (y-0.3)**2) - 0.1 )
!    endif


!    if ((x > 0.0).and.(x < 0.55).and.(y > 0.5)) then
!        dist = ( sqrt((x-0.4)**2 + (y-0.7)**2) - 0.1 )
!    endif

!    if ((x > 0.55).and.(x < 0.85)) then
!        dist = ( sqrt((x-0.7)**2 + (y-0.5)**2) - 0.1 )
!    endif


    dist = min(dist, sqrt((x-1.0)**2 + (y-0.3)**2) - 0.1 )
    dist = min(dist, sqrt((x-1.0)**2 + (y-0.7)**2) - 0.1 )

    dist = min(dist, sqrt((x-1.3)**2 + (y-0.5)**2) - 0.1 )

!    if ((x > 0.85).and.(x < 1.15).and.(y < 0.5)) then
!        dist = ( sqrt((x-1.0)**2 + (y-0.3)**2) - 0.1 )
!    endif


!    if ((x > 0.85).and.(x < 1.15).and.(y > 0.5)) then
!        dist = ( sqrt((x-1.0)**2 + (y-0.7)**2) - 0.1 )
!    endif

!    if ((x > 1.15).and.(x < 1.45)) then
!        dist = ( sqrt((x-1.3)**2 + (y-0.5)**2) - 0.1 )
!    endif

    dist = min(dist, sqrt((x-1.6)**2 + (y-0.3)**2) - 0.1 )
    dist = min(dist, sqrt((x-1.6)**2 + (y-0.7)**2) - 0.1 )

!    if ((x > 1.45).and.(y < 0.5)) then
!        dist = ( sqrt((x-1.6)**2 + (y-0.3)**2) - 0.1 )
!    endif
!
!    if ((x > 1.45).and.(y > 0.5)) then
!        dist = ( sqrt((x-1.6)**2 + (y-0.7)**2) - 0.1 )
!    endif


!    dist = 1.0E10
!    radius = 0.15; width = 0.2;
!    do i = 0,1
!      dist = min(dist, sqrt((x-0.4)**2 + (y- (width/2.0 + radius + i*(width+2.0*radius)) )**2) - radius )
!    enddo
!
!    do i = 0,1
!      dist = min(dist, sqrt((x-1.0)**2 + (y- (width/2.0 + radius + i*(width+2.0*radius)) )**2) - radius )
!    enddo
!
!    do i = 0,1
!      dist = min(dist, sqrt((x-1.6)**2 + (y- (width/2.0 + radius + i*(width+2.0*radius)) )**2) - radius )
!    enddo

!    dist = 1.0E10
!    radius = 0.075; width = 0.1;
!    do i = 0,3
!      dist = min(dist, sqrt((x-0.4)**2 + (y- (width/2.0 + radius + i*(width+2.0*radius)) )**2) - radius )
!    enddo
!    do i = 0,3
!      dist = min(dist, sqrt((x-1.0)**2 + (y- (width/2.0 + radius + i*(width+2.0*radius)) )**2) - radius )
!    enddo
!    do i = 0,3
!      dist = min(dist, sqrt((x-1.6)**2 + (y- (width/2.0 + radius + i*(width+2.0*radius)) )**2) - radius )
!    enddo


    !dist = ( sqrt((x-0.2*slen)**2 + (y-0.5*shig)**2) - 0.1 )
    !dist = sqrt((x-0.5*slen)**2/2.0 + (y-0.5*shig)**2) - 0.2

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


    end function dist
