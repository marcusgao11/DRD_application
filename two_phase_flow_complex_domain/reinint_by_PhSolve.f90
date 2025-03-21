  subroutine reinint_by_PhSolve()
  use global
  implicit none
  integer i,j

  integer i_dt
  character lr_bd_type_bp, top_type_ph_bp, bot_type_ph_bp
  real ph_after(-1:nr,-1:nz)
  real rmu_after(-1:nr,-1:nz)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!  calcuating the RHS of CH equation
    ! initialize by solve CH equation for several steps for input structure
    ph = bdr_ph
    do j=0,nz-1
    ph(-1,j) = ph(0,j)
    ph(nr,j) = ph(nr-1,j)
    enddo

    do i=-1,nr
    ph(i,-1) = ph(i,0)
    ph(i,nz) = ph(i,nz-1)
    end do

    ! store boundary conditions type and restore it after initialization
    lr_bd_type_bp = lr_bd_type
    top_type_ph_bp = top_type_ph
    bot_type_ph_bp = bot_type_ph

    lr_bd_type = 'n'
    top_type_ph = 'n'
    bot_type_ph = 'n'
    dt = 1.0E-4

    !=================================================================
    ! need give value for initialization
    ph_M_0 = 1.0
    ph_M_inner = 1.0 ! mobility in CH equation

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  setup for the scheme
    call setup()

    do j=0,nz-1
      do i=0,nr-1
          rmu(i,j) = -I_thickness*( (ph(i+1,j)+ph(i-1,j)-2*ph(i,j))/dr2 &
            + (ph(i,j+1)+ph(i,j-1)-2*ph(i,j))/dz2 )&
          & - ph(i,j)/I_thickness + ph(i,j)**3/I_thickness
      enddo
    enddo

    write(*,*) 'initialize boundary phi by CH equation'
    do i_dt=0,50
      call PhSolve_ori(rmu_after,ph_after,0.0)
      rmu = rmu_after
      ph = ph_after
    enddo

    bdr_ph = (ph + 1.0)/2.0

    ! truncate
    do j=0,nz-1
      do i=0,nr-1
          if (bdr_ph(i,j) >= 1.0) then
            bdr_ph(i,j) = 1.0
          else if (bdr_ph(i,j) <= 0.0) then
            bdr_ph(i,j) = 0.0
          endif
      enddo
    enddo

    ! store boundary conditions type and restore it after initialization
    lr_bd_type = lr_bd_type_bp
    top_type_ph = top_type_ph_bp
    bot_type_ph = bot_type_ph_bp
    ! end of initialize boundary phi
    !--------------------------------------------------------------------------

    do j=-1,nz
      do i=-1,nr
        if (bdr_ph(i,j) > 0.5) then
          bdr_ph_penalty(i,j) = 1.0
        else
          bdr_ph_penalty(i,j) = 0.0
        endif
      enddo
    enddo
    

  end subroutine

