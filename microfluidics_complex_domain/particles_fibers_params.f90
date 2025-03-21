module particles_params
  ! filbers

    real,parameter::PI=4.D0*DATAN(1.D0)

    integer, parameter :: bc = 1           ! boundary condition 0 free, 1 supporting, 2 hinged, 3 clamped
    integer, parameter :: nl = 1           ! length of chain
    real, parameter :: ks = 1000.0         ! K=k_s/k_b=50
    real, parameter :: kb = 1.0
    real, parameter :: R0 = 1.0            ! reference length
    real, parameter :: Theta0 = PI / 30.0 ! np.pi/30

    real,parameter::l0     = Theta0 * R0
    real,parameter::ntheta = Theta0 / PI ! Theta0 / np.pi

    integer,parameter::Np   = int(nl / ntheta) + 1
    real,parameter::radius   = 0.1


end module
