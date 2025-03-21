module particles_params
  ! filbers

    real,parameter::PI=4.D0*DATAN(1.D0)
!
!    integer, parameter :: bc = 1           ! boundary condition 0 free, 1 supporting, 2 hinged, 3 clamped
!    integer, parameter :: nl = 1           ! length of chain
!    real, parameter :: ks = 1000.0         ! K=k_s/k_b=50
!    real, parameter :: kb = 1.0
!    real, parameter :: R0 = 1.0            ! reference length
!    real, parameter :: Theta0 = PI / 30.0 ! np.pi/30
!
!    real,parameter::l0     = Theta0 * R0
!    real,parameter::ntheta = Theta0 / PI ! Theta0 / np.pi
!
!    integer,parameter::Np   = int(nl / ntheta) + 1
!
!    real,parameter::fpd_radius   = 0.05



    ! !=============================================================
    ! ! test cases 1
    ! real,parameter::swm_radius   = 0.01

    ! integer,parameter::Np_swm = 1

    ! real,parameter::swm_bp_dist  = 0.03 ! position determined by angle and dist


    ! real,parameter::k0_inter = 10.0 ! coefficient for interaction force

    ! ! active force
    ! real,parameter::active_force   = 10.0  !2.1250 !10 !0.
    ! real,parameter::active_force_ang   = 0.0 !PI/3


    ! ! new parameters
    ! real, parameter :: ks = 100.0         ! K=k_s/k_b=50
    ! real, parameter :: kb = 0.001
    ! real, parameter :: Theta0 = 0.0


    ! real,parameter::fpd_radius   = 0.019
    ! real,parameter::fpd_inter_dist   = 0.001


    ! real,parameter::l0     = 2*(fpd_radius+fpd_inter_dist)
    ! integer,parameter::Np   = 0 ! 9

    ! logical,parameter:: is_Phantoms_on = .true.

    ! ! parameters for tumbles
    ! logical,parameter:: is_Tumbles_on = .true.
    ! real,parameter::fpd_Alpha = 0.05


    !=============================================================
    ! test cases
    ! swimmer
    integer,parameter::Np_swm = 0 !40 !15 !25

    real,parameter::swm_radius   = 0.05/2.0

    real,parameter::swm_bp_dist  = 0.05 ! position determined by angle and dist


    ! ! active force
    ! real,parameter::active_force   = 0.1 !5 !20.0 !1.0  !2.1250 !10 !0.
    ! real,parameter::active_force_ang   = 0.0 !PI/3

    ! fiber
    ! new parameters
    real, parameter :: ks = 5 !100 !100.0  !1000 !5.0         ! K=k_s/k_b=50
    real, parameter :: kb = 0.0025*100 !0.05
    real, parameter :: Theta0 = 0.0

    integer,parameter::Np = 9

    real,parameter::fpd_radius   = 0.05/2.0
    real,parameter::fpd_inter_dist   = 0.00 !1

    real,parameter::l0 = 2*(fpd_radius+fpd_inter_dist)

    logical,parameter:: is_Phantoms_on = .true.

    ! parameters for tumbles
    logical,parameter:: is_Tumbles_on = .true. !.true.
    real,parameter::fpd_Alpha = 1/0.5 !1/0.05  ! 0.05

    ! interaction force and distance
    real,parameter::k0_inter = 2.5  !*10  !* 10 ! coefficient for interaction force
    real,parameter::k0_inter_bdr = 2.5*10
    real,parameter::fpd_swm_dist = fpd_radius + swm_radius + 0.01 !05 !* 10 ! coefficient for interaction force
    real,parameter::swm_swm_dist = swm_radius + swm_radius + 0.01 !05 !* 10 ! coefficient for interaction force
    real,parameter::fpd_fpd_dist = fpd_radius + fpd_radius + 0.01 !05 !* 10 ! coefficient for interaction force
    real,parameter::fpd_bdr_dist = fpd_radius + 0.02 ! froce applied by bdr


end module


