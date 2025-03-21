module particles
  use particles_params
  implicit none

    real points(0:Np-1,0:1)

    real pos_swm(0:Np_swm-1, 0:1)

    real tot_F(0:Np-1, 0:1)

    real tot_F_swm(0:Np_swm-1, 0:1)
    real active_F_swm(0:Np_swm-1, 0:1)
    real active_F_swm_ang(0:Np_swm-1)

    ! tumbles
    real swm_tP(0:Np_swm-1)

    !real :: tot_F(0:Np-1, 0:1)

    ! active force
    real active_force   != 0.1 !5 !20.0 !1.0  !2.1250 !10 !0.
    real active_force_ang   != 0.0 !PI/3

end module
