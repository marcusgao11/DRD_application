module global

    integer,parameter::nr=384
    integer,parameter::nz=192

    !integer,parameter::nr=256
    !integer,parameter::nz=128
    real,parameter::pi=3.14159265358979323846264338327950288419716939937510
!    real,parameter::pi=2*acos(0.0)

    real dt,dr,dz,dr2,dz2
    real::r(0:nr,0:nz),z(0:nr,0:nz)
    real angleB(0:nr),angleT(0:nr)

    real xld,re,ca,vs,xlsT,xlsB,xlsL,xlsR,uw,slen,shig,angle_bot, angle_top,slipratio
    real uwT,uwB
    real I_thickness,s_implicit
    real alpha_s   !!s_implicit: parameter to inforce stablity, TODO: alpha_s is not correct in calc, and = 0.0 for now
    real ggx,ggz,tol,alphaB(0:nr),alphaT(0:nr)
    real lambda_rho,lambda_eta

    real::u(-1:nr+1,-1:nz),v(-1:nr,0:nz)
    ! here p is used to calculate in order to enforce index's correct
    real ph(-1:nr,-1:nz),rmu(-1:nr,-1:nz)
    real ph_trunc(-1:nr,-1:nz)
    real conc(-1:nr,-1:nz)  ! concentration
    real ph_grad(-1:nr,-1:nz), D_inner(-1:nr,-1:nz)  ! store |\grad\phi|
    real C_0  ! reference concentration
    real D_0  ! mobility for diffusion equation
    !real c_bot, c_top
    real c_left


    real p(-1:nr,-1:nz)
    real rho(-1:nr,-1:nz),eta(-1:nr,-1:nz)


    real fmy(0:nz),wsaveyc(3*(nz+1)+100)!,wsaveyf(nz+15)
    real fmx(0:nr),wsavexc(3*(nr+1)+100), wsavexc_v(3*nr+100) !, !wsavexf(nr+15),

    real fmy_ph(0:nz-1),wsaveyc_ph(3*nz+100)!,wsaveyf(nz+15)
    real fmx_ph(0:nr-1),wsavexc_ph(3*nr+100)


!Laplace coefficient,Velocity, pressure, \phi for fft: periodic boundary conditions in two directions
!    real Ph_Lpp(0:nr,0:nz),Ph_Lpp_2(0:nr,0:nz),Ve_Lpp(0:nr,0:nz),Pre_Lpp(0:nr,0:nz),Ve_Lpp_2(0:nr,0:nz)


!ph Laplace coefficient,Velocity Laplace coefficient,pressure Laplace for cosine
    real Ph_Lcc(0:nr-1,0:nz-1), Ph_Lcc_2(0:nr-1,0:nz-1)
    real Ve_Lcc(0:nr,0:nz),Pre_Lcc(0:nr,0:nz)
    real Ve_Lcc_2(0:nr,0:nz)

    real Ph_Lcf(0:nr,0:nz),Ph_Lcf_2(0:nr,0:nz)


    real xa1(0:nr),xb1(0:nr),xc1(0:nr),ya1(0:nz),yb1(0:nz),yc1(0:nz)
    real xa2(0:nr),xb2(0:nr),xc2(0:nr),ya2(0:nz),yb2(0:nz),yc2(0:nz)

    real ua(1:nz,0:nr),ub(1:nz,0:nr),uc(1:nz,0:nr)
    real ua_2(1:nz,0:nr),ub_2(1:nz,0:nr),uc_2(1:nz,0:nr)
    real va(1:nz-1,0:nr-1),vb(1:nz-1,0:nr-1),vc(1:nz-1,0:nr-1)
    real va_2(1:nz-1,0:nr-1),vb_2(1:nz-1,0:nr-1),vc_2(1:nz-1,0:nr-1)
    character lr_bd_type, top_type_ph, bot_type_ph
    character top_type_u, bot_type_u

    real coef_bot, coef_top ! coefficient for relaxed boundary conditions of phi

    character only_cahn_hillarid

    real,parameter::penalty_k = 1.0E5


!!!!!!!         for CSR form store of the matrix for equation  of ph and mu
!    integer,parameter::NN=2*(nr+1)*(nz+1)
!    real aqin(1:7*NN)
!    integer iaqin(1:NN+1),jaqin(1:7*NN)

!    public Ph_Lcf,Ph_Lcf_2

    real shift_delta

    public shift_delta

    real uL(0:nz), uR(0:nz), UU ! velocity at inlet, and shared in pressure solver




    public nr,nz,dt,dr,dz,dr2,dz2,r,z
    public u,v,p,ph,rmu,conc, ph_grad, D_inner, ph_trunc
    public C_0, D_0, c_bot, c_top, c_left

    public xld,re,ca,vs,xlsT,xlsB,xlsL,xlsR,uw,slen,shig,angle,slipratio,I_thickness
    public uwT,uwB,tol,alpha_s,s_implicit,ggz,ggx
    public angleB,angleT,alphaB,alphaT,lambda_rho,lambda_eta

    public fmy,wsaveyc,fmx,wsavexc,wsavexc_v,Ph_Lcc,Ve_Lcc,Pre_Lcc,Ph_Lcc_2,Ve_Lcc_2!,Pre_Lf,wsavexf,wsaveyf
    public Ph_Lcf,Ph_Lcf_2
    public wsavexf!,wsaveyf
!    public Ph_Lpp,Ph_Lpp_2,Ve_Lpp,Ve_Lpp_2,Pre_Lpp


    public xa1,xb1,xc1,ya1,yb1,yc1
    public xa2,xb2,xc2,ya2,yb2,yc2

    public ua,ub,uc,va,vb,vc
    public ua_2,ub_2,uc_2,va_2,vb_2,vc_2

    public coef_bot, coef_top, only_cahn_hillarid

    public uL, uR, UU! velocity at inlet, and shared in pressure solver


end module
