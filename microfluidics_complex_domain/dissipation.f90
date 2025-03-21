  subroutine dissipation(ph_pre,dissip,mass,t)
  use global
  implicit none
  integer i,k,j
  real ph_pre(-1:nr,-1:nz),rmu_t(-1:nr+1,-1:nz+1)
  real wa(-1:nr+1,-1:nz+1),wb(-1:nr+1,-1:nz+1)
  real uua(-1:nr+1,-1:nz+1),uub(-1:nr+1,-1:nz+1)
  real vva(-1:nr+1,-1:nz+1),vvb(-1:nr+1,-1:nz+1)
  real e1,e2,e3,e4,cf,xlsinv,dissip,t,mass
  real Cx(0:nr),Cz(0:nz)

  real F,F1,F2,F3,Work !! free energy  work done to wall

Cx=1.0
Cz=1.0
Cx(0)=0.5
Cx(nr)=0.5
Cz(0)=0.5
Cz(nz)=0.5


!!!!!!!!!!!! for chemical potential
!    call dif2x(nr,nz,xa2,xb2,xc2,ph,wa)
!    call dif2y(nr,nz,ya2,yb2,yc2,ph,wb)

!    do j=0,nz
!       do i=0,nr
!          rmu_t(i,j) = -I_thickness**2*wa(i,j)-I_thickness**2*wb(i,j)+s_implicit*ph(i,j)&
!&- (1+s_implicit)*ph_pre(i,j) + ph_pre(i,j)**3
!       enddo
!    enddo
!
!   do i=0,nr
!       rmu_t(i,-1) = rmu_t(i,1)
!       rmu_t(i,nz+1) = rmu_t(i,nz-1)
!   enddo

!   do j=0,nz
!      rmu_t(-1,j) = rmu_t(1,j)
!      rmu_t(nr+1,j) = rmu_t(nr-1,j)
!   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! dissip using rmu
!    call  dif1x(nr,nz,xa1,xb1,xc1,rmu_t,wa)
!    call  dif1y(nr,nz,ya1,yb1,yc1,rmu_t,wb)

!    call  dif1x(nr,nz,xa1,xb1,xc1,rmu,wa)
!    call  dif1y(nr,nz,ya1,yb1,yc1,rmu,wb)
    e1=0.
    do k=0,nz-1
      do i=0,nr-1
         e1=e1+( ((rmu(i+1,k) - rmu(i-1,k))/2.0/dr) **2 + &
          & ((rmu(i,k+1) - rmu(i,k-1))/2.0/dz) **2 )
      enddo
    enddo
    e1=e1*xld*ca*dr*dz

    e2=0.
    do i=0,nr-1
       cf=sqrt(2.0)*cos(angleT(i)*pi/180.0)/3.0*pi/2./I_thickness  !!!! need \alpha here
       e2=e2 + ( (ph(i,nz)-ph(i,nz-1))/dz - cf*cos(pi* (ph(i,nz)+ph(i,nz-1))/2.0 /2.) )**2  !+alpha_s*(ph(i,nz)-ph_pre(i,nz)))**2

       cf=sqrt(2.0)*cos(angleB(i)*pi/180.0)/3.0*pi/2./I_thickness  !!!! need \alpha here
       e2=e2 + ( (ph(i,-1)-ph(i,0))/dz - cf*cos(pi* (ph(i,0)+ph(i,-1))/2.0 /2.) )**2  !+alpha_s*(ph(i,nz)-ph_pre(i,nz)))**2
    enddo

    e2=e2*vs*ca*I_thickness**2*dr

    call  dif1x(nr,nz,xa1,xb1,xc1,u,uua)
    call  dif1y(nr,nz,ya1,yb1,yc1,u,uub)
    call  dif1x(nr,nz,xa1,xb1,xc1,v,vva)
    !call  dif1y(nr,nz,ya1,yb1,yc1,v,vvb)


if(xlsT <= 1.0e-8) then
  do i=0,nr
     k=0
     uub(i,k)=(u(i,k+1)-u(i,k))/dz
     vvb(i,k)=(v(i,k+1)-v(i,k))/dz
     k=nz
     uub(i,k)=(u(i,k)-u(i,k-1))/dz
     vvb(i,k)=(v(i,k)-v(i,k-1))/dz
  enddo
endif


    e3=0.
    do k=0,nz
      do i=0,nr
!         e3 = e3 + Cx(i)*Cz(k)*( (uua(i,k)**2+vvb(i,k)**2)*2 + (uub(i,k)+vva(i,k))**2 )
         e3 = e3 + Cx(i)*Cz(k)*( (uua(i,k)**2+vvb(i,k)**2) + uub(i,k)**2 + vva(i,k)**2 )
      enddo
    enddo
    e3 = e3*dr*dz


    e4=0.
if(xlsB >= 1.0E-8) then
    do i=0,nr
       xlsinv=0.5*(1.0-ph(i,nz))/(slipratio*xlsT)+0.5*(1.0+ph(i,nz))/xlsT
       e4= e4 + Cx(i)*abs(u(i,nz)-uwT)**2*xlsinv
       xlsinv=0.5*(1.0-ph(i,0))/(slipratio*xlsB)+0.5*(1.0+ph(i,0))/xlsB
       e4= e4 + Cx(i)*abs(u(i,0)-uwB)**2*xlsinv
    end do
    e4=e4*dr
endif

    dissip = e1+e2+e3+e4

    write(23,110)t,e1,e2,e3,e4,dissip
!    print *,'dissipation'
!    write(6,110) t,e1,e2,e3,e4,dissip
110    format(6e14.7)


    e4=0.0
    do k=0,nz-1
       do i=0,nr-1
          e4 = e4 + ph(i,k)
       enddo
     enddo

     mass = e4*dr*dz
!     write(*,*)
!     print *,'mass',mass

!    call  dif1x(nr,nz,xa1,xb1,xc1,ph,wa)
!    call  dif1y(nr,nz,ya1,yb1,yc1,ph,wb)
!
!   F1=0.
!    do k=0,nz
!       do i=0,nr
!          F1 = F1 + Cx(i)*Cz(k)* (0.5*I_thickness**2*(wa(i,k)**2+wb(i,k)**2) + (ph(i,k)**2-1)**2/4)
!       enddo
!    enddo
!
!   F1 = F1*ca*dr*dz
!
!    F2=0.
!    do i=0,nr
!       cf= -ca*I_thickness*sqrt(2.0)*cos(angleT(i)*pi/180.0)/3.0
!       F2=F2 + Cx(i)*dr*cf*sin(pi*ph(i,nz)/2.)
!       cf= -ca*I_thickness*sqrt(2.0)*cos(angleB(i)*pi/180.0)/3.0
!       F2=F2 + Cx(i)*dr*cf*sin(pi*ph(i,0)/2.)
!    enddo
!
!    F2 = F2
!
!
!   F3=0.
!    do k=0,nz
!       do i=0,nr
!          F3 = F3 + Cx(i)*Cz(k)*(u(i,k)**2+v(i,k)**2)
!       enddo
!    enddo
!
!   F3 =F3*dr*dz*re/2.
!
!   F = F1+F2+F3
!
!   Work = 0.
!    do i=0,nr
!       Work=Work + Cx(i)*(u(i,nz)-uwT)*uwT/xlsT
!       Work=Work + Cx(i)*(u(i,0)-uwB)*uwB/xlsB
!    enddo
!   Work = Work*dr
!
!!     write(*,*)
!!   print *,'free energy',F1,F2,F3,F
!!   print *,'work done to wall by flow',Work
!!   print *,'free energy plus work per time',F+Work*dt
!     write(24,*) F1,F2,F3
!     write(25,*) F,Work*dt
  return

  end subroutine
