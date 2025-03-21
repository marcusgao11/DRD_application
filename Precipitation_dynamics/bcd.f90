!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   This bcd subroutine is for GNBC conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bcd(t)
    use global
    implicit none
    real t,vsd,cf,xlsinv
    integer i,k

    real aB(0:nr),bB(0:nr),cB(0:nr),fB(0:nr),ggB(0:nr)
    real aT(0:nr),bT(0:nr),cT(0:nr),fT(0:nr),ggT(0:nr)
    real w1(-1:nr+1,-1:nz+1),w2(-1:nr+1,-1:nz+1)

    write(6,*) 're= 7', re

!! Left BD & righte BD ghost points for velocity and ph
if (lr_bd_type == 'p') then
   do k=0,nz
       u(-1,k) = u(nr-1,k)
       u(nr+1,k) = u(1,k)
       v(-1,k) = v(nr-1,k)
       v(nr,k) = v(0,k)
   enddo

  do k=0,nz-1
   p(-1,k) = p(nr-1,k)
   p(nr,k) = p(0,k)

   ph(-1,k) = ph(nr-1,k)
   ph(nr,k) = ph(0,k)
  enddo
else  ! (lr_bd_type == 'n')
   do k=0,nz
       u(-1,k) = u(1,k)
       u(nr+1,k) = u(nr-1,k)
       v(-1,k) = v(0,k)
       v(nr,k) = v(nr-1,k)
   enddo

  do k=0,nz-1
   p(-1,k) = p(0,k)
   p(nr,k) = p(nr-1,k)

   ph(-1,k) = ph(0,k)
   ph(nr,k) = ph(nr-1,k)
  enddo
endif

 ! simple boundary conditions used at top and bottom
  do i=-1,nr+1
   u(i,-1) = -u(i,0)
   u(i,nz) = -u(i,nz-1)
  enddo

  do i=-1,nr
   v(i,0) = 0.0
   v(i,nz) = 0.0
  end do

  do i=-1,nr
    ph(i,-1) = ph(i,0)
    ph(i,nz) = ph(i,nz-1)
  enddo

! chemical potential
  do k=0,nz-1
     do i=0,nr-1
        rmu(i,k) = -I_thickness*( (ph(i+1,k)+ph(i-1,k)-2*ph(i,k))/dr2 + (ph(i,k+1)+ph(i,k-1)-2*ph(i,k))/dz2 )&
         & - ph(i,k)/I_thickness + ph(i,k)**3/I_thickness
     enddo
  enddo

if (lr_bd_type == 'p') then
  do k=0,nz-1
   rmu(-1,k) = rmu(nr-1,k)
   rmu(nr,k) = rmu(0,k)
  enddo

else
  do k=0,nz-1
   rmu(-1,k) = rmu(0,k)
   rmu(nr,k) = rmu(nr-1,k)
  enddo
endif

  do i=-1,nr
    rmu(i,-1) = rmu(i,0)
    rmu(i,nz) = rmu(i,nz-1)
  enddo



!! Calculte rho and eta
!     do j=0,nz
!        do i=-1,nr+1
!           if(ph(i,j) >= 1.0) then
!              rho(i,j)=lambda_rho
!              eta(i,j)=lambda_eta
!           elseif(ph(i,j) <= -1.0) then
!              rho(i,j)=1.0
!              eta(i,j)=1.0
!           else
!              rho(i,j) = 0.5*lambda_rho*(1.0+ph(i,j)) + 0.5*(1.0-ph(i,j))
!              eta(i,j) = 0.5*lambda_eta*(1.0+ph(i,j)) + 0.5*(1.0-ph(i,j))
!           endif
!        enddo
!    enddo
!
!    do i=0,nr
!    	j=0
!    	rho(i,j-1) = rho(i,j+2) - 3*rho(i,j+1) + 3*rho(i,j)
!    	eta(i,j-1) = eta(i,j+2) - 3*eta(i,j+1) + 3*eta(i,j)
!    	j=nz
!    	rho(i,j+1) = rho(i,j-2) - 3*rho(i,j-1) + 3*rho(i,j)
!    	eta(i,j+1) = eta(i,j-2) - 3*eta(i,j-1) + 3*eta(i,j)
!    enddo

!!!!!!!!!!!!! for chemical potential
!    call dif2x(nr,nz,xa2,xb2,xc2,ph,w1)
!    call dif2y(nr,nz,ya2,yb2,yc2,ph,w2)
!
!    do k=0,nz
!       do i=0,nr
!          rmu(i,k) = -I_thickness*w1(i,k)-I_thickness*w2(i,k)-ph(i,k)/I_thickness+ph(i,k)**3/I_thickness
!       enddo
!    enddo
!
!!! Top BD & Bottom BD ghost points for ph
!!!!!! Bottom BD
!    vsd=vs/xld/2.
!
!    do i=0,nr
!       cf=sqrt(2.0)*cos(angleB(i)*pi/180.0)/3.0*pi/2.
!       cB(i)=-2./dr2-2./dz2-vsd*dz
!       aB(i)=1.0/dr2
!       bB(i)=1.0/dr2
!       ggB(i)= -( (aB(i)*ph(i+1,0)+bB(i)*ph(i-1,0)-2.0*ph(i,0)/dr2)*I_thickness&
!& + ph(i,0)/I_thickness-ph(i,0)**3/I_thickness+(ph(i,1)-2.*ph(i,0))/dz2*I_thickness )
!       fB(i) = vsd*ph(i,1)/dz*I_thickness-2.*rmu(i,1)/dz2-vsd*dz*ggB(i) + vs/xld*cf*cos(ph(i,0)*pi/2.)
!    end do
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    aB(0)=2.0/dr2
!    bB(0)=2.0/dr2
!    aB(nr)=2.0/dr2
!    bB(nr)=2.0/dr2
!
!    call rfords(nr+1,fB,aB,cB,bB)
!
!    do i=0,nr
!       rmu(i,0) = fB(i)
!       ph(i,-1) = (ggB(i) - rmu(i,0))*dz**2/(I_thickness)
!    enddo
!
!!!! Top boundary
!
!    do i=0,nr
!       cf=sqrt(2.0)*cos(angleT(i)*pi/180.0)/3.0*pi/2.
!       cT(i)=-2./dr2-2./dz2-vsd*dz
!       aT(i)=1.0/dr2
!       bT(i)=1.0/dr2
!       ggT(i)= -( (aT(i)*ph(i+1,nz)+bT(i)*ph(i-1,nz)-2.0*ph(i,nz)/dr2)*I_thickness&
!&+ph(i,nz)/I_thickness-ph(i,nz)**3/I_thickness + (ph(i,nz-1)-2.*ph(i,nz))/dz2*I_thickness )
!       fT(i) = vsd*ph(i,nz-1)/dz*I_thickness-2.*rmu(i,nz-1)/dz2-vsd*dz*ggT(i)+vs/xld*cf*cos(ph(i,nz)*pi/2.)
!    end do
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    aT(0)=2.0/dr2
!    bT(0)=2.0/dr2
!    aT(nr)=2.0/dr2
!    bT(nr)=2.0/dr2
!
!    call rfords(nr+1,fT,aT,cT,bT)
!
!    do i=0,nr
!       rmu(i,nz) = fT(i)
!       ph(i,nz+1) = (ggT(i) - rmu(i,nz))*dz**2/I_thickness
!    enddo
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! BCDs for \partial_n \mu =0
!
!   do i=0,nr
!       rmu(i,-1) = rmu(i,1)
!       rmu(i,nz+1) = rmu(i,nz-1)
!   enddo
!
!   do k=0,nz
!      rmu(-1,k) = rmu(1,k)
!      rmu(nr+1,k) = rmu(nr-1,k)
!   enddo



! Calculte rho and eta
call interpolate_coef(ph)


! TODO:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ghost points for velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if(xlsB >= 1.0E-8) then
!      do i=0,nr
!         cf=sqrt(2.0)*cos(angleB(i)*pi/180.0)/3.0*pi/2.
!         xlsinv=0.5*(1.0+ph(i,0))/(xlsB)+0.5*(1.0-ph(i,0))/(slipratio*xlsB)
!         u(i,-1)=u(i,1) - 2.*dz*( (u(i,0)-uwB)*xlsinv - ca*( (ph(i,-1)-ph(i,1))/(2.*dz)*I_thickness&
!&-cf*cos(pi*ph(i,0)/2.) )*(xc1(i)*ph(i+1,0)+xa1(i)*ph(i-1,0)+xb1(i)*ph(i,0))/eta(i,0) )
!      end do
!      do i=0,nr
!         cf=sqrt(2.0)*cos(angleT(i)*pi/180.0)/3.0*pi/2.
!         xlsinv=0.5*(1.0+ph(i,nz))/(xlsT)+0.5*(1.0-ph(i,nz))/(slipratio*xlsT)
!         u(i,nz+1)=u(i,nz-1) - 2.*dz*( (u(i,nz)-uwT)*xlsinv - ca*((ph(i,nz+1)-ph(i,nz-1))/(2.*dz)*I_thickness&
!&-cf*cos(pi*ph(i,nz)/2.) )*(xc1(i)*ph(i+1,nz)+xa1(i)*ph(i-1,nz)+xb1(i)*ph(i,nz))/eta(i,nz) )
!      end do
!
!      do i=0,nr
!      	v(i,-1)=v(i,1)
!	      v(i,nz+1)=v(i,nz-1)
!      enddo
!endif

end subroutine
