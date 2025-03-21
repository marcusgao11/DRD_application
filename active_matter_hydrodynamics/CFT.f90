! periodic BD in x direction, homogeneous \_n phi=0
! for \phi

subroutine CFT(nr,nz,fO,fN,cof,wsavex,wsavey)
  implicit none
  integer nr,nz
  real fO(-1:nr+1,-1:nz+1),fN(-1:nr+1,-1:nz+1),cof(0:nr,0:nz)

  real wsavex(nr+15),wsavey(3*(nz+1)+100)

  integer i,k
  real tarry(1:nr+3,1:nz),tarrx(1:nz+3,1:nr),temp(-1:nz+1,-1:nr+1)


!  call VRFFTF(nr,nz,fO(0,0),tarry,nr+3,wsavey)

call  VCOST(nr+1,nz+1,fO(0,0),tarrx,nr+3,wsavey)

  do k=0,nr
     do i=0,nz
        temp(i,k)=fO(k,i)
     enddo
  enddo

!  call VRFFTF(nr,nz,temp(0,0),tarry,nr+3,wsavex)

call VRFFTF(nz,nr,temp(0,0),tarrx,nz+3,wsavex)

  do k=0,nr
     do i=0,nz
        temp(i,k) = temp(i,k)/cof(k,i)
     enddo
  enddo

 !  call VRFFTB(nr,nz,temp(0,0),tarry,nr+3,wsavex)
   call VRFFTB(nz,nr,temp(0,0),tarrx,nz+3,wsavex)

   do k=0,nz
      do i=0,nr
         fN(i,k) = temp(k,i)
      enddo
   enddo

call  VCOST(nr+1,nz+1,fN(0,0),tarrx,nr+3,wsavey)

  !  call VRFFTB(nr,nz,fN(0,0),tarry,nr+3,wsavey)
  do k=0,nz
     fN(nr,k) = fN(0,k)
     fN(-1,k) = fN(nr-1,k)
     fN(nr+1,k) = fN(1,k)
  enddo
   
  do i=0,nr
     fN(i,nz+1) = fN(i,nz-1)
     fN(i,-1) = fN(i,1)
  enddo
  end subroutine
