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

call  VCOST(nr,nz+1,fO(0,0),tarry,nr+3,wsavey)

  do k=0,nr-1
     do i=0,nz
        temp(i,k)=fO(k,i)
     enddo
  enddo

!  call VRFFTF(nr,nz,temp(0,0),tarry,nr+3,wsavex)

call VRFFTF(nz+1,nr,temp(0,0),tarrx,nz+3,wsavex)

  do k=0,nr-1
     do i=0,nz
        temp(i,k) = temp(i,k)/cof(k,i)
     enddo
  enddo

 !  call VRFFTB(nr,nz,temp(0,0),tarry,nr+3,wsavex)
   call VRFFTB(nz+1,nr,temp(0,0),tarrx,nz+3,wsavex)

   do k=0,nz
      do i=0,nr-1
         fN(i,k) = temp(k,i)
      enddo
   enddo

call  VCOST(nr,nz+1,fN(0,0),tarry,nr+3,wsavey)

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


  subroutine DCT(nr,nz,fO,fN,cof,wsavex,wsavey)
  implicit none
  integer nr,nz
  real fO(-1:nr+1,-1:nz+1),fN(-1:nr+1,-1:nz+1),cof(0:nr,0:nz)

  real wsavex(3*(nr+1)+100),wsavey(3*(nz+1)+100)
  integer i,k
  real tarrx(1:nr+3,1:nz),tarry(1:nz+3,1:nr),temp(-1:nz+1,-1:nr+1)!,tempp(-1:nr+1,-1:nz+1)


  call  VCOST(nr+1,nz+1,fO(0,0),tarrx,nr+3,wsavey)

  do k=0,nr
    do i=0,nz
       temp(i,k)=fO(k,i)
    enddo
  enddo

  call  VCOST(nz+1,nr+1,temp(0,0),tarry,nz+3,wsavex)


  do k=0,nz
    do i=0,nr
       fO(i,k)=temp(k,i)/cof(i,k)
    enddo
  enddo


  call  VCOST(nr+1,nz+1,fO(0,0),tarrx,nr+3,wsavey)

  do k=0,nr
    do i=0,nz
       temp(i,k)=fO(k,i)
    enddo
  enddo

  call  VCOST(nz+1,nr+1,temp(0,0),tarry,nz+3,wsavex)

  do k=0,nz
    do i=0,nr
       fN(i,k)=temp(k,i)
    enddo
  enddo

  do i=0,nr
    fN(i,-1)=fN(i,1)
    fN(i,nz+1)=fN(i,nz-1)
  enddo

  do k=0,nz
    fN(-1,k)=fN(1,k)
    fN(nr+1,k)=fN(nr-1,k)
  enddo


return
  end subroutine



  subroutine FFT(nr,nz,fO,fN,cof,wsavex,wsavey)
  implicit none
  integer nr,nz
  real fO(-1:nr+1,-1:nz+1),fN(-1:nr+1,-1:nz+1),cof(0:nr,0:nz)
  real wsavex(nr+15),wsavey(nz+15)

  integer i,k
  real tarrx(1:nr+3,1:nz),tarry(1:nz+3,1:nr),temp(-1:nz+1,-1:nr+1)


  call VRFFTF(nr,nz,fO(0,0),tarrx,nr+3,wsavey)

  do k=0,nr-1
     do i=0,nz-1
        temp(i,k)=fO(k,i)
     enddo
  enddo

  call VRFFTF(nz,nr,temp(0,0),tarry,nz+3,wsavex)

  do k=0,nr-1
     do i=0,nz-1
        temp(i,k) = temp(i,k)/cof(k,i)
     enddo
  enddo

  call VRFFTB(nz,nr,temp(0,0),tarry,nz+3,wsavex)


   do k=0,nz-1
      do i=0,nr-1
         fN(i,k) = temp(k,i)
      enddo
   enddo


  call VRFFTB(nr,nz,fN(0,0),tarrx,nr+3,wsavey)
   
  do i=0,nr
     fN(i,nz) = fN(i,0)
     fN(i,nz+1) = fN(i,1)
     fN(i,-1) = fN(i,nz-1)
  enddo

  do k=0,nz
     fN(nr,k) = fN(0,k)
     fN(nr+1,k) = fN(1,k)
     fN(-1,k) = fN(nr-1,k)
  enddo
 
  end subroutine


! nz+1 is the number of sequences need to transformed
! rhs(-1:nr+1,-1:nz+1), where rhs(0:nr-1,0:nz) contains the elments and return result is same dimension as rhs
   subroutine ESOLP_x(nr,nz,pa,pb,pc,b,p,wsavex)
   implicit none
   integer nr,nz
   real pa(1:nz+1,0:nr),pb(1:nz+1,0:nr),pc(1:nz+1,0:nr)     
   real b(-1:nr+1,-1:nz+1),bb(-1:nz+1,-1:nr+1),p(-1:nr+1,-1:nz+1),wsavex(nr+15)

   real wa(1:nz+1),wb(1:nz+1),wc(1:nz+1),wk(1:nz+1)
 
   real tarr(1:nz+3,1:nr)
   integer i,k


  do k=0,nr-1
     do i=0,nz
        bb(i,k) = b(k,i)
     enddo
  enddo

  call  VRFFTF(nz+1,nr,bb(0,0),tarr,nz+3,wsavex)

  do k=0,nr-1

     do i=1,nz+1
        wk(i)=bb(i-1,k)
     enddo

     do i=1,nz+1
        wa(i)=pa(i,k)
        wb(i)=pb(i,k)
        wc(i)=pc(i,k)
     enddo

     call rfords(nz+1,wk,wa,wb,wc)

     do i=0,nz
        bb(i,k)=wk(i+1)
     enddo

  enddo


   call  VRFFTB(nz+1,nr,bb(0,0),tarr,nz+3,wsavex)

  do k=0,nz
     do i=0,nr-1
        p(i,k) = bb(k,i)
     enddo
  enddo

  do k=0,nz
     p(nr,k)=p(0,k)
  enddo

   end subroutine
