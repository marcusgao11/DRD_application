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
