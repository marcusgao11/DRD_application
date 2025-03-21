   subroutine ESOLP_x(nr,nz,pa,pb,pc,b,p,wsavex)
   implicit none
   integer nr,nz
   real pa(1:nz+1,0:nr),pb(1:nz+1,0:nr),pc(1:nz+1,0:nr)     
   real tarr(1:nz+3,1:nr),wa(1:nz+1),wb(1:nz+1),wc(1:nz+1),wk(1:nz+1)
   real b(-1:nr+1,-1:nz+1),bb(-1:nz+1,-1:nr+1),p(-1:nr+1,-1:nz+1),wsavex(3*(nr+1)+100)
   integer i,k

   integer OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS


  do k=0,nr
     do i=0,nz
        bb(i,k) = b(k,i)
     enddo
  enddo

call  VCOST(nz+1,nr+1,bb(0,0),tarr,nz+3,wsavex)

!!$OMP PARALLEL PRIVATE(k,i,wk,wa,wb,wc) DEFAULT(SHARED) 
!!$OMP DO SCHEDULE(STATIC)

do k=0,nr
!print *,OMP_GET_THREAD_NUM()
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
!!$OMP END DO
!!$OMP END PARALLEL

   call  VCOST(nz+1,nr+1,bb(0,0),tarr,nz+3,wsavex)

  do k=0,nz
     do i=0,nr
        p(i,k) = bb(k,i)
     enddo
  enddo


   end subroutine
