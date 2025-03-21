subroutine PreSolve_init(nr,nz,plan,plan_back)
implicit none
    integer*8 plan
    integer*8 plan_back
    integer nr,nz
    include 'fftw3.f'
    double precision,allocatable:: in(:,:),out(:,:)

    allocate(in(nr,nz),out(nr,nz))

    call dfftw_plan_r2r_2d(plan,nr,nz,in,out,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE)
    call dfftw_plan_r2r_2d(plan_back,nr,nz,out,in,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE)

    deallocate(in,out)

end subroutine
