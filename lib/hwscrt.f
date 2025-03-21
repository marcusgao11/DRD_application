
      Subroutine hwscrt(a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     1                   elmbda,f,idimf,pertrb,ierror,w)
c
c
c       DOUBLE PRECISION VERSION  - LOCAL MODIFICATION
c
	implicit real*8(a-h,o-z)
c
c   
c
c
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                 *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c     * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c          subroutine hwscrt solves the standard five-point finite
c     difference approximation to the helmholtz equation in cartesian
c     coordinates:
c
c          (d/dx)(du/dx) + (d/dy)(du/dy) + lambda*u = f(x,y).
c
c
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     a,b
c       the range of x, i.e., a .le. x .le. b.  a must be less than b.
c
c     m
c       the number of panels into which the interval (a,b) is
c       subdivided.  hence, there will be m+1 grid points in the
c       x-direction given by x(i) = a+(i-1)dx for i = 1,2,...,m+1,
c       where dx = (b-a)/m is the panel width. m must be greater than 3.
c
c     mbdcnd
c       indicates the type of boundary conditions at x = a and x = b.
c
c       = 0  if the solution is periodic in x, i.e., u(i,j) = u(m+i,j).
c       = 1  if the solution is specified at x = a and x = b.
c       = 2  if the solution is specified at x = a and the derivative of
c            the solution with respect to x is specified at x = b.
c       = 3  if the derivative of the solution with respect to x is
c            specified at x = a and x = b.
c       = 4  if the derivative of the solution with respect to x is
c            specified at x = a and the solution is specified at x = b.
c
c     bda
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to x at x = a.
c       when mbdcnd = 3 or 4,
c
c            bda(j) = (d/dx)u(a,y(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bda is a dummy variable.
c
c     bdb
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to x at x = b.
c       when mbdcnd = 2 or 3,
c
c            bdb(j) = (d/dx)u(b,y(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value bdb is a dummy variable.
c
c     c,d
c       the range of y, i.e., c .le. y .le. d.  c must be less than d.
c
c     n
c       the number of panels into which the interval (c,d) is
c       subdivided.  hence, there will be n+1 grid points in the
c       y-direction given by y(j) = c+(j-1)dy for j = 1,2,...,n+1, where
c       dy = (d-c)/n is the panel width.  n must be greater than 3.
c
c     nbdcnd
c       indicates the type of boundary conditions at y = c and y = d.
c
c       = 0  if the solution is periodic in y, i.e., u(i,j) = u(i,n+j).
c       = 1  if the solution is specified at y = c and y = d.
c       = 2  if the solution is specified at y = c and the derivative of
c            the solution with respect to y is specified at y = d.
c       = 3  if the derivative of the solution with respect to y is
c            specified at y = c and y = d.
c       = 4  if the derivative of the solution with respect to y is
c            specified at y = c and the solution is specified at y = d.
c
c     bdc
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to y at y = c.
c       when nbdcnd = 3 or 4,
c
c            bdc(i) = (d/dy)u(x(i),c), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdc is a dummy variable.
c
c     bdd
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to y at y = d.
c       when nbdcnd = 2 or 3,
c
c            bdd(i) = (d/dy)u(x(i),d), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdd is a dummy variable.
c
c     elmbda
c       the constant lambda in the helmholtz equation.  if
c       lambda .gt. 0, a solution may not exist.  however, hwscrt will
c       attempt to find a solution.
c
c     f
c       a two-dimensional array which specifies the values of the right
c       side of the helmholtz equation and boundary values (if any).
c       for i = 2,3,...,m and j = 2,3,...,n
c
c            f(i,j) = f(x(i),y(j)).
c
c       on the boundaries f is defined by
c
c            mbdcnd     f(1,j)        f(m+1,j)
c            ------     ---------     --------
c
c              0        f(a,y(j))     f(a,y(j))
c              1        u(a,y(j))     u(b,y(j))
c              2        u(a,y(j))     f(b,y(j))     j = 1,2,...,n+1
c              3        f(a,y(j))     f(b,y(j))
c              4        f(a,y(j))     u(b,y(j))
c
c
c            nbdcnd     f(i,1)        f(i,n+1)
c            ------     ---------     --------
c
c              0        f(x(i),c)     f(x(i),c)
c              1        u(x(i),c)     u(x(i),d)
c              2        u(x(i),c)     f(x(i),d)     i = 1,2,...,m+1
c              3        f(x(i),c)     f(x(i),d)
c              4        f(x(i),c)     u(x(i),d)
c
c       f must be dimensioned at least (m+1)*(n+1).
c
c       note
c
c       if the table calls for both the solution u and the right side f
c       at  a corner then the solution must be specified.
c
c     idimf
c       the row (or first) dimension of the array f as it appears in the
c       program calling hwscrt.  this parameter is used to specify the
c       variable dimension of f.  idimf must be at least m+1  .
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.  w may require up to 4*(n+1) +
c       (13 + int(log2(n+1)))*(m+1) locations.  the actual number of
c       locations used is computed by hwscrt and is returned in location
c       w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c     f
c       contains the solution u(i,j) of the finite difference
c       approximation for the grid point (x(i),y(j)), i = 1,2,...,m+1,
c       j = 1,2,...,n+1  .
c
c     pertrb
c       if a combination of periodic or derivative boundary conditions
c       is specified for a poisson equation (lambda = 0), a solution may
c       not exist.  pertrb is a constant, calculated and subtracted from
c       f, which ensures that a solution exists.  hwscrt then computes
c       this solution, which is a least squares solution to the original
c       approximation.  this solution plus any constant is also a
c       solution.  hence, the solution is not unique.  the value of
c       pertrb should be small compared to the right side f.  otherwise,
c       a solution is obtained to an essentially different problem.
c       this comparison should always be made to insure that a
c       meaningful solution has been obtained.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for numbers 0 and 6, a solution is not attempted.
c
c       = 0  no error.
c       = 1  a .ge. b.
c       = 2  mbdcnd .lt. 0 or mbdcnd .gt. 4  .
c       = 3  c .ge. d.
c       = 4  n .le. 3
c       = 5  nbdcnd .lt. 0 or nbdcnd .gt. 4  .
c       = 6  lambda .gt. 0  .
c       = 7  idimf .lt. m+1  .
c       = 8  m .le. 3
c
c       since this is the only means of indicating a possibly incorrect
c       call to hwscrt, the user should test ierror after the call.
c
c     w
c       w(1) contains the required length of w.
c
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c
c     dimension of   bda(n+1),bdb(n+1),bdc(m+1),bdd(m+1),f(idimf,n+1),
c     arguments      w(see argument list)
c
c     latest         june 1, 1976
c     revision
c
c     subprograms    hwscrt,genbun,poisd2,poisn2,poisp2,cosgen,merge,
c     required       trix,tri3,pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        standardized september 1, 1973
c                    revised april 1, 1976
c
c     algorithm      the routine defines the finite difference
c                    equations, incorporates boundary data, and adjusts
c                    the right side of singular systems and then calls
c                    genbun to solve the system.
c
c     space          13110(octal) = 5704(decimal) locations on the ncar
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hwscrt is roughly proportional
c                    to m*n*log2(n), but also depends on the input
c                    parameters nbdcnd and mbdcnd.  some typical values
c                    are listed in the table below.
c                       the solution process employed results in a loss
c                    of no more than three significant digits for n and
c                    m as large as 64.  more detailed information about
c                    accuracy can be found in the documentation for
c                    subroutine genbun which is the routine that
c                    solves the finite difference equations.
c
c
c                       m(=n)    mbdcnd    nbdcnd    t(msecs)
c                       -----    ------    ------    --------
c
c                        32        0         0          31
c                        32        1         1          23
c                        32        3         3          36
c                        64        0         0         128
c                        64        1         1          96
c                        64        3         3         142
c
c     portability    american national standards institute fortran.
c                    all machine dependent constants are located in the
c                    function pimach.
c
c     reference      swarztrauber,p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
      dimension       f(idimf,1)
      dimension       bda(1)     ,bdb(1)     ,bdc(1)     ,bdd(1)     ,
     1                w(1)
c
c     check for invalid parameters.
c

      ierror = 0
      if (a .ge. b) ierror = 1
      if (mbdcnd.lt.0 .or. mbdcnd.gt.4) ierror = 2
      if (c .ge. d) ierror = 3
      if (n .le. 3) ierror = 4
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) ierror = 5
      if (idimf .lt. m+1) ierror = 7
      if (m .le. 3) ierror = 8
      if (ierror .ne. 0) return
      nperod = nbdcnd
      mperod = 0
      if (mbdcnd .gt. 0) mperod = 1
      deltax = (b-a)/dble(m)
      twdelx = 2.d0/deltax
      delxsq = 1.d0/deltax**2
      deltay = (d-c)/dble(n)
      twdely = 2.d0/deltay
      delysq = 1.d0/deltay**2
      np = nbdcnd+1
      np1 = n+1
      mp = mbdcnd+1
      mp1 = m+1
      nstart = 1
      nstop = n
      nskip = 1
      go to (104,101,102,103,104),np
  101 nstart = 2
      go to 104
  102 nstart = 2
  103 nstop = np1
      nskip = 2
  104 nunk = nstop-nstart+1
c
c     enter boundary data for x-boundaries.
c
      mstart = 1
      mstop = m
      mskip = 1
      go to (117,105,106,109,110),mp
  105 mstart = 2
      go to 107
  106 mstart = 2
      mstop = mp1
      mskip = 2
  107 do 108 j=nstart,nstop
         f(2,j) = f(2,j)-f(1,j)*delxsq
  108 continue
      go to 112
  109 mstop = mp1
      mskip = 2
  110 do 111 j=nstart,nstop
         f(1,j) = f(1,j)+bda(j)*twdelx
  111 continue
  112 go to (113,115),mskip
  113 do 114 j=nstart,nstop
         f(m,j) = f(m,j)-f(mp1,j)*delxsq
  114 continue
      go to 117
  115 do 116 j=nstart,nstop
         f(mp1,j) = f(mp1,j)-bdb(j)*twdelx
  116 continue
  117 munk = mstop-mstart+1
c
c     enter boundary data for y-boundaries.
c
      go to (127,118,118,120,120),np
  118 do 119 i=mstart,mstop
         f(i,2) = f(i,2)-f(i,1)*delysq
  119 continue
      go to 122
  120 do 121 i=mstart,mstop
         f(i,1) = f(i,1)+bdc(i)*twdely
  121 continue
  122 go to (123,125),nskip
  123 do 124 i=mstart,mstop
         f(i,n) = f(i,n)-f(i,np1)*delysq
  124 continue
      go to 127
  125 do 126 i=mstart,mstop
         f(i,np1) = f(i,np1)-bdd(i)*twdely
  126 continue
c
c    multiply right side by deltay**2.
c
  127 delysq = deltay*deltay
      do 129 i=mstart,mstop
         do 128 j=nstart,nstop
            f(i,j) = f(i,j)*delysq
  128    continue
  129 continue
c
c     define the a,b,c coefficients in w-array.
c
      id2 = munk
      id3 = id2+munk
      id4 = id3+munk
      s = delysq*delxsq
      st2 = 2.d0*s
      do 130 i=1,munk
         w(i) = s
         j = id2+i
         w(j) = -st2+elmbda*delysq
         j = id3+i
         w(j) = s
  130 continue
      if (mp .eq. 1) go to 131
      w(1) = 0.d0
      w(id4) = 0.d0
  131 continue
      go to (135,135,132,133,134),mp
  132 w(id2) = st2
      go to 135
  133 w(id2) = st2
  134 w(id3+1) = st2
  135 continue
      pertrb = 0.d0
      if (elmbda) 144,137,136
  136 ierror = 6
      go to 144
  137 if ((nbdcnd.eq.0 .or. nbdcnd.eq.3) .and.
     1    (mbdcnd.eq.0 .or. mbdcnd.eq.3)) go to 138
      go to 144
c
c     for singular problems must adjust data to insure that a solution
c     will exist.
c
  138 a1 = 1.
      a2 = 1.
      if (nbdcnd .eq. 3) a2 = 2.
      if (mbdcnd .eq. 3) a1 = 2.
      s1 = 0.d0
      msp1 = mstart+1
      mstm1 = mstop-1
      nsp1 = nstart+1
      nstm1 = nstop-1
      do 140 j=nsp1,nstm1
         s = 0.d0
         do 139 i=msp1,mstm1
            s = s+f(i,j)
  139    continue
         s1 = s1+s*a1+f(mstart,j)+f(mstop,j)
  140 continue
      s1 = a2*s1
      s = 0.d0
      do 141 i=msp1,mstm1
         s = s+f(i,nstart)+f(i,nstop)
  141 continue
      s1 = s1+s*a1+f(mstart,nstart)+f(mstart,nstop)+f(mstop,nstart)+
     1     f(mstop,nstop)
      s = (2.d0+dble(nunk-2)*a2)*(2.d0+dble(munk-2)*a1)
      pertrb = s1/s
      do 143 j=nstart,nstop
         do 142 i=mstart,mstop
            f(i,j) = f(i,j)-pertrb
  142    continue
  143 continue
      pertrb = pertrb/delysq
c
c     solve the equation.
c
  144 call genbun (nperod,nunk,mperod,munk,w(1),w(id2+1),w(id3+1),
     1             idimf,f(mstart,nstart),ierr1,w(id4+1))
      w(1) = w(id4+1)+3.d0*dble(munk)
c
c     fill in identical values when have periodic boundary conditions.
c
      if (nbdcnd .ne. 0) go to 146
      do 145 i=mstart,mstop
         f(i,np1) = f(i,1)
  145 continue
  146 if (mbdcnd .ne. 0) go to 148
      do 147 j=nstart,nstop
         f(mp1,j) = f(1,j)
  147 continue
      if (nbdcnd .eq. 0) f(mp1,np1) = f(1,np1)
  148 continue
      return
      end
      subroutine genbun (nperod,n,mperod,m,a,b,c,idimy,y,ierror,w)
c
c
      implicit real*8(a-h,o-z)
c
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                  *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c     * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c
c     subroutine genbun solves the linear system of equations
c
c          a(i)*x(i-1,j) + b(i)*x(i,j) + c(i)*x(i+1,j)
c
c          + x(i,j-1) - 2.*x(i,j) + x(i,j+1) = y(i,j)
c
c               for i = 1,2,...,m  and  j = 1,2,...,n.
c
c     the indices i+1 and i-1 are evaluated modulo m, i.e.,
c     x(0,j) = x(m,j) and x(m+1,j) = x(1,j), and x(i,0) may be equal to
c     0, x(i,2), or x(i,n) and x(i,n+1) may be equal to 0, x(i,n-1), or
c     x(i,1) depending on an input parameter.
c
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     nperod
c       indicates the values that x(i,0) and x(i,n+1) are assumed to
c       have.
c
c       = 0  if x(i,0) = x(i,n) and x(i,n+1) = x(i,1).
c       = 1  if x(i,0) = x(i,n+1) = 0  .
c       = 2  if x(i,0) = 0 and x(i,n+1) = x(i,n-1).
c       = 3  if x(i,0) = x(i,2) and x(i,n+1) = x(i,n-1).
c       = 4  if x(i,0) = x(i,2) and x(i,n+1) = 0.
c
c     n
c       the number of unknowns in the j-direction.  n must be greater
c       than 2.
c
c     mperod
c       = 0 if a(1) and c(m) are not zero
c       = 1 if a(1) = c(m) = 0
c
c     m
c       the number of unknowns in the i-direction.  m must be greater
c       than 2.
c
c     a,b,c
c       one-dimensional arrays of length m that specify the
c       coefficients in the linear equations given above.  if mperod = 0
c       the array elements must not depend upon the index i, but must be
c       constant.  specifically, the subroutine checks the following
c       condition
c
c             a(i) = c(1)
c             c(i) = c(1)
c             b(i) = b(1)
c
c       for i=1,2,...,m.
c
c     idimy
c       the row (or first) dimension of the two-dimensional array y as
c       it appears in the program calling genbun.  this parameter is
c       used to specify the variable dimension of y.  idimy must be at
c       least m.
c
c     y
c       a two-dimensional array that specifies the values of the right
c       side of the linear system of equations given above.  y must be
c       dimensioned at least m*n.
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.  w may require up to 4*n + (10 + int(log2(n)))*m
c       locations.  the actual number of locations used is computed by
c       genbun and is returned in location w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c     y
c       contains the solution x.
c
c     ierror
c       an error flag that indicates invalid input parameters  except
c       for number zero, a solution is not attempted.
c
c       = 0  no error.
c       = 1  m .le. 2  .
c       = 2  n .le. 2
c       = 3  idimy .lt. m
c       = 4  nperod .lt. 0 or nperod .gt. 4
c       = 5  mperod .lt. 0 or mperod .gt. 1
c       = 6  a(i) .ne. c(1) or c(i) .ne. c(1) or b(i) .ne. b(1) for
c            some i=1,2,...,m.
c       = 7  a(1) .ne. 0 or c(m) .ne. 0 and mperod = 1
c
c     w
c       w(1) contains the required length of w.
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   a(m),b(m),c(m),y(idimy,n),w(see parameter list)
c     arguments
c
c     latest         june 1, 1976
c     revision
c
c     subprograms    genbun,poisd2,poisn2,poisp2,cosgen,merge,trix,tri3,
c     required       pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        standardized april 1, 1973
c                    revised august 20,1973
c                    revised january 1, 1976
c
c     algorithm      the linear system is solved by a cyclic reduction
c                    algorithm described in the reference.
c
c     space          4944(decimal) = 11520(octal) locations on the ncar
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine genbun is roughly proportional
c                    to m*n*log2(n), but also depends on the input
c                    parameter nperod.  some typical values are listed
c                    in the table below.  more comprehensive timing
c                    charts may be found in the reference.
c                       to measure the accuracy of the algorithm a
c                    uniform random number generator was used to create
c                    a solution array x for the system given in the
c                    'purpose' with
c
c                       a(i) = c(i) = -0.5*b(i) = 1,       i=1,2,...,m
c
c                    and, when mperod = 1
c
c                       a(1) = c(m) = 0
c                       a(m) = c(1) = 2.
c
c                    the solution x was substituted into the given sys-
c                    tem and, using double precision, a right side y was
c                    computed.  using this array y subroutine genbun was
c                    called to produce an approximate solution z.  then
c                    the relative error, defined as
c
c                       e = max(abs(z(i,j)-x(i,j)))/max(abs(x(i,j)))
c
c                    where the two maxima are taken over all i=1,2,...,m
c                    and j=1,2,...,n, was computed.  the value of e is
c                    given in the table below for some typical values of
c                    m and n.
c
c
c                       m (=n)    mperod    nperod    t(msecs)    e
c                       ------    ------    ------    --------  ------
c
c                         31        0         0          36     6.e-14
c                         31        1         1          21     4.e-13
c                         31        1         3          41     3.e-13
c                         32        0         0          29     9.e-14
c                         32        1         1          32     3.e-13
c                         32        1         3          48     1.e-13
c                         33        0         0          36     9.e-14
c                         33        1         1          30     4.e-13
c                         33        1         3          34     1.e-13
c                         63        0         0         150     1.e-13
c                         63        1         1          91     1.e-12
c                         63        1         3         173     2.e-13
c                         64        0         0         122     1.e-13
c                         64        1         1         128     1.e-12
c                         64        1         3         199     6.e-13
c                         65        0         0         143     2.e-13
c                         65        1         1         120     1.e-12
c                         65        1         3         138     4.e-13
c
c     portability    american national standards institue fortran.
c                    all machine dependent constants are located in the
c                    function pimach.
c
c     required       cos
c     resident
c     routines
c
c     reference      sweet, r., 'a cyclic reduction algorithm for
c                    solving block tridiagonal systems of arbitrary
c                    dimensions,' siam j. on numer. anal.,
c                    14(sept., 1977), pp. 706-720.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
      dimension       y(idimy,1)
      dimension       w(1)       ,b(1)       ,a(1)       ,c(1)
      ierror = 0
      if (m .le. 2) ierror = 1
      if (n .le. 2) ierror = 2
      if (idimy .lt. m) ierror = 3
      if (nperod.lt.0 .or. nperod.gt.4) ierror = 4
      if (mperod.lt.0 .or. mperod.gt.1) ierror = 5
      if (mperod .eq. 1) go to 102
      do 101 i=2,m
         if (a(i) .ne. c(1)) go to 103
         if (c(i) .ne. c(1)) go to 103
         if (b(i) .ne. b(1)) go to 103
  101 continue
      go to 104
  102 if (a(1).ne.0. .or. c(m).ne.0.) ierror = 7
      go to 104
  103 ierror = 6
  104 if (ierror .ne. 0) return
      mp1 = m+1
      iwba = mp1
      iwbb = iwba+m
      iwbc = iwbb+m
      iwb2 = iwbc+m
      iwb3 = iwb2+m
      iww1 = iwb3+m
      iww2 = iww1+m
      iww3 = iww2+m
      iwd = iww3+m
      iwtcos = iwd+m
      iwp = iwtcos+4*n
      do 106 i=1,m
         k = iwba+i-1
         w(k) = -a(i)
         k = iwbc+i-1
         w(k) = -c(i)
         k = iwbb+i-1
         w(k) = 2.-b(i)
         do 105 j=1,n
            y(i,j) = -y(i,j)
  105    continue
  106 continue
      mp = mperod+1
      np = nperod+1
      go to (114,107),mp
  107 go to (108,109,110,111,123),np
  108 call poisp2 (m,n,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2),
     1             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),
     2             w(iwp))
      go to 112
  109 call poisd2 (m,n,1,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iww1),
     1             w(iwd),w(iwtcos),w(iwp))
      go to 112
  110 call poisn2 (m,n,1,2,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2),
     1             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),
     2             w(iwp))
      go to 112
  111 call poisn2 (m,n,1,1,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2),
     1             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),
     2             w(iwp))
  112 ipstor = w(iww1)
      irev = 2
      if (nperod .eq. 4) go to 124
  113 go to (127,133),mp
  114 continue
c
c     reorder unknowns when mp =0
c
      mh = (m+1)/2
      mhm1 = mh-1
      modd = 1
      if (mh*2 .eq. m) modd = 2
      do 119 j=1,n
         do 115 i=1,mhm1
            mhpi = mh+i
            mhmi = mh-i
            w(i) = y(mhmi,j)-y(mhpi,j)
            w(mhpi) = y(mhmi,j)+y(mhpi,j)
  115    continue
         w(mh) = 2.d0*y(mh,j)
         go to (117,116),modd
  116    w(m) = 2.d0*y(m,j)
  117    continue
         do 118 i=1,m
            y(i,j) = w(i)
  118    continue
  119 continue
      k = iwbc+mhm1-1
      i = iwba+mhm1
      w(k) = 0.d0
      w(i) = 0.d0
      w(k+1) = 2.d0*w(k+1)
      go to (120,121),modd
  120 continue
      k = iwbb+mhm1-1
      w(k) = w(k)-w(i-1)
      w(iwbc-1) = w(iwbc-1)+w(iwbb-1)
      go to 122
  121 w(iwbb-1) = w(k+1)
  122 continue
      go to 107
c
c     reverse columns when nperod = 4.
c
  123 irev = 1
      nby2 = n/2
  124 do 126 j=1,nby2
         mskip = n+1-j
         do 125 i=1,m
            a1 = y(i,j)
            y(i,j) = y(i,mskip)
            y(i,mskip) = a1
  125    continue
  126 continue
      go to (110,113),irev
  127 continue
      do 132 j=1,n
         do 128 i=1,mhm1
            mhmi = mh-i
            mhpi = mh+i
            w(mhmi) = .5d0*(y(mhpi,j)+y(i,j))
            w(mhpi) = .5d0*(y(mhpi,j)-y(i,j))
  128    continue
         w(mh) = .5d0*y(mh,j)
         go to (130,129),modd
  129    w(m) = .5d0*y(m,j)
  130    continue
         do 131 i=1,m
            y(i,j) = w(i)
  131    continue
  132 continue
  133 continue
c
c     return storage requirements for w array.
c
      w(1) = ipstor+iwp-1
      return
      end
      subroutine poisd2 (mr,nr,istag,ba,bb,bc,q,idimq,b,w,d,tcos,p)
      implicit real*8(a-h,o-z)
c
c     subroutine to solve poisson's equation for dirichlet boundary
c     conditions.
c
c     istag = 1 if the last diagonal block is the matrix a.
c     istag = 2 if the last diagonal block is the matrix a+i.
c
      dimension       q(idimq,1) ,ba(1)      ,bb(1)      ,bc(1)      ,
     1                tcos(1)    ,b(1)       ,d(1)       ,w(1)       ,
     2                p(1)
      m = mr
      n = nr
      jsh = 0
      fi = 1.d0/dble(istag)
      ip = -m
      ipstor = 0
      go to (101,102),istag
  101 kr = 0
      irreg = 1
      if (n .gt. 1) go to 106
      tcos(1) = 0.d0
      go to 103
  102 kr = 1
      jstsav = 1
      irreg = 2
      if (n .gt. 1) go to 106
      tcos(1) = -1.d0
  103 do 104 i=1,m
         b(i) = q(i,1)
  104 continue
      call trix (1,0,m,ba,bb,bc,b,tcos,d,w)
      do 105 i=1,m
         q(i,1) = b(i)
  105 continue
      go to 183
  106 lr = 0
      do 107 i=1,m
         p(i) = 0.d0
  107 continue
      nun = n
      jst = 1
      jsp = n
c
c     irreg = 1 when no irregularities have occurred, otherwise it is 2.
c
  108 l = 2*jst
      nodd = 2-2*((nun+1)/2)+nun
c
c     nodd = 1 when nun is odd, otherwise it is 2.
c
      go to (110,109),nodd
  109 jsp = jsp-l
      go to 111
  110 jsp = jsp-jst
      if (irreg .ne. 1) jsp = jsp-l
  111 continue
c
c     regular reduction
c
      call cosgen (jst,1,0.5d0,0.0d0,tcos)
      if (l .gt. jsp) go to 118
      do 117 j=l,jsp,l
         jm1 = j-jsh
         jp1 = j+jsh
         jm2 = j-jst
         jp2 = j+jst
         jm3 = jm2-jsh
         jp3 = jp2+jsh
         if (jst .ne. 1) go to 113
         do 112 i=1,m
            b(i) = 2.d0*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  112    continue
         go to 115
  113    do 114 i=1,m
            t = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = t+q(i,j)-q(i,jm3)-q(i,jp3)
            q(i,j) = t
  114    continue
  115    continue
         call trix (jst,0,m,ba,bb,bc,b,tcos,d,w)
         do 116 i=1,m
            q(i,j) = q(i,j)+b(i)
  116    continue
  117 continue
c
c     reduction for last unknown
c
  118 go to (119,136),nodd
  119 go to (152,120),irreg
c
c     odd number of unknowns
c
  120 jsp = jsp+l
      j = jsp
      jm1 = j-jsh
      jp1 = j+jsh
      jm2 = j-jst
      jp2 = j+jst
      jm3 = jm2-jsh
      go to (123,121),istag
  121 continue
      if (jst .ne. 1) go to 123
      do 122 i=1,m
         b(i) = q(i,j)
         q(i,j) = 0.d0
  122 continue
      go to 130
  123 go to (124,126),noddpr
  124 do 125 i=1,m
         ip1 = ip+i
         b(i) = .5d0*(q(i,jm2)-q(i,jm1)-q(i,jm3))+p(ip1)+q(i,j)
  125 continue
      go to 128
  126 do 127 i=1,m
      b(i) = .5d0*(q(i,jm2)-q(i,jm1)-q(i,jm3))+q(i,jp2)-q(i,jp1)+q(i,j)
  127 continue
  128 do 129 i=1,m
         q(i,j) = .5d0*(q(i,j)-q(i,jm1)-q(i,jp1))
  129 continue
  130 call trix (jst,0,m,ba,bb,bc,b,tcos,d,w)
      ip = ip+m
      ipstor = max0(ipstor,ip+m)
      do 131 i=1,m
         ip1 = ip+i
         p(ip1) = q(i,j)+b(i)
         b(i) = q(i,jp2)+p(ip1)
  131 continue
      if (lr .ne. 0) go to 133
      do 132 i=1,jst
         krpi = kr+i
         tcos(krpi) = tcos(i)
  132 continue
      go to 134
  133 continue
      call cosgen (lr,jstsav,0.d0,fi,tcos(jst+1))
      call merge (tcos,0,jst,jst,lr,kr)
  134 continue
      call cosgen (kr,jstsav,0.0d0,fi,tcos)
      call trix (kr,kr,m,ba,bb,bc,b,tcos,d,w)
      do 135 i=1,m
         ip1 = ip+i
         q(i,j) = q(i,jm2)+b(i)+p(ip1)
  135 continue
      lr = kr
      kr = kr+l
      go to 152
c
c     even number of unknowns
c
  136 jsp = jsp+l
      j = jsp
      jm1 = j-jsh
      jp1 = j+jsh
      jm2 = j-jst
      jp2 = j+jst
      jm3 = jm2-jsh
      go to (137,138),irreg
  137 continue
      jstsav = jst
      ideg = jst
      kr = l
      go to 139
  138 call cosgen (kr,jstsav,0.0d00,fi,tcos)
      call cosgen (lr,jstsav,0.0d00,fi,tcos(kr+1))
      ideg = kr
      kr = kr+jst
  139 if (jst .ne. 1) go to 141
      irreg = 2
      do 140 i=1,m
         b(i) = q(i,j)
         q(i,j) = q(i,jm2)
  140 continue
      go to 150
  141 do 142 i=1,m
         b(i) = q(i,j)+.5d00*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  142 continue
      go to (143,145),irreg
  143 do 144 i=1,m
         q(i,j) = q(i,jm2)+.5d00*(q(i,j)-q(i,jm1)-q(i,jp1))
  144 continue
      irreg = 2
      go to 150
  145 continue
      go to (146,148),noddpr
  146 do 147 i=1,m
         ip1 = ip+i
         q(i,j) = q(i,jm2)+p(ip1)
  147 continue
      ip = ip-m
      go to 150
  148 do 149 i=1,m
         q(i,j) = q(i,jm2)+q(i,j)-q(i,jm1)
  149 continue
  150 call trix (ideg,lr,m,ba,bb,bc,b,tcos,d,w)
      do 151 i=1,m
         q(i,j) = q(i,j)+b(i)
  151 continue
  152 nun = nun/2
      noddpr = nodd
      jsh = jst
      jst = 2*jst
      if (nun .ge. 2) go to 108
c
c     start solution.
c
      j = jsp
      do 153 i=1,m
         b(i) = q(i,j)
  153 continue
      go to (154,155),irreg
  154 continue
      call cosgen (jst,1,0.5d00,0.0d00,tcos)
      ideg = jst
      go to 156
  155 kr = lr+jst
      call cosgen (kr,jstsav,0.0d0,fi,tcos)
      call cosgen (lr,jstsav,0.0d0,fi,tcos(kr+1))
      ideg = kr
  156 continue
      call trix (ideg,lr,m,ba,bb,bc,b,tcos,d,w)
      jm1 = j-jsh
      jp1 = j+jsh
      go to (157,159),irreg
  157 do 158 i=1,m
         q(i,j) = .5d0*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
  158 continue
      go to 164
  159 go to (160,162),noddpr
  160 do 161 i=1,m
         ip1 = ip+i
         q(i,j) = p(ip1)+b(i)
  161 continue
      ip = ip-m
      go to 164
  162 do 163 i=1,m
         q(i,j) = q(i,j)-q(i,jm1)+b(i)
  163 continue
  164 continue
c
c     start back substitution.
c
      jst = jst/2
      jsh = jst/2
      nun = 2*nun
      if (nun .gt. n) go to 183
      do 182 j=jst,n,l
         jm1 = j-jsh
         jp1 = j+jsh
         jm2 = j-jst
         jp2 = j+jst
         if (j .gt. jst) go to 166
         do 165 i=1,m
            b(i) = q(i,j)+q(i,jp2)
  165    continue
         go to 170
  166    if (jp2 .le. n) go to 168
         do 167 i=1,m
            b(i) = q(i,j)+q(i,jm2)
  167    continue
         if (jst .lt. jstsav) irreg = 1
         go to (170,171),irreg
  168    do 169 i=1,m
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  169    continue
  170    continue
         call cosgen (jst,1,0.5d0,0.0d0,tcos)
         ideg = jst
         jdeg = 0
         go to 172
  171    if (j+l .gt. n) lr = lr-jst
         kr = jst+lr
         call cosgen (kr,jstsav,0.0d0,fi,tcos)
         call cosgen (lr,jstsav,0.0d0,fi,tcos(kr+1))
         ideg = kr
         jdeg = lr
  172    continue
         call trix (ideg,jdeg,m,ba,bb,bc,b,tcos,d,w)
         if (jst .gt. 1) go to 174
         do 173 i=1,m
            q(i,j) = b(i)
  173    continue
         go to 182
  174    if (jp2 .gt. n) go to 177
  175    do 176 i=1,m
            q(i,j) = .5d0*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
  176    continue
         go to 182
  177    go to (175,178),irreg
  178    if (j+jsh .gt. n) go to 180
         do 179 i=1,m
            ip1 = ip+i
            q(i,j) = b(i)+p(ip1)
  179    continue
         ip = ip-m
         go to 182
  180    do 181 i=1,m
            q(i,j) = b(i)+q(i,j)-q(i,jm1)
  181    continue
  182 continue
      l = l/2
      go to 164
  183 continue
c
c     return storage requirements for p vectors.
c
      w(1) = ipstor
      return
      end
      subroutine poisn2 (m,n,istag,mixbnd,a,bb,c,q,idimq,b,b2,b3,w,w2,
     1                   w3,d,tcos,p)
      implicit real*8(a-h,o-z)
c
c     subroutine to solve poisson's equation with neumann boundary
c     conditions.
c
c     istag = 1 if the last diagonal block is a.
c     istag = 2 if the last diagonal block is a-i.
c     mixbnd = 1 if have neumann boundary conditions at both boundaries.
c     mixbnd = 2 if have neumann boundary conditions at bottom and
c     dirichlet condition at top.  (for this case, must have istag = 1.)
c
      dimension       a(1)       ,bb(1)      ,c(1)       ,q(idimq,1) ,
     1                b(1)       ,b2(1)      ,b3(1)      ,w(1)       ,
     2                w2(1)      ,w3(1)      ,d(1)       ,tcos(1)    ,
     3                k(4)       ,p(1)
      equivalence     (k(1),k1)  ,(k(2),k2)  ,(k(3),k3)  ,(k(4),k4)
      fistag = 3-istag
      fnum = 1.d0/dble(istag)
      fden = 0.5d0*dble(istag-1)
      mr = m
      ip = -mr
      ipstor = 0
      i2r = 1
      jr = 2
      nr = n
      nlast = n
      kr = 1
      lr = 0
      go to (101,103),istag
  101 continue
      do 102 i=1,mr
         q(i,n) = .5d0*q(i,n)
  102 continue
      go to (103,104),mixbnd
  103 if (n .le. 3) go to 155
  104 continue
      jr = 2*i2r
      nrod = 1
      if ((nr/2)*2 .eq. nr) nrod = 0
      go to (105,106),mixbnd
  105 jstart = 1
      go to 107
  106 jstart = jr
      nrod = 1-nrod
  107 continue
      jstop = nlast-jr
      if (nrod .eq. 0) jstop = jstop-i2r
      call cosgen (i2r,1,0.5d0,0.0d0,tcos)
      i2rby2 = i2r/2
      if (jstop .ge. jstart) go to 108
      j = jr
      go to 116
  108 continue
c
c     regular reduction.
c
      do 115 j=jstart,jstop,jr
         jp1 = j+i2rby2
         jp2 = j+i2r
         jp3 = jp2+i2rby2
         jm1 = j-i2rby2
         jm2 = j-i2r
         jm3 = jm2-i2rby2
         if (j .ne. 1) go to 109
         jm1 = jp1
         jm2 = jp2
         jm3 = jp3
  109    continue
         if (i2r .ne. 1) go to 111
         if (j .eq. 1) jm2 = jp2
         do 110 i=1,mr
            b(i) = 2.d0*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  110    continue
         go to 113
  111    continue
         do 112 i=1,mr
            fi = q(i,j)
            q(i,j) = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = fi+q(i,j)-q(i,jm3)-q(i,jp3)
  112    continue
  113    continue
         call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
         do 114 i=1,mr
            q(i,j) = q(i,j)+b(i)
  114    continue
c
c     end of reduction for regular unknowns.
c
  115 continue
c
c     begin special reduction for last unknown.
c
      j = jstop+jr
  116 nlast = j
      jm1 = j-i2rby2
      jm2 = j-i2r
      jm3 = jm2-i2rby2
      if (nrod .eq. 0) go to 128
c
c     odd number of unknowns
c
      if (i2r .ne. 1) go to 118
      do 117 i=1,mr
         b(i) = fistag*q(i,j)
         q(i,j) = q(i,jm2)
  117 continue
      go to 126
  118 do 119 i=1,mr
         b(i) = q(i,j)+.5d0*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  119 continue
      if (nrodpr .ne. 0) go to 121
      do 120 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)
  120 continue
      ip = ip-mr
      go to 123
  121 continue
      do 122 i=1,mr
         q(i,j) = q(i,j)-q(i,jm1)+q(i,jm2)
  122 continue
  123 if (lr .eq. 0) go to 124
      call cosgen (lr,1,0.5d0,fden,tcos(kr+1))
      go to 126
  124 continue
      do 125 i=1,mr
         b(i) = fistag*b(i)
  125 continue
  126 continue
      call cosgen (kr,1,0.5d0,fden,tcos)
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 127 i=1,mr
         q(i,j) = q(i,j)+b(i)
  127 continue
      kr = kr+i2r
      go to 151
  128 continue
c
c     even number of unknowns
c
      jp1 = j+i2rby2
      jp2 = j+i2r
      if (i2r .ne. 1) go to 135
      do 129 i=1,mr
         b(i) = q(i,j)
  129 continue
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      ip = 0
      ipstor = mr
      go to (133,130),istag
  130 do 131 i=1,mr
         p(i) = b(i)
         b(i) = b(i)+q(i,n)
  131 continue
      tcos(1) = 1.d0
      tcos(2) = 0.d0
      call trix (1,1,mr,a,bb,c,b,tcos,d,w)
      do 132 i=1,mr
         q(i,j) = q(i,jm2)+p(i)+b(i)
  132 continue
      go to 150
  133 continue
      do 134 i=1,mr
         p(i) = b(i)
         q(i,j) = q(i,jm2)+2.d0*q(i,jp2)+3.d0*b(i)
  134 continue
      go to 150
  135 continue
      do 136 i=1,mr
         b(i) = q(i,j)+.5d0*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  136 continue
      if (nrodpr .ne. 0) go to 138
      do 137 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  137 continue
      go to 140
  138 continue
      do 139 i=1,mr
         b(i) = b(i)+q(i,jp2)-q(i,jp1)
  139 continue
  140 continue
      call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
      ip = ip+mr
      ipstor = max0(ipstor,ip+mr)
      do 141 i=1,mr
         ii = ip+i
         p(ii) = b(i)+.5d0*(q(i,j)-q(i,jm1)-q(i,jp1))
         b(i) = p(ii)+q(i,jp2)
  141 continue
      if (lr .eq. 0) go to 142
      call cosgen (lr,1,0.5d0,fden,tcos(i2r+1))
      call merge (tcos,0,i2r,i2r,lr,kr)
      go to 144
  142 do 143 i=1,i2r
         ii = kr+i
         tcos(ii) = tcos(i)
  143 continue
  144 call cosgen (kr,1,0.5d0,fden,tcos)
      if (lr .ne. 0) go to 145
      go to (146,145),istag
  145 continue
      call trix (kr,kr,mr,a,bb,c,b,tcos,d,w)
      go to 148
  146 continue
      do 147 i=1,mr
         b(i) = fistag*b(i)
  147 continue
  148 continue
      do 149 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)+b(i)
  149 continue
  150 continue
      lr = kr
      kr = kr+jr
  151 continue
      go to (152,153),mixbnd
  152 nr = (nlast-1)/jr+1
      if (nr .le. 3) go to 155
      go to 154
  153 nr = nlast/jr
      if (nr .le. 1) go to 192
  154 i2r = jr
      nrodpr = nrod
      go to 104
  155 continue
c
c      begin solution
c
      j = 1+jr
      jm1 = j-i2r
      jp1 = j+i2r
      jm2 = nlast-i2r
      if (nr .eq. 2) go to 184
      if (lr .ne. 0) go to 170
      if (n .ne. 3) go to 161
c
c     case n = 3.
c
      go to (156,168),istag
  156 continue
      do 157 i=1,mr
         b(i) = q(i,2)
  157 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 158 i=1,mr
         q(i,2) = b(i)
         b(i) = 4.d0*b(i)+q(i,1)+2.d0*q(i,3)
  158 continue
      tcos(1) = -2.d0
      tcos(2) = 2.d0
      i1 = 2
      i2 = 0
      call trix (i1,i2,mr,a,bb,c,b,tcos,d,w)
      do 159 i=1,mr
         q(i,2) = q(i,2)+b(i)
         b(i) = q(i,1)+2.d0*q(i,2)
  159 continue
      tcos(1) = 0.d0
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 160 i=1,mr
         q(i,1) = b(i)
  160 continue
      jr = 1
      i2r = 0
      go to 194
c
c     case n = 2**p+1
c
  161 continue
      go to (162,170),istag
  162 continue
      do 163 i=1,mr
         b(i) = q(i,j)+.5d0*q(i,1)-q(i,jm1)+q(i,nlast)-q(i,jm2)
  163 continue
      call cosgen (jr,1,0.5d0,0.0d0,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 164 i=1,mr
         q(i,j) = .5d0*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
         b(i) = q(i,1)+2.d0*q(i,nlast)+4.d0*q(i,j)
  164 continue
      jr2 = 2*jr
      call cosgen (jr,1,0.0d0,0.0d0,tcos)
      do 165 i=1,jr
         i1 = jr+i
         i2 = jr+1-i
         tcos(i1) = -tcos(i2)
  165 continue
      call trix (jr2,0,mr,a,bb,c,b,tcos,d,w)
      do 166 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+2.d0*q(i,j)
  166 continue
      call cosgen (jr,1,0.5d0,0.0d0,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 167 i=1,mr
         q(i,1) = .5d0*q(i,1)-q(i,jm1)+b(i)
  167 continue
      go to 194
c
c     case of general n with nr = 3 .
c
  168 do 169 i=1,mr
         b(i) = q(i,2)
         q(i,2) = 0.d0
         b2(i) = q(i,3)
         b3(i) = q(i,1)
  169 continue
      jr = 1
      i2r = 0
      j = 2
      go to 177
  170 continue
      do 171 i=1,mr
         b(i) = .5d0*q(i,1)-q(i,jm1)+q(i,j)
  171 continue
      if (nrod .ne. 0) go to 173
      do 172 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  172 continue
      go to 175
  173 do 174 i=1,mr
         b(i) = b(i)+q(i,nlast)-q(i,jm2)
  174 continue
  175 continue
      do 176 i=1,mr
         t = .5d0*(q(i,j)-q(i,jm1)-q(i,jp1))
         q(i,j) = t
         b2(i) = q(i,nlast)+t
         b3(i) = q(i,1)+2.*t
  176 continue
  177 continue
      k1 = kr+2*jr-1
      k2 = kr+jr
      tcos(k1+1) = -2.
      k4 = k1+3-istag
      call cosgen (k2+istag-2,1,0.0d0,fnum,tcos(k4))
      k4 = k1+k2+1
      call cosgen (jr-1,1,0.0d0,1.0d0,tcos(k4))
      call merge (tcos,k1,k2,k1+k2,jr-1,0)
      k3 = k1+k2+lr
      call cosgen (jr,1,0.5d0,0.0d0,tcos(k3+1))
      k4 = k3+jr+1
      call cosgen (kr,1,0.5d0,fden,tcos(k4))
      call merge (tcos,k3,jr,k3+jr,kr,k1)
      if (lr .eq. 0) go to 178
      call cosgen (lr,1,0.5d0,fden,tcos(k4))
      call merge (tcos,k3,jr,k3+jr,lr,k3-lr)
      call cosgen (kr,1,0.5d0,fden,tcos(k4))
  178 k3 = kr
      k4 = kr
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 179 i=1,mr
         b(i) = b(i)+b2(i)+b3(i)
  179 continue
      tcos(1) = 2.d0
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 180 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+2.*q(i,j)
  180 continue
      call cosgen (jr,1,0.5d0,0.0d0,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      if (jr .ne. 1) go to 182
      do 181 i=1,mr
         q(i,1) = b(i)
  181 continue
      go to 194
  182 continue
      do 183 i=1,mr
         q(i,1) = .5d0*q(i,1)-q(i,jm1)+b(i)
  183 continue
      go to 194
  184 continue
      if (n .ne. 2) go to 188
c
c     case  n = 2
c
      do 185 i=1,mr
         b(i) = q(i,1)
  185 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 186 i=1,mr
         q(i,1) = b(i)
         b(i) = 2.d0*(q(i,2)+b(i))*fistag
  186 continue
      tcos(1) = -fistag
      tcos(2) = 2.d0
      call trix (2,0,mr,a,bb,c,b,tcos,d,w)
      do 187 i=1,mr
         q(i,1) = q(i,1)+b(i)
  187 continue
      jr = 1
      i2r = 0
      go to 194
  188 continue
c
c     case of general n and nr = 2 .
c
      do 189 i=1,mr
         ii = ip+i
         b3(i) = 0.
         b(i) = q(i,1)+2.*p(ii)
         q(i,1) = .5d0*q(i,1)-q(i,jm1)
         b2(i) = 2.d0*(q(i,1)+q(i,nlast))
  189 continue
      k1 = kr+jr-1
      tcos(k1+1) = -2.d0
      k4 = k1+3-istag
      call cosgen (kr+istag-2,1,0.0d0,fnum,tcos(k4))
      k4 = k1+kr+1
      call cosgen (jr-1,1,0.0d0,1.0d0,tcos(k4))
      call merge (tcos,k1,kr,k1+kr,jr-1,0)
      call cosgen (kr,1,0.5d0,fden,tcos(k1+1))
      k2 = kr
      k4 = k1+k2+1
      call cosgen (lr,1,0.5d0,fden,tcos(k4))
      k3 = lr
      k4 = 0
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 190 i=1,mr
         b(i) = b(i)+b2(i)
  190 continue
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 191 i=1,mr
         q(i,1) = q(i,1)+b(i)
  191 continue
      go to 194
  192 do 193 i=1,mr
         b(i) = q(i,nlast)
  193 continue
      go to 196
  194 continue
c
c     start back substitution.
c
      j = nlast-jr
      do 195 i=1,mr
         b(i) = q(i,nlast)+q(i,j)
  195 continue
  196 jm2 = nlast-i2r
      if (jr .ne. 1) go to 198
      do 197 i=1,mr
         q(i,nlast) = 0.
  197 continue
      go to 202
  198 continue
      if (nrod .ne. 0) go to 200
      do 199 i=1,mr
         ii = ip+i
         q(i,nlast) = p(ii)
  199 continue
      ip = ip-mr
      go to 202
  200 do 201 i=1,mr
         q(i,nlast) = q(i,nlast)-q(i,jm2)
  201 continue
  202 continue
      call cosgen (kr,1,0.5d0,fden,tcos)
      call cosgen (lr,1,0.5d0,fden,tcos(kr+1))
      if (lr .ne. 0) go to 204
      do 203 i=1,mr
         b(i) = fistag*b(i)
  203 continue
  204 continue
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 205 i=1,mr
         q(i,nlast) = q(i,nlast)+b(i)
  205 continue
      nlastp = nlast
  206 continue
      jstep = jr
      jr = i2r
      i2r = i2r/2
      if (jr .eq. 0) go to 222
      go to (207,208),mixbnd
  207 jstart = 1+jr
      go to 209
  208 jstart = jr
  209 continue
      kr = kr-jr
      if (nlast+jr .gt. n) go to 210
      kr = kr-jr
      nlast = nlast+jr
      jstop = nlast-jstep
      go to 211
  210 continue
      jstop = nlast-jr
  211 continue
      lr = kr-jr
      call cosgen (jr,1,0.5d0,0.0d0,tcos)
      do 221 j=jstart,jstop,jstep
         jm2 = j-jr
         jp2 = j+jr
         if (j .ne. jr) go to 213
         do 212 i=1,mr
            b(i) = q(i,j)+q(i,jp2)
  212    continue
         go to 215
  213    continue
         do 214 i=1,mr
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  214    continue
  215    continue
         if (jr .ne. 1) go to 217
         do 216 i=1,mr
            q(i,j) = 0.
  216    continue
         go to 219
  217    continue
         jm1 = j-i2r
         jp1 = j+i2r
         do 218 i=1,mr
            q(i,j) = .5d0*(q(i,j)-q(i,jm1)-q(i,jp1))
  218    continue
  219    continue
         call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
         do 220 i=1,mr
            q(i,j) = q(i,j)+b(i)
  220    continue
  221 continue
      nrod = 1
      if (nlast+i2r .le. n) nrod = 0
      if (nlastp .ne. nlast) go to 194
      go to 206
  222 continue
c
c     return storage requirements for p vectors.
c
      w(1) = ipstor
      return
      end
      subroutine poisp2 (m,n,a,bb,c,q,idimq,b,b2,b3,w,w2,w3,d,tcos,p)
      implicit real*8(a-h,o-z)
c
c     subroutine to solve poisson equation with periodic boundary
c     conditions.
c
      dimension       a(1)       ,bb(1)      ,c(1)       ,q(idimq,1) ,
     1                b(1)       ,b2(1)      ,b3(1)      ,w(1)       ,
     2                w2(1)      ,w3(1)      ,d(1)       ,tcos(1)    ,
     3                p(1)
      mr = m
      nr = (n+1)/2
      nrm1 = nr-1
      if (2*nr .ne. n) go to 107
c
c     even number of unknowns
c
      do 102 j=1,nrm1
         nrmj = nr-j
         nrpj = nr+j
         do 101 i=1,mr
            s = q(i,nrmj)-q(i,nrpj)
            t = q(i,nrmj)+q(i,nrpj)
            q(i,nrmj) = s
            q(i,nrpj) = t
  101    continue
  102 continue
      do 103 i=1,mr
         q(i,nr) = 2.d0*q(i,nr)
         q(i,n) = 2.d0*q(i,n)
  103 continue
      call poisd2 (mr,nrm1,1,a,bb,c,q,idimq,b,w,d,tcos,p)
      ipstor = w(1)
      call poisn2 (mr,nr+1,1,1,a,bb,c,q(1,nr),idimq,b,b2,b3,w,w2,w3,d,
     1             tcos,p)
      ipstor = max0(ipstor,int(w(1)))
      do 105 j=1,nrm1
         nrmj = nr-j
         nrpj = nr+j
         do 104 i=1,mr
            s = .5d0*(q(i,nrpj)+q(i,nrmj))
            t = .5d0*(q(i,nrpj)-q(i,nrmj))
            q(i,nrmj) = s
            q(i,nrpj) = t
  104    continue
  105 continue
      do 106 i=1,mr
         q(i,nr) = .5d0*q(i,nr)
         q(i,n) = .5d0*q(i,n)
  106 continue
      go to 118
  107 continue
c
c     odd  number of unknowns
c
      do 109 j=1,nrm1
         nrpj = n+1-j
         do 108 i=1,mr
            s = q(i,j)-q(i,nrpj)
            t = q(i,j)+q(i,nrpj)
            q(i,j) = s
            q(i,nrpj) = t
  108    continue
  109 continue
      do 110 i=1,mr
         q(i,nr) = 2.d0*q(i,nr)
  110 continue
      lh = nrm1/2
      do 112 j=1,lh
         nrmj = nr-j
         do 111 i=1,mr
            s = q(i,j)
            q(i,j) = q(i,nrmj)
            q(i,nrmj) = s
  111    continue
  112 continue
      call poisd2 (mr,nrm1,2,a,bb,c,q,idimq,b,w,d,tcos,p)
      ipstor = w(1)
      call poisn2 (mr,nr,2,1,a,bb,c,q(1,nr),idimq,b,b2,b3,w,w2,w3,d,
     1             tcos,p)
      ipstor = max0(ipstor,int(w(1)))
      do 114 j=1,nrm1
         nrpj = nr+j
         do 113 i=1,mr
            s = .5d0*(q(i,nrpj)+q(i,j))
            t = .5d0*(q(i,nrpj)-q(i,j))
            q(i,nrpj) = t
            q(i,j) = s
  113    continue
  114 continue
      do 115 i=1,mr
         q(i,nr) = .5d0*q(i,nr)
  115 continue
      do 117 j=1,lh
         nrmj = nr-j
         do 116 i=1,mr
            s = q(i,j)
            q(i,j) = q(i,nrmj)
            q(i,nrmj) = s
  116    continue
  117 continue
  118 continue
c
c     return storage requirements for p vectors.
c
      w(1) = ipstor
      return
      end
      subroutine tri3 (m,a,b,c,k,y1,y2,y3,tcos,d,w1,w2,w3)
      implicit real*8(a-h,o-z)
      dimension       a(1)       ,b(1)       ,c(1)       ,k(4)       ,
     1                tcos(1)    ,y1(1)      ,y2(1)      ,y3(1)      ,
     2                d(1)       ,w1(1)      ,w2(1)      ,w3(1)
c
c     subroutine to solve three linear systems whose common coefficient
c     matrix is a rational function in the matrix given by
c
c                  tridiagonal (...,a(i),b(i),c(i),...)
c
      mm1 = m-1
      k1 = k(1)
      k2 = k(2)
      k3 = k(3)
      k4 = k(4)
      f1 = k1+1
      f2 = k2+1
      f3 = k3+1
      f4 = k4+1
      k2k3k4 = k2+k3+k4
      if (k2k3k4 .eq. 0) go to 101
      l1 = f1/f2
      l2 = f1/f3
      l3 = f1/f4
      lint1 = 1
      lint2 = 1
      lint3 = 1
      kint1 = k1
      kint2 = kint1+k2
      kint3 = kint2+k3
  101 continue
      do 115 n=1,k1
         x = tcos(n)
         if (k2k3k4 .eq. 0) go to 107
         if (n .ne. l1) go to 103
         do 102 i=1,m
            w1(i) = y1(i)
  102    continue
  103    if (n .ne. l2) go to 105
         do 104 i=1,m
            w2(i) = y2(i)
  104    continue
  105    if (n .ne. l3) go to 107
         do 106 i=1,m
            w3(i) = y3(i)
  106    continue
  107    continue
         z = 1.d0/(b(1)-x)
         d(1) = c(1)*z
         y1(1) = y1(1)*z
         y2(1) = y2(1)*z
         y3(1) = y3(1)*z
         do 108 i=2,m
            z = 1.d0/(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y1(i) = (y1(i)-a(i)*y1(i-1))*z
            y2(i) = (y2(i)-a(i)*y2(i-1))*z
            y3(i) = (y3(i)-a(i)*y3(i-1))*z
  108    continue
         do 109 ip=1,mm1
            i = m-ip
            y1(i) = y1(i)-d(i)*y1(i+1)
            y2(i) = y2(i)-d(i)*y2(i+1)
            y3(i) = y3(i)-d(i)*y3(i+1)
  109    continue
         if (k2k3k4 .eq. 0) go to 115
         if (n .ne. l1) go to 111
         i = lint1+kint1
         xx = x-tcos(i)
         do 110 i=1,m
            y1(i) = xx*y1(i)+w1(i)
  110    continue
         lint1 = lint1+1
         l1 = (dble(lint1)*f1)/f2
  111    if (n .ne. l2) go to 113
         i = lint2+kint2
         xx = x-tcos(i)
         do 112 i=1,m
            y2(i) = xx*y2(i)+w2(i)
  112    continue
         lint2 = lint2+1
         l2 = (dble(lint2)*f1)/f3
  113    if (n .ne. l3) go to 115
         i = lint3+kint3
         xx = x-tcos(i)
         do 114 i=1,m
            y3(i) = xx*y3(i)+w3(i)
  114    continue
         lint3 = lint3+1
         l3 = (dble(lint3)*f1)/f4
  115 continue
      return
      end
      subroutine trix (idegbr,idegcr,m,a,b,c,y,tcos,d,w)
      implicit real*8(a-h,o-z)
c
c     subroutine to solve a system of linear equations where the
c     coefficient matrix is a rational function in the matrix given by
c     tridiagonal  ( . . . , a(i), b(i), c(i), . . . ).
c
      dimension       a(1)       ,b(1)       ,c(1)       ,y(1)       ,
     1                tcos(1)    ,d(1)       ,w(1)
      mm1 = m-1
      fb = idegbr+1
      fc = idegcr+1
      l = fb/fc
      lint = 1
      do 108 k=1,idegbr
         x = tcos(k)
         if (k .ne. l) go to 102
         i = idegbr+lint
         xx = x-tcos(i)
         do 101 i=1,m
            w(i) = y(i)
            y(i) = xx*y(i)
  101    continue
  102    continue
         z = 1.d0/(b(1)-x)
         d(1) = c(1)*z
         y(1) = y(1)*z
         do 103 i=2,mm1
            z = 1.d0/(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y(i) = (y(i)-a(i)*y(i-1))*z
  103    continue
         z = b(m)-x-a(m)*d(mm1)
         if (z .ne. 0.d0) go to 104
         y(m) = 0.d0
         go to 105
  104    y(m) = (y(m)-a(m)*y(mm1))/z
  105    continue
         do 106 ip=1,mm1
            i = m-ip
            y(i) = y(i)-d(i)*y(i+1)
  106    continue
         if (k .ne. l) go to 108
         do 107 i=1,m
            y(i) = y(i)+w(i)
  107    continue
         lint = lint+1
         l = (dble(lint)*fb)/fc
  108 continue
      return
      end
      subroutine cosgen (n,ijump,fnum,fden,a)
      implicit real*8(a-h,o-z)
      dimension       a(1)
c
c
c     this subroutine computes required cosine values in ascending
c     order.  when ijump .gt. 1 the routine computes values
c
c        2*cos(j*pi/l) , j=1,2,...,l and j .ne. 0(mod n/ijump+1)
c
c     where l = ijump*(n/ijump+1).
c
c
c     when ijump = 1 it computes
c
c            2*cos((j-fnum)*pi/(n+fden)) ,  j=1, 2, ... ,n
c
c     where
c        fnum = 0.5, fden = 0.0,  for regular reduction values
c        fnum = 0.0, fden = 1.0, for b-r and c-r when istag = 1
c        fnum = 0.0, fden = 0.5, for b-r and c-r when istag = 2
c        fnum = 0.5, fden = 0.5, for b-r and c-r when istag = 2
c                                in poisn2 only.
c
c
      pi = pimach(dum)
      if (n .eq. 0) go to 105
      if (ijump .eq. 1) go to 103
      k3 = n/ijump+1
      k4 = k3-1
      pibyn = pi/dble(n+ijump)
      do 102 k=1,ijump
         k1 = (k-1)*k3
         k5 = (k-1)*k4
         do 101 i=1,k4
            x = k1+i
            k2 = k5+i
            a(k2) = -2.d0*dcos(x*pibyn)
  101    continue
  102 continue
      go to 105
  103 continue
      np1 = n+1
      y = pi/(dble(n)+fden)
      do 104 i=1,n
         x = dble(np1-i)-fnum
         a(i) = 2.d0*dcos(x*y)
  104 continue
  105 continue
      return
      end
      subroutine merge (tcos,i1,m1,i2,m2,i3)
      implicit real*8(a-h,o-z)
      dimension       tcos(1)
c
c
c     this subroutine merges two ascending strings of numbers in the
c     array tcos.  the first string is of length m1 and starts at
c     tcos(i1+1).  the second string is of length m2 and starts at
c     tcos(i2+1).  the merged string goes into tcos(i3+1).
c
c
      j1 = 1
      j2 = 1
      j = i3
      if (m1 .eq. 0) go to 107
      if (m2 .eq. 0) go to 104
  101 j = j+1
      l = j1+i1
      x = tcos(l)
      l = j2+i2
      y = tcos(l)
      if (x-y) 102,102,103
  102 tcos(j) = x
      j1 = j1+1
      if (j1 .gt. m1) go to 106
      go to 101
  103 tcos(j) = y
      j2 = j2+1
      if (j2 .le. m2) go to 101
      if (j1 .gt. m1) go to 109
  104 k = j-j1+1
      do 105 j=j1,m1
         m = k+j
         l = j+i1
         tcos(m) = tcos(l)
  105 continue
      go to 109
  106 continue
      if (j2 .gt. m2) go to 109
  107 k = j-j2+1
      do 108 j=j2,m2
         m = k+j
         l = j+i2
         tcos(m) = tcos(l)
  108 continue
  109 continue
      return
      end
c      double precision function pimach (dum)
c      implicit real*8(a-h,o-z)
c
c     this subprogram supplies the value of the constant pi correct to
c     machine precision where
c
c     pi=3.1415926535897932384626433832795028841971693993751058209749446
c
c      pimach = 3.14159265358979d00
c      return
c      end
c

