* -*- fortran -*-
*     GetBreak
      double precision FUNCTION GETBREAK()
*
*     Get breakdown parameter tolerance; for the test routine,
*     set to machine precision.
*
      double precision EPS
      double precision,external::DLAMCH
*
      EPS = DLAMCH('EPS')
      GETBREAK = EPS**2
*
      RETURN
*
      END FUNCTION



