!*******************************************************************************
!                              INTEL CONFIDENTIAL
!   Copyright(C) 2005-2008 Intel Corporation. All Rights Reserved.
!   The source code contained  or  described herein and all documents related to
!   the source code ("Material") are owned by Intel Corporation or its suppliers
!   or licensors.  Title to the  Material remains with  Intel Corporation or its
!   suppliers and licensors. The Material contains trade secrets and proprietary
!   and  confidential  information of  Intel or its suppliers and licensors. The
!   Material  is  protected  by  worldwide  copyright  and trade secret laws and
!   treaty  provisions. No part of the Material may be used, copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way without Intel's prior express written permission.
!   No license  under any  patent, copyright, trade secret or other intellectual
!   property right is granted to or conferred upon you by disclosure or delivery
!   of the Materials,  either expressly, by implication, inducement, estoppel or
!   otherwise.  Any  license  under  such  intellectual property  rights must be
!   express and approved by Intel in writing.
C*******************************************************************************
C  Content:
C  Intel MKL RCI (P)FGMRES ((Preconditioned) Flexible Generalized Minimal
C                                                       RESidual method) example
C*******************************************************************************

C---------------------------------------------------------------------------
C  Example program for solving non-symmetric indefinite system of equations
C  Simplest case: no preconditioning and no user-defined stopping tests
C---------------------------------------------------------------------------

      SUBROUTINE FGMRES(N,N_RESTRT,COMPUTED_SOLUTION,RHS,tol,
     $ MatrixVector,PSolve)

C      INCLUDE "mkl_rci.fi"

      INTEGER N,N_RESTRT
      INTEGER SIZE
      PARAMETER (SIZE=128)

      real tol
      external MatrixVector,PSolve

C---------------------------------------------------------------------------
C Define arrays for the upper triangle of the coefficient matrix
C Compressed sparse row storage is used for sparse representation
C---------------------------------------------------------------------------

C      INTEGER IA(6)
C      DATA IA /1,3,6,9,12,14/
C      INTEGER JA(13)
C      DATA JA    /  1,        3,
C     1              1,   2,        4,
C     2						 2,   3,        5,
C     3 								3,   4,   5,
C     4 									  4,   5  /
C      DOUBLE PRECISION A(13)
C      DATA A     / 1.0,     -1.0,
C     1 		  	   -1.0, 1.0,     -1.0,
C     2 				      1.0,-2.0,      1.0,
C     3 				          -1.0, 2.0,-1.0,
C     4 				               -1.0,-3.0 /

C---------------------------------------------------------------------------
C Allocate storage for the ?par parameters and the solution/rhs vectors
C---------------------------------------------------------------------------

      INTEGER IPAR(SIZE)
      DOUBLE PRECISION DPAR(SIZE)
      DOUBLE PRECISION TMP(N*(2*N_RESTRT+1)+(N_RESTRT*(N_RESTRT+9))/2+1)
      DOUBLE PRECISION RHS(N)
      DOUBLE PRECISION COMPUTED_SOLUTION(N)

C      DOUBLE PRECISION EXPECTED_SOLUTION(N)
C      DATA EXPECTED_SOLUTION /-1.0,1.0,0.0,1.0,-1.0/

C---------------------------------------------------------------------------
C Some additional variables to use with the RCI (P)FGMRES solver
C---------------------------------------------------------------------------
      INTEGER ITERCOUNT
      INTEGER RCI_REQUEST, I

C---------------------------------------------------------------------------
C Initialize variables and the right hand side through matrix-vector product
C---------------------------------------------------------------------------
C      CALL MKL_DCSRGEMV('N', N, A, IA, JA, EXPECTED_SOLUTION, RHS)

C---------------------------------------------------------------------------
C Initialize the initial guess
C---------------------------------------------------------------------------

C      DO I=1,N
C         COMPUTED_SOLUTION(I)=1.0
C      ENDDO

C---------------------------------------------------------------------------
C Initialize the solver
C---------------------------------------------------------------------------


      CALL DFGMRES_INIT(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR,
     1 DPAR, TMP)

      IF (RCI_REQUEST.NE.0) GOTO 999

C---------------------------------------------------------------------------
C Set the desired parameters:
C LOGICAL parameters:
C do residual stopping test
C do not request for the user defined stopping test
C do the check of the norm of the next generated vector automatically
C DOUBLE PRECISION parameters
C set the relative tolerance to 1.0D-3 instead of default value 1.0D-6
C---------------------------------------------------------------------------

      IPAR(9)=1
      IPAR(8)=0
      IPAR(10)=0
      IPAR(11)=1
      IPAR(12)=1
      DPAR(1)=tol
      IPAR(15)=N_RESTRT

C---------------------------------------------------------------------------
C Check the correctness and consistency of the newly set parameters
C---------------------------------------------------------------------------
      CALL DFGMRES_CHECK(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST,
     1 IPAR, DPAR, TMP)
      IF (RCI_REQUEST.NE.0) GOTO 999


C---------------------------------------------------------------------------
C Print the info about the RCI FGMRES method
C---------------------------------------------------------------------------

C      PRINT *, ''
C      PRINT *,'Some info about the current run of RCI FGMRES method:'
C      PRINT *, ''
C      IF (IPAR(8).NE.0) THEN
C         WRITE(*,'(A,I1,A,A)') 'As IPAR(8)=',IPAR(8),', the automatic',
C     1 ' test for the maximal number of iterations will be'
C         PRINT *,'performed'
C      ELSE
C      	WRITE(*,'(A,I1,A,A)') 'As IPAR(8)=',IPAR(8),', the automatic',
C     1 ' test for the maximal number of iterations will be'
C      	PRINT *,'skipped'
C      ENDIF
C      PRINT *,'+++'
C      IF (IPAR(9).NE.0) THEN
C      	WRITE(*,'(A,I1,A,A)') 'As IPAR(9)=',IPAR(9),', the automatic',
C     1 ' residual test will be performed'
C      ELSE
C      	WRITE(*,'(A,I1,A,A)') 'As IPAR(9)=',IPAR(9),', the automatic',
C     1 ' residual test will be skipped'
C      ENDIF
C      PRINT *,'+++'
C      IF (IPAR(10).NE.0) THEN
C      	WRITE(*,'(A,I1,A,A)') 'As IPAR(10)=',IPAR(10),', the',
C     1 ' user-defined stopping test will be requested via'
C      	PRINT *,'RCI_REQUEST=2'
C      ELSE
C      	WRITE(*,'(A,I1,A,A)') 'As IPAR(10)=',IPAR(10),', the',
C     1 ' user-defined stopping test will not be requested, thus,'
C      	PRINT *,'RCI_REQUEST will not take the value 2'
C      ENDIF
C      PRINT *,'+++'
C      IF (IPAR(11).NE.0) THEN
C      	WRITE(*,'(A,I1,A,A)') 'As IPAR(11)=',IPAR(11),', the',
C     1 ' Preconditioned FGMRES iterations will be performed, thus,'
C      	WRITE(*,'(A,A)') 'the preconditioner action will be requested',
C     1 ' via RCI_REQUEST=3'
C      ELSE
C      	WRITE(*,'(A,I1,A,A)') 'As IPAR(11)=',IPAR(11),', the',
C     1 ' Preconditioned FGMRES iterations will not be performed,'
C      	WRITE(*,'(A)') 'thus, RCI_REQUEST will not take the value 3'
C      ENDIF
C      PRINT *,'+++'
C      IF (IPAR(12).NE.0) THEN
C      	WRITE(*,'(A,I1,A,A)') 'As IPAR(12)=',IPAR(12),', the automatic',
C     1 ' test for the norm of the next generated vector is'
C      	WRITE(*,'(A,A)') 'not equal to zero up to rounding and',
C     1 ' computational errors will be performed,'
C      	PRINT *,'thus, RCI_REQUEST will not take the value 4'
C      ELSE
C      	WRITE(*,'(A,I1,A,A)') 'As IPAR(12)=',IPAR(12),', the automatic',
C     1 ' test for the norm of the next generated vector is'
C      	WRITE(*,'(A,A)') 'not equal to zero up to rounding and',
C     1 ' computational errors will be skipped,'
C      	WRITE(*,'(A,A)') 'thus, the user-defined test will be requested',
C     1 ' via RCI_REQUEST=4'
C      ENDIF
C      PRINT *,'+++'
C---------------------------------------------------------------------------
C Compute the solution by RCI (P)FGMRES solver without preconditioning
C Reverse Communication starts here
C---------------------------------------------------------------------------


1     CALL DFGMRES(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR,
     1 DPAR, TMP)

C---------------------------------------------------------------------------
C If RCI_REQUEST=0, then the solution was found with the required precision
C---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.0) GOTO 3


C---------------------------------------------------------------------------
C If RCI_REQUEST=1, then compute the vector A*TMP(IPAR(22))
C and put the result in vector TMP(IPAR(23))
C---------------------------------------------------------------------------

      IF (RCI_REQUEST.EQ.1) THEN
      	CALL MatrixVector(TMP(IPAR(22)), TMP(IPAR(23)))
      	GOTO 1
      ENDIF

C---------------------------------------------------------------------------
C preconditioner by FFT
C---------------------------------------------------------------------------

      IF (RCI_REQUEST.EQ.3) THEN
      	CALL PSolve(TMP(IPAR(22)), TMP(IPAR(23)))
      	GOTO 1
      ELSE
      	GOTO 999
      ENDIF



C---------------------------------------------------------------------------
C Reverse Communication ends here
C Get the current iteration number and the FGMRES solution (DO NOT FORGET to
C call DFGMRES_GET routine as COMPUTED_SOLUTION is still containing
C the initial guess!)
C---------------------------------------------------------------------------
3     CALL DFGMRES_GET(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR,
     1 DPAR, TMP, ITERCOUNT)



C---------------------------------------------------------------------------
C Print solution vector: COMPUTED_SOLUTION(N) and
C the number of iterations: ITERCOUNT
C---------------------------------------------------------------------------
C      PRINT *, ''
C      PRINT *,' The system has been SUCCESSFULLY solved'
C      PRINT *, ''
C      PRINT *,' The following solution has been obtained:'
C      DO I=1,N
C         WRITE(*,'(A18,I1,A2,E10.3)') 'COMPUTED_SOLUTION(',I,')=',
C     1 COMPUTED_SOLUTION(I)
C      ENDDO
C      PRINT *, ''
C      PRINT *,' The expected solution is:'
C      DO I=1,N
C         WRITE(*,'(A18,I1,A2,E10.3)') 'EXPECTED_SOLUTION(',I,')=',
C     1 EXPECTED_SOLUTION(I)
C      ENDDO
C      PRINT *, ''
C      PRINT *,' Number of iterations: ',ITERCOUNT
      WRITE(*,*) 'FGMRES', ITERCOUNT
      GOTO 1000

999   PRINT *,'The solver has returned the ERROR code ', RCI_REQUEST
      STOP

1000  CONTINUE
      END
