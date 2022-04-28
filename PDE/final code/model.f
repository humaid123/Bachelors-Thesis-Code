
C-----------------------------------------------------------------------
      SUBROUTINE F_time(T, X, U, UX, UXX, FVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE DEFINES THE RIGHT HAND SIDE VECTOR OF THE 
C       NPDE DIMENSIONAL PARABOLIC PARTIAL DIFFERENTIAL EQUATION 
C                        UT = F(T, X, U, UX, UXX). 
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
         DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
      DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
      DOUBLE PRECISION        UXX(NPDE)
C                               UXX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SECOND SPATIAL DERIVATIVE OF THE 
C                               SOLUTION AT THE POINT (T,X).
C
C OUTPUT:
      DOUBLE PRECISION        FVAL(NPDE)
C                               FVAL(1:NPDE) IS THE RIGHT HAND SIDE
C                               VECTOR F(T, X, U, UX, UXX) OF THE PDE.
  !DOUBLE PRECISION PI
  !PI = ACOS(-1.D0)
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE 
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      DOUBLE PRECISION mu, beta, alpha, gamma, N
      DOUBLE PRECISION dS, dE, dI, dR

      integer :: nfev
      common /nfev/ nfev
      nfev = nfev + 1 ! count the number of function evaluations
      ! write(12, *) t, nfev

      alpha = (1.0/8.0)
      beta = 0.90
      gamma = 0.06
      mu = (0.01/365.0)
      N = 37d6 

      if (T > 30d0) then
            beta = 0.035
      end if

      call dSfunc(X, dS)
      call dEfunc(X, dE)
      call dIfunc(X, dI)
      call dRfunc(X, dR)

      FVAL(1) = dS*Uxx(1) + MU* N - MU*U(1) - (beta/N)*U(3)*U(1)
      FVAL(2) = dE*Uxx(2)+(beta/N)*U(3)*U(1)-alpha*U(2)-mu*U(2)
      FVAL(3) = dI*Uxx(3) + alpha*U(2) - gamma*U(3) - mu*U(3)
      FVAL(4) = dR*Uxx(4) + gamma*U(3) - mu*U(4)
C
      RETURN
      END



C-----------------------------------------------------------------------
      SUBROUTINE F_state(T, X, U, UX, UXX, FVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE DEFINES THE RIGHT HAND SIDE VECTOR OF THE 
C       NPDE DIMENSIONAL PARABOLIC PARTIAL DIFFERENTIAL EQUATION 
C                        UT = F(T, X, U, UX, UXX). 
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
            INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
            DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
      DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
      DOUBLE PRECISION        UXX(NPDE)
C                               UXX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SECOND SPATIAL DERIVATIVE OF THE 
C                               SOLUTION AT THE POINT (T,X).
C
C OUTPUT:
      DOUBLE PRECISION        FVAL(NPDE)
C                               FVAL(1:NPDE) IS THE RIGHT HAND SIDE
C                               VECTOR F(T, X, U, UX, UXX) OF THE PDE.
      !DOUBLE PRECISION PI
      !PI = ACOS(-1.D0)
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE 
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C

      DOUBLE PRECISION mu, beta, alpha, gamma, N
      DOUBLE PRECISION dS, dE, dI, dR
      integer :: measures_implemented
      common /measures_implemented/ measures_implemented


      integer :: nfev
      common /nfev/ nfev
      nfev = nfev + 1 ! count the number of function evaluations

      if (measures_implemented == 1) then
            beta = 0.035
      else
            beta = 0.9
      end if

      alpha = (1.0/8.0)
      gamma = 0.06
      mu = (0.01/365.0)
      N = 37d6

      call dSfunc(X, dS)
      call dEfunc(X, dE)
      call dIfunc(X, dI)
      call dRfunc(X, dR)

      FVAL(1) = dS*Uxx(1) + MU* N - MU*U(1) - (beta/N)*U(3)*U(1)
      FVAL(2) = dE*Uxx(2)+(beta/N)*U(3)*U(1)-alpha*U(2)-mu*U(2)
      FVAL(3) = dI*Uxx(3) + alpha*U(2) - gamma*U(3) - mu*U(3)
      FVAL(4) = dR*Uxx(4) + gamma*U(3) - mu*U(4)
C
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE F1(T, X, U, UX, UXX, FVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE DEFINES THE RIGHT HAND SIDE VECTOR OF THE 
C       NPDE DIMENSIONAL PARABOLIC PARTIAL DIFFERENTIAL EQUATION 
C                        UT = F(T, X, U, UX, UXX). 
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
            INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
            DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
      DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
      DOUBLE PRECISION        UXX(NPDE)
C                               UXX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SECOND SPATIAL DERIVATIVE OF THE 
C                               SOLUTION AT THE POINT (T,X).
C
C OUTPUT:
      DOUBLE PRECISION        FVAL(NPDE)
C                               FVAL(1:NPDE) IS THE RIGHT HAND SIDE
C                               VECTOR F(T, X, U, UX, UXX) OF THE PDE.
      !DOUBLE PRECISION PI
      !PI = ACOS(-1.D0)
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE 
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C

      DOUBLE PRECISION mu, beta, alpha, gamma, N
      DOUBLE PRECISION dS, dE, dI, dR
      integer :: nfev
      common /nfev/ nfev
      nfev = nfev + 1 ! count the number of function evaluations
      ! write(12, *) t, nfev
      
      beta = 0.90d0
      alpha = (1.0/8.0)
      gamma = 0.06
      mu = (0.01/365.0)
      N = 37d6

      call dSfunc(X, dS)
      call dEfunc(X, dE)
      call dIfunc(X, dI)
      call dRfunc(X, dR)

      FVAL(1) = dS*Uxx(1) + MU* N - MU*U(1) - (beta/N)*U(3)*U(1)
      FVAL(2) = dE*Uxx(2)+(beta/N)*U(3)*U(1)-alpha*U(2)-mu*U(2)
      FVAL(3) = dI*Uxx(3) + alpha*U(2) - gamma*U(3) - mu*U(3)
      FVAL(4) = dR*Uxx(4) + gamma*U(3) - mu*U(4)
C
      RETURN
      END


C-----------------------------------------------------------------------
      SUBROUTINE F2(T, X, U, UX, UXX, FVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE DEFINES THE RIGHT HAND SIDE VECTOR OF THE 
C       NPDE DIMENSIONAL PARABOLIC PARTIAL DIFFERENTIAL EQUATION 
C                        UT = F(T, X, U, UX, UXX). 
C		IT DOES SO FOR THE CASE AFTER MEASURES ARE IMPLEMENTED
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
         DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
      DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
      DOUBLE PRECISION        UXX(NPDE)
C                               UXX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SECOND SPATIAL DERIVATIVE OF THE 
C                               SOLUTION AT THE POINT (T,X).
C
C OUTPUT:
      DOUBLE PRECISION        FVAL(NPDE)
C                               FVAL(1:NPDE) IS THE RIGHT HAND SIDE
C                               VECTOR F(T, X, U, UX, UXX) OF THE PDE.
  !DOUBLE PRECISION PI
  !PI = ACOS(-1.D0)
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE 
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C

      DOUBLE PRECISION mu, beta, alpha, gamma, N
      DOUBLE PRECISION dS, dE, dI, dR
      integer :: nfev
      common /nfev/ nfev
      nfev = nfev + 1 ! count the number of function evaluations
      ! write(12, *) t, nfev
      
      alpha = (1.0/8.0)
      beta = 0.035
      gamma = 0.06
      mu = (0.01/365.0)
      N = 37d6

      call dSfunc(X, dS)
      call dEfunc(X, dE)
      call dIfunc(X, dI)
      call dRfunc(X, dR)

      FVAL(1) = dS*Uxx(1)+MU*N-MU*U(1)-(beta/N)*U(3)*U(1)
      FVAL(2) = dE*Uxx(2)+(beta/N)*U(3)*U(1)-alpha*U(2)-mu*U(2)
      FVAL(3) = dI*Uxx(3) + alpha*U(2) - gamma*U(3) - mu*U(3)
      FVAL(4) = dR*Uxx(4) + gamma*U(3) - mu*U(4)
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE dSfunc(X, dS)
      DOUBLE PRECISION X
      DOUBLE PRECISION dS
        
      ! OLD values => I change the values to make the growth faster
      ! DOUBLE PRECISION, parameter :: maxDs = 0.05
      ! DOUBLE PRECISION, parameter :: minDs = 0.001!0.035
      ! DOUBLE PRECISION, parameter :: var1 = 8
      ! DOUBLE PRECISION, parameter :: var2 = 0.5

      ! NEW from DOUBLE PRECISION, parameter :: maxDs = 0.1
      ! NEW from !DOUBLE PRECISION, parameter :: minDs = 0.05
      ! worked with this => DOUBLE PRECISION, parameter :: maxDs = 0.8
      ! worked with this => DOUBLE PRECISION, parameter :: minDs = 0.1
      DOUBLE PRECISION, parameter :: maxDs = 0.8
      DOUBLE PRECISION, parameter :: minDs = 0.01
      DOUBLE PRECISION, parameter :: var1 = 8
      DOUBLE PRECISION, parameter :: var2 = 0.5

    !Continuous definition
      ! NEW => dS = (maxDs-minDs) * exp(-10 * (sqrt(x**2)-1)**2) + minDs

      dS = (maxDs-minDs) * exp(-10 * (sqrt(x**2)-1d0)**2) + minDs

      RETURN
      END

      SUBROUTINE dEfunc(X, dE)
      DOUBLE PRECISION X
      DOUBLE PRECISION dE
      call dSfunc(X, dE)
      RETURN
      END

      SUBROUTINE dIfunc(X, dI)
      DOUBLE PRECISION X
      DOUBLE PRECISION dI
      call dSfunc(X, dI)
      dI = dI/10.0
      RETURN
      END

      SUBROUTINE dRfunc(X, dR)
      DOUBLE PRECISION X
      DOUBLE PRECISION dR
      call dSfunc(X, dR)
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE DERIVF(T, X, U, UX, UXX, DFDU, DFDUX, DFDUXX, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE IS USED TO DEFINE THE INFORMATION ABOUT THE 
C       PDE REQUIRED TO FORM THE ANALYTIC JACOBIAN MATRIX FOR THE DAE
C       OR ODE SYSTEM. ASSUMING THE PDE IS OF THE FORM
C                        UT = F(T, X, U, UX, UXX)
C       THIS ROUTINE RETURNS THE JACOBIANS D(F)/D(U), D(F)/D(UX), AND
C       D(F)/D(UXX).
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
      INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
      DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
      DOUBLE PRECISION        UXX(NPDE)
C                               UXX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SECOND SPATIAL DERIVATIVE OF THE 
C                               SOLUTION AT THE POINT (T,X).
C
C OUTPUT:
      DOUBLE PRECISION        DFDU(NPDE,NPDE)
C                               DFDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DFDUX(NPDE,NPDE)
C                               DFDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DFDUXX(NPDE,NPDE)
C                               DFDUXX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SECOND SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.

C-----------------------------------------------------------------------
C
C     ASSIGN DFDU(1:NPDE,1:NPDE), DFDUX(1:NPDE,1:NPDE), AND
C     DFDUXX(1:NPDE,1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE 
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C

    !Stub function to allow compilation
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE BNDXA(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C       LEFT SPATIAL END POINT X = XA.
C                           B(T, U, UX) = 0
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
      INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,XA).
C
      DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,XA).
C
C OUTPUT:
      DOUBLE PRECISION        BVAL(NPDE)
C                               BVAL(1:NPDE) IS THE BOUNDARY CONTIDITION
C                               AT THE LEFT BOUNDARY POINT.


      BVAL(1)= UX(1) !- #FILL# - Value at the left boundary
      BVAL(2) = UX(2)
      BVAL(3) = UX(3)
      BVAL(4) = UX(4)
C
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE BNDXB(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C       RIGHT SPATIAL END POINT X = XB.
C                           B(T, U, UX) = 0
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
      INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,XB).
C
      DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,XB).
C
C OUTPUT:
      DOUBLE PRECISION        BVAL(NPDE)
C                               BVAL(1:NPDE) IS THE BOUNDARY CONTIDITION
C                               AT THE RIGHT BOUNDARY POINT.



      BVAL(1)= UX(1) !- #FILL# - Value at the left boundary
      BVAL(2) = UX(2)
      BVAL(3) = UX(3)
      BVAL(4) = UX(4)

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE DIFBXA(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY 
C       CONDITIONS AT THE LEFT SPATIAL END POINT X = XA. FOR THE 
C       BOUNDARY CONDITION EQUATION
C                              B(T, U, UX) = 0
C       THE PARTIAL DERIVATIVES DB/DU, DB/DUX, AND DB/DT ARE SUPPLIED 
C       BY THIS ROUTINE.
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
      INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C                                 U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
      DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
C OUTPUT:
      DOUBLE PRECISION        DBDU(NPDE,NPDE)
C                               DBDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DBDUX(NPDE,NPDE)
C                               DBDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DBDT(NPDE)
C                               DBDT(I) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO TIME T.

C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF 
C     U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C

    !Stub function to allow compilation

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE DIFBXB(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY 
C       CONDITIONS AT THE RIGHT SPATIAL END POINT 1 = XB. FOR THE 
C       BOUNDARY CONDITION EQUATION
C                              B(T, U, UX) = 0
C       THE PARTIAL DERIVATIVES DB/DU, DB/DUX, AND DB/DT ARE SUPPLIED 
C       BY THIS ROUTINE.
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
      INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
      DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
C OUTPUT:
      DOUBLE PRECISION        DBDU(NPDE,NPDE)
C                               DBDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DBDUX(NPDE,NPDE)
C                               DBDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE 
C                               UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DBDT(NPDE)
C                               DBDT(I) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO TIME T.

C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF 
C     U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
    !Stub function to allow compilation
C
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE UINIT(X, U, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE IS USED TO RETURN THE NPDE-VECTOR OF INITIAL 
C       CONDITIONS OF THE UNKNOWN FUNCTION AT THE INITIAL TIME T = T0 
C       AT THE SPATIAL COORDINATE X.
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C  INPUT:
      DOUBLE PRECISION        X
C                               THE SPATIAL COORDINATE.
C
      INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
C OUTPUT:
      DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS VECTOR OF INITIAL VALUES OF
C                               THE UNKNOWN FUNCTION AT T = T0 AND THE
C                               GIVEN VALUE OF X.

C-----------------------------------------------------------------------
C
C     ASSIGN U(1:NPDE) THE INITIAL VALUES OF U(T0,X).
C
! C-----------------------------------------------------------------------
!       DOUBLE PRECISION PI

!       PI = ACOS(-1.d0)
! C-----------------------------------------------------------------------

    !Discontinuous definition
    !IF (X .LT. 0.d0) THEN
        !U(1) = 0.96 * exp(-10 * ((X+1.d0) ** 2))
    !    U(3) = 0.2 * exp(-10 * ((X+1.d0) ** 2))
    !    U(1) = 1.0 - U(3)
    !ELSE 
        !U(1) = 0.96 * exp(-10 * ((X-1.d0) ** 2))
    !    U(3) = 0.0!0.2 * exp(-10 * ((X-1.d0) ** 2))
    !    U(1) = 1.0 - U(3)
    !END IF
    !U(1) = 0.96 * exp(-10 * (X**2))
    !U(2) = 0.0d0
    !U(3) = 0.04 * exp(-10 * (X**2))
    !U(4) = 0.0d0
C

    !Continuous defintion
      ! OLD DEFINITION OF THE initial conditions
      ! U(3) = 0.2 * exp(-10 * ((X+1.d0) ** 2))
      ! U(1) = N - U(3)
      ! U(2) = 0.0d0
      ! U(4) = 0.0d0
      
      DOUBLE PRECISION N
      ! NEW changed N from 10d0
      N = 37d6 ! make N much greater than the actual values S, E, I, R

      ! NEW changed from - U(3) = 2 * exp(-10 * ((X+1.d0) ** 2))
      ! this seemed good => U(3) = 100 * exp(-10 * ((X+1.d0) ** 2))
      ! U(3) = 10 * exp(-10 * ((X+1.d0) ** 2))
      U(3) = 100 * exp(-(X ** 2))
      U(1) = N - U(3)
      U(2) = 0.0d0
      U(4) = 0.0d0
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE TRUU(T, X, U, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C     THIS FUNCTION PROVIDES THE EXACT SOLUTION OF THE PDE.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
      INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
C OUTPUT:
      DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE EXACT SOLUTION AT THE 
C                               POINT (T,X).
C-----------------------------------------------------------------------
      DOUBLE PRECISION PI

      PI = ACOS(-1.d0)
C-----------------------------------------------------------------------

      U(1) = 0! #FILL# - I don't think this does anything
C
      RETURN
      END

c-----------------------------------------------------------------------
      subroutine header(nout)
c-----------------------------------------------------------------------
c purpose:
c       this subroutine writes a header describing the npde dimensional
c       parabolic partial differential equation
c                        ut = f(t, x, u, ux, uxx).
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
      integer                 nout
c                               nout is the output unit number.
c-----------------------------------------------------------------------
c constants:
      double precision        t0
      parameter              (t0 = 0.0d0)
c
      double precision        xa
      parameter              (xa = 0.0d0)
c
      double precision        xb
      parameter              (xb = 1.0d0)
c-----------------------------------------------------------------------
c
c	write(nout,95) '#FILL#Name:'
c	write(nout,95) 'pde: '
c	write(nout,95) '   u_t =  #FILL#Equation,  '
c	write(nout,95) 'domain:'
c	write(nout,96) '   t0 =', t0, ' < t,'
c	write(nout,96) '   xa =', xa, ' <= x <= xb =', xb, ','
c
c	return
c	95 format(a)
c	96 format(a,e13.5,a,e13.5,a,e13.5,a,e13.5,a)
      end



