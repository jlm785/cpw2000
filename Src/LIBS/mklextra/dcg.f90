!------------------------------------------------------------!
! This file is distributed as part of the cpw2000 code and   !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the cpw2000        !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the cpw2000 code is not yet written         !
!                                                            !
! The cpw2000 code is hosted on GitHub:                      !
!                                                            !
! https://github.com/jlm785/cpw2000                          !
!------------------------------------------------------------!

!>  RCI (reverse communication interface) to solve a linear problem A x = b,
!>  for a symmetric matrix A.  If A is not positive definite it may not work...
!>
!>  The cg_rc subroutine of of John Burkardt was modified to make it
!>  compatible with the dcg subroutine of intel MKL library,
!>  and hopefuly thread safe.  However it does not check for convergence,
!>  so use it with care.
!>
!>  \author       John Burkardt, Jose Luis Martins
!>  \version      5.06
!>  \date         12 January 2013, 16 January 2023.
!>  \copyright    GNU Public License v2

subroutine dcg ( n, x, b, job, ipar, dpar, tmp )

!*****************************************************************************80
!
!! cg_rc() is a reverse communication conjugate gradient routine.
!
!  Discussion:
!
!    This routine seeks a solution of the linear system A*x=b
!    where b is a given right hand side vector, A is an n by n
!    symmetric positive definite matrix, and x is an unknown vector
!    to be determined.
!
!    Under the assumptions that the matrix A is large and sparse,
!    the conjugate gradient method may provide a solution when
!    a direct approach would be impractical because of excessive
!    requirements of storage or even of time.
!
!    The conjugate gradient method presented here does not require the
!    user to store the matrix A in a particular way.  Instead, it only
!    supposes that the user has a way of calculating
!      y = alpha * A * x + b * y
!    and of solving the preconditioned linear system
!      M * x = b
!    where M is some preconditioning matrix, which might be merely
!    the identity matrix, or a diagonal matrix containing the
!    diagonal entries of A.
!
!    This routine was extracted from the "templates" package.
!    There, it was not intended for direct access by a user;
!    instead, a higher routine called "cg()" was called once by
!    the user.  The cg() routine then made repeated calls to
!    cgrevcom() before returning the result to the user.
!
!    The reverse communication feature of cgrevcom() makes it, by itself,
!    a very powerful function.  It allows the user to handle issues of
!    storage and implementation that would otherwise have to be
!    mediated in a fixed way by the function argument list.  Therefore,
!    this version of cgrecom() has been extracted from the templates
!    library and documented as a stand-alone procedure.
!
!    The user sets the value of JOB to 1 before the first call,
!    indicating the beginning of the computation, and to the value of
!    2 thereafter, indicating a continuation call.
!    The output value of JOB is set by cgrevcom(), which
!    will return with an output value of JOB that requests a particular
!    new action from the user.
!
!    THIS WAS CHANGED, INITIALIZATION IS IN DGC-INIT, PLUS FIRST CALL TO DGC
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994,
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!  Parameters:
!
!    Input, integer N, the dimension of the matrix.
!
!    Input, real(REAL64) B(N), the right hand side vector.
!
!    Input/output, real(REAL64) X(N).  On first call, the user
!    should store an initial guess for the solution in X.  On return with
!    JOB = 4, X contains the latest solution estimate.
!
!    Input/output, real(REAL64) R(N)->tmp_4, Z(N)->tmp_3, P(N)->tmp_1, Q(N)->tmp_2,
!    information used by the program during the calculation.  The user
!    does not need to initialize these vectors.  However, specific
!    return values of JOB may require the user to carry out some computation
!    using data in some of these vectors.
!
!    Input/output, integer JOB, communicates the task to be done.
!    The user needs to set the input value of JOB to 1, before the first call,
!    and then to 2 for every subsequent call for the given problem.
!    The output value of JOB indicates the requested user action.
!    * JOB = 1, compute Q = A * P;
!    * JOB = 2: solve M*Z=R, where M is the preconditioning matrix;
!    * JOB = 3: compute R = R - A * X;
!    * JOB = 4: check the residual R for convergence.
!               If satisfactory, terminate the iteration.
!               If too many iterations were taken, terminate the iteration.
!
!    THIS HAS BEEN CHANGED     WARNING
!
!    * JOB = 1, compute Q = A * P;
!    * JOB = 2: check the residual R for convergence.
!               If satisfactory, terminate the iteration.
!               If too many iterations were taken, terminate the iteration.
!    * JOB = 3: solve M*Z=R, where M is the preconditioning matrix;



! dpar(5) can be used to check for convergence. 26 January 2023. JLM


  implicit none


  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  n                               !<  matrix dimension
  real(REAL64), intent(in)           ::  b(n)                            !<  right hand-side

! input and output

  real(REAL64), intent(inout)        ::  x(n)                            !<  guess of solution

  integer, intent(inout)             ::  job                             !<  task to be performed on output, 1 on input on first call

  real(REAL64), intent(inout)        ::  tmp(n,4)                        !<  tmp(:,1)=p, input for A p;  tmp(:,2)=q=A p;   tmp(:,3)=z=Pr;   tmp(:,4)=r, error estimate.

  integer, intent(inout)             ::  ipar(128)                       !<  integer parameters.   ipar(12:128) can be used internally.  ipar(4) is the iteration number.
  real(REAL64), intent(inout)        ::  dpar(128)                       !<  real parameters.  dpar(9:128) can be used internally.  dpar(5) is the square norm of the error and can be used to check for convergence.

! "saved" variables,

  integer               ::  iter                                         !  iteration number
  integer               ::  ilbl                                         !  state of algorithm
  real(REAL64)          ::  rho                                          !  rho = r.z
  real(REAL64)          ::  rho_old                                      !  old value of rho

! local variables

  real(REAL64)    ::  alpha
  real(REAL64)    ::  beta
  real(REAL64)    ::  pdotq

!  Memory

   iter = ipar(20)
   ilbl = ipar(21)

   rho = dpar(20)
   rho_old = dpar(21)

! initialization

  if( ilbl == 1 ) then

  tmp(1:n, 3 ) = b(1:n)
  tmp(1:n, 1 ) = x(1:n)

  job = 1
  ilbl = 2

!  Ask the user to compute the initial residual.
  elseif ( ilbl == 2 ) then

    iter = 1

    tmp(1:n, 3 ) = tmp(1:n, 3 ) - tmp(1:n, 2 )

    job = 3
    ilbl = 3
!
!  Compute the direction.
!  Ask the user to compute ALPHA.
!  Save A*P to Q.
!
  else if ( ilbl == 3 ) then

    rho = dot_product ( tmp(1:n, 3 ), tmp(1:n, 4) )

    if ( 1 < iter ) then
      beta = rho / rho_old
      tmp(1:n, 4 ) = tmp(1:n, 4 ) + beta * tmp(1:n, 1 )
    end if

    tmp(1:n, 1 ) = tmp(1:n, 4 )

    job = 1
    ilbl = 4
!
!  Compute current solution vector.
!  Ask the user to check the stopping criterion.
!
  else if ( ilbl == 4 ) then

    pdotq = dot_product ( tmp(1:n, 1 ), tmp(1:n, 2 ) )
    alpha = rho / pdotq
    x(1:n) = x(1:n) + alpha * tmp(1:n, 1 )
    tmp(1:n, 3 ) = tmp(1:n, 3 ) - alpha * tmp(1:n, 2 )

    dpar(5) = dot_product ( tmp(1:n, 3 ), tmp(1:n, 3 ) )

    job = 2
    ilbl = 5
!
!  Begin the next step.
!  Ask for a preconditioner solve.
!
  else if ( ilbl == 5 ) then

    rho_old = rho
    iter = iter + 1

    job = 3
    ilbl = 3

  end if

  ipar(20) = iter
  ipar(21) = ilbl

  dpar(20) = rho
  dpar(21) = rho_old

! for compatibility with MKL

  ipar(4) = iter

  return

end subroutine dcg

