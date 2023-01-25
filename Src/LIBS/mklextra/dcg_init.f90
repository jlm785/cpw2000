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

subroutine dcg_init ( n, x, b, job, ipar, dpar, tmp )

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  n                               !<  matrix dimension
  real(REAL64), intent(in)           ::  b(n)                            !<  right hand-side

! input and output

  real(REAL64), intent(inout)        ::  x(n)                            !<  guess of solution

  integer, intent(inout)             ::  job                             !<  task to be performed on output, 1 on input on first call

  real(REAL64), intent(inout)        ::  tmp(n,4)                        !<  tmp(:,1)=p, input for A p;  tmp(:,2)=q=A p;   tmp(:,3)=r, error estimate;   tmp(:,4)=z=Pr.

  integer, intent(inout)             ::  ipar(128)                       !<  integer parameters.   ipar(12:128) can be used internally
  real(REAL64), intent(inout)        ::  dpar(128)                       !<  real parameters.  dpar(9:128) can be used internally

! "saved" variables,

  integer               ::  iter                                         !  iteration number
  integer               ::  ilbl                                         !  state of algorithm

  iter = 0
  ilbl = 1

  ipar(20) = iter
  ipar(21) = ilbl

! get rid of "unused" messages while debugging

  job = 0
  dpar(31) = b(1)
  tmp(1, 1 ) = x(1)

! for compatibility with MKL not used here...

  ipar( 1) = n
  ipar( 2) = 6
  ipar( 3) = 1
  ipar( 4) = 0
  ipar( 5) = min(150,n)
  ipar( 6) = 1
  ipar( 7) = 1
  ipar( 8) = 1
  ipar( 9) = 0
  ipar(10) = 1
  ipar(11) = 1

  dpar(1) = 1.0E-6_REAL64
  dpar(2) = 0.0_REAL64
  dpar(3) = 0.0_REAL64
  dpar(4) = 0.0_REAL64
  dpar(5) = 0.0_REAL64
  dpar(6) = 0.0_REAL64
  dpar(7) = 0.0_REAL64

  return

end subroutine dcg_init

