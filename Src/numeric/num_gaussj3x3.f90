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

!>  Solves the 3x3 linear system a x = b by Gauss-Jordan elimination
!>  with partial (row) pivoting.  The solution x is returned in b.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         3 October 2026.
!>  \copyright    GNU Public License v2

subroutine num_gaussj3x3(a, b)

! Written 3 October 2026, replaces gaussj (general n, full pivoting,
! also computed the inverse) in the eqst tool. JLM+claude

  implicit none
  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)        ::  a(3,3)                             !<  matrix of the linear system

! input and output

  real(REAL64), intent(inout)     ::  b(3)                               !<  on input right hand side, on output the solution

! local

  real(REAL64)         ::  aw(3,3)                                       !  working copy of a
  real(REAL64)         ::  big
  real(REAL64)         ::  pivinv
  real(REAL64)         ::  factor
  real(REAL64)         ::  tmp
  integer              ::  ipiv                                          !  pivot row

! constants

  real(REAL64), parameter    :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer              ::  i, j, k


  aw(:,:) = a(:,:)

  do k = 1,3

!   finds the pivot in column k

    ipiv = k
    big = abs(aw(k,k))
    do i = k+1,3
      if(abs(aw(i,k)) > big) then
        ipiv = i
        big = abs(aw(i,k))
      endif
    enddo

    if(big == ZERO) stop 'num_gaussj3x3: singular matrix'

!   interchanges rows k and ipiv

    if(ipiv /= k) then
      do j = 1,3
        tmp = aw(k,j)
        aw(k,j) = aw(ipiv,j)
        aw(ipiv,j) = tmp
      enddo
      tmp = b(k)
      b(k) = b(ipiv)
      b(ipiv) = tmp
    endif

!   normalizes row k

    pivinv = UM / aw(k,k)
    do j = k,3
      aw(k,j) = aw(k,j)*pivinv
    enddo
    b(k) = b(k)*pivinv

!   eliminates column k from the other rows

    do i = 1,3
      if(i /= k) then
        factor = aw(i,k)
        do j = k,3
          aw(i,j) = aw(i,j) - factor*aw(k,j)
        enddo
        b(i) = b(i) - factor*b(k)
      endif
    enddo

  enddo

  return

end subroutine num_gaussj3x3
