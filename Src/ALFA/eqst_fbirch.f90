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

!>  Birch equation of state.
!>  Least squares fit of the energy for fixed vzero.
!>  The value of the function is the sum of the squares of the residuals,
!>  minimal for the optimal vzero.
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         September 1989, 3 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_fbirch(qbirch, vzero, npt, xar, yscl, ezero, bzero, bprim, mxdnpt)

! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! f1 and f2 calculated here instead of calling eqst_f1 and eqst_f2. 3 October 2026. JLM+claude
! Data passed as arguments instead of common block. 3 October 2026. JLM+claude
! num_gaussj3x3 instead of eqst_gaussj. 3 October 2026. JLM+claude
! Transformed into subroutine, function had side effects. 5 October 2026. JLM

  implicit none
  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)         ::  mxdnpt                                 !<  array dimension of number of points
  real(REAL64), intent(in)    ::  vzero                                  !<  equilibrium volume
  integer, intent(in)         ::  npt                                    !<  number of calculated energies
  real(REAL64), intent(in)    ::  xar(mxdnpt)                            !<  volume
  real(REAL64), intent(in)    ::  yscl(mxdnpt)                           !<  scaled energy

! output

  real(REAL64), intent(out)   ::  qbirch                                 !<  quality of fit (least squares)

  real(REAL64), intent(out)   ::  ezero                                  !<  scaled energy at vzero
  real(REAL64), intent(out)   ::  bzero                                  !<  scaled bulk modulus
  real(REAL64), intent(out)   ::  bprim                                  !<  pressure derivative of bulk modulus

! other variables

  real(REAL64)                ::  a(3,3), b(3), c(3)
  real(REAL64)                ::  csi, sumx, yfit
  real(REAL64)                ::  x23                                    !  (vzero/v)**(2/3)
  real(REAL64)                ::  f1(npt), f2(npt)                       !  the two functions of the fit at each point

! constants

  real(REAL64), parameter    :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer                     ::  i, j, k

  do j = 1,3
    b(j) = ZERO
    do k = 1,3
      a(k,j) = ZERO
    enddo
  enddo

! the Birch fit is linear in the functions 1, f1 and f2 of (vzero/v)**(2/3)
!   f1 = (9/8) (1 - x23)^2  and  f2 = (3/4) (1 - x23)^3

  do i = 1,npt
    x23 = (vzero/xar(i))**((2*UM)/(3*UM))
    f1(i) = 9*(1 + x23*(-2+x23)) / (8*UM)
    f2(i) = 3*(1 + x23*(-3 + x23*(3-x23))) / (4*UM)
  enddo

  do i = 1,npt
    c(1) = UM
    c(2) = f1(i)
    c(3) = f2(i)
    do j = 1,3
      b(j) = b(j) + yscl(i)*c(j)
      do k = 1,3
        a(k,j) = a(k,j) + c(k)*c(j)
      enddo
    enddo
  enddo

  call num_gaussj3x3(a, b)

  ezero = b(1)
  bzero = b(2)/vzero
  csi = b(3)/b(2)
  bprim = 4*(UM-csi/(3*UM))

  sumx = ZERO
  do i = 1,npt
    yfit = b(1) + b(2)*f1(i) + b(3)*f2(i)
    sumx = sumx + (yscl(i)-yfit)*(yscl(i)-yfit)
  enddo
  qbirch = sumx

  return

end subroutine eqst_fbirch
