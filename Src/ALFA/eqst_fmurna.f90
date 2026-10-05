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

!>  Murnaghan equation of state.
!>  Least squares fit of y = a + b x + c x**d for fixed d.
!>  The value of the function is zero for the optimal d.
!>
!>  \author       David M. Wood, Jose Luis Martins
!>  \version      5.13
!>  \date         early eighties, 3 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_fmurna(emurna, d, npt, xscl, yscl, a, b, c, mxdnpt)

! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! Data passed as arguments instead of common blocks. 3 October 2026. JLM+claude

  implicit none
  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)         ::  mxdnpt                                 !<  array dimension of number of points
  real(REAL64), intent(in)    ::  d                                      !<  exponent of the fit
  integer, intent(in)         ::  npt                                    !<  number of calculated energies
  real(REAL64), intent(in)    ::  xscl(mxdnpt)                           !<  scaled volume
  real(REAL64), intent(in)    ::  yscl(mxdnpt)                           !<  scaled energy

! output

  real(REAL64), intent(out)   ::  emurna                                 !<  error, should be zero for optimal d

  real(REAL64), intent(out)   ::  a                                      !<  fit parameter, y = a + b x + c x**d
  real(REAL64), intent(out)   ::  b                                      !<  fit parameter, y = a + b x + c x**d
  real(REAL64), intent(out)   ::  c                                      !<  fit parameter, y = a + b x + c x**d

! local

  real(REAL64)              ::  sx, sy,sxd, sxy, sx2, sx1pd
  real(REAL64)              ::  sxdy, sx2d, sxdp1l
  real(REAL64)              ::  sxdlxy, sxdlx, sx2dlx
  real(REAL64)              ::  x, alx, xd, y
  real(REAL64)              ::  alpha, beta, gama
  real(REAL64)              ::  u, w
  real(REAL64)              ::  fn                                       !  real(npt)

! parameters

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer    ::  i


  fn = UM*npt

  sx = ZERO
  sy = ZERO
  sxd = ZERO
  sxy = ZERO
  sx2 = ZERO
  sx1pd = ZERO
  sxdy = ZERO
  sx2d = ZERO
  sxdlxy = ZERO
  sxdlx = ZERO
  sx2dlx = ZERO
  sxdp1l = ZERO
  do i = 1,npt
    x = xscl(i)
    alx = log(x)
    xd = x**d
    y = yscl(i)
    sx = sx+x
    sy = sy+y
    sxd = sxd+xd
    sxy = sxy+x*y
    sx2 = sx2+x*x
    sx1pd = sx1pd+x*xd
    sxdy = sxdy+y*xd
    sx2d = sx2d+xd*xd
    sxdlxy = sxdlxy+alx*y*xd
    sxdlx = sxdlx+alx*xd
    sx2dlx = sx2dlx+xd*xd*alx
    sxdp1l = sxdp1l+x*alx*xd
  enddo
  alpha = -sxd/sx2d
  beta = -sx1pd/sx2d
  gama = sxdy/sx2d
  u = -(sx+alpha*sx1pd)/(sx2+beta*sx1pd)
  w = (sxy-gama*sx1pd)/(sx2+beta*sx1pd)
  a = -(sxd*(gama+beta*w)+sx*w-sy)/(fn+u*sx+(alpha+beta*u)*sxd)
  b = u*a+w
  c = (alpha+beta*u)*a+beta*w+gama

  emurna = sxdlxy-a*sxdlx-b*sxdp1l-c*sx2dlx

  return

end subroutine eqst_fmurna
