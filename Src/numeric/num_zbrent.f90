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

!>  Finds a zero of a function in the interval ax, bx, where f(ax) and
!>  f(bx) have opposite signs.
!>
!>  Brent method: combination of bisection, linear (secant) and inverse
!>  quadratic interpolation.
!>
!>  Translated from the function zeroin of
!>    G. E. Forsythe, M. A. Malcolm and C. B. Moler,
!>    Computer Methods for Mathematical Computations,
!>    (Prentice-Hall, 1977), section 7.2, (netlib.org/fmm/zeroin.f)
!>  which is a translation of the algol 60 procedure zero of
!>    R. P. Brent, Algorithms for Minimization without Derivatives,
!>    (Prentice-Hall, 1973), chapter 4.
!>    R. P. Brent, Computer Journal 14, 422 (1971).
!>
!>  The zero is found within a tolerance 4*macheps*abs(x) + tol.
!>
!>  This is a reverse communication interface:
!>    call it first with iflag = 0.
!>    On return, if iflag = 1, calculate fval = f(xval) and call it again
!>    without changing anything else.
!>    On return, if iflag = 0, xval is the zero.
!>
!>    iflag = 0
!>    do
!>      call num_zbrent(xval, fval, ax, bx, tol, iflag)
!>      if(iflag == 0) exit
!>      fval = f(xval)
!>    enddo
!>
!>  \author       R. P. Brent + claude + JLM
!>  \version      5.13
!>  \date         3 October 2026.
!>  \copyright    GNU Public License v2

subroutine num_zbrent(xval, fval, ax, bx, tol, iflag)

! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! Rewritten from zeroin with a reverse communication interface. 3 October 2026. JLM+claude
! Renamed num_zbrent, pure numerical routine. 3 October 2026. JLM+claude

  implicit none
  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)    ::  ax                                     !<  one end of the interval
  real(REAL64), intent(in)    ::  bx                                     !<  other end of the interval
  real(REAL64), intent(in)    ::  tol                                    !<  absolute tolerance for the zero
  real(REAL64), intent(in)    ::  fval                                   !<  f(xval), when iflag = 1 on input

! input and output

  integer, intent(inout)      ::  iflag                                  !<  0 start (input) or finished (output), 1 fval requested (output) or computed (input)

! output

  real(REAL64), intent(out)   ::  xval                                   !<  point where f is requested, or the zero when finished

! local, saved between calls

  real(REAL64), save   ::  a, b, c                                       !  b best estimate, zero between b and c, a previous b
  real(REAL64), save   ::  fa, fb, fc
  real(REAL64), save   ::  d, e                                          !  last step and the step before
  integer, save        ::  iter
  integer, save        ::  istate                                        !  which evaluation of f was requested

! local

  real(REAL64)      ::  p, q, r, s
  real(REAL64)      ::  xm
  real(REAL64)      ::  tol1
  logical           ::  lnewc                                            !  the zero is between a and b
  logical           ::  lbisect

! parameters

  integer, parameter        ::  ITMAX = 200                              !  safety, not in the original
  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter   ::  EPS = epsilon(UM)                        !  relative machine precision


! start, request f(ax)

  if(iflag == 0) then
    a = ax
    b = bx
    istate = 1
    xval = a
    iflag = 1

    return

  endif

! processes the requested value of f

  if(istate == 1) then

!   request f(bx)

    fa = fval
    istate = 2
    xval = b

    return

  elseif(istate == 2) then

    fb = fval

    if(fa*fb > ZERO) stop 'root must be bracketed for num_zbrent'

    iter = 0
    lnewc = .TRUE.

  else

    fb = fval
    lnewc = (fb*(fc/abs(fc)) > ZERO)

  endif

  iter = iter + 1
  if(iter > ITMAX) stop 'num_zbrent exceeded maximum iterations'

! begin step

  if(lnewc) then
    c = a
    fc = fa
    d = b - a
    e = d
  endif

  if(abs(fc) < abs(fb)) then
    a = b
    b = c
    c = a
    fa = fb
    fb = fc
    fc = fa
  endif

! convergence test

  tol1 = 2*EPS*abs(b) + tol/2
  xm = (c - b)/2

  if(abs(xm) <= tol1 .or. fb == ZERO) then
    xval = b
    iflag = 0

    return

  endif

! is bisection necessary

  lbisect = .TRUE.

  if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then

    if(a == c) then

!     linear interpolation

      s = fb/fa
      p = 2*xm*s
      q = UM - s

    else

!     inverse quadratic interpolation

      q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2*xm*q*(q - r) - (b - a)*(r - UM))
      q = (q - UM)*(r - UM)*(s - UM)

    endif

!   adjust signs

    if(p > ZERO) q = -q
    p = abs(p)

!   is interpolation acceptable

    if(2*p < 3*xm*q - abs(tol1*q) .and. p < abs(e*q/2)) then
      e = d
      d = p/q
      lbisect = .FALSE.
    endif

  endif

  if(lbisect) then
    d = xm
    e = d
  endif

! complete step

  a = b
  fa = fb
  if(abs(d) > tol1) then
    b = b + d
  else
    b = b + sign(tol1, xm)
  endif

! request f(b)

  istate = 3
  xval = b

  return

end subroutine num_zbrent
