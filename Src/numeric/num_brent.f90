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

!>  Finds the minimum of a function in the interval ax, bx.
!>
!>  Brent method: combination of golden section search and
!>  successive parabolic interpolation.
!>
!>  Translated from the function fmin of
!>    G. E. Forsythe, M. A. Malcolm and C. B. Moler,
!>    Computer Methods for Mathematical Computations,
!>    (Prentice-Hall, 1977), section 8.2, (netlib.org/fmm/fmin.f)
!>  which is a translation of the algol 60 procedure localmin of
!>    R. P. Brent, Algorithms for Minimization without Derivatives,
!>    (Prentice-Hall, 1973), chapter 5.
!>
!>  The minimum is found with an error less than 3*sqrt(macheps)*abs(x) + tol.
!>  If f is not unimodal in the interval it may find a local minimum.
!>
!>  This is a reverse communication interface:
!>    call it first with iflag = 0.
!>    On return, if iflag = 1, calculate fval = f(xval) and call it again
!>    without changing anything else.
!>    On return, if iflag = 0, xval is the minimum.
!>
!>    iflag = 0
!>    do
!>      call num_brent(xval, fval, ax, bx, tol, iflag)
!>      if(iflag == 0) exit
!>      fval = f(xval)
!>    enddo
!>
!>  \author       R. P. Brent + claude + JLM
!>  \version      5.13
!>  \date         3 October 2026.
!>  \copyright    GNU Public License v2

subroutine num_brent(xval, fval, ax, bx, tol, iflag)

! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! Rewritten from fmin with a reverse communication interface. 3 October 2026. JLM+claude
! Renamed num_brent, pure numerical routine. 3 October 2026. JLM+claude

  implicit none
  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)    ::  ax                                     !<  left end of the interval
  real(REAL64), intent(in)    ::  bx                                     !<  right end of the interval
  real(REAL64), intent(in)    ::  tol                                    !<  absolute tolerance for the minimum
  real(REAL64), intent(in)    ::  fval                                   !<  f(xval), when iflag = 1 on input

! input and output

  integer, intent(inout)      ::  iflag                                  !<  0 start (input) or finished (output), 1 fval requested (output) or computed (input)

! output

  real(REAL64), intent(out)   ::  xval                                   !<  point where f is requested, or the minimum when finished

! local, saved between calls

  real(REAL64), save   ::  a, b                                          !  the minimum is in the interval a, b
  real(REAL64), save   ::  x, fx                                         !  point with the lowest value of f
  real(REAL64), save   ::  w, fw                                         !  point with the second lowest value of f
  real(REAL64), save   ::  v, fv                                         !  previous value of w
  real(REAL64), save   ::  u                                             !  point where f was requested
  real(REAL64), save   ::  d, e                                          !  last step and the step before
  integer, save        ::  iter
  integer, save        ::  istate                                        !  which evaluation of f was requested

! local

  real(REAL64)      ::  p, q, r
  real(REAL64)      ::  fu
  real(REAL64)      ::  xm
  real(REAL64)      ::  tol1, tol2
  logical           ::  lgolden

! parameters

  integer, parameter        ::  ITMAX = 200                              !  safety, not in the original
  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter   ::  CGOLD = (3 - sqrt(5*UM)) / 2             !  squared inverse of the golden ratio
  real(REAL64), parameter   ::  EPS = sqrt(epsilon(UM))                  !  square root of the relative machine precision


! start, request f(x) at the golden section point

  if(iflag == 0) then
    a = ax
    b = bx
    v = a + CGOLD*(b - a)
    w = v
    x = v
    d = ZERO
    e = ZERO
    istate = 1
    xval = x
    iflag = 1

    return

  endif

! processes the requested value of f

  if(istate == 1) then

    fx = fval
    fv = fx
    fw = fx
    iter = 0

  else

!   update a, b, v, w, and x

    fu = fval

    if(fu <= fx) then
      if(u >= x) then
        a = x
      else
        b = x
      endif
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
    else
      if(u < x) then
        a = u
      else
        b = u
      endif
      if(fu <= fw .or. w == x) then
        v = w
        fv = fw
        w = u
        fw = fu
      elseif(fu <= fv .or. v == x .or. v == w) then
        v = u
        fv = fu
      endif
    endif

  endif

  iter = iter + 1
  if(iter > ITMAX) stop 'num_brent exceeded maximum iterations'

! check stopping criterion

  xm = (a + b)/2
  tol1 = EPS*abs(x) + tol/3
  tol2 = 2*tol1

  if(abs(x - xm) <= tol2 - (b - a)/2) then
    xval = x
    iflag = 0

    return

  endif

! is golden-section necessary

  lgolden = .TRUE.

  if(abs(e) > tol1) then

!   fit parabola

    r = (x - w)*(fx - fv)
    q = (x - v)*(fx - fw)
    p = (x - v)*q - (x - w)*r
    q = 2*(q - r)
    if(q > ZERO) p = -p
    q = abs(q)
    r = e
    e = d

!   is parabola acceptable

    if(abs(p) < abs(q*r/2) .and. p > q*(a - x) .and. p < q*(b - x)) then

!     a parabolic interpolation step

      d = p/q
      u = x + d

!     f must not be evaluated too close to ax or bx

      if(u - a < tol2 .or. b - u < tol2) d = sign(tol1, xm - x)

      lgolden = .FALSE.

    endif

  endif

! a golden-section step

  if(lgolden) then
    if(x >= xm) then
      e = a - x
    else
      e = b - x
    endif
    d = CGOLD*e
  endif

! f must not be evaluated too close to x

  if(abs(d) >= tol1) then
    u = x + d
  else
    u = x + sign(tol1, d)
  endif

! request f(u)

  istate = 2
  xval = u

  return

end subroutine num_brent
