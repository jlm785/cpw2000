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

!>  Interpolates from one grid to an uniform grid using 
!>  Neville's algorithm for Lagrange interpolation.
!>  Based on the Numerical Recipes algorithm (recurrence on differences).

  subroutine grid_interp(xin,fin,nin,xgmin,xgmax,ng,fg,nordp1,dymax)

! Written 23 April 2018 by J. L. Martins
! copyright  Jose Luis Martins/INESC-MN

! version 4.98

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          ::  iowrite = 6                       !  default tape for writing

! input

  integer, intent(in)               ::  nin                              !<  number of points in the input (old) grid
  real(REAL64), intent(in)          ::  xin(nin)                         !<  points on the input grid, x(i) > x(i-1)
  real(REAL64), intent(in)          ::  fin(nin)                         !<  f(x(i))

  real(REAL64), intent(in)          ::  xgmin, xgmax                     !<  first and last points on the uniform grid
  integer, intent(in)               ::  ng                               !<  number of points in the uniform grid

  integer, intent(in)               ::  nordp1                           !<  norder of Lagrange interpolation + 1

! output

  real(REAL64), intent(out)         ::  fg(ng)                           !<  value of the function on the uniform grid 
  real(REAL64), intent(out)         ::  dymax                            !<  estimate of the accuracy (from nordp1-1 interpolation)

! allocatable arrays

  real(REAL64), allocatable         ::  xg(:)                            !  uniform grid

  real(REAL64), allocatable         ::  cmj(:), dmj(:)                   !  difference arrays of Neville's algorithm

! local variables

  real(REAL64)   ::  xdif
  integer        ::  near, jmin
  real(REAL64)   ::  dif, difmin
  integer        ::  j0, ns
  real(REAL64)   ::  xnum, dm, dp, dy

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-15_REAL64

 
! counters

  integer     ::  i, j, m


! check grid sizes and if original grid is ordered

  if(nin < 3) then
    write(iowrite,*) '  grid_interp,  nin = ', nin

    stop

  endif

  if(ng < 2) then
    write(iowrite,*) '  grid_interp,  ng = ', ng

    stop

  endif

  if(nordp1 < 2) then
    write(iowrite,*) '  grid_interp,  nordp1 = ', nordp1

    stop

  endif

  if(nordp1 > nin) then
    write(iowrite,*) '  grid_interp,  nordp1, nin = ', nordp1, nin

    stop

  endif

  if(nordp1 > 10) then
    write(iowrite,*) '  grid_interp,  nordp1 = ', nordp1
    write(iowrite,*) '  Do you know what you are doing? '
  endif

  xdif = EPS*abs(xin(nin)-xin(1))
  do i = 1,nin-1
    if(xin(i+1) < xin(i) + xdif) then
      write(iowrite,*) '  grid_interp,  i,xin(i+1),xin(i) = ', i,xin(i+1),xin(i)

      stop

    endif
  enddo

  if(xgmax < xgmin) then
      write(iowrite,*) '  grid_interp,  xgmin,xgmax = ', xgmin,xgmax

      stop

  endif

  if(xgmin < xin(1)-xdif .or. xgmax > xin(nin)+xdif) then
    write(iowrite,*) '  grid_interp,  extrapolation not allowed'
    write(iowrite,*) '  xin(1),xin(nin),xgmin,xgmax = ', xin(1),xin(nin),xgmin,xgmax

    stop

  endif

! allocate arrays

  allocate(xg(ng))
  allocate(cmj(nordp1),dmj(nordp1))

  xg(1) = xgmin
  xdif = (xgmax - xgmin)/(ng-1)
  do i = 2,ng
    xg(i) = xgmin + (i-1)*xdif
  enddo

  near = 1
  dymax = ZERO
  do i = 1,ng

    difmin = abs(xg(i)-xin(near))
    jmin = near

    do j = near,nin
      dif = xin(j) - xg(i)
      jmin = j
          
      if(dif > ZERO) exit

    enddo

    near = jmin

!   effective 0 for Lagrange interpolation

    j0 = max(0,near -1 - nordp1/2)
    j0 = min(j0,nin-nordp1)

    do j = 1,nordp1
      cmj(j) = fin(j0+j)
      dmj(j) = fin(j0+j)
    enddo

    ns = near-j0-1
    fg(i) = fin(j0+ns)
    ns = ns - 1

    do m = 1,nordp1-1

      do j = 1,nordp1-m
        dm = xin(j0+j) - xg(i)
        dp = xin(j0+j+m) - xg(i)
        xnum = (dmj(j) - cmj(j+1)) / (dp - dm)
        dmj(j) = dp*xnum
        cmj(j) = dm*xnum
      enddo

      if (2*ns < nordp1 - m)then
        dy = cmj(ns+1)
      else
        dy = dmj(ns)
        ns = ns - 1
      endif

      fg(i) = fg(i) + dy

    enddo

    dymax = max(dymax,abs(dy))

  enddo

  deallocate(xg)
  deallocate(cmj,dmj)

  return
  end subroutine grid_interp
