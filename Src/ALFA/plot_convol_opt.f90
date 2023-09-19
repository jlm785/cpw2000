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

!>  performs a convolution with a square function
!>  The width of the square function optimized to
!>  be flat at the center of the material
!>
!>  \author       Jose Luis Martins
!>  \version      5.07
!>  \date         11 September 2023.
!>  \copyright    GNU Public License v2

subroutine plot_convol_opt(n3, xave, ave, dave, yave, lfound, lopt,      &
           nrepeat, rbottom, rtop)

! Written 11 September 2023. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)         ::  n3                                     !<  number of points
  real(REAL64), intent(in)    ::  xave                                   !<  initial width of square function
  real(REAL64), intent(in)    ::  ave(n3)                                !<  original function

  logical, intent(in)         ::  lopt                                   !<  if true the width is optimized
  integer, intent(in)         ::  nrepeat                                !<  number of repeat units for the material
  real(REAL64), intent(in)    ::  rbottom                                !<  bottom lattice coordinate of the material
  real(REAL64), intent(in)    ::  rtop                                   !<  top lattice coordinate of the material

! output

  real(REAL64), intent(out)   ::  dave(n3)                               !<  function convoluted with square function
  real(REAL64), intent(out)   ::  yave                                   !<  averaged function in the center of the material.
  logical, intent(out)        ::  lfound                                 !<  if true the optimized width is within the search interval, if not optimized is true.


! local variables

  integer                ::  jbottom, jtop
  integer                ::  jj, jmin, jmax

  real(REAL64)           ::  fac
  real(REAL64)           ::  xsum, xsumsq
  real(REAL64)           ::  sigma

  real(REAL64)           ::  xtry, xopt

  integer                ::  kmax
  real(REAL64)           ::  xrange

! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer         ::  j, k, kk


! hard coded values.  Change if you know what you are doing

  kmax = 201
  xrange = 0.6

  call plot_convol(n3, xave, ave, dave)

  lfound = .TRUE.

! calculates relevant interval in the center of the material
! discard 1/4 on each side up to a maximum of 2 repeat units.

  jbottom = nint(rbottom*n3)
  jtop = nint(rtop*n3)
  if(jtop <= jbottom) jtop = jtop+n3

  fac = min(UM/4, (2*UM)/nrepeat)

  jmin = jbottom + nint((jtop-jbottom)*fac)
  jmax = jtop - nint((jtop-jbottom)*fac)

  if(lopt) then

    xopt = xave

    xsum = ZERO
    do j = jmin,jmax
      jj = mod(j-1,n3)+1
      xsum = xsum + dave(jj)
    enddo
    xsum = xsum / (jmax-jmin+1)

    xsumsq = ZERO
    do j = jmin,jmax
      jj = mod(j-1,n3)+1
      xsumsq = xsumsq + (dave(jj)-xsum)*(dave(jj)-xsum)
    enddo
    xsumsq = xsumsq / (jmax-jmin+1)
    sigma = sqrt(xsumsq)

!   Stupid but safe minimization.  Function is V-shaped and fast to calculate
!   when compared to initial calculation of ave(j)

    kk = 1
    sigma = sigma + UM
    do k = 1,kmax
      xtry = (UM-xrange/2)*xave + (k-1)*(xrange*xave/(kmax-1))

      call plot_convol(n3, xtry, ave, dave)

      xsum = ZERO
      do j = jmin,jmax
        jj = mod(j-1,n3)+1
        xsum = xsum + dave(jj)
      enddo
      xsum = xsum / (jmax-jmin+1)

      xsumsq = ZERO
      do j = jmin,jmax
        jj = mod(j-1,n3)+1
        xsumsq = xsumsq + (dave(jj)-xsum)*(dave(jj)-xsum)
      enddo
      xsumsq = xsumsq / (jmax-jmin+1)

      if(sqrt(xsumsq) < sigma) then
        xopt = xtry
        sigma = sqrt(xsumsq)
        kk = k
      endif

    enddo

!   Warn when minimum is expected outside search boundary

    if(kk == 1 .or. kk == kmax) lfound = .FALSE.

    call plot_convol(n3, xopt, ave, dave)

  endif

  yave = ZERO
  do j = jmin,jmax
    jj = mod(j-1,n3)+1
    yave = yave + dave(jj)
  enddo
  yave = yave / (jmax-jmin+1)

  return

end subroutine plot_convol_opt

