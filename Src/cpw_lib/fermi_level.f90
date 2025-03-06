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

!>  Finds the Fermi level and computes occupations numbers.
!>  If the temperature is small compared with the range
!>  of energies, it is considered as zero.
!>
!>  \author       Jose Luis Martins, Renata Wentzcovitch
!>  \version      5.11
!>  \date         June 16, 1987, March 5 2025.
!>  \copyright    GNU Public License v2

subroutine fermi_level(el, ztot, teleck,                                 &
      nrk, wgk, nband,                                                   &
      frac, efermi, eband, elects, bandwid, penngap,                     &
      mxdnrk, mxdbnd)

! Adapted from Sverre Froyen plane wave program.
! written june 16, 1987.jlm
! modified august 3, 1987.jlm
! modified 22 march 1999. jlm
! temperature added 25 nov 90. RMM
! modified for temperature 11 may 99.jlm
! Modified, documentation, January 2020. JLM
! Modified, indentation, ztot=0. March 5 2025. JLM


  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)                ::  mxdnrk                          !<  size of k-points
  integer, intent(in)                ::  mxdbnd                          !<  size of number of bands

  real(REAL64), intent(in)           ::  el(mxdbnd*mxdnrk)               !<  eigenvalues in Hartree for all the k-points,
  real(REAL64), intent(in)           ::  ztot                            !<  total valence charge density
  real(REAL64), intent(in)           ::  teleck                          !<  temperature (Kelvin) of electron system.
  integer, intent(in)                ::  nrk                             !<  number of k-points for integration in the irreducible wedge of the brillouin zone
  real(REAL64), intent(in)           ::  wgk(mxdnrk)                     !<  weight in the integration of k-point
  integer, intent(in)                ::  nband(mxdnrk)                   !<  number of bands for each k-points

! output

  real(REAL64), intent(out)          ::  frac(mxdbnd*mxdnrk)             !<  fractional ocupation of level j
  real(REAL64), intent(out)          ::  efermi                          !<  Fermi energy (or highest occupied state)
  real(REAL64), intent(out)          ::  eband                           !<  band energy (Hartree)
  real(REAL64), intent(out)          ::  elects                          !<  electronic temperature*entropy (hartree)
  real(REAL64), intent(out)          ::  bandwid                         !<  occupied band width estimate (Hartree)
  real(REAL64), intent(out)          ::  penngap                         !<  Penn gap estimate (Hartree)

! local allocatable arrays

  integer, allocatable               ::  jrk(:)
  integer, allocatable               ::  ind(:)

! local variables

  integer           ::  iel, ifrm, imx, jmin, jmax
  integer           ::  maxit                                            !  maximum number of iterations for root finding
  real(REAL64)      ::  sumt, sumw, sw2, fr, tempau
  real(REAL64)      ::  fl, fh, xl, xh, ff
  real(REAL64)      ::  f, df, dff, dx, dxold, efold, arg
  logical           ::  lexit


! counters

  integer   ::    i, j

! constants

  real(REAL64), parameter ::  ZERO = 0.0_REAL64 , UM = 1.0_REAL64
  real(REAL64), parameter ::  EPS = 0.000001_REAL64
  real(REAL64), parameter ::  SMALL = EPS*EPS
  real(REAL64), parameter ::  EV = 27.2116_REAL64
  real(REAL64), parameter ::  TAUTOK = 11604.9_REAL64 * EV


  allocate(jrk(mxdnrk*mxdbnd))
  allocate(ind(mxdnrk*mxdbnd))

  iel = 0
  do i=1,nrk
    do j=1,nband(i)
      iel = iel + 1
      jrk(iel) = i
    enddo
  enddo

! initialize array frac.

  do i=1,iel
    frac(i) = ZERO
  enddo

! sort the energy levels

  call sort(iel,el,ind)

! finds approximate level

  sumt = ztot/(2*UM)
  sumw = ZERO
  do i=1,iel
    ifrm = i
    sumw = sumw + wgk(jrk(ind(i)))

    if (sumw > sumt-SMALL) exit

  enddo

! finds degeneracies. ifrm identifies where the fermi level is

  imx = ifrm - 1
  if(ifrm > 1) then
    do i=1,imx
      jmin = ifrm - i

      if (abs(el(ind(jmin))-el(ind(ifrm))) > EPS) exit

    enddo
  else
    jmin = 0
  endif
  imx = ifrm + 1
  if(ifrm < iel) then
    do i=imx,iel
      jmax = i

      if (abs(el(ind(jmax))-el(ind(ifrm))) > EPS) exit

    enddo
    jmax = jmax - 1
  else
    jmax = ifrm
  endif

  tempau = teleck / TAUTOK

  if( tempau < EPS*abs(el(ind(iel))-el(ind(1))) ) then

!   zero temperature

!   the levels between jmin+1 and jmax are degenerate and at efermi

    sumw = ZERO
    do i=1,jmin
      sumw = sumw + wgk(jrk(ind(i)))
    enddo
    do i=1,jmin
      frac(ind(i)) = UM
    enddo
    imx = jmin + 1
    sw2 = ZERO
    do i=imx,jmax
      sw2 = sw2 + wgk(jrk(ind(i)))
    enddo
    if(abs(sw2) > SMALL) then
      fr = (ztot/(2*UM) - sumw)/sw2
      do i=imx,jmax
        frac(ind(i)) = fr
      enddo
    endif

!   fermi level

    efermi = el(ind(ifrm))

    elects = ZERO

  else

!   finite temperature

!   find out fermi level (rtsafe from numerical recipes)

    maxit = 300

!   bracket roots

    do i=1,jmin
      xl = el(ind(jmin - i + 1))
      fl = - ztot/(2*UM)
      do j=1,iel
        arg = (el(ind(j)) - xl) / tempau
        call fermi_dirac(arg,ff,dff)
        fl = fl + wgk(jrk(ind(j))) * ff
      enddo

      if (fl < ZERO) exit

    enddo

    do i = jmax,iel
      xh = el(ind(i))
      fh = - ztot/(2*UM)
      do j=1,iel
        arg = (el(ind(j)) - xh) / tempau
        call fermi_dirac(arg,ff,dff)
        fh = fh + wgk(jrk(ind(j))) * ff
      enddo

      if (fh > ZERO) exit

    enddo

!   finds root

    efermi = (xl + xh) / (2*UM)
    dxold = abs(xh - xl)
    dx = dxold

    f = - ztot/(2*UM)
    df = ZERO
    do j=1,iel
      arg = (el(ind(j)) - efermi) / tempau
      call fermi_dirac(arg,ff,dff)
      dff = -dff / tempau
      f = f + wgk(jrk(ind(j))) * ff
      df = df + wgk(jrk(ind(j))) * dff
    enddo

!   loop over a maximum number of iterations

    lexit = .FALSE.
    do i = 1,maxit
      if(((efermi - xh)*df - f)*((efermi - xl)*df - f) >= ZERO      &
              .or. abs(2*f) > abs(dxold*df)) then

!        bisection

         dxold = dx
         dx = (xh - xl) / (2*UM)
         efermi = xl + dx
         if (xl == efermi) lexit = .TRUE.

      else

!       newton

        dxold = dx
        dx = f / df
        efold = efermi
        efermi = efermi - dx
        if (efold == efermi) lexit = .TRUE.
      endif

      if (abs(dx) < SMALL) lexit = .TRUE.

      if(lexit) exit

      f = - ztot/(2*UM)
      df = ZERO
      do j=1,iel
        arg = (el(ind(j)) - efermi) / tempau
        call fermi_dirac(arg,ff,dff)
        dff = -dff / tempau
        f = f + wgk(jrk(ind(j))) * ff
        df = df + wgk(jrk(ind(j))) * dff
      enddo

      if (f < ZERO) then
          xl = efermi
      else
          xh = efermi
      endif
    enddo

    if(.not. lexit) then

      write(6,*)
      write(6,*) '   *** STOPPED in fermi_level  '
      write(6,*) ' unable to find it....'

      stop

    endif

!   fermi level was found

    do j=1,iel
      arg = (el(ind(j)) - efermi) / tempau
      call fermi_dirac(arg, ff, dff)
      frac(ind(j)) = ff
    enddo

    elects = ZERO

    iel = 0

    do i=1,nrk
      do j=1,nband(i)
        iel = iel + 1
        if (frac(iel) > SMALL) then
          elects = elects - 2*tempau * wgk(i)*frac(iel)*log(frac(iel))
        endif
        if ((UM - frac(iel)) > SMALL) then
          elects = elects - 2*tempau * wgk(i)*(UM-frac(iel))*log(UM-frac(iel))
        endif
      enddo
    enddo
  endif

! calculates band energy

  eband = ZERO

  iel = 0
  do i=1,nrk
    do j=1,nband(i)
      iel = iel + 1
      eband = eband + 2*wgk(i)*frac(iel)*el(iel)
    enddo
  enddo

! estimates occupied band width
! beware of ztot=0

  bandwid = efermi - el(ind(1))
  if(bandwid < EPS) then
    if(iel > 1) then
      bandwid = el(ind(2)) - el(ind(1))
    else
      bandwid = 0.1
    endif
  endif

! estimates the penn gap

  penngap = ZERO
  jmax = nint(ztot/(2*UM))
  if(abs(ztot - real(2*jmax)) < EPS .and. jmax > 0) then
    imx = 0
    iel = 0
    do i=1,nrk
      penngap = penngap + el(iel+jmax+1) - el(iel+jmax)
      if(el(iel+jmax+1) - efermi < -SMALL) imx = imx + 1
      if(efermi - el(iel+jmax) < -SMALL) imx = imx + 1
      iel = iel + nband(i)
    enddo
    if(imx == 0) then
      penngap = penngap / real(nrk)
    else
      penngap = ZERO
    endif
  endif

  deallocate(jrk)
  deallocate(ind)

  return

end subroutine fermi_level
