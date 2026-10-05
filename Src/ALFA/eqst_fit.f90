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

!>  Fits the Murnaghan or Birch equation of state to the
!>  energies E(V) of one structure and prints the equilibrium
!>  volume, energy, bulk modulus and its pressure derivative.
!>
!>  The energies are scaled to make the least squares fit more stable.
!>
!>  \author       David M. Wood, Jose Luis Martins
!>  \version      5.13
!>  \date         early eighties, 5 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_fit(ftype, npt, volfac, xar, yar, iprint,                &
    ezero, vzero, bzero, bprim,                                          &
    mxdnpt)

! Extracted from the main program eqst.f90. 5 October 2026. JLM+claude

  implicit none

  integer, parameter  ::  REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                   ::  mxdnpt                       !<  array dimension of number of points

  character(len=5), intent(in)          ::  ftype                        !<  type of equation of state, MURNA or BIRCH
  integer, intent(in)                   ::  npt                          !<  number of calculated energies
  real(REAL64), intent(in)              ::  volfac                       !<  volume is volfac*alatt**3, if zero volume, if negative area of the base (epitaxial)
  real(REAL64), intent(in)              ::  xar(mxdnpt)                  !<  cell volume (bohr^3)
  real(REAL64), intent(in)              ::  yar(mxdnpt)                  !<  energy (Rydberg)
  integer, intent(in)                   ::  iprint                       !<  0: prints only the results, 1: also prints the details of the fit

! output

  real(REAL64), intent(out)             ::  ezero                        !<  energy at equilibrium (Rydberg)
  real(REAL64), intent(out)             ::  vzero                        !<  equilibrium volume (bohr^3)
  real(REAL64), intent(out)             ::  bzero                        !<  bulk modulus (Rydberg/bohr^3)
  real(REAL64), intent(out)             ::  bprim                        !<  pressure derivative of the bulk modulus

! local variables

  real(REAL64)              ::  xscl(mxdnpt)                             !  scaled volume (same as xar)
  real(REAL64)              ::  yscl(mxdnpt)                             !  scaled energy
  real(REAL64)              ::  alat(mxdnpt)                             !  lattice constant (or cell height) for printing

  real(REAL64)              ::  fn                                       !  real(npt)

  real(REAL64)              ::  a                                        !  Murnaghan fit, y = a + b x + c x**dd
  real(REAL64)              ::  b                                        !  Murnaghan fit, y = a + b x + c x**dd
  real(REAL64)              ::  c                                        !  Murnaghan fit, y = a + b x + c x**dd
  real(REAL64)              ::  dd                                       !  Murnaghan fit, y = a + b x + c x**dd

  real(REAL64)              ::  fval                                     !  value of the function for num_zbrent and num_brent
  integer                   ::  iflag                                    !  reverse communication flag

  real(REAL64)              ::  y                                        !  energy of a point
  real(REAL64)              ::  yav                                      !  average energy
  real(REAL64)              ::  y2av                                     !  average of energy squared
  real(REAL64)              ::  vary                                     !  variance of the energy
  real(REAL64)              ::  sqvry                                    !  square root of the variance
  real(REAL64)              ::  below                                    !  lower end of the search interval
  real(REAL64)              ::  above                                    !  upper end of the search interval
  real(REAL64)              ::  emin                                     !  lowest scaled energy
  integer                   ::  nmin                                     !  point with the lowest energy
  real(REAL64)              ::  vmin                                     !  volume of the point with the lowest energy
  real(REAL64)              ::  sq                                       !  sum of the squares of the fit error
  real(REAL64)              ::  qq                                       !  fitted scaled energy
  real(REAL64)              ::  dy                                       !  fit error
  real(REAL64)              ::  percnt                                   !  fit error in percentage

  real(REAL64)              ::  goodns                                   !  quality of fit
  real(REAL64)              ::  qual                                     !  typical error of the energy
  real(REAL64)              ::  aeq                                      !  equilibrium lattice constant
  real(REAL64)              ::  bulkmd                                   !  bulk modulus in GPa

! constants

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter   ::  BOHR = 0.5291772105_REAL64
  real(REAL64), parameter   ::  EPS = 1.0E-8_REAL64
  real(REAL64), parameter   ::  AUTOGPA = 29421.0158_REAL64              !  Hartree/bohr^3 to GPa
  real(REAL64), parameter   ::  GPA = AUTOGPA/2                          !  Rydberg/bohr^3 to GPa

  real(REAL64), parameter   ::  BPMIN = 2*UM                             !  lower bound for B' in the Murnaghan fit
  real(REAL64), parameter   ::  BPMAX = 7*UM                             !  upper bound for B' in the Murnaghan fit
  real(REAL64), parameter   ::  EPSZB = 1.0E-4_REAL64                    !  tolerance (in dd) for num_zbrent
  real(REAL64), parameter   ::  EPSMU = 1.0E-6_REAL64                    !  tolerance (in volume) for num_brent

! counters

  integer                   ::  ii, n


! lattice constant (or height of the cell in the epitaxial case) for printing

  do n = 1,npt
    if(volfac > EPS) then
      alat(n) = (xar(n)/volfac)**(UM/3)
    elseif(volfac < -EPS) then
      alat(n) = xar(n)/abs(volfac)
    else
      alat(n) = ZERO
    endif
  enddo

! evaluate ybar,y2bar for use in the variance
! and average of each variable, in order to scale the array y
! to make the least squares algorithm more stable

  fn = UM*npt

  yav = ZERO
  y2av = ZERO
  do ii = 1,npt
    y = yar(ii)
    yav = yav+y
    y2av = y2av+y*y
  enddo

  vary = (y2av-yav*yav/fn)/fn
  yav = yav/fn
  sqvry = sqrt(vary)

  if(iprint > 0) then
    write(6,*)
    write(6,'(2x,a5," fit from ",i3," (x,y) data pairs")') ftype, npt
    write(6,'( 5x,"mean of y = ",f10.5,"  variance = ",f10.5,5x,         &
        &   "---y array are re-scaled")') yav, sqvry
    write(6,*)
    write(6,*)
  endif

  do ii = 1,npt
    xscl(ii) = xar(ii)
    yscl(ii) = (yar(ii)-yav)/sqvry
  enddo

! fit of the equation of state

  if(ftype == 'murna' .or. ftype == 'MURNA') then

!   dd is evaluated from a non linear fit
!   limits for the search of optimum bprim
!   eqst_fmurna(below)*eqst_fmurna(above) < 0

    below  =  UM-BPMAX
    above  =  UM-BPMIN

    iflag = 0
    do
      call num_zbrent(dd, fval, below, above, EPSZB, iflag)
      if(iflag == 0) exit
      call eqst_fmurna(fval, dd, npt, xscl, yscl, a, b, c, mxdnpt)
    enddo

!   a, b, c for the final value of dd

    call eqst_fmurna(fval, dd, npt, xscl, yscl, a, b, c, mxdnpt)

  else

!   estimates the equilibrium volume

    nmin = 1
    emin = yscl(1)
    do n = 1,npt
      if(yscl(n) < emin) then
        nmin = n
        emin = yscl(n)
      endif
    enddo
    vmin = xar(nmin)
    below = 0.75_REAL64*vmin
    above = 1.5_REAL64*vmin
    do n = 1,npt
      if(xar(n) > below .and. xar(n) < vmin-EPS) below = xar(n)
      if(xar(n) < above .and. xar(n) > vmin+EPS) above = xar(n)
    enddo

    iflag = 0
    do
      call num_brent(vzero, fval, below, above, EPSMU, iflag)
      if(iflag == 0) exit
      call eqst_fbirch(fval, vzero, npt, xar, yscl, ezero, bzero, bprim, mxdnpt)
    enddo

!   ezero, bzero, bprim for the final value of vzero

    call eqst_fbirch(fval, vzero, npt, xar, yscl, ezero, bzero, bprim, mxdnpt)
  endif

! use sq for sum of (y-yscl)**2

  if(iprint > 0) then
    write(6,'("  cols are x(in), y(in), alat, y(in,scaled), y(fit,scaled), and % diff")')
    write(6,*)
  endif

  sq = ZERO
  do n = 1,npt
    if(ftype == 'murna' .or. ftype == 'MURNA') then
      qq = a + b*xscl(n) + c*xscl(n)**dd
    else
      call eqst_eofv(qq, xscl(n), ezero, vzero, bzero, bprim)
    endif
    dy = yscl(n)-qq
    percnt = 100.0_REAL64*dy
    sq = sq+dy*dy
    if(iprint > 0) then
      write(6,'(2x,f10.5,4x,2(f10.5,2x),2x,2(f10.5,2x),4x,               &
          &     f10.5,3x,f8.3)') xar(n), yar(n)/2, alat(n), yscl(n), qq, percnt
    endif
  enddo

  if(iprint > 0) then
    if(npt > 4) then
      goodns = sqrt(sq/real(npt-4,REAL64))

      write(6,*)
      write(6,'("  goodness of fit = scaled variance = (mn sqr)/(npt-4) = ")')
      write(6,*)
      write(6,'(3x,f12.8)') goodns
      write(6,*)

      qual = sqvry*goodns/sqrt(real(npt,REAL64))
      write(6,*)
      write(6,'("   typical error in y   ",e14.6)') qual
      write(6,*)
    else
      write(6,*)
      write(6,'("  number of data pairs read in  =  4  ==>  perfect fit")')
      write(6,*)
    endif
  endif

  if(ftype == 'murna' .or. ftype == 'MURNA') then

!   now find the coefficients a,b,c,d for the original, unscaled y variable

    a = yav+a*sqvry
    b = sqvry*b
    c = sqvry*c
    if(iprint > 0) then
      write(6,'(5x,"for unscaled quantities, y = a+bx+c*x**d where ")')
      write(6,'(4x,"a,b,c,d= ",1pe18.10,3e18.10)') a,b,c,dd
      write(6,*)
    endif
    bprim = UM-dd
    bzero = b*bprim
    vzero = (-c*dd/b)**(UM/bprim)
    ezero = a-b*bprim*vzero/dd
  else
    ezero = yav+ezero*sqvry
    bzero = bzero*sqvry
  endif

! now calculate the actual bulk modulus in gpa

  bulkmd = GPA*bzero

  if(abs(volfac) < EPS) then

    aeq = vzero**(UM/(3*UM))

    write(6,*)
    write(6,'("      Veq(a.u.)=",f9.3,"   E0(Hartree)=",f10.5)')         &
             vzero, ezero/2
    write(6,'("   B0(GPa)=",f8.2,"   B0PRIM= ",f6.2)') bulkmd, bprim
    write(6,*)
    write(6,'("  Cubic root of volume in Angstroms: ",f10.5)') aeq*BOHR
    write(6,*)

  elseif(volfac > ZERO) then

    aeq = (vzero/abs(volfac))**(UM/(3*UM))

    write(6,*)
    write(6,'("   aeq(a.u.)=",f8.4,"   Veq(a.u.)=",f9.3,                 &
        &     "   E0(Hartree)=",f10.5)') aeq, vzero, ezero/2
    write(6,'("   B0(GPa)=",f8.2,"   B0PRIM= ",f6.2)') bulkmd, bprim
    write(6,*)
    write(6,'("  Lattice constant in Angstroms: ",f10.5)') aeq*BOHR
    write(6,*)

  else

    aeq = vzero/abs(volfac)

    write(6,*)
    write(6,'("   ceq(a.u.)=",f8.4,"   Veq(a.u.)=",f9.3,                 &
        &     "   E0(Hartree)=",f10.5)') aeq, vzero, ezero/2
    write(6,'("   Elastic constant(GPa)=",f8.2,"   B0PRIM= ",f6.2)') bulkmd, bprim
    write(6,*)
    write(6,'("  Perpendicular lattice constant in Angstroms: ",f10.5)') aeq*BOHR
    write(6,*)

  endif

  return

end subroutine eqst_fit
