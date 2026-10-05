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

!>    Fits an equation of state (Birch or Murnaghan)
!>    to a few E(V) data points for several structures.
!>    Different ways of ploting the results are available.

!>    Murnaghan equation of state:
!>    F. D. Murnaghan, Proc. Natl. Acad. Sci. 30, 244 (1944).
!>    Given a set of ordered pairs (x(i),y(i)),i = 1,npt
!>    where npt is the total number of points provided to the fit.
!>    Fit is to the form
!>    y = total energy per primitive cell in Hartree
!>      =  a + b*x + c*x**d
!>    where x  =  volume of primitive cell in cubic atomic units

!>    Birch equation of state:
!>    F. Birch J. Geophys. Res. 83, 1257 (1978).
!>    The funcyion y is different...
!>
!>  \author       David M. Wood, Jose Luis Martins
!>  \version      5.13
!>  \date         early eighties. 4 October 2026.
!>  \copyright    GNU Public License v2


program eqst

!   The Murnaghan fit was written by D. Wood
!   The Birch fit was added by jlm. Sep 89
!   Changed to hartree jlm. jun 04
!   Changed to f90 and gnuplot. June 2022. JLM
!   Interactive prompts, prints lattice constant in angstroms.  19 April 2024. JLM
!   Read data in separate subroutine. Epitaxial situation.
!
!   Major rewriting, subroutines separated and prefixed by eqst_ or num_. 3 October 2026. JLM+claude
!   num_brent and num_zbrent as reverse communication interfaces.  3 October 2026. claude
!   Number of formula units nform, used only to compare structures. 4 October 2026. JLM+claude
!   fbirch and fmurna replaced by subroutines. 5 Octiber 2026. JLM


  implicit none

  integer, parameter  ::  REAL64 = selected_real_kind(12)


  integer                   ::  mxdstr                                   !<  maximum number of structures
  integer                   ::  mxdnpt                                   !<  maximum number of energy points (per structure)

  character(len=5)          ::  ftype                                    !<  type of equation of state, MURNA or BIRCH
  integer                   ::  nstr                                     !<  number of structures

! fom comparing several structures

  real(REAL64), allocatable ::  xarall(:,:)                              !<  cell volume (bohr^3)
  real(REAL64), allocatable ::  yarall(:,:)                              !<  energy in Rydberg :-(
  integer, allocatable      ::  nptall(:)                                !<  number of points for each structure
  integer, allocatable      ::  nform(:)                                 !<  number of formula units in the cell
  real(REAL64), allocatable ::  volfacall(:)                             !<  volume is volfac*alatt**3 (see eqst_read_data)
  character(len=10), allocatable  ::  label(:)                           !<  label of the structure

! present structure

  integer                   ::  npt                                      !<  number of calculated energies
  real(REAL64)              ::  fn                                       !<  real(npt)
  real(REAL64)              ::  volfac                                   !<  volume is volfac*alatt**3
  real(REAL64), allocatable ::  xar(:)                                   !<  cell volume (bohr^3)
  real(REAL64), allocatable ::  yar(:)                                   !<  energy (Rydberg)
  real(REAL64), allocatable ::  xscl(:)                                  !<  scaled volume (same as xar)
  real(REAL64), allocatable ::  yscl(:)                                  !<  scaled energy
  real(REAL64), allocatable ::  alat(:)                                  !<  lattice constant (or cell height) for printing

! results for each structure, for the plots

  real(REAL64), allocatable ::  wvzero(:)                                !<  equilibrium volume (bohr^3)
  real(REAL64), allocatable ::  wezero(:)                                !<  energy at equilibrium (Rydberg)
  real(REAL64), allocatable ::  wbzero(:)                                !<  bulk modulus (Rydberg/bohr^3)
  real(REAL64), allocatable ::  wbprim(:)                                !<  pressure derivative of the bulk modulus
  real(REAL64), allocatable ::  xarz(:,:)                                !<  cell volume (bohr^3)
  real(REAL64), allocatable ::  yarz(:,:)                                !<  energy (Rydberg)
  integer, allocatable      ::  nptz(:)                                  !<  number of calculated energies

! fits

  real(REAL64)              ::  a                                        !<  Murnaghan fit, y = a + b x + c x**dd
  real(REAL64)              ::  b                                        !<  Murnaghan fit, y = a + b x + c x**dd
  real(REAL64)              ::  c                                        !<  Murnaghan fit, y = a + b x + c x**dd
  real(REAL64)              ::  dd                                       !<  Murnaghan fit, y = a + b x + c x**dd

  real(REAL64)              ::  fval                                     !<  value of the function for num_zbrent and num_brent
  integer                   ::  iflag                                    !<  reverse communication flag

  real(REAL64)              ::  y                                        !<  energy of a point
  real(REAL64)              ::  yav                                      !<  average energy
  real(REAL64)              ::  y2av                                     !<  average of energy squared
  real(REAL64)              ::  vary                                     !<  variance of the energy
  real(REAL64)              ::  sqvry                                    !<  square root of the variance
  real(REAL64)              ::  below                                    !<  lower end of the search interval
  real(REAL64)              ::  above                                    !<  upper end of the search interval
  real(REAL64)              ::  emin                                     !<  lowest scaled energy
  integer                   ::  nmin                                     !<  point with the lowest energy
  real(REAL64)              ::  vmin                                     !<  volume of the point with the lowest energy
  real(REAL64)              ::  sq                                       !<  sum of the squares of the fit error
  real(REAL64)              ::  qq                                       !<  fitted scaled energy
  real(REAL64)              ::  dy                                       !<  fit error
  real(REAL64)              ::  percnt                                   !<  fit error in percentage

  real(REAL64)              ::  goodns                                   !<  quality of fit
  real(REAL64)              ::  qual                                     !<  typical error of the energy
  real(REAL64)              ::  aeq                                      !<  equilibrium lattice constant
  real(REAL64)              ::  bprim                                    !<  B'
  real(REAL64)              ::  bzero                                    !<  B(0)
  real(REAL64)              ::  bulkmd                                   !<  Bulk modulus
  real(REAL64)              ::  ezero                                    !<  E(aeq)
  real(REAL64)              ::  vzero                                    !<  equilibrium volume

  real(REAL64)              ::  pmax                                     !<  maximum pressure for the plots (GPa)
  integer                   ::  npts                                     !<  number of points in the plots

  character(len=15)         ::  filename                                 !<  writes input for future reuse
  integer                   ::  ioreplay                                 !<  tape number of filename

  character(len=200)        ::  datafile                                 !<  file with the input data
  integer                   ::  iodata                                   !<  tape number of datafile

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

  integer                   ::  ns, ii, jj, n


! file stuff

  filename = 'replay_eqst.dat'
  ioreplay = 15

  iodata = 14

! opens the data file and finds the dimensions

  call eqst_size_data(iodata, datafile, mxdstr, mxdnpt)

! allocations

  allocate(xarall(mxdnpt,mxdstr))
  allocate(yarall(mxdnpt,mxdstr))
  allocate(nptall(mxdstr))
  allocate(nform(mxdstr))
  allocate(volfacall(mxdstr))
  allocate(label(mxdstr))

  allocate(xar(mxdnpt))
  allocate(yar(mxdnpt))
  allocate(xscl(mxdnpt))
  allocate(yscl(mxdnpt))
  allocate(alat(mxdnpt))

  allocate(wvzero(mxdstr))
  allocate(wezero(mxdstr))
  allocate(wbzero(mxdstr))
  allocate(wbprim(mxdstr))
  allocate(xarz(mxdnpt,mxdstr))
  allocate(yarz(mxdnpt,mxdstr))
  allocate(nptz(mxdstr))

  xarall(:,:) = ZERO
  yarall(:,:) = ZERO

  call eqst_read_data(iodata, ioreplay, filename, ftype, nstr,           &
      volfacall, nform, label, nptall, xarall, yarall,                   &
      mxdstr, mxdnpt)


!   evaluate ybar,y2bar for use in the variance
!   and average of each variable, in order to scale the array y
!   to make the least squares algorithm more stable

  do ns = 1,nstr

    xar(:) = xarall(:,ns)
    yar(:) = yarall(:,ns)
    npt = nptall(ns)
    volfac = volfacall(ns)

!   lattice constant (or height of the cell in the epitaxial case) for printing

    do n = 1,npt
      if(volfac > EPS) then
        alat(n) = (xar(n)/volfac)**(UM/3)
      elseif(volfac < -EPS) then
        alat(n) = xar(n)/abs(volfac)
      else
        alat(n) = ZERO
      endif
    enddo

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

    write(6,*)
    write(6,'(2x,a5," fit from ",i3," (x,y) data pairs")') ftype, npt
    write(6,'( 5x,"mean of y = ",f10.5,"  variance = ",f10.5,5x,         &
        &   "---y array are re-scaled")') yav, sqvry
    write(6,*)
    write(6,*)

    do ii = 1,npt
      xscl(ii) = xar(ii)
      yscl(ii) = (yar(ii)-yav)/sqvry
    enddo

!   fit of the equation of state

    if(ftype == 'murna' .or. ftype == 'MURNA') then

!     dd is evaluated from a non linear fit
!     limits for the search of optimum bprim
!     eqst_fmurna(below)*eqst_fmurna(above) < 0

      below  =  UM-BPMAX
      above  =  UM-BPMIN

      iflag = 0
      do
        call num_zbrent(dd, fval, below, above, EPSZB, iflag)
        if(iflag == 0) exit
        call eqst_fmurna(fval, dd, npt, xscl, yscl, a, b, c, mxdnpt)
      enddo

!     a, b, c for the final value of dd

      call eqst_fmurna(fval, dd, npt, xscl, yscl, a, b, c, mxdnpt)

    else

!     estimates the equilibrium volume

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

!     ezero, bzero, bprim for the final value of vzero

      call eqst_fbirch(fval, vzero, npt, xar, yscl, ezero, bzero, bprim, mxdnpt)
    endif

    write(6,'("  cols are x(in), y(in), alat, y(in,scaled), y(fit,scaled), and % diff")')
    write(6,*)


!   use sq for sum of (y-yscl)**2

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
      write(6,'(2x,f10.5,4x,2(f10.5,2x),2x,2(f10.5,2x),4x,               &
          &     f10.5,3x,f8.3)') xar(n), yar(n)/2, alat(n), yscl(n), qq, percnt
    enddo

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

    if(ftype == 'murna' .or. ftype == 'MURNA') then

!     now find the coefficients a,b,c,d for the original, unscaled y variabl

      a = yav+a*sqvry
      b = sqvry*b
      c = sqvry*c
      write(6,'(5x,"for unscaled quantities, y = a+bx+c*x**d where ")')
      write(6,'(4x,"a,b,c,d= ",1pe18.10,3e18.10)') a,b,c,dd
      write(6,*)
      bprim = UM-dd
      bzero = b*bprim
      vzero = (-c*dd/b)**(UM/bprim)
      ezero = a-b*bprim*vzero/dd
    else
      ezero = yav+ezero*sqvry
      bzero = bzero*sqvry
    endif

!   now calculate the actual bulk modulus in gpa
    bulkmd = GPA*bzero
!     if(volfac > ZERO) then
!       aeq = (vzero/volfac)**(1./3.)
!     else
!       aeq = vzero**(1./3.)
!     endif

    if(abs(volfac) < EPS) then

      aeq = vzero**(UM/(3*UM))

      write(6,*)
      write(6,'("      Veq(a.u.)=",f9.3,"   E0(Hartree)=",f10.5)')       &
               vzero, ezero/2
      write(6,'("   B0(GPa)=",f8.2,"   B0PRIM= ",f6.2)') bulkmd, bprim
      write(6,*)
      write(6,'("  Cubic root of volume in Angstroms: ",f10.5)') aeq*BOHR
      write(6,*)

    elseif(volfac > ZERO) then

      aeq = (vzero/abs(volfac))**(UM/(3*UM))

      write(6,*)
      write(6,'("   aeq(a.u.)=",f8.4,"   Veq(a.u.)=",f9.3,               &
          &     "   E0(Hartree)=",f10.5)') aeq, vzero, ezero/2
      write(6,'("   B0(GPa)=",f8.2,"   B0PRIM= ",f6.2)') bulkmd, bprim
      write(6,*)
      write(6,'("  Lattice constant in Angstroms: ",f10.5)') aeq*BOHR
      write(6,*)

    else

      aeq = vzero/abs(volfac)

      write(6,*)
      write(6,'("   ceq(a.u.)=",f8.4,"   Veq(a.u.)=",f9.3,               &
          &     "   E0(Hartree)=",f10.5)') aeq, vzero, ezero/2
      write(6,'("   Elastic constant(GPa)=",f8.2,"   B0PRIM= ",f6.2)') bulkmd, bprim
      write(6,*)
      write(6,'("  Perpendicular lattice constant in Angstroms: ",f10.5)') aeq*BOHR
      write(6,*)

    endif

!   stores for plot subroutines, where the structures are compared.
!   Energies and volumes per formula unit (bzero and bprim are intensive)

    wvzero(ns)  =  vzero / nform(ns)
    wezero(ns)  =  ezero / nform(ns)
    wbzero(ns)  =  bzero
    wbprim(ns)  =  bprim
    nptz(ns)  =  npt
    do jj = 1,npt
      xarz(jj,ns) = xar(jj) / nform(ns)
      yarz(jj,ns) = yar(jj) / nform(ns)
    enddo
  enddo

  if(nstr > 1) then
    read(iodata,*) pmax,npts

    write(6,*)
    write(6,'("  Maximum pressure (GPa) and number of points for the plots: ",f10.3,i8)') pmax, npts
    write(6,*)

    if(pmax < ZERO .or. pmax > 5000*UM) then
      pmax = 100*UM
      write(6,*)
      write(6,*) '   Unreasonable value of pmax. Set to 100 GPa'
      write(6,*)
    endif
    if(npts < 10 .or. npts > 5000) then
      npts = 50
      write(6,*)
      write(6,*) '   Unreasonable value of npts. Set to 50 GPa'
      write(6,*)
    endif
    write(ioreplay,'(3x,f12.3,3x,i8,10x,                                 &
        &     "maximum pressure (GPa), number of points")') pmax, npts

    call eqst_gpeqst(wezero, wvzero, wbzero, wbprim, nstr, label,        &
        ftype, pmax, npts, mxdstr)

  endif

  close(unit = iodata)
  close(unit = ioreplay)

  call eqst_gpeofv(xarz, yarz, nptz, wezero, wvzero, wbzero,             &
      wbprim, nstr, label, ftype,                                        &
      mxdstr, mxdnpt)

  deallocate(xarall, yarall, nptall, nform, volfacall, label)
  deallocate(xar, yar, xscl, yscl, alat)
  deallocate(wvzero, wezero, wbzero, wbprim)
  deallocate(xarz, yarz, nptz)

  stop

end program eqst
