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
!>  \date         early eighties. 5 October 2026.
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
!   Fit of each structure in subroutine eqst_fit. 5 October 2026. JLM+claude


  implicit none

  integer, parameter  ::  REAL64 = selected_real_kind(12)


  integer                   ::  mxdstr                                   !<  maximum number of structures
  integer                   ::  mxdnpt                                   !<  maximum number of energy points (per structure)

  character(len=5)          ::  ftype                                    !<  type of equation of state, MURNA or BIRCH
  integer                   ::  nstr                                     !<  number of structures

! for comparing several structures

  real(REAL64), allocatable ::  xarall(:,:)                              !<  cell volume (bohr^3)
  real(REAL64), allocatable ::  yarall(:,:)                              !<  energy in Rydberg :-(
  integer, allocatable      ::  nptall(:)                                !<  number of points for each structure
  integer, allocatable      ::  nform(:)                                 !<  number of formula units in the cell
  real(REAL64), allocatable ::  volfacall(:)                             !<  volume is volfac*alatt**3 (see eqst_read_data)
  character(len=10), allocatable  ::  label(:)                           !<  label of the structure

! present structure

  integer                   ::  npt                                      !<  number of calculated energies
  real(REAL64), allocatable ::  xar(:)                                   !<  cell volume (bohr^3)
  real(REAL64), allocatable ::  yar(:)                                   !<  energy (Rydberg)

! results for each structure, for the plots

  real(REAL64), allocatable ::  wvzero(:)                                !<  equilibrium volume (bohr^3)
  real(REAL64), allocatable ::  wezero(:)                                !<  energy at equilibrium (Rydberg)
  real(REAL64), allocatable ::  wbzero(:)                                !<  bulk modulus (Rydberg/bohr^3)
  real(REAL64), allocatable ::  wbprim(:)                                !<  pressure derivative of the bulk modulus
  real(REAL64), allocatable ::  xarz(:,:)                                !<  cell volume (bohr^3)
  real(REAL64), allocatable ::  yarz(:,:)                                !<  energy (Rydberg)
  integer, allocatable      ::  nptz(:)                                  !<  number of calculated energies

! results of the fit

  real(REAL64)              ::  bprim                                    !<  B'
  real(REAL64)              ::  bzero                                    !<  B(0)
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

! counters

  integer                   ::  ns, jj


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


  do ns = 1,nstr

    xar(:) = xarall(:,ns)
    yar(:) = yarall(:,ns)
    npt = nptall(ns)

!   fit of the equation of state and printing of the results

    call eqst_fit(ftype, npt, volfacall(ns), xar, yar, 1,                &
        ezero, vzero, bzero, bprim,                                      &
        mxdnpt)

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
  deallocate(xar, yar)
  deallocate(wvzero, wezero, wbzero, wbprim)
  deallocate(xarz, yarz, nptz)

  stop

end program eqst
