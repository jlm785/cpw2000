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

!>  Reads the energy as a function of lattice constant (or volume)
!>  for several structures from the data file.
!>  Writes the data to a file for reuse.
!>
!>  The line with the volume factor is
!>    volfac  [nform]  [label]
!>  where nform, the number of formula units in the cell, is optional
!>  (default 1) for compatibility with older files.  The second item
!>  is taken as nform only if it is a positive integer (digits only),
!>  so a label cannot be made only of digits.
!>
!>  \author       David M. Wood, Jose Luis Martins
!>  \version      5.13
!>  \date         31 March 2026, 4 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_read_data(iodata, ioreplay, filename, ftype, nstr,       &
        volfacall, nform, label, nptall, xarall, yarall,                 &
        mxdstr, mxdnpt)

! extracted from the main program. 31 March 2026. JLM
! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! Reads from file opened by eqst_size_data instead of standard input. 3 October 2026. JLM+claude
! Optional number of formula units, nform. 4 October 2026. JLM+claude


  implicit none

  integer, parameter  ::  REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                   ::  mxdstr                       !<  maximum number of structures
  integer, intent(in)                   ::  mxdnpt                       !<  maximum number of energy points (per structure)

  integer, intent(in)                   ::  iodata                       !<  tape number of the (open) data file

  character(len=*), intent(in)          ::  filename                     !<  writes input for future reuse
  integer, intent(in)                   ::  ioreplay                     !<  tape number of filename

!output

  character(len=5), intent(out)         ::  ftype                        !<  type of function

  integer, intent(out)                  ::  nstr                         !<  number of structures
  real(REAL64), intent(out)             ::  volfacall(mxdstr)            !<  volume is volfac*alatt**3 if zero, volume is read, if negative it is the are of the base in epitaxial
  character(len=10), intent(out)        ::  label(mxdstr)                !<  label of the structure

  integer, intent(out)                  ::  nptall(mxdstr)               !<  number of calculated energies
  integer, intent(out)                  ::  nform(mxdstr)                !<  number of formula units in the cell

  real(REAL64), intent(out)             ::  xarall(mxdnpt,mxdstr)        !<  Cell volume
  real(REAL64), intent(out)             ::  yarall(mxdnpt,mxdstr)        !<  Energy in Rydberg:-(

! local variables

  character(len=100)        ::  fline

  real(REAL64)              ::  alat(mxdnpt)
  integer                   ::  ioerr
  character(len=100)        ::  word                                     !  second item in the volfac line
  real(REAL64)              ::  vtmp                                     !  volfac read again

! constants

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter   ::  EPS = 1.0E-8_REAL64

! counters

  integer                   ::  ns, ii, mm




  open(unit = ioreplay, file = adjustl(trim(filename)),                  &
       status  ='UNKNOWN', form='FORMATTED')

! reads the number of structures and type of fit

  read(iodata,*) ftype
  if(ftype /= 'MURNA' .and. ftype /= 'murna' .and.                       &
     ftype /= 'BIRCH' .and. ftype /= 'birch') then

    write(6,*)
    write(6,*) '  Unknown equation of state ',ftype
    write(6,*) '  Will assume Murnaghan equation of state'
    write(6,*)
    ftype = 'MURNA'

  endif
  write(ioreplay,'(2x,a5,10x,"Type of equation of state")') ftype

  read(iodata,*) nstr
  if(nstr > mxdstr) then
    write(6,*)
    write(6,*) '  STOPPED in eqst_read_data:  nstr = ',nstr,' > mxdstr = ',mxdstr
    write(6,*)

    stop

  endif
  write(ioreplay,'(2x,i5,10x,"Number of structures")') nstr

  write(6,*)
  write(6,'("  Equation of state: ",a5,"     number of structures: ",i5)') ftype, nstr
  write(6,*)

! loop over structures

  do ns = 1,nstr

!   read number of points

    read(iodata,*) nptall(ns)
    if(nptall(ns) > mxdnpt) then
      write(6,*)
      write(6,*) '  STOPPED in eqst_read_data:  npt = ',nptall(ns),' > mxdnpt = ',mxdnpt
      write(6,*)

      stop

    endif
    if(nptall(ns) < 4) then
      write(6,*)
      write(6,*) '  STOPPED in eqst_read_data:  npt = ',nptall(ns),' < 4 '
      write(6,*)

      stop

    endif
    write(ioreplay,'(2x,i5,10x,"Number of points for structure",i5)') nptall(ns), ns

!   volfacall is the volume of primitive cell in (lattice constant)**3,
!   e.g. for fcc bravais lattice, volfacall = .25
!   if volfacall is zero volume instead of lattice constants are read
!   if it is negative onle the height of cell is being changed

    do ii = 1,100
      fline(ii:ii) = ' '
    enddo

    read(iodata,'(a100)') fline

    read(fline,*,iostat=ioerr) volfacall(ns)
    if(ioerr /= 0) then
      write(6,*)
      write(6,*) '  STOPPED in eqst_read_data:  cannot read volfac for structure ',ns
      write(6,*) '  from line: ',trim(fline)
      write(6,*)

      stop

    endif

!   defaults for old files

    nform(ns) = 1
    if(ns < 10) then
      write(label(ns),'("struct-N.",i1)') ns
    elseif(ns < 100) then
      write(label(ns),'("structN.",i2)') ns
    else
      label(ns) = '          '
    endif

!   second item is either nform (digits only) or the label

    read(fline,*,iostat=ioerr) vtmp, word
    if(ioerr == 0) then
      if(verify(trim(word),'0123456789') == 0) then
        read(word,*) nform(ns)
        read(fline,*,iostat=ioerr) vtmp, word, label(ns)
      else
        label(ns) = word(1:10)
      endif
    endif

    if(nform(ns) < 1) then
      write(6,*)
      write(6,*) '  STOPPED in eqst_read_data:  number of formula units ',nform(ns)
      write(6,*) '  for structure ',ns
      write(6,*)

      stop

    endif

    write(ioreplay,'(3x,f15.6,3x,i5,3x,a10,10x,                          &
        &     "volfac, nform, label")') volfacall(ns), nform(ns), label(ns)

    write(6,'("  Structure ",i3,2x,a10,"   volfac = ",f12.6,             &
        &     "   formula units: ",i3,"   number of points: ",i5)')      &
        ns, label(ns), volfacall(ns), nform(ns), nptall(ns)

!   reads in the data points. the lattice constant is
!   in atomic units, the energy in Hartree and converted to Rydberg.
!   For angstroms exchange commented lines
!   If volfac is positive the data lines have lattice constants,
!   otherwise they have primitive cell volumes.

    do mm = 1,nptall(ns)

      read (iodata,*) alat(mm),yarall(mm,ns)

      if(volfacall(ns) < EPS) then
        xarall(mm,ns) = alat(mm)
!       xar(mm) = alat(mm)/0.529177**3
        write(ioreplay,'(2(3x,f17.8),10x,"volume energy for point ",     &
            &     2i5)') alat(mm), yarall(mm,ns), ns, mm
      else
        xarall(mm,ns) = volfacall(ns)*(alat(mm))**3
        write(ioreplay,'(2(3x,f17.8),10x,                                &
            &     "lattice const. energy for point ",2i5)')              &
            alat(mm), yarall(mm,ns), ns, mm
      endif
!     back to Rydberg :-(
      yarall(mm,ns)  =  2*yarall(mm,ns)

    enddo

  enddo

  write(6,*)

end subroutine eqst_read_data
