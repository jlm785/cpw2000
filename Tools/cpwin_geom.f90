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

!>  Performs the analysis of the crystal geometry
!>  from the cpw.in, cpw.out or other esdf file
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         22 April 2021, 28 October 2026.
!>  \copyright    GNU Public License v2

program cpwin_geom

! Added call to sym_space_group_name. 23 January 2022. JLM
! Interactive, analysis in voronoi_analysis. 28 October 2026. JLM+claude
! Choice of file (cpw.in, cpw.out, other). 28 October 2026. JLM+claude

  use esdf

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! dimensions

  integer                            ::  mxdtyp                          !<  array dimension of types of atoms
  integer                            ::  mxdatm                          !<  array dimension of number of atoms of a given type

! crystal variables

  real(REAL64)                       ::  alatt                           !<  lattice constant

  real(REAL64)                       ::  adot(3,3)                       !<  metric in direct space in atomic units (Bohr radius)

  integer                            ::  ntype                           !<  number of types of atoms
  integer, allocatable               ::  natom(:)                        !<  number of atoms of type i
  real(REAL64), allocatable          ::  rat(:,:,:)                      !<  lattice coordinates of atom j of type i
  real(REAL64), allocatable          ::  atmass(:)                       !<  atomic mass of atoms of type i
  character(len=2), allocatable      ::  nameat(:)                       !<  chemical symbol for the type i

! other variables

  integer               ::  ipr
  character(len=200)    ::  fcpwin                                       !<  filename with input data (new style)
  logical               ::  lgeom, lgeom2
  logical               ::  lex                                          !<  file exists
  integer               ::  ichoice, iblank

  integer               ::  ioreplay = 35                                !<  hard coded tape unit for replay
  character(len=15)     ::  fname = 'replay_geom.dat'                    !<  hard coded filename for replay

  ipr = 3

  open(unit = ioreplay, file = fname,  status = 'UNKNOWN',               &
                    form = 'FORMATTED')

! choice of file with the crystal structure

  write(6,*)
  write(6,*) '  Which file has the crystal structure?'
  write(6,*)
  write(6,*) '  1) cpw.in'
  write(6,*) '  2) cpw.out (optimized structure after a relaxation)'
  write(6,*) '  3) enter filename'
  write(6,*)

  read(5,*) ichoice
  write(ioreplay,*) ichoice,'   file choice'

  if(ichoice == 1) then
    fcpwin = 'cpw.in'
  elseif(ichoice == 2) then
    fcpwin = 'cpw.out'
  elseif(ichoice == 3) then
    write(6,*)
    write(6,*) '  Enter file name'
    write(6,*)
    read(5,'(a)') fcpwin
    fcpwin = adjustl(fcpwin)
    iblank = index(fcpwin,' ')
    if(iblank > 0) fcpwin(iblank:) = ' '
    write(ioreplay,'(1x,a,"   filename")') trim(fcpwin)
  else
    fcpwin = 'cpw.in'
    write(6,*)
    write(6,*) '  Wrong answer, using cpw.in'
    write(6,*)
  endif

  inquire(file = trim(fcpwin), exist = lex)

  if(.not. lex) then
    write(6,*)
    write(6,*) '  STOPPED in cpwin_geom:  file ',trim(fcpwin),' does not exist'
    write(6,*)

    stop

  endif

  call size_mxdtyp_mxdatm_esdf(ipr, trim(fcpwin), lgeom,                 &
              mxdtyp, mxdatm)

  if(lgeom) then

    allocate(natom(mxdtyp))
    allocate(rat(3,mxdatm,mxdtyp))
    allocate(atmass(mxdtyp))
    allocate(nameat(mxdtyp))

    call esdf_init(trim(fcpwin))

    call read_esdf_crystal(ipr,                                          &
        adot, ntype, natom, nameat, rat,                                 &
        atmass, alatt, lgeom2,                                           &
        mxdtyp, mxdatm)

    call esdf_close

    if(lgeom2) then

      call voronoi_analysis(ioreplay,                                    &
          adot, ntype, natom, nameat, rat,                               &
          mxdtyp, mxdatm)

    else
      write(6,*)
      write(6,*) '  could not get enough information from ',trim(fcpwin)
      write(6,*) '  to figure out the geometry'
      write(6,*)
    endif

    deallocate(natom)
    deallocate(rat)
    deallocate(atmass)
    deallocate(nameat)

  else
    write(6,*)
    write(6,*) '  could not get enough information from ',trim(fcpwin)
    write(6,*) '  to find mxdtyp, mxdatm '
    write(6,*)
  endif

  close(unit = ioreplay)

  stop
end program cpwin_geom
