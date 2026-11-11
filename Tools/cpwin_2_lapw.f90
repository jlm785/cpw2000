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

!>  Reads the crystal structure, exchange-correlation functional and
!>  k-point grid from a cpw.in (or cpw.out, or other esdf file)
!>  and writes input files for the all-electron LAPW codes
!>  Elk (elk.in) and Fleur (fleur_inpgen.in, input of inpgen).
!>
!>  Muffin-tin radii, basis cutoffs and local orbitals are left
!>  to the defaults of each code.  For a precise comparison they
!>  should be set (and converged) by hand.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         11 November 2026.
!>  \copyright    GNU Public License v2

program cpwin_2_lapw

! Written 11 November 2026. claude

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

! other input variables

  character(len=4)                   ::  author                          !<  type of xc wanted (CA, PBE, R2SC,...)
  real(REAL64)                       ::  emax                            !<  kinetic energy cutoff of plane wave expansion (Hartree).
  integer                            ::  nbandin                         !<  target for number of bands
  integer                            ::  nx, ny, nz                      !<  size of the integration mesh in k-space
  real(REAL64)                       ::  sx, sy, sz                      !<  offset of the integration mesh (usually 0.5)

! lattice vectors

  real(REAL64)                       ::  avec(3,3)                       !<  primitive lattice vectors (columns), adot = transpose(avec)*avec
  real(REAL64)                       ::  bvec(3,3)                       !<  reciprocal lattice vectors

! xc in the LAPW codes

  character(len=40)                  ::  elkxc                           !<  xctype of Elk
  character(len=200)                 ::  fleurxc                         !<  xctyp of Fleur inpgen
  logical                            ::  lmgga                           !<  meta-GGA
  logical                            ::  lxcok                           !<  functional known

! other variables

  integer               ::  ipr
  character(len=200)    ::  fcpwin                                       !  filename with input data (new style)
  character(len=200)    ::  sppath                                       !  directory of the Elk species files
  logical               ::  lgeom, lgeom2, lbz
  logical               ::  lex                                          !  file exists
  logical               ::  lnonrel                                      !  non-relativistic calculation
  integer               ::  ichoice, iblank
  character(len=1)      ::  yesno
  integer               ::  ioreplay = 35                                !  hard coded tape unit for replay
  character(len=16)     ::  fname = 'replay_lapw.dat'                    !  hard coded filename for replay

  integer               ::  ioelk = 36                                   !  tape unit for elk.in
  integer               ::  iofleur = 37                                 !  tape unit for fleur_inpgen.in

  ipr = 1

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
    write(6,*) '  STOPPED in cpwin_2_lapw:  file ',trim(fcpwin),' does not exist'
    write(6,*)

    stop

  endif

  call size_mxdtyp_mxdatm_esdf(ipr, trim(fcpwin), lgeom,                 &
              mxdtyp, mxdatm)

  if(.not. lgeom) then
    write(6,*)
    write(6,*) '  STOPPED in cpwin_2_lapw:  could not get the number of atoms'
    write(6,*) '  from ',trim(fcpwin)
    write(6,*)

    stop

  endif

  allocate(natom(mxdtyp))
  allocate(rat(3,mxdatm,mxdtyp))
  allocate(atmass(mxdtyp))
  allocate(nameat(mxdtyp))

  call esdf_init(trim(fcpwin))

  call read_esdf_crystal(ipr,                                            &
      adot, ntype, natom, nameat, rat,                                   &
      atmass, alatt, lgeom2,                                             &
      mxdtyp, mxdatm)

  call read_esdf_bz(ipr,                                                 &
      emax, nbandin, nx,ny,nz, sx,sy,sz, lbz)

  author = 'CA  '
  author = esdf_string('XC.Authors', author)

  call esdf_close

  if(.not. lgeom2) then
    write(6,*)
    write(6,*) '  STOPPED in cpwin_2_lapw:  could not get the geometry'
    write(6,*) '  from ',trim(fcpwin)
    write(6,*)

    stop

  endif

  if(.not. lbz) then
    nx = 4
    ny = 4
    nz = 4
    sx = 0.5_REAL64
    sy = 0.5_REAL64
    sz = 0.5_REAL64
    write(6,*)
    write(6,*) '  WARNING:  no k-point grid in ',trim(fcpwin),', using 4x4x4 shifted'
    write(6,*)
  endif

! usual lattice vectors (cartesian, bohr) from the metric of cpw2000,
! adot is already rounded to rationals or square roots of rationals
! by read_esdf_crystal, so the symmetry is exact

  call adot_to_avec_sym(adot, avec, bvec)

! functional

  call cpwin_2_lapw_xc(author, elkxc, fleurxc, lmgga, lxcok)

  if(.not. lxcok) then
    write(6,*)
    write(6,*) '  WARNING:  XC.Authors = ',author,' has no equivalent'
    write(6,*) '  in the LAPW codes, edit the xc lines by hand'
    write(6,*)
  endif

! questions

  write(6,*)
  write(6,*) '  Do you want a non-relativistic calculation (y/n)?'
  write(6,*) '  (cpw2000 pseudopotentials generated with "nrl")'
  write(6,*)

  read(5,*) yesno
  write(ioreplay,*) yesno,'   non-relativistic'
  lnonrel = (yesno == 'y' .or. yesno == 'Y')

  write(6,*)
  write(6,*) '  Enter the directory of the Elk species files'
  write(6,*) '  (e.g. ~/elk/species/,  "." for the current directory)'
  write(6,*)

  read(5,'(a)') sppath
  sppath = adjustl(sppath)
  if(len_trim(sppath) == 0) sppath = './'
  if(sppath(len_trim(sppath):len_trim(sppath)) /= '/') sppath = trim(sppath)//'/'
  write(ioreplay,'(1x,a,"   Elk species directory")') trim(sppath)

  close(unit = ioreplay)

! writes the files

  open(unit = ioelk, file = 'elk.in', status = 'UNKNOWN', form = 'FORMATTED')

  call cpwin_2_lapw_elk(ioelk, trim(fcpwin), author, elkxc, lmgga,       &
      lnonrel, trim(sppath),                                             &
      avec, ntype, natom, nameat, rat,                                   &
      nx, ny, nz, sx, sy, sz,                                            &
      mxdtyp, mxdatm)

  close(unit = ioelk)

  open(unit = iofleur, file = 'fleur_inpgen.in', status = 'UNKNOWN',     &
      form = 'FORMATTED')

  call cpwin_2_lapw_fleur(iofleur, trim(fcpwin), author, fleurxc,        &
      avec, ntype, natom, nameat, rat,                                   &
      nx, ny, nz, sx, sy, sz,                                            &
      mxdtyp, mxdatm)

  close(unit = iofleur)

! notes

  write(6,*)
  write(6,*) '  Wrote elk.in and fleur_inpgen.in   (XC.Authors = ',author,')'
  write(6,*)
  write(6,*) '  Elk:    copy or link the species files, choose the quality'
  write(6,*) '          (e.g. vhighq; it sets autokpt, keep ngridk after it)'
  write(6,*) '          and set rmtall to the radius used in Fleur.'
  write(6,*) '  Fleur:  inpgen -f fleur_inpgen.in   then edit inp.xml:'
  write(6,*) '          muffin-tin radii and cutoffs (inpgen defaults),'
  write(6,*) '          check kpts.xml against the cpw2000 k-points'
  if(lnonrel) then
    write(6,*) '          non-relativistic:  <expertModes ... warp_factor="10000.0"/>'
    write(6,*) '          (Fleur develop branch)'
  endif
  if(lmgga) then
    write(6,*)
    write(6,*) '  meta-GGA:  Fleur needs the develop branch (generalized Kohn-Sham).'
    write(6,*) '             Elk replaces the orbital tau by a kinetic energy'
    write(6,*) '             functional (partial deorbitalisation), use with care.'
  endif
  write(6,*)

  deallocate(natom)
  deallocate(rat)
  deallocate(atmass)
  deallocate(nameat)

end program cpwin_2_lapw




!>  Equivalent of the cpw2000 XC.Authors in Elk (xctype)
!>  and in the Fleur input generator (xctyp)
!>
!>  \author       claude
!>  \version      5.13
!>  \date         11 November 2026.
!>  \copyright    GNU Public License v2

subroutine cpwin_2_lapw_xc(author, elkxc, fleurxc, lmgga, lxcok)

! Written 11 November 2026. claude

  implicit none

! input

  character(len=4), intent(in)       ::  author                          !<  type of xc wanted (CA, PBE, R2SC,...)

! output

  character(len=40), intent(out)     ::  elkxc                           !<  xctype of Elk
  character(len=200), intent(out)    ::  fleurxc                         !<  xctyp of Fleur inpgen
  logical, intent(out)               ::  lmgga                           !<  meta-GGA
  logical, intent(out)               ::  lxcok                           !<  functional known

! local

  character(len=4)                   ::  xcbase                          !  meta-GGA functional used in xc_mgga
  character(len=4)                   ::  tausrc                          !  source of tau
  character(len=*), parameter        ::  AUX = ', AuxExchID: 117, AuxCorID: 130'

! functions

  logical                            ::  chrsameinfo                     !  strings are the same irrespective of case or blanks


  lmgga = .FALSE.
  lxcok = .TRUE.

  call xc_author_tau(author, xcbase, tausrc)

  if(tausrc /= 'NONE' .and. tausrc /= 'PSI ') then
    write(6,*)
    write(6,*) '  WARNING:  ',author,' is a deorbitalised version of ',xcbase
    write(6,*) '  the LAPW inputs use ',xcbase,' with the orbital tau'
    write(6,*)
  endif

! the Fleur meta-GGA inputs use an auxiliary GGA (libxc 117, 130)
! for the radial basis functions, as testing/inputfiles/extra/MetaGGA_SCAN

  if(chrsameinfo(author, 'CA') .or. chrsameinfo(author, 'PZ')) then
    elkxc = '2'
    fleurxc = 'pz'
  elseif(chrsameinfo(author, 'PW92')) then
    elkxc = '3'
    fleurxc = 'LibXC: Exch: lda_x, Cor: lda_c_pw'
  elseif(chrsameinfo(author, 'VWN')) then
    elkxc = '100 1 7'
    fleurxc = 'LibXC: Exch: lda_x, Cor: lda_c_vwn'
  elseif(chrsameinfo(author, 'WI')) then
    elkxc = '100 1 2'
    fleurxc = 'LibXC: Exch: lda_x, Cor: lda_c_wigner'
  elseif(chrsameinfo(author, 'PBE')) then
    elkxc = '20'
    fleurxc = 'pbe'
  elseif(chrsameinfo(xcbase, 'R2SC')) then
    lmgga = .TRUE.
    elkxc = '100 497 498'
    fleurxc = 'LibXC: Exch: mgga_x_r2scan, Cor: mgga_c_r2scan'//AUX
  elseif(chrsameinfo(xcbase, 'TASK')) then
    lmgga = .TRUE.
    elkxc = '100 707 12'
    fleurxc = 'LibXC: Exch: mgga_x_task, Cor: lda_c_pw'//AUX
  elseif(chrsameinfo(xcbase, 'LAK')) then
    lmgga = .TRUE.
    elkxc = '100 342 345'
    fleurxc = 'LibXC: Exch: mgga_x_lak, Cor: mgga_c_lak'//AUX
  elseif(chrsameinfo(author, 'TBL') .or. chrsameinfo(author, 'TB09')) then
    lmgga = .TRUE.
    elkxc = '100 208 12'
    fleurxc = 'LibXC: Exch: mgga_x_tb09, Cor: lda_c_pw'//AUX
    write(6,*)
    write(6,*) '  WARNING:  TB09 is a potential-only functional'
    write(6,*)
  else
    lxcok = .FALSE.
    elkxc = '2'
    fleurxc = 'pz'
  endif

  return

end subroutine cpwin_2_lapw_xc




!>  Writes an Elk input file (elk.in) for a ground state calculation
!>
!>  \author       claude
!>  \version      5.13
!>  \date         11 November 2026.
!>  \copyright    GNU Public License v2

subroutine cpwin_2_lapw_elk(io, fcpwin, author, elkxc, lmgga,            &
    lnonrel, sppath,                                                     &
    avec, ntype, natom, nameat, rat,                                     &
    nx, ny, nz, sx, sy, sz,                                              &
    mxdtyp, mxdatm)

! Written 11 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type

  integer, intent(in)                ::  io                              !<  tape number of elk.in
  character(len=*), intent(in)       ::  fcpwin                          !<  name of the cpw2000 input file
  character(len=4), intent(in)       ::  author                          !<  type of xc in cpw2000
  character(len=40), intent(in)      ::  elkxc                           !<  xctype of Elk
  logical, intent(in)                ::  lmgga                           !<  meta-GGA
  logical, intent(in)                ::  lnonrel                         !<  non-relativistic calculation
  character(len=*), intent(in)       ::  sppath                          !<  directory of the Elk species files

  real(REAL64), intent(in)           ::  avec(3,3)                       !<  primitive lattice vectors (columns), bohr
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i

  integer, intent(in)                ::  nx, ny, nz                      !<  size of the integration mesh in k-space
  real(REAL64), intent(in)           ::  sx, sy, sz                      !<  offset of the integration mesh

! local

  integer               ::  nt, ja, j
  character(len=40)     ::  spfile

  write(io,'("! from ",a," of cpw2000 by cpwin_2_lapw, XC.Authors = ",a)') fcpwin, author
  write(io,'("! set the quality (e.g. vhighq .true.) and rmtall by hand")')
  write(io,*)
  write(io,'("tasks")')
  write(io,'("  0")')
  write(io,*)
  write(io,'("xctype")')
  write(io,'("  ",a)') trim(elkxc)
  write(io,*)
  if(lmgga) then
    write(io,'("metagga")')
    write(io,'("  .true.")')
    write(io,*)
    write(io,'("mixtype")')
    write(io,'("  1")')
    write(io,*)
  endif
  if(lnonrel) then
    write(io,'("! non-relativistic (scaled speed of light)")')
    write(io,'("solscf")')
    write(io,'("  1.0e4")')
    write(io,*)
  endif
  write(io,'("! k-point grid of cpw2000, after any quality block (autokpt)")')
  write(io,'("autokpt")')
  write(io,'("  .false.")')
  write(io,*)
  write(io,'("ngridk")')
  write(io,'(3i5)') nx, ny, nz
  write(io,*)
  write(io,'("vkloff")')
  write(io,'(3f8.4)') sx, sy, sz
  write(io,*)
  write(io,'("avec")')
  do j = 1,3
    write(io,'(3f20.12)') avec(1,j), avec(2,j), avec(3,j)
  enddo
  write(io,*)
  write(io,'("scale")')
  write(io,'("  1.0")')
  write(io,*)
  write(io,'("sppath")')
  write(io,'("  ''",a,"''")') sppath
  write(io,*)
  write(io,'("atoms")')
  write(io,'(i4,32x,": nspecies")') ntype
  do nt = 1,ntype
    spfile = "'"//trim(adjustl(nameat(nt)))//".in'"
    write(io,'(2x,a,t38,": spfname")') trim(spfile)
    write(io,'(i4,32x,": natoms; atposl below")') natom(nt)
    do ja = 1,natom(nt)
      write(io,'(3f18.12)') rat(1,ja,nt), rat(2,ja,nt), rat(3,ja,nt)
    enddo
  enddo

  return

end subroutine cpwin_2_lapw_elk




!>  Writes an input file for the Fleur input generator (inpgen)
!>
!>  \author       claude
!>  \version      5.13
!>  \date         11 November 2026.
!>  \copyright    GNU Public License v2

subroutine cpwin_2_lapw_fleur(io, fcpwin, author, fleurxc,               &
    avec, ntype, natom, nameat, rat,                                     &
    nx, ny, nz, sx, sy, sz,                                              &
    mxdtyp, mxdatm)

! Written 11 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type

  integer, intent(in)                ::  io                              !<  tape number of the inpgen file
  character(len=*), intent(in)       ::  fcpwin                          !<  name of the cpw2000 input file
  character(len=4), intent(in)       ::  author                          !<  type of xc in cpw2000
  character(len=200), intent(in)     ::  fleurxc                         !<  xctyp of Fleur inpgen

  real(REAL64), intent(in)           ::  avec(3,3)                       !<  primitive lattice vectors (columns), bohr
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i

  integer, intent(in)                ::  nx, ny, nz                      !<  size of the integration mesh in k-space
  real(REAL64), intent(in)           ::  sx, sy, sz                      !<  offset of the integration mesh

! local

  integer               ::  nt, ja, j, iz, ntot
  real(REAL64)          ::  zid
  logical               ::  lsame                                        !  same element in several types
  integer               ::  nt2

  real(REAL64), parameter    ::  EPS = 1.0E-6_REAL64

! types with the same element need an identifier Z.id

  lsame = .FALSE.
  do nt = 1,ntype
    do nt2 = nt+1,ntype
      if(nameat(nt) == nameat(nt2)) lsame = .TRUE.
    enddo
  enddo
  if(lsame .and. ntype > 9) then
    write(6,*)
    write(6,*) '  WARNING:  more than 9 types, check the Z.id labels in fleur_inpgen.in'
    write(6,*)
  endif

  ntot = 0
  do nt = 1,ntype
    ntot = ntot + natom(nt)
  enddo

  write(io,'(a," of cpw2000 (cpwin_2_lapw), XC.Authors = ",a)') fcpwin, author
  write(io,'("&input cartesian=F /")')
  do j = 1,3
    write(io,'(3f20.12)') avec(1,j), avec(2,j), avec(3,j)
  enddo
  write(io,'("  1.0")')
  write(io,'("  1.0  1.0  1.0")')
  write(io,*)
  write(io,'(i5)') ntot
  do nt = 1,ntype
    call p_tbl_charge(nameat(nt), iz)
    zid = real(iz,REAL64)
    if(lsame) zid = zid + real(nt,REAL64) / 10
    do ja = 1,natom(nt)
      write(io,'(f8.1,3f18.12)') zid, rat(1,ja,nt), rat(2,ja,nt), rat(3,ja,nt)
    enddo
  enddo
  write(io,'("&exco xctyp=''",a,"'' /")') trim(fleurxc)

! Fleur grid with gamma=T contains Gamma, otherwise Monkhorst-Pack

  if(abs(sx) < EPS .and. abs(sy) < EPS .and. abs(sz) < EPS) then
    write(io,'("&kpt div1=",i0," div2=",i0," div3=",i0," gamma=T /")') nx, ny, nz
  else
    write(io,'("&kpt div1=",i0," div2=",i0," div3=",i0," /")') nx, ny, nz
    if(abs(sx-0.5) > EPS .or. abs(sy-0.5) > EPS .or. abs(sz-0.5) > EPS) then
      write(6,*)
      write(6,*) '  WARNING:  k-grid offset is not 0 or 1/2, check kpts.xml'
      write(6,*)
    endif
  endif

  return

end subroutine cpwin_2_lapw_fleur
