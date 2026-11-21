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

!>  Filter that reads cpw.in (or another esdf file), chooses a general
!>  Monkhorst-Pack mesh and writes cpw.out, a verbatim copy of the input
!>  where only the kgrid_Monkhorst_Pack block is replaced (the old block
!>  is kept as comments).
!>
!>
!>  \author       claude
!>  \version      5.13
!>  \date         20 November 2026.
!>  \copyright    GNU Public License v2

program cpwin_kmesh

! Written 20 November 2026. claude

! The mesh with supercell  A_i = sum_j mpgrid(i,j) a_j  and shift s
! averages exp(i k.R) to  w(R) = cos(2 pi s.m)  if  R = sum_i m_i A_i,
! and to zero otherwise.  The integration error of a periodic function
! f(k) = sum_R f_R exp(i k.R)  is  sum_{R/=0} f_R w(R).
!
! Only meshes invariant under the point group of the parent lattice
! (the lattice of highest symmetry within a strain delta) are
! considered, so that a change of symmetry during a relaxation does not
! change the effective mesh.  The quality of a mesh is R_cc, the length
! of the shortest lattice vectors (a star of the parent group) whose sum
! of w(R) does not vanish (Chadi and Cohen).  For an unshifted mesh it is
! the length of the shortest supercell vector, twice the kgrid_cutoff of
! siesta (Moreno and Soler).
!
! The chosen mesh has R_cc >= target and the smallest number of
! irreducible k-points for the symmetry of the crystal.  If the
! Rede.Superlattice block is present (with Rede.NumberOfLatticePlanes)
! the block with the same k-points in the primitive cell of the
! superlattice is also given, for the bulk reference calculations.
!
! usage:   cpwin_kmesh.exe  [file  [R_cc target in bohr  [delta]]]
!
! defaults: file = cpw.in, the target is asked, delta = 0.05.




  use esdf

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! dimensions

  integer                            ::  mxdtyp                          !<  array dimension of types of atoms
  integer                            ::  mxdatm                          !<  array dimension of number of atoms of a given type

! crystal

  real(REAL64)                       ::  alatt                           !<  lattice constant
  real(REAL64)                       ::  adot(3,3)                       !<  metric in direct space in atomic units (Bohr radius)
  integer                            ::  ntype                           !<  number of types of atoms
  integer, allocatable               ::  natom(:)                        !<  number of atoms of type i
  real(REAL64), allocatable          ::  rat(:,:,:)                      !<  lattice coordinates of atom j of type i
  real(REAL64), allocatable          ::  atmass(:)                       !<  atomic mass of atoms of type i
  character(len=2), allocatable      ::  nameat(:)                       !<  chemical symbol for the type i

! Brillouin zone

  real(REAL64)                       ::  emax                            !<  kinetic energy cutoff of plane wave expansion (Hartree).
  integer                            ::  nbandin                         !<  target for number of bands
  integer                            ::  mpgrid0(3,3)                    !<  present block, A_i = sum_j mpgrid(i,j) a_j
  real(REAL64)                       ::  mpshift0(3)                     !<  present shift
  logical                            ::  lbz                             !<  kgrid_Monkhorst_Pack was read

! superlattice

  character(len=250)                 ::  meta_cpw2000                    !<  metadata of rede
  logical                            ::  lrede                           !<  Rede metadata present
  integer                            ::  isl(3,3)                        !<  isl(:,i) line i of Rede.Superlattice
  integer                            ::  islr(3,3)                       !<  line i of Rede.Superlattice, supercell vector i in primitive cell coordinates
  real(REAL64)                       ::  adotpc(3,3)                     !<  metric of the primitive cell of the superlattice

! symmetry

  real(REAL64)                       ::  symtol                          !<  tolerance of the crystal symmetry
  integer                            ::  ntrans                          !<  number of operations of the crystal point group
  integer                            ::  mtrx(3,3,48)                    !<  operations in reciprocal lattice coordinates
  real(REAL64)                       ::  tnp(3,48)                       !<  fractional translations
  integer                            ::  geff                            !<  largest orbit of a k-point (with time reversal)
  integer                            ::  nh                              !<  number of operations of the parent lattice
  integer                            ::  dh(3,3,48)                      !<  operations of the parent lattice, direct lattice coordinates
  real(REAL64)                       ::  strain                          !<  largest strain to the parent lattice
  integer                            ::  nhb                             !<  number of operations of the parent of the primitive cell
  integer                            ::  dhb(3,3,48)                     !<  operations of the parent of the primitive cell
  real(REAL64)                       ::  strainb                         !<  largest strain to that parent

! present mesh

  logical                            ::  linv0, lfound0, linvb0
  real(REAL64)                       ::  mpshc0(3), mpshc1(3), mpshcb(3)
  real(REAL64)                       ::  rcc0, sumw0, rni0, rmax0
  integer                            ::  nk0, nirr0
  integer                            ::  mpgridr(3,3), mpc0(3,3)
  real(REAL64)                       ::  mpshiftr(3)

! new mesh

  real(REAL64)                       ::  rtarget                         !<  target R_cc (bohr)
  real(REAL64)                       ::  delta                           !<  largest strain to the parent lattice
  logical                            ::  lfound1, lbulk1
  integer                            ::  mpgrid1(3,3), nk1, nirr1, mpc1(3,3)
  real(REAL64)                       ::  mpshift1(3), rcc1
  logical                            ::  lfoundb
  integer                            ::  mpgridb(3,3), nkb, nirrb, mpcb(3,3)
  real(REAL64)                       ::  mpshiftb(3), rccb

! output comments

  integer, parameter                 ::  MXCOM = 30
  integer                            ::  ncom
  character(len=120)                 ::  com(MXCOM)

! other variables

  integer               ::  ipr, nargs, ioerr, idet
  character(len=200)    ::  fcpwin, carg
  logical               ::  lgeom, lgeom2, lex, linvc
  real(REAL64)          ::  adotc(3,3)
  real(REAL64), allocatable     ::  ratc(:,:,:)
  real(REAL64)          ::  slinv(3,3), xdet
  integer               ::  dc(3,3), adjm(3,3)
  character(len=40)     ::  family

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counters

  integer   ::  i, j, k, l, n


  ipr = 1

! arguments

  fcpwin = 'cpw.in'
  rtarget = -1.0_REAL64
  delta = 0.05_REAL64

  nargs = command_argument_count()
  if(nargs >= 1) call get_command_argument(1, fcpwin)
  if(nargs >= 2) then
    call get_command_argument(2, carg)
    read(carg,*,iostat=ioerr) rtarget
    if(ioerr /= 0) rtarget = -1.0_REAL64
  endif
  if(nargs >= 3) then
    call get_command_argument(3, carg)
    read(carg,*,iostat=ioerr) delta
    if(ioerr /= 0 .or. delta <= ZERO) delta = 0.05_REAL64
  endif

  inquire(file = trim(fcpwin), exist = lex)

  if(.not. lex) then
    write(6,*)
    write(6,*) '  STOPPED in cpwin_kmesh:  file ',trim(fcpwin),' does not exist'
    write(6,*)

    stop

  endif

! reads the file

  call size_mxdtyp_mxdatm_esdf(ipr, trim(fcpwin), lgeom,                 &
              mxdtyp, mxdatm)

  if(.not. lgeom) then
    write(6,*)
    write(6,*) '  STOPPED in cpwin_kmesh:  could not get the number of atoms'
    write(6,*) '  from ',trim(fcpwin)
    write(6,*)

    stop

  endif

  allocate(natom(mxdtyp))
  allocate(rat(3,mxdatm,mxdtyp))
  allocate(ratc(3,mxdatm,mxdtyp))
  allocate(atmass(mxdtyp))
  allocate(nameat(mxdtyp))

  call esdf_init(trim(fcpwin))

  call read_esdf_crystal(ipr,                                            &
      adot, ntype, natom, nameat, rat,                                   &
      atmass, alatt, lgeom2,                                             &
      mxdtyp, mxdatm)

  call read_esdf_bz(ipr,                                                 &
      emax, nbandin, mpgrid0, mpshift0, lbz)

  call read_esdf_rede(ipr, meta_cpw2000, lrede)

  symtol = 1.0E-5_REAL64
  symtol = esdf_double('SymmTolerance', symtol)

  call esdf_close

  if(.not. lgeom2) then
    write(6,*)
    write(6,*) '  STOPPED in cpwin_kmesh:  could not get the geometry'
    write(6,*) '  from ',trim(fcpwin)
    write(6,*)

    stop

  endif

! superlattice matrix, same order as read_esdf_rede

  if(lrede) then
    read(meta_cpw2000,'(3(3i4,2x))',iostat=ioerr)                         &
          (isl(i,1),i=1,3), (isl(i,2),i=1,3), (isl(i,3),i=1,3)
    if(ioerr /= 0) lrede = .FALSE.
    do i = 1,3
    do j = 1,3
      islr(i,j) = isl(j,i)
    enddo
    enddo
    call cpwin_kmesh_adj3(islr, adjm, idet)
    if(idet == 0) then
      lrede = .FALSE.
      write(6,*)
      write(6,*) '  WARNING:  the Rede.Superlattice matrix is singular, ignored'
      write(6,*)
    endif
  endif

  write(6,*)
  write(6,'("  cpwin_kmesh:  ",a)') trim(fcpwin)
  write(6,*)

! symmetry of the crystal

  do j = 1,3
  do i = 1,3
    adotc(i,j) = adot(i,j)
  enddo
  enddo
  do k = 1,mxdtyp
  do j = 1,mxdatm
  do i = 1,3
    ratc(i,j,k) = rat(i,j,k)
  enddo
  enddo
  enddo

  call sym_identify(1, 0, symtol,                                        &
      ntrans, mtrx, tnp,                                                 &
      adotc, ntype, natom, ratc,                                         &
      mxdtyp, mxdatm)

  linvc = .FALSE.
  do n = 1,ntrans
    if(mtrx(1,1,n) == -1 .and. mtrx(2,2,n) == -1 .and. mtrx(3,3,n) == -1 .and.  &
       mtrx(1,2,n) == 0 .and. mtrx(1,3,n) == 0 .and. mtrx(2,1,n) == 0 .and.     &
       mtrx(2,3,n) == 0 .and. mtrx(3,1,n) == 0 .and. mtrx(3,2,n) == 0) linvc = .TRUE.
  enddo
  geff = ntrans
  if(.not. linvc) geff = 2*ntrans

! parent lattice

  call cpwin_kmesh_parent(adot, delta, nh, dh, strain)

  call cpwin_kmesh_family(nh, family)

  write(6,'("  Crystal point group:  ",i3," operations (symmetry tolerance",    &
        & es9.2,")")') ntrans, symtol
  write(6,'("  Parent lattice:       ",a,", ",i3," operations,",                &
        & " largest strain",f8.4," (delta =",f6.3,")")')                         &
        trim(family), nh, strain, delta

! the crystal operations must be operations of the parent lattice,
! direct lattice coordinates  D = (mtrx^-1)^T

  do n = 1,ntrans
    call cpwin_kmesh_adj3(mtrx(:,:,n), adjm, idet)
    do i = 1,3
    do j = 1,3
      dc(i,j) = adjm(j,i)*idet
    enddo
    enddo
    lex = .FALSE.
    do k = 1,nh
      l = 0
      do i = 1,3
      do j = 1,3
        if(dc(i,j) /= dh(i,j,k)) l = 1
      enddo
      enddo
      if(l == 0) lex = .TRUE.
    enddo
    if(.not. lex) then
      write(6,*)
      write(6,*) '  WARNING:  an operation of the crystal is not an operation',  &
                 ' of the parent lattice'
      write(6,*) '  (delta smaller than the symmetry tolerance?)'
      write(6,*)
      exit
    endif
  enddo

! parent lattice of the primitive cell of the superlattice

  adotpc = adot
  nhb = 0
  if(lrede) then
    call cpwin_kmesh_inv3(real(islr,REAL64), slinv, xdet)
    do i = 1,3
    do j = 1,3
      adotpc(i,j) = ZERO
      do k = 1,3
      do l = 1,3
        adotpc(i,j) = adotpc(i,j) + slinv(i,k)*adot(k,l)*slinv(j,l)
      enddo
      enddo
    enddo
    enddo
    call cpwin_kmesh_parent(adotpc, delta, nhb, dhb, strainb)
    call cpwin_kmesh_family(nhb, family)
    write(6,'("  Rede.Superlattice:    ",3(3i4,3x))') ((islr(i,j),j=1,3),i=1,3)
    write(6,'("  Primitive cell:       ",a,", ",i3," operations,",              &
        & " largest strain",f8.4)') trim(family), nhb, strainb
  endif
  write(6,*)

! present mesh

  rmax0 = ZERO
  if(lbz) then

    call cpwin_kmesh_reduce(adot, mpgrid0, mpshift0, mpgridr, mpshiftr)
    rmax0 = 1.0E10_REAL64
    do i = 1,3
      xdet = ZERO
      do k = 1,3
      do l = 1,3
        xdet = xdet + mpgridr(i,k)*adot(k,l)*mpgridr(i,l)
      enddo
      enddo
      rmax0 = min(rmax0,sqrt(xdet))
    enddo
    rmax0 = 2.5_REAL64*rmax0
    if(rtarget > ZERO) rmax0 = max(rmax0, 1.5_REAL64*rtarget)

    call cpwin_kmesh_rcc(adot, mpgrid0, mpshift0, nh, dh, rmax0,          &
        linv0, lfound0, rcc0, sumw0, rni0)

    call cpwin_kmesh_adj3(mpgrid0, adjm, nk0)
    nk0 = abs(nk0)

    call cpwin_kmesh_nirr(adot, mpgrid0, mpshift0, ntrans, mtrx, nirr0)

    write(6,*) '  Present kgrid_Monkhorst_Pack:'
    do i = 1,3
      write(6,'(8x,3(2x,i6),3x,f12.6)') (mpgrid0(i,j),j=1,3), mpshift0(i)
    enddo
    if(lfound0) then
      write(6,'(5x,i6," k-points,",i5," irreducible,  R_cc =",f9.2," bohr")') &
          nk0, nirr0, rcc0
    else
      write(6,'(5x,i6," k-points,",i5," irreducible,  R_cc >",f9.2," bohr")') &
          nk0, nirr0, rcc0
    endif
    if(linv0) then
      write(6,'(5x,"invariant under the parent lattice")')
    else
      if(rni0 < rmax0) then
        write(6,'(5x,"NOT invariant under the parent lattice, first star",   &
          & " with different weights at",f9.2," bohr")') rni0
      else
        write(6,'(5x,"NOT invariant under the parent lattice")')
      endif
    endif
    if(lrede) then
      call cpwin_kmesh_bulk(mpgrid0, mpshift0, islr, adotpc, nhb, dhb,  &
          mpc0, mpshc0, linvb0)
      write(6,'(5x,"same k-points in the primitive cell:",3(3i4,3x),3f6.2)') &
          ((mpc0(i,j),j=1,3),i=1,3), (mpshc0(i),i=1,3)
      if(linvb0) then
        write(6,'(5x,"(invariant under the parent of the primitive cell)")')
      else
        write(6,'(5x,"(NOT invariant under the parent of the primitive cell)")')
      endif
    endif
    write(6,*)

  else

    write(6,*) '  No kgrid_Monkhorst_Pack block in ',trim(fcpwin)
    write(6,*)

  endif

! target

  if(rtarget <= ZERO) then
    write(6,*) '  Enter the target R_cc in bohr'
    if(lbz) write(6,'("  (present mesh R_cc =",f9.2,")")') rcc0
    read(5,*) rtarget
    if(rtarget <= ZERO) then
      write(6,*)
      write(6,*) '  STOPPED in cpwin_kmesh:  the target must be positive'
      write(6,*)

      stop

    endif
  endif

  write(6,'("  Target R_cc =",f9.2," bohr.  Largest R_cc for each number",   &
     & " of irreducible k-points:")') rtarget
  write(6,*)

  call cpwin_kmesh_search(adot, nh, dh, ntrans, mtrx, geff,              &
      lrede, islr, adotpc, nhb, dhb, rtarget,                            &
      lfound1, mpgrid1, mpshift1, nk1, nirr1, rcc1, lbulk1,              &
      lfoundb, mpgridb, mpshiftb, nkb, nirrb, rccb)

  if(.not. lfound1) then
    write(6,*)
    write(6,*) '  STOPPED in cpwin_kmesh:  no mesh found that reaches the target'
    write(6,*) '  cpw.out not written'
    write(6,*)

    stop

  endif

  write(6,*)
  write(6,*) '  Chosen mesh:'
  write(6,*)
  write(6,'("%block kgrid_Monkhorst_Pack")')
  do i = 1,3
    write(6,'(8x,3(2x,i6),3x,f12.6)')  (mpgrid1(i,j),j=1,3), mpshift1(i)
  enddo
  write(6,'("%endblock kgrid_Monkhorst_Pack")')
  write(6,*)
  write(6,'(5x,i6," k-points,",i5," irreducible,  R_cc =",f9.2," bohr")')   &
      nk1, nirr1, rcc1

! comments for cpw.out

  ncom = 0
  ncom = ncom + 1
  write(com(ncom),'("# kgrid_Monkhorst_Pack chosen by cpwin_kmesh,",        &
       & " target R_cc =",f9.2," bohr")') rtarget
  ncom = ncom + 1
  write(com(ncom),'("#   R_cc =",f9.2," bohr,",i6," k-points,",i5,           &
       & " irreducible")') rcc1, nk1, nirr1
  ncom = ncom + 1
  call cpwin_kmesh_family(nh, family)
  write(com(ncom),'("#   invariant under the ",i2," operations of the",      &
       & " parent lattice (",a,", strain",f7.4,")")') nh, trim(family), strain
  ncom = ncom + 1
  com(ncom) = '#   R_cc: length of the shortest lattice vectors whose mesh average does not cancel'

  if(lrede) then
    call cpwin_kmesh_bulk(mpgrid1, mpshift1, islr, adotpc, nhb, dhb,    &
        mpc1, mpshc1, lbulk1)
    write(6,*)
    write(6,*) '  The same k-points in the primitive cell of Rede.Superlattice:'
    write(6,*)
    do i = 1,3
      write(6,'(8x,3(2x,i6),3x,f12.6)')  (mpc1(i,j),j=1,3), mpshc1(i)
    enddo
    if(lbulk1) then
      write(6,'(5x,"(invariant under the parent of the primitive cell)")')
    else
      write(6,'(5x,"(NOT invariant under the parent of the primitive cell:",  &
          & " for the same k-points the bulk")')
      write(6,'(5x," calculation needs UseSymmetry .false. or the",          &
          & " supercell lattice)")')
    endif
    ncom = ncom + 1
    com(ncom) = '#   same k-points in the primitive cell of Rede.Superlattice:'
    do i = 1,3
      ncom = ncom + 1
      write(com(ncom),'("#",7x,3(2x,i6),3x,f12.6)') (mpc1(i,j),j=1,3), mpshc1(i)
    enddo
    if(lbulk1) then
      ncom = ncom + 1
      com(ncom) = '#   (invariant under the parent lattice of the primitive cell)'
    else
      ncom = ncom + 1
      com(ncom) = '#   (NOT invariant under the parent lattice of the primitive cell)'
    endif

    if(.not. lbulk1 .and. lfoundb) then
      call cpwin_kmesh_bulk(mpgridb, mpshiftb, islr, adotpc, nhb, dhb,  &
          mpcb, mpshcb, linvb0)
      write(6,*)
      write(6,*) '  Cheapest mesh with an invariant block in the primitive cell:'
      write(6,*)
      do i = 1,3
        write(6,'(8x,3(2x,i6),3x,f12.6,8x,3(2x,i6),3x,f6.2)')                &
             (mpgridb(i,j),j=1,3), mpshiftb(i), (mpcb(i,j),j=1,3), mpshcb(i)
      enddo
      write(6,'(5x,i6," k-points,",i5," irreducible,  R_cc =",f9.2," bohr")') &
          nkb, nirrb, rccb
      ncom = ncom + 1
      write(com(ncom),'("#   cheapest alternative with an invariant block",  &
          & " in the primitive cell:",i5," irreducible")') nirrb
      do i = 1,3
        ncom = ncom + 1
        write(com(ncom),'("#",7x,3(2x,i6),3x,f12.6,4x,"(",3i5,f6.2," )")')    &
             (mpgridb(i,j),j=1,3), mpshiftb(i), (mpcb(i,j),j=1,3), mpshcb(i)
      enddo
    endif
  endif

  if(lbz) then
    ncom = ncom + 1
    if(lfound0) then
      write(com(ncom),'("# previous block:  R_cc =",f9.2," bohr,",i6,          &
         & " k-points,",i5," irreducible")') rcc0, nk0, nirr0
    else
      write(com(ncom),'("# previous block:  R_cc >",f9.2," bohr,",i6,          &
         & " k-points,",i5," irreducible")') rcc0, nk0, nirr0
    endif
    ncom = ncom + 1
    if(linv0) then
      com(ncom) = '#   invariant under the parent lattice'
    else
      com(ncom) = '#   NOT invariant under the parent lattice'
    endif
  endif

  call cpwin_kmesh_write(trim(fcpwin), 'cpw.out', mpgrid1, mpshift1,     &
      ncom, com, MXCOM)

  write(6,*)
  write(6,*) '  cpw.out written (copy of ',trim(fcpwin),' with the new block)'
  write(6,*)

  deallocate(natom, rat, ratc, atmass, nameat)

  stop

end program cpwin_kmesh




!>  Point group of the parent lattice: the largest group of lattice
!>  operations (finite order, integer in the Buerger reduced basis)
!>  whose strain, the largest |eigenvalue| of (O^T O - I)/2 with
!>  O the cartesian matrix of the operation, is below delta.

subroutine cpwin_kmesh_parent(adot, delta, nh, dh, strain)

! Written 20 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  delta                           !<  largest strain between the lattice and its parent

! output

  integer, intent(out)               ::  nh                              !<  number of operations of the point group of the parent lattice
  integer, intent(out)               ::  dh(3,3,48)                      !<  operations in direct lattice coordinates, n -> dh n for R = sum_j n_j a_j
  real(REAL64), intent(out)          ::  strain                          !<  largest strain of those operations

! local

  integer, parameter         ::  NE = 2                                  !  matrix elements -NE,...,NE in the reduced basis
  integer, parameter         ::  MXCAND = 400
  integer, parameter         ::  MXGRP = 24                              !  largest proper crystallographic point group

  real(REAL64)     ::  adotb(3,3)
  integer          ::  mtb(3,3), mtbinv(3,3), idet
  real(REAL64)     ::  avec(3,3), ainv(3,3), xdet
  real(REAL64)     ::  o(3,3), e(3,3), dev
  integer          ::  g(3,3), gp(3,3), gt(3,3)
  integer          ::  ncand, ntot, ic, idig, iord
  integer          ::  gcand(3,3,MXCAND)
  real(REAL64)     ::  dcand(MXCAND)
  integer, allocatable          ::  mult(:,:)
  integer          ::  ngrp, ngbest, iden
  integer          ::  igrp(MXGRP), igbest(MXGRP)
  real(REAL64)     ::  dgrp, dgbest
  logical          ::  lok, lnew

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer   ::  i, j, k, l, m, n, ip, iq


  call metric_buerger(adot, adotb, mtb)

! exact lattice vectors of the reduced metric (Cholesky, avec^T avec = adotb)

  do j = 1,3
  do i = 1,3
    avec(i,j) = ZERO
  enddo
  enddo
  avec(1,1) = sqrt(adotb(1,1))
  avec(1,2) = adotb(1,2) / avec(1,1)
  avec(1,3) = adotb(1,3) / avec(1,1)
  avec(2,2) = sqrt(adotb(2,2) - avec(1,2)*avec(1,2))
  avec(2,3) = (adotb(2,3) - avec(1,2)*avec(1,3)) / avec(2,2)
  avec(3,3) = sqrt(adotb(3,3) - avec(1,3)*avec(1,3) - avec(2,3)*avec(2,3))

  call cpwin_kmesh_inv3(avec, ainv, xdet)

! proper operations of finite order with small strain

  ncand = 0
  ntot = (2*NE+1)**9
  do ic = 0,ntot-1
    idig = ic
    do j = 1,3
    do i = 1,3
      g(i,j) = mod(idig,2*NE+1) - NE
      idig = idig / (2*NE+1)
    enddo
    enddo

    idet = g(1,1)*(g(2,2)*g(3,3) - g(3,2)*g(2,3))                        &
         + g(2,1)*(g(3,2)*g(1,3) - g(1,2)*g(3,3))                        &
         + g(3,1)*(g(1,2)*g(2,3) - g(2,2)*g(1,3))
    if(idet /= 1) cycle

!   O = avec g avec^-1,  strain = max |O^T O - 1| / 2

    do j = 1,3
    do i = 1,3
      o(i,j) = ZERO
      do k = 1,3
      do l = 1,3
        o(i,j) = o(i,j) + avec(i,k)*g(k,l)*ainv(l,j)
      enddo
      enddo
    enddo
    enddo
    do j = 1,3
    do i = 1,3
      e(i,j) = o(1,i)*o(1,j) + o(2,i)*o(2,j) + o(3,i)*o(3,j)
    enddo
    e(j,j) = e(j,j) - UM
    enddo
    call cpwin_kmesh_symeig_max(e, dev)
    dev = dev / 2
    if(dev > delta) cycle

!   finite order (1, 2, 3, 4 or 6)

    gp = g
    iord = 0
    do n = 1,6
      lok = .TRUE.
      do j = 1,3
      do i = 1,3
        if(i == j .and. gp(i,j) /= 1) lok = .FALSE.
        if(i /= j .and. gp(i,j) /= 0) lok = .FALSE.
      enddo
      enddo
      if(lok) then
        iord = n
        exit
      endif
      gt = matmul(gp,g)
      gp = gt
    enddo
    if(iord == 0 .or. iord == 5) cycle

    if(ncand == MXCAND) then
      write(6,*)
      write(6,*) '  STOPPED in cpwin_kmesh_parent:  too many candidates'
      write(6,*) '  decrease delta'
      write(6,*)

      stop

    endif

    ncand = ncand + 1
    gcand(:,:,ncand) = g
    dcand(ncand) = dev

  enddo

! multiplication table

  allocate(mult(ncand,ncand))

  iden = 0
  do ip = 1,ncand
    lok = .TRUE.
    do j = 1,3
    do i = 1,3
      if(i == j .and. gcand(i,j,ip) /= 1) lok = .FALSE.
      if(i /= j .and. gcand(i,j,ip) /= 0) lok = .FALSE.
    enddo
    enddo
    if(lok) iden = ip
  enddo

  do ip = 1,ncand
  do iq = 1,ncand
    gp = matmul(gcand(:,:,ip),gcand(:,:,iq))
    mult(ip,iq) = 0
    do m = 1,ncand
      lok = .TRUE.
      do j = 1,3
      do i = 1,3
        if(gp(i,j) /= gcand(i,j,m)) lok = .FALSE.
      enddo
      enddo
      if(lok) then
        mult(ip,iq) = m
        exit
      endif
    enddo
  enddo
  enddo

! largest group generated by two candidates (all proper crystallographic
! point groups have at most two generators), ties: smallest strain

  ngbest = 1
  igbest(1) = iden
  dgbest = ZERO

  do ip = 1,ncand
  do iq = ip,ncand

    ngrp = 1
    igrp(1) = iden
    if(ip /= iden) then
      ngrp = ngrp + 1
      igrp(ngrp) = ip
    endif
    if(iq /= iden .and. iq /= ip) then
      ngrp = ngrp + 1
      igrp(ngrp) = iq
    endif

    lok = .TRUE.
    do
      lnew = .FALSE.
      do k = 1,ngrp
      do l = 1,ngrp
        m = mult(igrp(k),igrp(l))
        if(m == 0) then
          lok = .FALSE.
          exit
        endif
        if(.not. any(igrp(1:ngrp) == m)) then
          if(ngrp == MXGRP) then
            lok = .FALSE.
            exit
          endif
          ngrp = ngrp + 1
          igrp(ngrp) = m
          lnew = .TRUE.
        endif
      enddo
      if(.not. lok) exit
      enddo
      if(.not. lok .or. .not. lnew) exit
    enddo

    if(lok) then
      dgrp = ZERO
      do k = 1,ngrp
        dgrp = max(dgrp,dcand(igrp(k)))
      enddo
      if(ngrp > ngbest .or. (ngrp == ngbest .and. dgrp < dgbest)) then
        ngbest = ngrp
        igbest(1:ngrp) = igrp(1:ngrp)
        dgbest = dgrp
      endif
    endif

  enddo
  enddo

  deallocate(mult)

! adds the inversion and goes back to the input basis, n = mtb n'

  call cpwin_kmesh_adj3(mtb, mtbinv, idet)

  nh = 2*ngbest
  do k = 1,ngbest
    g = gcand(:,:,igbest(k))
    gp = matmul(mtb,matmul(g,mtbinv))
    dh(:,:,k) = gp
    dh(:,:,ngbest+k) = -gp
  enddo

  strain = dgbest

  return

end subroutine cpwin_kmesh_parent




!>  Name of the crystal family from the order of the holohedry

subroutine cpwin_kmesh_family(nh, family)

! Written 20 November 2026. claude

  implicit none

  integer, intent(in)                ::  nh                              !<  order of the holohedry
  character(len=40), intent(out)     ::  family                          !<  crystal family

  select case(nh)
    case(48)
      family = 'cubic'
    case(24)
      family = 'hexagonal'
    case(16)
      family = 'tetragonal'
    case(12)
      family = 'rhombohedral'
    case(8)
      family = 'orthorhombic'
    case(4)
      family = 'monoclinic'
    case(2)
      family = 'triclinic'
    case default
      family = 'unknown'
  end select

  return

end subroutine cpwin_kmesh_family




!>  Checks if the supercell lattice of the mesh is invariant under the
!>  operations dh, and returns the operations in supercell coordinates,
!>  xh = M^-T dh M^T  (R = sum_i m_i A_i -> m' = xh m).

subroutine cpwin_kmesh_lsinv(mpgrid, nh, dh, linv, xh)

! Written 20 November 2026. claude

  implicit none

! input

  integer, intent(in)                ::  mpgrid(3,3)                     !<  supercell of the mesh,  A_i = sum_j mpgrid(i,j) a_j
  integer, intent(in)                ::  nh                              !<  number of operations
  integer, intent(in)                ::  dh(3,3,48)                      !<  operations in direct lattice coordinates

! output

  logical, intent(out)               ::  linv                            !<  supercell lattice invariant
  integer, intent(out)               ::  xh(3,3,48)                      !<  operations in supercell coordinates

! local

  integer          ::  mt(3,3), adjt(3,3), idet, y(3,3)

! counters

  integer   ::  i, j, k

  do j = 1,3
  do i = 1,3
    mt(i,j) = mpgrid(j,i)
  enddo
  enddo

  call cpwin_kmesh_adj3(mt, adjt, idet)

  linv = .FALSE.
  if(idet == 0) return

  do k = 1,nh
    y = matmul(adjt,matmul(dh(:,:,k),mt))
    do j = 1,3
    do i = 1,3
      if(mod(y(i,j),idet) /= 0) return
      xh(i,j,k) = y(i,j) / idet
    enddo
    enddo
  enddo

  linv = .TRUE.

  return

end subroutine cpwin_kmesh_lsinv




!>  Checks if the shift (in units of the reciprocal vectors of the
!>  supercell) is invariant,  xh^T s - s  integer for all operations.

subroutine cpwin_kmesh_shiftinv(nh, xh, mpshift, linv)

! Written 20 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nh                              !<  number of operations
  integer, intent(in)                ::  xh(3,3,48)                      !<  operations in supercell coordinates
  real(REAL64), intent(in)           ::  mpshift(3)                      !<  shift of the mesh

! output

  logical, intent(out)               ::  linv                            !<  shift invariant

! local

  real(REAL64)     ::  v

! parameters

  real(REAL64), parameter    ::  EPS = 1.0E-6_REAL64

! counters

  integer   ::  j, k

  linv = .FALSE.
  do k = 1,nh
    do j = 1,3
      v = xh(1,j,k)*mpshift(1) + xh(2,j,k)*mpshift(2) + xh(3,j,k)*mpshift(3) - mpshift(j)
      if(abs(v - nint(v)) > EPS) return
    enddo
  enddo
  linv = .TRUE.

  return

end subroutine cpwin_kmesh_shiftinv




!>  Quality of a mesh:  R_cc, the radius of the shortest star (parent
!>  group) of lattice vectors whose sum of w(R) = <exp(i k.R)> does not
!>  vanish, and the radius of the shortest star where w(R) is not
!>  constant (only for a mesh that is not invariant).

subroutine cpwin_kmesh_rcc(adot, mpgrid, mpshift, nh, dh, rmax,          &
    linv, lfound, rcc, sumw, rni)

! Written 20 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  mpgrid(3,3)                     !<  supercell of the mesh,  A_i = sum_j mpgrid(i,j) a_j
  real(REAL64), intent(in)           ::  mpshift(3)                      !<  shift of the mesh
  integer, intent(in)                ::  nh                              !<  number of operations of the parent lattice
  integer, intent(in)                ::  dh(3,3,48)                      !<  operations in direct lattice coordinates
  real(REAL64), intent(in)           ::  rmax                            !<  largest radius searched

! output

  logical, intent(out)               ::  linv                            !<  mesh invariant under the parent group
  logical, intent(out)               ::  lfound                          !<  a non-vanishing star was found within rmax
  real(REAL64), intent(out)          ::  rcc                             !<  radius of the first non-vanishing star (rmax if not found)
  real(REAL64), intent(out)          ::  sumw                            !<  sum of w(R) in that star
  real(REAL64), intent(out)          ::  rni                             !<  radius of the first star with different w(R) (huge if none)

! local

  logical          ::  llat, lshift
  integer          ::  xh(3,3,48), ops(3,3,48)
  integer          ::  mt(3,3), adjt(3,3), idet
  real(REAL64)     ::  gb(3,3), gbinv(3,3), xdet
  integer          ::  mm(3)
  logical, allocatable  ::  visited(:,:,:)
  integer          ::  p(3), q(3), y(3), mem(3,48), nm
  real(REAL64)     ::  r2, rq, rstar, sw, w, wmin, wmax
  logical          ::  lnew, lin

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-6_REAL64
  real(REAL64), parameter    ::  PI = 3.14159265358979323846_REAL64

! counters

  integer   ::  i, j, k, l, h, i1, i2, i3


  call cpwin_kmesh_lsinv(mpgrid, nh, dh, llat, xh)
  lshift = .FALSE.
  if(llat) call cpwin_kmesh_shiftinv(nh, xh, mpshift, lshift)
  linv = llat .and. lshift

  do j = 1,3
  do i = 1,3
    mt(i,j) = mpgrid(j,i)
  enddo
  enddo
  call cpwin_kmesh_adj3(mt, adjt, idet)

! enumerates the supercell lattice if invariant, otherwise the lattice

  if(llat) then
    do j = 1,3
    do i = 1,3
      gb(i,j) = ZERO
      do k = 1,3
      do l = 1,3
        gb(i,j) = gb(i,j) + mpgrid(i,k)*adot(k,l)*mpgrid(j,l)
      enddo
      enddo
    enddo
    enddo
    ops(:,:,1:nh) = xh(:,:,1:nh)
  else
    gb = adot
    ops(:,:,1:nh) = dh(:,:,1:nh)
  endif

  call cpwin_kmesh_inv3(gb, gbinv, xdet)
  do i = 1,3
    mm(i) = int(1.25_REAL64*rmax*sqrt(abs(gbinv(i,i)))) + 1
  enddo

  allocate(visited(-mm(1):mm(1),-mm(2):mm(2),-mm(3):mm(3)))
  visited = .FALSE.

  lfound = .FALSE.
  rcc = rmax
  sumw = ZERO
  rni = 1.0E30_REAL64

  do i3 = -mm(3),mm(3)
  do i2 = -mm(2),mm(2)
  do i1 = -mm(1),mm(1)
    if(i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
    if(visited(i1,i2,i3)) cycle
    p(1) = i1
    p(2) = i2
    p(3) = i3
    r2 = ZERO
    do j = 1,3
    do i = 1,3
      r2 = r2 + p(i)*gb(i,j)*p(j)
    enddo
    enddo
    if(r2 > rmax*rmax) cycle

!   the star of p

    nm = 0
    sw = ZERO
    wmin = 1.0E30_REAL64
    wmax = -1.0E30_REAL64
    rstar = 1.0E30_REAL64
    do h = 1,nh
      do i = 1,3
        q(i) = ops(i,1,h)*p(1) + ops(i,2,h)*p(2) + ops(i,3,h)*p(3)
      enddo
      lnew = .TRUE.
      do k = 1,nm
        if(q(1) == mem(1,k) .and. q(2) == mem(2,k) .and. q(3) == mem(3,k)) then
          lnew = .FALSE.
          exit
        endif
      enddo
      if(.not. lnew) cycle
      nm = nm + 1
      mem(:,nm) = q
      lin = .TRUE.
      do i = 1,3
        if(abs(q(i)) > mm(i)) lin = .FALSE.
      enddo
      if(lin) visited(q(1),q(2),q(3)) = .TRUE.

      if(llat) then
        w = cos(2*PI*(mpshift(1)*q(1) + mpshift(2)*q(2) + mpshift(3)*q(3)))
      else
        do i = 1,3
          y(i) = adjt(i,1)*q(1) + adjt(i,2)*q(2) + adjt(i,3)*q(3)
        enddo
        if(mod(y(1),idet) == 0 .and. mod(y(2),idet) == 0 .and. mod(y(3),idet) == 0) then
          w = cos(2*PI*(mpshift(1)*(y(1)/idet) + mpshift(2)*(y(2)/idet) +           &
                        mpshift(3)*(y(3)/idet)))
        else
          w = ZERO
        endif
      endif

      sw = sw + w
      wmin = min(wmin,w)
      wmax = max(wmax,w)
      rq = ZERO
      do j = 1,3
      do i = 1,3
        rq = rq + q(i)*gb(i,j)*q(j)
      enddo
      enddo
      rstar = min(rstar,sqrt(rq))
    enddo

    if(abs(sw) > EPS .and. rstar < rcc) then
      rcc = rstar
      sumw = sw
      lfound = .TRUE.
    endif
    if(wmax - wmin > EPS .and. rstar < rni) rni = rstar

  enddo
  enddo
  enddo

  deallocate(visited)

  return

end subroutine cpwin_kmesh_rcc




!>  Buerger reduced basis of the supercell of the mesh (same mesh),
!>  A'_j = sum_i A_i mtb(i,j),  s'_j = sum_i mtb(i,j) s_i  (modulo 1),
!>  then rows ordered for a positive determinant and a large diagonal.

subroutine cpwin_kmesh_reduce(adot, mpgrid, mpshift, mpgridr, mpshiftr)

! Written 20 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  mpgrid(3,3)                     !<  supercell of the mesh,  A_i = sum_j mpgrid(i,j) a_j
  real(REAL64), intent(in)           ::  mpshift(3)                      !<  shift of the mesh

! output

  integer, intent(out)               ::  mpgridr(3,3)                    !<  reduced supercell of the mesh
  real(REAL64), intent(out)          ::  mpshiftr(3)                     !<  shift in the reduced basis

! local

  real(REAL64)     ::  gs(3,3), gsb(3,3)
  integer          ::  mtb(3,3)

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-8_REAL64

! counters

  integer   ::  i, j, k, l

  do j = 1,3
  do i = 1,3
    gs(i,j) = ZERO
    do k = 1,3
    do l = 1,3
      gs(i,j) = gs(i,j) + mpgrid(i,k)*adot(k,l)*mpgrid(j,l)
    enddo
    enddo
  enddo
  enddo

  call metric_buerger(gs, gsb, mtb)

  do j = 1,3
    do l = 1,3
      mpgridr(j,l) = 0
      do i = 1,3
        mpgridr(j,l) = mpgridr(j,l) + mtb(i,j)*mpgrid(i,l)
      enddo
    enddo
    mpshiftr(j) = ZERO
    do i = 1,3
      mpshiftr(j) = mpshiftr(j) + mtb(i,j)*mpshift(i)
    enddo
    mpshiftr(j) = mpshiftr(j) - floor(mpshiftr(j) + EPS)
    if(abs(mpshiftr(j)) < EPS) mpshiftr(j) = ZERO
  enddo

  call cpwin_kmesh_pretty(mpgridr, mpshiftr)

  return

end subroutine cpwin_kmesh_reduce




!>  Number of irreducible k-points (as in cpw, with time reversal)

subroutine cpwin_kmesh_nirr(adot, mpgrid, mpshift, ntrans, mtrx, nirr)

! Written 20 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  mpgrid(3,3)                     !<  supercell of the mesh,  A_i = sum_j mpgrid(i,j) a_j
  real(REAL64), intent(in)           ::  mpshift(3)                      !<  shift of the mesh
  integer, intent(in)                ::  ntrans                          !<  number of operations of the crystal point group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  operations in reciprocal lattice coordinates

! output

  integer, intent(out)               ::  nirr                            !<  number of irreducible k-points

! local

  integer          ::  nx, ny, nz, kmat(3,3)
  real(REAL64)     ::  sx, sy, sz

  call kgrid_box(mpgrid, mpshift, nx, ny, nz, sx, sy, sz, kmat)

  call size_mxdnrk(nx,ny,nz, sx,sy,sz, kmat, adot, ntrans, mtrx, nirr)

  return

end subroutine cpwin_kmesh_nirr




!>  Same k-points in the primitive cell of the superlattice (reduced
!>  basis), and invariance under the parent lattice of the primitive cell.

subroutine cpwin_kmesh_bulk(mpgrid, mpshift, islr, adotpc, nhb, dhb,     &
    mpc, mpshiftc, linvb)

! Written 20 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mpgrid(3,3)                     !<  supercell of the mesh in superlattice coordinates
  real(REAL64), intent(in)           ::  mpshift(3)                      !<  shift of the mesh
  integer, intent(in)                ::  islr(3,3)                       !<  line i: superlattice vector i in primitive cell coordinates
  real(REAL64), intent(in)           ::  adotpc(3,3)                     !<  metric of the primitive cell
  integer, intent(in)                ::  nhb                             !<  number of operations of the parent of the primitive cell
  integer, intent(in)                ::  dhb(3,3,48)                     !<  operations in primitive cell coordinates

! output

  integer, intent(out)               ::  mpc(3,3)                        !<  supercell of the mesh in primitive cell coordinates
  real(REAL64), intent(out)          ::  mpshiftc(3)                     !<  shift in that basis
  logical, intent(out)               ::  linvb                           !<  invariant under the parent of the primitive cell

! local

  integer          ::  xh(3,3,48), mpc0(3,3)
  logical          ::  llat, lshift

  mpc0 = matmul(mpgrid,islr)

  call cpwin_kmesh_reduce(adotpc, mpc0, mpshift, mpc, mpshiftc)

  call cpwin_kmesh_lsinv(mpc, nhb, dhb, llat, xh)
  lshift = .FALSE.
  if(llat) call cpwin_kmesh_shiftinv(nhb, xh, mpshiftc, lshift)
  linvb = llat .and. lshift

  return

end subroutine cpwin_kmesh_bulk




!>  Equivalent block with the rows permuted and with signs changed so that
!>  the determinant is positive and the diagonal as large as possible
!>  (the same mesh, a row -A_i has the shift -s_i).

subroutine cpwin_kmesh_pretty(mpgrid, mpshift)

! Written 20 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input and output

  integer, intent(inout)             ::  mpgrid(3,3)                     !<  supercell of the mesh
  real(REAL64), intent(inout)        ::  mpshift(3)                      !<  shift of the mesh

! local

  integer          ::  iperm(3,6), m(3,3), mbest(3,3), adj(3,3), idet
  real(REAL64)     ::  s(3), sbest(3)
  integer          ::  npos, nsum, nneg, nposb, nsumb, nnegb, isg
  logical          ::  lbetter

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-8_REAL64

! counters

  integer   ::  i, j, ip, is

  data iperm / 1,2,3,  1,3,2,  2,1,3,  2,3,1,  3,1,2,  3,2,1 /

  nposb = -1
  nsumb = 0
  nnegb = 0
  mbest = mpgrid
  sbest = mpshift

  do ip = 1,6
  do is = 0,7
    do i = 1,3
      isg = 1
      if(mod(is/(2**(i-1)),2) == 1) isg = -1
      do j = 1,3
        m(i,j) = isg*mpgrid(iperm(i,ip),j)
      enddo
      s(i) = isg*mpshift(iperm(i,ip))
      s(i) = s(i) - floor(s(i) + EPS)
      if(abs(s(i)) < EPS) s(i) = ZERO
    enddo

    call cpwin_kmesh_adj3(m, adj, idet)
    if(idet <= 0) cycle

    npos = 0
    nsum = 0
    nneg = 0
    do i = 1,3
      if(m(i,i) > 0) npos = npos + 1
      nsum = nsum + m(i,i)
      do j = 1,3
        if(m(i,j) < 0) nneg = nneg + 1
      enddo
    enddo

    lbetter = .FALSE.
    if(npos > nposb) then
      lbetter = .TRUE.
    elseif(npos == nposb) then
      if(nsum > nsumb) then
        lbetter = .TRUE.
      elseif(nsum == nsumb .and. nneg < nnegb) then
        lbetter = .TRUE.
      endif
    endif

    if(lbetter) then
      nposb = npos
      nsumb = nsum
      nnegb = nneg
      mbest = m
      sbest = s
    endif
  enddo
  enddo

  mpgrid = mbest
  mpshift = sbest

  return

end subroutine cpwin_kmesh_pretty




!>  Search of invariant meshes (Hermite normal forms of the supercell,
!>  shifts 0 or 1/2) with R_cc >= target and the smallest number of
!>  irreducible k-points.  Prints the Pareto front (best R_cc for each
!>  number of irreducible k-points).

subroutine cpwin_kmesh_search(adot, nh, dh, ntrans, mtrx, geff,          &
    lrede, islr, adotpc, nhb, dhb, rtarget,                              &
    lfound1, mpgrid1, mpshift1, nk1, nirr1, rcc1, lbulk1,                &
    lfoundb, mpgridb, mpshiftb, nkb, nirrb, rccb)

! Written 20 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  nh                              !<  number of operations of the parent lattice
  integer, intent(in)                ::  dh(3,3,48)                      !<  operations of the parent lattice (direct lattice coordinates)
  integer, intent(in)                ::  ntrans                          !<  number of operations of the crystal point group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  operations of the crystal (reciprocal lattice coordinates)
  integer, intent(in)                ::  geff                            !<  largest orbit of a k-point
  logical, intent(in)                ::  lrede                           !<  superlattice information present
  integer, intent(in)                ::  islr(3,3)                       !<  line i: superlattice vector i in primitive cell coordinates
  real(REAL64), intent(in)           ::  adotpc(3,3)                     !<  metric of the primitive cell of the superlattice
  integer, intent(in)                ::  nhb                             !<  number of operations of the parent of the primitive cell
  integer, intent(in)                ::  dhb(3,3,48)                     !<  operations in primitive cell coordinates
  real(REAL64), intent(in)           ::  rtarget                         !<  target R_cc

! output

  logical, intent(out)               ::  lfound1                         !<  a mesh was found
  integer, intent(out)               ::  mpgrid1(3,3)                    !<  chosen mesh
  real(REAL64), intent(out)          ::  mpshift1(3)                     !<  chosen shift
  integer, intent(out)               ::  nk1                             !<  number of k-points
  integer, intent(out)               ::  nirr1                           !<  number of irreducible k-points
  real(REAL64), intent(out)          ::  rcc1                            !<  R_cc of the chosen mesh
  logical, intent(out)               ::  lbulk1                          !<  invariant block in the primitive cell
  logical, intent(out)               ::  lfoundb                         !<  a mesh with invariant block in the primitive cell was found
  integer, intent(out)               ::  mpgridb(3,3)                    !<  cheapest such mesh
  real(REAL64), intent(out)          ::  mpshiftb(3)                     !<  its shift
  integer, intent(out)               ::  nkb                             !<  number of k-points
  integer, intent(out)               ::  nirrb                           !<  number of irreducible k-points
  real(REAL64), intent(out)          ::  rccb                            !<  R_cc

! local

  integer, parameter         ::  NCAP = 4096                             !  largest number of k-points
  integer, parameter         ::  MXDSHAPE = 32                           !  largest determinant of a stored shape
  integer, parameter         ::  MXSHAPE = 5000                          !  largest number of stored shapes

  integer          ::  nfull                                             !  complete search up to nfull k-points
  integer          ::  nmaxs                                             !  search stops after nmaxs k-points
  integer          ::  nshape
  integer, allocatable          ::  ishape(:,:,:)
  integer, allocatable          ::  dshape(:)
  real(REAL64)     ::  rmax, rrun
  integer          ::  m(3,3), xh(3,3,48)
  integer          ::  ia, ic, ifz, ib, id, ie, iq, nsc, nr, nk, nlast
  logical          ::  llat

! Pareto front, best mesh for each number of irreducible points

  real(REAL64), allocatable     ::  bestr(:)
  integer, allocatable          ::  bestm(:,:,:)
  real(REAL64), allocatable     ::  bests(:,:)
  integer, allocatable          ::  bestnk(:)
  logical, allocatable          ::  bestlb(:)
  character(len=1) ::  cgt, cb

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-6_REAL64

! counters

  integer   ::  i, j, n

  nfull = 128
  if(nh < 8) nfull = 48

  rmax = 1.5_REAL64*rtarget
  nshape = 0
  nmaxs = NCAP

  allocate(ishape(3,3,MXSHAPE), dshape(MXSHAPE))
  allocate(bestr(NCAP), bestm(3,3,NCAP), bests(3,NCAP), bestnk(NCAP), bestlb(NCAP))
  bestr = -UM
  bestlb = .FALSE.

  lfound1 = .FALSE.
  lfoundb = .FALSE.
  lbulk1 = .FALSE.
  nirr1 = 0
  nirrb = 0
  nk1 = 0
  nkb = 0
  rcc1 = ZERO
  rccb = ZERO
  mpgrid1 = 0
  mpgridb = 0
  mpshift1 = ZERO
  mpshiftb = ZERO

  do nk = 1,NCAP

    if(nk > nmaxs) exit

    if(nk <= nfull) then

!     all Hermite normal forms (a 0 0 / b c 0 / d e f), 0 <= b,d < a, 0 <= e < c

      do ia = 1,nk
        if(mod(nk,ia) /= 0) cycle
        do ic = 1,nk/ia
          if(mod(nk/ia,ic) /= 0) cycle
          ifz = nk/(ia*ic)
          do ib = 0,ia-1
          do id = 0,ia-1
          do ie = 0,ic-1
            m(1,1) = ia
            m(1,2) = 0
            m(1,3) = 0
            m(2,1) = ib
            m(2,2) = ic
            m(2,3) = 0
            m(3,1) = id
            m(3,2) = ie
            m(3,3) = ifz
            call cpwin_kmesh_search_try()
          enddo
          enddo
          enddo
        enddo
      enddo

    else

!     scaled shapes  n Q

      do iq = 1,nshape
        if(mod(nk,dshape(iq)) /= 0) cycle
        nr = nk / dshape(iq)
        nsc = nint(real(nr,REAL64)**(UM/3))
        if(nsc < 2 .or. nsc*nsc*nsc /= nr) cycle
        m = nsc*ishape(:,:,iq)
        call cpwin_kmesh_search_try()
      enddo

    endif

  enddo

! Pareto front up to the chosen mesh (complete up to there)

  write(6,'("   irred      N    R_cc(bohr)      kgrid_Monkhorst_Pack (lines)",  &
       &   "               shift")')

  nlast = NCAP
  if(lfound1) nlast = nirr1
  rrun = ZERO
  do n = 1,nlast
    if(bestr(n) > rrun + EPS) then
      rrun = bestr(n)
      cgt = ' '
      if(bestr(n) >= rmax) cgt = '>'
      cb = ' '
      if(bestlb(n)) cb = 'B'
      write(6,'(2x,i6,i7,3x,a1,f9.2,4x,3(3i4,2x),2x,3f5.2,2x,a1)')        &
          n, bestnk(n), cgt, bestr(n), ((bestm(i,j,n),j=1,3),i=1,3),       &
          (bests(i,n),i=1,3), cb
    endif
  enddo
  if(lrede) then
    write(6,*)
    write(6,*) '  B: the same k-points in the primitive cell form an invariant mesh'
  endif

  deallocate(ishape, dshape)
  deallocate(bestr, bestm, bests, bestnk, bestlb)

  return

  contains


!   tries all invariant shifts for the supercell m

    subroutine cpwin_kmesh_search_try()

      implicit none

      integer          ::  mr(3,3), mpc(3,3)
      real(REAL64)     ::  s(3), sr(3), sc(3)
      logical          ::  linv, lfound, lbulk, lbetter, lgood
      real(REAL64)     ::  rcc, sumw, rni, rdom, r2
      integer          ::  nirr, nlow, nz, nz1, k, ix, it
      integer          ::  t(3)

      call cpwin_kmesh_lsinv(m, nh, dh, llat, xh)
      if(.not. llat) return

      if(nk <= MXDSHAPE .and. nshape < MXSHAPE) then
        nshape = nshape + 1
        ishape(:,:,nshape) = m
        dshape(nshape) = nk
      endif

      do it = 0,7
        t(1) = mod(it,2)
        t(2) = mod(it/2,2)
        t(3) = mod(it/4,2)

!       shift invariance, xh^T t - t even

        lgood = .TRUE.
        do k = 1,nh
          do j = 1,3
            ix = xh(1,j,k)*t(1) + xh(2,j,k)*t(2) + xh(3,j,k)*t(3) - t(j)
            if(mod(ix,2) /= 0) lgood = .FALSE.
          enddo
          if(.not. lgood) exit
        enddo
        if(.not. lgood) cycle

        do i = 1,3
          s(i) = t(i) / 2.0_REAL64
        enddo

        call cpwin_kmesh_reduce(adot, m, s, mr, sr)

        if(nh == 2) then

!         only +-1:  stars {R,-R} never cancel, R_cc is the shortest
!         supercell vector (the shortest vector of the reduced basis)

          rcc = rmax
          do i = 1,3
            r2 = ZERO
            do k = 1,3
            do ix = 1,3
              r2 = r2 + mr(i,k)*adot(k,ix)*mr(i,ix)
            enddo
            enddo
            rcc = min(rcc,sqrt(r2))
          enddo

        else

          call cpwin_kmesh_rcc(adot, mr, sr, nh, dh, rmax,               &
              linv, lfound, rcc, sumw, rni)

        endif

!       lower bound of the number of irreducible points

        nlow = (nk + geff - 1) / geff

!       dominated (not on the Pareto front and below the target)

        rdom = -UM
        do k = 1,nlow
          rdom = max(rdom,bestr(k))
        enddo
        if(rcc < rtarget .and. rcc <= rdom + EPS) cycle

        if(rcc >= rtarget .and. lfound1 .and. nlow > nirr1) then
          if(.not. lrede) cycle
          if(lfoundb .and. nlow > nirrb) cycle
        endif

        call cpwin_kmesh_nirr(adot, mr, sr, ntrans, mtrx, nirr)

        lbulk = .FALSE.
        if(lrede) call cpwin_kmesh_bulk(mr, sr, islr, adotpc, nhb, dhb, mpc, sc, lbulk)

!       Pareto front

        if(rcc > bestr(nirr) + EPS .or.                                  &
           (abs(rcc - bestr(nirr)) <= EPS .and. lbulk .and. .not. bestlb(nirr))) then
          bestr(nirr) = rcc
          bestm(:,:,nirr) = mr
          bests(:,nirr) = sr
          bestnk(nirr) = nk
          bestlb(nirr) = lbulk
        endif

        if(rcc < rtarget) cycle

        nz = 0
        nz1 = 0
        do i = 1,3
          if(abs(sr(i)) > EPS) nz = nz + 1
          if(abs(mpshift1(i)) > EPS) nz1 = nz1 + 1
        enddo

!       smallest nirr, then invariant in the primitive cell, larger R_cc,
!       fewer points, fewer shifts

        lbetter = .FALSE.
        if(.not. lfound1) then
          lbetter = .TRUE.
        elseif(nirr < nirr1) then
          lbetter = .TRUE.
        elseif(nirr == nirr1) then
          if(lbulk .and. .not. lbulk1) then
            lbetter = .TRUE.
          elseif(lbulk .eqv. lbulk1) then
            if(rcc > rcc1 + EPS) then
              lbetter = .TRUE.
            elseif(abs(rcc - rcc1) <= EPS) then
              if(nk < nk1) then
                lbetter = .TRUE.
              elseif(nk == nk1 .and. nz < nz1) then
                lbetter = .TRUE.
              endif
            endif
          endif
        endif

        if(lbetter) then
          lfound1 = .TRUE.
          mpgrid1 = mr
          mpshift1 = sr
          nk1 = nk
          nirr1 = nirr
          rcc1 = rcc
          lbulk1 = lbulk
        endif

        if(lbulk) then
          lbetter = .FALSE.
          if(.not. lfoundb) then
            lbetter = .TRUE.
          elseif(nirr < nirrb) then
            lbetter = .TRUE.
          elseif(nirr == nirrb .and. rcc > rccb + EPS) then
            lbetter = .TRUE.
          endif
          if(lbetter) then
            lfoundb = .TRUE.
            mpgridb = mr
            mpshiftb = sr
            nkb = nk
            nirrb = nirr
            rccb = rcc
          endif
        endif

!       new end of the search (no mesh with more points can have fewer irreducible points)

        nmaxs = min(NCAP, nirr1*geff)
        if(lrede) then
          if(lfoundb) then
            nmaxs = min(NCAP, max(nirr1,nirrb)*geff)
          else
            nmaxs = min(NCAP, 4*nirr1*geff)
          endif
        endif

      enddo

      return

    end subroutine cpwin_kmesh_search_try

end subroutine cpwin_kmesh_search




!>  Writes a verbatim copy of the input file where only the block
!>  kgrid_Monkhorst_Pack is replaced (the old block becomes comments).
!>  If there is no block the new one is added at the end.

subroutine cpwin_kmesh_write(fin, fout, mpgrid, mpshift, ncom, com, mxcom)

! Written 20 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len=*), intent(in)       ::  fin                             !<  input file
  character(len=*), intent(in)       ::  fout                            !<  output file
  integer, intent(in)                ::  mpgrid(3,3)                     !<  new block
  real(REAL64), intent(in)           ::  mpshift(3)                      !<  new shift
  integer, intent(in)                ::  mxcom                           !<  dimension of com
  integer, intent(in)                ::  ncom                            !<  number of comment lines
  character(len=*), intent(in)       ::  com(mxcom)                      !<  comment lines

! local

  integer                            ::  io
  character(len=:), allocatable      ::  line
  character(len=:), allocatable      ::  lines(:)
  integer, allocatable               ::  lens(:)
  integer                            ::  nl, lmax, ios
  integer                            ::  istart, iend
  character(len=80)                  ::  t1, t2

! counters

  integer   ::  i, j

  io = 21

! reads all the lines (the input may be the output file)

  open(unit = io, file = fin, status = 'OLD', form = 'FORMATTED', action = 'READ')

  nl = 0
  lmax = 1
  do
    call cpwin_kmesh_getline(io, line, ios)
    if(ios /= 0) exit
    nl = nl + 1
    lmax = max(lmax,len(line))
  enddo

  rewind(io)

  allocate(character(len=lmax) :: lines(nl))
  allocate(lens(nl))

  do i = 1,nl
    call cpwin_kmesh_getline(io, line, ios)
    lens(i) = len(line)
    lines(i) = line
  enddo

  close(unit = io)

! finds the block

  istart = 0
  iend = 0
  do i = 1,nl
    call cpwin_kmesh_tokens(lines(i)(1:lens(i)), t1, t2)
    if(istart == 0) then
      if(t1 == '%block' .and. t2 == 'kgridmonkhorstpack') istart = i
    elseif(iend == 0) then
      if(t1 == '%endblock') iend = i
    endif
  enddo
  if(iend == 0) istart = 0

! writes

  open(unit = io, file = fout, status = 'REPLACE', form = 'FORMATTED')

  if(istart > 0) then
    do i = 1,istart-1
      write(io,'(a)') lines(i)(1:lens(i))
    enddo
  else
    do i = 1,nl
      write(io,'(a)') lines(i)(1:lens(i))
    enddo
    write(io,*)
  endif

  do i = 1,ncom
    write(io,'(a)') trim(com(i))
  enddo

  if(istart > 0) then
    do i = istart,iend
      write(io,'("# ",a)') lines(i)(1:lens(i))
    enddo
  endif

  write(io,'("%block kgrid_Monkhorst_Pack")')
  do i = 1,3
    write(io,'(8x,3(2x,i6),3x,f12.6)')  (mpgrid(i,j),j=1,3), mpshift(i)
  enddo
  write(io,'("%endblock kgrid_Monkhorst_Pack")')

  if(istart > 0) then
    do i = iend+1,nl
      write(io,'(a)') lines(i)(1:lens(i))
    enddo
  endif

  close(unit = io)

  deallocate(lines, lens)

  return

  contains


!   reads a line of any length (non-advancing input)

    subroutine cpwin_kmesh_getline(io, line, ios)

      implicit none

      integer, intent(in)                            ::  io                  !<  tape number
      character(len=:), allocatable, intent(out)     ::  line                !<  the line, without the end of record
      integer, intent(out)                           ::  ios                 !<  0 if a line was read

      character(len=256)    ::  buf
      integer               ::  nsz

      line = ''
      do
        read(io,'(a)',advance='NO',size=nsz,iostat=ios) buf
        if(nsz > 0) line = line // buf(1:nsz)
        if(is_iostat_eor(ios)) then
          ios = 0
          return
        endif
        if(is_iostat_end(ios)) then
          if(len(line) > 0) ios = 0
          return
        endif
        if(ios /= 0) return
      enddo

    end subroutine cpwin_kmesh_getline

end subroutine cpwin_kmesh_write




!>  First two tokens of a line, reduced as the esdf labels
!>  (comments removed, lower case, without '.', '_' and '-')

subroutine cpwin_kmesh_tokens(line, t1, t2)

! Written 20 November 2026. claude

  implicit none

  character(len=*), intent(in)       ::  line                            !<  line of the file
  character(len=*), intent(out)      ::  t1                              !<  first token
  character(len=*), intent(out)      ::  t2                              !<  second token

  integer          ::  n, i, k, ic, ntok
  character(len=1) ::  c
  logical          ::  lsep, lin

  t1 = ' '
  t2 = ' '

  n = len(line)
  do i = 1,len(line)
    if(line(i:i) == '#' .or. line(i:i) == ';' .or. line(i:i) == '!') then
      n = i - 1
      exit
    endif
  enddo

  ntok = 0
  lin = .FALSE.
  k = 0
  do i = 1,n
    c = line(i:i)
    lsep = (c == ' ' .or. c == char(9) .or. c == '=' .or. c == ':')
    if(lsep) then
      if(lin) then
        lin = .FALSE.
        if(ntok == 2) exit
      endif
      cycle
    endif
    if(.not. lin) then
      lin = .TRUE.
      ntok = ntok + 1
      k = 0
    endif
    if(c == '.' .or. c == '_' .or. c == '-') cycle
    ic = ichar(c)
    if(ic >= ichar('A') .and. ic <= ichar('Z')) c = char(ic - ichar('A') + ichar('a'))
    k = k + 1
    if(ntok == 1 .and. k <= len(t1)) t1(k:k) = c
    if(ntok == 2 .and. k <= len(t2)) t2(k:k) = c
  enddo

  return

end subroutine cpwin_kmesh_tokens




!>  Largest absolute value of the eigenvalues of a real symmetric 3x3
!>  matrix (closed form, trigonometric solution of the cubic)

subroutine cpwin_kmesh_symeig_max(e, emax)

! Written 20 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  real(REAL64), intent(in)           ::  e(3,3)                          !<  symmetric matrix
  real(REAL64), intent(out)          ::  emax                            !<  largest |eigenvalue|

  real(REAL64)     ::  p1, p2, p, q, r, phi, b(3,3), detb
  real(REAL64)     ::  ev1, ev2, ev3

  real(REAL64), parameter    ::  UM = 1.0_REAL64
  real(REAL64), parameter    ::  PI = 3.14159265358979323846_REAL64

  integer   ::  i, j

  p1 = e(1,2)*e(1,2) + e(1,3)*e(1,3) + e(2,3)*e(2,3)
  q = (e(1,1) + e(2,2) + e(3,3)) / 3
  p2 = (e(1,1)-q)*(e(1,1)-q) + (e(2,2)-q)*(e(2,2)-q) + (e(3,3)-q)*(e(3,3)-q) + 2*p1
  p = sqrt(p2/6)

  if(p < 1.0E-14_REAL64) then
    emax = abs(q)
    return
  endif

  do j = 1,3
  do i = 1,3
    b(i,j) = e(i,j) / p
  enddo
  b(j,j) = b(j,j) - q/p
  enddo
  detb = b(1,1)*(b(2,2)*b(3,3) - b(2,3)*b(3,2))                          &
       - b(1,2)*(b(2,1)*b(3,3) - b(2,3)*b(3,1))                          &
       + b(1,3)*(b(2,1)*b(3,2) - b(2,2)*b(3,1))
  r = detb / 2
  r = max(-UM,min(UM,r))
  phi = acos(r) / 3

  ev1 = q + 2*p*cos(phi)
  ev3 = q + 2*p*cos(phi + 2*PI/3)
  ev2 = 3*q - ev1 - ev3

  emax = max(abs(ev1),abs(ev2),abs(ev3))

  return

end subroutine cpwin_kmesh_symeig_max




!>  Adjugate and determinant of an integer 3x3 matrix,  m adj = det I

subroutine cpwin_kmesh_adj3(m, adj, idet)

! Written 20 November 2026. claude

  implicit none

  integer, intent(in)                ::  m(3,3)                          !<  matrix
  integer, intent(out)               ::  adj(3,3)                        !<  adjugate
  integer, intent(out)               ::  idet                            !<  determinant

  adj(1,1) = m(2,2)*m(3,3) - m(2,3)*m(3,2)
  adj(1,2) = m(1,3)*m(3,2) - m(1,2)*m(3,3)
  adj(1,3) = m(1,2)*m(2,3) - m(1,3)*m(2,2)
  adj(2,1) = m(2,3)*m(3,1) - m(2,1)*m(3,3)
  adj(2,2) = m(1,1)*m(3,3) - m(1,3)*m(3,1)
  adj(2,3) = m(1,3)*m(2,1) - m(1,1)*m(2,3)
  adj(3,1) = m(2,1)*m(3,2) - m(2,2)*m(3,1)
  adj(3,2) = m(1,2)*m(3,1) - m(1,1)*m(3,2)
  adj(3,3) = m(1,1)*m(2,2) - m(1,2)*m(2,1)

  idet = m(1,1)*adj(1,1) + m(1,2)*adj(2,1) + m(1,3)*adj(3,1)

  return

end subroutine cpwin_kmesh_adj3




!>  Inverse and determinant of a real 3x3 matrix

subroutine cpwin_kmesh_inv3(a, ainv, det)

! Written 20 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  real(REAL64), intent(in)           ::  a(3,3)                          !<  matrix
  real(REAL64), intent(out)          ::  ainv(3,3)                       !<  inverse
  real(REAL64), intent(out)          ::  det                             !<  determinant

  integer   ::  i, j

  ainv(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
  ainv(1,2) = a(1,3)*a(3,2) - a(1,2)*a(3,3)
  ainv(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
  ainv(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
  ainv(2,2) = a(1,1)*a(3,3) - a(1,3)*a(3,1)
  ainv(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
  ainv(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
  ainv(3,2) = a(1,2)*a(3,1) - a(1,1)*a(3,2)
  ainv(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

  det = a(1,1)*ainv(1,1) + a(1,2)*ainv(2,1) + a(1,3)*ainv(3,1)

  do j = 1,3
  do i = 1,3
    ainv(i,j) = ainv(i,j) / det
  enddo
  enddo

  return

end subroutine cpwin_kmesh_inv3
