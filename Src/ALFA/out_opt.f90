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

!>  This subroutine calculates the bands on an uniform grid
!>  and the oscillator strengths for later processing.

  subroutine out_opt(diag_type, lworkers,                                &
    title, subtitle,                                                     &
    emax, flgdal, flgpsd, iguess, epspsi, icmax, ztot,                   &
    adot, ntype, natom, rat, ntrans, mtrx,                               &
    ng, kgv, phase, conj,                                                &
    ns, inds, kmax, indv, ek,                                            &
    sfact, icmplx,                                                       &
    veff,                                                                &
    nqnl, delqnl, vkb, nkb,                                              &
    latorb, norbat, nqwf, delqwf, wvfao, lorb,                           &     
    mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)          

! Adapted from out_dos, CLR
! Modified, latorb, 7 June 2020. JLM
! Modified, vmax, vmin, 27 November 2020. JLM
! Modified allk, details of workers, 7 December 2020. JLM

! copyright  Carlos Loia Reis/Jose Luis Martins/INESC-MN

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: REAL32 = selected_real_kind(6)


! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxdlao                          !< array dimension of orbital per atom type

  character(len=50), intent(in)      ::  title                           !<  title for plots
  character(len=140), intent(in)     ::  subtitle                        !<  subtitle for plots

  character(len=4), intent(in)       ::  diag_type                       !<  selects diagonalization, 'pw  ','ao  ','aojc'
  logical, intent(in)                ::  lworkers                        !<  use lworkers

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4)                   ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  integer, intent(in)                ::  iguess                          !<  if guess eigenvectors are available, iguess = 1, otherwise iguess = 0
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase

  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|
  integer, intent(in)                ::  indv(mxdcub)                    !<  kgv(i,indv(jadd)) is the g-vector associated with jadd. jadd is defined by the g-vector components and kmax
  real(REAL64), intent(in)           ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j

  complex(REAL64), intent(in)        ::  sfact(mxdtyp,mxdnst)            !<  structure factor
  integer, intent(in)                ::  icmplx                          !<  indicates if the structure factor is complex

  complex(REAL64), intent(in)        ::  veff(mxdnst)                    !<  ionic potential (local+Hartree+XC) for the prototype g-vector in star j

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)        ::   vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)    !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<   KB pseudo.  normalization for atom k, ang. mom. l

  logical, intent(in)                ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  integer, intent(in)                ::  norbat(mxdtyp)                  !< number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !< number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !< step used in the wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !< wavefunction for atom k, ang. mom. l (normalized to vcell)
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !< angular momentum of orbital n of atom k


! local allocatable arrays

  integer, allocatable               ::  kmap(:,:,:,:)                   !  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
  integer, allocatable               ::  nband(:)                        !  number of bands for each k-points
  integer, allocatable               ::  indk(:,:)                       !  index of the six k-points neighbouring k-point i
  real(REAL64),allocatable           ::  rk(:,:)                         !  component in lattice coordinates of the k-point in the mesh
  real(REAL64),allocatable           ::  wgk(:)                          !  weight in the integration of k-point

! allocatable arrays with larger scope

  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi(:,:)                        !  component j of eigenvector i (guess on input)
  complex(REAL64), allocatable       ::  hpsi(:,:)                       !  H | psi> 
  real(REAL64), allocatable          ::  ekpsi(:)                        !  kinetic energy of eigenvector i. (hartree)
  real(REAL64), allocatable          ::  ekpsi_so(:)                     !  kinetic energy of eigenvector i. (hartree)
  real(REAL64), allocatable          ::  vscr(:)                         !  screened potential in the fft real space mesh

  real(REAL64), allocatable          ::  ei_so(:)                        !  spin-orbit eigenvalue (hartree)
  complex(REAL64), allocatable       ::  psi_so(:,:)                     !  component j of eigenvector i (guess on input)

! oscillator strength stuff 

  complex(REAL64), allocatable       ::  h0(:,:)                         !  <Psi|H|Psi> without spin-orbit
  complex(REAL64), allocatable       ::  dh0drk(:,:,:)                   !  d <Psi|H|Psi> d k
  complex(REAL64), allocatable       ::  d2h0drk2(:,:,:,:)               !  d^2 <Psi|H|Psi> d k^2 (not allocated/computed)
  
  complex(REAL64), allocatable       ::  hso0(:,:)                       !  <Psi|H|Psi> with spin-orbit
  complex(REAL64), allocatable       ::  dhso0drk(:,:,:)                 !  d <Psi|H|Psi> d k
  complex(REAL64), allocatable       ::  d2hso0drk2(:,:,:,:)             !  d^2 <Psi|H|Psi> d k^2 (not allocated/computed)

  real(REAL64), allocatable          ::  e_of_k(:,:)                     !  band energies of k-point
  real(REAL64), allocatable          ::  e_of_k_so(:,:)                  !  spin-orbit band energies of k-point

! allocatable single precision arrays

  complex(REAL32), allocatable       ::  dhdrk_32(:,:,:)                 !  d <Psi|H|Psi> d k
  complex(REAL32), allocatable       ::  dhsodrk_32(:,:,:)               !  d <Psi|H|Psi> d k
  
! local variables

  integer                            ::  mxdscr                          !  array dimension for screening potential

  integer                            ::  mxddim                          !  array dimension for the hamiltonian
  integer                            ::  mxdbnd                          !  array dimension for the number of bands
  integer                            ::  mxdwrk                          !  array dimension for fft transform workspace
  integer                            ::  mxdnrk                          !  size of k-points

  integer                            ::  mtxd                            !  dimension of the hamiltonian
  integer                            ::  neig                            !  number of eigenvectors required (maybe modified on output)
  real(REAL64)                       ::  rkpt(3)                         !  j-th component in lattice coordinates of the k-point
  integer                            ::  kmscr(7)                        !  max value of kgv(i,n) used for the potential fft mesh
  integer                            ::  idshift                         !  shift of the fft mesh, used /= 0 only in highly banked memory.

  real(REAL64)                       ::  vmax, vmin                      !  maximum and minimum values of vscr

  integer                            ::  irk,nrka
  character(len=5)                   ::  labelk
  integer                            ::  ipr,nrk,nd
  logical                            ::  lfile
  integer                            ::  nbandi,nx,ny,nz
  real(REAL64)                       ::  sx,sy,sz
  integer                            ::  nsfft(3)

  integer                            ::  io11, io66, io67, io68, iodos   !  tape numbers

  integer                            ::  ncond, nval                     !  number of conduction, valence bands

  integer                            ::  ioerr                           !  error in opening/reading/writing files
  character(len=12)                  ::  filemesh
  character(len=13)                  ::  filedhdrk 
  character(len=16)                  ::  filedhdrkso 
  character(len=12)                  ::  fileband
  character(len=12)                  ::  filedos 

  logical                            ::  lex
  integer                            ::  nocc
  logical                            ::  lmyjob

  integer                            ::  identif                         !  identifier

  logical                            ::  lkpg                            !  If true use the previous G-vectors (same mtxd and isort)
  
  real(REAL64)                       ::  t1, t2
  
  integer                            ::  nder

  integer                            ::  icmax                           !  maximum value of outer iteration
  integer                            ::  ifail                           !  if ifail=0 the ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.

! workers and restart  
  
  integer                            ::  irec_err
  integer                            ::  irk_rd
  integer                            ::  iworker, nworker, cworker
  integer                            ::  ir_size, pp_flag


  io11 = 11

  io66 = 66
  io67 = 67
  io68 = 68
  iodos = 21

  filemesh = 'DOS_MESH.DAT'
  filedhdrk = 'opt_dhdrk.dat'
  filedhdrkso = 'opt_dhdrk_so.dat'
  fileband = 'tmp_band.dat'
  filedos = 'dos_file.dat'

! calculates local potential in fft mesh


  if(flgdal == 'DUAL') then
    kmscr(1) = kmax(1)/2 + 2
    kmscr(2) = kmax(2)/2 + 2
    kmscr(3) = kmax(3)/2 + 2
  else
    kmscr(1) = kmax(1)
    kmscr(2) = kmax(2)
    kmscr(3) = kmax(3)
  endif

  call size_fft(kmscr,nsfft,mxdscr,mxdwrk)

  allocate(vscr(mxdscr))

  ipr = 1

  idshift=0

  call pot_local(ipr, vscr, vmax, vmin, veff, kmscr, idshift,            &
    ng, kgv, phase, conj, ns, inds,                                      &
    mxdscr, mxdgve, mxdnst)

  ipr = 2

  lfile = .false.
  open(unit=io11, file=filemesh, status='old', iostat=ioerr, form = 'formatted')

  if(ioerr == 0) lfile = .true.

  if(lfile) then
    read(io11,*) nbandi, nx,ny,nz, sx,sy,sz

    close(unit=io11)

    mxdbnd = nbandi

  else

    mxdbnd = 2*nint(ztot) + 4

    nbandi = mxdbnd
    nx = 8
    ny = 8
    nz = 8
    sx = 0.0
    sy = 0.0
    sz = 0.0
  endif

 
  call size_mxdnrk(nx,ny,nz, sx,sy,sz, adot, ntrans, mtrx,               &
  mxdnrk)

  allocate(kmap(3,nx,ny,nz))
  allocate(nband(mxdnrk))
  allocate(indk(6,mxdnrk))
  allocate(rk(3,mxdnrk))
  allocate(wgk(mxdnrk))

  call int_pnt(nbandi, nx,ny,nz, sx,sy,sz, ipr,                          &
    adot, ntrans, mtrx,                                                  &
    nrk, rk, wgk, nband, indk, kmap,                                     &
    mxdnrk, mxdbnd)

  neig = nbandi
  
  nval = nint(0.5*ztot)
  ncond = neig - nval

! finds mxddim

  if(lworkers) then
    write(6,*) "Using Workers and Restart Capablities"

    write(6,*) 'how many workers ?'
    read(5,*) nworker

    write(6,*) 'current worker ?'
    read(5,*) iworker
  else
    nworker = 1
    iworker = 1
  endif

  mxddim = 1
  do irk = 1,nrk
    cworker = mod(irk-1,nworker) + 1
    if(cworker== iworker) then
      rkpt(1) = rk(1,irk)
      rkpt(2) = rk(2,irk)
      rkpt(3) = rk(3,irk)
      call size_mtxd(emax,rkpt,adot,ng,kgv,nd)
      if(nd > mxddim) mxddim = nd
    endif
  enddo

! mxddim = int(1.01*mxddim)

! allocates arrays

  allocate(ei(mxdbnd))
  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(psi(mxddim,mxdbnd))
  allocate(hpsi(mxddim,mxdbnd))
  allocate(ekpsi(mxdbnd))
  allocate(ekpsi_so(2*mxdbnd))

  allocate(ei_so(2*mxdbnd))
  allocate(psi_so(2*mxddim,2*mxdbnd))

! loop over k-points

  allocate(h0(mxdbnd,mxdbnd))
  allocate(dh0drk(mxdbnd,mxdbnd,3))

  allocate(hso0(2*mxdbnd,2*mxdbnd))
  allocate(dhso0drk(2*mxdbnd,2*mxdbnd,3))

  allocate(dhdrk_32(neig,neig,3))
  allocate(dhsodrk_32(2*neig,2*neig,3))


! if run by a human do not restart

  if(.not. lworkers) then
    inquire(unit = io66, exist = lex)
    if(lex) close(unit = io66)
    call execute_command_line("rm " // fileband // " 2> /dev/null ")
    inquire(unit = io67, exist = lex)
    if(lex) close(unit = io67)
    call execute_command_line("rm " // filedhdrk // " 2> /dev/null ")
    inquire(unit = io68, exist = lex)
    if(lex) close(unit = io68)
    call execute_command_line("rm " // filedhdrkso // " 2> /dev/null ")
  endif

  inquire(iolength = ir_size) irk_rd, ei(:), ei_so(:)
  open(unit = io66, file = fileband, access="direct", recl=ir_size)

  inquire(iolength = ir_size) irk_rd, dhdrk_32
  open(unit = io67, file = filedhdrk, access="direct", recl=ir_size)

  inquire(iolength = ir_size) irk_rd, dhsodrk_32
  open(unit = io68, file = filedhdrkso, access="direct", recl=ir_size)

  call zeelap(t1)

  do irk=1,nrk

    cworker = mod(irk-1,nworker) +1
    lmyjob = .TRUE.

    if (cworker == iworker) then
      irk_rd = -10
      read(io67,rec=irk, iostat=irec_err) irk_rd

      if(irk_rd ==irk) then
        write(*,'("not computing k-point #",i5)') irk
        lmyjob = .FALSE.
      endif
    else
      write(*,'("not computing k-point #",2i5)') irk, irk_rd
      lmyjob = .FALSE.
    endif

    if(lmyjob) then

      rkpt(1) = rk(1,irk)
      rkpt(2) = rk(2,irk)
      rkpt(3) = rk(3,irk)

      call size_mtxd(emax,rkpt,adot,ng,kgv,nd)

      if(nd > mxddim) then
        write(6,*) "no resize to new mtxd implemented. aborting"         !  should never occur, just paranoid
        stop
      endif

      lkpg = .FALSE.
      ipr = 0

      nocc = neig

      call h_kb_dia_all(diag_type, emax, rkpt, neig, nocc,               &
        flgpsd, ipr, ifail, icmax, iguess, epspsi,                       &
        ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                  &
        sfact, veff, icmplx,                                             &
        nqnl, delqnl, vkb, nkb,                                          &
        ntype, natom, rat, adot,                                         &
        mtxd, hdiag, isort, qmod, ekpg, lkpg,                            &
        psi, hpsi, ei,                                                   &
        vscr, kmscr,                                                     &
        latorb, norbat, nqwf, delqwf, wvfao, lorb,                       &
        mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,          &
        mxdbnd, mxdscr, mxdlao)

      nder = 1
      allocate(d2h0drk2(1,1,3,3))

      call kdotp_matrix(mtxd, neig, psi, ei, rkpt, isort, nder,          &
        h0, dh0drk, d2h0drk2,                                            &
        ng, kgv,                                                         &
        ntype, natom, rat, adot,                                         &
        nqnl, delqnl, vkb, nkb,                                          &
        mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

      deallocate(d2h0drk2)

      call kinetic_energy(neig, mtxd, ekpg, psi, ekpsi,                  &
        mxddim,mxdbnd)

      ipr = 1
      nrka = -1
      call print_eig(ipr, irk, labelk, nrka, rkpt,                       &
        mtxd, icmplx, neig, psi,                                         &
        adot, ei, ekpsi, isort, kgv,                                     &
        mxddim, mxdbnd, mxdgve)

      call spin_orbit_perturb(rkpt, mtxd, isort,                         &
        neig, psi, ei, ei_so, psi_so, .false.,                           &
        ng, kgv,                                                         &
        nqnl, delqnl, vkb, nkb,                                          &
        ntype, natom, rat, adot,                                         &
        mxdtyp,mxdatm,mxdlqp,mxddim,mxdbnd,mxdgve)

      if(ipr == 2) then
        call kinetic_energy_so(neig, mtxd, ekpg, psi_so, ekpsi_so,       &
          mxddim, mxdbnd)
      endif

      call print_eig_so(ipr, irk, labelk, nrka, rkpt,                    &
        mtxd, neig, psi_so,                                              &
        adot, ei_so, ekpsi_so, isort, kgv,                               &
        mxddim, mxdbnd, mxdgve)

      nder = 1
      allocate(d2hso0drk2(1,1,3,3))

      call kdotp_matrix_so_pert(mtxd, neig, psi, ei, rkpt, isort, nder,  &
        hso0, dhso0drk, d2hso0drk2,                                      &
        ng, kgv,                                                         &
        ntype, natom, rat, adot,                                         &
        nqnl, delqnl, vkb, nkb,                                          &
        mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)
     
      call kdotp_matrix_so_convert(neig, hso0, dhso0drk, d2hso0drk2,     &
        nder,                                                            &
        mxdbnd)

      deallocate(d2hso0drk2)

      write(6,'( "iworker #",i5, "   writing in irk # "                  &
              & ,i5, "   of ", i5)') iworker, irk,nrk
      write(6,*)

      dhdrk_32(:,:,:) = dh0drk(:,:,:)
      dhsodrk_32(:,:,:) = dhso0drk(:,:,:)

      write(io66,rec=irk, iostat=irec_err) irk, ei(:), ei_so(:)
      write(io67,rec=irk, iostat=irec_err) irk, dhdrk_32
      write(io68,rec=irk, iostat=irec_err) irk, dhsodrk_32

    endif

  enddo

! checks if all points have been calculated

  do irk = 1, nrk
    irk_rd = -10
    read(io67,rec=irk, iostat=irec_err) irk_rd

    if(irk_rd /= irk) then
      pp_flag = 0

      exit

    endif

    pp_flag = 1

  enddo

  call zeelap(t2)

  write(6,*)
  write(6,'(" elapsed time (s):", 2f14.3)') (t2-t1)


  deallocate(vscr)

  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)

  deallocate(psi)
  deallocate(hpsi)

  deallocate(ekpsi)
  deallocate(ekpsi_so)

  deallocate(psi_so)

  deallocate(h0)
  deallocate(dh0drk)
  deallocate(hso0)
  deallocate(dhso0drk)
  deallocate(dhdrk_32)
  deallocate(dhsodrk_32)

  allocate(e_of_k(neig,nrk))
  allocate(e_of_k_so(2*neig,nrk))

  if (pp_flag ==1) then
    do irk = 1,nrk

      read(io66,rec=irk, iostat=irec_err) irk_rd, ei(:), ei_so(:)

      e_of_k(1:neig,irk) = ei(1:neig)
      e_of_k_so(1:2*neig,irk) = ei_so(1:2*neig)

    enddo
  endif

  deallocate(ei)
  deallocate(ei_so)

  if (pp_flag ==1) then

    write(6,*)
    write(6,*) '  Writing DOS files (and identifier)!'
    write(6,*)

    call get_identifier(identif)


    call out_dos_write(filedos, iodos, title, subtitle,                  &
      .TRUE. , .TRUE. ,identif,                                          &
      nrk, nx, ny, nz, ztot, adot, ntrans, mtrx,                         &
      nband, rk, wgk, indk, kmap, e_of_k, e_of_k_so,                     &
      mxdnrk, mxdbnd)

    irk = nrk + 1
    write(io67,rec=irk, iostat=irec_err) identif
    write(io68,rec=irk, iostat=irec_err) identif

  else
    write(6,*) "this worker is done. run again when all workers are done to see results."
  endif
  
  close(io66)
  close(io67)
  close(io68)
 
! if run by a human clean up temporary files

  if(.not. lworkers) then
    call execute_command_line("rm " // fileband // " 2> /dev/null ")
  endif

  deallocate(e_of_k)
  deallocate(e_of_k_so)

  deallocate(kmap)
  deallocate(nband)
  deallocate(indk)
  deallocate(rk)
  deallocate(wgk)

  return
  end subroutine out_opt
