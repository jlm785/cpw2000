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

!>  Calculates the energies on an uniform grid
!>  for later processing by the density of states or optical program
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.03
!>  \date         8 may 2004, 7 December 2021.
!>  \copyright    GNU Public License v2

  subroutine out_dos(diag_type, lworkers, lproj, lso,                    &
  title, subtitle,                                                       &
  emax, flgdal, flgpsd,                                                  &
  epspsi, icmax, ztot,                                                   &
  adot, ntype, natom, rat, ntrans, mtrx,                                 &
  ng, kgv, phase, conj,                                                  &
  ns, inds, kmax, indv, ek,                                              &
  sfact, icmplx,                                                         &
  veff,                                                                  &
  nqnl, delqnl, vkb, nkb,                                                &
  latorb, norbat, nqwf, delqwf, wvfao, lorb,                             &
  mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)


! version 4.42. 8 may 2004. jlm
! modified 11 february 2008 (read from file). JLM
! modified December 18, 2013 (f90, new interface). JLM
! modified, dimensions vkb, May 23, 2014. JLM
! workers and restart capabilities for version 4.93. CLR
! adpated for version 4.96 March 2020. CLR
! Modified, latorb, 7 June 2020.  JLM
! Modified, lworkers, diag_type, icmax, 8-14 June 2020. JLM
! Modified, mxddim, title, cleaning, 21-23 September 2020. JLM
! Modified, vmax, vmin, 27 November 2020. JLM
! Modified, project wf on atomic orbitals, July 2021. CLR
! Modified, lproj, lso in input. JLM
! Modified, bug correction, Jan 2022. CLR.

! copyright  Jose Luis Martins/Carlos Loia Reis/INESC-MN

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  character(len=4), intent(in)       ::  diag_type                       !<  selects diagonalization, 'pw  ','ao  ','aojc'
  logical, intent(in)                ::  lworkers                        !<  use lworkers
  logical, intent(in)                ::  lproj                           !<  projection in atomic orbitals
  logical, intent(in)                ::  lso                             !<  include spin-orbit in projection

  character(len=50), intent(in)      ::  title                           !<  title for plots
  character(len=140), intent(in)     ::  subtitle                        !<  subtitle for plots

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4)                   ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors
  integer, intent(in)                ::  icmax                           !<  maximum number of iterations for diagonalization
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
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<   KB pseudo.  normalization for atom k, ang. mom. l

  logical, intent(in)                ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  wavefunction for atom k, ang. mom. l (normalized to vcell)
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k

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

  real(REAL64), allocatable          ::  e_of_k(:,:)                     !  band energies of k-point in plot
  real(REAL64), allocatable          ::  e_of_k_so(:,:)                  !  spin-orbit band energies of k-point in plot


! local variables

  integer                            ::  mxdscr                          !  array dimension for screening potential

  integer                            ::  ifail                           !  if ifail=0 h_kb_dia_all was successfull. Otherwise ifail indicates the number of correct digits.
  integer                            ::  iguess                          !  if guess eigenvectors are available, iguess = 1, otherwise iguess = 0

  logical                            ::  lkpg                            !  If true use the previous G-vectors (same mtxd and isort)
  logical                            ::  lpsiso                          !  If true calculates the spin-orbit perturbed wave-functions

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
  integer                            ::  nbandi, nx,ny,nz
  real(REAL64)                       ::  sx,sy,sz
  integer                            ::  nsfft(3)

  integer                            ::  io, io11, io67                  !  tape numbers
  integer                            ::  ioerr

  character(len=12)                  ::  filemesh
  character(len=11)                  ::  filename
  character(len=12)                  ::  filedos

! workers and restart

  integer                            ::  irec_err
  integer                            ::  irk_rd
  integer                            ::  iworker, nworker, cworker
  integer                            ::  ir_size, pp_flag
  logical                            ::  lmyjob

  logical                            ::  lex
  integer                            ::  nocc

  integer                            ::  identif                         !  identifier

  real(REAL64)                       ::  t1, t2

! work in progress stuff clr

! variables for local orbitals

  integer                            ::  mxdorb                          !  array dimension of number of local orbitals
  integer                            ::  nbaslcao                        !  number of atomic orbitals
  complex(REAL64), allocatable       ::  baslcao(:,:)                    !  atomic orbitals in plane-wave basis
  complex(REAL64), allocatable       ::  baslcao_aux(:,:)                !  atomic orbitals in plane-wave basis
  integer, allocatable               ::  infolcao(:,:)                   !  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)
  integer, allocatable               ::  infolcao_aux(:,:)               !  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)
  real(REAL64), allocatable          ::  basxpsi(:,:,:)                  !  |<bas|psi>|^2 for each k
  complex(REAL64), allocatable       ::  prod(:,:)                       !  <bas|psi>

  integer, allocatable               ::  infolcao_so(:,:)                !  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)
  complex(REAL64), allocatable       ::  prod_so(:,:)                    !  <bas|psi>
  real(REAL64), allocatable          ::  basxpsi_so(:,:,:)               !  |<bas|psi>|^2 for each k
  complex(REAL64), allocatable       ::  psi_in(:,:)                     !  component j of eigenvector i (guess on input)


! orbital information with Lowdin orthogonalization

  complex(REAL64), allocatable       ::  S(:,:)
  complex(REAL64), allocatable       ::  S12(:,:)
  complex(REAL64), allocatable       ::  S12_inv(:,:)
  complex(REAL64), allocatable       ::  Swrk(:,:)
  real(REAL64),    allocatable       ::  ev_wrk(:)

!      constants

  real(REAL64), parameter            :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter         :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter         :: C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer      :: n, i, j
  real(REAL64) :: sq


!------------------------------------------------------------------


  io11 = 11
  io67 = 67
  io = 12

  filemesh = 'DOS_MESH.DAT'
  filename = 'dos_rec.dat'
  filedos  = 'dos_file.dat'

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
  idshift = 0

  call pot_local(ipr, vscr, vmax, vmin, veff, kmscr, idshift,            &
  ng, kgv, phase, conj, ns, inds,                                        &
  mxdscr, mxdgve, mxdnst)

  ipr = 2

  lfile = .false.
  open(unit=io11, file=filemesh, status='old', IOSTAT=ioerr, form = 'formatted')
  if(ioerr == 0) lfile = .true.

  if(lfile) then
    read(io11,*) nbandi, nx,ny,nz, sx,sy,sz
    close(unit=io11)
    mxdbnd = nbandi
  else
    mxdbnd = nint(ztot) + 4
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
  adot,                                                                  &
  ntrans, mtrx,                                                          &
  nrk, rk, wgk, nband, indk, kmap,                                       &
  mxdnrk, mxdbnd)

  neig = nbandi

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

! allocates arrays related to atomic orbitals

  if (lproj) then

    call size_nbaslcao(ntype,natom,norbat,lorb,nbaslcao,                   &
    mxdtyp,mxdlao)

    mxdorb = nbaslcao

    allocate(baslcao(mxddim,mxdorb))
    allocate(baslcao_aux(mxddim,mxdorb))
    allocate(infolcao(5,mxdorb))
    allocate(infolcao_aux(5,mxdorb))

    allocate(basxpsi(mxdorb,mxdorb,nrk))  !! check dimensions not consistent with out_band_info_fold :(
!           allocate(basxpsi(mxdorb,mxdbnd,nrk2))

    allocate(S(mxdorb,mxdorb))
    allocate(S12(mxdorb,mxdorb))
    allocate(S12_inv(mxdorb,mxdorb))
    allocate(Swrk(mxdorb,mxdorb))
    allocate(ev_wrk(mxdorb))

    allocate(prod(mxdorb,mxdbnd))

    if (lso) then
      allocate(infolcao_so(5,2*mxdorb))
      allocate(prod_so(2*mxdorb,2*mxdbnd))
      allocate(psi_in(2*mxddim,2*mxdorb))
      allocate(basxpsi_so(2*mxdorb,2*mxdorb,nrk))  !! check dimensions
      lpsiso = .TRUE.
      allocate(psi_so(2*mxddim,2*mxdbnd))
    else
      lpsiso = .FALSE.
      allocate(psi_so(1,1))
    endif    
    else
      lpsiso = .FALSE.
      allocate(psi_so(1,1))
  endif


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

! loop over k-points

  irk=1
  if (lproj) then
    if (lso) then
      inquire(iolength = ir_size) irk_rd, ei(:), ei_so(:), basxpsi_so(:,:,irk)
    else
      inquire(iolength = ir_size) irk_rd, ei(:), ei_so(:), basxpsi(:,:,irk)
    endif
  else
    inquire(iolength = ir_size) irk_rd, ei(:), ei_so(:)
  endif

! if run by a human do not restart

  if(.not. lworkers) then
    inquire(unit = io67, exist = lex)
    if(lex) close(unit = io67)
    call execute_command_line("rm " // filename // " 2> /dev/null ")
  endif

  open(unit = io67, file = filename, access="direct", recl=ir_size)

  pp_flag = 0

  call zeelap(t1)

  do irk=1,nrk

    cworker = mod(irk-1,nworker) + 1
    lmyjob = .TRUE.
    if(cworker == iworker) then
      irk_rd = -10
      read(io67,rec=irk, iostat=irec_err) irk_rd
      if(irk_rd == irk) then
        write(6,'("not computing k-point #",i5)') irk
        lmyjob = .FALSE.
      endif
    else
      write(6,'("not computing k-point #",i5)') irk
      lmyjob = .FALSE.
    endif

    if(lmyjob) then

      rkpt(1) = rk(1,irk)
      rkpt(2) = rk(2,irk)
      rkpt(3) = rk(3,irk)

      lkpg = .FALSE.
      ipr = 0

      nocc = neig
      iguess = 0

      call h_kb_dia_all(diag_type, emax, rkpt, neig, nocc,               &
      flgpsd, ipr, ifail, icmax, iguess, epspsi,                         &
      ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                    &
      sfact, veff, icmplx,                                               &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      mtxd, hdiag, isort, qmod, ekpg, lkpg,                              &
      psi, hpsi, ei,                                                     &
      vscr, kmscr,                                                       &
      latorb, norbat, nqwf, delqwf, wvfao,lorb,                          &
      mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,            &
      mxdbnd, mxdscr, mxdlao)

!     atomic orbital stuff

      if (lproj) then

        call atomic_orbital_c16(rkpt,mtxd,isort,1,                       &
        nbaslcao,baslcao_aux,infolcao,                                   &
        ng,kgv,                                                          &
        norbat,nqwf,delqwf,wvfao,lorb,                                   &
        ntype,natom,rat,adot,                                            &
        mxdtyp,mxdatm,mxdlqp,mxddim,mxdorb,mxdgve,mxdlao)

        call zgemm('C','N',nbaslcao,nbaslcao,mtxd,C_UM,baslcao_aux,      &
        mxddim,baslcao_aux,mxddim,C_ZERO,S,mxdorb)

        call  GetS12(S,S12,S12_inv,Swrk,ev_wrk,nbaslcao)

        call zgemm('n','n',mtxd,nbaslcao,nbaslcao,C_UM,baslcao_aux,      &
        mxddim, S12,mxdorb,C_ZERO,baslcao,mxddim)

        call zgemm('C','N',nbaslcao,nbaslcao,mtxd,C_UM,baslcao,          &
        mxddim, baslcao,mxddim,C_ZERO,S,mxdorb)

        do n = 1,nbaslcao
          sq = ZERO
          do j = 1,mtxd
            sq = sq + real(baslcao(j,n)*conjg(baslcao(j,n)),REAL64)
          enddo
          sq = UM/sqrt(sq)
          do j = 1,mtxd
            baslcao(j,n) = sq * baslcao(j,n)
          enddo
        enddo

        call zgemm('C','N',nbaslcao,neig,mtxd,C_UM,baslcao,mxddim,       &
        psi,mxddim,C_ZERO,prod,mxdorb)

!!!       needs to be done with worker data when tested

        do n = 1,neig
        do j = 1,nbaslcao
          basxpsi(j,n,irk) = real(prod(j,n)*conjg(prod(j,n)),REAL64)
        enddo
        enddo
!!!
      endif

      call kinetic_energy(neig,mtxd,ekpg,psi,ekpsi,                      &
      mxddim,mxdbnd)

      ipr = 1
      nrka = -1
      call print_eig(ipr,irk,labelk,nrka,rkpt,                           &
      mtxd,icmplx,neig,psi,                                              &
      adot,ei,ekpsi,isort,kgv,                                           &
      mxddim,mxdbnd,mxdgve)

      call spin_orbit_perturb(rkpt,mtxd,isort,                           &
      neig,psi,ei,ei_so,psi_so,lpsiso,                                   &
      ng,kgv,                                                            &
      nqnl,delqnl,vkb,nkb,                                               &
      ntype,natom,rat,adot,                                              &
      mxdtyp,mxdatm,mxdlqp,mxddim,mxdbnd,mxdgve)


!    new adapted from ao_h_and_s_spin_orbit

      if (lproj .and. lso) then
        do i=1,nbaslcao
        infolcao_so(:,2*i)   = infolcao(:,i)
        infolcao_so(:,2*i-1) = infolcao(:,i)
        enddo


        do i=1,nbaslcao
        do j=1,mtxd
          psi_in(2*j-1,2*i         ) = baslcao(j,i)
          psi_in(2*j  ,2*i         ) = C_ZERO
          psi_in(2*j-1,2*i-1       ) = C_ZERO
          psi_in(2*j  ,2*i-1       ) = baslcao(j,i)
        enddo
        enddo

        call zgemm('c','n',2*nbaslcao,2*neig,2*mtxd,                       &
        C_UM,psi_in,2*mxddim,psi_so,2*mxddim,C_ZERO,prod_so,2*mxdorb)

!!!       needs to be done with worker data when tested

        do n = 1,2*neig
        do j = 1,2*nbaslcao
          basxpsi_so(j,n,irk) = real(prod_so(j,n)*conjg(prod_so(j,n)),     &
          REAL64)
        enddo
        enddo
!!!
      endif


      if(ipr == 2) then
        call kinetic_energy_so(neig,mtxd,ekpg,psi_so,ekpsi_so,           &
        mxddim,mxdbnd)
      endif

      call print_eig_so(ipr,irk,labelk,nrka,rkpt,                        &
      mtxd,neig,psi_so,                                                  &
      adot,ei_so,ekpsi_so,isort,kgv,                                     &
      mxddim,mxdbnd,mxdgve)


      write(6,'( "iworker #",i5, "   writing in irk # "                  &
                 ,i5, "   of ", i5)') iworker, irk,nrk
      write(6,*)

      if (lproj) then
        if (lso) then
          write(io67, rec = irk) irk, ei(:), ei_so(:), basxpsi_so(:,:,irk)
        else
          write(io67, rec = irk) irk, ei(:), ei_so(:), basxpsi(:,:,irk)
        endif
      else
        write(io67, rec = irk) irk, ei(:), ei_so(:)
      endif

    endif

  enddo

! end loop over k-points

  do irk = 1, nrk
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

  allocate(e_of_k(neig,nrk))
  allocate(e_of_k_so(2*neig,nrk))

  if (pp_flag ==1) then
    do irk = 1,nrk

      if (lproj) then
        if (lso) then
          read(io67,rec=irk, iostat=irec_err) irk_rd, ei(:), ei_so(:), basxpsi_so(:,:,irk)
        else
          read(io67,rec=irk, iostat=irec_err) irk_rd, ei(:), ei_so(:), basxpsi(:,:,irk)
        endif
      else
        read(io67,rec=irk, iostat=irec_err) irk_rd, ei(:), ei_so(:)
      endif

      e_of_k(1:neig,irk) = ei(1:neig)
      e_of_k_so(1:2*neig,irk) = ei_so(1:2*neig)

    enddo
  endif

  deallocate(ei)
  deallocate(ei_so)

  if(pp_flag == 1) then

    write(6,*)
    write(6,*) '  Generating DOS files!'
    write(6,*)

    identif = 0

    call out_dos_write(filedos, io, title, subtitle,                     &
      .TRUE. , .TRUE. ,identif,                                          &
      nrk, nx, ny, nz, ztot, adot, ntrans, mtrx,                         &
      nband, rk, wgk, indk, kmap, e_of_k, e_of_k_so,                     &
      mxdnrk, mxdbnd)

!  this needs to be done differently .

    if (lproj) then
      open(unit=190,file="dos_file_proj.dat", form="unformatted")
      if (lso) then
        write(190) nrk
        write(190) 2*neig
        write(190) 2*nbaslcao
        do irk=1,nrk
!          write(*,*) basxpsi_so(:,:,irk)
!          stop
          write(190) basxpsi_so(:,:,irk)
        enddo
      else
        write(190) nrk
        write(190) neig
        write(190) nbaslcao
        do irk=1,nrk
!          write(*,*) basxpsi(:,:,irk)
!          stop
          write(190) basxpsi(:,:,irk)
        enddo
      endif
      close(unit=190)
    endif

  else
    close(io67)
    write(6,*) "this worker is done. run again when all workers are done to see results."
  endif

! if run by a human clean up temporary files

  if(.not. lworkers) then
    call execute_command_line("rm " // filename // " 2> /dev/null ")
  endif

  deallocate(kmap)
  deallocate(nband)
  deallocate(indk)
  deallocate(rk)
  deallocate(wgk)

  deallocate(e_of_k)
  deallocate(e_of_k_so)
  
  if (lproj) then

    deallocate(baslcao)
    deallocate(baslcao_aux)
    deallocate(infolcao)
    deallocate(infolcao_aux)

    deallocate(basxpsi)

    deallocate(S)
    deallocate(S12)
    deallocate(S12_inv)
    deallocate(Swrk)
    deallocate(ev_wrk)

    deallocate(prod)

    if (lso) then
      deallocate(infolcao_so)
      deallocate(prod_so)
      deallocate(psi_in)
      deallocate(basxpsi_so)  
      deallocate(psi_so)
    else
      deallocate(psi_so)
    endif    
    else
      deallocate(psi_so)
  endif
  
  
  
  

  return
  end subroutine out_dos
