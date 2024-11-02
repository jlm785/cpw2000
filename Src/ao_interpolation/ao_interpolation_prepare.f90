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

!>  Calculates the Hamiltonian and overlap for a grid of k-points in
!>  the Brillouin zone with non-orthogonal atomic orbitals.
!>
!>  Having computed H(k) and S(k) these quantities are dft transformed
!>  to a Wigner-Seitz cell in direct space H(R), S(R).
!>  With the latter it is then possible to obtain H(k') and S(k') in
!>  an arbitrary k' point with idft.
!>
!>  \author       Carlos Loia Reis
!>  \version      5.11
!>  \date         2014. 1 November 2024.
!>  \copyright    GNU Public License v2


subroutine ao_interpolation_prepare(ioreplay, noiData,                   &
        emax, flgdal, flgpsd, iguess, epspsi, icmax, ztot,               &
        adot, ntype, natom, rat, ntrans, mtrx,                           &
        ng, kgv, phase, conj,                                            &
        ns, inds, kmax, indv, ek,                                        &
        sfact, icmplx,                                                   &
        veff,                                                            &
        nqnl,delqnl, vkb, nkb,                                           &
        norbat, nqwf, delqwf, wvfao, lorb,                               &
        mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)

! Written by Carlos Loia Reis adapting previous code.
! First stable version in 4.54.1 2014.
! Modified 25 November 2015. JLM
! Modified, mxdlao, December 1 2015.  JLM.
! Calculation and Interpolation of Hamiltonian derivatives and
! spectral weights in unfolding. November 2018. CLR.
! Fixed a major bug in IrredBz and
! cleaned up unfolding routine, October 2019. CLR.
! Modified Wigner-Seitz cell routine, October 2019. CLR.
! Modified May 2020, for 4.96. CLR
! Modified documentation May 2020, nocc, icmax, June 12 2020.jlm
! Modified, vmax, vmin, 27 November 2020. JLM
! Modified, indentation, more comments. 2 October 2024. JLM
! Modified, removed IrredBZM, only used here, 28 October 2024. JLM
! Array dimensions, ao_h_and_s API, eigenvalue mapping. 1 November 2024. JLM


  use NonOrthoInterp

  implicit none

  integer, parameter            :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

  type(noiData_t)                    ::  noiData                         !<  see NonOrthoInterp

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4)                   ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  integer, intent(in)                ::  iguess                          !<  if guess eigenvectors are available, iguess = 1, otherwise iguess = 0
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

  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  (1/q**l) * wavefunction for atom k, ang. mom. l (unnormalized to vcell)

! allocatable arrays with larger scope

  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi(:,:)                        !  component j of eigenvector i (guess on input)
  complex(REAL64), allocatable       ::  hpsi(:,:)                       !  component j of eigenvector i (guess on input)

  real(REAL64), allocatable          ::  ekpsi(:)                        !  kinetic energy of eigenvector i. (hartree)
  real(REAL64), allocatable          ::  vscr(:)                         !  screened potential in the fft real space mesh

  real(REAL64), allocatable          ::  ei_so(:)                        !  spin-orbit eigenvalue (hartree)
  complex(REAL64), allocatable       ::  psi_so(:,:)                     !  component j of eigenvector i (guess on input)

  real(REAL64), allocatable          ::  ei_ao(:)                        !  atomic orbital eigenvalue (hartree)
  complex(REAL64), allocatable       ::  psi_ao(:,:)                     !  component j of atomic orbital eigenvector i (guess on input)
  complex(REAL64), allocatable       ::  hpsi_ao(:,:)                    !  H | psi >

  real(REAL64), allocatable          ::  ev_pw(:,:), ev_pw_irred(:,:), ev_interp(:)
  complex(REAL64), allocatable       ::  Hao(:,:), S(:,:), dh0drk_AO(:,:,:)

  real(REAL64), allocatable:: rk_grid(:,:)

  integer, allocatable               ::  imatch(:)                       !  points to the old state that is most similar to the new state
  real(REAL64), allocatable          ::  xover(:)                        !  overlap of new state with old state

! local variables

  integer                            ::  mxdorb                          !  array dimension for orbitals
  integer                            ::  mxdbnd                          !  array dimension for pw bands
  integer                            ::  mxdscr                          !  array dimension for screening potential

  integer                            ::  mxddim                          !  array dimension for the hamiltonian
  integer                            ::  mxdwrk                          !  array dimension for fft transform workspace

  integer                            ::  mtxd                            !  dimension of the hamiltonian
  real(REAL64)                       ::  rkpt(3)                         !  j-th component in lattice coordinates of the k-point
  integer                            ::  kmscr(7)                        !  max value of kgv(i,n) used for the potential fft mesh
  integer                            ::  idshift                         !  shift of the fft mesh, used /= 0 only in highly banked memory.

  real(REAL64)                       ::  vmax, vmin                      !  maximum and minimum values of vscr

  integer                            ::  irk
  integer                            ::  idk

  integer                            ::  ipr,nd
  integer                            ::  nsfft(3)
  integer                            ::  iband, jband, idir

  integer                            ::  ifail                           !  if ifail=0 the ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.

  logical                            ::  lkpg                            !  If true use the previous G-vectors (same mtxd and isort)
  logical                            ::  lpsiso                          !  If true calculates the spin-orbit perturbed wave-functions


!--------------------------------------------------------------------
  integer           ::  nk1,nk2,nk3
  integer           ::  ws_n1, ws_n2, ws_n3

  integer           ::  nkpt, lmax
  integer           ::  norb_ao                                          !  number of orbitals for ao diagonalization
  integer           ::  neig                                             !  number of orbitals for pw diagonalization
  integer           ::  nocc

! Arrays replacing IrredBZ

  integer, allocatable               ::  kmap(:,:,:,:)                   !  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
  integer, allocatable               ::  wrk_nband(:)                    !  number of bands for each k-points
  integer, allocatable               ::  indk(:,:)                       !  index of the six k-points neighbouring k-point i
  real(REAL64), allocatable          ::  rk_irred(:,:)                   !  component in lattice coordinates of the k-point in the mesh
  real(REAL64), allocatable          ::  w_mesh(:)                       !  weight in the integration of k-point
  integer                            ::  nrk_irred                       !  number of irreducible k-points


  integer           ::  mxdpnt

  logical           ::  lso                                              !  if true uses spin-orbit
  logical           ::  loptical                                         !  if true calculates optical matrix elements
  integer           ::  nspin                                            !  1 for no spin-orbit, 2 for spin orbit calculation

  character(len=1)  ::  yesno

  integer           ::  inoso, irem

! parameters

  real(REAL64), parameter :: ZERO = 0.0_REAL64
  real(REAL64), parameter :: UM = 1.0_REAL64

! counters

  integer           ::  i, j, k
  integer           ::  ix,iy,iz


!--------------------------------------------------------------------
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

  call size_fft(kmscr, nsfft, mxdscr, mxdwrk)

  allocate(vscr(mxdscr))

  ipr = 1

  idshift = 0
  call pot_local(ipr, vscr, vmax, vmin, veff, kmscr, idshift,            &
      ng, kgv, phase, conj, ns, inds,                                    &
      mxdscr, mxdgve, mxdnst)

  ipr = 2


!-------------------------------------------------------------------------
! initializes interpolation routines

  write(6,*)
  write(6,*) '  initializing ao interpolation'
  write(6,*)

  write(6,*) '  please enter grid of k-points'
  write(6,*) '  recomended is 4 4 4 minimum'

  read(5,*) nk1,nk2,nk3    ! this is to be removed
  write(ioreplay,*) nk1,nk2,nk3,'     k-point grid'

  write(6,'( "  interpolation in a mesh of", i3," x ",i3," x ",i3, " k-points"  )') nk1,nk2,nk3

  write(6,*) '  do you want to modify Wigner-Seitz cell search size ?'
  read(5,*) yesno
  write(ioreplay,*) yesno,'     ws_search_size_yesno'

  if(yesno == 'y' .or. yesno == 'Y') then
     write(6,*) '  please enter Wigner Seitz cell search size'
     write(6,*) '  recomended is 2 2 2 minimum'

     read(5,*) ws_n1,ws_n2,ws_n3
     write(ioreplay,*) ws_n1,ws_n2,ws_n3,'     ws_search_size'
  else
   ws_n1 = 2
   ws_n2 = 2
   ws_n3 = 2
  endif

! Allocates the arrays and calls int_pnt
! It is uniform nk1 nk2 nk3 grid with zero shift.

  mxdpnt = nk1*nk2*nk3

  allocate(kmap(3,nk1,nk2,nk3))
  allocate(wrk_nband(mxdpnt))
  allocate(indk(6,mxdpnt))
  allocate(rk_irred(3,mxdpnt))
  allocate(w_mesh(mxdpnt))

  call int_pnt(1, nk1,nk2,nk3, ZERO,ZERO,ZERO, 2,                        &
      adot,                                                              &
      ntrans, mtrx,                                                      &
      nrk_irred, rk_irred, w_mesh, wrk_nband, indk, kmap,                &
      mxdpnt, 1)

  deallocate(wrk_nband)
  deallocate(indk)
  deallocate(w_mesh)

!----------------------------------------------------

  write(6,*)
  write(6,*) '  do you want a calculation including Spin-Orbit (y/n) ?'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'     spin-orbit'

  if(yesno == 'y' .or. yesno == 'Y') then
    lso = .TRUE.
  else
    lso = .FALSE.
  endif

  if(lso) then
    write(6,*)
    write(6,*) '  calculation will be done using SO perturbation.'
    write(6,*)
  endif

!----------------------------------------------------

  write(6,*)
  write(6,*) '  do you want a calculation including Hamiltonian Derivatives (y/n) ?'
  write(6,*) '  this will be needed to compute the dielectric function '
  write(6,*) '  using interpolation (it is memory intensive...).'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'     optical'

  if(yesno == 'y' .or. yesno == 'Y') then
    loptical = .TRUE.
  else
    loptical = .FALSE.
  endif
  if(loptical) then
    write(6,*)
    write(6,*) '  calculation will also include dh0drk computation.'
    write(6,*)
  endif

! allocates and fills the uniform grid

  allocate(rk_grid(3,nk1*nk2*nk3))

  irk = 1
    do ix = 0,nk1-1
      do iy = 0,nk2-1
        do iz = 0,nk3-1
        rk_grid(1,irk) = (UM*ix)/nk1
        rk_grid(2,irk) = (UM*iy)/nk2
        rk_grid(3,irk) = (UM*iz)/nk3
        irk = irk+1
      enddo
    enddo
  enddo

!--------------------------------------------------------------------

! finds mxddim

  mxddim = 1
  do irk = 1,nrk_irred
    rkpt(1) = rk_irred(1,irk)
    rkpt(2) = rk_irred(2,irk)
    rkpt(3) = rk_irred(3,irk)
    call size_mtxd(emax, rkpt, adot, ng, kgv, nd)
    if(nd > mxddim) mxddim = nd
  enddo

!--------------------------------------------------------------------

! finds initial value of norb_ao

  norb_ao = 0
  lmax = 0
  do k = 1,ntype
    do j = 1,norbat(k)
      norb_ao = norb_ao + (2*lorb(j,k)+1)*natom(k)
      lmax = max(lorb(j,k),lmax)
    enddo
  enddo

  nkpt = nk1*nk2*nk3

! Finds neig, mxdorb, mxdbnd  THE RECIPE OF neig COULD BE IMPROVED

  neig = min(int(ztot)+4,norb_ao+10)

  mxdorb = norb_ao
  mxdbnd = neig

  write(6,*)
  write(6,*) '  number k-points:', nkpt
  write(6,*) '  number of atomic orbitals:', norb_ao
  write(6,*) '  number of eigenvalues calculated with plane waves:', neig
  write(6,*) '  total charge:', nint(ztot)
  write(6,*)

! allocates arrays

  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(ekpsi(mxdorb))

  allocate(ei(mxdbnd))
  allocate(psi(mxddim,mxdbnd))
  allocate(hpsi(mxddim,mxdbnd))

  allocate(ei_ao(mxdorb))
  allocate(psi_ao(mxddim,mxdorb))
  allocate(hpsi_ao(mxddim,mxdorb))

  if(lso) then
    allocate(ei_so(2*mxdorb))
    allocate(psi_so(2*mxddim,2*mxdorb))
  endif

!--------------------------------------------------------------------

! Initialize structures and more allocations

  if(lso) then
    nspin = 2
  else
    nspin = 1
  endif

  call NonOrthoInterpInit(noiData, lso, loptical, nk1,nk2,nk3,           &
      ws_n1,ws_n2,ws_n3, adot, nkpt, nspin*norb_ao, nspin*neig)

  allocate(ev_pw(nspin*norb_ao,nkpt))
  allocate(ev_pw_irred(nspin*norb_ao,nkpt))
  allocate(ev_interp(nspin*norb_ao))
  allocate(Hao(nspin*norb_ao,nspin*norb_ao))
  allocate(S(nspin*norb_ao,nspin*norb_ao))
  allocate(dh0drk_AO(nspin*norb_ao,nspin*norb_ao,3))

  allocate(imatch(mxdorb))
  allocate(xover(mxdorb))



!-------------------------------------------------------------------------

  if (lso) then

    write(6,*)
    write(6,*) '  Computing full precision eigenvalues ev(i) including Spin-Orbit'
    write(6,*)

  else

    write(6,*)
    write(6,*) '  Computing full precision eigenvalues ev(i) without spin orbit'
    write(6,*)

  endif

  do irk = 1,nrk_irred

    rkpt(1) = rk_irred(1,irk)
    rkpt(2) = rk_irred(2,irk)
    rkpt(3) = rk_irred(3,irk)

    write(6,'(2i5,3f8.5)') nrk_irred, irk, rkpt(1),rkpt(2),rkpt(3)

    lkpg = .FALSE.
    ipr = 0

    nocc = norb_ao

!   THERE IS A LOT OF DUPLICATION HERE

    call h_kb_dia_all('ao  ', emax, rkpt, norb_ao, nocc, flgpsd,         &
        ipr, ifail, icmax, iguess, epspsi,                               &
        ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                  &
        sfact, veff, icmplx,                                             &
        nqnl, delqnl, vkb, nkb,                                          &
        ntype, natom, rat, adot,                                         &
        mtxd, hdiag, isort, qmod, ekpg, .FALSE.,                         &
        psi_ao, hpsi_ao, ei_ao,                                          &
        vscr, kmscr,                                                     &
        .TRUE., norbat, nqwf, delqwf, wvfao, lorb,                       &
        mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,          &
        mxdorb, mxdscr, mxdlao)

    nocc = neig

    call h_kb_dia_all('pw  ', emax, rkpt, neig, nocc, flgpsd,            &
        ipr, ifail, icmax, iguess, epspsi,                               &
        ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                  &
        sfact, veff, icmplx,                                             &
        nqnl, delqnl, vkb, nkb,                                          &
        ntype, natom, rat, adot,                                         &
        mtxd, hdiag, isort, qmod, ekpg, .FALSE.,                         &
        psi, hpsi, ei,                                                   &
        vscr, kmscr,                                                     &
        .TRUE., norbat, nqwf, delqwf, wvfao, lorb,                       &
        mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,          &
        mxdbnd, mxdscr, mxdlao)

    call ao_match_state(imatch, xover, mtxd,                             &
        neig, psi, ei, norb_ao, psi_ao, ei_ao,                           &
        mxddim, mxdbnd, mxdorb)


    if(lso) then

      lpsiso = .false.

      call spin_orbit_perturb(rkpt, mtxd, isort,                         &
          neig, psi, ei, ei_so, psi_so, lpsiso,                          &
          ng, kgv,                                                       &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdorb, mxdgve)

      do i = 1,2*norb_ao
        inoso = (i+1)/2
        irem = mod(i+1,2)
        if(imatch(inoso) > 0) then
          ev_pw_irred(i,irk) = ei_so(2*imatch(inoso)+irem-1)
        else
          ev_pw_irred(i,irk) = ei_ao(inoso)
        endif
      enddo

    else

      do i = 1,norb_ao
        if(imatch(i) > 0) then
          ev_pw_irred(i,irk) = ei(imatch(i))
        else
          ev_pw_irred(i,irk) = ei_ao(i)
        endif
      enddo

    endif

  enddo

! now unpack the eigenvalues

  if(lso) then

    irk = 1
    do i = 1,nk1
    do j = 1,nk2
    do k = 1,nk3
      idk = iabs(kmap(1,i,j,k))
      do iband = 1,2*norb_ao
        ev_pw(iband,irk) = ev_pw_irred(iband,idk)
      enddo
      irk = irk+1
    enddo
    enddo
    enddo

  else

    irk = 1
    do i = 1,nk1
    do j = 1,nk2
    do k = 1,nk3
      idk = iabs(kmap(1,i,j,k))
      do iband = 1,norb_ao
        ev_pw(iband,irk) = ev_pw_irred(iband,idk)
      enddo
      irk = irk+1
    enddo
    enddo
    enddo

  endif

  deallocate(kmap)
  deallocate(rk_irred)


!-------------------------------------------------------------------------

! Now let's get to business...

  if (lso) then

    write(6,*)
    write(6,*) '  Computing atomic orbitals: |psi> and overlaps: S=<psi|psi>, H=<psi|H|psi>'
    write(6,*) '  SO version'
    write(6,*)

    do irk=1, nkpt

      rkpt(1) = rk_grid(1,irk)
      rkpt(2) = rk_grid(2,irk)
      rkpt(3) = rk_grid(3,irk)

      write(*,'(2i5,3f8.5)') nkpt,irk, rkpt(1),rkpt(2),rkpt(3)

      call ao_h_and_s_spin_orbit(emax, rkpt, norb_ao, flgpsd,            &
          ng, kgv,                                                       &
          veff, icmplx,                                                  &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
!          mtxd, hdiag, isort, qmod, ekpg,                                &
          norbat, nqwf, delqwf, wvfao, lorb,                             &
          psi_ao, hpsi_ao,                                               &
          Hao, S, dh0drk_AO,                                             &
          vscr, kmscr,                                                   &
          mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxddim, mxdorb,        &
          mxdscr, mxdlao)

      call NonOrthoInterpSetGridData(noiData, irk, Hao, S, ev_pw(:,irk))

      if (loptical) then

        do idir=1, 3
          do iband=1, 2*norb_ao
          do jband=1, 2*norb_ao
            noiData%fiData%dh0drk_GridK(iband,jband,idir,irk) = dh0drk_AO(iband,jband,idir)
          enddo
          enddo
        enddo

      endif

    enddo

  else

    write(6,*)
    write(6,*) '  Computing atomic orbitals: |psi> and overlaps: S=<psi|psi>, H=<psi|H|psi>'
    write(6,*)

    do irk=1, nkpt

      rkpt(1) = rk_grid(1,irk)
      rkpt(2) = rk_grid(2,irk)
      rkpt(3) = rk_grid(3,irk)

      write(*,'(2i5,3f8.5)') nkpt,irk, rkpt(1),rkpt(2),rkpt(3)

!     hdiag etc... are not needed here.  Should simplify call

      call ao_h_and_s(emax, rkpt, norb_ao, flgpsd,                        &
          ng, kgv,                                                        &
          veff, icmplx,                                                   &
          nqnl, delqnl, vkb, nkb,                                         &
          ntype, natom, rat, adot,                                        &
          norbat, nqwf, delqwf, wvfao, lorb,                              &
          psi_ao, hpsi_ao,                                                &
          Hao, S, dh0drk_AO,                                              &
          vscr, kmscr,                                                    &
          mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxddim, mxdorb,         &
          mxdscr, mxdlao)

      call NonOrthoInterpSetGridData(noiData, irk, Hao, S, ev_pw(:,irk))

      if (loptical) then

        do idir=1, 3
          do iband=1, norb_ao
          do jband=1, norb_ao
            noiData%fiData%dh0drk_GridK(iband,jband,idir,irk) = dh0drk_AO(iband,jband,idir)
          enddo
          enddo
        enddo

      endif

    enddo

  endif

!-------------------------------------------------------------------------

  write(6,*)
  write(6,*) '  Finishing Setup'
  call NonOrthoInterpFinishSetup(noiData)
  write(6,*) '  Writing to file for later use.'
  call NonOrthoInterpWriteToFile(noiData)
  write(6,*)

  deallocate(vscr)
  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)
  deallocate(ekpsi)

  deallocate(ei)
  deallocate(psi)
  deallocate(hpsi)

  deallocate(ei_ao)
  deallocate(psi_ao)
  deallocate(hpsi_ao)

  deallocate(imatch)
  deallocate(xover)

  if(lso) then
    deallocate(ei_so)
    deallocate(psi_so)
  endif

  deallocate(rk_grid)

  deallocate(ev_pw)
  deallocate(ev_pw_irred)
  deallocate(ev_interp)
  deallocate(Hao)
  deallocate(S)
  deallocate(dh0drk_AO)


  return

end subroutine ao_interpolation_prepare
