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

!>  This subroutine calculates the band structure along a path
!>  It uses the kdotp aproximation
!>  Files band.agr, band_so.agr, band.gp and band_so.gp,
!>  for later ploting with gnuplot and xmgrace are written
!>  Circuit for band structure is defined in BAND_LINES.DAT
!>  It also writes kdotp_matrix.dat and kdotp_matrix_so.dat for further processing
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         8 may 2004, 13 August 2025.
!>  \copyright    GNU Public License v2

subroutine out_band_kdotp_2nd(title, subtitle,                           &
      emax, flgdal, flgpsd, epspsi, icmax, ztot, efermi,                 &
      adot, ntype, natom, rat,                                           &
      ng, kgv, phase, conj,                                              &
      ns, inds, kmax, indv, ek,                                          &
      sfact, icmplx,                                                     &
      veff,                                                              &
      nqnl, delqnl, vkb, nkb,                                            &
      latorb, norbat, nqwf, delqwf, wvfao, lorb,                         &
      mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)

! version 4.42. 8 may 2004. jlm
! modified 11 february 2008 (read from file). JLM
! modified September 5, 2012 (f90, circuit). jlm
! modified October 18, 2013 (band index permutations). jlm
! modified January 8, 2013 (new interface). jlm
! modified, vkb dimensions, March 31, 2014. jlm
! modified (title,subtitle) xmgrace plot. 4 August 2012. JLM
! modified 1 August 2019. Introduce changes in out_band.f90. JLM
! Modified February 2020. Documentation, ifail, icmax,
! Modified h_kb_dia_all, 8 June 2020, icmax 14 June 2020. JLM
! Modified, vmax, vmin, 27 November 2020. JLM
! Modified, efermi, 29 November 2021. JLM
! Modified, iguess, indentation, 11 November 2023. JLM
! Modified, ztot in out_band_circuit_size. 26 July 2024. JLM
! Modified, rk in out_band_eref, 13 August 2025. JLM

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

  character(len=50), intent(in)      ::  title                           !<  title for plots
  character(len=140), intent(in)     ::  subtitle                        !<  subtitle for plots

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4), intent(in)       ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors
  integer, intent(in)                ::  icmax                           !<  maximum number of iterations for diagonalization
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  efermi                          !  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

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
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

  logical, intent(in)                ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  wavefunction for atom k, ang. mom. l
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k

! allocatable arrays for Brillouin zone path

  integer                            ::  nlines                          !  number of lines in reciprocal space
  integer, allocatable               ::  nkstep(:)                       !  number of steps in line
  logical, allocatable               ::  ljump(:)                        !  indicates if the new line contains a jump from the preceeding
  integer                            ::  nvert                           !  number of vertical lines in plot
  real(REAL64), allocatable          ::  xcvert(:)                       !  x coordinate of vertical line
  real(REAL64), allocatable          ::  xk(:)                           !  x coordinate of k-point in plot
  real(REAL64), allocatable          ::  rk(:,:)                         !  x coordinate of k-point in plot
  real(REAL64), allocatable          ::  e_of_k(:,:)                     !  band energies of k-point in plot
  real(REAL64), allocatable          ::  e_of_k_so(:,:)                  !  spin-orbit band energies of k-point in plot
  character(len=6), allocatable      ::  label(:)                        !  label of symmetry k-points
  real(REAL64), allocatable          ::  xklab(:)                        !  x coordinate of label

! variables for match_state

  complex(REAL64),  allocatable      ::  psiold(:,:)                     !  old eigenvectors
  integer, allocatable               ::  imatch(:)                       !  new matching
  integer, allocatable               ::  isold(:)                        !  old isort
  integer, allocatable               ::  iperm(:), ipermold(:)           !  permutation of the eigenvalues

  complex(REAL64),  allocatable      ::  psiold_so(:,:)                  !  old eigenvectors
  integer, allocatable               ::  imatch_so(:)                    !  new matching
  integer, allocatable               ::  iperm_so(:), ipermold_so(:)     !  permutation of the eigenvalues

! allocatable arrays with larger scope

  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  ei_so(:)                        !  spin-orbit eigenvalue (hartree)
  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  integer, allocatable               ::  isort_so(:)                     !  g-vector associated with row/column i of spin-orbit hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi(:,:)                        !  component j of eigenvector i
  real(REAL64), allocatable          ::  ekpsi(:)                        !  kinetic energy of eigenvector i. (hartree)

  real(REAL64), allocatable          ::  ei0(:)                          !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  hdiag0(:)                       !  hamiltonian diagonal
  integer, allocatable               ::  isort0(:)                       !  g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod0(:)                        !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg0(:)                        !  kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi0(:,:)                       !  component j of eigenvector i
  complex(REAL64), allocatable       ::  hpsi0(:,:)                      !  H | psi > for component j of eigenvector i
  real(REAL64), allocatable          ::  ekpsi0(:)                       !  kinetic energy of eigenvector i. (hartree)

  real(REAL64), allocatable          ::  vscr(:)                         !  screened potential in the fft real space mesh
  complex(REAL64), allocatable       ::  psi_so(:,:)                     !  component j of eigenvector i (guess on input)
  complex(REAL64), allocatable       ::  psi_so0(:,:)                    !  component j of eigenvector i (guess on input)
  real(REAL64), allocatable          ::  ekpsi_so(:)                     !  kinetic energy of eigenvector i. (hartree)

! allocatable arrays for k.p


  complex(REAL64), allocatable       ::  h0(:,:)                         ! <Psi|H|Psi> without spin-orbit
  complex(REAL64), allocatable       ::  dh0drk(:,:,:)                   ! d <Psi|H|Psi> d k
  complex(REAL64), allocatable       ::  d2h0drk2(:,:,:,:)               ! d^2 <Psi|H|Psi> d k^2

  complex(REAL64), allocatable       ::  hso0(:,:)                       ! <Psi|H_so|Psi> with spin-orbit
  complex(REAL64), allocatable       ::  dhso0drk(:,:,:)                 ! d <Psi|H_so|Psi> d k
  complex(REAL64), allocatable       ::  d2hso0drk2(:,:,:,:)             ! d^2 <Psi|H_so|Psi> d k^2

! local variables

  integer           ::  mxdscr         !  array dimension for screening potential

  integer           ::  iguess         !  if guess eigenvectors are available, iguess = 1, otherwise iguess = 0
  integer           ::  ifail          !  if ifail=0 the ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.

  integer           ::  mxddim         !  array dimension for the hamiltonian
  integer           ::  mxdbnd         !  array dimension for the number of bands
  integer           ::  mxdwrk         !  array dimension for fft transform workspace

  integer           ::  mtxd           !  dimension of the hamiltonian
  integer           ::  mtxd0          !  dimension of the hamiltonian
  integer           ::  neig           !  number of eigenvectors required (maybe modified on output)
  real(REAL64)      ::  rkpt(3)        !  j-th component in lattice coordinates of the k-point
  integer           ::  kmscr(7)       !  max value of kgv(i,n) used for the potential fft mesh
  integer           ::  idshift        !  shift of the fft mesh, used /= 0 only in highly banked memory.

  real(REAL64)      ::  vmax, vmin     !  maximum and minimum values of vscr

  real(REAL64)      ::  eref           !  reference energy for plot
  integer           ::  nocc           !  number of occupied states (different color) or recycled
  integer           ::  nstyle         !  choice of plot style

  integer           ::  irk,nrka
  integer           ::  iotape
  character(len=5)  ::  labelk
  integer           ::  ipr,nrk2,nd
  integer           ::  nsfft(3)

  real(REAL64)      ::  rk0(3)

  integer           ::  natot          !  total number of atoms
  integer           ::  neltot         !  total number of electrons
  real(REAL64)      ::  emidgap        !  rough estimate of the mid-gap

  real(REAL64)      ::  xsum
  integer           ::  mxdold         ! psiold mxddim

  integer           ::  nder

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64
  complex(REAL64), parameter  :: C_ZERO = cmplx(ZERO,ZERO,REAL64)

!  real(REAL64), parameter  :: XSC = 0.95_REAL64                     !  criteria for reduced vector size
  real(REAL64), parameter  :: XSC = 1.00_REAL64                     !  criteria for reduced vector size

! counters

  integer    ::  j, n, m


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

  iotape = 13
  call out_band_circuit_size('BAND_LINES.DAT', iotape, 1, adot, ztot,    &
                   neig, nrk2, nlines, nvert)

  allocate(xk(nrk2))
  allocate(rk(3,nrk2))
  allocate(xcvert(nvert))
  allocate(ljump(nlines))
  allocate(nkstep(nlines))
  allocate(label(nvert+nlines))
  allocate(xklab(nvert+nlines))

  call out_band_get_circuit('BAND_LINES.DAT', iotape, 1, adot,           &
                   xk, rk, xcvert, ljump, nkstep, label, xklab,          &
                   neig, nrk2, nlines, nvert)


  allocate(e_of_k(neig,nrk2))
  allocate(e_of_k_so(2*neig,nrk2))


! finds mxddim, mxdbnd

  mxdbnd = neig


  rk0(1) = 0.0
  rk0(2) = 0.0
  rk0(3) = 0.0

  call size_mtxd(emax, rk0, adot, ng, kgv, mtxd0)

  mxddim = mtxd0

  do irk=1,nrk2

!   loop over k-points

    do j=1,3
      rkpt(j) = rk(j,irk)
    enddo

    call size_mtxd(emax, rkpt, adot, ng, kgv, nd)

    if(nd > mxddim) mxddim = nd

  enddo

! allocates arrays


  allocate(ei(mxdbnd))
  allocate(ei_so(2*mxdbnd))
  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(psi(mxddim,mxdbnd))
  allocate(ekpsi(mxdbnd))

  allocate(isort_so(2*mxddim))
  allocate(psi_so(2*mxddim,2*mxdbnd))
  allocate(ekpsi_so(2*mxdbnd))

  allocate(imatch(mxdbnd))
  allocate(isold(mxddim))
  allocate(iperm(mxdbnd),ipermold(mxdbnd))

  allocate(imatch_so(2*mxdbnd))
  allocate(iperm_so(2*mxdbnd),ipermold_so(2*mxdbnd))

  allocate(ei0(mxddim))
  allocate(hdiag0(mxddim))
  allocate(isort0(mxddim))
  allocate(qmod0(mxddim))
  allocate(ekpg0(mxddim))
  allocate(psi0(mxddim,mxdbnd))
  allocate(hpsi0(mxddim,mxdbnd))
  allocate(ekpsi0(mxdbnd))

! Calculates reference k-point

  iguess = 0
  ipr = 0

  nocc = neig

  write(6,*)
  write(6,*) '   Calculating reference state for k.p'
  write(6,*)

  call h_kb_dia_all('pw  ', emax, rk0, neig, nocc,                       &
      flgpsd, ipr, ifail, icmax, iguess, epspsi,                         &
      ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                    &
      sfact, veff, icmplx,                                               &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      mtxd0, hdiag0, isort0, qmod0, ekpg0, .FALSE.,                      &
      psi0, hpsi0, ei0,                                                  &
      vscr, kmscr,                                                       &
      latorb, norbat, nqwf, delqwf, wvfao, lorb,                         &
      mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,            &
      mxdbnd, mxdscr, mxdlao)


  ipr = 1
  nrka = -1
  call print_eig(ipr, 0, labelk, nrka, rk0,                              &
      mtxd0, icmplx, neig, psi0,                                         &
      adot, ei0, ekpsi, isort0, kgv,                                     &
      mxddim, mxdbnd, mxdgve)

! prepares k.p matrices

  allocate(h0(mxdbnd,mxdbnd))
  allocate(dh0drk(mxdbnd,mxdbnd,3))
  allocate(d2h0drk2(mxdbnd,mxdbnd,3,3))

  nder = 2

  call kdotp_matrix(mtxd0, neig, psi0, ei0, rk0, isort0, nder,           &
      h0, dh0drk, d2h0drk2,                                              &
      ng, kgv,                                                           &
      ntype, natom, rat, adot,                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

  natot = 0
  do n = 1,ntype
    natot = natot + natom(n)
  enddo
  neltot = nint(ztot)
  j = neltot/2
  emidgap = real((h0(j,j) + h0(j+1,j+1)) / 2,REAL64)

! writes information for further processing

  call kdotp_silvaco_out('matrix_kdotp.dat', 23, neig, adot, rk0,        &
      h0,dh0drk, d2h0drk2, natot, neltot, emidgap, .FALSE.,              &
      mxdbnd)


  allocate(hso0(2*mxdbnd,2*mxdbnd))
  allocate(dhso0drk(2*mxdbnd,2*mxdbnd,3))
  allocate(d2hso0drk2(2*mxdbnd,2*mxdbnd,3,3))

  nder = 2

  call kdotp_matrix_so_pert(mtxd0, neig, psi0, ei0, rk0, isort0, nder,   &
      hso0, dhso0drk, d2hso0drk2,                                        &
      ng, kgv,                                                           &
      ntype, natom, rat, adot,                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

  call kdotp_matrix_so_convert(neig, hso0, dhso0drk, d2hso0drk2, nder,   &
           mxdbnd)

  j = neltot
  emidgap = real((hso0(j,j) + hso0(j+1,j+1)) / 2,REAL64)

  call kdotp_silvaco_out('matrix_kdotp_so.dat', 23, 2*neig, adot, rk0,   &
      hso0, dhso0drk, d2hso0drk2, natot, neltot, emidgap, .TRUE.,        &
      2*mxdbnd)


  allocate(psi_so0(2*mxddim,2*mxdbnd))

  do n = 1,neig
  do m = 1,mxddim
    psi_so0(2*m-1,2*n-1) = psi0(m,n)
    psi_so0(2*m  ,2*n-1) = C_ZERO
    psi_so0(2*m-1,2*n  ) = C_ZERO
    psi_so0(2*m  ,2*n  ) = psi0(m,n)
  enddo
  enddo

! allocate psiold here.  Defines mtxdold such that it
! has most of the spectral weight for highest eigenvalue

  xsum = ZERO
  do j = 1,mtxd0
    xsum = xsum + real(psi(j,neig)*conjg(psi(j,neig)),REAL64)
    mxdold = j
    if(xsum > XSC) exit
  enddo

  allocate(psiold(mxdold,mxdbnd))
  allocate(psiold_so(2*mxdold,2*mxdbnd))

  do irk=1,nrk2

!   loop over k-points

    do j=1,3
      rkpt(j) = rk(j,irk)
    enddo

    write(6,*)
    write(6,'("  Working on k-point ",i7,"  out of ", i7)') irk, nrk2
    write(6,*)


    call kdotp_diag(emax, rkpt, neig,                                    &
        rk0, mtxd0, isort0, psi0, h0, dh0drk, d2h0drk2,                  &
        psi, ei, mtxd, isort, qmod, ekpg,                                &
        ng, kgv, adot,                                                   &
        mxdgve, mxddim, mxdbnd)

    call kinetic_energy(neig, mtxd, ekpg, psi, ekpsi,                    &
    mxddim, mxdbnd)

    ipr = 1
    nrka = -1
    call print_eig(ipr, irk, labelk, nrka, rkpt,                         &
        mtxd, icmplx, neig, psi,                                         &
        adot, ei, ekpsi, isort, kgv,                                     &
        mxddim, mxdbnd, mxdgve)


    call kdotp_diag_so(emax, rkpt, neig,                                 &
        rk0, mtxd0, isort0, psi_so0, hso0, dhso0drk, d2hso0drk2,         &
        psi_so, ei_so, mtxd, isort, qmod, ekpg,                          &
        ng, kgv, adot,                                                   &
        mxdgve, mxddim, mxdbnd)

    ipr = 1
    nrka = -1

    if(ipr == 2) then
      call kinetic_energy_so(neig, mtxd, ekpg, psi_so, ekpsi_so,         &
      mxddim, mxdbnd)
    endif

    call print_eig_so(ipr, irk, labelk, nrka, rkpt,                      &
        mtxd, neig, psi_so,                                              &
        adot, ei_so, ekpsi_so, isort, kgv,                               &
        mxddim, mxdbnd, mxdgve)


    call out_band_match(irk, nlines, ljump, nkstep,                     &
         neig, icmplx, xsc,                                             &
         mtxd, psi, isort, psi_so, ei, ei_so,                           &
         e_of_k, e_of_k_so, nrk2,                                       &
         mxddim, mxdbnd)


  enddo

! writes the output files for xmgrace and gnuplot

  iotape = 15
  nstyle = 2

  call out_band_eref(neig, nrk2, rk, ztot, efermi, 2, 1, e_of_k, eref, nocc)

  call out_band_gnuplot('band_kp.gp', iotape,                            &
         neig, nrk2, xk, e_of_k, eref,                                   &
         nvert, xcvert, nlines, ljump, nkstep, label, xklab)

  call out_band_xmgrace('band_kp.agr', iotape,                           &
         title, subtitle, nstyle,                                        &
         neig, nrk2, xk, e_of_k, eref, nocc,                             &
         nvert, xcvert, nlines, ljump, nkstep, label, xklab)

  call out_band_eref(neig, nrk2, rk, ztot, efermi, 1, 1, e_of_k_so, eref, nocc)

  call out_band_gnuplot('band_so_kp.gp', iotape,                         &
      2*neig, nrk2, xk, e_of_k_so, eref,                                 &
      nvert, xcvert, nlines, ljump, nkstep, label, xklab)

  call out_band_xmgrace('band_so_kp.agr', iotape,                        &
     title, subtitle, nstyle,                                            &
     2*neig, nrk2, xk, e_of_k_so, eref, nocc,                            &
     nvert, xcvert, nlines, ljump, nkstep, label, xklab)


  deallocate(psiold)
  deallocate(imatch)
  deallocate(isold)
  deallocate(iperm,ipermold)

  deallocate(psiold_so)
  deallocate(imatch_so)
  deallocate(iperm_so,ipermold_so)

  deallocate(xk)
  deallocate(rk)
  deallocate(xcvert)
  deallocate(ljump)
  deallocate(nkstep)
  deallocate(label)
  deallocate(xklab)

  deallocate(e_of_k)

  deallocate(vscr)

  deallocate(ei)
  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)
  deallocate(psi)
  deallocate(ekpsi)

  deallocate(ei0)
  deallocate(hdiag0)
  deallocate(isort0)
  deallocate(qmod0)
  deallocate(ekpg0)
  deallocate(psi0)
  deallocate(ekpsi0)

  deallocate(psi_so)

  return

end subroutine out_band_kdotp_2nd
