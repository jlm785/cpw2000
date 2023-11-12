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

!>   Prepares the reference energies and wavefunctions
!>   For the GLK interpolation
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.09
!>  \date         23 September 2020, 12 November 2023.
!>  \copyright    GNU Public License v2

subroutine out_glk_prepare(diag_type, io66,                              &
      nrk3, rk_ref,                                                      &
      emax, neig, flgpsd,                                                &
      epspsi, icmax,                                                     &
      adot, ntype, natom, rat,                                           &
      ng, kgv, phase, conj,                                              &
      ns, inds, kmax, indv, ek,                                          &
      sfact, icmplx,                                                     &
      veff,                                                              &
      nqnl, delqnl, vkb, nkb,                                            &
      vscr, kmscr,                                                       &
      latorb, norbat, nqwf, delqwf, wvfao, lorb,                         &
      mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxddim,            &
      mxdbnd, mxdscr, mxdlao)


! Extracted from out_dos_glk, 23 September 2020. JLM
! Name change, indentation. 12 November 2023.

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdscr                          !<  array dimension for vscr

  character(len=4), intent(in)       ::  diag_type                       !<  selects diagonalization, 'pw  ','ao  ','aojc'

  integer, intent(in)                ::  io66                            !<  tape number

  integer, intent(in)                ::  nrk3                            !<  number of reference k-points
  real(REAL64), intent(in)           ::  rk_ref(3,nrk3)                  !<  reference k-points

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential

! input and output

  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors
  integer, intent(in)                ::  icmax                           !<  maximum number of iterations for diagonalization

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
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<   KB pseudo.  normalization for atom k, ang. mom. l

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh and fft mesh size
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh

  logical, intent(in)                ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  wavefunction for atom k, ang. mom. l (normalized to vcell)
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k

! input and output

  integer, intent(inout)             ::  neig                            !<  number of eigenvectors (requested on input, modified by


! allocatable arrays for eigensolutions

  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi(:,:)                        !  component j of eigenvector i (guess on input)
  complex(REAL64), allocatable       ::  hpsi(:,:)                       !  H | psi>
  real(REAL64), allocatable          ::  ekpsi(:)                        !  kinetic energy of eigenvector i. (hartree)

! local variables

  integer                            ::  mtxd                            !  dimension of the hamiltonian

  integer                            ::  iguess                          !   if guess eigenvectors are available, iguess = 1, otherwise iguess = 0

  integer                            ::  ipr
  integer                            ::  nrka
  character(len=5)                   ::  labelk

  real(REAL64)                       ::  rkpt(3)                         !  k-vector

  logical                            ::  lkpg                            !  If true use the previous G-vectors (same mtxd and isort)
  integer                            ::  nocc
  integer                            ::  ifail                           !  if ifail=0 ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.

! parameters

  real(REAL64), parameter            ::  UM = 1.0_REAL64

! counters

  integer                            ::  irk



! allocates arrays

  allocate(ei(mxdbnd))
  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(psi(mxddim,mxdbnd))
  allocate(hpsi(mxddim,mxdbnd))
  allocate(ekpsi(mxdbnd))

  iguess = 0

  lkpg = .FALSE.
  ipr = 0

  nocc = neig

  do irk = 1,nrk3

    rkpt(1) = rk_ref(1,irk)
    rkpt(2) = rk_ref(2,irk)
    rkpt(3) = rk_ref(3,irk)

    call h_kb_dia_all(diag_type, emax, rkpt, neig, nocc,                 &
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


    call kinetic_energy(neig, mtxd, ekpg, psi, ekpsi,                    &
        mxddim, mxdbnd)


    ipr = 1
    nrka = -1
    call print_eig(ipr, irk, labelk, nrka, rkpt,                         &
        mtxd, icmplx, neig, psi,                                         &
        adot, ei, ekpsi, isort, kgv,                                     &
        mxddim, mxdbnd, mxdgve)

    write(io66,rec=irk)  irk, mtxd, rkpt(1:3), psi(:,:), isort(:)

  enddo

  deallocate(ei)
  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)
  deallocate(psi)
  deallocate(hpsi)
  deallocate(ekpsi)

  return

end subroutine out_glk_prepare
