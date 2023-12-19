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

!>  Iterative diagonalization of the spin-hamiltonian.
!>  starting from the non-spin wave-functions.
!>  First it does 1st order perturbation theory,
!>  then does jacobi-ritz iteration
!>
!>  \author       José Luís Martins
!>  \version      5.09
!>  \date         18 December 2023.
!>  \copyright    GNU Public License v2

subroutine diag_improve_psi_spin(rkpt, mtxd, neig, njac, nritz, tol,     &
    ei, psi,                                                             &
    ei_pert, ei_sp, psi_sp, hpsi_sp,                                     &
    ng, kgv,                                                             &
    ekpg, isort, vscr_sp, kmscr, nsp,                                    &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve, mxdscr, mxdnsp)


! Written 18 December 2023 from early out_mass_berry code.  JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors
  integer, intent(in)                ::  mxdscr                          !< array dimension of vscr_sp
  integer, intent(in)                ::  mxdnsp                          !<  array dimension for number of spin components (1,2,4)

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  number of eigenvectors (without spin)

  integer, intent(in)                ::  njac                            !< number of jacobi iterations
  integer, intent(in)                ::  nritz                           !< number of ritz iterations
!  real(REAL64), intent(in)           ::  epsdeg                          !< two eigenvalues are considered quasi degenerate if |e_i - e_j| < eps
  real(REAL64), intent(in)           ::  tol                             !<  criteria for convergence

  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalues without spin
  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  eigenvectors without spin

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates

  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  length of k+g-vector of row/column i
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  real(REAL64), intent(in)           ::  vscr_sp(mxdscr,mxdnsp)          !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh
  integer, intent(in)                ::  nsp                             !<  number of spin components ox xc-potential (1,2,4)

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

! output

  real(REAL64), intent(out)          ::  ei_pert(2*mxdbnd)               !<  spin-orbit eigenvalues by perturbation theory
  real(REAL64), intent(out)          ::  ei_sp(2*mxdbnd)                 !<  spin-orbit eigenvalues
  complex(REAL64), intent(out)       ::  psi_sp(2*mxddim,2*mxdbnd)       !<  vectors in the spin-orbit form | psi_sp >
  complex(REAL64), intent(out)       ::  hpsi_sp(2*mxddim,2*mxdbnd)      !<  H | psi_sp >

! allocatable local arrays

  complex(REAL64), allocatable       ::  anlsp(:,:)
  real(REAL64), allocatable          ::  xnlkbsp(:)

! main local variables

  integer       ::  mxdanl                                               !  array dimension of number of projectors
  integer       ::  mxdasp                                               !  array dimension of number of projectors with spin-orbit
  integer       ::  nanl                                                 !  half of number of projectors without spin
  integer       ::  nanlso, nanlsp                                       !  number of projectors witn spin

  integer       ::  nout

  logical       ::  lnewanl                                              ! indicates that anlsp has been recalculated (not used in default implementation)


  lnewanl = .TRUE.

  call spin_orbit_perturb(rkpt, mtxd, isort,                           &
      neig, psi, ei, ei_pert, psi_sp, lnewanl,                         &
      ng, kgv,                                                         &
      nqnl, delqnl, vkb, nkb,                                          &
      ntype, natom, rat, adot,                                         &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlsp,        &
    mxdtyp)

  mxdanl = nanl
  mxdasp = nanlsp

  allocate(anlsp(2*mxddim,mxdasp))
  allocate(xnlkbsp(mxdasp))

  call proj_nl_kb_so_c16(rkpt, mtxd, isort,                            &
      nanlsp,                                                          &
      ng, kgv,                                                         &
      nqnl, delqnl, vkb, nkb,                                          &
      ntype, natom, rat, adot,                                         &
      anlsp, xnlkbsp,                                                  &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdasp, mxdgve)

  call diag_jacobi_ritz_spin_c16(mtxd, 2*neig, ei_sp,                  &
      psi_sp, hpsi_sp, nout,                                           &
      njac, nritz, tol,                                                &
      ng, kgv,                                                         &
      ekpg, isort, vscr_sp, kmscr, nsp,                                &
      anlsp, xnlkbsp, nanlsp, lnewanl,                                 &
      mxddim, 2*mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

  deallocate(anlsp)
  deallocate(xnlkbsp)

  return

end subroutine diag_improve_psi_spin
