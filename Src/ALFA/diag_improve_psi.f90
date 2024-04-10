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

subroutine diag_improve_psi(rkpt, mtxd, neig, njac, nritz, tol,          &
    ei, psi, hpsi,                                                       &
    ng, kgv,                                                             &
    ekpg, isort, hdiag, vscr, kmscr,                                     &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve, mxdscr)


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
  integer, intent(in)                ::  mxdscr                          !< array dimension of vscr

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  number of eigenvectors (without spin)

  integer, intent(in)                ::  njac                            !< number of jacobi iterations
  integer, intent(in)                ::  nritz                           !< number of ritz iterations
!  real(REAL64), intent(in)           ::  epsdeg                          !< two eigenvalues are considered quasi degenerate if |e_i - e_j| < eps
  real(REAL64), intent(in)           ::  tol                             !<  criteria for convergence

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates

  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  length of k+g-vector of row/column i
  real(REAL64), intent(in)           ::  hdiag(mxddim)                   !<  hamiltonian diagonal
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

! output

  real(REAL64), intent(out)          ::  ei(mxdbnd)                      !<  spin-orbit eigenvalues
  complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)             !<  H | psi >

! input and output

  complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)              !<  eigenvectors | psi >

! allocatable local arrays

  complex(REAL64), allocatable       ::  anlga(:,:)                      !  KB projectors
  real(REAL64), allocatable          ::  xnlkb(:)                        !  KB normalization

! main local variables

  integer       ::  mxdanl                                               !  array dimension of number of projectors
  integer       ::  nanl                                                 !  half of number of projectors without spin
  integer       ::  nanlso, nanlsp                                       !  number of projectors witn spin

  integer       ::  nout

  logical       ::  lnewanl                                              ! indicates that anlga has been recalculated (not used in default implementation)


  lnewanl = .TRUE.

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlsp,        &
    mxdtyp)

  mxdanl = nanl

  allocate(anlga(mxddim,mxdanl))
  allocate(xnlkb(mxdanl))

  call proj_nl_kb_c16(rkpt, mtxd, isort, nanl,                         &
      ng, kgv,                                                         &
      nqnl, delqnl, vkb, nkb,                                          &
      ntype, natom, rat, adot,                                         &
      anlga, xnlkb,                                                    &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)

  call diag_jacobi_ritz_c16(mtxd, neig, ei,                            &
      psi, hpsi, nout,                                                 &
      njac, nritz, tol,                                                &
      ng, kgv,                                                         &
      ekpg, isort, hdiag, vscr, kmscr,                                 &
      anlga, xnlkb, nanl, lnewanl,                                     &
      mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

  deallocate(anlga)
  deallocate(xnlkb)

  return

end subroutine diag_improve_psi
