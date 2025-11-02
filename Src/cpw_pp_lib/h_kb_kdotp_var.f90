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

!>  This subroutine calculates the Luttinger-Kohn
!>  variational k.p energies and eigenvectors.
!>
!>  \author       José Luís Martins
!>  \version      5.12
!>  \date         7 February 2014. 2 November 2025.
!>  \copyright    GNU Public License v2

subroutine h_kb_kdotp_var(emax, rkpt, neig, mtxd0, isort0, psi0,         &
    psi, ei, mtxd, isort, qmod, ekpg, lkpg,                              &
    ng, kgv,                                                             &
    vscr, kmscr, nqnl, delqnl, vkb, nkb,                                 &
    ntype, natom, rat, adot,                                             &
    mxdtyp, mxdatm, mxdgve, mxddim, mxdbnd, mxdlqp, mxdscr)


! Written 7 February 2014. JLM
! modified, documentation, lkpg, February 2020. JLM
! Modified, qmod-->ekpg in hk_psi. 13 February 2021. JLM
! Modified, nanlspin, 30 November 2023. JLM
! Modified, psi_convert. 2 November 2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  integer, intent(in)                ::  neig                            !<  number of eigenvectors (requested on input, modified by degeneracies on output)

  integer, intent(in)                ::  mtxd0                           !<  dimension of the reference hamiltonian
  integer, intent(in)                ::  isort0(mxddim)                  !<  g-vector associated with row/column i of reference hamiltonian
  complex(REAL64), intent(in)        ::  psi0(mxddim,mxdbnd)             !<  component j of reference eigenvector

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)           ::  vkb(mxdlqp,0:3,-1:1,mxdtyp)     !<  kb nonlocal pseudo. for atom k, ang. mom. l. (not normalized to vcell, Hartree)
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<   kb pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh

  logical, intent(in)                ::  lkpg                            !<  If true use the previous G-vectors (same mtxd and isort)

! output

  complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)              !<  component j of eigenvector i (guess on input)
  real(REAL64), intent(out)          ::  ei(mxdbnd)                      !<  eigenvalues (Hartree)
  integer, intent(out)               ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(out)               ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  real(REAL64), intent(out)          ::  qmod(mxddim)                    !<  length of k+g-vector of row/column i
  real(REAL64), intent(out)          ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

! allocatable arrays

  complex(REAL64), allocatable       ::  psi_lk(:,:)                     !  component j of eigenvector i
  complex(REAL64), allocatable       ::  hpsi_lk(:,:)                    !  component j of h*eigenvector i
  complex(REAL64), allocatable       ::  vec(:,:)                        !  overlap matrix S
  complex(REAL64), allocatable       ::  hred(:,:)                       !  reduced hamiltonian
  complex(REAL64), allocatable       ::  anlga(:,:)                      !  KB projectors without spin-orbit
  real(REAL64), allocatable          ::  xnlkb(:)                        !  KB normalization without spin-orbit

! local variables

  integer             ::  nanl, nanlso, nanlspin                         !  number of KB projectors
  integer             ::  mxdanl                                         !  array dimension of number of projectors

  logical             ::  lnewanl                                        !  indicates that anlga has been recalculated (not used in default implementation)

  integer             ::  info

! constants

  real(REAL64), parameter     :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  :: C_UM = cmplx(UM,ZERO,REAL64)
  complex(REAL64), parameter  :: C_ZERO = cmplx(ZERO,ZERO,REAL64)


  call hamilt_struct(emax, rkpt, mtxd, isort, qmod, ekpg, lkpg,          &
      ng, kgv, adot,                                                     &
      mxdgve, mxddim)

  allocate(psi_lk(mxddim,mxdbnd))
  allocate(hpsi_lk(mxddim,mxdbnd))

  call psi_convert(neig, 1, mtxd0, isort0, psi0, mtxd, isort, psi_lk,    &
      mxddim, mxdbnd)

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlspin,        &
       mxdtyp)

  mxdanl = nanl

  allocate(anlga(mxddim,mxdanl))
  allocate(xnlkb(mxdanl))

  call proj_nl_kb_c16(rkpt, mtxd, isort, nanl,                           &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      anlga, xnlkb,                                                      &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)

  lnewanl = .TRUE.

  call hk_psi_c16(mtxd, neig, psi_lk, hpsi_lk, lnewanl,                  &
      ng, kgv,                                                           &
      ekpg, isort, vscr, kmscr,                                          &
      anlga, xnlkb, nanl,                                                &
      mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

  allocate(vec(neig,neig))
  allocate(hred(neig,neig))

  call zgemm('c','n', neig, neig, mtxd, C_UM, psi_lk, mxddim,            &
      hpsi_lk, mxddim, C_ZERO, hred, neig)


  call diag_c16(neig, hred, ei, vec, neig, info)

  if(info /= 0) stop


  call zgemm('n','n', mtxd, neig, neig, C_UM, psi_lk, mxddim,            &
      vec, neig,C_ZERO, psi, mxddim)

  deallocate(vec)
  deallocate(hred)

  deallocate(psi_lk)
  deallocate(hpsi_lk)

  deallocate(anlga)
  deallocate(xnlkb)

  return
end subroutine h_kb_kdotp_var
