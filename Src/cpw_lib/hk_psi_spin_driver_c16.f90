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

!>  Calculates the product of the hamiltonian times neig spin-wavevectors,
!>  either the Kohn-Sham hamiltonian (hk_psi_spin_c16) or the generalized
!>  Kohn-Sham hamiltonian of a meta-GGA (hk_psi_spin_mgga_c16).
!>  vtau does not depend on spin (vtau_sp has only the first component).
!>
!>  \author       claude
!>  \version      5.13
!>  \date         8 October 2026.
!>  \copyright    GNU Public License v2

subroutine hk_psi_spin_driver_c16(lgks, mtxd, neig, psi_sp, hpsi_sp,     &
    lnewanl,                                                             &
    ng, kgv, rkpt, adot,                                                 &
    ekpg, isort, vscr_sp, vtaumsh, kmscr, nsp,                           &
    anlsp, xnlkbsp, nanlsp,                                              &
    mxddim, mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

! Written 8 October 2026 from hk_psi_driver_c16. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdasp                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr_sp
  integer, intent(in)                ::  mxdnsp                          !<  array dimension for number of spin components (1,2,4)

  logical, intent(in)                ::  lgks                            !<  generalized Kohn-Sham meta-GGA (vtau term)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size, does not include spin)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  real(REAL64), intent(in)           ::  vscr_sp(mxdscr,mxdnsp)          !<  screened potential in the fft real space mesh (1, sigma_z, sigma_x, sigma_y)
  real(REAL64), intent(in)           ::  vtaumsh(mxdscr)                 !<  d (rho eps_xc) / d tau in the fft real space mesh (only used if lgks)
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size
  integer, intent(in)                ::  nsp                             !<  number of spin components of the potential (1,2,4)

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nanlsp                          !<  number of projectors
  complex(REAL64), intent(in)        ::  anlsp(2*mxddim,mxdasp)          !<  Kleinman-Bylander projectors
  real(REAL64), intent(in)           ::  xnlkbsp(mxdasp)                 !<  Kleinman-Bylander normalization
  complex(REAL64), intent(in)        ::  psi_sp(2*mxddim,mxdbnd)         !<  spin-wavevectors

! input and output

  logical, intent(inout)             ::  lnewanl                         !<  indicates that anlsp has been recalculated

! output

  complex(REAL64), intent(out)       ::  hpsi_sp(2*mxddim,mxdbnd)        !<  |hpsi_sp> =  H |psi_sp>

! local allocatable arrays

  real(REAL64), allocatable          ::  vtau_sp(:,:)                    !  vtau with spin components, only the first is not zero

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64


  if(lgks) then

    allocate(vtau_sp(mxdscr,mxdnsp))
    vtau_sp(:,:) = ZERO
    vtau_sp(:,1) = vtaumsh(:)

    call hk_psi_spin_mgga_c16(mtxd, neig, psi_sp, hpsi_sp, lnewanl,      &
        ng, kgv, rkpt, adot,                                             &
        ekpg, isort, vscr_sp, vtau_sp, kmscr, nsp,                       &
        anlsp, xnlkbsp, nanlsp,                                          &
        mxddim, mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

    deallocate(vtau_sp)

  else

    call hk_psi_spin_c16(mtxd, neig, psi_sp, hpsi_sp, lnewanl,           &
        ng, kgv,                                                         &
        ekpg, isort, vscr_sp, kmscr, nsp,                                &
        anlsp, xnlkbsp, nanlsp,                                          &
        mxddim, mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

  endif

  return

end subroutine hk_psi_spin_driver_c16
