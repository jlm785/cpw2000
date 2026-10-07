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

!>  Calculates the product of the hamiltonian times neig wavevectors,
!>  either the Kohn-Sham hamiltonian (hk_psi_c16) or the generalized
!>  Kohn-Sham hamiltonian of a meta-GGA (hk_psi_mgga_c16).
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         6 October 2026.
!>  \copyright    GNU Public License v2

subroutine hk_psi_driver_c16(lgks, mtxd, neig, psi, hpsi, lnewanl,       &
    ng, kgv, rkpt, adot,                                                 &
    ekpg, isort, vscr, vtau, kmscr,                                      &
    anlga, xnlkb, nanl,                                                  &
    mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

! Written 6 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr

  logical, intent(in)                ::  lgks                            !<  generalized Kohn-Sham meta-GGA (vtau term)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  real(REAL64), intent(in)           ::  vtau(mxdscr)                    !<  d (rho eps_xc) / d tau in the fft real space mesh (only used if lgks)
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nanl                            !<  number of projectors
  complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)            !<  Kleinman-Bylander projectors
  real(REAL64), intent(in)           ::  xnlkb(mxdanl)                   !<  Kleinman-Bylander normalization
  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  wavevector

! input and output

  logical, intent(inout)             ::  lnewanl                         !<  indicates that anlga has been recalculated

! output

  complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)             !<  |hpsi> =  H |psi>


  if(lgks) then

    call hk_psi_mgga_c16(mtxd, neig, psi, hpsi, lnewanl,                 &
        ng, kgv, rkpt, adot,                                             &
        ekpg, isort, vscr, vtau, kmscr,                                  &
        anlga, xnlkb, nanl,                                              &
        mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

  else

    call hk_psi_c16(mtxd, neig, psi, hpsi, lnewanl,                      &
        ng, kgv,                                                         &
        ekpg, isort, vscr, kmscr,                                        &
        anlga, xnlkb, nanl,                                              &
        mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

  endif

  return

end subroutine hk_psi_driver_c16
