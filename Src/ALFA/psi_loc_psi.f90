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

!>  Calculates the matrix <Psi|V_loc|Psi> for a local
!>  pseudopotential V_loc for neig wavevectors.  complex*16 version
!>  Works for spin-polarized and no spin.
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         February 18, 2014. 4 November 2025.
!>  \copyright    GNU Public License v2

subroutine psi_loc_psi(mtxd, neig, psi, hloc, nspin, nsp,                &
    ng, kgv, isort,                                                      &
    vscr, kmscr,                                                         &
    mxddim, mxdbnd, mxdgve, mxdscr, mxdnsp)

! Written February 18, 2014, from hk_psi_c16.   jlm
! See that file for historical record.
! Modified 8 November 2015. Compatibility new libpw. JLM
! Modified, documentation, 29 February 2020. JLM
! Added spin polarization option, 4 November 2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves (not counting spin)
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands (not counting spin)
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr
  integer, intent(in)                ::  mxdnsp                          !<  array dimension for number of spin components (1,2,4)

  integer, intent(in)                ::  nspin                           !<  spin components (1:no spin or 2:spin present)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  real(REAL64), intent(in)           ::  vscr(mxdscr,mxdnsp)             !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh
  integer, intent(in)                ::  nsp                             !<  number of spin components ox xc-potential (1,2,4)

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  complex(REAL64), intent(in)        ::  psi(nspin*mxddim,nspin*mxdbnd)  !<  wavevector

! output

  complex(REAL64), intent(out)       ::  hloc(nspin*mxdbnd,nspin*mxdbnd)  !<  <Psi|V_loc|Psi>

! local variables

  complex(REAL64), allocatable     ::  hpsi(:,:)

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)


  allocate(hpsi(nspin*mxddim,neig))

  call hk_psi_loc_c16(mtxd, neig, psi, hpsi, .FALSE., nspin, nsp,        &
      ng, kgv,                                                           &
      isort, vscr, kmscr,                                                &
      mxddim, mxdbnd, mxdgve, mxdscr, mxdnsp)

  call zgemm('c', 'n', neig, neig, nspin*mtxd, C_UM, psi, nspin*mxddim,  &
             hpsi, nspin*mxddim, C_ZERO, hloc, nspin*mxdbnd)

  deallocate(hpsi)

  return

end subroutine psi_loc_psi
