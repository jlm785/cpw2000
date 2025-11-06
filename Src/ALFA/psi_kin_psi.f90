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

!>  Calculates the matrix <Psi|K|Psi> for the kinetic energy
!>  operator for neig wavevectors.  complex*16 version
!>  Works for spin-polarized and no spin.
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         February 18, 2014. 4 November 2025.
!>  \copyright    GNU Public License v2

subroutine psi_kin_psi(mtxd, neig, psi, hkin, ekpg, nspin,               &
    mxddim, mxdbnd)

! Written February 18, 2014, from hk_psi_c16.   jlm
! See that file for historical record.
! Modified, documentation, 21 January 2020. JLM
! Modified, qmod-->ekpg in hk_psi. 13 February 2021. JLM
! Added spin polarization option, 4 November 2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves (not counting spin)
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands (including spin)

  integer, intent(in)                ::  nspin                           !<  spin components (1:no spin or 2:spin present)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (not counting spin)
  integer, intent(in)                ::  neig                            !<  wavefunction dimension (including spin)
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

  complex(REAL64), intent(in)        ::  psi(nspin*mxddim,mxdbnd)        !<  wavevectors

! output

  complex(REAL64), intent(out)       ::  hkin(mxdbnd,nspin)              !<  <Psi|K|Psi>

! local variables

  complex(REAL64), allocatable     ::  hpsi(:,:)

! constants

  real(REAL64), parameter     :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)


  allocate(hpsi(nspin*mxddim,neig))

  call hk_psi_kin_c16(mtxd, neig, psi, hpsi, ekpg, .FALSE., nspin,       &
      mxddim, mxdbnd)

  call zgemm('c','n', neig, neig, nspin*mtxd, C_UM, psi, nspin*mxddim,   &
      hpsi, nspin*mxddim, C_ZERO, hkin, mxdbnd)

  deallocate(hpsi)

  return

end subroutine psi_kin_psi
