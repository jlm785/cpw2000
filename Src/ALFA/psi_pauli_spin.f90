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

!>  Calculates the expectation value of the Pauli matrices.
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         28 January 2025.
!>  \copyright    GNU Public License v2

subroutine psi_pauli_psi(mtxd, neig, psi_sp, pauli,                      &
    mxddim, mxdbnd)

! Written 28 January 2025 based on psi_orient_spin. JLM



  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  number of states


  complex(REAL64), intent(in)        ::  psi_sp(2*mxddim, 2*mxdbnd)      !<  |psi> spin-wavefunctions

! output

  complex(REAL64), intent(out)       ::  pauli(3, 2*mxdbnd)             !<  <psi|sigma_xyz|psi>

! local allocatable variables

  complex(REAL64), allocatable       ::  psi_tmp(:)                     !  temporary |psi>

! parameters

  real(REAL64), parameter       ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter    ::  C_UM = cmplx(UM,ZERO,REAL64)
  complex(REAL64), parameter    ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter    ::  C_I = cmplx(ZERO,UM,REAL64)

! counters

  integer    ::  i, n

! external functions

  complex(REAL64),external      :: zdotc


  allocate(psi_tmp(2*mtxd))

! loop over states

  do n = 1, neig

! sigma_x

    do i = 1,mtxd
      psi_tmp(2*i-1) =  psi_sp(2*i  ,n)
      psi_tmp(2*i  ) =  psi_sp(2*i-1,n)
    enddo

    pauli(1,n) = zdotc(2*mtxd, psi_sp(:,n), 1, psi_tmp, 1)

! sigma_y

    do i = 1,mtxd
      psi_tmp(2*i-1) = -C_I*psi_sp(2*i  ,n)
      psi_tmp(2*i  ) =  C_I*psi_sp(2*i-1,n)
    enddo

    pauli(2,n) = zdotc(2*mtxd, psi_sp(:,n), 1, psi_tmp, 1)

! sigma_z

    do i = 1,mtxd
      psi_tmp(2*i-1) =  psi_sp(2*i-1,n)
      psi_tmp(2*i  ) = -psi_sp(2*i  ,n)
    enddo

    pauli(3,n) = zdotc(2*mtxd, psi_sp(:,n), 1, psi_tmp, 1)

  enddo

  deallocate(psi_tmp)

  return

end subroutine psi_pauli_psi
