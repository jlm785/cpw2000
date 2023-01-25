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

!>  Aligns a set of wave-vectors at one k-point with another
!>  set at another k-point, band by band.
!>  Applies the same phase to another vector.
!>  Does not align if the wave-vectors are not similar.
!>  Check berry_psi_align_svd for aligning subspaces.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         24 January 2023.
!>  \copyright    GNU Public License v2

subroutine berry_psi_align_two(neig ,mtxd, psi_ref, psi, hpsi,           &
     mxddim, mxdbnd)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddim                          !<  array dimension for the hamiltonian
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands

  integer, intent(in)                ::  neig                            !<  number of eigenvectors (requested on input, modified by degeneracies on output)
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  complex(REAL64), intent(in)        ::  psi_ref(mxddim,mxdbnd)          !<  component j of eigenvector i (reference)

! input and output

  complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)              !<  component j of eigenvector i (only change of global phase)
  complex(REAL64), intent(inout)     ::  hpsi(mxddim,mxdbnd)             !<  auxiliary vector (only change of global phase)

! local varaibles

  complex(REAL64)         ::  phase
  real(REAL64)            ::  xabs

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter ::  C_UM = cmplx(UM,ZERO,REAL64)


  complex(REAL64), external   ::  zdotc

! counters

  integer ::  n

  do n = 1,neig
    phase = zdotc(mtxd, psi_ref(:,n),1,psi(:,n),1)
    xabs = abs(phase)
    if(xabs > 0.7) then
      phase = conjg(phase) / xabs
      psi(:,n) = phase*psi(:,n)
      hpsi(:,n) = phase*hpsi(:,n)
    endif
  enddo

  return

end subroutine berry_psi_align_two
