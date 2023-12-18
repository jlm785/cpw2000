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

!>  Calculates the derivatives of V_NL|psi>  + K |psi> for a separable non-local
!>  pseudopotential with respect to the k-vector.
!>
!>  spin-wavefunction version.
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      5.09
!>  \date         15 December 2023.
!>  \copyright    GNU Public License v2


subroutine berry_dhdk_psi_spin(rkpt, adot, mtxd, neig,                   &
    psi_sp, dhdkpsi_sp,                                                  &
    kgv, isort,                                                          &
    nanlsp, anlsp, xnlkbsp, danlspdrk,                                   &
    mxddim, mxdbnd, mxdgve, mxdasp)

! adapted from non-spin version, 15 December 2023. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors
  integer, intent(in)                ::  mxdasp                          !<  array dimension of number of projectors

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  neig                            !<  wavefunction dimension

  complex(REAL64), intent(in)        ::  psi_sp(2*mxddim, mxdbnd)        !<  |psi_sp> (in principle eigen-functions)

  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates
  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction

  integer, intent(in)                ::  nanlsp                          !<  half of number of projectors without spin
  complex(REAL64), intent(in)        ::  anlsp(2*mxddim,mxdasp)          !<  KB projectors without spin-orbit
  real(REAL64), intent(in)           ::  xnlkbsp(mxdasp)                 !<  KB normalization without spin-orbit
  complex(REAL64), intent(in)        ::  danlspdrk(2*mxddim,mxdasp,3)    !<  d anlsp / d rkpt

! output

  complex(REAL64), intent(out)       ::  dhdkpsi_sp(2*mxddim,mxdbnd,3)   !<  (d H /d k) |psi_sp>

! local allocatable variables

  real(REAL64), allocatable          ::  qcontra(:,:)                    !  contravariant rkpt+kgv
  real(REAL64), allocatable          ::  bdot(:,:)                       !  metric in reciprocal space.  It is allocatable to keep subroutine thread-safe

! local variables

  real(REAL64)      ::  qk1, qk2, qk3

  real(REAL64)      ::  vcell

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer    ::  j, m, n

! calculates the non local part

  call berry_dvnldk_psi(2*mtxd, neig, psi_sp, dhdkpsi_sp,                &
    nanlsp, anlsp, xnlkbsp, danlspdrk,                                   &
    2*mxddim, mxdbnd, mxdasp)

! kinetic energy operator derivative.  momentum operator in contravariant coordinates

  allocate(qcontra(3,mtxd))
  allocate(bdot(3,3))

  call adot_to_bdot(adot,vcell,bdot)

  do m = 1,mtxd
    qk1 = rkpt(1) + UM*(kgv(1,isort(m)))
    qk2 = rkpt(2) + UM*(kgv(2,isort(m)))
    qk3 = rkpt(3) + UM*(kgv(3,isort(m)))
    do j = 1,3
      qcontra(j,m) = bdot(j,1)*qk1 + bdot(j,2)*qk2 + bdot(j,3)*qk3
    enddo
  enddo

  do j = 1,3
    do n = 1,neig
      do m = 1,mtxd
        dhdkpsi_sp(2*m-1,n,j) = dhdkpsi_sp(2*m-1,n,j) + qcontra(j,m)*psi_sp(2*m-1,n)
        dhdkpsi_sp(2*m  ,n,j) = dhdkpsi_sp(2*m  ,n,j) + qcontra(j,m)*psi_sp(2*m  ,n)
      enddo
    enddo
  enddo

  deallocate(qcontra)
  deallocate(bdot)

  return

end subroutine berry_dhdk_psi_spin

