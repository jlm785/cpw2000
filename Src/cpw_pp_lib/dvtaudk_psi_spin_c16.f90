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

!>  Spinor version of dvtaudk_psi_c16.  As vtau does not depend on
!>  spin the two components (interleaved, 2*i-1 and 2*i) are treated
!>  as independent wave-functions.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         8 October 2026.
!>  \copyright    GNU Public License v2

subroutine dvtaudk_psi_spin_c16(mtxd, neig, psi_sp, vtaupsi_sp,         &
    dvtaupsi_sp,                                                         &
    rkpt, adot, isort, kgv, vtaumsh, kmscr,                              &
    mxddim, mxdbnd, mxdgve, mxdscr)

! Written 8 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vtaumsh

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size, does not include spin)
  integer, intent(in)                ::  neig                            !<  number of spin-wavefunctions
  complex(REAL64), intent(in)        ::  psi_sp(2*mxddim,mxdbnd)         !<  spin-wavevectors

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  real(REAL64), intent(in)           ::  vtaumsh(mxdscr)                 !<  d (rho eps_xc) / d tau in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

! output

  complex(REAL64), intent(out)       ::  vtaupsi_sp(2*mxddim,mxdbnd)     !<  vtau |psi_sp>
  complex(REAL64), intent(out)       ::  dvtaupsi_sp(2*mxddim,mxdbnd,3)  !<  (d H_tau / d k_l) |psi_sp>, k in lattice coordinates

! local allocatable arrays

  complex(REAL64), allocatable       ::  psi(:,:)                        !  spin components as independent wave-functions
  complex(REAL64), allocatable       ::  vtaupsi(:,:)
  complex(REAL64), allocatable       ::  dvtaupsi(:,:,:)

! counters

  integer   ::   m, n, l


  allocate(psi(mtxd,2*neig))
  allocate(vtaupsi(mtxd,2*neig))
  allocate(dvtaupsi(mtxd,2*neig,3))

  do n = 1,neig
    do m = 1,mtxd
      psi(m,n     ) = psi_sp(2*m-1,n)
      psi(m,n+neig) = psi_sp(2*m  ,n)
    enddo
  enddo

  call dvtaudk_psi_c16(mtxd, 2*neig, psi, vtaupsi, dvtaupsi,             &
      rkpt, adot, isort, kgv, vtaumsh, kmscr,                            &
      mtxd, 2*neig, mxdgve, mxdscr)

  do n = 1,neig
    do m = 1,mtxd
      vtaupsi_sp(2*m-1,n) = vtaupsi(m,n     )
      vtaupsi_sp(2*m  ,n) = vtaupsi(m,n+neig)
    enddo
  enddo
  do l = 1,3
    do n = 1,neig
      do m = 1,mtxd
        dvtaupsi_sp(2*m-1,n,l) = dvtaupsi(m,n     ,l)
        dvtaupsi_sp(2*m  ,n,l) = dvtaupsi(m,n+neig,l)
      enddo
    enddo
  enddo

  deallocate(psi)
  deallocate(vtaupsi)
  deallocate(dvtaupsi)

  return

end subroutine dvtaudk_psi_spin_c16
