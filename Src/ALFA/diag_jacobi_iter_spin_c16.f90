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

!>  Performs a succession of Jacobi relaxations
!>  For the spinor representation.
!>
!>  \author       Jose Luis Martins
!>  \version      5.09
!>  \date         13 December 2023.
!>  \copyright    GNU Public License v2

subroutine diag_jacobi_iter_spin_c16(mtxd, neig, psi_sp, hpsi_sp, njac,  &
      dpsi_sp, hdpsi_sp,                                                 &
      ng, kgv,                                                           &
      ekpg, isort, vscr_sp, kmscr, nsp,                                  &
      anlsp, xnlkbsp, nanlsp, lnewanl,                                   &
      mxddim, mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

! Adapted from the non spin version. 13 December 2023. JLM


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdasp                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr_sp
  integer, intent(in)                ::  mxdnsp                          !<  array dimension for number of spin components (1,2,4)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions

  integer, intent(in)                ::  njac                            !<  number of jacobi iterations

  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  length of k+g-vector of row/column i
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  real(REAL64), intent(in)           ::  vscr_sp(mxdscr,mxdnsp)          !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh
  integer, intent(in)                ::  nsp                             !<  number of spin components ox xc-potential (1,2,4)

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nanlsp                            !<  number of projectors
  complex(REAL64), intent(in)        ::  anlsp(2*mxddim,mxdasp)            !<  Kleinman-Bylander projectors
  real(REAL64), intent(in)           ::  xnlkbsp(mxdasp)                   !<  Kleinman-Bylander normalization

  complex(REAL64), intent(in)        ::  psi_sp(2*mxddim,mxdbnd)              !<  |psi> to be improved

! output

  complex(REAL64), intent(out)       ::  hpsi_sp(2*mxddim,mxdbnd)             !<  H |psi_sp>
  complex(REAL64), intent(out)       ::  dpsi_sp(2*mxddim,mxdbnd)             !<  correction to |psi_sp>.  Used has work vector in subroutine
  complex(REAL64), intent(out)       ::  hdpsi_sp(2*mxddim,mxdbnd)            !<  H |dpsi_sp>

! input and output

  logical, intent(inout)             ::  lnewanl                         !<  indicates that anlsp has been recalculated (not used in default implementation)

! local variables

  real(REAL64)            ::  xn
!  complex(REAL64)         ::  cn

! local allocatable arrays

  complex(REAL64), allocatable       ::  xerror(:,:)                     !  error vector  H |psi_sp) - |psi_sp><psi_sp| H |psi_sp> or respective jacobian relaxed vector
  complex(REAL64), allocatable       ::  hxerror(:,:)                    !  H | relaxed xerror >
  real(REAL64), allocatable          ::  eg(:)                           !  <psi_sp| H |psi_sp>

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! external functions

  complex(REAL64),external   :: zdotc

! counters

  integer   ::   n, j


  allocate(eg(mxdbnd))
  allocate(xerror(2*mxddim,mxdbnd))
  allocate(hxerror(2*mxddim,mxdbnd))

! apply H to initial eigenvectors

  call hk_psi_spin_c16(mtxd, neig, psi_sp, hpsi_sp, lnewanl,             &
      ng, kgv,                                                           &
      ekpg, isort, vscr_sp, kmscr, nsp,                                  &
      anlsp, xnlkbsp, nanlsp,                                            &
      mxddim, mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

  do n = 1,neig
    call zcopy(2*mtxd, psi_sp(1,n), 1, dpsi_sp(1,n), 1)
  enddo

  do n = 1,neig
    call zcopy(2*mtxd, hpsi_sp(1,n), 1, hdpsi_sp(1,n), 1)
  enddo

! start relaxation

  do j = 1,njac

!   xerror is the error vector or jacobian relaxed error

    do n = 1,neig
      eg(n) = real(zdotc(2*mtxd, dpsi_sp(1, n), 1, hdpsi_sp(1, n), 1), REAL64)
      call zcopy(2*mtxd, hdpsi_sp(1, n), 1, xerror(1, n), 1)
      call zaxpy(2*mtxd, cmplx(-eg(n), ZERO, REAL64),                    &
                           dpsi_sp(1, n), 1, xerror(1, n), 1)
    enddo

    call diag_rq_jac_spin_c16(xerror, eg, ekpg, mtxd, neig,              &
    mxddim, mxdbnd)


    call hk_psi_spin_c16(mtxd, neig, xerror, hxerror, lnewanl,           &
      ng, kgv,                                                           &
      ekpg, isort, vscr_sp, kmscr, nsp,                                  &
      anlsp, xnlkbsp, nanlsp,                                            &
      mxddim, mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

    do n = 1,neig

!     adds

      call zaxpy(2*mtxd, C_UM, xerror(1, n), 1, dpsi_sp(1, n), 1)
      call zaxpy(2*mtxd, C_UM, hxerror(1, n), 1, hdpsi_sp(1, n), 1)
      xn = real(zdotc(2*mtxd, dpsi_sp(1, n), 1, dpsi_sp(1, n), 1), REAL64)
      xn = UM/sqrt(xn)
      call zscal(2*mtxd, cmplx(xn, ZERO, REAL64), dpsi_sp(1, n), 1)
      call zscal(2*mtxd, cmplx(xn, ZERO, REAL64), hdpsi_sp(1, n), 1)

    enddo

  enddo

  deallocate(eg)
  deallocate(xerror)
  deallocate(hxerror)

  return

end subroutine diag_jacobi_iter_spin_c16
