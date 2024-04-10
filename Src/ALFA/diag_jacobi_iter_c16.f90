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

!> Performs a succession of Jacobi relaxations
!>
!>  \author       Jose Luis Martins
!>  \version      5.09
!>  \date         May 3 2019. 10 December 2023.
!>  \copyright    GNU Public License v2

subroutine diag_jacobi_iter_c16(mtxd, neig, psi, hpsi, njac,             &
      dpsi, hdpsi,                                                       &
      ng, kgv,                                                           &
      ekpg, isort, hdiag, vscr, kmscr,                                   &
      anlga, xnlkb, nanl, lnewanl,                                       &
      mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

! Adapted May 3 2019. JLM
! Adapted for 5.0X, 10 December 2023. JLM
! Prefix diag_, 17 March 2024. JLM


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions

  integer, intent(in)                ::  njac                            !<  number of jacobi iterations

  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  length of k+g-vector of row/column i
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  real(REAL64), intent(in)           ::  hdiag(mxddim)                   !<  hamiltonian diagonal

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nanl                            !<  number of projectors
  complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)            !<  Kleinman-Bylander projectors
  real(REAL64), intent(in)           ::  xnlkb(mxdanl)                   !<  Kleinman-Bylander normalization

  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  |psi> to be improved

! output

  complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)             !<  H |psi>
  complex(REAL64), intent(out)       ::  dpsi(mxddim,mxdbnd)             !<  correction to |psi>.  Used has work vector in subroutine
  complex(REAL64), intent(out)       ::  hdpsi(mxddim,mxdbnd)            !<  H |dpsi>

! input and output

  logical, intent(inout)             ::  lnewanl                         !<  indicates that anlga has been recalculated (not used in default implementation)

! local variables

  real(REAL64)            ::  xn
!  complex(REAL64)         ::  cn

! local allocatable arrays

  complex(REAL64), allocatable       ::  xerror(:,:)                     !  error vector  H |psi) - |psi><psi| H |psi> or respective jacobian relaxed vector
  complex(REAL64), allocatable       ::  hxerror(:,:)                    !  H | relaxed xerror >
  real(REAL64), allocatable          ::  eg(:)                           !  <psi| H |psi>

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! external functions

  complex(REAL64),external   :: zdotc

! counters

  integer   ::   n, j

  allocate(eg(mxdbnd))
  allocate(xerror(mxddim,mxdbnd))
  allocate(hxerror(mxddim,mxdbnd))


! apply H to initial eigenvectors

  call hk_psi_c16(mtxd, neig, psi, hpsi, lnewanl,                        &
      ng, kgv,                                                           &
      ekpg, isort, vscr, kmscr,                                          &
      anlga, xnlkb, nanl,                                                &
      mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

  do n = 1,neig
    call zcopy(mtxd, psi(1,n), 1, dpsi(1,n), 1)
  enddo

  do n = 1,neig
    call zcopy(mtxd, hpsi(1,n), 1, hdpsi(1,n), 1)
  enddo

! start relaxation

   WRITE(6,*)
   WRITE(6,*)

  do j = 1,njac

!   xerror is the error vector or jacobian relaxed error

    do n = 1,neig
      eg(n) = real(zdotc(mtxd, dpsi(1, n), 1, hdpsi(1, n), 1), REAL64)
      call zcopy(mtxd, hdpsi(1, n), 1, xerror(1, n), 1)
      call zaxpy(mtxd, cmplx(-eg(n), ZERO, REAL64),                      &
                           dpsi(1, n), 1, xerror(1, n), 1)


       XN = REAL(ZDOTC(MTXD, XERROR(1, N), 1, XERROR(1, N), 1), REAL64)

       WRITE(6,'("  J,N,XN,EG = ",2I5,3X,F15.8,3X,F15.8)') J,N,XN,EG(N)*27.212


    enddo

   WRITE(6,*)


    call rq_jac_c16(xerror, eg, hdiag, mtxd, neig,                       &
    mxddim, mxdbnd)

    call hk_psi_c16(mtxd, neig, xerror, hxerror, lnewanl,                &
        ng, kgv,                                                         &
        ekpg, isort, vscr, kmscr,                                        &
        anlga, xnlkb, nanl,                                              &
        mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)


    do n = 1,neig

!     adds

      call zaxpy(mtxd, C_UM, xerror(1, n), 1, dpsi(1, n), 1)
      call zaxpy(mtxd, C_UM, hxerror(1, n), 1, hdpsi(1, n), 1)
      xn = real(zdotc(mtxd, dpsi(1, n), 1, dpsi(1, n), 1), REAL64)

      WRITE(6,'("  XN NORM",F12.8 )') XN

      xn = UM/sqrt(xn)
      call zscal(mtxd, cmplx(xn, ZERO, REAL64), dpsi(1, n), 1)
      call zscal(mtxd, cmplx(xn, ZERO, REAL64), hdpsi(1, n), 1)

    enddo

  enddo

  return

end subroutine diag_jacobi_iter_c16
