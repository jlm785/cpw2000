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

!>  Finds improved guess eigen-vector by Jacobian relaxation and
!>  diagonalizes hamiltonian in the initial+improved subspace.
!>
!>  \author       Jose Luis Martins
!>  \version      5.09
!>  \date         April 2019. 10 December 2023.
!>  \copyright    GNU Public License v2

subroutine diag_jacobi_ritz_c16(mtxd, neig, eg,                          &
     psi, hpsi, nout,                                                    &
     njac, nritz, epsdeg,                                                &
     ng, kgv,                                                            &
     ekpg, isort, hdiag, vscr, kmscr,                                    &
     anlga, xnlkb, nanl, lnewanl,                                        &
     mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

! Adapted April-May 2019. JLM
! Adapted fo 5.0X, 10 December 2023.


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !< array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !< array dimension for number of bands
  integer, intent(in)                ::  mxdanl                          !< array dimension of number of projectors
  integer, intent(in)                ::  mxdgve                          !< array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !< array dimension of vscr

  integer, intent(in)                ::  mtxd                            !< wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !< number of wavefunctions

  integer, intent(in)                ::  njac                            !< number of jacobi iterations
  integer, intent(in)                ::  nritz                           !< number of ritz iterations
  real(REAL64), intent(in)           ::  epsdeg                          !< two eigenvalues are considered quasi degenerate if |e_i - e_j| < eps

  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !< length of k+g-vector of row/column i
  integer, intent(in)                ::  isort(mxddim)                   !< g-vector associated with row/column i of hamiltonian
  real(REAL64), intent(in)           ::  hdiag(mxddim)                   !< hamiltonian diagonal

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !< screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !< max value of kgv(i,n) used for the potential fft mesh

  integer, intent(in)                ::  ng                              !< total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !< i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nanl                            !< number of projectors
  complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)            !< Kleinman-Bylander projectors
  real(REAL64), intent(in)           ::  xnlkb(mxdanl)                   !< Kleinman-Bylander normalization

! output

  complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)             !< H |psi>
  real(REAL64), intent(out)          ::  eg(mxdbnd)                      !< <psi| H |psi>
  integer, intent(out)               ::  nout                            !< number of eigenvectors on output  neig <= nout <= mxdbnd, if negative error in the diagonalization.

! input and output

  complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)              !< |psi> on input, improved |psi> on output
  logical, intent(inout)             ::  lnewanl                         !< indicates that anlga has been recalculated (not used in default implementation)

! local allocatable arrays

  complex(REAL64), allocatable       ::  bas(:,:)
  complex(REAL64),allocatable        ::  hbas(:,:)

  integer, allocatable               ::  irow(:)                         !  if irow(i)=0 the vector is linearly dependent on the others and the corresponding xvec is zero

! local variables

  integer           ::  nbas
  integer           ::  mxdorb
  integer           ::  ired

! counters

  integer     ::  j, n

! work arrays.

  allocate(bas(mxddim,2*mxdbnd))
  allocate(hbas(mxddim,2*mxdbnd))

  do j = 1,nritz

    call zlacpy(' ', mtxd, neig, psi, mxddim, bas, mxddim)

! Jacobi relaxation

    call diag_jacobi_iter_c16(mtxd, neig, bas, hbas, njac,               &
        bas(1, mxdbnd+1), hbas(1, mxdbnd+1),                             &
        ng, kgv,                                                         &
        ekpg, isort, hdiag, vscr, kmscr,                                 &
        anlga, xnlkb, nanl, lnewanl,                                     &
        mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)
  if(mxdbnd > neig) then
    do n = 1,neig
      call zcopy(mtxd, bas(1,mxdbnd+n), 1, bas(1,neig+n), 1)
      call zcopy(mtxd, hbas(1,mxdbnd+n), 1, hbas(1,neig+n), 1)
    enddo
  endif

! orthogonalization.  Could use SVD...

    allocate(irow(2*neig))

    call grsch_loop_c16(.TRUE., bas, hbas, mtxd, neig, 2*neig, irow, mxddim)

    call grsch_shift_c16(bas, mtxd, neig, 2*neig, irow, ired, mxddim)

    call grsch_shift_c16(hbas, mtxd, neig, 2*neig, irow, ired, mxddim)

    deallocate(irow)

    nbas = 2*neig - ired
    mxdorb = 2*mxdbnd

! diagonalization

!     call diag_psi_hsq_psi_c16(mtxd, neig, eg, psi, hpsi, nout,           &
!         nbas, bas, hbas, elambda, epsdeg,                                &
!         mxddim, mxdbnd, mxdorb)

    call diag_psi_h_psi_c16(mtxd, neig, eg, psi, hpsi, nout, .TRUE.,     &
        nbas, bas, hbas,  epsdeg,                                        &
        mxddim, mxdbnd, mxdorb)

  enddo

  deallocate(hbas)
  deallocate(bas)

  return

end subroutine diag_jacobi_ritz_c16
