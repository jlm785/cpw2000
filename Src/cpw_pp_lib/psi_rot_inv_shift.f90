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

!>  Copies the wave-function coeficients from one k-point
!>  to another k-point obtained by a rotation, followed by
!>  a faculative inversion and shift by a reciprocal lattice vector.
!>
!>  k(j) = sum_i mtrx_n(j,i)*k_ref(i) + kgshift(j) in lattice coordinates
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         21 November 2022.
!>  \copyright    GNU Public License v2

subroutine psi_rot_inv_shift(mtrx_n, tnp_n, kmap2, kgshift, neig,        &
    rkpt_ref, mtxd_ref,isort_ref,psi_ref,                                &
    rkpt, mtxd, isort, psi,                                              &
    ng, kgv,                                                             &
    mxdgve, mxddim, mxdbnd)

! Written 9 November 2022.  See related psi_translate.  JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxddim                          !<  array dimension for the hamiltonian
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands

  integer, intent(in)                ::  mtrx_n(3,3)                     !<  rotation matrix (in reciprocal lattice coordinates)
  real(REAL64), intent(in)           ::  tnp_n(3)                        !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector
  integer, intent(in)                ::  kmap2                           !<  if kmap2=1 the k-point was inverted (time reversal)
  integer, intent(in)                ::  kgshift(3)                      !<  reciprocal lattice vector shift between k-vectors

  integer, intent(in)                ::  neig                            !<  number of eigenvectors (requested on input, modified by degeneracies on output)
  real(REAL64), intent(in)           ::  rkpt_ref(3)                     !<  reference k-point (used for check)
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point (used for check)

  integer, intent(in)                ::  mtxd_ref                        !<  dimension of the hamiltonian (reference)
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian

  complex(REAL64), intent(in)        ::  psi_ref(mxddim,mxdbnd)          !<  component j of eigenvector i (reference)
  integer, intent(in)                ::  isort_ref(mxddim)               !<  g-vector associated with row/column i of hamiltonian (reference)
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

! output

  complex(REAL64), intent(out)       ::  psi(mxddim,mxdbnd)              !<  component j of eigenvector i (converted)

! local variables

  integer    ::  n1m, n2m, n3m
  integer    ::  kgv_rot(3)
  real(REAL64)     ::  xp
  real(REAL64)     ::  diff


! local allocatable arrays

  integer, allocatable               ::  isofkg(:,:,:)                   !  reverse mapping of rotated kgv(isort_ref)
  integer, allocatable               ::  ib(:)                           !  reverse mapping
  complex(REAL64), allocatable       ::  phase(:)                        !  phase from partial translation

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter ::  C_UM = cmplx(UM,ZERO,REAL64)
  complex(REAL64), parameter ::  C_I = cmplx(ZERO,UM,REAL64)
  real(REAL64), parameter    ::  EPS = 1.0E-8_REAL64

! counters

  integer ::  i, j, k, n, m


! The basis vectors should be the same for kgshift=0, but their ordering
! as described of the isort may not be the same.   So it is always necessary
! to make them compatible.

! paranoid check.  Guards against confusion in the relevant symmetry operation.

  do j = 1,3
    diff = ZERO
    do k = 1,3
      diff = diff + mtrx_n(j,k)*rkpt_ref(k)
    enddo
    if(kmap2 == 1) then
      diff = -diff - rkpt(j)
    else
      diff = diff - rkpt(j)
    endif
    diff = diff + kgshift(j)
    if(abs(diff) > EPS) then
      write(6,*)
      write(6,*) '   stopped in psi_rot_inv_shift'
      write(6,*) '   k-points not compatible with symmetry'
      write(6,*)

      stop

    endif
  enddo

! initializes output

  do n = 1,neig
  do i = 1,mtxd
    psi(i,n) = C_ZERO
  enddo
  enddo

! constructs the inverse mapping

  n1m = 0
  n2m = 0
  n3m = 0
  do m = 1,mtxd_ref

    do j = 1,3
      kgv_rot(j) = 0
      do k = 1,3
        kgv_rot(j) = kgv_rot(j) + mtrx_n(j,k)*kgv(k,isort_ref(m))
      enddo
    enddo
    if(abs(kgv_rot(1)) > n1m) n1m = abs(kgv_rot(1))
    if(abs(kgv_rot(2)) > n2m) n2m = abs(kgv_rot(2))
    if(abs(kgv_rot(3)) > n3m) n3m = abs(kgv_rot(3))
  enddo
  do i = 1,mtxd
    if(abs(kgv(1,isort(i))) > n1m) n1m = abs(kgv(1,isort(i)))
    if(abs(kgv(2,isort(i))) > n2m) n2m = abs(kgv(2,isort(i)))
    if(abs(kgv(3,isort(i))) > n3m) n3m = abs(kgv(3,isort(i)))
  enddo

  n1m = n1m + iabs(kgshift(1))
  n2m = n2m + iabs(kgshift(2))
  n3m = n3m + iabs(kgshift(3))


! constructs the inverse mapping

  allocate(isofkg(-n1m:n1m,-n2m:n2m,-n3m:n3m))

  isofkg(:,:,:) = 0
  do m = 1,mtxd_ref
!   finds rotated G-vector

    do j = 1,3
      kgv_rot(j) = 0
      do k = 1,3
        kgv_rot(j) = kgv_rot(j) + mtrx_n(j,k)*kgv(k,isort_ref(m))
      enddo
    enddo

    if(kmap2 == 1) then
      isofkg(-kgv_rot(1)-kgshift(1),-kgv_rot(2)-kgshift(2),-kgv_rot(3)-kgshift(3) ) = m
    else
      isofkg( kgv_rot(1)-kgshift(1), kgv_rot(2)-kgshift(2), kgv_rot(3)-kgshift(3) ) = m
    endif
  enddo

! mapping and phase

  allocate(ib(mtxd))
  allocate(phase(mtxd))

  do i = 1,mtxd

    ib(i) = isofkg( kgv(1,isort(i)),  kgv(2,isort(i)),  kgv(3,isort(i)) )

    if(ib(i) == 0) then
!     this should not occur, unless input is inconsistent.
      phase(i) = C_ZERO
    else
      xp = ZERO
      do k = 1,3
        xp = xp + kgv(k,isort(i))*tnp_n(k)
      enddo
      if(kmap2 == 1) xp = -xp
      phase(i) = C_UM*cos(xp) - C_I*sin(xp)
    endif

  enddo

! rotated functions

  do i = 1,mtxd
    if(ib(i) /= 0) then
      do n = 1,neig
        psi(i,n) = psi_ref(ib(i),n)*phase(i)
      enddo
    endif
  enddo

! possible inversion

  if(kmap2 == 1) then
    do n = 1,neig
    do i = 1,mtxd
      psi(i,n) = conjg(psi(i,n))
    enddo
    enddo
  endif

  deallocate(ib)
  deallocate(phase)

  deallocate(isofkg)

  return

end subroutine psi_rot_inv_shift
