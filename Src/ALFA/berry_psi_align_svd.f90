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

!>  Aligns a subspace of wave-vectors at one k-point with another
!>  subspace at another k-point.
!>  See appendix B of Vanderbilt's Berry Phase book.
!>  Can also be used to multiply by inverse overlap matrix.
!>  Check berry_psi_align_one for aligning band by band.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         9 January 2023.
!>  \copyright    GNU Public License v2

subroutine berry_psi_align_svd(neig, mtxd, psi_ref, psi, lsigma,         &
     mxddim, mxdbnd)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddim                          !<  array dimension for the hamiltonian
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands

  integer, intent(in)                ::  neig                            !<  number of eigenvectors (requested on input, modified by degeneracies on output)
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  complex(REAL64), intent(in)        ::  psi_ref(mxddim,mxdbnd)          !<  component j of eigenvector i (reference)

  logical, intent(in)                ::  lsigma                          !<  multiply by 1/singval (inverse of overlap) or not (optimal alignment).  Default should be .FALSE.

! input and output

  complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)              !<  component j of eigenvector i (only change of global phase)

! local allocatable arrays

  complex(REAL64), allocatable       ::  ovlp(:,:)                       !  overlap matrix

  real(REAL64), allocatable          ::  singval(:)                      !  singular values
  complex(REAL64), allocatable       ::  u(:,:)                          !  left singular vectors (not used)
  complex(REAL64), allocatable       ::  vt(:,:)                         !  hermitian of v (right singular vectors)

  complex(REAL64), allocatable       ::  work(:)
  real(REAL64), allocatable          ::  rwork(:)
  integer, allocatable               ::  iwork(:)

  complex(REAL64), allocatable       ::  psi_tmp(:,:)

! local variables

  integer                 ::  info, lwork, lrwork

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter ::  C_UM = cmplx(UM,ZERO,REAL64)
  complex(REAL64), parameter ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  real(REAL64), parameter    ::  EPS = 1.0E-2_REAL64

! counters

  integer      ::  n, i, j

  COMPLEX(REAL64), EXTERNAL   ::  ZDOTC



  allocate(ovlp(neig,neig))

  call zgemm('c', 'n', neig, neig, mtxd, C_UM, psi_ref, mxddim, psi,     &
                   mxddim, C_ZERO, ovlp, neig)

! diagonal should be "positive"

  do n = 1,neig
    if(real(ovlp(n,n), REAL64) < -EPS) then
      ovlp(:,n) = - ovlp(:,n)
      psi(:,n) = -psi(:,n)
    endif
  enddo


  allocate(singval(neig))
  allocate(u(neig,neig))
  allocate(vt(neig,neig))
  allocate(work(1))
  allocate(rwork(1))
  allocate(iwork(8*neig))

  call zgesdd('A', neig, neig, ovlp, neig, singval, u, neig, vt ,neig,   &
         work, -1, rwork, iwork, info)

  if(info /= 0) then
    write(6,*) '    berry_psi_align_svd:  error reported by zgesdd (1), info = ',info

    stop

  endif

  lwork = nint(real(work(1)))
  deallocate(work)
  allocate(work(lwork))
  lrwork = 5*neig*neig + 7*neig
  deallocate(rwork)
  allocate(rwork(lrwork))

  call zgesdd('A', neig, neig, ovlp, neig, singval, u, neig, vt, neig,   &
        work, lwork, rwork, iwork, info)


  if(info /= 0) then
    write(6,*) '    berry_psi_align_svd:  error reported by zgesdd(2), info = ',info

    stop

  endif

  deallocate(work)
  deallocate(iwork)
  deallocate(rwork)

! difference between inverse overlap and optimal alignment
! uncomment and add lsigma to input parameters if you want this feature.

  if(lsigma) then
    do i = 1,neig
    do j = 1,neig
      vt(i,j) = (UM/singval(i))*vt(i,j)
    enddo
    enddo
  endif

  deallocate(singval)

! recycles ovlp

  call zgemm('n', 'n', neig, neig, neig, C_UM, u, neig, vt, neig,        &
                   C_ZERO, ovlp, neig)

  deallocate(u)
  deallocate(vt)

  allocate(psi_tmp(mxddim,mxdbnd))

  psi_tmp(:,:) = C_ZERO

  call zgemm('n', 'n', mtxd, neig, neig, C_UM, psi, mxddim, ovlp, neig,  &
                   C_ZERO, psi_tmp, mxddim)

  psi(:,:) = psi_tmp(:,:)

  deallocate(psi_tmp)
  deallocate(ovlp)

  return

end subroutine berry_psi_align_svd
