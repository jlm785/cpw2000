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

!>  interface with lapack.
!>  diagonalizes a real symmetric matrix
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         29 January 2023, 28 October 2026.
!>  \copyright    GNU Public License v2

subroutine diag_r8(neig, ham, ev, vec, mxdbnd, info)

! Adapted from diag_c16
! dsyev instead of dsyevd (MKL bug). 28 October 2026. JLM+claude

  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  neig                            !<  number of bands (matrix size) ( <= mxdbnd)
  real(REAL64), intent(in)           ::  ham(mxdbnd,mxdbnd)              !<  <Psi_i|H|Psi_j>

! output

  real(REAL64), intent(out)          ::  ev(mxdbnd)                      !<  eigenvalues
  real(REAL64), intent(out)          ::  vec(mxdbnd,mxdbnd)              !<  eigenvector

  integer, intent(out)               ::  info                            !<  if info /=0 subroutine returned with error

! local variables and arrays

  integer    ::  lwork
  real(REAL64), allocatable    :: work(:)

! counters

  integer i,j

  do i=1,neig
  do j=1,neig
    vec(j,i) = ham(j,i)
  enddo
  enddo

! finds the dimension of work array

  allocate(work(2))

  lwork = -1
  call dsyev( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork, info )

  if( info /= 0) then
     write(6,*)
     write(6,*)'    ERROR    diag_r8 FAILED A, info = ',info
     write(6,*)

     deallocate(work)

     return

  endif

  lwork = max(1,3*neig-1,int( work( 1 ) ))

  deallocate(work)

  allocate(work(lwork))

  call dsyev( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork, info )

  deallocate(work)

  if( info /= 0) then
     write(6,*)
     write(6,*)'    ERROR    diag_r8 FAILED B, info = ',info
     write(6,*)

     return

  endif

  return

end subroutine diag_r8
