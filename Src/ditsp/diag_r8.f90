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
!>  \date         29 January 2023.
!>  \copyright    GNU Public License v2

subroutine diag_r8(neig, ham, ev, vec, mxdbnd, info)

! Adapted from diag_c16

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

  integer    ::  lwork,liwork
  real(REAL64), allocatable    :: work(:)
  integer, allocatable            :: iwork(:)

! counters

  integer i,j

  do i=1,neig
  do j=1,neig
    vec(j,i) = ham(j,i)
  enddo
  enddo

! finds the dimension of work arrays

  allocate(work(2),iwork(2))

  lwork = -1
  liwork = -1
  call dsyevd( 'V', 'L', neig, vec, mxdbnd , ev, work,lwork,iwork,liwork, info )

  if( info /= 0) then
     write(6,*)
     write(6,*)'    ERROR    diag_r8 FAILED A, info = ',info
     write(6,*)

     return

  endif

  lwork = int( work( 1 ) )
  liwork = iwork( 1 )

  deallocate(work,iwork)

  allocate(work(lwork),iwork(liwork))

  call dsyevd( 'V', 'L', neig, vec, mxdbnd , ev, work,lwork,iwork,liwork, info )

  deallocate(work,iwork)

  if( info /= 0) then
     write(6,*)
     write(6,*)'    ERROR    diag_r8 FAILED B, info = ',info
     write(6,*)

     return

  endif

  return

end subroutine diag_r8
