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

!>  Interface with lapack.
!>  diagonalizes the Hamiltonian matrix.
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         June 2012, 28 October 2026.
!>  \copyright    GNU Public License v2

subroutine diag_c16(neig, ham, ev, vec, mxdbnd, info)

! Written June 2012. jlm
! Modified, documentation, January 2020. JLM
! Modified 7 January 2014, style. jlm
! Indentation. 28 September 2026. JLM+claude
! zheev instead of zheevd (MKL bug). 28 October 2026. JLM+claude

  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  neig                            !<  number of bands (matrix size) ( <= mxdbnd)
  complex(REAL64), intent(in)        ::  ham(mxdbnd,mxdbnd)              !<  <Psi_i|H|Psi_j>

! output

  real(REAL64), intent(out)          ::  ev(mxdbnd)                      !<  eigenvalues
  complex(REAL64), intent(out)       ::  vec(mxdbnd,mxdbnd)              !<  eigenvector

  integer, intent(out)               ::  info                            !<  if info /=0 subroutine returned with error

! local variables and arrays

  integer    ::  lwork
  complex(REAL64), allocatable    :: work(:)
  real(REAL64), allocatable       :: rwork(:)

! counters

  integer i,j

  do i=1,neig
  do j=1,neig
    vec(j,i) = ham(j,i)
  enddo
  enddo

! finds the dimension of work array

  allocate(work(2))
  allocate(rwork(max(1,3*neig-2)))

  lwork = -1
  call zheev( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork, rwork, info )

  if( info /= 0) then
     write(6,*)
     write(6,*)'    ERROR    diag_c16 FAILED A, info = ',info
     write(6,*)

     deallocate(work,rwork)

     return

  endif

  lwork = max(1,2*neig-1,int( work( 1 ) ))

  deallocate(work)

  allocate(work(lwork))

  call zheev( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork, rwork, info )

  deallocate(work,rwork)

  if( info /= 0) then
     write(6,*)
     write(6,*)'    ERROR    diag_c16 FAILED B, info = ',info
     write(6,*)

     return

  endif

  return

end subroutine diag_c16
