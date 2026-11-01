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
!>  Uses dsyevd (divide and conquer), the fastest for all eigenvectors.
!>  The code for dsyev (QR, slower) and dsyevr (MRRR) is kept in
!>  comments at the end.  With MKL, dsyevd and dsyevr may divide by
!>  zero legitimately, so they stop if the code is compiled with
!>  trapping of floating point exceptions (ifx -fpe0).
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         29 January 2023, 1 November 2026.
!>  \copyright    GNU Public License v2

subroutine diag_r8(neig, ham, ev, vec, mxdbnd, info)

! Adapted from diag_c16
! dsyev instead of dsyevd (MKL bug). 28 October 2026. JLM+claude
! dsyevd again, dsyev and dsyevr in comments. 1 November 2026. JLM+claude

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
  real(REAL64), allocatable       :: work(:)
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
  call dsyevd( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork,            &
      iwork, liwork, info )

  if( info /= 0) then
     write(6,*)
     write(6,*)'    ERROR    diag_r8 FAILED A, info = ',info
     write(6,*)

     deallocate(work,iwork)

     return

  endif

  lwork = int( work( 1 ) )
  liwork = iwork( 1 )

  deallocate(work,iwork)

  allocate(work(lwork),iwork(liwork))

  call dsyevd( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork,            &
      iwork, liwork, info )

  deallocate(work,iwork)

  if( info /= 0) then
     write(6,*)
     write(6,*)'    ERROR    diag_r8 FAILED B, info = ',info
     write(6,*)

     return

  endif

  return

! Alternative with dsyev.  Replace the code from "finds the dimension
! of work arrays" to the end, there is no iwork.
!
!  allocate(work(2))
!
!  lwork = -1
!  call dsyev( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork, info )
!
!  if( info /= 0) then
!     write(6,*)
!     write(6,*)'    ERROR    diag_r8 FAILED A, info = ',info
!     write(6,*)
!
!     deallocate(work)
!
!     return
!
!  endif
!
!  lwork = max(1,3*neig-1,int( work( 1 ) ))
!
!  deallocate(work)
!
!  allocate(work(lwork))
!
!  call dsyev( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork, info )
!
!  deallocate(work)
!
!  if( info /= 0) then
!     write(6,*)
!     write(6,*)'    ERROR    diag_r8 FAILED B, info = ',info
!     write(6,*)
!
!     return
!
!  endif
!
!  return

! Alternative with dsyevr.  It destroys the input matrix, so ham is copied
! to a and not to vec, and the eigenvectors are returned in vec.
! Extra declarations:
!
!  real(REAL64), allocatable       :: a(:,:)                          !  copy of ham, destroyed by dsyevr
!  integer, allocatable            :: isuppz(:)                       !  support of eigenvectors
!
!  integer         ::  m                                              !  number of eigenvalues found
!  real(REAL64)    ::  vl, vu, abstol
!  integer         ::  il, iu
!
!  real(REAL64), parameter  :: ZERO = 0.0_REAL64
!
! Code (replaces everything after the declarations):
!
!  allocate(a(neig,neig))
!
!  do i=1,neig
!  do j=1,neig
!    a(j,i) = ham(j,i)
!  enddo
!  enddo
!
!  allocate(isuppz(2*max(1,neig)))
!
!! all eigenvalues, vl, vu, il, iu not referenced
!
!  vl = ZERO
!  vu = ZERO
!  il = 1
!  iu = neig
!  abstol = ZERO
!
!  allocate(work(2),iwork(2))
!
!  lwork = -1
!  liwork = -1
!  call dsyevr( 'V', 'A', 'L', neig, a, neig, vl, vu, il, iu, abstol,     &
!      m, ev, vec, mxdbnd, isuppz, work, lwork, iwork, liwork, info )
!
!  if( info /= 0) then
!     write(6,*)
!     write(6,*)'    ERROR    diag_r8 FAILED A, info = ',info
!     write(6,*)
!
!     deallocate(work,iwork)
!     deallocate(a,isuppz)
!
!     return
!
!  endif
!
!  lwork = max(1,26*neig,int( work( 1 ) ))
!  liwork = max(1,10*neig,iwork( 1 ))
!
!  deallocate(work,iwork)
!
!  allocate(work(lwork),iwork(liwork))
!
!  call dsyevr( 'V', 'A', 'L', neig, a, neig, vl, vu, il, iu, abstol,     &
!      m, ev, vec, mxdbnd, isuppz, work, lwork, iwork, liwork, info )
!
!  deallocate(work,iwork)
!  deallocate(a,isuppz)
!
!  if( info /= 0) then
!     write(6,*)
!     write(6,*)'    ERROR    diag_r8 FAILED B, info = ',info
!     write(6,*)
!
!     return
!
!  endif
!
!  if( m /= neig) then
!     write(6,*)
!     write(6,*)'    ERROR    diag_r8 FAILED C, found ',m,' of ',neig,' eigenvalues'
!     write(6,*)
!
!     info = -100
!
!     return
!
!  endif
!
!  return

end subroutine diag_r8
