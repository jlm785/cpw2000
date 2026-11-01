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
!>  Uses zheevd (divide and conquer), the fastest for all eigenvectors.
!>  The code for zheev (QR, slower) and zheevr (MRRR) is kept in
!>  comments at the end.  With MKL, zheevd and zheevr may divide by
!>  zero legitimately, so they stop if the code is compiled with
!>  trapping of floating point exceptions (ifx -fpe0).
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         June 2012, 1 November 2026.
!>  \copyright    GNU Public License v2

subroutine diag_c16(neig, ham, ev, vec, mxdbnd, info)

! Written June 2012. jlm
! Modified, documentation, January 2020. JLM
! Modified 7 January 2014, style. jlm
! Indentation. 28 September 2026. JLM+claude
! zheev instead of zheevd (MKL bug). 28 October 2026. JLM+claude
! zheevd again, zheev and zheevr in comments. 1 November 2026. JLM+claude

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

  integer    ::  lwork,liwork,lrwork
  complex(REAL64), allocatable    :: work(:)
  real(REAL64), allocatable       :: rwork(:)
  integer, allocatable            :: iwork(:)

! counters

  integer i,j

  do i=1,neig
  do j=1,neig
    vec(j,i) = ham(j,i)
  enddo
  enddo

! finds the dimension of work arrays

  allocate(work(2),rwork(2),iwork(2))

  lwork = -1
  liwork = -1
  lrwork = -1
  call zheevd( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork, rwork,     &
      lrwork, iwork, liwork, info )

  if( info /= 0) then
     write(6,*)
     write(6,*)'    ERROR    diag_c16 FAILED A, info = ',info
     write(6,*)

     deallocate(work,rwork,iwork)

     return

  endif

  lwork = int( work( 1 ) )
  lrwork = int( rwork( 1 ) )
  liwork = iwork( 1 )

  deallocate(work,rwork,iwork)

  allocate(work(lwork),rwork(lrwork),iwork(liwork))

  call zheevd( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork, rwork,     &
      lrwork, iwork, liwork, info )

  deallocate(work,rwork,iwork)

  if( info /= 0) then
     write(6,*)
     write(6,*)'    ERROR    diag_c16 FAILED B, info = ',info
     write(6,*)

     return

  endif

  return

! Alternative with zheev.  Replace the code from "finds the dimension
! of work arrays" to the end, rwork has a fixed size and there is no iwork.
!
!  allocate(work(2))
!  allocate(rwork(max(1,3*neig-2)))
!
!  lwork = -1
!  call zheev( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork, rwork, info )
!
!  if( info /= 0) then
!     write(6,*)
!     write(6,*)'    ERROR    diag_c16 FAILED A, info = ',info
!     write(6,*)
!
!     deallocate(work,rwork)
!
!     return
!
!  endif
!
!  lwork = max(1,2*neig-1,int( work( 1 ) ))
!
!  deallocate(work)
!
!  allocate(work(lwork))
!
!  call zheev( 'V', 'L', neig, vec, mxdbnd , ev, work, lwork, rwork, info )
!
!  deallocate(work,rwork)
!
!  if( info /= 0) then
!     write(6,*)
!     write(6,*)'    ERROR    diag_c16 FAILED B, info = ',info
!     write(6,*)
!
!     return
!
!  endif
!
!  return

! Alternative with zheevr.  It destroys the input matrix, so ham is copied
! to a and not to vec, and the eigenvectors are returned in vec.
! Extra declarations:
!
!  complex(REAL64), allocatable    :: a(:,:)                          !  copy of ham, destroyed by zheevr
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
!  allocate(work(2),rwork(2),iwork(2))
!
!  lwork = -1
!  liwork = -1
!  lrwork = -1
!  call zheevr( 'V', 'A', 'L', neig, a, neig, vl, vu, il, iu, abstol,     &
!      m, ev, vec, mxdbnd, isuppz, work, lwork, rwork, lrwork,            &
!      iwork, liwork, info )
!
!  if( info /= 0) then
!     write(6,*)
!     write(6,*)'    ERROR    diag_c16 FAILED A, info = ',info
!     write(6,*)
!
!     deallocate(work,rwork,iwork)
!     deallocate(a,isuppz)
!
!     return
!
!  endif
!
!  lwork = max(1,2*neig,int( work( 1 ) ))
!  lrwork = max(1,24*neig,int( rwork( 1 ) ))
!  liwork = max(1,10*neig,iwork( 1 ))
!
!  deallocate(work,rwork,iwork)
!
!  allocate(work(lwork),rwork(lrwork),iwork(liwork))
!
!  call zheevr( 'V', 'A', 'L', neig, a, neig, vl, vu, il, iu, abstol,     &
!      m, ev, vec, mxdbnd, isuppz, work, lwork, rwork, lrwork,            &
!      iwork, liwork, info )
!
!  deallocate(work,rwork,iwork)
!  deallocate(a,isuppz)
!
!  if( info /= 0) then
!     write(6,*)
!     write(6,*)'    ERROR    diag_c16 FAILED B, info = ',info
!     write(6,*)
!
!     return
!
!  endif
!
!  if( m /= neig) then
!     write(6,*)
!     write(6,*)'    ERROR    diag_c16 FAILED C, found ',m,' of ',neig,' eigenvalues'
!     write(6,*)
!
!     info = -100
!
!     return
!
!  endif
!
!  return

end subroutine diag_c16
