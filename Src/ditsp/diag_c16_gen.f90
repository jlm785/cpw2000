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

!>     interface with lapack.
!>     diagonalizes the Hamiltonian matrix in a non-orthogonal basis

       subroutine diag_c16_gen(neig,ham,sm,ev,vec,mxdbnd,info)

!      written June 2012. jlm
!      modified by CLR for the generalized eigenvalue problem, March 2020
!      copyright inesc-mn/jose luis martins/Carlos Loia Reis
!

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

      
!      input

       integer, intent(in)               ::  mxdbnd                      !<  array dimension for number of bands

       integer, intent(in)               ::  neig                        !<  number of bands

       complex(REAL64), intent(in)       ::  ham(mxdbnd,mxdbnd)          !<  <Psi_i|H|Psi_j>
       complex(REAL64), intent(in)       ::  sm(mxdbnd,mxdbnd)           !<  <Psi_i|Psi_j>

!      output

       real(REAL64), intent(out)         ::  ev(mxdbnd)                  !<  eigenvalues
       complex(REAL64), intent(out)      ::  vec(mxdbnd,mxdbnd)          !<  eigenvector

       integer, intent(out)               ::  info                       !<  if info /=0 subroutine returned with error

!      local variables

       complex(REAL64), allocatable    :: work(:)
       complex(REAL64), allocatable    :: S(:,:)

       real(REAL64), allocatable    :: rwork(:)
       integer, allocatable    :: iwork(:)

       integer    ::  lwork,liwork,lrwork
       integer itype

!      counters

       integer i,j
        
       allocate(S(mxdbnd,mxdbnd))

       do i=1,neig
       do j=1,neig
         vec(j,i) = ham(j,i)
         S(j,i) = sm(j,i)
       enddo
       enddo

!      finds the dimension of work arrays

       allocate(work(2),rwork(2),iwork(2))

       lwork = -1
       liwork = -1
       lrwork = -1
       
       itype = 1
       
       call zhegvd (itype, 'V', 'L', neig, vec, mxdbnd, S, mxdbnd, ev,   &
     &      work, lwork, rwork, lrwork, iwork, liwork, info )

       if( info /= 0) then
          write(6,*)
          write(6,*)'    ERROR    diag_c16_gen FAILED A, info = ',info
          write(6,*)

          return

       endif
       
       lwork = int( work( 1 ) ) 
       lrwork = int( rwork( 1 ) ) 
       liwork = iwork( 1 )

       deallocate(work,rwork,iwork)

       allocate(work(lwork),rwork(lrwork),iwork(liwork))

      
       itype = 1

      call zhegvd (itype, 'V', 'L', neig, vec, mxdbnd, S, mxdbnd, ev,    &
     &     work, lwork, rwork, lrwork, iwork, liwork, info )

       if( info /= 0) then
          write(6,*)
          write(6,*)'     ERROR   diag_c16_gen  FAILED    FAILED ', info
          write(6,*)

          return
          
       endif

       deallocate(S)

       return

       end subroutine diag_c16_gen
