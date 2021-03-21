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

!>     Given |psi> and H|psi> estimates the error in the 
!>     eigenvector estimates of the maximum error: max | H |psi> - |psi><psi|H|psi> |^2

       subroutine ditsp_error(xmax, neig, mtxd, psi, hpsi,               &
     & mxddim, mxdbnd)
  
!      Written 11 June 2020 from previous code. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.95


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands
       
       integer, intent(in)                ::  neig                       !<  number of eigenvectors
       integer, intent(in)                ::  mtxd                       !<  dimension of the hamiltonian

       complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)         !<  |psi> 
       complex(REAL64), intent(in)        ::  hpsi(mxddim,mxdbnd)        !<  H |psi> 

!      output

       real(REAL64), intent(out)          ::  xmax                       !<  maximum error: | H |psi> - |psi><psi|H|psi> |^2

!      local allocatable array

       complex(REAL64),allocatable        ::  xerror(:)

!      local variables

       real(REAL64)  ::  xm
       real(REAL64)  ::  egn

!      counters

       integer       ::  n

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64

!      external functions

       complex(REAL64),external   :: zdotc


       allocate(xerror(mtxd))

       xmax = ZERO
       do n = 1,neig
         egn = real(zdotc(mtxd,psi(:,n),1,hpsi(:,n),1),REAL64)
         call zcopy(mtxd,hpsi(:,n),1,xerror,1)
         call zaxpy(mtxd,cmplx(-egn,ZERO,REAL64),psi(:,n),1,xerror,1)
         xm = real(zdotc(mtxd,xerror,1,xerror,1),REAL64)
         if(xmax < xm) xmax = xm
       enddo

       deallocate(xerror)

       return
       end subroutine ditsp_error
