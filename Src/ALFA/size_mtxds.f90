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

!>     calculates the small matrix size for the iterative
!>     diagonalization subroutines.
!>     works correctly only if the elements of the diagonal
!>     of the hamiltonian are essentially in increasing order.

       subroutine size_mtxds(ipr,hdiag,neig,mtxd,mtxds)

!      Based on sizdit.f written 4 November 1989, and
!      modified 25 March 1999. jlm
!      Written 20 December 1989. jlm
!      Modified, 3n+20, September 10, 2015. JLM
!      Modified, documentation, January 2020. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  ipr                        !<  prints message if ipr=1
       integer, intent(in)                ::  mtxd                       !<  dimension of the hamiltonian
       real(REAL64), intent(in)           ::  hdiag(mtxd)                !<  diagonal of the hamiltonian
       integer, intent(in)                ::  neig                       !<  number of eigenvectors required

!      output

       integer, intent(out)               ::  mtxds                      !<  dimension of the small matrix

!      local variables

       real(REAL64)  ::  href

!      counters

       integer    ::  i

!      parameters

       real(REAL64), parameter ::  EPS = 0.000001_REAL64
       integer, parameter      ::  NFULL = 250                           !  hard coded parameter for full diagonalization


       if(mtxd < neig) then
         write(6,*)
         write(6,'("   stopped in size_mtxds:     matrix size is ",i6,   &
     &        " and you want to calculate",i6," states")') mtxd,neig
         write(6,'("   check cutoff energy!!")')

         stop

       endif
      
       if(mtxd < NFULL) then
         mtxds = mtxd
       elseif(3*neig+20 >= mtxd) then
         mtxds = mtxd
       else
         href = abs(2*hdiag(neig)-hdiag(1)) + EPS
         do i = 3*neig+20,mtxd
           mtxds = i
           if(hdiag(i) > href .and. (hdiag(i+1)-hdiag(i)) > EPS) exit
         enddo    
       endif

       if(ipr > 1) write(6,'("  small matrix size = ",i5)') mtxds

       return
       end subroutine size_mtxds
